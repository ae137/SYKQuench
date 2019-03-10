/* main_quench.cpp: Numerical simulation of non-equilibrium dynamics of certain Majorana fermions
 * 
 * Copyright (C) 2017 Andreas Eberlein
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <sstream>
#include <boost/math/special_functions.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;

typedef complex<double> cmplx;

inline double pow2(double x) {
    return x * x;
}

inline cmplx pow2(cmplx z) {
	return z * z;
}

inline cmplx pow3(cmplx z) {
	return z * z * z;
}

inline cmplx pow4(cmplx z) {
	const cmplx z2 = pow2(z);
	return z2 * z2;
}

inline double GenerateWeight(const int i, const int i_begin, const int i_end) {
	/***************************************************************************
	 * Generates the weights for numerical integration using Trapezoidal or
	 * Simpson rule for different numbers of subintervals
	 * i		Index numbering points (running over [i_begin, i_end])
	 * i_begin	Index of first point
	 * i_end	Index of last point
	 **************************************************************************/
	return (i == i_begin || i == i_end) ? 0.5 : 1.;
}

inline cmplx iGgtr_random_hopping_model_Re(const double t, const double J2) {
	/***************************************************************************
	 * Compute real part of bare greater propgator i G^>_0:
	 * t		first time argument
	 * J2		Hopping in random hopping model
	 **************************************************************************/
	cmplx result;
	if (t == 0 || J2 == 0) { result.real(0.5); result.imag(0.); }
	else {
		result.real(boost::math::cyl_bessel_j(1, 2.*J2*t) / (2.*J2*t));
		result.imag(0.);		// Note: For the imaginary part, we don't have a library function at the moment
	}
	return result;
}

cmplx iGgtr_random_hopping_model_Numerical(const double t, const double J2, const int N) {
	/***************************************************************************
	 * Compute the Green's function iG^>(t) for the random hopping model from the
	 * spectral function
	 * t		time
	 * J2		random hopping amplitude
	 * N		Number of subintervals in the interval [0, 2J2]
	 **************************************************************************/
	const double h = (2. * J2) / N;
	cmplx sum = 0;
	for (int i = 0; i <= N; ++i) {
		const double omega = i * h;
		const cmplx exponential_factor = {cos(omega * t), -sin(omega * t)};
		sum += exponential_factor * sqrt(1 - pow2(omega / (2*J2))) * GenerateWeight(i, 0, N);
	}
	return sum * h / (M_PI * J2);
}

class Coupling {
private:
	double J_init, J_final;
	int Quench_time;
public:
	Coupling(const double j_init, const double j_final, const int quench_time) : J_init(j_init), J_final(j_final), Quench_time(quench_time) {}
	double operator()(const int time) const { return time > Quench_time ? J_final : J_init; }
	double get_J_init() const { return J_init; }
	double get_J_final() const { return J_final; }
};

cmplx iGK(const int m, const double omega, const int N, const double h, const boost::numeric::ublas::matrix<cmplx>& iG_gtr, const int index_max = 0) {
	const double dtau = 2. * h;
	cmplx result = 0;
	const int max_index = (index_max == 0) ? min(m, N - m - 1) : index_max;
	for (int n = -max_index; n <= max_index; ++n) {
		const double tau = dtau * n;
		const cmplx exp_factor = {cos(omega * tau), sin(omega * tau)};
		result += exp_factor * (iG_gtr(m + n, m - n) - iG_gtr(m - n, m + n)) * GenerateWeight(n, -max_index, max_index);
	}
	return result * dtau;
}

cmplx iGR(const int m, const double omega, const int N, const double h, const boost::numeric::ublas::matrix<cmplx>& iG_gtr, const int index_max = 0) {
	const double dtau = 2. * h;
	cmplx result = 0;
	const int max_index = (index_max == 0) ? min(m, N - m - 1) : index_max;
	for (int n = 0; n <= max_index; ++n) {
		const double tau = dtau * n;
		const cmplx exp_factor = {cos(omega * tau), sin(omega * tau)};
		result += exp_factor * (iG_gtr(m + n, m - n) + iG_gtr(m - n, m + n)) * GenerateWeight(n, 0, max_index);
	}
	return result * dtau;
}

cmplx Energy(const int i, const double h, const boost::numeric::ublas::matrix<cmplx>& iG_gtr, const Coupling J2, const Coupling J4) {
	cmplx I = {0, 1};
	cmplx integral_1 = 0;
	cmplx integral_2 = 0;
	for (int n = 0; n <= i; ++n) {
		integral_1 += J2(n) * (pow2(iG_gtr(i, n)) - pow2(iG_gtr(n, i))) * GenerateWeight(n, 0, i);
		integral_2 += J4(n) * (pow4(iG_gtr(i, n)) - pow4(iG_gtr(n, i))) * GenerateWeight(n, 0, i);
	}
	return h * I * (-0.5 * J2(i) * integral_1 - 0.25 * J4(i) * integral_2);
}

// For use in Kadanoff-Baym equation:
inline cmplx iSigma_gtr(const int i, const int j, const boost::numeric::ublas::matrix<cmplx>& iG_gtr, const Coupling J2, const Coupling J4) {
	return J2(i) * J2(j) * iG_gtr(i, j) + J4(i) * J4(j) * pow3(iG_gtr(i, j));
}

// For use in Dyson equation
inline cmplx iSigma_gtr(const cmplx iG_gtr_t, const double J2, const double J4) {
	return pow2(J2) * iG_gtr_t + pow2(J4) * pow3(iG_gtr_t);
}

cmplx f_rhs_t(const int i, const int j, const boost::numeric::ublas::matrix<cmplx>& iG_gtr,
			  const boost::numeric::ublas::matrix<cmplx>& iS_gtr, const double h) {
	/***************************************************************************
	 * Computes right hand side of ODE for d/dt iG^>(t, t')
	 * 
	 * i		Index for time t
	 * j		Index for time t'
	 * iG_gtr	Reference to matrix containing the data for i G^>
	 * iS_gtr	Reference to matrix containing the data for i Sigma^>
	 * h		Step width for integrations
	 **************************************************************************/
			
	cmplx integrals_first_term = 0;
	cmplx integrals_second_term= 0;
			
	for (int m = 0; m <= i; ++m) {
		integrals_first_term += -(iS_gtr(i, m) + iS_gtr(m, i)) * iG_gtr(m, j) * GenerateWeight(m, 0, i);
	}
	
	if (j > 0) {
		for (int m = 0; m <= j; ++m) {
			integrals_second_term += iS_gtr(i, m) * (iG_gtr(m, j) + iG_gtr(j, m)) * GenerateWeight(m, 0, j);
		}
	}
	
	return (integrals_first_term + integrals_second_term) * h;
}

cmplx f_rhs_tp(const int i, const int j, const boost::numeric::ublas::matrix<cmplx>& iG_gtr,
			  const boost::numeric::ublas::matrix<cmplx>& iS_gtr, const double h) {
	/***************************************************************************
	 * Computes right hand side of ODE for d/dt' iG^>(t, t')
	 * 
	 * i		Index for time t
	 * j		Index for time t'
	 * iG_gtr	Reference to matrix containing the data for i G^>
	 * iS_gtr	Reference to matrix containing the data for i Sigma^>
	 * h		Step width for integrations
	 **************************************************************************/
	cmplx integrals_first_term = 0;
	cmplx integrals_second_term= 0;
	
	if (i > 0) {
		for (int m = 0; m <= i; ++m) {
			integrals_first_term += (iG_gtr(i, m) + iG_gtr(m, i)) * iS_gtr(m, j) * GenerateWeight(m, 0, i);
		}
	}
	
	for (int m = 0; m <= j; ++m) {
		integrals_second_term += -iG_gtr(i, m) * (iS_gtr(m, j) + iS_gtr(j, m)) * GenerateWeight(m, 0, j);
	}
	
	return (integrals_first_term + integrals_second_term) * h;
}

cmplx FourierTransform(const double omega, const int N, const boost::numeric::ublas::vector<double>& t, const boost::numeric::ublas::vector<cmplx>& Data) {
	/***************************************************************************
	 * Computes the Fourier transform of Data at frequency omega
	 *                     F(Data)(omega)
	 * 
	 * omega		Frequency
	 * N			Number of time points in t
	 * t			Vector with times
	 * Data			Data that should be Fourier transformed
	 **************************************************************************/
	const double h = t(1) - t(0);
	cmplx result = 0;
	for (int n = 0; n < N; ++n) {
		const cmplx exp_factor = {cos(omega * t(n)), sin(omega * t(n))};
		result += exp_factor * Data(n) * GenerateWeight(n, 0, N-1);
	}
	return result * h;
}

cmplx InvFourierTransform(const double t, const int N, boost::numeric::ublas::vector<double> omega, boost::numeric::ublas::vector<cmplx> Data) {
	/***************************************************************************
	 * Computes the inverse Fourier transform of Data at time t
	 *                     F(Data)(t)
	 * 
	 * t			Frequency
	 * N			Number of frequency points in omega
	 * omega		Vector with frequencies
	 * Data			Data that should be Fourier transformed
	 **************************************************************************/
	const double h = omega(1) - omega(0);
	cmplx result = 0;
	for (int n = 0; n < N; ++n) {
		const cmplx exp_factor = {cos(omega(n) * t), -sin(omega(n) * t)};
		result += exp_factor * Data(n) * GenerateWeight(n, 0, N-1);
	}
	return result * h / (2. * M_PI);
}

double HeavisideTheta(const double x) {
	if (x > 1e-14) { return 1; }
	else if (x < -1e-14) { return 0; }
	else { return 0.5; }
}

inline double nF(const double x, const double beta = 1000) {
	return 1. / (1. + exp(beta * x));
}

double f(double beta_guess, const int m_T, const double omega_max, 
         const boost::numeric::ublas::vector<double>& omega_vec, 
         const boost::numeric::ublas::matrix<cmplx>& iGR_val, 
         const boost::numeric::ublas::matrix<cmplx>& iGK_val) {
	// Helper function for determining beta by fitting tanh(0.5 beta omega) to iG^K / A.
	double f_beta = 0;
	for (unsigned int i = 0; i < omega_vec.size(); ++i) {
		if (abs(omega_vec(i)) < omega_max) {
			f_beta += pow2(tanh(0.5 * beta_guess * omega_vec(i)) - iGK_val(m_T, i).real() / (2 * iGR_val(m_T, i).real()));
		}
	}
	f_beta *= 0.5;
	return f_beta;
}

/***************************************************************************
 * Main function of main_quench
 * 
 * Main function for simulation of Kadanoff-Baym equations for the non-equilibrium dynamics of
 * Majorana fermions after a quantum quench
 * 
 * Program start:
 * ./main_quench N h beta_init J2_init J2_final J4_init J4_final type
 * 
 * Command line parameters:
 *      N           Number of points in time grid (should be odd)
 *      h           Time step for numerical solutions of differential equations
 *      beta_init   Temperature of the system before the quench
 *      J2_init     Value of quadratic coupling before quench
 *      J2_final    Value of quadratic coupling after quench
 *      J4_init     Value of quartic coupling before quench
 *      J4_final    Value of quartic coupling after quench
 *      type        String describing the quench type
 * 
 * Return value:
 *      None
 **************************************************************************/
int main(int argc, char **argv)
{
    // Programmaufruf
    // ./main_quench N h beta_init J2_init J2_final J4_init J4_final type

	// Setting up the initial data
    // Number of point in grid, should be odd so that 0 is part of the grid
	const int N = atoi(argv[1]);							// 5001;
    // Time step for integrations and differential equations
	const double h = atof(argv[2]);							// 4e-2;
	
	boost::numeric::ublas::vector<double> t(N);
    // This index for the time right before the quench
	const int quench_time = N/2;												
    // Note: t = 0 is at t(N/2) and t(N-1) = -t(0) for N odd
	for (int i = 0; i < N; ++i) t(i) = (-N/2 + i) * h;
	
	// Inverse initial temperature
	const double beta_init = atof(argv[3]); 				// 50;			
	
	const double J2_init = atof(argv[4]); 					// 0.03125;
	const double J2_final = atof(argv[5]);					// 0;
	const double J4_init = atof(argv[6]);					// 1;
	const double J4_final = atof(argv[7]);					// 1;
	string type_of_calculation = argv[8];					// "J2J4_J4";
    
	const Coupling J2(J2_init, J2_final, quench_time);
	const Coupling J4(J4_init, J4_final, quench_time);
	
	boost::numeric::ublas::matrix<cmplx> iG_init(N, N);
	boost::numeric::ublas::matrix<cmplx> iG_gtr(N, N);
	boost::numeric::ublas::matrix<cmplx> iS_gtr(N, N);
	
	ostringstream outfile_details;
	outfile_details << "_Type_" << type_of_calculation << "_J2init_" << J2.get_J_init() << "_J2final_" << J2.get_J_final()
			<< "_J4init_" << J4.get_J_init() << "_J4final_" << J4.get_J_final() << "_beta_" << beta_init << "_N_" << N << "_h_" << h << ".out";
			
	string outfile_energy_name = "diG_dt_Energy";
	outfile_energy_name.append(outfile_details.str());
	ofstream outfile_energy(outfile_energy_name, ios_base::out);
	
	string outfile_diagnostics_name = "Diagnostics";
	outfile_diagnostics_name.append(outfile_details.str());
	ofstream outfile_diagnostics(outfile_diagnostics_name, ios_base::out);
	
	chrono::steady_clock::time_point time_begin = chrono::steady_clock::now();
	
	// Computation of initial conditions
	// Solution of Dyson equation:
	const int N_Dyson = 2*(N-1) + 1;
	boost::numeric::ublas::vector<double> t_Dyson(N_Dyson);
	boost::numeric::ublas::vector<cmplx> iG_gtr_t(N_Dyson);
	boost::numeric::ublas::vector<cmplx> iG_gtr_t_temp(N_Dyson);
	boost::numeric::ublas::vector<cmplx> iSigma_R_t(N_Dyson);
	for (int i = 0; i < N_Dyson; ++i) { t_Dyson(i) = 2.*t(0) + h * i; }
	
	const int N_omega_Dyson = N;
	const double omega_max_Dyson = 6. * max(J2_init, J4_init);
	const double omega_min_Dyson = -omega_max_Dyson;
	boost::numeric::ublas::vector<double> omega_vec_Dyson(N_omega_Dyson);
	boost::numeric::ublas::vector<cmplx> iG_R_omega(N_omega_Dyson);
	boost::numeric::ublas::vector<cmplx> iG_gtr_omega(N_omega_Dyson);
	for (int i = 0; i < N_omega_Dyson; ++i) {
		omega_vec_Dyson(i) = omega_min_Dyson + (omega_max_Dyson - omega_min_Dyson) * i / (N_omega_Dyson - 1);
	}
	
	for (int l = 0; l < N_Dyson; ++l) {
// 		 iG_gtr_t(l) = iG_gtr_t_temp(l) = iGgtr_random_hopping_model_Re(t_Dyson(l), max(J2.get_J_init(), 0.01));
		 iG_gtr_t(l) = iG_gtr_t_temp(l) = iGgtr_random_hopping_model_Numerical(t_Dyson(l), max(J2.get_J_init(), 0.01), N);
	}
	
	int m, n;
	
	
	// Numerical solution of Dyson equation in order to obtain initial conditions:
	const int max_iter_Dyson = 100;
	int iter_Dyson = 0;
	double norm_Dyson = numeric_limits<double>::max();
	for (iter_Dyson = 0; iter_Dyson < max_iter_Dyson; ++iter_Dyson) {
		for (int i = 0; i < N_Dyson; ++i) {
			iSigma_R_t(i) = HeavisideTheta(t_Dyson(i)) * (iSigma_gtr(iG_gtr_t(i), J2.get_J_init(), J4.get_J_init()) + iSigma_gtr(iG_gtr_t(N_Dyson-1 - i), J2.get_J_init(), J4.get_J_init()));
		}
		// Fourier transformation and computation of retarded Green's function:
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) private(m)
#endif
		for (m = 0; m < N_omega_Dyson; ++m) {
			cmplx iSigma_R_omega = FourierTransform(omega_vec_Dyson(m), N_Dyson, t_Dyson, iSigma_R_t);
			cmplx i_omega = {0, omega_vec_Dyson(m)};
			iG_R_omega(m) = -1. / (i_omega - iSigma_R_omega);
			iG_gtr_omega(m) = (1 - nF(omega_vec_Dyson(m), beta_init)) * 2. * iG_R_omega(m).real();
		}
		
		// Inverse Fourier transformation:
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) private(m)
#endif
		for (m = 0; m < N_Dyson; ++m) {
			const double h_mix = 0.1;
			iG_gtr_t(m) = h_mix * InvFourierTransform(t_Dyson(m), N_omega_Dyson, omega_vec_Dyson, iG_gtr_omega) + (1 - h_mix) * iG_gtr_t(m);
		}
		
		if (iG_gtr_t(N_Dyson/2).real() < 0) {				// Make sure that we pick the solution where G^>(t, t) is correct (= 1/2)
// 			cerr << "# Switched sign of iG_gtr" << endl;
			for (m = 0; m < N_Dyson; ++m) { iG_gtr_t(m) *= -1; }
		}
	
		
		norm_Dyson = 0;
		// Computation of difference norm:
		for (int i = 0; i < N_Dyson; ++i) {
			norm_Dyson += abs(iG_gtr_t(i) - iG_gtr_t_temp(i));
			iG_gtr_t_temp(i) = iG_gtr_t(i);
		}
		norm_Dyson /= N_Dyson;
		
// 		cerr << iter_Dyson << '\t' << norm_Dyson << endl;
		
		if (norm_Dyson < 1e-7) { break; }		
	}
	
	if (norm_Dyson > 1e-5 && iter_Dyson == (max_iter_Dyson-1)) {
		cerr << "Iterative solution of Dyson equation failed to converge!" << endl;
		exit(1);
	}
	

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) collapse(2) private(m, n)
#endif
	for (m = 0; m <= quench_time; ++m) {
		for (n = 0; n <= quench_time; ++n) {
			iG_init(m, n) = iG_gtr_t(static_cast<int>(round((t(m) - t(n) - t_Dyson(0)) / h)));		// iGgtr_random_hopping_model_Re(t(m) - t(n), J2_init);
			iG_gtr(m, n) = iG_init(m, n);
			iS_gtr(m, n) = iSigma_gtr(m, n, iG_gtr, J2, J4);
		}
	}
	
	
	/***************************************************************************
	 * First step:
	 * - Trapezoidal rule for time integrals (which are 1d-integrals)
	 * - Euler method with predictor-corrector scheme for ODE
	 * - Initial condition is Green's function for random hopping model
	 * 
	 * Remark: t = 0 (i. e. t(N/2)) is the time step right before the quench
	 **************************************************************************/
	
	// Solution of differential equations using first-order explicit Euler scheme
	for (int i = quench_time; i < N-1; ++i) {		// quench_time = N/2
		iG_gtr(i+1, i+1) = 0.5;
		iS_gtr(i+1, i+1) = iSigma_gtr(i+1, i+1, iG_gtr, J2, J4);
		
		boost::numeric::ublas::vector<cmplx> iG_gtr_predictor_ip1_j(i+1);		// Stores the predictor values of iG_gtr(i+1, j)
		boost::numeric::ublas::vector<cmplx> iG_gtr_predictor_j_ip1(i+1);		// Stores the predictor values of iG_gtr(j, i+1)
		
		// Predictor step:
		// Note: In the predictor step, we don't access the elements (i+1, j) and (j, i+1) of iG_gtr on the right hand side!
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
		for (int j = 0; j <= i; ++j) {
		// Time step in the t-direction:
		// i -> t, j -> t'
			iG_gtr(i+1, j) = iG_gtr(i, j) + h * f_rhs_t(i, j, iG_gtr, iS_gtr, h);
			iG_gtr_predictor_ip1_j(j) = iG_gtr(i+1, j);
			iS_gtr(i+1, j) = iSigma_gtr(i+1, j, iG_gtr, J2, J4);	// Now stores self-energy for predictor
		}
		
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
		for (int j = 0; j <= i; ++j) {
		// Time step in the t'-direction:
		// i -> t', j -> t
			iG_gtr(j, i+1) = iG_gtr(j, i) + h * f_rhs_tp(j, i, iG_gtr, iS_gtr, h);
			iG_gtr_predictor_j_ip1(j) = iG_gtr(j, i+1);
			iS_gtr(j, i+1) = iSigma_gtr(j, i+1, iG_gtr, J2, J4);	// Now stores self-energy for predictor
		}
		
		boost::numeric::ublas::vector<cmplx> iG_gtr_temp_ip1_j(i+1);		// Stores iG_gtr(i+1, j) from after the last iteration
		boost::numeric::ublas::vector<cmplx> iG_gtr_temp_j_ip1(i+1);		// Stores iG_gtr(j, i+1) from after the last iteration
		
	// Corrector code starts here
		// Note: iG_gtr and iS_gtr at beginning of iteration contain the values of the predictor
		// Corrector step:
		for (int pc_iter = 0; pc_iter < 10; ++pc_iter) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
			for (int j = 0; j <= i; ++j) {
			// Time step in the t-direction:
			// i -> t, j -> t'
				iG_gtr_temp_ip1_j(j) = iG_gtr(i, j) + h * f_rhs_t(i+1, j, iG_gtr, iS_gtr, h);
			}
		
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
			for (int j = 0; j <= i; ++j) {
			// Time step in the t'-direction:
			// i -> t', j -> t
				iG_gtr_temp_j_ip1(j) = iG_gtr(j, i) + h * f_rhs_tp(j, i+1, iG_gtr, iS_gtr, h);
			}
			
			double norm_new = 0;
			// Determine how much the values change during one iteration:
			boost::numeric::ublas::vector<cmplx> iG_gtr_old_ip1_j(i+1);
			boost::numeric::ublas::vector<cmplx> iG_gtr_old_j_ip1(i+1);
			for (int j = 0; j <= i; ++j) {
				iG_gtr_old_ip1_j(j) = iG_gtr(i+1, j);
				iG_gtr_old_j_ip1(j) = iG_gtr(j, i+1);
				
				norm_new += abs(iG_gtr_temp_ip1_j(j) - iG_gtr_old_ip1_j(j)) + abs(iG_gtr_temp_j_ip1(j) - iG_gtr_old_j_ip1(j));
			}
			
			// Output some diagnostic information about iterative solution of self-consistency condition
			outfile_diagnostics << i << '\t' << pc_iter << '\t' << norm_new << endl;
			
			// Update propagator and self-energies:
			for (int j = 0; j <= i; ++j) {
				iG_gtr(i+1, j) = iG_gtr_temp_ip1_j(j);
				iS_gtr(i+1, j) = iSigma_gtr(i+1, j, iG_gtr, J2, J4);
				
				iG_gtr(j, i+1) = iG_gtr_temp_j_ip1(j);
				iS_gtr(j, i+1) = iSigma_gtr(j, i+1, iG_gtr, J2, J4);
			}
			if (norm_new < 1e-3) { break; }
		}
		
		for (int j = 0; j <= i; ++j) {
			iG_gtr(i+1, j) = 0.5 * (iG_gtr_predictor_ip1_j(j) + iG_gtr(i+1, j));
			iG_gtr(j, i+1) = 0.5 * (iG_gtr_predictor_j_ip1(j) + iG_gtr(j, i+1));
			
			iS_gtr(i+1, j) = iSigma_gtr(i+1, j, iG_gtr, J2, J4);
			iS_gtr(j, i+1) = iSigma_gtr(j, i+1, iG_gtr, J2, J4);
		}
		
		// Compute and output diGgtr_dt for diagnostic purposes:
		cmplx diGgtr_dt = (iG_gtr(i + 1, i) - iG_gtr(i - 1, i)) / (2*h);
		outfile_diagnostics << i << '\t' << -1 << '\t' << diGgtr_dt.imag() << '\t' << -diGgtr_dt.real() << endl;
			
	// Corrector code ends here
	}
	
	chrono::steady_clock::time_point time_end= chrono::steady_clock::now();
	
	// Output of results for time dependence of iG^>(t, t'):
	string outfile_results_name = "Results";
	outfile_results_name.append(outfile_details.str());
	ofstream outfile_t_tp(outfile_results_name, ios_base::out);	
	// Output of parameters:
	outfile_t_tp << "# h =\t" << h << endl;
	outfile_t_tp << "# N =\t" << N << endl;
	
	outfile_t_tp << "# J2init =\t" << J2.get_J_init() << endl;
	outfile_t_tp << "# J2final =\t" << J2.get_J_final() << endl;
	outfile_t_tp << "# J4init =\t" << J4.get_J_init() << endl;
	outfile_t_tp << "# J4final =\t" << J4.get_J_final() << endl;
	outfile_t_tp << "# beta_init =\t" << beta_init << endl;
	outfile_t_tp << "# Type of calculation:\t" << type_of_calculation << endl << endl;
	
	outfile_t_tp << "# t (1)" << '\t' << "t' (2)" << '\t' << "T (3)" << '\t' << "tau (4)" << '\t' << "Re iG_gtr (5)" << '\t' << "Im iG_gtr (6)" << '\t' << "Re iG_R (7)" << '\t' << "Im iG_R (8)" << endl;

	const int output_increment = 5;
	for (int i = 0; i < N; i += output_increment) {
		for (int j = 0; j < N; j += output_increment) {
			const double t_val = t(i);
			const double tp_val = t(j);
			const double T = 0.5 * (t(i) + t(j));
			const double tau = t(i) - t(j);
			const cmplx iG_R_ij = (tau >= 0) ? iG_gtr(i, j) + iG_gtr(j, i) : 0;
			outfile_t_tp << t_val << '\t' << tp_val << '\t' << T << '\t' << tau << '\t' << iG_gtr(i, j).real() << '\t' << iG_gtr(i, j).imag() << '\t' << iG_R_ij.real() << '\t' << iG_R_ij.imag() << endl;
		}
	}
	outfile_t_tp.close();
	
	// Output time that it took to compute the results
	outfile_diagnostics << endl;
	outfile_diagnostics << "# Time difference = " << chrono::duration_cast<chrono::seconds>(time_end - time_begin).count() << " s" << endl;
	outfile_diagnostics << "# Time difference = " << chrono::duration_cast<chrono::microseconds>(time_end - time_begin).count() << " mu s" << endl;
	
	outfile_diagnostics.close();
	
	boost::numeric::ublas::vector<cmplx> E_t(N);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
	for (int i = 1; i < N-1; ++i) {
		E_t(i) = Energy(i, h, iG_gtr, J2, J4);
	}
	
	outfile_energy << "# i (1)" << '\t' << "t(i) (2)" << '\t' << "diGgtr_dt.imag() (3)" << '\t' << "-diGgtr_dt.real() (4)" << '\t' << "E_t(i).real() (5)" << '\t' << "E_t(i).imag() (6)" << '\t' << "J2(i) (7)" << '\t' << "J4(i) (8)" << endl;
	
	// Output of energy, first test:
	for (int i = 1; i < N-1; ++i) {
		cmplx diGgtr_dt = (iG_gtr(i + 1, i) - iG_gtr(i - 1, i)) / (2*h);
		outfile_energy << i << '\t' << t(i) << '\t' << diGgtr_dt.imag() << '\t' << -diGgtr_dt.real() << '\t' << E_t(i).real() << '\t' << E_t(i).imag() << '\t' << J2(i) << '\t' << J4(i) << endl;
	}
	
	outfile_energy.close();
	
	// Computation and output of spectral function:
	
	string outfile_G_name = "iGR_iGK";
	outfile_G_name.append(outfile_details.str());
	ofstream outfile_G(outfile_G_name, ios_base::out);
	
	// Output of parameters:
	outfile_G << "# h =\t" << h << endl;
	outfile_G << "# N =\t" << N << endl;
	
	outfile_G << "# J2init =\t" << J2.get_J_init() << endl;
	outfile_G << "# J2final =\t" << J2.get_J_final() << endl;
	outfile_G << "# J4init =\t" << J4.get_J_init() << endl;
	outfile_G << "# J4final =\t" << J4.get_J_final() << endl;
	outfile_G << "# Type of calculation:\t" << type_of_calculation << endl << endl;
	
	outfile_G << "# T (1)" << '\t' << "omega (2)" << '\t' << "Re iGR(omega) = A(omega)/2 (3)" << '\t' << "Im iGR(omega) (4)" << '\t' << "Re iGK(omega) (5)" << '\t' << "Im iGK(omega) (6)";
	outfile_G << "Re Sigma^R(omega) (7)" << '\t' << "Im Sigma^R(omega) (8)" << '\t' << "beta_eff(T) (9)" << endl;
	
	const int N_omega = N;
	const double omega_max = 6. * max(max(J2_init, J2_final), max(J4_init, J4_final));
	const double omega_min = -omega_max;
	boost::numeric::ublas::vector<double> omega_vec(N_omega);
	for (int i = 0; i < N_omega; ++i) {
		omega_vec(i) = omega_min + (omega_max - omega_min) * i / (N_omega - 1);
	}
	
	boost::numeric::ublas::matrix<cmplx> iGK_val(N, N_omega);
	boost::numeric::ublas::matrix<cmplx> iGR_val(N, N_omega);
	
	int m_T, i_w;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1) collapse(2) private(m_T, i_w)
#endif
	for (m_T = N/4; m_T < 3*N/4; ++m_T) {
		for (i_w = 0; i_w < N_omega; ++i_w) {
			iGK_val(m_T, i_w) = iGK(m_T, omega_vec(i_w), N, h, iG_gtr, N/4);
			iGR_val(m_T, i_w) = iGR(m_T, omega_vec(i_w), N, h, iG_gtr, N/4);
		}
	}
	
	for (m_T = N/4; m_T < 3*N/4; m_T += output_increment) {
		for (i_w = 0; i_w < N_omega; i_w += (abs(omega_vec(i_w)) < max(J2_final, J4_final) ? 1 : output_increment)) {
			outfile_G << t(m_T) << '\t' << omega_vec(i_w) << '\t' << iGR_val(m_T, i_w).real() << '\t' << iGR_val(m_T, i_w).imag() << '\t';
			outfile_G << iGK_val(m_T, i_w).real() << '\t' << iGK_val(m_T, i_w).imag() << '\t';
			outfile_G << omega_vec(i_w) + (1. / iGR_val(m_T, i_w)).imag() << '\t' << -(1. / iGR_val(m_T, i_w)).real();
			
			if (omega_vec(i_w) == 0) {
// 				const double h_deriv = omega_vec(i_w + 1) - omega_vec(i_w);
// 				double beta_eff_T = -(iGK_val(m_T, i_w + 2).real() / iGR_val(m_T, i_w + 2).real()) + 8 * (iGK_val(m_T, i_w + 1).real() / iGR_val(m_T, i_w + 1).real())
// 							- 8 * (iGK_val(m_T, i_w - 1).real() / iGR_val(m_T, i_w - 1).real()) + (iGK_val(m_T, i_w - 2).real() / iGR_val(m_T, i_w - 2).real());
// 				beta_eff_T /= 12. * h_deriv;
				
				const double h_deriv = omega_vec(i_w + 2) - omega_vec(i_w);
				const int step = 2;
				double beta_eff_T = -(iGK_val(m_T, i_w + 2*step).real() / iGR_val(m_T, i_w + 2*step).real()) + 8 * (iGK_val(m_T, i_w + step).real() / iGR_val(m_T, i_w + step).real())
							- 8 * (iGK_val(m_T, i_w - step).real() / iGR_val(m_T, i_w - step).real()) + (iGK_val(m_T, i_w - 2*step).real() / iGR_val(m_T, i_w - 2*step).real());
				beta_eff_T /= 12. * h_deriv;
				
				outfile_G << '\t' << beta_eff_T;
				
				// Bestimmung von beta durch Fit von tanh(beta omega/2) über einen breiteren Frequenzbereich:
				// Minimierung einer Kostenfunktion durch Golden Section Search:
				const double rho = 0.5 * (3 - sqrt(5));
				
				const double omega_max_beta = 0.1 * max(J2_final, J4_final);
				
				const double beta_a_val = 0.00001;
				const double beta_d_val = beta_init;
				
				const double beta_b_val = beta_a_val + rho * (beta_d_val - beta_a_val);
				const double beta_c_val = beta_d_val - rho * (beta_d_val - beta_a_val);
				
				pair<double, double> a = {beta_a_val, f(beta_a_val, m_T, omega_max_beta, omega_vec, iGR_val, iGK_val)};
				pair<double, double> b = {beta_b_val, f(beta_b_val, m_T, omega_max_beta, omega_vec, iGR_val, iGK_val)};
				pair<double, double> c = {beta_c_val, f(beta_c_val, m_T, omega_max_beta, omega_vec, iGR_val, iGK_val)};
				pair<double, double> d = {beta_d_val, f(beta_d_val, m_T, omega_max_beta, omega_vec, iGR_val, iGK_val)};
				
				int iter_beta = 0;
				const int max_iter_beta = 100;
				const double epsilon = 1e-3;
				
				while (d.first - a.first > epsilon && iter_beta < max_iter_beta) {
					if (b.second < c.second) {
						const pair<double, double> temp = c;
						c = b;
						d = temp;
						const double b_val = a.first + rho * (d.first - a.first);
						b = {b_val, f(b_val, m_T, omega_max_beta, omega_vec, iGR_val, iGK_val)};
					} else {
						const pair<double, double> temp = b;
						b = c;
						a = temp;
						const double c_val = d.first - rho * (d.first - a.first);
						c = {c_val, f(c_val, m_T, omega_max_beta, omega_vec, iGR_val, iGK_val)};
					}
					++iter_beta;
				}
				
				outfile_G << '\t' << 0.5 * (a.first + d.first);
				
				
				
				
			}
			outfile_G << endl;
		}
		outfile_G << endl;
	}
	
	outfile_G.close();
	
// 	// Check Kramers-Kronig consistency of retarded Green's function:
// 	boost::numeric::ublas::vector<cmplx> iGR_vec(N_omega);
// 	boost::numeric::ublas::vector<cmplx> iGR_KK_vec(N_omega);
// 	for (int i = 0; i < N_omega; ++i) {
// 		iGR_vec(i) = iGR(3*N/4, omega_vec(i), N, h, iG_gtr);
// 	}
// 	
// 	for (int i = 0; i < N_omega; ++i) {
// 		double iGR_KK_real_i = 0;
// 		double iGR_KK_imag_i = 0;
// 		for (int j = 0; j < N_omega; ++j) {
// 			iGR_KK_real_i += (i != j) ? iGR_vec(j).imag() / (M_PI * (omega_vec(j) - omega_vec(i))) : 0;
// 			iGR_KK_imag_i += (i != j) ? -iGR_vec(j).real() / (M_PI * (omega_vec(j) - omega_vec(i))) : 0;
// 		}
// 		iGR_KK_vec(i) = {iGR_KK_real_i, iGR_KK_imag_i};
// 		iGR_KK_vec(i) *= (omega_vec(1) - omega_vec(0));
// 	}
// 	
// 	for (int i = 0; i < N_omega; ++i) {
// 		cout << omega_vec(i) << '\t' << iGR_KK_vec(i).real() << '\t' << iGR_KK_vec(i).imag() << endl;
// 	}
	
	// Create simple output file in order to assure that copying of data works (in load leveler script)
	ofstream outfileTest("job_finished", ios_base::out);
	outfileTest << endl;
	outfileTest.close();
}
