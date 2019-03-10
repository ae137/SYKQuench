# Makefile for main_quench

# Macro for compiler name
CC = g++

# Macro for compiler options
CFLAGS = -Wall -g -O3 -DHAVE_INLINE -march=native --std=c++11 -fopenmp

# Macro for linker options
LFLAGS = -Wall -g -O3 -march=native -lm --std=c++11 -fopenmp

all: main_quench

main_quench: main_quench.o
	$(CC) $(LFLAGS) -o main_quench main_quench.o

main_quench.o: main_quench.cpp
	$(CC) $(CFLAGS) -c main_quench.cpp

clean:
	rm *.o
	rm main_quench
