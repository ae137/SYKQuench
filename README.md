# Non-equilibrium dynamics of Majorana fermions with all-to-all interactions

This repository contains code for the numerical simulation of the non-equilibrium dynamics of Majorana fermions with all-to-all interactions. The non-equilibrium dynamics is caused by a quantum quench and described using the Kadanoff-Baym equations.

The code in this repository was used for obtaining the results that were published in "Quantum quench of the Sachdev-Ye-Kitaev Model" ([https://arxiv.org/abs/1706.07803](https://arxiv.org/abs/1706.07803)). The paper also contains more information.

### Note:
The code in this repository is research code and not under active development and without any support. It is quite old already and nowadays I would presumably implement a few things differently. Nevertheless, I hope publishing this code is helpful to some people for implementing their own numerical solver of Kadanoff-Baym equations for Majorana fermions.


### Dependencies:
The source code is header-only and includes several headers from

    * C++ standard headers
    * Boost library (for example libboost1.67)

### License
GPLv3, see LICENSE for more information
