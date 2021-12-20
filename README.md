# KCAP BOSS DR12 correlation function wedges module

This repository provides a build system and Python/CosmoSIS interfaces to compute the predictions for the BOSS DR12 correlation function wedges as used in [Sánchez et al. 2017](https://arxiv.org/abs/1607.03147), as well as [Tröster et al. 2020](https://arxiv.org/abs/1909.11006), and the KiDS-1000 3x2pt analysis in [Heymans, Tröster et al. 2021](https://arxiv.org/abs/2007.15632). The code is essentially a wrapper around the original CosmoMC implementation. 
Please cite these works if you use this module in your work.

## Installation
```
mkdir build
cd build
cmake ..
make boss_module
```
On macOS, the compilers need to be set to some recent GCC, e.g., `CC=gcc-8 CXX=g++-8`.
