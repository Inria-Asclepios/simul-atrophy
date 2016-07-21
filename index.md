---
layout: index
---
### Simul@trophy
Simul@atrophy is a C++ based simulator tool that can simulate specified volume changes in different parts of the image.
The tool was developed to generate synthetic time-series structural MRIs with specified ground truth atrophy.
At its core, it implements a biophysical model of brain deformation with atrophy in Alzheimer's disease detailed in [1].
It can also simulate realistic variation of intensity.
The whole pipeline and the simulation of variation of intensity are detailed in [2].

The core components of the software are implemented in C++.
It also contains a number of helper scripts (written in python) to preprocess input images and generate desired map of atrophy (or volume changes).

The source code is available at https://github.com/Inria-Asclepios/simul-atrophy. 
Please cite the following papers if you use this software:

[1] Khanal, B., Lorenzi, M., Ayache, N., & Pennec, X. (2016). A biophysical model of brain deformation to simulate and analyze longitudinal MRIs of patients with Alzheimer’s disease . NeuroImage , 134, 35–52. http://doi.org/10.1016/j.neuroimage.2016.03.061

[2] Khanal, B., Ayache, N., & Pennec, X. (2016). Simulating Realistic Synthetic Longitudinal Brain MRIs with known Volume Changes. (under Review).

## Requirements
1. Cmake > 2.8
2. [PETSc](https://www.mcs.anl.gov/petsc/index.html) version 3.6.3.1
**Use `--with-clanguage=cxx` option** when configuring PETSc.
3. [ITK](https://itk.org/) version 4.9

## Installation
First make sure that:

* You have **right versions** of PETSc and ITK.
* That you configured PETSc with `--with-clanguage=cxx` option.
* Set `PETSC_DIR`, `PETSC_ARCH` and `ITK_DIR`.
For example, I have my `.bashrc` with:
```
ITK_DIR=/home/bkhanal/Documents/softwares/ITK-build
export ITK_DIR
PETSC_DIR=/home/bkhanal/Documents/softwares/petsc
export PETSC_DIR
PETSC_ARCH=arch-linux2-cxx-debug
export PETSC_ARCH
```

Now, build `simul@trophy` with the following steps. *I assume that you are inside a `curr_dir`, replace `curr_dir` with a directory of your choice.*

```bash
curr_dir$ git clone https://github.com/Inria-Asclepios/simul-atrophy.git
curr_dir$ cd simul-atrophy
simul-atrophy$ mkdir build && cd build
```
Use `cmake` to build the files.
```bash
build$ cmake ..
build$ make
```

If you have ITK installed with `c++11`, then use the following options with cmake:

```bash
build$ cmake -DUSE_CXX11:bool=true ..
build$ make
```

This should install the software to your machine.

## Quickstart
The main simulator executable is `build/src/simul_atrophy`.
Run with `-h` option to see in detail how to use it.

Otherwise, there is also a python script at 'simul-atrophy/scripts/simul_atrophy.py' which you can use.
Once again, use `-h` option to see the details.


## Example:
Coming soon here.
