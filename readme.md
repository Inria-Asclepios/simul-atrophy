# Simul@trophy
## Simulating Realistic Longitudinal Images with Atrophy

### Requirements
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


