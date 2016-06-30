# Simul@trophy
## Simulating Realistic Longitudinal Images with Atrophy

### Requirements
1. Cmake > 2.8
2. [PETSc](https://www.mcs.anl.gov/petsc/index.html) version 3.6.3.1
3. [ITK](https://itk.org/) version 4.9

First install right versions of PETSc and ITK.

## Quickstart
Clone this repository and create a `build` directory inside it.
Get into the `build` directory and run `cmake ..` followed by `make all`

The main simulator executable is `build/src/AdLemMain`.
Run with `-h` option to see the arguments that must be given.

If you have python, it might be easier to use `scripts/runAdLemModel.py` instead of the executable.
Once again, use `-h` option to see the details.


