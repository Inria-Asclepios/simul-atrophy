# Simul@trophy
## Simulating Realistic Longitudinal Images with Atrophy

### Requirements
1. Cmake > 2.8
2. [PETSc](https://www.mcs.anl.gov/petsc/index.html) version 3.12.4
**Use `--with-clanguage=cxx` option** when configuring PETSc.
**Consider  --with-64-bit-indices=yes if using larger images.
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

### A basic example
One very basic example of how to use simul@trophy is shown below on a relative small image so that this can be run quickly on a normal laptop/desktop.
Note that if your image sizes are too large, usually greater than 120^3 you must use cluster computing for running the model due to memory limitations on normal desktops.
Here I have used a dummy object instead of brain, but you can create your own test set by downsampling brain images too.
Any file format readable by ITK will work; the example here uses `.mha` format.

The input files and the results of this basic example can be found at the directory `basicExample`.

There are three main input files:

1. Input atrophy map: `test1Atrophy1.mha`

1. Input segmentation file: `bMask.mha`

1. Input image file: `bTest1.mha`

The input parameters $\mu$ and $\lambda$ are set to 1 for both tissue and CSF regions.

The example is run for two steps.

```
$build/src/simul_atrophy -parameters 1,1,1,1 -boundary_condition dirichlet_at_skull --relax_ic_in_csf -atrophyFile basicExample/test1Atrophy1.mha -maskFile basicExample/bMask1.mha -imageFile basicExample/bTest1.mha --invert_field_to_warp -numOfTimeSteps 2 -resPath ./basicExample/ -resultsFilenamesPrefix basicTest1_
```

If everything is all right you should be getting the following output in your display terminal:
```
Model will be run for 2 time steps
dmda of size: (21,41,21)
grid with the spacings hx, hy, hz: (1.000000,1.000000,1.000000)
dirichlet_at_skull=>Skull velocity = 0.
Incompressibility constraint relaxed where brain mask has label 1 
Incompressibility constraint relaxation coefficient: 1/lambda.
solver uses discretization for constant viscosity case!
 muBrain=1.000000, muCsf=1.000000, lambdaBrain=1.000000, lambdaCsf=1.000000
First Lame Parameter (lambda) is a scalar.

 RHS will be taken as (mu + lambda)grad(a), i.e. with Lame parameters.

 computing the operator for linear solve with 8 point stencil for divergence

 Relax IC with div(u) + kp = 0 on relax cells

 Displacement field inversion: tolerance not reached in 0 voxels 


 Displacement field inversion: tolerance not reached in 0 voxels 
```

The output files created will be placed in the folder `basicExample`.
The repository already contains the results of same example with prefix `basicTest_`.
Now, after running this example you should get new files with the prefix `basicTest1_`.


Instead of using the `simul_atrophy` executable directly, you can also a run the model by using the python script at 'simul-atrophy/scripts/simul_atrophy.py'.
Please see `-h` option to see the details of the script.


