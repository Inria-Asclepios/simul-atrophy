#include"ADLinearElastic.hxx"
#include<petscsys.h>
#include<petscksp.h>
#include<petscdm.h>

static char help[] = "Solves 2D inhomogeneous Laplacian using multigrid.\n\n";
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  KSP            ksp;
  DM             da;
//  const char     *bcTypes[2] = {"dirichlet","neumann"};
  PetscErrorCode ierr;
//  PetscInt       bc;
  Vec            b,x;
  ADLinearElastic model1(10,10,1.0);
//  ADLinearElastic *model_ptr = &model1;

  PetscInitialize(&argc,&argv,(char*)0,help);

  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = DMDACreate2d(PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE,DMDA_STENCIL_STAR,-3,-3,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0,&da);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(da,0,1,0,1,0,0);CHKERRQ(ierr);
  ierr = DMDASetFieldName(da,0,"Pressure");CHKERRQ(ierr);

//  ierr        = PetscOptionsBegin(PETSC_COMM_WORLD, "", "Options for the inhomogeneous Poisson equation", "DMqq");
//  user.rho    = 1.0;
//  ierr        = PetscOptionsReal("-rho", "The conductivity", "ex29.c", user.rho, &user.rho, NULL);CHKERRQ(ierr);
//  user.nu     = 0.1;
//  ierr        = PetscOptionsReal("-nu", "The width of the Gaussian source", "ex29.c", user.nu, &user.nu, NULL);CHKERRQ(ierr);
//  bc          = (PetscInt)DIRICHLET;
//  ierr        = PetscOptionsEList("-bc_type","Type of boundary condition","ex29.c",bcTypes,2,bcTypes[0],&bc,NULL);CHKERRQ(ierr);
//  user.bcType = (BCType)bc;
//  ierr        = PetscOptionsEnd();

  ierr = KSPSetComputeRHS(ksp,ComputeRHS,&model1);CHKERRQ(ierr);
  ierr = KSPSetComputeOperators(ksp,ComputeMatrix,&model1);CHKERRQ(ierr);
  ierr = KSPSetDM(ksp,da);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,NULL,NULL);CHKERRQ(ierr);
  ierr = KSPGetSolution(ksp,&x);CHKERRQ(ierr);
  ierr = KSPGetRhs(ksp,&b);CHKERRQ(ierr);

  ierr = DMDestroy(&da);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = PetscFinalize();

  return 0;
}
