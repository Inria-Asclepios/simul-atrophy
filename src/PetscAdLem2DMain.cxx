#include"AdLem2D.hxx"
#include"PetscAdLemTaras2D.hxx"

#include<petscsys.h>


static char help[] = "Solves 2D inhomogeneous Laplacian using multigrid.\n\n";
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
    PetscInitialize(&argc,&argv,(char*)0,help);

    {AdLem2D model(5,5);
    AdLem2D *model_ptr;
    model_ptr = &model;

    PetscAdLemTaras2D solver1(model_ptr);
    solver1.solveModel();}
    PetscErrorCode ierr;

    ierr = PetscFinalize();CHKERRXX(ierr);

    return 0;
}

