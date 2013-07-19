#include"AdLem2D.hxx"
#include"PetscAdLemTaras2D.hxx"

#include<petscsys.h>


static char help[] = "Solves 2D inhomogeneous Laplacian using multigrid.\n\n";
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
    int xn = 16;
    int yn = 16;
    AdLem2D model(xn,yn,20,10,1,1);
    AdLem2D *model_ptr;
    model_ptr = &model;

    PetscInitialize(&argc,&argv,(char*)0,help);
    { /*This scope is necessary because the object PetscAdLemTaras2D will free all the
        memory allocated to petsc-objects within its methods only when its destructor
        is called. For this call to happen, the object should be gone out of the scope.
        If PetscFinalize() is called before the object go out of scope, it will result
        in error because the PetscFinalize() will see pets-objects allocated but not freed!*/

        PetscAdLemTaras2D solver1(model_ptr);
        std::string wFileName("/user/bkhanal/home/works/ADDeformationModeling/src/matlab/lin_sys");
        solver1.solveModel(true,wFileName);
//        solver1.writeLinearSystemMatlabFile("/user/bkhanal/home/works/ADDeformationModeling/src/matlab/lin_sys");
    }
    PetscErrorCode ierr;

    ierr = PetscFinalize();CHKERRXX(ierr);

    return 0;
}

