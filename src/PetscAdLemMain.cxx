//#include"AdLem2D.hxx"
//#include"PetscAdLemTaras2D.hxx"

#include"AdLem3D.hxx"
#include"PetscAdLemTaras3D.hxx"
#include<petscsys.h>

#include<iostream>
#include<fstream>

static char help[] = "Solves AdLem model.\n\n";

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
    int xn = 16;
    int yn = 16;
    int zn = 16;
//    AdLem2D model(xn,yn,20,10,1,1);
//    AdLem2D *model_ptr;
    AdLem3D model(xn,yn,zn,20,10,1,1);
    AdLem3D *model_ptr;

    model_ptr = &model;

    PetscInitialize(&argc,&argv,(char*)0,help);
    { /*This scope is necessary because the object PetscAdLemTaras2D will free all the
        memory allocated to petsc-objects within its methods only when its destructor
        is called. For this call to happen, the object should be gone out of the scope.
        If PetscFinalize() is called before the object go out of scope, it will result
        in error because the PetscFinalize() will see pets-objects allocated but not freed!*/

//        PetscAdLemTaras2D solver1(model_ptr);
        PetscAdLemTaras3D solver1(model_ptr);
        std::string wFileName("/user/bkhanal/home/works/AdLemModel/results/lin_sys");
        solver1.solveModel(true,wFileName);
    }
    PetscErrorCode ierr;

    ierr = PetscFinalize();CHKERRQ(ierr);
    std::ofstream size_file;
    size_file.open("/user/bkhanal/home/works/AdLemModel/results/size_lin_sys");
    size_file<<xn+1<<" "<<yn+1<<" "<<zn+1;
    size_file.close();

    return 0;
}

