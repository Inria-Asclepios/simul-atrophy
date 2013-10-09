//#include"AdLem2D.hxx"
//#include"PetscAdLemTaras2D.hxx"

#include"AdLem3D.hxx"
//#include"PetscAdLemTaras3D.hxx"
#include<petscsys.h>

#include<iostream>
#include<iomanip>
#include<fstream>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

static char help[] = "Solves AdLem model.\n\n";

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
    std::string dataPath("/user/bkhanal/home/works/AdLemModel/data/Template/");
    std::string resultsPath("/user/bkhanal/home/works/AdLemModel/results/");
//    std::string imVelocityFile("/user/bkhanal/home/works/AdLemModel/results/vel.mha");
//    std::string imPressureFile("/user/bkhanal/home/works/AdLemModel/results/press.mha");

    PetscInitialize(&argc,&argv,(char*)0,help);
    { /*This scope is necessary because the object PetscAdLemTaras2D will free all the
        memory allocated to petsc-objects within its methods only when its destructor
        is called. For this call to happen, the object should be gone out of the scope.
        If PetscFinalize() is called before the object go out of scope, it will result
        in error because the PetscFinalize() will see pets-objects allocated but not freed!*/

        AdLem3D AdLemModel; //xn,yn,zn,1,1,1,1);
        unsigned int origin[3] = {95,110,95};
        unsigned int size[3] = {40,40,40};
//        AdLemModel.createAtrophy(size);
//        AdLemModel.writeAtrophyToFile(resultsPath + "atrophyOrig.mha");
        AdLemModel.setAtrophy(dataPath+"divergence.nii");

        /*for(unsigned int i=0; i<3;++i) {
            origin[i] = 0;
            size[i] = 18;
        }*/
        AdLemModel.setDomainRegion(origin,size);
        //Must call modifyAtrophy after setDomainRegion, if you want
        //to use valid atrophy for this new region!
        AdLemModel.scaleAtorphy(-1);
        AdLemModel.modifyAtrophy();
//        AdLemModel.writeAtrophyToFile(resultsPath + "atrophyModified.mha");
        AdLemModel.writeAtrophyToFile(resultsPath + "divergenceModified.nii");

        AdLemModel.setLameParameters(2,2,1,1);

        //Test if atrophy for sub region is correct:
        /*std::cout<<std::endl;
        for(unsigned int k=0; k<size[2]; ++k) {
            for(unsigned int j=0; j<size[1]; ++j) {
                for(unsigned int i=0; i<size[0]; ++i) {
                    std::cout<<std::setw(16)<<AdLemModel.dataAt("atrophy",i,j,k);
                }
                std::cout<<std::endl;
            }
            std::cout<<std::endl<<std::endl;
        }*/
        //mu, lambda value has effect on no. of iterations!!
        AdLemModel.solveModel(); //it's not just the ratios!! Figure out how it affects the velocity!!
        AdLemModel.writeSolution(resultsPath);

    }
    PetscErrorCode ierr;
    ierr = PetscFinalize();CHKERRQ(ierr);
    return 0;
}

