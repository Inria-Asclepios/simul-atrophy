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
    std::string divFileName;
    std::string resultsPath;

    PetscInitialize(&argc,&argv,(char*)0,help);
    {
        PetscErrorCode ierr;
        PetscBool optionFlag = PETSC_FALSE;
        char optionString[PETSC_MAX_PATH_LEN];
        ierr = PetscOptionsGetString(NULL,"-divFile",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(optionFlag) divFileName = optionString;
        ierr = PetscOptionsGetString(NULL,"-resPath",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(!optionFlag) {
            std::cout<<"Must provide a valid path with -resPath option: e.g. -resPath ~/results"<<std::endl;
            return(EXIT_FAILURE);
        } else {
            resultsPath = optionString;
        }
        AdLem3D AdLemModel; //xn,yn,zn,1,1,1,1);
        unsigned int origin[3] = {0,0,0};
        unsigned int size[3] = {10,10,10};
        AdLemModel.createAtrophy(size);
        AdLemModel.setDomainRegion(origin,size,true);

        //        AdLemModel.writeAtrophyToFile(resultsPath + "atrophyOrig.mha");


        //AdLemModel.setAtrophy(divFileName);
        //AdLemModel.setAtrophy(divFileName);
        //origin[0] = 50;  origin[1] = 50;  origin[2] = 50;
        //size[0] = 182;   size[1] = 218;       size[2] = 182;
        //AdLemModel.setDomainRegion(origin,size);
        //AdLemModel.setDomainRegion(origin,size,true);
        //Must call modifyAtrophy after setDomainRegion, if you want
        //to use valid atrophy for this new region!
        //        AdLemModel.scaleAtorphy(-1);

        //        AdLemModel.modifyAtrophy();
        //        AdLemModel.writeAtrophyToFile(resultsPath + "atrophyModified.mha");

        AdLemModel.setLameParameters(1,1,1,1);

        //mu, lambda value has effect on no. of iterations!!
        AdLemModel.solveModel(); //it's not just the ratios!! Figure out how it affects the velocity!!
        AdLemModel.writeSolution(resultsPath);
        //        AdLemModel.writeSolution(resultsPath,true);
        AdLemModel.writeResidual(resultsPath);
    }
    PetscErrorCode ierr;
    ierr = PetscFinalize();CHKERRQ(ierr);
    return 0;
}

