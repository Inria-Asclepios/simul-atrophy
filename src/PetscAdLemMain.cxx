#include"AdLem3D.hxx"

#include<iostream>
#include<iomanip>
#include<fstream>

#include<petscsys.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

static char help[] = "Solves AdLem model.\n\n";

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
    std::string divFileName, maskFileName;
    std::string resultsPath;

    std::vector<double> wallVelocities(18);
    /*0,1,2,     //south wall
       3,4,5,     //west wall
       6,7,8,     //north wall
       9,10,11,     //east wall
       12,13,14,     //front wall
       15,16,17    //back wall*/
    unsigned int wallPos = 6;
        wallVelocities.at(wallPos) = 1;

    PetscInitialize(&argc,&argv,(char*)0,help);
    {
        //---------------*** Get the atrophy, mask file and result directory ***--------------//
        PetscErrorCode ierr;
        PetscBool optionFlag = PETSC_FALSE;
        char optionString[PETSC_MAX_PATH_LEN];
        ierr = PetscOptionsGetString(NULL,"-divFile",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(optionFlag) divFileName = optionString;
        ierr = PetscOptionsGetString(NULL,"-maskFile",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(optionFlag) maskFileName = optionString;
        ierr = PetscOptionsGetString(NULL,"-resPath",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(!optionFlag) {
            std::cout<<"Must provide a valid path with -resPath option: e.g. -resPath ~/results"<<std::endl;
            return(EXIT_FAILURE);
        } else {
            resultsPath = optionString;
        }

        //------------------------*** Set up the model parameters ***-----------------------//
        AdLem3D AdLemModel; //xn,yn,zn,1,1,1,1);
        AdLemModel.setWallVelocities(wallVelocities);
        AdLemModel.setLameParameters(1,1);
        //        AdLemModel.setLameParameters(1,1,false,10,20);
//        AdLemModel.setBrainMask(maskFileName,1,1); //1 is the CSF label, 1 coeff for p dof.
        AdLemModel.setBrainMask(maskFileName,1,0); //don't relase IC

        //-------------------*** Set the computational region***-------------------//
        unsigned int origin[3] = {0,0,0};
        unsigned int size[3] = {10,10,10};
        //origin[0] = 50;  origin[1] = 50;  origin[2] = 50;
        //size[0] = 182;   size[1] = 218;       size[2] = 182;
        //AdLemModel.setDomainRegion(origin,size);
        AdLemModel.setDomainRegion(origin,size,true);

        //-------------------------*** Set up the atrophy map ***--------------------------//
        AdLemModel.setAtrophy(divFileName);
        //AdLemModel.writeAtrophyToFile(resultsPath + "atrophyOrig.mha");
        //        AdLemModel.scaleAtorphy(-1);
        /*AdLemModel.modifyAtrophy(1,0);  //CSF region, set zero atrophy
        AdLemModel.modifyAtrophy(0,0);  //non-brain region, set zero atrophy
        AdLemModel.writeAtrophyToFile(resultsPath + "atrophyModified.mha");*/

        //-----------------------------*** Solve the system ***----------------------------//
        AdLemModel.solveModel();

        //-------------------------***  Write the results ***-------------------------------//
        AdLemModel.writeSolution(resultsPath);
        //        AdLemModel.writeSolution(resultsPath,true,true);
        AdLemModel.writeResidual(resultsPath);

    }
    PetscErrorCode ierr;
    ierr = PetscFinalize();CHKERRQ(ierr);
    return 0;
}

