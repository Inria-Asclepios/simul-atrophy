#include"AdLem3D.hxx"
#include "GlobalConstants.hxx"

#include <iostream>
#include <sstream>
#include <petscsys.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkWarpImageFilter.h>
#include <itkMultiplyImageFilter.h>

static char help[] = "Solves AdLem model.\n\n";

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
    std::string atrophyFileName, maskFileName, lambdaFileName;
    std::string resultsPath;
    int stepIndex;

    std::vector<double> wallVelocities(18);
    /*0,1,2,     //south wall
       3,4,5,     //west wall
       6,7,8,     //north wall
       9,10,11,     //east wall
       12,13,14,     //front wall
       15,16,17    //back wall*/
    //    unsigned int wallPos = 6;
    //        wallVelocities.at(wallPos) = 1;

    PetscInitialize(&argc,&argv,(char*)0,help);
    {
        //---------------*** Get the atrophy, mask file and result directory ***--------------//
        PetscErrorCode ierr;
        PetscBool optionFlag = PETSC_FALSE;
        char optionString[PETSC_MAX_PATH_LEN];
        ierr = PetscOptionsGetString(NULL,"-atrophyFile",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(optionFlag) atrophyFileName = optionString;
        ierr = PetscOptionsGetString(NULL,"-lambdaFile",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(optionFlag) lambdaFileName = optionString;
        ierr = PetscOptionsGetString(NULL,"-maskFile",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(optionFlag) maskFileName = optionString;
        ierr = PetscOptionsGetString(NULL,"-resPath",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(!optionFlag) {
            std::cout<<"Must provide a valid path with -resPath option: e.g. -resPath ~/results"<<std::endl;
            return(EXIT_FAILURE);
        } else {
            resultsPath = optionString;
        }
        ierr = PetscOptionsGetInt(NULL,"-stepIndex",&stepIndex,&optionFlag);CHKERRQ(ierr);
        if(!optionFlag) {
            std::cerr<<"MUST provide a valid step index option: e.g. -stepIndex 1"<<std::endl;
            return(EXIT_FAILURE);
        }

        /*std::cout<<"atrophy file: "<<atrophyFileName<<std::endl;
        std::cout<<"mask file: "<<maskFileName<<std::endl;
        std::cout<<"results path: "<<resultsPath<<std::endl;*/

        //------------------------*** Set up the model parameters ***-----------------------//
        AdLem3D AdLemModel; //xn,yn,zn,1,1,1,1);
        AdLemModel.setWallVelocities(wallVelocities);
//        AdLemModel.setLameParameters(true,false);
        AdLemModel.setLameParameters(true,true,1,1,1,1,lambdaFileName);

        AdLemModel.setBrainMask(maskFileName,constants::CSF_LABEL,1,true,constants::NBR_LABEL); //1 coeff for p dof.
//        AdLemModel.setBrainMask(maskFileName,constants::CSF_LABEL,1,false,constants::NBR_LABEL); //1 coeff for p dof.

        //        AdLemModel.setBrainMask(maskFileName,constants::CSF_LABEL,0); //don't relase IC

        //-------------------*** Set the computational region***-------------------//
        //        unsigned int origin[3] = {0,0,0};
        //        unsigned int size[3] = {10,10,10};
        //AdLemModel.setDomainRegion(origin,size);
        AdLemModel.setDomainRegionFullImage();

        //-------------------------*** Set up the atrophy map ***--------------------------//
        AdLemModel.setAtrophy(atrophyFileName);
        AdLemModel.modifyAtrophy(constants::CSF_LABEL,0);  //CSF region, set zero atrophy
        AdLemModel.modifyAtrophy(constants::NBR_LABEL,0);  //non-brain region, set zero atrophy

        std::stringstream fileIndex;
        fileIndex << stepIndex;
        AdLemModel.writeAtrophyToFile(resultsPath + "step" + fileIndex.str() + "AtrophyModified.nii.gz");
        AdLemModel.solveModel();
        AdLemModel.writeSolution(resultsPath + "step" + fileIndex.str());
        //        AdLemModel.writeSolution(resultsPath,true,true);
        AdLemModel.writeResidual(resultsPath + "step" + fileIndex.str());
    }
    PetscErrorCode ierr;
    ierr = PetscFinalize();CHKERRQ(ierr);
    return 0;
}

