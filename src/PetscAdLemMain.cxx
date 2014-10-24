#include"AdLem3D.hxx"
#include "GlobalConstants.hxx"

#include <iostream>
#include <sstream>
#include <petscsys.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkWarpImageFilter.h>
#include "InverseDisplacementImageFilter.h"
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>

#include <itkComposeDisplacementFieldsImageFilter.h>

#include <itkAbsoluteValueDifferenceImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkMultiplyImageFilter.h>

static char help[] = "Solves AdLem model.\n\n";

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
    std::string     atrophyFileName, maskFileName, lambdaFileName;
    std::string     baselineImageFileName;  //used only when debug priority is highest.
    std::string     resultsPath;            // Directory where all the results will be stored.
    std::string     resultsFilenamesPrefix;           // Prefix for all the filenames of the results to be stored in the resultsPath.
    bool            useTensorLambda;
    int             numOfTimeSteps;
    bool            isMaskChanged;

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

        ierr = PetscOptionsGetString(NULL,"-maskFile",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(optionFlag) maskFileName = optionString;

        ierr = PetscOptionsGetString(NULL,"-imageFile",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(optionFlag) baselineImageFileName = optionString;

        ierr = PetscOptionsGetString(NULL,"-useTensorLambda",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(!optionFlag) {
            std::cerr<<"Must provide true or false for the option -useTensorLambda\n";
            return(EXIT_FAILURE);
        } else {
            if(strcmp(optionString,"true")==0)  useTensorLambda = true;
            else useTensorLambda = false;
        }

        ierr = PetscOptionsGetInt(NULL,"-numOfTimeSteps",&numOfTimeSteps,&optionFlag);CHKERRQ(ierr);
        if(!optionFlag) {
            std::cout<<"Using default number of steps: 1 since -numOfTimeSteps option was not used.\n";
            numOfTimeSteps = 1;
        } else {
            std::cout<<"Computing for "<<numOfTimeSteps<<" time steps"<<std::endl;
        }

        if (useTensorLambda) {
            ierr = PetscOptionsGetString(NULL,"-lambdaFile",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
            if(useTensorLambda) {
                if(!optionFlag) {
                    std::cerr<<"Must provide valid tensor image filename since -useTensorLambda is set true.\n";
                    return(EXIT_FAILURE);
                } else {
                    lambdaFileName = optionString;
                }
            }
        }

        ierr = PetscOptionsGetString(NULL,"-resPath",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(!optionFlag) {
            std::cout<<"Must provide a valid path with -resPath option: e.g. -resPath ~/results"<<std::endl;
            return(EXIT_FAILURE);
        } else {
            resultsPath = optionString;
        }

        ierr = PetscOptionsGetString(NULL,"-resultsFilenamesPrefix",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(!optionFlag) {
            std::cerr<<"MUST provide a valid prefix for output results filenames: e.g. -resultsFilenamesPrefix step1"<<std::endl;
            return(EXIT_FAILURE);
        } else {
            resultsFilenamesPrefix = optionString;
        }

        /*std::cout<<"atrophy file: "<<atrophyFileName<<std::endl;
        std::cout<<"mask file: "<<maskFileName<<std::endl;
        std::cout<<"results path: "<<resultsPath<<std::endl;*/

        //------------------------*** Set up the model parameters ***-----------------------//
        AdLem3D AdLemModel; //xn,yn,zn,1,1,1,1);
        AdLemModel.setWallVelocities(wallVelocities);
        if(useTensorLambda) {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n using tensor image to get lambda\n");
            AdLemModel.setLameParameters(true,true,1,1,1,1,lambdaFileName);
        }
        else {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n using scalar lambda\n");
            AdLemModel.setLameParameters(true,false);
        }

        //-------******** Read baseline image ****** ---------------//
        //TODO: Do this only when debug priority is set highest.
        AdLem3D::ScalarImageType::Pointer baselineImage = AdLem3D::ScalarImageType::New();
        {
            AdLem3D::ScalarImageReaderType::Pointer   imageReader = AdLem3D::ScalarImageReaderType::New();
            imageReader->SetFileName(baselineImageFileName);
            imageReader->Update();
            baselineImage = imageReader->GetOutput();
        }

        double k = 1;  //compressibility coefficient k for CSF region.
        //---------*********** Read and set baseline brainMask ****************---------//
        AdLem3D::ScalarImageType::Pointer baselineBrainMask = AdLem3D::ScalarImageType::New();
        {
            AdLem3D::ScalarImageReaderType::Pointer   imageReader = AdLem3D::ScalarImageReaderType::New();
            imageReader->SetFileName(maskFileName);
            imageReader->Update();
            baselineBrainMask = imageReader->GetOutput();
        }
        AdLemModel.setBrainMask(baselineBrainMask,maskLabels::CSF,k,true,maskLabels::NBR); //1 coeff for p dof.
        //        AdLemModel.setBrainMask(baselineBrainMask,constants::CSF_LABEL,1,false,constants::NBR_LABEL); //1 coeff for p dof.

        //        AdLemModel.setBrainMask(baselineBrainMask,constants::CSF_LABEL,0); //don't relase IC

        //-------------------*** Set the computational region***-------------------//
        //        unsigned int origin[3] = {0,0,0};
        //        unsigned int size[3] = {10,10,10};
        //AdLemModel.setDomainRegion(origin,size);
        AdLemModel.setDomainRegionFullImage();

        //-------------------------*** Read and set baseline atrophy map ***--------------------------//
        AdLem3D::ScalarImageType::Pointer baselineAtrophy = AdLem3D::ScalarImageType::New();
        {
            AdLem3D::ScalarImageReaderType::Pointer   imageReader = AdLem3D::ScalarImageReaderType::New();
            imageReader->SetFileName(atrophyFileName);
            imageReader->Update();
            baselineAtrophy = imageReader->GetOutput();
        }
        AdLemModel.setAtrophy(baselineAtrophy);

        //Define itk types required for the warping of the mask and atrophy map:
        typedef InverseDisplacementImageFilter<AdLem3D::VectorImageType> FPInverseType;
        typedef itk::WarpImageFilter<AdLem3D::ScalarImageType,AdLem3D::ScalarImageType,AdLem3D::VectorImageType> WarpFilterType;
        typedef itk::NearestNeighborInterpolateImageFunction<AdLem3D::ScalarImageType> InterpolatorFilterNnType;
        typedef itk::AbsoluteValueDifferenceImageFilter<AdLem3D::ScalarImageType,AdLem3D::ScalarImageType,AdLem3D::ScalarImageType> DiffImageFilterType;
        typedef itk::StatisticsImageFilter<AdLem3D::ScalarImageType> StatisticsImageFilterType;
        typedef itk::ComposeDisplacementFieldsImageFilter<AdLem3D::VectorImageType, AdLem3D::VectorImageType> VectorComposerType;

        AdLem3D::VectorImageType::Pointer composedDisplacementField; //declared outside loop because we need this for two different iteration steps.
        isMaskChanged = true;
        for (int t=1; t<=numOfTimeSteps; ++t) {
            //-------------- Get the string for the current time step and add it to the prefix of all the files to be saved -----//
            std::stringstream timeStep;
            timeStep << t;
            std::string stepString("T"+timeStep.str());

            //------------- Modify atrophy map to adapt to the provided mask. ----------------//
            AdLemModel.modifyAtrophy(maskLabels::CSF,0);  //CSF region, set zero atrophy
            AdLemModel.modifyAtrophy(maskLabels::NBR,0);  //non-brain region, set zero atrophy
            AdLemModel.writeAtrophyToFile(resultsPath + resultsFilenamesPrefix + stepString + "AtrophyModified.nii.gz");

            //---------------------------*** Solve the system of equations ****-------------------//
            AdLemModel.solveModel(isMaskChanged);
            //--- Must check here to recompute the operator ---//
            //TODO: the above comment.
            //This requires:
            //1. Modifying other classes to support this recomputation of the operator.

            //----------------------------**** Write the solutions and residuals ***----------------//
            //TODO: Add debug flag and write these solutions only for the case of debug.
            //or set level of detail and write different files based on this level.
            //E.g. Files to be saved:
            //1 => save final composed disp. field, final warped image.
            //2 => all of 1 + intermediate warped images.
            //3 => all of 2 + intermediate velocity fields.
            //4 => all of 3 + intermediate residual files.
            // can add more as required!
            AdLemModel.writeSolution(resultsPath + resultsFilenamesPrefix + stepString);
            //        AdLemModel.writeSolution(resultsPath + resultsFilenamesPrefix + stepString,true,true);
            AdLemModel.writeResidual(resultsPath + resultsFilenamesPrefix + stepString);

            //--------------- Compose velocity fields before inversion and the consequent warping ---------------//
            if(t == 1) {
                composedDisplacementField = AdLemModel.getVelocityImage();
            }
            else {
                VectorComposerType::Pointer vectorComposer = VectorComposerType::New();
                vectorComposer->SetDisplacementField(AdLemModel.getVelocityImage());
                vectorComposer->SetWarpingField(composedDisplacementField);
                vectorComposer->Update();
                composedDisplacementField = vectorComposer->GetOutput();
            }

            if(numOfTimeSteps > 1) {

                //------------*** Invert the current displacement field to create warping field *** -------------//
                FPInverseType::Pointer inverter1 = FPInverseType::New();
                inverter1->SetInput(composedDisplacementField);
                inverter1->SetErrorTolerance(1e-1);
                inverter1->SetMaximumNumberOfIterations(50);
                inverter1->Update();
                std::cout<<"tolerance not reached for "<<inverter1->GetNumberOfErrorToleranceFailures()<<" pixels"<<std::endl;
                AdLem3D::VectorImageType::Pointer warperField = inverter1->GetOutput();

                //---------------- *** Warp baseline brain mask with an itk warpFilter, nearest neighbor *** -------------------//
                //------------- *** Using the inverted composed field. ** --------------//
                WarpFilterType::Pointer warper = WarpFilterType::New();
                warper->SetDisplacementField(warperField);
                warper->SetInput(baselineBrainMask);
                warper->SetOutputSpacing(baselineBrainMask->GetSpacing());
                warper->SetOutputOrigin(baselineBrainMask->GetOrigin());
                warper->SetOutputDirection(baselineBrainMask->GetDirection());

                InterpolatorFilterNnType::Pointer nnInterpolatorFilter = InterpolatorFilterNnType::New();
                warper->SetInterpolator(nnInterpolatorFilter);
                warper->Update();

                //-------------------- ** Compare warped mask with the previous mask ** -----------------------//
                DiffImageFilterType::Pointer diffImageFilter = DiffImageFilterType::New();
                diffImageFilter->SetInput1(AdLemModel.getBrainMaskImage());
                diffImageFilter->SetInput2(warper->GetOutput());
                diffImageFilter->Update();

                StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
                statisticsImageFilter->SetInput(diffImageFilter->GetOutput());
                statisticsImageFilter->Update();
                if(statisticsImageFilter->GetSum() < 1)
                    isMaskChanged = false;
                else {
                    isMaskChanged = true;
                    AdLemModel.writeBrainMaskToFile(resultsPath + resultsFilenamesPrefix + stepString+"Mask.nii.gz");
                    AdLemModel.setBrainMask(warper->GetOutput(),maskLabels::CSF,k,true,maskLabels::NBR);
                }

                //------------------- ******* Warp baseline atrophy with an itk WarpFilter, linear interpolation *******-------//
                //------------- *** using the inverted composed field *** ---------------//
                WarpFilterType::Pointer warper1 = WarpFilterType::New();
                warper1->SetDisplacementField(warperField);
                warper1->SetInput(baselineAtrophy);
                warper1->SetOutputSpacing(baselineAtrophy->GetSpacing());
                warper1->SetOutputOrigin(baselineAtrophy->GetOrigin());
                warper1->SetOutputDirection(baselineAtrophy->GetDirection());
                warper1->Update();
                AdLemModel.setAtrophy(warper1->GetOutput());
                //WHAT DOES IT MEAN for an atrophy present in some voxels to be
                //substituted by 0 as a result of these voxels being converted to CSF voxels due to NN interpolation of the mask ??
                //Or, WHAT ABOUT the dilution of atrophy value due to linear interpolation because of zero atrophy contribution from
                //the neighboring CSF or NBR regions ?

                //-------------- **** Warp the baseline image with the inverted composed field with linear interpolation ****------//
                //TODO: Do this and save only when the debug priority level is set to maximum.
                WarpFilterType::Pointer warper2 = WarpFilterType::New();
                warper2->SetDisplacementField(warperField);
                warper2->SetInput(baselineImage);
                warper2->SetOutputSpacing(baselineImage->GetSpacing());
                warper2->SetOutputOrigin(baselineImage->GetOrigin());
                warper2->SetOutputDirection(baselineImage->GetDirection());
                warper2->Update();
                AdLem3D::ScalarImageWriterType::Pointer imageWriter = AdLem3D::ScalarImageWriterType::New();
                //Write with the stepString lying at the end of the filename so that other tool can combine all the images
                //into 4D by using the number 1, 2, 3 ... at the end.
                imageWriter->SetFileName(resultsPath + resultsFilenamesPrefix + "WarpedImage" + stepString+ ".nii.gz");
                imageWriter->SetInput(warper2->GetOutput());
                imageWriter->Update();
            }
        }
        AdLem3D::VectorImageWriterType::Pointer   displacementWriter = AdLem3D::VectorImageWriterType::New();
        displacementWriter->SetFileName(resultsPath+resultsFilenamesPrefix+"ComposedField.nii.gz");
        displacementWriter->SetInput(composedDisplacementField);
        displacementWriter->Update();

    }
    PetscErrorCode ierr;
    ierr = PetscFinalize();CHKERRQ(ierr);
    return 0;
}

