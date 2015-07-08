#include "AdLem3D.hxx"
#include "GlobalConstants.hxx"

#include <iostream>
#include <sstream>
#include <petscsys.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkWarpImageFilter.h>
//#include "InverseDisplacementImageFilter.h"
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include "itkLabelImageGenericInterpolateImageFunction.h"

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
    // Options for results storage: Velocity and divergence will always be written to file. For others options must be passed.
    bool            writePressure, writeForce, writeResidual;

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
            numOfTimeSteps = 1;
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n Using default number of steps: 1 since -numOfTimeSteps option was not used.\n");
        } else {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n Computing for %d time steps", numOfTimeSteps);
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

        ierr = PetscOptionsGetString(NULL,"-writePressure",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(!optionFlag) {
            writePressure = false;
        } else {
            if(strcmp(optionString,"true")==0)  writePressure = true;
            else writePressure = false;
        }

        ierr = PetscOptionsGetString(NULL,"-writeForce",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(!optionFlag) {
            writeForce = false;
        } else {
            if(strcmp(optionString,"true")==0)  writeForce = true;
            else writeForce = false;
        }

        ierr = PetscOptionsGetString(NULL,"-writeResidual",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(!optionFlag) {
            writeResidual = false;
        } else {
            if(strcmp(optionString,"true")==0)  writeResidual = true;
            else writeResidual = false;
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
        //typedef InverseDisplacementImageFilter<AdLem3D::VectorImageType> FPInverseType;
        typedef itk::WarpImageFilter<AdLem3D::ScalarImageType,AdLem3D::ScalarImageType,AdLem3D::VectorImageType> WarpFilterType;
        typedef itk::NearestNeighborInterpolateImageFunction<AdLem3D::ScalarImageType> InterpolatorFilterNnType;
	typedef itk::LabelImageGenericInterpolateImageFunction<AdLem3D::ScalarImageType, itk::LinearInterpolateImageFunction> InterpolatorGllType; //General Label interpolator with linear interpolation.
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
	    // ---------- do the modification after the first step. That means I expect the atrophy map to be valid
	    // ---------- when input by the user. i.e. only GM/WM has atrophy and 0 on CSF and NBR regions. 

            //---------------------------*** Solve the system of equations ****-------------------//
            AdLemModel.solveModel(isMaskChanged);
            //----------------------------**** Write the solutions and residuals ***----------------//
            AdLemModel.writeVelocityImage(resultsPath+resultsFilenamesPrefix+stepString+"vel.nii.gz");
            AdLemModel.writeDivergenceImage(resultsPath+resultsFilenamesPrefix+stepString+"div.nii.gz");
            if (writeForce) AdLemModel.writeForceImage(resultsPath+resultsFilenamesPrefix+stepString+"force.nii.gz");
            if (writePressure) AdLemModel.writePressureImage(resultsPath+resultsFilenamesPrefix+stepString+"press.nii.gz");
            if (writeResidual) AdLemModel.writeResidual(resultsPath + resultsFilenamesPrefix + stepString);

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
//                FPInverseType::Pointer inverter1 = FPInverseType::New();
//                inverter1->SetInput(composedDisplacementField);
//                inverter1->SetErrorTolerance(1e-1);
//                inverter1->SetMaximumNumberOfIterations(50);
//                inverter1->Update();
//                std::cout<<"tolerance not reached for "<<inverter1->GetNumberOfErrorToleranceFailures()<<" pixels"<<std::endl;
//                AdLem3D::VectorImageType::Pointer warperField = inverter1->GetOutput();
                //----------------*** Let's not invert the field, rather assume the atrophy is provided to be negative
                //---------------- so that we can use the field obtained from the model itself as being already inverted.//
                AdLem3D::VectorImageType::Pointer warperField = composedDisplacementField;

                //---------------- *** Warp baseline brain mask with an itk warpFilter, nearest neighbor *** -------------------//
                //------------- *** Using the inverted composed field. ** --------------//
                WarpFilterType::Pointer warper = WarpFilterType::New();
                warper->SetDisplacementField(warperField);
                warper->SetInput(baselineBrainMask);
                warper->SetOutputSpacing(baselineBrainMask->GetSpacing());
                warper->SetOutputOrigin(baselineBrainMask->GetOrigin());
                warper->SetOutputDirection(baselineBrainMask->GetDirection());

                InterpolatorFilterNnType::Pointer nnInterpolatorFilter = InterpolatorFilterNnType::New();
		InterpolatorGllType::Pointer gllInterpolator = InterpolatorGllType::New();
                warper->SetInterpolator(nnInterpolatorFilter);
		//warper->SetInterpolator(gllInterpolator);
                warper->Update();

                //-------------------- ** Compare warped mask with the previous mask ** -----------------------//
                DiffImageFilterType::Pointer diffImageFilter = DiffImageFilterType::New();
                diffImageFilter->SetInput1(AdLemModel.getBrainMaskImage());
                diffImageFilter->SetInput2(warper->GetOutput());
                diffImageFilter->Update();

                StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
                statisticsImageFilter->SetInput(diffImageFilter->GetOutput());
                statisticsImageFilter->Update();
                if(statisticsImageFilter->GetSum() < 1) {
                    isMaskChanged = false;
                } else {
                    isMaskChanged = true;
                    AdLemModel.setBrainMask(warper->GetOutput(),maskLabels::CSF,k,true,maskLabels::NBR);
                    AdLemModel.writeBrainMaskToFile(resultsPath + resultsFilenamesPrefix + stepString+"Mask.nii.gz");
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
                //Atrophy present at the newly created CSF regions are redistributed to the nearest GM/WM tissues voxels.
	        AdLemModel.modifyAtrophy(maskLabels::CSF,0,true); //redistributeatrophy; then set zero atrophy in CSF region.
		//AdLemModel.modifyAtrophy(maskLabels::CSF,0,false); //set zero atrophy in CSF region without redistribution.
		AdLemModel.modifyAtrophy(maskLabels::NBR,0);  //set zero atrophy at non-brain region.
		AdLemModel.writeAtrophyToFile(resultsPath + resultsFilenamesPrefix + stepString + "AtrophyModified.nii.gz");

                //-------------- **** Warp the baseline image with the inverted composed field with linear interpolation ****------//
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

                //-------------- **** Warp the baseline image with b-spline interpolation too. ****------//
		typedef itk::BSplineInterpolateImageFunction<AdLem3D::ScalarImageType> InterpolatorFilterType;
		InterpolatorFilterType::Pointer interpolatorFilter = InterpolatorFilterType::New();
		interpolatorFilter->SetSplineOrder(3);
                WarpFilterType::Pointer warper3 = WarpFilterType::New();
                warper3->SetDisplacementField(warperField);
		warper3->SetInterpolator(interpolatorFilter);
                warper3->SetInput(baselineImage);
                warper3->SetOutputSpacing(baselineImage->GetSpacing());
                warper3->SetOutputOrigin(baselineImage->GetOrigin());
                warper3->SetOutputDirection(baselineImage->GetDirection());
                warper3->Update();
                AdLem3D::ScalarImageWriterType::Pointer imageWriter1 = AdLem3D::ScalarImageWriterType::New();
                //Write with the stepString lying at the end of the filename so that other tool can combine all the images
                //into 4D by using the number 1, 2, 3 ... at the end.
                imageWriter1->SetFileName(resultsPath + resultsFilenamesPrefix + "WarpedImageBspline" + stepString+ ".nii.gz");
                imageWriter1->SetInput(warper3->GetOutput());
                imageWriter1->Update();

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

