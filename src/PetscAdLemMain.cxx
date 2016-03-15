#include "AdLem3D.h"
#include "GlobalConstants.h"

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
#include <itkBSplineInterpolateImageFunction.h>
//#include "itkLabelImageGenericInterpolateImageFunction.h"

#include <itkComposeDisplacementFieldsImageFilter.h>

#include <itkAbsoluteValueDifferenceImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkMultiplyImageFilter.h>

static char help[] = "Solves AdLem model. Equations solved: "
    " --------------------------------\n"
    " mu div(u) + grad(p)	= (mu + lambda) grad(a)\n"
    " div(u) + k.p		= -a\n"
    " --------------------------------\n"
    " k = 0 in brain tissue. Value of k is either 0 in CSF or non-zero, can be chosen by the user.\n"
    " Choosing non-zero k relaxes the Incompressibility Constraint (IC) in CSF.\n"
    "Arguments: \n"
    "-muFile			: Lame parameter mu image file. Use this if mu is spatially varying\n"
    "         If piecewise constant in tissue and in CSF, use -parameters instead.\n\n"
    "-parameters		: Lame parameters separated by comma WITHOUT SPACE in the "
    "following order		:\n"
    "             muBrain,muCsf,lambdaBrain,lambdaCsf. E.g. -parameters "
    "2.4,4.5,1.2,3\n\n"
    "--div12pt_stencil           : If given, uses 12 point Stencil for divergence in Staggered grid. This stencil is compatible "
    "with the way most atrophy measurement tools compute divergence from the transformation field in image space. Thus use this "
    "option if you want the divergence of the velocity field obtained from the model has to exactly match the divergenc/atrophy "
    "you prescribe when computed with following stencil: (vx(i+1,j)-vx(i-1,j))/(2*hx) + (vy(i,j+1)-vy(i,j-1)/(2*hy)). \n"
    "If NOT provided, the output velocity field will have divergence compatible to the one computed internally in the staggered grid."
    "--no_lame_in_rhs           : If given, momentum equation will be mu div(u) + grad(p) = grad(a). So lambda is irrelevant and mu only in LHS."
    "-boundary_condition	: Possible values: dirichlet_at_walls dirichlet_at_skull.\n\n"
    "--relax_ic_in_csf		: If given, relaxes IC, that is non-zero k. If not given, modifies"
    "\natrophy map distributing uniform volume change in CSF to compensate global volume change in"
    "\n other regions. \n\n"
    "-relax_ic_coeff		: value of k. If not given when using --relax_ic_in_csf, reciprocal of lambda is used.\n\n"
    "--zero_vel_at_falx		: If given, makes Falx cerebri labeled voxels rigid and won't move.\n\n"
    "--sliding_at_falx		: If given, sets sliding boundary condition at Falx cerebri. Must also provide -falx_zero_vel_dir if this option is provided.\n\n"
    "-falx_zero_vel_dir	        : Possible values: 0, 1 or 2. Relevant only when --sliding_at_falx used. Sets the vel in the given dir as zero in the voxels labels as Falx Cerebri in the given maskFile. \n\n"
    "-atrophyFile		: Filename of a valid existing atrophy file that prescribes desired volume change.\n\n"
    "-maskFile			: Segmentation file that segments the image into CSF, tissue and non-brain regions.\n\n"
    "-imageFile			: Input image filename.\n\n"
    "-domainRegion		: Origin (in image coordinate => integer values) and size (in image coord => integer values) "
    "of the image region selected as computational domain.\n"
    "    x y z sx sy sz e.g. '0 0 0 30 40 50' Selects the region with origin at (0, 0, 0) and size (30, 40, 50) \n"
    "    If not provided uses full image regions.\n\n"
    "--invert_field_to_warp	: If given, inverts the obtained displacement field to warp the baseline image. This means the output field from the model is considered to be taking a point in baseline to follow-up. "
    "Otherwise the field  is assumed to be taking a point in follow-up to baseline and hence when warping the baseline image does not invert the field to perform warping..\n\n"
    "--useTensorLambda		: true or false. If true must provide a DTI image for lame parameter lambda.\n\n"
    "-lambdaFile		: filename of the DTI lambda-value image. Used when -useTensorlambda is true.\n\n"
    "-numOfTimeSteps		: number of time-steps to run the model.\n\n"
    "-resPath			: Path where all the results are to be placed.\n\n"
    "-resultsFilenamesPrefix	: Prefix to be added to all output files.\n\n"
    "--writePressure		: If given, writes the pressure image file output.\n\n"
    "--writeForce		: If given, writes the force image file output.\n\n"
    "--writeResidual		: If given, writes the residual image file output.\n\n"
    ;

struct UserOptions {
    std::string atrophyFileName, maskFileName, lambdaFileName, muFileName;
    std::string baselineImageFileName;	//used only when debug priority is highest.

    unsigned int	domainOrigin[3], domainSize[3]; //Currently not taken from commaind line!
    bool		isDomainFullSize;

    std::string boundaryCondition;
    bool	div12ptStencil, noLameInRhs;
    float	lameParas[4];	//muBrain, muCsf, lambdaBrain, lambdaCsf
    bool	relaxIcInCsf, zeroVelAtFalx, slidingAtFalx;
    float	relaxIcCoeff;	//compressibility coefficient k for CSF region.
    int		falxZeroVelDir; //Component of the velocity to be set to zero in the Falx sliding boundary condition.
    bool        useTensorLambda, isMuConstant, invertFieldToWarp;
    int         numOfTimeSteps;

    std::string resultsPath;    // Directory where all the results will be stored.
    std::string resultsFilenamesPrefix;	// Prefix for all the filenames of the results to be stored in the resultsPath.
    bool        writePressure, writeForce, writeResidual;
};


#undef __FUNCT__
#define __FUNCT_	"opsParser"
int opsParser(UserOptions &ops) {
/*
  Parse the options provided by the user and set relevant UserOptions variables.
*/
    PetscErrorCode	ierr;
    PetscBool	optionFlag = PETSC_FALSE;
    char		optionString[PETSC_MAX_PATH_LEN];
    ierr = PetscOptionsGetString(NULL,"-muFile",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    if(optionFlag) {	//Use image option only when mu is not even piecewise constant.
	ops.muFileName	   = optionString;
	ops.isMuConstant   = false;
    }
    else{//Provide parameters from options when mu is piecewise consant.
	ops.muFileName = "";
	ops.isMuConstant = true;
    }

    ierr = PetscOptionsGetString(NULL,"-parameters",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    if(!optionFlag){ //If parameters not given.
	if(ops.muFileName.empty()) throw "Must provide -parameters when -muFile not given\n";
	//Set all paras to 1 but this isn't supposed to be used by the program.
	//rather values should be used from muImageFileName.
	for(int i=0; i<4; ++i){
	    ops.lameParas[i] = 1.;
	}
    }
    else {//Set values that will be used instead of the image.
	std::stringstream parStream(optionString);
	char dummy; //for comma
	for(int i=0; i<4; ++i){
	    parStream >> ops.lameParas[i] >> dummy;
	}
    }

    // --------- Set divergence Stencil option
    ierr = PetscOptionsGetString(NULL,"--div12pt_stencil",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    ops.div12ptStencil = (bool)optionFlag;

    // --------- Set momentum equatin RHS option
    ierr = PetscOptionsGetString(NULL,"--no_lame_in_rhs",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    ops.noLameInRhs = (bool)optionFlag;

    // ---------- Set boundary condition option
    ierr = PetscOptionsGetString(NULL,"-boundary_condition",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    if(!optionFlag) throw "Must provide a valid boundary condition. e.g. dirichlet_at_skull";
    ops.boundaryCondition = optionString;

    ierr = PetscOptionsGetString(NULL,"--relax_ic_in_csf",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    ops.relaxIcInCsf = (bool)optionFlag;

    PetscReal optionReal;
    ierr = PetscOptionsGetReal(NULL, "-relax_ic_coeff", &optionReal, &optionFlag);CHKERRQ(ierr);
    if(optionFlag) ops.relaxIcCoeff = (float)optionReal;
    else ops.relaxIcCoeff = 0.;

    ierr = PetscOptionsGetString(NULL,"--zero_vel_at_falx",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    ops.zeroVelAtFalx = (bool)optionFlag;

    ierr = PetscOptionsGetString(NULL,"--sliding_at_falx",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    ops.slidingAtFalx = (bool)optionFlag;

    if(ops.slidingAtFalx)
    {
	if(ops.zeroVelAtFalx) throw "--zero_vel_at_falx and --sliding_at_falx are mutually exlclusive";
	else
	{
	    PetscInt optionInt;
	    ierr = PetscOptionsGetInt(NULL, "-falx_zero_vel_dir", &optionInt, &optionFlag);CHKERRQ(ierr);
	    if(optionFlag) ops.falxZeroVelDir = (int)optionInt;
	    else throw "must provide -falx_zero_vel_dir when --sliding_at_falx is used.";
	}
    }

    // ---------- Set input image files, computational region and other options.
    ierr = PetscOptionsGetString(NULL,"-atrophyFile",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    if(optionFlag) ops.atrophyFileName = optionString;

    ierr = PetscOptionsGetString(NULL,"-maskFile",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    if(optionFlag) ops.maskFileName = optionString;

    ierr = PetscOptionsGetString(NULL,"-imageFile",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    if(optionFlag) ops.baselineImageFileName = optionString;

    ierr = PetscOptionsGetString(NULL,"-domainRegion",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    if(optionFlag) {
	ops.isDomainFullSize = false;
	std::stringstream regionStream(optionString);
	for(int i=0; i<3; ++i) regionStream >> ops.domainOrigin[i];
	for(int i=0; i<3; ++i) regionStream >> ops.domainSize[i];
    }else ops.isDomainFullSize = true;

    ierr = PetscOptionsGetString(NULL,"--invert_field_to_warp",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    ops.invertFieldToWarp = (bool)optionFlag;

    ierr = PetscOptionsGetString(NULL,"--useTensorLambda",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    ops.useTensorLambda = (bool)optionFlag;

    ierr = PetscOptionsGetInt(NULL,"-numOfTimeSteps",&ops.numOfTimeSteps,&optionFlag);CHKERRQ(ierr);
    if(!optionFlag) {
	ops.numOfTimeSteps = 1;
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n Using default number of steps: 1 since -numOfTimeSteps option was not used.\n");
    } else {
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n Model will be run for %d time steps\n", ops.numOfTimeSteps);
    }

    ierr = PetscOptionsGetString(NULL,"-lambdaFile",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    if (ops.useTensorLambda) {
	if(!optionFlag) throw "Must provide valid tensor image filename when using --useTensorLambda.\n";
	ops.lambdaFileName = optionString;
    }
    else
	if(optionFlag) throw "-lambdaFile option can be used only when --useTensorLambda is provided.\n";

    // ---------- Set output path and prefixes.
    ierr = PetscOptionsGetString(NULL,"-resPath",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    if(!optionFlag) throw "Must provide a valid path with -resPath option: e.g. -resPath ~/results";
    ops.resultsPath = optionString;

    ierr = PetscOptionsGetString(NULL,"-resultsFilenamesPrefix",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    if(!optionFlag) throw "Must provide a prefix for output filenames: e.g. -resultsFilenamesPrefix step1";
    ops.resultsFilenamesPrefix = optionString;

    // ---------- Set the choice of the output files to be written.
    ierr = PetscOptionsGetString(NULL,"--writePressure",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    ops.writePressure = (bool)optionFlag;
    ierr = PetscOptionsGetString(NULL,"--writeForce",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    ops.writeForce = (bool)optionFlag;
    ierr = PetscOptionsGetString(NULL,"--writeResidual",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
    ops.writeResidual = (bool)optionFlag;
    return 0;

}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{

    const unsigned int DIM = 3;
    std::vector<double> wallVelocities(18);
    //TODO: Set from the user when dirichlet_at_walls boundary condition is used.
    /*0,1,2,		//south wall
      3,4,5,			//west wall
      6,7,8,			//north wall
      9,10,11,			//east wall
      12,13,14,		//front wall
      15,16,17			//back wall*/
    //    unsigned int wallPos		 = 6;
    //        wallVelocities.at(wallPos) = 1;

    //Define types:
    typedef AdLem3D<DIM>::ScalarImageType	ScalarImageType;
    typedef AdLem3D<DIM>::IntegerImageType	IntegerImageType;
    typedef AdLem3D<DIM>::VectorImageType	VectorImageType;
    typedef AdLem3D<DIM>::ScalarImageReaderType ScalarImageReaderType;
    typedef AdLem3D<DIM>::ScalarImageWriterType ScalarImageWriterType;
    typedef AdLem3D<DIM>::VectorImageWriterType VectorImageWriterType;

    PetscInitialize(&argc,&argv,(char*)0,help);
    {
	// ---------- Get user inputs
	UserOptions	ops;
	try {
	    opsParser(ops);
	}
	catch (const char* msg){
	    std::cerr<<msg<<std::endl;
	    return EXIT_FAILURE;
	}
        // ---------- Read input baseline image
        ScalarImageType::Pointer baselineImage = ScalarImageType::New();
        {
            ScalarImageReaderType::Pointer   imageReader = ScalarImageReaderType::New();
            imageReader->SetFileName(ops.baselineImageFileName);
            imageReader->Update();
            baselineImage = imageReader->GetOutput();
        }
        // ---------- Read input baseline brainMask image
        IntegerImageType::Pointer baselineBrainMask = IntegerImageType::New();
        {
            AdLem3D<DIM>::IntegerImageReaderType::Pointer   imageReader = AdLem3D<DIM>::IntegerImageReaderType::New();
            imageReader->SetFileName(ops.maskFileName);
            imageReader->Update();
            baselineBrainMask = imageReader->GetOutput();
        }
        // ---------- Read input atrophy map
        ScalarImageType::Pointer baselineAtrophy = ScalarImageType::New();
        {
            ScalarImageReaderType::Pointer   imageReader = ScalarImageReaderType::New();
            imageReader->SetFileName(ops.atrophyFileName);
            imageReader->Update();
            baselineAtrophy = imageReader->GetOutput();
        }

	// ---------- Set up output prefix with proper path
	std::string filesPref(ops.resultsPath+ops.resultsFilenamesPrefix);

	AdLem3D<DIM>		AdLemModel;
	try { // ---------- Set up the model parameters
	    AdLemModel.setBoundaryConditions(ops.boundaryCondition, ops.relaxIcInCsf, ops.relaxIcCoeff, ops.zeroVelAtFalx,
					     ops.slidingAtFalx, ops.falxZeroVelDir);
	    if(AdLemModel.getBcType() == AdLem3D<DIM>::DIRICHLET_AT_WALLS)
		AdLemModel.setWallVelocities(wallVelocities);
	    AdLemModel.setLameParameters(
		ops.isMuConstant, ops.useTensorLambda, ops.lameParas[0], ops.lameParas[1], ops.lameParas[2],
		ops.lameParas[3], ops.lambdaFileName, ops.muFileName);
	    AdLemModel.setBrainMask(baselineBrainMask, maskLabels::NBR, maskLabels::CSF, maskLabels::FALX_CEREBRI);
	    // ---------- Set up the atrophy map
	    AdLemModel.setAtrophy(baselineAtrophy);

	    // ---------- Set the computational region (Can be set only after setting all required images!)
	    if(ops.isDomainFullSize)
		AdLemModel.setDomainRegionFullImage();
	    else
		AdLemModel.setDomainRegion(ops.domainOrigin, ops.domainSize);
	    if (!ops.relaxIcInCsf) {
		AdLemModel.prescribeUniformExpansionInCsf();
		baselineAtrophy = AdLemModel.getAtrophyImage();
		AdLemModel.writeAtrophyToFile(filesPref + "T0AtrophyModified.nii.gz");

	    }
	} catch(const char* msg) {
	    std::cerr<<msg<<std::endl;
	    return EXIT_FAILURE;
	}

        // ---------- Define itk types required for the warping of the mask and atrophy map:
        typedef InverseDisplacementImageFilter<VectorImageType> FPInverseType;
        typedef itk::WarpImageFilter<ScalarImageType,ScalarImageType,VectorImageType> WarpFilterType;
	typedef itk::WarpImageFilter<IntegerImageType,IntegerImageType,VectorImageType> IntegerWarpFilterType;
        typedef itk::NearestNeighborInterpolateImageFunction<IntegerImageType> InterpolatorFilterNnType;
	//typedef itk::LabelImageGenericInterpolateImageFunction<ScalarImageType, itk::LinearInterpolateImageFunction> InterpolatorGllType; //General Label interpolator with linear interpolation.
        typedef itk::AbsoluteValueDifferenceImageFilter<IntegerImageType,IntegerImageType,ScalarImageType> AbsDiffImageFilterType;
        typedef itk::StatisticsImageFilter<ScalarImageType> StatisticsImageFilterType;
        typedef itk::ComposeDisplacementFieldsImageFilter<VectorImageType, VectorImageType> VectorComposerType;
        VectorImageType::Pointer composedDisplacementField; //declared outside loop because we need this for two different iteration steps.


	bool isMaskChanged(true);	//tracker flag to see if the brain mask is changed or not after the previous warp and NN interpolation.
        for (int t=1; t<=ops.numOfTimeSteps; ++t) {
            //-------------- Get the string for the current time step and add it to the prefix of all the files to be saved -----//
            std::stringstream	timeStep;
            timeStep << t;
            std::string		stepString("T"+timeStep.str());
            //------------- Modify atrophy map to adapt to the provided mask. ----------------//
	    // ---------- do the modification after the first step. That means I expect the atrophy map to be valid
	    // ---------- when input by the user. i.e. only GM/WM has atrophy and 0 on CSF and NBR regions.
            // ---------- Solve the system of equations
            AdLemModel.solveModel(ops.noLameInRhs, ops.div12ptStencil, isMaskChanged);
            // ---------- Write the solutions and residuals
            AdLemModel.writeVelocityImage(filesPref+stepString+"vel.nii.gz");
	    if(!ops.div12ptStencil) //Div computation from within Adlem3d supported only for 9 point div stencil.
		AdLemModel.writeDivergenceImage(filesPref+stepString+"div.nii.gz");
            if (ops.writeForce) AdLemModel.writeForceImage(filesPref+stepString+"force.nii.gz");
            if (ops.writePressure) AdLemModel.writePressureImage(filesPref+stepString+"press.nii.gz");
            if (ops.writeResidual) AdLemModel.writeResidual(filesPref+stepString);
	    VectorImageType::Pointer currentDisplacementField = AdLemModel.getVelocityImage();
	    if(ops.invertFieldToWarp)
	    {// Invert the current displacement field to create warping field
		FPInverseType::Pointer inverter = FPInverseType::New();
		inverter->SetInput(AdLemModel.getVelocityImage());
		inverter->SetErrorTolerance(1e-1);
		inverter->SetMaximumNumberOfIterations(50);
		inverter->Update();
		PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n Displacement field inversion: tolerance not reached in %d voxels \n\n",
					inverter->GetNumberOfErrorToleranceFailures());
		currentDisplacementField = inverter->GetOutput();
	    }
            if(t == 1) composedDisplacementField = currentDisplacementField;
            else
	    { // Compose the velocity field.
                VectorComposerType::Pointer vectorComposer = VectorComposerType::New();
                vectorComposer->SetDisplacementField(currentDisplacementField);
                vectorComposer->SetWarpingField(composedDisplacementField);
                vectorComposer->Update();
                composedDisplacementField = vectorComposer->GetOutput();
            }
	    // ---------- Warp the baseline image with the composed field with BSpline interpolation
	    typedef itk::BSplineInterpolateImageFunction<ScalarImageType> InterpolatorFilterType;
	    InterpolatorFilterType::Pointer interpolatorFilter = InterpolatorFilterType::New();
	    interpolatorFilter->SetSplineOrder(3);
	    WarpFilterType::Pointer baselineWarper = WarpFilterType::New();
	    baselineWarper->SetDisplacementField(composedDisplacementField);
	    baselineWarper->SetInterpolator(interpolatorFilter);
	    baselineWarper->SetInput(baselineImage);
	    baselineWarper->SetOutputSpacing(baselineImage->GetSpacing());
	    baselineWarper->SetOutputOrigin(baselineImage->GetOrigin());
	    baselineWarper->SetOutputDirection(baselineImage->GetDirection());
	    baselineWarper->Update();
	    ScalarImageWriterType::Pointer imageWriter = ScalarImageWriterType::New();
	    imageWriter->SetFileName(filesPref + "WarpedImageBspline" + stepString+ ".nii.gz"); //step at the end facilitate external tools to combine images later into 4D.
	    imageWriter->SetInput(baselineWarper->GetOutput());
	    imageWriter->Update();

            if(ops.numOfTimeSteps > 1)
	    { // Prepare brain mask and atrophy map for next step by warping them with current composed displacement field.
                // ---------- Warp baseline brain mask with an itk warpFilter, nearest neighbor.
                IntegerWarpFilterType::Pointer brainMaskWarper = IntegerWarpFilterType::New();
                brainMaskWarper->SetDisplacementField(composedDisplacementField);
                brainMaskWarper->SetInput(baselineBrainMask);
                brainMaskWarper->SetOutputSpacing(baselineBrainMask->GetSpacing());
                brainMaskWarper->SetOutputOrigin(baselineBrainMask->GetOrigin());
                brainMaskWarper->SetOutputDirection(baselineBrainMask->GetDirection());

                InterpolatorFilterNnType::Pointer nnInterpolatorFilter = InterpolatorFilterNnType::New();
                brainMaskWarper->SetInterpolator(nnInterpolatorFilter);
                brainMaskWarper->Update();

                // ---------- Compare warped mask with the previous mask
                AbsDiffImageFilterType::Pointer absDiffImageFilter = AbsDiffImageFilterType::New();
                absDiffImageFilter->SetInput1(AdLemModel.getBrainMaskImage());
                absDiffImageFilter->SetInput2(brainMaskWarper->GetOutput());
                absDiffImageFilter->Update();

                StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
                statisticsImageFilter->SetInput(absDiffImageFilter->GetOutput());
                statisticsImageFilter->Update();
                if(statisticsImageFilter->GetSum() < 1) {
                    isMaskChanged = false;
                } else {
                    isMaskChanged = true;
                    AdLemModel.setBrainMask(brainMaskWarper->GetOutput(), maskLabels::NBR, maskLabels::CSF, maskLabels::FALX_CEREBRI);
                    AdLemModel.writeBrainMaskToFile(filesPref+stepString+"Mask.nii.gz");
                }

                // ---------- Warp baseline atrophy with an itk WarpFilter, linear interpolation; using composed field.
                WarpFilterType::Pointer atrophyWarper = WarpFilterType::New();
                atrophyWarper->SetDisplacementField(composedDisplacementField);
                atrophyWarper->SetInput(baselineAtrophy);
                atrophyWarper->SetOutputSpacing(baselineAtrophy->GetSpacing());
                atrophyWarper->SetOutputOrigin(baselineAtrophy->GetOrigin());
                atrophyWarper->SetOutputDirection(baselineAtrophy->GetDirection());
                atrophyWarper->Update();
                AdLemModel.setAtrophy(atrophyWarper->GetOutput());
		//AdLemModel.writeAtrophyToFile(filesPref+stepString+"AtrophyWarpedNotModified.nii.gz"); //Useful to see
		// how i) warping  ii) modifying affects the total atrophy in the image.
                //Atrophy present at the newly created CSF regions are redistributed to the nearest GM/WM tissues voxels.
		// And in CSF put the values as the ops.relaxIcInCsf dictates.
	        AdLemModel.modifyAtrophy(maskLabels::CSF, 0, true, ops.relaxIcInCsf);
		//AdLemModel.modifyAtrophy(maskLabels::CSF,0,false, ops.relaxIcInCsf); //no redistribution.
		AdLemModel.modifyAtrophy(maskLabels::NBR,0);  //set zero atrophy at non-brain region., don't change values elsewhere.
		AdLemModel.writeAtrophyToFile(filesPref+stepString+"AtrophyModified.nii.gz");

            }
        }
	if(ops.numOfTimeSteps > 1) //Write composed field only if num_of_time_steps > 1
	{
	    VectorImageWriterType::Pointer   displacementWriter = VectorImageWriterType::New();
	    displacementWriter->SetFileName(filesPref+"ComposedField.nii.gz");
	    displacementWriter->SetInput(composedDisplacementField);
	    displacementWriter->Update();
	}
    }
    PetscErrorCode ierr;
    ierr = PetscFinalize();CHKERRQ(ierr);
    return 0;
}
