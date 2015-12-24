#include "AdLem3D.hxx"
#include "GlobalConstants.hxx"
#include"PetscAdLemTaras3D.hxx"
#include <itkImageRegionIteratorWithIndex.h>
#include <itkMultiplyImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkLabelStatisticsImageFilter.h>
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkExtractImageFilter.h"
#include<iostream>

#undef __FUNCT__
#define __FUNCT__ "AdLem3D"
AdLem3D::AdLem3D():mWallVelocities(18)
{
    mIsBrainMaskSet	  = false;
    mIsAtrophySet	  = false;
    mIsMuImageSet	  = false;
    mIsLambdaImageSet	  = false;

    mIsBcSet		  = false;
    mPetscSolverTarasUsed = false;
    mRelaxIcPressureCoeff = 0;  //This default changed only when setting brain mask.

    // number of times the solver is called.
    mNumOfSolveCalls		= 0;
    //results allocation track variables
    mVelocityAllocated		= false;
    mPressureAllocated		= false;
    mForceAllocated		= false;
    mDivergenceAllocated	= false;

    mVelocityLatest		= false;
    mPressureLatest		= false;
    mForceLatest		= false;
    mDivergenceLatest		= false;
}

#undef __FUNCT__
#define __FUNCT__ "~AdLem3D"
AdLem3D::~AdLem3D()
{
    if(mPetscSolverTarasUsed) delete mPetscSolverTaras;
}


#undef __FUNCT__
#define __FUNCT__ "setBoundaryConditions"
void AdLem3D::setBoundaryConditions(const std::string& boundaryCondition, bool relaxIcInCsf, float relaxIcPressureCoeff) {
    if((boundaryCondition.compare("dirichlet_at_walls") != 0) &&
	   (boundaryCondition.compare("dirichlet_at_skull") != 0))
	    throw "Invalid boundary condition.";
    else {
	if(boundaryCondition.compare("dirichlet_at_skull") == 0)
	    mBc = DIRICHLET_AT_SKULL;
	else
	    mBc = DIRICHLET_AT_WALLS;
    }
    if(!relaxIcInCsf)
	if(fabs(relaxIcPressureCoeff) > 1e-6)
	    throw "cannot set non-zero relaxIcPressureCoeff when relaxIcIncsf is false!";
    mRelaxIcInCsf = relaxIcInCsf;
    mRelaxIcPressureCoeff = relaxIcPressureCoeff;
}

#undef __FUNCT__
#define __FUNCT__ "setWallVelocities"
void AdLem3D::setWallVelocities(std::vector<double>& wallVelocities) {
    if(mBc != DIRICHLET_AT_WALLS)
	throw "setWallVelocities can be called only when dirichlet_at_walls boundary condition is set.";
    mWallVelocities = wallVelocities;
}

#undef __FUNCT__
#define __FUNCT__ "setLameParameters"
void AdLem3D::setLameParameters(bool isMuConstant, bool useTensorLambda,
                                double muBrain, double muCsf,
                                double lambdaBrain, double lambdaCsf,
				std::string lambdaImageFile, std::string muImageFile)
{
    mIsMuConstant	= isMuConstant;
    mUseMuImage		=  !isMuConstant; //Currently constant mu means don't use mu image.
    mUseTensorLambda	= useTensorLambda;
    mMuBrain		= muBrain;
    mMuCsf		= muCsf;
    mLambdaBrain	= lambdaBrain;
    mLambdaCsf		= lambdaCsf;
    if(!isMuConstant) { // That is if not even piecewise constant, use image.
        ScalarImageReaderType::Pointer   scalarImageReader = ScalarImageReaderType::New();
        scalarImageReader->SetFileName(muImageFile);
        scalarImageReader->Update();
        setMu(scalarImageReader->GetOutput());
    }
    if(useTensorLambda) {
        TensorImageReaderType::Pointer   tensorImageReader = TensorImageReaderType::New();
        tensorImageReader->SetFileName(lambdaImageFile);
        tensorImageReader->Update();
        setLambda(tensorImageReader->GetOutput());
    }
}

#undef __FUNCT__
#define __FUNCT__ "setMu"
void AdLem3D::setMu(ScalarImageType::Pointer inputMu)
{
    mMu = inputMu;
    mIsMuImageSet = true;
}


#undef __FUNCT__
#define __FUNCT__ "setLambda"
void AdLem3D::setLambda(TensorImageType::Pointer inputLambda)
{
    mLambda = inputLambda;
    mIsLambdaImageSet = true;
}


#undef __FUNCT__
#define __FUNCT__ "setBrainMask"
int AdLem3D::setBrainMask(std::string maskImageFile, int skullLabel, int relaxIcLabel)
{
    if(!mIsBcSet)
	throw "Boundary conditions must be set before setting brain mask.";
    IntegerImageReaderType::Pointer   imageReader = IntegerImageReaderType::New();
    imageReader->SetFileName(maskImageFile);
    imageReader->Update();
    setBrainMask(imageReader->GetOutput(),skullLabel, relaxIcLabel);
    return 0;
}

#undef __FUNCT__
#define __FUNCT__ "setBrainMask"
int AdLem3D::setBrainMask(IntegerImageType::Pointer brainMask, int skullLabel, int relaxIcLabel)
{
    if(mIsBcSet)
	throw "Boundary conditions must be set before setting brain mask.";
    mBrainMask			= brainMask;
    mRelaxIcLabel		= relaxIcLabel;
    mSkullLabel			= skullLabel;
    mIsBrainMaskSet		= true;
    return 0;
}

#undef __FUNCT__
#define __FUNCT__ "getBrainMaskImage"
AdLem3D::IntegerImageType::Pointer AdLem3D::getBrainMaskImage()
{
    return mBrainMask;
}

#undef __FUNCT__
#define __FUNCT__ "setDomainRegionFullImage"
void AdLem3D::setDomainRegionFullImage()
{
    if(!mIsBrainMaskSet || !mIsAtrophySet)
        throw "must first set brain mask and atrophy image before setting domainRegion.\n";
    if(mUseMuImage && (!mIsMuImageSet))
	throw "for the given choices, mu image is expected to be set. Domain region can be set only after setting this image!\n";
    if(mUseTensorLambda && (!mIsLambdaImageSet))
	throw "for the given choices, lambda image is expected to be set. Domain region can be set only after setting this image!\n";
    mDomainRegion = mBrainMask->GetLargestPossibleRegion();
    // No need to extract when using the full image.
}

#undef __FUNCT__
#define __FUNCT__ "setDomainRegion"
void AdLem3D::setDomainRegion(unsigned int origin[], unsigned int size[])
{
    // ---------- Set the region
    ScalarImageType::IndexType domainOrigin;
    ScalarImageType::SizeType domainSize;
    for(int i=0;i<3;++i) {     //FIXME: ADD error-guard to ensure that size
        domainSize.SetElement(i,size[i]);     //size!!
        domainOrigin.SetElement(i,origin[i]); //provided is less or equal to the input
    }
    mDomainRegion.SetIndex(domainOrigin);
    mDomainRegion.SetSize(domainSize);
    //    std::cout<<"domain size: "<<mDomainRegion.GetSize()<<std::endl;

    // ---------- Extract the region from all the images and reset them to contain only this
    // ---------- region.
    typedef itk::ExtractImageFilter<IntegerImageType, IntegerImageType> ExtractIntegerImageFilterType;
    typedef itk::ExtractImageFilter<ScalarImageType, ScalarImageType> ExtractScalarImageFilterType;
    typedef itk::ExtractImageFilter<TensorImageType, TensorImageType> ExtractTensorImageFilterType;

    if(mIsBrainMaskSet) {
	ExtractIntegerImageFilterType::Pointer imageExtracter = ExtractIntegerImageFilterType::New();
	imageExtracter->SetInput(mBrainMask);
	imageExtracter->SetExtractionRegion(mDomainRegion);
	imageExtracter->Update();
	mBrainMask = imageExtracter->GetOutput();
    }else
        throw "must set brain mask image before setting domain region.\n";

    if(mIsAtrophySet){
	ExtractScalarImageFilterType::Pointer imageExtracter = ExtractScalarImageFilterType::New();
	imageExtracter->SetInput(mAtrophy);
	imageExtracter->SetExtractionRegion(mDomainRegion);
	imageExtracter->Update();
	mAtrophy = imageExtracter->GetOutput();
	//FIXME: 1. Do I need disconnectPipeline() here ? 2. DirectionCollapse ?
    }else
	throw "must set atrophy image before setting domain region.\n";

    if(mUseMuImage){
	if(mIsMuImageSet){
	    //std::cout<<"in muImage extract image"<<std::endl;
	    ExtractScalarImageFilterType::Pointer imageExtracter = ExtractScalarImageFilterType::New();
	    imageExtracter->SetInput(mMu);
	    imageExtracter->SetExtractionRegion(mDomainRegion);
	    imageExtracter->Update();
	    mMu = imageExtracter->GetOutput();
	}else
	    throw "for the given choices, mu image is expected to be set. Domain region can be set only after setting this image!\n";
    }

    if(mUseTensorLambda){
	if(mIsLambdaImageSet){
	    //std::cout<<"in lamda extract image"<<std::endl;
	    ExtractTensorImageFilterType::Pointer imageExtracter = ExtractTensorImageFilterType::New();
	    imageExtracter->SetInput(mLambda);
	    imageExtracter->SetExtractionRegion(mDomainRegion);
	    imageExtracter->Update();
	    mLambda = imageExtracter->GetOutput();
	}else
	    throw "for the given choices, lambda image is expected to be set. Domain region can be set only after setting this image!\n";
    }
}

#undef __FUNCT__
#define __FUNCT__ "getBcType"
AdLem3D::bcType AdLem3D::getBcType() const { return mBc; }

#undef __FUNCT__
#define __FUNCT__ "isMuConstant"
bool AdLem3D::isMuConstant() const
{
    return mIsMuConstant;
}

#undef __FUNCT__
#define __FUNCT__ "isLambdaTensor"
bool AdLem3D::isLambdaTensor() const
{
    return mUseTensorLambda;
}


#undef __FUNCT__
#define __FUNCT__ "getWallVelocities"
void AdLem3D::getWallVelocities(std::vector<double>& wallVelocities) {
    wallVelocities = mWallVelocities;
}

#undef __FUNCT__
#define __FUNCT__ "relaxIcInCsf"
bool AdLem3D::relaxIcInCsf()
{
    return mRelaxIcInCsf;
}


#undef __FUNCT__
#define __FUNCT__ "getRelaxIcPressureCoeff"
float AdLem3D::getRelaxIcPressureCoeff()
{
    return mRelaxIcPressureCoeff;
}

#undef __FUNCT__
#define __FUNCT__ "getSkullLabel"
int AdLem3D::getSkullLabel()
{
    return mSkullLabel;
}

#undef __FUNCT__
#define __FUNCT__ "muAt"
double AdLem3D::muAt(int x, int y, int z) const
{
    if (mIsMuConstant) {
	if (brainMaskAt(x, y, z) == maskLabels::CSF)
	    return mMuCsf;
	else
	    return mMuBrain;
    }
    else {
	ScalarImageType::IndexType pos;
	pos.SetElement(0, mMu->GetLargestPossibleRegion().GetIndex()[0] + x);
        pos.SetElement(1, mMu->GetLargestPossibleRegion().GetIndex()[1] + y);
        pos.SetElement(2, mMu->GetLargestPossibleRegion().GetIndex()[2] + z);
        return (mMu->GetPixel(pos));
    }
}

#undef __FUNCT__
#define __FUNCT__ "lambdaAt"
double AdLem3D::lambdaAt(int x, int y, int z,
                         unsigned int Li, unsigned int Lj) const
{
    if (mUseTensorLambda) {
        TensorImageType::IndexType pos;
        pos.SetElement(0, mLambda->GetLargestPossibleRegion().GetIndex()[0] + x);
        pos.SetElement(1, mLambda->GetLargestPossibleRegion().GetIndex()[1] + y);
        pos.SetElement(2, mLambda->GetLargestPossibleRegion().GetIndex()[2] + z);
        return (mLambda->GetPixel(pos)(Li,Lj));
    }
    // If the model is initialized for scalar lambda then it is same as
    // lamda times identity matrix. So all non-diagonal elements are zero
    // while diagonal elements are the scalar value.
    if (Li == Lj) {
	if (brainMaskAt(x, y, z) == maskLabels::CSF){
	    return mLambdaCsf;
	}
	else
	    return mLambdaBrain;

    }
    else
        return 0.0;
}

#undef __FUNCT__
#define __FUNCT__ "aAt"
double AdLem3D::aAt(int x, int y, int z) const
{
    ScalarImageType::IndexType pos;
    pos.SetElement(0, mAtrophy->GetLargestPossibleRegion().GetIndex()[0] + x);
    pos.SetElement(1, mAtrophy->GetLargestPossibleRegion().GetIndex()[1] + y);
    pos.SetElement(2, mAtrophy->GetLargestPossibleRegion().GetIndex()[2] + z);
    return(mAtrophy->GetPixel(pos));

}

#undef __FUNCT__
#define __FUNCT__ "brainMaskAt"
double AdLem3D::brainMaskAt(int x, int y, int z) const
{
    ScalarImageType::IndexType pos;
    pos.SetElement(0, mBrainMask->GetLargestPossibleRegion().GetIndex()[0] + x);
    pos.SetElement(1, mBrainMask->GetLargestPossibleRegion().GetIndex()[1] + y);
    pos.SetElement(2, mBrainMask->GetLargestPossibleRegion().GetIndex()[2] + z);

    return(mBrainMask->GetPixel(pos));

}

#undef __FUNCT__
#define __FUNCT__ "getRelaxIcLabel"
int AdLem3D::getRelaxIcLabel() const
{
/*
  Return a label of the region in brainMask for which the Incompressibility Constraint (IC) is to be relaxed.
  That is the regions where div(u) - kp = 0 with non-zero k.
*/
    return mRelaxIcLabel;
}

#undef __FUNCT__
#define __FUNCT__ "dataAt"
double AdLem3D::dataAt(std::string dType, int x, int y, int z, unsigned int Mi, unsigned int Mj)
{
    if (dType.compare("mu") == 0)
        return muAt(x,y,z);
    else if (dType.compare("lambda") == 0)
        return lambdaAt(x,y,z,Mi,Mj);
    else if (dType.compare("atrophy") == 0)
        return aAt(x,y,z);
    else
        std::cout<<"invalid option: "<<dType<<" : for function dataAt"<<std::endl;
    return 0;
}

#undef __FUNCT__
#define __FUNCT__ "getXnum"
int AdLem3D::getXnum() const
{
    return mDomainRegion.GetSize()[0];
}

#undef __FUNCT__
#define __FUNCT__ "getYnum"
int AdLem3D::getYnum() const
{
    return mDomainRegion.GetSize()[1];
}

#undef __FUNCT__
#define __FUNCT__ "getZnum"
int AdLem3D::getZnum() const
{
    return mDomainRegion.GetSize()[2];
}

#undef __FUNCT__
#define __FUNCT__ "getXspacing"
double AdLem3D::getXspacing() const
{
    return mBrainMask->GetSpacing()[0];
}

#undef __FUNCT__
#define __FUNCT__ "getYspacing"
double AdLem3D::getYspacing() const
{
    return mBrainMask->GetSpacing()[1];
}

#undef __FUNCT__
#define __FUNCT__ "getZspacing"
double AdLem3D::getZspacing() const
{
    return mBrainMask->GetSpacing()[2];
}

#undef __FUNCT__
#define __FUNCT__ "solveModel"
void AdLem3D::solveModel(bool operatorChanged)
{
    if(!mPetscSolverTarasUsed) {
        mPetscSolverTarasUsed = true;
        mPetscSolverTaras = new PetscAdLemTaras3D(this,false);
    }
    mPetscSolverTaras->solveModel(operatorChanged);
    updateStateAfterSolveCall();

}
#undef __FUNCT__
#define __FUNCT__ "updateStateAfterSolveCall"
void AdLem3D::updateStateAfterSolveCall()
{
    ++mNumOfSolveCalls;
    mVelocityLatest = false;
    mPressureLatest = false;
    mDivergenceLatest = false;
    mForceLatest = false;
}

#undef __FUNCT__
#define __FUNCT__ "setAtrophy"
void AdLem3D::setAtrophy(std::string atrophyImageFile)
{
    ScalarImageReaderType::Pointer   scalarReader = ScalarImageReaderType::New();
    scalarReader->SetFileName(atrophyImageFile);
    scalarReader->Update();
    setAtrophy(scalarReader->GetOutput());
}


#undef __FUNCT__
#define __FUNCT__ "setAtrophy"
void AdLem3D::setAtrophy(ScalarImageType::Pointer inputAtrophy)
{
    mAtrophy = inputAtrophy;
    mIsAtrophySet = true;
}

#undef __FUNCT__
#define __FUNCT__ "scaleAtrophy"
void AdLem3D::scaleAtrophy(double factor)
{
    typedef itk::MultiplyImageFilter<ScalarImageType, ScalarImageType, ScalarImageType> MultiplyImageFilterType;
    MultiplyImageFilterType::Pointer multiplyImageFilter = MultiplyImageFilterType::New();
    multiplyImageFilter->SetInput(mAtrophy);
    multiplyImageFilter->SetConstant(factor);
    multiplyImageFilter->Update();
    mAtrophy=multiplyImageFilter->GetOutput();
}

#undef __FUNCT__
#define __FUNCT__ "isAtrophyValid"
//total sum close to zero, i.e. < sumMaxValue.
//FIXME: all boundary voxels has zero ? But what to do with CSF voxels that touch the skull!
// This could be important when relaxIcIncsf is false, where compatibility condition is important.
bool AdLem3D::isAtrophySumZero(double sumMaxValue) {
    //Check if sum is zero:
    itk::ImageRegionIterator<ScalarImageType> it(mAtrophy, mAtrophy->GetLargestPossibleRegion());
    it.GoToBegin();
    double sum = 0;
    while(!it.IsAtEnd()) {
        sum+=it.Get();
        ++it;
    }
    if(std::abs(sum) > sumMaxValue) return false;

    //Check whether the boundary should have zero atrophy.
    //FIXME
    return true;

}

#undef __FUNCT__
#define __FUNCT__ "prescribeUniformExpansionInCsf"
//FIXME: Need to make it work when domainRegion is not the full image!
//i.e. Exclude boundary voxels from being considered for computing atrophy and expansion.
void AdLem3D::prescribeUniformExpansionInCsf(){//Need to check this function!
/*
    Set atrophy_in_csf  = - total_atrophy_elsewhere/num_of_csf_voxels.
      // 	//And CHECK how Dirichlet condition at the skull creates issues because CSF voxels touching
      // 	//skull voxels possibly won't get the desired expansion. But hopefully it shouldn't cause problems!
      */
    typedef itk::MaskImageFilter< ScalarImageType, IntegerImageType, ScalarImageType> MaskFilterType;
    typedef itk::StatisticsImageFilter<ScalarImageType> StatisticsImageFilterType;
    typedef itk::LabelStatisticsImageFilter<ScalarImageType, IntegerImageType> LabelStatisticsImageFilterType;
    // ---------- First set atrophy in CSF region to zero before taking total sum of atrophy.
    {
	MaskFilterType::Pointer maskFilter = MaskFilterType::New();
	maskFilter->SetInput(mAtrophy);
	maskFilter->SetMaskImage(mBrainMask);
	maskFilter->SetMaskingValue(maskLabels::CSF);
	maskFilter->SetOutsideValue(0.);
	maskFilter->Update();
	setAtrophy(maskFilter->GetOutput());
    }

    StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
    statisticsImageFilter->SetInput(mAtrophy);
    statisticsImageFilter->Update();
    float total_atrophy = statisticsImageFilter->GetSum();
    //std::cout<<"total atrophy = "<<total_atrophy<<std::endl;

    LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
    labelStatisticsImageFilter->SetInput(mAtrophy);
    labelStatisticsImageFilter->SetLabelInput(mBrainMask);
    labelStatisticsImageFilter->Update();
    unsigned int csf_voxel_count = labelStatisticsImageFilter->GetCount(maskLabels::CSF);
    //std::cout<<"total number of voxels with csf labels = "<<csf_voxel_count<<std::endl;

    float expansion = -total_atrophy / (float)csf_voxel_count;
    //std::cout<<"csf expansion values = "<<expansion<<std::endl;

    MaskFilterType::Pointer maskFilter = MaskFilterType::New();
    maskFilter->SetInput(mAtrophy);
    maskFilter->SetMaskImage(mBrainMask);
    maskFilter->SetMaskingValue(maskLabels::CSF);
    maskFilter->SetOutsideValue(expansion);
    maskFilter->Update();
    setAtrophy(maskFilter->GetOutput());
}

#undef __FUNCT__
#define __FUNCT__ "modifyAtrophy"
void AdLem3D::modifyAtrophy(int maskLabel, double maskValue, bool redistributeAtrophy, bool relaxIcInCsf) {
    redistributeAtrophy = false; //FIXME: seems there is bug when small domain regions are selected, so
    // not going in there right now, need to figure out later.
    if(redistributeAtrophy) {
	int rad = 1; 	// redistribution only on 3X3 neigborhood at present.
	typedef itk::ConstShapedNeighborhoodIterator< ScalarImageType > ConstShapedNeighborhoodIteratorType;
	typedef itk::ConstNeighborhoodIterator<	IntegerImageType > ConstNeighborhoodIteratorType;
	typedef itk::ShapedNeighborhoodIterator< ScalarImageType > ShapedNeighborhoodIteratorType;

	ConstShapedNeighborhoodIteratorType::RadiusType radius;
	radius.Fill(rad);

	typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< ScalarImageType > FaceCalculatorType;
	FaceCalculatorType faceCalculator;
	FaceCalculatorType::FaceListType faceList;
	faceList = faceCalculator( mAtrophy, mAtrophy->GetRequestedRegion(), radius );
	FaceCalculatorType::FaceListType::iterator fit;
	for ( fit=faceList.begin(); fit != faceList.end(); ++fit) {
	    ConstNeighborhoodIteratorType maskIt( radius, mBrainMask, *fit); //iterate over brain mask to find CSF voxel and its neighboring tissue voxels.
	    ShapedNeighborhoodIteratorType atrophyIt (radius, mAtrophy, *fit); //iterate over atrophy map to change neighborhood values if required.
	    for (maskIt.GoToBegin(), atrophyIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt, ++atrophyIt) {
		double eps = 1e-6;
		if ( (int)maskIt.GetCenterPixel() == maskLabels::CSF) {
		    ScalarImageType::PixelType atrophyCurrentPixel = atrophyIt.GetCenterPixel();
		    if (std::abs(atrophyCurrentPixel) > eps) {
			int adjacentTissueVolume = 0;
			atrophyIt.ClearActiveList();
			for (int z = -rad; z <= rad; z++) { //activateIndex is protected so must create offset!
			    for (int y = -rad; y <= rad; y++) {
				for (int x = -rad; x <= rad; x++) {
				    ConstShapedNeighborhoodIteratorType::OffsetType off;
				    off[0] = x; off[1] = y; off[2] = z;
				    int neighbMaskPixel = maskIt.GetPixel(off);
				    if (neighbMaskPixel == maskLabels::GM || neighbMaskPixel == maskLabels::WM) {
					++adjacentTissueVolume;
					atrophyIt.ActivateOffset(off);
				    }
				}
			    }
			}
			if (adjacentTissueVolume > 0) {
			    double inc = atrophyCurrentPixel/adjacentTissueVolume;
			    ShapedNeighborhoodIteratorType::Iterator neighbTissueIt;
			    for (neighbTissueIt = atrophyIt.Begin(); neighbTissueIt != atrophyIt.End(); neighbTissueIt++)
				neighbTissueIt.Set(neighbTissueIt.Get() + inc);
			}
		    }
		}
	    }
	}
    }
    if(relaxIcInCsf) { //If IC is relaxed, set atrophy values to maskValue in regions with label maskLabel.
        typedef itk::MaskImageFilter< ScalarImageType, IntegerImageType,
				      ScalarImageType> MaskFilterType;
        MaskFilterType::Pointer maskFilter = MaskFilterType::New();
        maskFilter->SetInput(mAtrophy);
        maskFilter->SetMaskImage(mBrainMask);
        maskFilter->SetMaskingValue(maskLabel); //The naming conventions are tricky :)
        maskFilter->SetOutsideValue(maskValue); //with itk's convention for mask label and values.
        maskFilter->Update();
        mAtrophy = maskFilter->GetOutput();
    } else //not relaxing IC => need to have uniformcsfexpansion!
	prescribeUniformExpansionInCsf();
}

#undef __FUNCT__
#define __FUNCT__ "getAtrophyImage"
AdLem3D::ScalarImageType::Pointer AdLem3D::getAtrophyImage()
{
    return mAtrophy;
}

#undef __FUNCT__
#define __FUNCT__ "writeAtrophyToFile"
void AdLem3D::writeAtrophyToFile(std::string fileName) {
    ScalarImageWriterType::Pointer writer = ScalarImageWriterType::New();
    writer->SetFileName(fileName);
    writer->SetInput(mAtrophy);
    writer->Update();
}

#undef __FUNCT__
#define __FUNCT__ "writeBrainMaskToFile"
void AdLem3D::writeBrainMaskToFile(std::string fileName) {
    IntegerImageWriterType::Pointer writer = IntegerImageWriterType::New();
    writer->SetFileName(fileName);
    writer->SetInput(mBrainMask);
    writer->Update();
}

#undef __FUNCT__
#define __FUNCT__ "getVelocityImage"
AdLem3D::VectorImageType::Pointer AdLem3D::getVelocityImage()
{
    if(mNumOfSolveCalls == 0) {
        std::cerr<<"the model is not solved yet, first solve the system to get the solution."<<std::endl;
        return NULL;
    }
    if(!mVelocityAllocated) {
        createVelocityImage();
    }
    if(!mVelocityLatest) {
	updateImages("velocity");
    }
    return mVelocity;
}

#undef __FUNCT__
#define __FUNCT__ "getPressureImage"
AdLem3D::ScalarImageType::Pointer AdLem3D::getPressureImage()
{
    if(mNumOfSolveCalls == 0) {
        std::cerr<<"the model is not solved yet, first solve the system to get the solution."<<std::endl;
        return NULL;
    }
    if(!mPressureAllocated) {
        createPressureImage();
    }
    if(!mPressureLatest) {
	updateImages("pressure");
    }
    return mPressure;
}

#undef __FUNCT__
#define __FUNCT__ "getDivergenceImage"
AdLem3D::ScalarImageType::Pointer AdLem3D::getDivergenceImage()
{
    if(mNumOfSolveCalls == 0) {
        std::cerr<<"the model is not solved yet, first solve the system to get the solution."<<std::endl;
        return NULL;
    }
    if(!mDivergenceAllocated) {
        createDivergenceImage();
    }
    if(!mDivergenceLatest) {
	updateImages("divergence");
    }
    return mDivergence;
}

#undef __FUNCT__
#define __FUNCT__ "getForceImage"
AdLem3D::VectorImageType::Pointer AdLem3D::getForceImage()
{
    if(mNumOfSolveCalls == 0) {
        std::cerr<<"the model is not solved yet, first solve the system to get the solution."<<std::endl;
        return NULL;
    }
    if(!mForceAllocated) {
        createForceImage();
    }
    if(!mForceLatest) {
	updateImages("force");
    }
    return mForce;
}


#undef __FUNCT__
#define __FUNCT__ "writeVelcoityImage"
void AdLem3D::writeVelocityImage(std::string fileName)
{
    VectorImageWriterType::Pointer   velocityWriter = VectorImageWriterType::New();
    velocityWriter->SetFileName(fileName);
    velocityWriter->SetInput(getVelocityImage());
    velocityWriter->Update();
}

#undef __FUNCT__
#define __FUNCT__ "writePressureImage"
void AdLem3D::writePressureImage(std::string fileName)
{
    ScalarImageWriterType::Pointer   pressureWriter = ScalarImageWriterType::New();
    pressureWriter->SetFileName(fileName);
    pressureWriter->SetInput(getPressureImage());
    pressureWriter->Update();
}

#undef __FUNCT__
#define __FUNCT__ "writeDivergenceImage"
void AdLem3D::writeDivergenceImage(std::string fileName)
{
    ScalarImageWriterType::Pointer divWriter = ScalarImageWriterType::New();
    divWriter->SetFileName(fileName);
    divWriter->SetInput(getDivergenceImage());
    divWriter->Update();
}

#undef __FUNCT__
#define __FUNCT__ "writeForceImage"
void AdLem3D::writeForceImage(std::string fileName)
{
    VectorImageWriterType::Pointer   forceWriter = VectorImageWriterType::New();
    forceWriter->SetFileName(fileName);
    forceWriter->SetInput(getForceImage());
    forceWriter->Update();
}


#undef __FUNCT__
#define __FUNCT__ "writeSolutionForMatlab"
//Writes in matlab format.
void AdLem3D::writeSolutionForMatlab(std::string resultsPath, bool writeSystemMatrix)
{

    std::string matSolFileName(resultsPath+"sol");
    std::string matSysFileName(resultsPath+"sys");
    std::string matSizeSysFileName(resultsPath+"size_lin_sys");
    mPetscSolverTaras->writeToMatFile(matSolFileName,writeSystemMatrix,matSysFileName);
    std::ofstream size_file;
    size_file.open(matSizeSysFileName.c_str());
    size_file<<mDomainRegion.GetSize()[0]+1<<" "<<mDomainRegion.GetSize()[1]+1<<" "
	     <<mDomainRegion.GetSize()[2]+1;
    size_file.close();
}

#undef __FUNCT__
#define __FUNCT__ "writeResidual"
void AdLem3D::writeResidual(std::string resultsPath)
{
    mPetscSolverTaras->writeResidual(resultsPath);
}

#undef __FUNCT__
#define __FUNCT__ "createVelocityImage"
void AdLem3D::createVelocityImage()
{
    mVelocity = VectorImageType::New();
    mVelocity->SetRegions(mAtrophy->GetLargestPossibleRegion());
    mVelocity->SetOrigin(mAtrophy->GetOrigin());
    mVelocity->SetSpacing(mAtrophy->GetSpacing());
    mVelocity->SetDirection(mAtrophy->GetDirection());
    mVelocity->Allocate();

    mVelocityAllocated = true;
}

#undef __FUNCT__
#define __FUNCT__ "createPressureImage"
void AdLem3D::createPressureImage()
{
    mPressure = ScalarImageType::New();
    mPressure->SetRegions(mAtrophy->GetLargestPossibleRegion());
    mPressure->SetOrigin(mAtrophy->GetOrigin());
    mPressure->SetSpacing(mAtrophy->GetSpacing());
    mPressure->SetDirection(mAtrophy->GetDirection());
    mPressure->Allocate();

    mPressureAllocated = true;
}

#undef __FUNCT__
#define __FUNCT__ "createForceImage"
void AdLem3D::createForceImage()
{
    mForce = VectorImageType::New();
    mForce->SetRegions(mAtrophy->GetLargestPossibleRegion());
    mForce->SetOrigin(mAtrophy->GetOrigin());
    mForce->SetSpacing(mAtrophy->GetSpacing());
    mForce->SetDirection(mAtrophy->GetDirection());
    mForce->Allocate();

    mForceAllocated = true;
}

#undef __FUNCT__
#define __FUNCT__ "createDivergenceImage"
void AdLem3D::createDivergenceImage()
{
    mDivergence = ScalarImageType::New();
    mDivergence->SetRegions(mAtrophy->GetLargestPossibleRegion());
    mDivergence->SetOrigin(mAtrophy->GetOrigin());
    mDivergence->SetSpacing(mAtrophy->GetSpacing());
    mDivergence->SetDirection(mAtrophy->GetDirection());
    mDivergence->Allocate();

    mDivergenceAllocated = true;
}

#undef __FUNCT__
#define __FUNCT__ "updateImages"
void AdLem3D::updateImages(const std::string& whichImage){
    if (whichImage.compare("pressure") == 0 || whichImage.compare("divergence") == 0) {
	ScalarImageType::Pointer img;
	bool isPressure; //Create bool var because perhaps faster to check bool than doing
	//string compare below inside the loop.
	if(whichImage.compare("pressure")==0){
	    isPressure = true;
	    img = mPressure;
	}else{
	    isPressure = false;
	    img = mDivergence;
	}
	ScalarImageType::RegionType roi = img->GetLargestPossibleRegion();
	//std::cout<<whichImage<<" image region: "<<roi<<std::endl;
	typedef itk::ImageRegionIterator<ScalarImageType> ImageIteratorType;
	ImageIteratorType it(img, roi);
	ScalarImageType::PixelType pix;
	unsigned int pos[3];
	it.GoToBegin();
	// printf("From (0, 0, 0) to (%lu, %lu, %lu)\n", roi.GetSize()[2], roi.GetSize()[1], roi.GetSize()[0]);
	for(unsigned int k = 0; k<roi.GetSize()[2]; ++k){
	    pos[2] = k;
	    for(unsigned int j = 0; j<roi.GetSize()[1]; ++j){
		pos[1] = j;
		for(unsigned int i = 0; i<roi.GetSize()[0]; ++i){
		    pos[0] = i;
		    if(isPressure)
			pix = mPetscSolverTaras->getSolPressureAt(pos);
		    else
			pix = mPetscSolverTaras->getDivergenceAt(pos);
		    it.Set(pix);
		    ++it;
		}
	    }
	}
	if(isPressure)
	    mPressureLatest = true;
	else
	    mDivergenceLatest = true;
    }else if (whichImage.compare("velocity") == 0 || whichImage.compare("force") == 0){
	VectorImageType::Pointer img;
	bool isVelocity; //Create bool var because perhaps faster to check bool than doing
	//string compare below inside the loop.
	if(whichImage.compare("velocity")==0){
	    isVelocity = true;
	    img = mVelocity;
	}else{
	    isVelocity = false;
	    img = mForce;
	}
	VectorImageType::RegionType roi = img->GetLargestPossibleRegion();
	typedef itk::ImageRegionIterator<VectorImageType> VectorIteratorType;
	//std::cout<<whichImage<<"  image region: "<<roi<<std::endl;
	VectorIteratorType it(img, roi);
	VectorImageType::PixelType pix;
	unsigned int pos[3];
	it.GoToBegin();
	// printf("From (0, 0, 0) to (%lu, %lu, %lu)\n", roi.GetSize()[2], roi.GetSize()[1], roi.GetSize()[0]);
	for(unsigned int k = 0; k<roi.GetSize()[2]; ++k){
	    pos[2] = k;
	    for(unsigned int j = 0; j<roi.GetSize()[1]; ++j){
		pos[1] = j;
		for(unsigned int i = 0; i<roi.GetSize()[0]; ++i){
		    pos[0] = i;
		    for(int cc = 0; cc<3; ++cc) {
			if(isVelocity)
			    pix[cc] = mPetscSolverTaras->getSolVelocityAt(pos, cc);
			else
			    pix[cc] = mPetscSolverTaras->getRhsAt(pos,cc);
                }
                it.Set(pix);
                ++it;
		}
	    }
	}
	if(isVelocity)
	    mVelocityLatest = true;
	else
	    mForceLatest = true;
    }else
	std::cout<<"invalid image type string: "<<whichImage<<" : for function updateImages"<<std::endl; //FIXME: Exception handling!
}
