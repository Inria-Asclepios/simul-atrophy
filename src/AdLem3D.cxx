#include "AdLem3D.hxx"
#include"PetscAdLemTaras3D.hxx"
#include <itkImageRegionIteratorWithIndex.h>
#include <itkMultiplyImageFilter.h>
#include <itkMaskImageFilter.h>
#include<iostream>

#undef __FUNCT__
#define __FUNCT__ "AdLem3D"
AdLem3D::AdLem3D():mWallVelocities(18)
{
    //Initialize with Dirichlet boundary condition, no other boundary condition for now.
    mBc= AdLem3D::DIRICHLET;
    mIsBrainMaskSet = false;
    mPetscSolverTarasUsed = false;
    mRelaxIcPressureCoeff = 0;  //This gets non-zero only if brain mask is set
}

#undef __FUNCT__
#define __FUNCT__ "~AdLem3D"
AdLem3D::~AdLem3D()
{
    if(mPetscSolverTarasUsed) delete mPetscSolverTaras;
}


#undef __FUNCT__
#define __FUNCT__ "setWallVelocities"
void AdLem3D::setWallVelocities(std::vector<double>& wallVelocities) {
    mWallVelocities = wallVelocities;
}

#undef __FUNCT__
#define __FUNCT__ "setLameParameters"
void AdLem3D::setLameParameters(double muCsf, double lambdaCsf,
                                bool isMuConstant,
                                double muRatio, double lambdaRatio)
{
    mMuCsf = muCsf;
    mLambdaCsf = lambdaCsf;
    if(!isMuConstant) {
        mIsMuConstant = false;
        mMuGm = muCsf*muRatio;      mMuWm = muCsf*muRatio;
        mLambdaGm = lambdaCsf*lambdaRatio;      mLambdaWm = lambdaCsf*lambdaRatio;
    } else {
        mIsMuConstant = true;
        mMuGm = muCsf;          mMuWm = muCsf;
        mLambdaGm = lambdaCsf;   mLambdaWm = lambdaCsf;
    }
}

#undef __FUNCT__
#define __FUNCT__ "setBrainMask"
void AdLem3D::setBrainMask(std::string maskImageFile, int relaxIcLabel, int relaxIcPressureCoeff, bool setSkullVelToZero, int skullLabel)
{

    IntegerImageReaderType::Pointer   imageReader = IntegerImageReaderType::New();
    imageReader->SetFileName(maskImageFile);
    imageReader->Update();
    setBrainMask(imageReader->GetOutput(),relaxIcLabel,relaxIcPressureCoeff, setSkullVelToZero, skullLabel);
}

#undef __FUNCT__
#define __FUNCT__ "setBrainMask"
void AdLem3D::setBrainMask(IntegerImageType::Pointer brainMask, int relaxIcLabel, int relaxIcPressureCoeff, bool setSkullVelToZero, int skullLabel)
{
    mBrainMask = brainMask;
    mRelaxIcLabel = relaxIcLabel;
    mRelaxIcPressureCoeff = relaxIcPressureCoeff;
    mSetSkullVelToZero = setSkullVelToZero;
    mSkullLabel = skullLabel;
    mIsBrainMaskSet = true;
}

#undef __FUNCT__
#define __FUNCT__ "setDomainRegionFullImage"
void AdLem3D::setDomainRegionFullImage()
{
    if(!mIsBrainMaskSet) {
        std::cerr<<"must first set the brain mask image to obtain the largest possible region"<<std::endl;
    }
    mDomainRegion = mBrainMask->GetLargestPossibleRegion();
}

#undef __FUNCT__
#define __FUNCT__ "setDomainRegion"
void AdLem3D::setDomainRegion(unsigned int origin[], unsigned int size[])
{
        ScalarImageType::IndexType domainOrigin;
        ScalarImageType::SizeType domainSize;
        for(int i=0;i<3;++i) {     //ADD error-guard to ensure that size
            domainSize.SetElement(i,size[i]);     //size!!
            domainOrigin.SetElement(i,origin[i]); //provided is less or equal to the input
        }
        mDomainRegion.SetIndex(domainOrigin);
        mDomainRegion.SetSize(domainSize);
    //    std::cout<<"domain size: "<<mDomainRegion.GetSize()<<std::endl;
}

#undef __FUNCT__
#define __FUNCT__ "getBcType"
AdLem3D::bcType AdLem3D::getBcType() const { return mBc; }

bool AdLem3D::isMuConstant() const
{
    return mIsMuConstant;
}

#undef __FUNCT__
#define __FUNCT__ "getWallVelocities"
void AdLem3D::getWallVelocities(std::vector<double>& wallVelocities) {
    wallVelocities = mWallVelocities;
}


#undef __FUNCT__
#define __FUNCT__ "getRelaxIcPressureCoeff"
int AdLem3D::getRelaxIcPressureCoeff()
{
    return mRelaxIcPressureCoeff;
}

#undef __FUNCT__
#define __FUNCT__ "shouldSkullVelSetToZero"
bool AdLem3D::shouldSkullVelSetToZero()
{
    return mSetSkullVelToZero;
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
    /*if ( (x > (mXnum/2. - 4)) && (x < (mXnum/2. + 4))
         && (y > (mYnum/2. - 4)) && (y < (mYnum/2. + 4))
         && (z > (mZnum/2. - 4)) && (z < (mZnum/2. + 4))) {
        return mMuGm;
    }*/
    return mMuCsf;
}

#undef __FUNCT__
#define __FUNCT__ "lambdaAt"
double AdLem3D::lambdaAt(int x, int y, int z) const
{
    /*if ( (x > (mXnum/2. - 4)) && (x < (mXnum/2. + 4))
         && (y > (mYnum/2. - 4)) && (y < (mYnum/2. + 4))
         && (z > (mZnum/2. - 4)) && (z < (mZnum/2. + 4))) {
        return mLambdaGm;
    }
    return mLambdaCsf;*/
    return mLambdaCsf;
}

#undef __FUNCT__
#define __FUNCT__ "aAt"
double AdLem3D::aAt(int x, int y, int z) const
{
    ScalarImageType::IndexType pos;
    pos.SetElement(0, mDomainRegion.GetIndex()[0] + x);
    pos.SetElement(1, mDomainRegion.GetIndex()[1] + y);
    pos.SetElement(2, mDomainRegion.GetIndex()[2] + z);

    return(mAtrophy->GetPixel(pos));

}

#undef __FUNCT__
#define __FUNCT__ "brainMaskAt"
int AdLem3D::brainMaskAt(int x, int y, int z) const
{
    ScalarImageType::IndexType pos;
    pos.SetElement(0, mDomainRegion.GetIndex()[0] + x);
    pos.SetElement(1, mDomainRegion.GetIndex()[1] + y);
    pos.SetElement(2, mDomainRegion.GetIndex()[2] + z);

    return(mBrainMask->GetPixel(pos));

}

#undef __FUNCT__
#define __FUNCT__ "getRelaxIcLabel"
int AdLem3D::getRelaxIcLabel() const
{
    return mRelaxIcLabel;
}

#undef __FUNCT__
#define __FUNCT__ "dataAt"
double AdLem3D::dataAt(std::string dType, int x, int y, int z)
{
    if (dType.compare("mu") == 0)
        return muAt(x,y,z);
    else if (dType.compare("lambda") == 0)
        return lambdaAt(x,y,z);
    else if (dType.compare("atrophy") == 0)
        return aAt(x,y,z);
    else
        std::cout<<"invalid option: "<<dType<<" : for funciton dataAt"<<std::endl;
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
#define __FUNCT__ "solveModel"
void AdLem3D::solveModel()
{
    if(!mPetscSolverTarasUsed) {
        mPetscSolverTarasUsed = true;
        mPetscSolverTaras = new PetscAdLemTaras3D(this,false);
    }
    mPetscSolverTaras->solveModel();
    createResultImages();
}

#undef __FUNCT__
#define __FUNCT__ "setAtrophy"
void AdLem3D::setAtrophy(std::string atrophyImageFile)
{
    ScalarReaderType::Pointer   scalarReader = ScalarReaderType::New();
    scalarReader->SetFileName(atrophyImageFile);
    scalarReader->Update();
    setAtrophy(scalarReader->GetOutput());
}


#undef __FUNCT__
#define __FUNCT__ "setAtrophy"
void AdLem3D::setAtrophy(ScalarImageType::Pointer inputAtrophy)
{
    mAtrophy = inputAtrophy;
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
#define __FUNCT__ "createAtrophy"
void AdLem3D::createAtrophy(unsigned int size[3])  //Must be removed and put to another class!!
{
    ScalarImageType::IndexType  start;
    start.Fill(0);
    ScalarImageType::SizeType   sizeAtrophy;
    for(int i=0; i<3; ++i)  sizeAtrophy.SetElement(i,size[i]);

    ScalarImageType::RegionType region(start,sizeAtrophy);
    mAtrophy = ScalarImageType::New();
    mAtrophy->SetRegions(region);
    mAtrophy->Allocate();
    mAtrophy->FillBuffer(0);

    unsigned int    xn = size[0];
    unsigned int    yn = size[1];
    unsigned int    zn = size[2];
    unsigned int    numOfCentreVoxels = 8;
    double          posAtrophy = 0.4;
    //Center mass loss: in centreVoxels:
    for(unsigned int i = (xn/2)-1; i<= xn/2; ++i) {
        for(unsigned int j = (yn/2)-1; j<= yn/2; ++j) {
            for(unsigned int k = (zn/2)-1; k<= zn/2; ++k) {
                ScalarImageType::IndexType  pos;
                pos[0] = i;
                pos[1] = j;
                pos[2] = k;
                mAtrophy->SetPixel(pos,posAtrophy);
            }
        }
    }
    unsigned int    numOfBorderVoxels = xn*yn*zn - (xn-2)*(yn-2)*(zn-2);
    double  negAtrophy = -1.* (posAtrophy*numOfCentreVoxels)/
            (double)(xn*yn*zn - numOfCentreVoxels - numOfBorderVoxels);

    //Uniform creation everywhere else except the borders:
    for(unsigned int i=1; i<xn-1; ++i) {
        for(unsigned int j=1; j<yn-1; ++j) {
            for(unsigned int k=1; k<zn-1; ++k) {
                ScalarImageType::IndexType pos;
                pos[0] = i;
                pos[1] = j;
                pos[2] = k;
                /*if(!((i>= (xn/2)-1) && (i<=xn/2) &&
                        (j>= (yn/2)-1) && (j<=yn/2) &&
                        (k>= (zn/2)-1) && (k<=zn/2)))*/
                if(mAtrophy->GetPixel(pos) != posAtrophy) //prolly not a good idea to have inequality for double type!
                    mAtrophy->SetPixel(pos,negAtrophy); //so use the above commented condition instead ?
            }
        }
    }

    double pixelSum = 0;
    for(unsigned int i=1; i<xn-1; ++i) {
        for(unsigned int j=1; j<yn-1; ++j) {
            for(unsigned int k=1; k<zn-1; ++k) {
                ScalarImageType::IndexType pos;
                pos[0] = i;
                pos[1] = j;
                pos[2] = k;
                pixelSum += mAtrophy->GetPixel(pos);
            }
        }
    }
    //    std::cout<<"Sum of created Atrophy: "<<pixelSum<<std::endl;
    if(!isAtrophyValid(0.000009)) {
        std::cout<<"yaha, strange, should be valid!!"<<std::endl;
    }
}

//Valid atrophy has:
//total sum close to zero, i.e. < sumMaxValue.
//all boundary voxels has zero.
#undef __FUNCT__
#define __FUNCT__ "isAtrophyValid"
bool AdLem3D::isAtrophyValid(double sumMaxValue) {
    //Check if sum is zero:
    itk::ImageRegionIterator<ScalarImageType> it(mAtrophy,mDomainRegion);
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
#define __FUNCT__ "modifyAtrophy"
void AdLem3D::modifyAtrophy(int maskLabel, double maskValue, bool makeSumZero) {
    if(!makeSumZero) {
        typedef itk::MaskImageFilter< ScalarImageType, IntegerImageType,
                ScalarImageType> MaskFilterType;
        MaskFilterType::Pointer maskFilter = MaskFilterType::New();
        maskFilter->SetInput(mAtrophy);
        maskFilter->SetMaskImage(mBrainMask);
        maskFilter->SetMaskingValue(maskLabel); //The naming conventions are tricky :)
        maskFilter->SetOutsideValue(maskValue); //with itk's convention for mask label and values.
        maskFilter->Update();
        mAtrophy = maskFilter->GetOutput();
    } else {
        ScalarImageType::RegionType innerRegion;
        ScalarImageType::SizeType offset;
        ScalarImageType::SizeType Size(mDomainRegion.GetSize());
        offset.Fill(2);
        innerRegion.SetIndex(mDomainRegion.GetIndex() + offset);
        offset.Fill(4);
        innerRegion.SetSize(Size - offset);
        itk::ImageRegionIterator<ScalarImageType> itInner(mAtrophy,innerRegion);
        double aSum = 0;
        for(itInner.GoToBegin(); !itInner.IsAtEnd(); ++itInner)
            aSum += itInner.Get();
        unsigned int numOfPixels = (Size[0]-2)*(Size[1]-2) * (Size[2]-2)
                -(Size[0]-4)*(Size[1]-4) * (Size[2]-4);
        aSum/=(-1*(double)numOfPixels);

        itk::ImageRegionIteratorWithIndex<ScalarImageType> itOuter(mAtrophy,mDomainRegion);
        for(itOuter.GoToBegin(); !itOuter.IsAtEnd(); ++itOuter) {
            ScalarImageType::IndexType pos = itOuter.GetIndex();
            if(pos[0] == mDomainRegion.GetIndex()[0] ||
                    pos[1] == mDomainRegion.GetIndex()[1] ||
                    pos[2] == mDomainRegion.GetIndex()[2] ||
                    pos[0] == mDomainRegion.GetUpperIndex()[0] ||
                    pos[1] == mDomainRegion.GetUpperIndex()[1] ||
                    pos[2] == mDomainRegion.GetUpperIndex()[2]) {
                itOuter.Set(0);
            } else if (pos[0] == (mDomainRegion.GetIndex()[0] + 1) ||
                       pos[1] == (mDomainRegion.GetIndex()[1] + 1) ||
                       pos[2] == (mDomainRegion.GetIndex()[2] + 1) ||
                       pos[0] == (mDomainRegion.GetUpperIndex()[0] - 1) ||
                       pos[1] == (mDomainRegion.GetUpperIndex()[1] - 1) ||
                       pos[2] == (mDomainRegion.GetUpperIndex()[2] - 1)
                       ) {
                itOuter.Set(aSum);
            }
        }
    }
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
    ScalarWriterType::Pointer writer = ScalarWriterType::New();
    writer->SetFileName(fileName);
    writer->SetInput(mAtrophy);
    writer->Update();
}

#undef __FUNCT__
#define __FUNCT__ "getVelocityImage"
AdLem3D::VectorImageType::Pointer AdLem3D::getVelocityImage()
{
    return mVelocity;
}

#undef __FUNCT__
#define __FUNCT__ "getPressureImage"
AdLem3D::ScalarImageType::Pointer AdLem3D::getPressureImage()
{
    return mPressure;
}



#undef __FUNCT__
#define __FUNCT__ "writeSolution"
void AdLem3D::writeSolution(std::string resultsPath, bool inMatlabFormat, bool inMatlabFormatSystemMatrix)
{

    if(inMatlabFormat) {
        std::string matSolFileName(resultsPath+"sol");
        std::string matSysFileName(resultsPath+"sys");
        std::string matSizeSysFileName(resultsPath+"size_lin_sys");
        mPetscSolverTaras->writeToMatFile(matSolFileName,inMatlabFormatSystemMatrix,matSysFileName);
        //        mPetscSolverTaras->writeToMatFile(matSolFileName,true,matSysFileName);
        std::ofstream size_file;
        size_file.open(matSizeSysFileName.c_str());
        size_file<<mDomainRegion.GetSize()[0]+1<<" "<<mDomainRegion.GetSize()[1]+1<<" "
                <<mDomainRegion.GetSize()[2]+1;
        size_file.close();
    }
    std::string velocityFileName(resultsPath+"vel.nii.gz");
    std::string pressureFileName(resultsPath+"press.nii.gz");
    std::string divergenceFileName(resultsPath+"div.nii.gz");

    VectorWriterType::Pointer   vectorWriter = VectorWriterType::New();
    vectorWriter->SetFileName(velocityFileName);
    vectorWriter->SetInput(mVelocity);
    vectorWriter->Update();

    ScalarWriterType::Pointer   scalarWriter = ScalarWriterType::New();
    scalarWriter->SetFileName(pressureFileName);
    scalarWriter->SetInput(mPressure);
    scalarWriter->Update();

    ScalarWriterType::Pointer divWriter = ScalarWriterType::New();
    divWriter->SetFileName(divergenceFileName);
    divWriter->SetInput(mDivergence);
    divWriter->Update();
}

#undef __FUNCT__
#define __FUNCT__ "writeResidual"
void AdLem3D::writeResidual(std::string resultsPath)
{
    mPetscSolverTaras->writeResidual(resultsPath);
}

#undef __FUNCT__
#define __FUNCT__ "createImageOf"
void AdLem3D::createResultImages()
{
    itk::Index<3>           outputImageStart;       //should be 0!
    ScalarImageType::RegionType domainRegion;
    for (unsigned int i=0; i<3; ++i) outputImageStart.SetElement(i,0);
    domainRegion.SetSize(mDomainRegion.GetSize());
    domainRegion.SetIndex(outputImageStart);

    mVelocity = VectorImageType::New();
    mVelocity->SetRegions(domainRegion);
    mVelocity->SetOrigin(mAtrophy->GetOrigin());
    mVelocity->SetSpacing(mAtrophy->GetSpacing());
    mVelocity->SetDirection(mAtrophy->GetDirection());
    mVelocity->Allocate();

    mPressure = ScalarImageType::New();
    mPressure->SetRegions(domainRegion);
    mPressure->SetOrigin(mAtrophy->GetOrigin());
    mPressure->SetSpacing(mAtrophy->GetSpacing());
    mPressure->SetDirection(mAtrophy->GetDirection());
    mPressure->Allocate();


    mDivergence = ScalarImageType::New();
    mDivergence->SetRegions(domainRegion);
    mDivergence->SetOrigin(mAtrophy->GetOrigin());
    mDivergence->SetSpacing(mAtrophy->GetSpacing());
    mDivergence->SetDirection(mAtrophy->GetDirection());
    mDivergence->Allocate();


    typedef itk::ImageRegionIterator<VectorImageType> VectorIteratorType;
    typedef itk::ImageRegionIterator<ScalarImageType> ScalarIteratorType;

    VectorIteratorType velocityIterator(mVelocity,mVelocity->GetLargestPossibleRegion());
    ScalarIteratorType pressureIterator(mPressure,mPressure->GetLargestPossibleRegion());
    ScalarIteratorType divergenceIterator(mDivergence,mDivergence->GetLargestPossibleRegion());
    VectorImageType::PixelType velocityPixel;
    ScalarImageType::PixelType pressurePixel;
    ScalarImageType::PixelType divergencePixel;

    unsigned int pos[3];
    unsigned int k = 0;
    velocityIterator.GoToBegin();
    pressureIterator.GoToBegin();
    divergenceIterator.GoToBegin();
    while(k<mDomainRegion.GetSize()[2] && !velocityIterator.IsAtEnd()) {
        pos[2] = k;
        unsigned int j = 0;
        while(j<mDomainRegion.GetSize()[1] && !velocityIterator.IsAtEnd()) {
            pos[1] = j;
            unsigned int i = 0;
            while(i<mDomainRegion.GetSize()[0] && !velocityIterator.IsAtEnd()) {
                pos[0] = i;
                for(int cc = 0; cc<3; ++cc) {
                    velocityPixel[cc] = mPetscSolverTaras->getSolVelocityAt(pos,cc);
                }
                velocityIterator.Set(velocityPixel);
                ++velocityIterator;
                pressurePixel = mPetscSolverTaras->getSolPressureAt(pos);
                pressureIterator.Set(pressurePixel);
                ++pressureIterator;
                divergencePixel = mPetscSolverTaras->getDivergenceAt(pos);
                divergenceIterator.Set(divergencePixel);
                ++divergenceIterator;
                ++i;
            }
            ++j;
        }
        ++k;
    }
}
