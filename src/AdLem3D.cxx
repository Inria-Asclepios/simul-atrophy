#include "AdLem3D.hxx"
#include"PetscAdLemTaras3D.hxx"
#include <itkImageRegionIteratorWithIndex.h>
#include <itkMultiplyImageFilter.h>
#include<iostream>

//Initialize with Dirichlet boundary condition, no other boundary condition for now.
AdLem3D::AdLem3D()
{
//    FIXME: Initialize the region!!
    //Ask whether full Image will be used or not, if not set a bool variable for
    //region as unInitialized?
    mBc= AdLem3D::DIRICHLET;
    mPetscSolverTarasUsed = false;
}

AdLem3D::~AdLem3D()
{
    if(mPetscSolverTarasUsed) delete mPetscSolverTaras;
}

#undef __FUNCT__
#define __FUNCT__ "setDomainRegion"
void AdLem3D::setDomainRegion(unsigned int origin[], unsigned int size[], bool fullImageSize)
{
    if(fullImageSize) {
        mDomainRegion = mAtrophy->GetLargestPossibleRegion();
    } else {
        ScalarImageType::IndexType domainOrigin;
        ScalarImageType::SizeType domainSize;
        for(int i=0;i<3;++i) {     //ADD error-guard to ensure that size
            domainSize.SetElement(i,size[i]);     //size!!
            domainOrigin.SetElement(i,origin[i]); //provided is less or equal to the input
        }
        mDomainRegion.SetIndex(domainOrigin);
        mDomainRegion.SetSize(domainSize);
    }
}

#undef __FUNCT__
#define __FUNCT__ "solveModel"
void AdLem3D::solveModel()
{
    mPetscSolverTaras = new PetscAdLemTaras3D(this,false);
    mPetscSolverTaras->solveModel();
    mPetscSolverTarasUsed = true;
}

//does not guarantee that this is a valid atrophy!
void AdLem3D::setAtrophy(std::string atrophyImageFile)
{
    ScalarReaderType::Pointer   scalarReader = ScalarReaderType::New();
    scalarReader->SetFileName(atrophyImageFile);
    scalarReader->Update();
    setAtrophy(scalarReader->GetOutput());
}

//does not guarantee that this is a valid atrophy!
void AdLem3D::setAtrophy(ScalarImageType::Pointer inputAtrophy)
{
    mAtrophy = inputAtrophy;
    //Make it valid or not (for DIRICHLET bc)! i.e.
    //1. Set all border pixels to 0;
    //2. Modify all second-last border pixels such that sum of all pixels = 0;
//    modifyAtrophy();
}

void AdLem3D::scaleAtorphy(double factor)
{
    typedef itk::MultiplyImageFilter<ScalarImageType, ScalarImageType, ScalarImageType> MultiplyImageFilterType;
    MultiplyImageFilterType::Pointer multiplyImageFilter = MultiplyImageFilterType::New();
    multiplyImageFilter->SetInput(mAtrophy);
    multiplyImageFilter->SetConstant(factor);
    multiplyImageFilter->Update();
    mAtrophy=multiplyImageFilter->GetOutput();
}

void AdLem3D::createAtrophy(unsigned int size[3])
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
                if(mAtrophy->GetPixel(pos) != posAtrophy)
                    mAtrophy->SetPixel(pos,negAtrophy);
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

void AdLem3D::writeAtrophyToFile(std::string fileName) {
    ScalarWriterType::Pointer writer = ScalarWriterType::New();
    writer->SetFileName(fileName);
    writer->SetInput(mAtrophy);
    writer->Update();
}

void AdLem3D::modifyAtrophy() {
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
//    std::cout<<"num of pixels: "<<numOfPixels<<std::endl;
    aSum/=(-1*(double)numOfPixels);
//    std::cout<<"asum = "<<aSum<<std::endl;

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
    if(!isAtrophyValid(0.000009)) {
        std::cout<<"strange, should be valid!!"<<std::endl;
    }
}

void AdLem3D::setBrainMask(std::string maskImageFile)
{

    ScalarReaderType::Pointer   scalarReader;
    scalarReader->SetFileName(maskImageFile);
    scalarReader->Update();
    mBrainMask = scalarReader->GetOutput();
}

void AdLem3D::setLameParameters(double muCsf, double lambdaCsf,
                                double muRatio, double lambdaRatio)
{
    mMuCsf = muCsf;
    mLambdaCsf = lambdaCsf;
    mMuGm = muCsf*muRatio;      mMuWm = muCsf*muRatio;
    mLambdaGm = lambdaCsf*lambdaRatio;      mLambdaWm = lambdaCsf*lambdaRatio;
}

void AdLem3D::writeSolution(std::string resultsPath)
{
    std::string velocityFileName(resultsPath+"vel.mha");
    std::string pressureFileName(resultsPath+"press.mha");
    std::string matSolFileName(resultsPath+"sol");
    std::string matSysFileName(resultsPath+"sys");
    std::string matSizeSysFileName(resultsPath+"size_lin_sys");
    mPetscSolverTaras->writeToMatFile(matSolFileName,false,matSysFileName);
    std::ofstream size_file;
    size_file.open(matSizeSysFileName.c_str());
    size_file<<mDomainRegion.GetSize()[0]+1<<" "<<mDomainRegion.GetSize()[1]+1<<" "
            <<mDomainRegion.GetSize()[2]+1;
    size_file.close();

    itk::Index<3>           outputImageStart;       //should be 0!
    ScalarImageType::RegionType domainRegion;
    for (unsigned int i=0; i<3; ++i) outputImageStart.SetElement(i,0);
    domainRegion.SetSize(mDomainRegion.GetSize());
    domainRegion.SetIndex(outputImageStart);

    ScalarWriterType::Pointer   scalarWriter = ScalarWriterType::New();
    VectorWriterType::Pointer   vectorWriter = VectorWriterType::New();

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

    typedef itk::ImageRegionIterator<VectorImageType> VectorIteratorType;
    typedef itk::ImageRegionIterator<ScalarImageType> ScalarIteratorType;

    VectorIteratorType velocityIterator(mVelocity,mVelocity->GetLargestPossibleRegion());
    ScalarIteratorType pressureIterator(mPressure,mPressure->GetLargestPossibleRegion());
    VectorImageType::PixelType velocityPixel;
    ScalarImageType::PixelType pressurePixel;

    unsigned int pos[3];
    unsigned int k = 0;
    velocityIterator.GoToBegin();
    pressureIterator.GoToBegin();
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
                ++i;
            }
            ++j;
        }
        ++k;
    }

    vectorWriter->SetFileName(velocityFileName);
    vectorWriter->SetInput(mVelocity);
    vectorWriter->Update();

    scalarWriter->SetFileName(pressureFileName);
    scalarWriter->SetInput(mPressure);
    scalarWriter->Update();

}


long double AdLem3D::muAt(int x, int y, int z) const
{
    /*if ( (x > (mXnum/2. - 4)) && (x < (mXnum/2. + 4))
         && (y > (mYnum/2. - 4)) && (y < (mYnum/2. + 4))
         && (z > (mZnum/2. - 4)) && (z < (mZnum/2. + 4))) {
        return mMuGm;
    }*/
    return mMuCsf;
}

long double AdLem3D::lambdaAt(int x, int y, int z) const
{
    /*if ( (x > (mXnum/2. - 4)) && (x < (mXnum/2. + 4))
         && (y > (mYnum/2. - 4)) && (y < (mYnum/2. + 4))
         && (z > (mZnum/2. - 4)) && (z < (mZnum/2. + 4))) {
        return mLambdaGm;
    }
    return mLambdaCsf;*/
    return mLambdaCsf;
}

long double AdLem3D::aAt(int x, int y, int z) const
{
    /*if ( (x > (mXnum/2. - 1)) && (x < (mXnum/2. + 1))
         && (y > (mYnum/2. - 1)) && (y < (mYnum/2. + 1))
         && (z > (mZnum/2. - 1)) && (z < (mZnum/2. + 1))) {
        return -0.9;
    } else {
        if (x == 3 && y == 3 && z == 3)
            return 0.6;
        else {
            if (x==mXnum-4 && y==3 && z==3)
                return 0.3;
            else
                return 0;
        }
    }*/
    ScalarImageType::IndexType pos;
    pos.SetElement(0, mDomainRegion.GetIndex()[0] + x);
    pos.SetElement(1, mDomainRegion.GetIndex()[1] + y);
    pos.SetElement(2, mDomainRegion.GetIndex()[2] + z);

    return(mAtrophy->GetPixel(pos));

}


long double AdLem3D::dataAt(std::string dType, int x, int y, int z)
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

int AdLem3D::getXnum() const
{
    return mDomainRegion.GetSize()[0];
}

int AdLem3D::getYnum() const
{
    return mDomainRegion.GetSize()[1];
}

int AdLem3D::getZnum() const
{
    return mDomainRegion.GetSize()[2];
}

AdLem3D::bcType AdLem3D::getBcType() const { return mBc; }
