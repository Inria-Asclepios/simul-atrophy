#include "HelmHoltzDecomposer.hxx"

#include <itkRegionOfInterestImageFilter.h>
#include <itkSubtractImageFilter.h>

#include <itkDerivativeImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkComposeImageFilter.h>

HelmHoltzDecomposer::HelmHoltzDecomposer():mPoissonSolver("HoltzPoisson")
{
    //    mV = VectorImageType::New(); Let's see if this is needed or not, prolly not!
}

HelmHoltzDecomposer::~HelmHoltzDecomposer()
{

}

void HelmHoltzDecomposer::setV(std::string inputFileName, bool fullImage, unsigned int origin[], unsigned int size[])
{
    VectorReaderPointerType vectorReader = VectorReaderType::New();
    vectorReader->SetFileName(inputFileName);
    vectorReader->Update();
    if(fullImage)
        setV(vectorReader->GetOutput(),true,NULL,NULL);
    else
        setV(vectorReader->GetOutput(),false,origin,size);
}

void HelmHoltzDecomposer::setV(VectorImagePointerType v, bool fullImage, unsigned int origin[], unsigned int size[])
{
//    std::cout<<"largest Region of v: "<<v->GetLargestPossibleRegion()<<std::endl;
//    std::cout<<"origin of v: "<<v->GetOrigin()<<std::endl;
    if(fullImage)
        mV = v;
    else {
        ScalarImageType::RegionType domainRegion;
        ScalarImageType::IndexType domainOrigin;
        ScalarImageType::SizeType  domainSize;
        for(int i=0;i<3;++i) {     //ADD error-guard to ensure that
            domainSize.SetElement(i,size[i]);     //size!!
            domainOrigin.SetElement(i,origin[i]); //provided is less or equal to the input
        }
        domainRegion.SetIndex(domainOrigin);
        domainRegion.SetSize(domainSize);

        typedef itk::RegionOfInterestImageFilter< VectorImageType, VectorImageType > RegionOfInterestType;
        RegionOfInterestType::Pointer cropper = RegionOfInterestType::New();
        cropper->SetInput(v);
        cropper->SetRegionOfInterest(domainRegion);
        cropper->Update();
        mV = cropper->GetOutput();
    }
    computeDivV();  //HH decomposer will always need div(v). Let's do it
    //as soon as we have the input V!! This way, if v changes as an input,
    //mDivV will also change!!
}


void HelmHoltzDecomposer::decompose()
{
    mPoissonSolver.setCoeff(1,"");
    mPoissonSolver.setRhs(mDivV);
    mPoissonSolver.setDomainRegionFullImage();
    mPoissonSolver.solveModel();
    mP = mPoissonSolver.getSolution();
    computeGradP();
    //mCurlA = mV - mGradP;
    typedef itk::SubtractImageFilter<VectorImageType, VectorImageType, VectorImageType> subtractorType;
    subtractorType::Pointer subtractor = subtractorType::New();
    subtractor->SetInput1(mV);
    subtractor->SetInput2(mGradP);
    subtractor->Update();
    mCurlA = subtractor->GetOutput();
}

HelmHoltzDecomposer::ScalarImagePointerType HelmHoltzDecomposer::getDivV()
{
    return mDivV;
}

HelmHoltzDecomposer::ScalarImagePointerType HelmHoltzDecomposer::getP()
{
    return mP;
}

//HelmHoltzDecomposer::VectorImagePointerType HelmHoltzDecomposer::getGradP()
HelmHoltzDecomposer::VectorImagePointerType HelmHoltzDecomposer::getGradP()
{
    return mGradP;
}

HelmHoltzDecomposer::VectorImagePointerType HelmHoltzDecomposer::getCurlA()
{
    return mCurlA;
}

void HelmHoltzDecomposer::computeDivV()
{
    //----------------Extract vector field components-------------------
    typedef typename itk::VectorIndexSelectionCastImageFilter< VectorImageType,ScalarImageType> IndexSelectFilterType;
    IndexSelectFilterType::Pointer componentExtractor0 = IndexSelectFilterType::New();
    componentExtractor0->SetIndex(0);
    IndexSelectFilterType::Pointer componentExtractor1 = IndexSelectFilterType::New();
    componentExtractor1->SetIndex(1);
    IndexSelectFilterType::Pointer componentExtractor2 = IndexSelectFilterType::New();
    componentExtractor2->SetIndex(2);

    componentExtractor0->SetInput(mV);
    componentExtractor1->SetInput(mV);
    componentExtractor2->SetInput(mV);

    //-----------------Compute partial derivatives----------------------
    typedef typename itk::DerivativeImageFilter<ScalarImageType, ScalarImageType > DerivativeFilterType;
    DerivativeFilterType::Pointer derivative0=DerivativeFilterType::New();
    DerivativeFilterType::Pointer derivative1=DerivativeFilterType::New();
    DerivativeFilterType::Pointer derivative2=DerivativeFilterType::New();

    derivative0->SetOrder(1);
    derivative1->SetOrder(1);
    derivative2->SetOrder(1);
    derivative0->SetDirection(0);
    derivative1->SetDirection(1);
    derivative2->SetDirection(2);

    derivative0->SetInput(componentExtractor0->GetOutput());
    derivative1->SetInput(componentExtractor1->GetOutput());
    derivative2->SetInput(componentExtractor2->GetOutput());

    //------------------ Sum partial derivatives------------------------
    typedef itk::AddImageFilter<ScalarImageType, ScalarImageType, ScalarImageType> AddFilterType;
    AddFilterType::Pointer addFilter = AddFilterType::New();
    addFilter->SetInput1( derivative0->GetOutput() );
    addFilter->SetInput2( derivative1->GetOutput() );
    addFilter->Update();

    AddFilterType::Pointer addFilter1 = AddFilterType::New();
    addFilter1->SetInput1( derivative2->GetOutput() );
    addFilter1->SetInput2(addFilter->GetOutput());
    addFilter1->Update();

    mDivV = addFilter1->GetOutput();

}

void HelmHoltzDecomposer::computeGradP()
{
    //-----------------Compute partial derivatives----------------------
    typedef typename itk::DerivativeImageFilter<ScalarImageType, ScalarImageType > DerivativeFilterType;
    DerivativeFilterType::Pointer derivative0=DerivativeFilterType::New();
    DerivativeFilterType::Pointer derivative1=DerivativeFilterType::New();
    DerivativeFilterType::Pointer derivative2=DerivativeFilterType::New();

    derivative0->SetOrder(1);
    derivative1->SetOrder(1);
    derivative2->SetOrder(1);
    derivative0->SetDirection(0);
    derivative1->SetDirection(1);
    derivative2->SetDirection(2);

    derivative0->SetInput(mP);
    derivative1->SetInput(mP);
    derivative2->SetInput(mP);

    typedef typename itk::ComposeImageFilter<ScalarImageType,VectorImageType> composeFilterType;
    composeFilterType::Pointer composeFilter = composeFilterType::New();
    composeFilter->SetInput(0,derivative0->GetOutput());
    composeFilter->SetInput(1,derivative1->GetOutput());
    composeFilter->SetInput(2,derivative2->GetOutput());
    composeFilter->Update();
    mGradP = composeFilter->GetOutput();
}

