#ifndef HELMHOLTZDECOMPOSER_HXX
#define HELMHOLTZDECOMPOSER_HXX

/*HH decomposition:             v = grad(p) + curl(A);
--------------------------------------------------------
Bishesh Khanal, Asclepios, INRIA Sophia Antipolis.
--------------------------------------------------------

Input: a 3D vector field: v
---------------------------------------------------------
Output:
grad(p): The irrotational component.
p:
curl(A): The divergence free component.
*/

#include "PoissonProblem.hxx"

#include <iostream>
#include <itkImage.hxx>
#include <itkImageFileReader.hxx>
#include <itkImageFileWriter.hxx>

class HelmHoltzDecomposer {
public:
    //----------******** Typedefinitions ******* ----------------//
    typedef typename itk::Image<double,3>                   ScalarImageType;
    typedef typename ScalarImageType::Pointer               ScalarImagePointerType;
    /*typedef typename ScalarImageType::PixelType  ScalarImagePixelType;
    typedef typename ScalarImageType::IndexType  ScalarImageIndexType;
    typedef typename ScalarImageType::SizeType   ScalarImageSizeType;
    typedef typename ScalarImageType::RegionType ScalarImageRegionType;*/

    typedef typename itk::ImageFileReader<ScalarImageType>  ScalarReaderType;
    typedef typename ScalarReaderType::Pointer              ScalarReaderPointerType;
    typedef typename itk::ImageFileWriter<ScalarImageType>  ScalarWriterType;
    typedef typename ScalarWriterType::Pointer              ScalarWriterPointerType;

    typedef typename itk::Vector<double,3>                  VectorType;
    typedef typename itk::Image<VectorType,3>               VectorImageType;
    typedef typename VectorImageType::Pointer               VectorImagePointerType;

    typedef typename itk::ImageFileReader<VectorImageType>  VectorReaderType;
    typedef typename VectorReaderType::Pointer              VectorReaderPointerType;
    typedef typename itk::ImageFileWriter<VectorImageType>  VectorWriterType;
    typedef typename VectorWriterType::Pointer              VectorWriterPointerType;


    //-------------------- User visible Methods -------------------------//
    HelmHoltzDecomposer();
    ~HelmHoltzDecomposer();

    void setV(std::string inputFileName, bool fullImage, unsigned int origin[3], unsigned int size[3]);
    void setV(VectorImagePointerType v, bool fullImage, unsigned int origin[3], unsigned int size[3]);

    void decompose();

    ScalarImagePointerType getDivV();
    ScalarImagePointerType getP();
    VectorImagePointerType getGradP();
    VectorImagePointerType getCurlA();

private:
    VectorImagePointerType      mV;
    ScalarImagePointerType      mDivV;
    ScalarImagePointerType      mP;
    VectorImagePointerType      mGradP;
    VectorImagePointerType      mCurlA;

    PoissonProblem<double>      mPoissonSolver;

    void computeDivV();
    void computeGradP();
};

#endif // HELMHOLTZDECOMPOSER_HXX
