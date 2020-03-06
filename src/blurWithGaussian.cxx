#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMaskImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMaskNeighborhoodOperatorImageFilter.h"
#include "itkSobelOperator.h"
#include <itkGaussianOperator.h>
#include "itkCastImageFilter.h"



int main(int argc, char *argv[])
{
    const int IMAGE_DIMENSION = 3;
    typedef itk::Image<float, IMAGE_DIMENSION>  InputImageType;
    typedef itk::Image<unsigned char, InputImageType::ImageDimension>  UnsignedCharImageType;
    typedef itk::Image<float, InputImageType::ImageDimension>  FloatImageType;

    if (argc!= 6) {
        std::cerr<<"usage: blurWithGaussian inputImageFilename maskImageFilename radius(i.e. kernelWidth) variance outputImageFileName"<<std::endl;
        return EXIT_FAILURE;
    }
    std::string inImageFile(argv[1]);
    std::string maskImageFile(argv[2]);
    int         inRadius(atoi(argv[3]));
    float       variance(atof(argv[4]));
    std::string outImageFile(argv[5]);

    InputImageType::Pointer inputImage = InputImageType::New();
    {
        typedef itk::ImageFileReader< InputImageType > InputReaderType;
        InputReaderType::Pointer inputReader = InputReaderType::New();
        inputReader->SetFileName(inImageFile);
        inputReader->Update();
        inputImage = inputReader->GetOutput();
    }

    UnsignedCharImageType::Pointer mask = UnsignedCharImageType::New();

    {
        typedef itk::ImageFileReader< UnsignedCharImageType > MaskReaderType;
        MaskReaderType::Pointer maskReader = MaskReaderType::New();
        maskReader->SetFileName(maskImageFile);
        maskReader->Update();
        mask = maskReader->GetOutput();
    }

    typedef itk::GaussianOperator<float, InputImageType::ImageDimension> GaussianOperatorType;
    itk::Size<InputImageType::ImageDimension> radius;
    radius.Fill(inRadius); // a radius of 1x1 creates a 3x3 operator

    GaussianOperatorType gaussianOperator;
    gaussianOperator.CreateToRadius(radius);
    gaussianOperator.SetVariance(variance);

    // Visualize mask image
    typedef itk::MaskNeighborhoodOperatorImageFilter< InputImageType, UnsignedCharImageType, FloatImageType, float> MaskNeighborhoodOperatorImageFilterType;
    MaskNeighborhoodOperatorImageFilterType::Pointer maskNeighborhoodOperatorImageFilter =
            MaskNeighborhoodOperatorImageFilterType::New();
    maskNeighborhoodOperatorImageFilter->SetMaskImage(mask);

    FloatImageType::Pointer outputImage = inputImage;
    for (unsigned int i = 0; i<InputImageType::ImageDimension; ++i) {
        maskNeighborhoodOperatorImageFilter->SetInput(outputImage);
        gaussianOperator.SetDirection(i);
        gaussianOperator.CreateDirectional();
        maskNeighborhoodOperatorImageFilter->SetOperator(gaussianOperator);
        maskNeighborhoodOperatorImageFilter->Update();
        outputImage = maskNeighborhoodOperatorImageFilter->GetOutput();
    }
    std::cout << "Size: " << gaussianOperator.GetSize() << std::endl;
    std::cout << gaussianOperator << std::endl;
    for(unsigned int i = 0; i < 9; i++)
    {
      std::cout << gaussianOperator.GetOffset(i) << " " << gaussianOperator.GetElement(i) << std::endl;
    }


//    typedef itk::RescaleIntensityImageFilter< FloatImageType, UnsignedCharImageType > RescaleFilterType;
//    RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
//    //  rescaleFilter->SetInput(maskNeighborhoodOperatorImageFilter->GetOutput());
//    rescaleFilter->SetInput(imageTmp);
//    rescaleFilter->Update();

    typedef  itk::ImageFileWriter< InputImageType  > WriterType;
//    typedef  itk::ImageFileWriter< UnsignedCharImageType  > WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outImageFile);
//    writer->SetInput(rescaleFilter->GetOutput());
    writer->SetInput(outputImage);
    writer->Update();

    return EXIT_SUCCESS;
}
