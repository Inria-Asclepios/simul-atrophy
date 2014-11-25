#include <algorithm>

#include <stdio.h>

//#include "itkCastImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkExtractImageFilter.h"
//#include "itkLabelStatisticsImageFilter.h"
#include "itkVersion.h"

#include <string>
#include <vector>

int main(int argc, char *argv[])
{
    const unsigned int ImageDimension = 3;

    if( argc < 5 || argc > 6 )
    {
        std::cout << "Extract a sub-region from image using the bounding"
                     " box from a label image, with optional padding radius."
                  << std::endl << "Usage : " << argv[0] << " inputImage outputImage labelMaskImage [label=1] [padRadius=0]"
                  << std::endl;
        if( argc >= 2 &&
                ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
        {
            return EXIT_SUCCESS;
        }
        return EXIT_FAILURE;
    }

    // ** Get in the inputs **//
    std::string inputImageFile(argv[1]), outputImageFile(argv[2]), labelMaskImageFile(argv[3]);
    const unsigned int label = (argc >= 5) ? atoi(argv[4]) : 1;
    const unsigned int padWidth = (argc >= 6) ? atoi(argv[5]) : 0;

    // --------------------------------------------------------------------//
    typedef itk::DiffusionTensor3D< float > PixelType;
    //  typedef itk::Image<itk::DiffusionTensor3D<double>, 3>       TensorImageType;
    typedef itk::Image<PixelType, ImageDimension> ImageType;
    typedef itk::ImageFileReader<ImageType>       ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(inputImageFile);
    reader->Update();

    typedef itk::Image<unsigned short, ImageDimension> ShortImageType;
    typedef itk::ImageFileReader<ShortImageType> ShortImageReaderType;
    typename ShortImageReaderType::Pointer shortReader = ShortImageReaderType::New();
    shortReader->SetFileName(labelMaskImageFile);
    shortReader->Update();

    typename ShortImageType::IndexType minIndex, maxIndex, currIndex;
    typename ShortImageType::SizeType  size;
    minIndex.Fill(100000);   maxIndex.Fill(0);
    size.Fill(0);
    itk::ImageRegionConstIteratorWithIndex< ShortImageType > it(shortReader->GetOutput(), shortReader->GetOutput()->GetRequestedRegion());
    it.GoToBegin();
    while( !it.IsAtEnd() ) {
        if(it.Get() == label) {
            currIndex = it.GetIndex();
            for(unsigned int i=0; i<ImageDimension; ++i) {
                if (minIndex[i] > currIndex[i]) minIndex.SetElement(i, currIndex[i]);
                if (maxIndex[i] < currIndex[i]) maxIndex.SetElement(i, currIndex[i]);
            }
        }
        ++it;
    }

    for(unsigned int i = 0; i<ImageDimension; ++i)
        size[i] = maxIndex[i] - minIndex[i] + 1;

    typename ImageType::RegionType region;
    region.SetIndex(minIndex);
    region.SetSize(size);

    std::cout << "bounding box of label=" << label
              << " : " << region << std::endl;

    region.PadByRadius(padWidth);

    std::cout << "padding radius = " << padWidth
              << " : " << region << std::endl;

    region.Crop(reader->GetOutput()->GetBufferedRegion() );

    std::cout << "crop with original image region " << reader->GetOutput()->GetBufferedRegion()
              << " : " << region << std::endl;

    std::cout << "final cropped region: " << region << std::endl;

    typedef itk::ExtractImageFilter<ImageType, ImageType> CropperType;
    typename CropperType::Pointer cropper = CropperType::New();
    cropper->SetInput(reader->GetOutput() );
    cropper->SetExtractionRegion(region);
    cropper->SetDirectionCollapseToSubmatrix();
    cropper->Update();

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput(cropper->GetOutput() );
    writer->SetFileName(outputImageFile);
    writer->Update();

    return EXIT_SUCCESS;
}

