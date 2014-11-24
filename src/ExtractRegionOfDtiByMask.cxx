#include <algorithm>

#include <stdio.h>

//#include "itkCastImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkExtractImageFilter.h"
//#include "itkLabelStatisticsImageFilter.h"

#include <string>
#include <vector>

int main(int argc, char *argv[])
{
    const unsigned int ImageDimension = 3;

    if( argc < 6 || argc > 7 )
    {
        std::cout << "Extract a sub-region from image using the bounding"
                     " box from a label image, with optional padding radius."
                  << std::endl << "Usage : " << argv[0] << " ImageDimension "
                  << "inputImage outputImage labelMaskImage [label=1] [padRadius=0]"
                  << std::endl;
        if( argc >= 2 &&
                ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
        {
            return EXIT_SUCCESS;
        }
        return EXIT_FAILURE;
    }
    if( argc < 6 || argc > 7 )
    {
        std::cout << "Extract a sub-region from DTI using the bounding"
                     " box from a label image, with optional padding radius."
                  << std::endl << "Usage : " << argv[0] << " ImageDimension "
                  << "inputDiffusionTensorImage outputDiffusionTensorImage labelMaskImage [label=1] [padRadius=0]"
                  << std::endl;
        if( argc >= 2 &&
                ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
        {
            return EXIT_SUCCESS;
        }
        return EXIT_FAILURE;
    }


    typedef itk::DiffusionTensor3D< float > PixelType;
    //  typedef itk::Image<itk::DiffusionTensor3D<double>, 3>       TensorImageType;
    typedef itk::Image<PixelType, ImageDimension> ImageType;
    typedef itk::ImageFileReader<ImageType>       ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(argv[2]);
    reader->Update();

//    typename ImageType::RegionType region;
    //    typename ImageType::RegionType::SizeType size;
    //    typename ImageType::RegionType::IndexType index;

    typedef itk::Image<unsigned short, ImageDimension> ShortImageType;
    typedef itk::ImageFileReader<ShortImageType> ShortImageReaderType;
    typename ShortImageReaderType::Pointer shortReader = ShortImageReaderType::New();
    shortReader->SetFileName(argv[4]);
    shortReader->Update();

    const unsigned int label = (argc >= 6) ? atoi(argv[5]) : 1;
    typename ShortImageType::IndexType minIndex, maxIndex, currIndex;
    typename ShortImageType::SizeType  size;
    minIndex.Fill(100000);   maxIndex.Fill(0);
    size.Fill(0);
    itk::ImageRegionConstIteratorWithIndex< ShortImageType > it(shortReader->GetOutput(), shortReader->GetOutput()->GetRequestedRegion());
    it.Begin();
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

    for(int i = 0; i<ImageDimension; ++i)
        size[i] = maxIndex[i] - minIndex[i] + 1;

    //    region = stats->GetRegion(label);
    typename ImageType::RegionType region;
    region.SetIndex(minIndex);
    region.SetSize(size);

    std::cout << "bounding box of label=" << label
              << " : " << region << std::endl;

    const unsigned int padWidth = (argc >= 7) ? atoi(argv[6]) : 0;

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
    writer->SetFileName(argv[3]);
    writer->Update();

    return EXIT_SUCCESS;
}

