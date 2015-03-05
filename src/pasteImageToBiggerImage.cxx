#include <algorithm>

#include <stdio.h>

//#include "itkCastImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkExtractImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkVersion.h"

#include <string>
#include <vector>

int main(int argc, char *argv[])
{
    const unsigned int ImageDimension = 3;

    if( argc < 5 )
    {
        std::cout << "Paste smalle image to a bigger image."
                     " The paste region is a bounding box computed from the given labelMask image."
                  << std::endl << "Usage : " << argv[0] << " inputSmallImage inputBigImage labelMaskImage outputImage "
                  << std::endl;
        if( argc >= 2 &&
                ( std::string( argv[1] ) == std::string("--help") || std::string( argv[1] ) == std::string("-h") ) )
        {
            return EXIT_SUCCESS;
        }
        return EXIT_FAILURE;
    }

    // ** Get in the inputs **//
    std::string inputSmallImageFile(argv[1]), inputBigImageFile(argv[2]);
    std::string labelMaskImageFile(argv[3]), outputImageFile(argv[4]);
    const unsigned int label = 1;

    // --------------------------------------------------------------------//
    typedef float PixelType;
    typedef itk::Image<PixelType, ImageDimension> ImageType;
    typedef itk::ImageFileReader<ImageType>       ReaderType;
    ReaderType::Pointer readerSmallImage = ReaderType::New();
    readerSmallImage->SetFileName(inputSmallImageFile);
    readerSmallImage->Update();

    ReaderType::Pointer readerBigImage = ReaderType::New();
    readerBigImage->SetFileName(inputBigImageFile);
    readerBigImage->Update();

    typedef itk::Image<unsigned short, ImageDimension> ShortImageType;
    typedef itk::ImageFileReader<ShortImageType> ShortImageReaderType;
    ShortImageReaderType::Pointer shortReader = ShortImageReaderType::New();
    shortReader->SetFileName(labelMaskImageFile);
    shortReader->Update();

    ShortImageType::IndexType minIndex, maxIndex, currIndex;
    //ShortImageType::SizeType  size;
    minIndex.Fill(100000);   maxIndex.Fill(0);
    //size.Fill(0);
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

    // for(unsigned int i = 0; i<ImageDimension; ++i)
    //     size[i] = maxIndex[i] - minIndex[i] + 1;

    // ImageType::RegionType region;
    // region.SetIndex(minIndex);
    // region.SetSize(size);

    // std::cout << "bounding box of label=" << label
    //           << " : " << region << std::endl;

    typedef itk::PasteImageFilter <ImageType, ImageType >
    PasteImageFilterType;
  // The SetDestinationIndex() method prescribes where in the first
  // input to start pasting data from the second input.
  // The SetSourceRegion method prescribes the section of the second
  // image to paste into the first.
    ImageType::IndexType destinationIndex;
    destinationIndex[0] = 10;
    destinationIndex[1] = 10;
    PasteImageFilterType::Pointer pasteFilter = PasteImageFilterType::New ();
    pasteFilter->SetSourceImage(readerSmallImage->GetOutput());
    pasteFilter->SetDestinationImage(readerBigImage->GetOutput());
    pasteFilter->SetSourceRegion(readerSmallImage->GetOutput()->GetLargestPossibleRegion());
    pasteFilter->SetDestinationIndex(minIndex);


    // region.Crop(reader->GetOutput()->GetBufferedRegion() );

    // std::cout << "crop with original image region " << reader->GetOutput()->GetBufferedRegion()
    //           << " : " << region << std::endl;

    // std::cout << "final cropped region: " << region << std::endl;

    // typedef itk::ExtractImageFilter<ImageType, ImageType> CropperType;
    // typename CropperType::Pointer cropper = CropperType::New();
    // cropper->SetInput(reader->GetOutput() );
    // cropper->SetExtractionRegion(region);
    // cropper->SetDirectionCollapseToSubmatrix();
    // cropper->Update();

    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput(pasteFilter->GetOutput() );
    writer->SetFileName(outputImageFile);
    writer->Update();

    return EXIT_SUCCESS;
}

