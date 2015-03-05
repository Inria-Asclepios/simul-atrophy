#include <iostream>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkSignedDanielssonDistanceMapImageFilter.h>

int main (int argc, char** argv) {
    if(argc<3) {
        std::cout<<"usage: "<<argv[0]<<" inputBinaryImage outputImage"<<std::endl;
        return EXIT_FAILURE;
    }
    std::string inFileName(argv[1]), outFileName(argv[2]);
    typedef typename itk::Image< int, 3 > ImageType;
    typedef typename itk::ImageFileReader< ImageType > ReaderType;
    typedef typename itk::ImageFileWriter< ImageType > WriterType;

    typedef typename itk::SignedDanielssonDistanceMapImageFilter< ImageType, ImageType > DistanceFilterType;

    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(inFileName);

    DistanceFilterType::Pointer distanceFilter = DistanceFilterType::New();
    distanceFilter->SetInput(reader->GetOutput());

    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outFileName);
    writer->SetInput(distanceFilter->GetDistanceMap());
    writer->Update();

    return(EXIT_FAILURE);

}

