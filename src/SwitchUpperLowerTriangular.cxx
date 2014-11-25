/*

SwitchFormatBetweenUpperAndLowerTriangular.cxx
Bishesh Khanal
Asclepios, INRIA Sophia Antipolis

ITK, FSL uses Upper Triangular format: [xx xy xz yy yz zz]
TTK, MedInria uses LowerTriangular format: [xx yx yy zx zy zz]

Switching between these two is done simply by swapping the 3rd and 4th element
of the vector in which the values are stored for all pixels.

*/

#include <iostream>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkDiffusionTensor3D.h>
#include <itkImageRegionIterator.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    //---------------------  Read input image ----------------------//
    if (argc < 3) {
        std::cout<<"usage: "<<argv[0]<<" inputDtiFile outputDtiFile"<<std::endl;
        return EXIT_FAILURE;
    }
    std::string inFile, outFile;

    typedef itk::Image<itk::DiffusionTensor3D<double>, 3> ImageType;
    typedef itk::ImageFileReader<ImageType>       ImageReaderType;

    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(inFile);
    reader->Update();
    ImageType::Pointer tensorImage = reader->GetOutput();

    // Test how the values are stored:
    //    std::cout<<"origin: "<<tensorImage->GetOrigin()<<std::endl;
    //    std::cout<<"direction: "<<tensorImage->GetDirection()<<std::endl;


    itk::ImageRegionIterator< ImageType > it(tensorImage, tensorImage->GetLargestPossibleRegion());
    it.GoToBegin();
    typename ImageType::PixelType currPixel;
    while (!it.IsAtEnd()) {
        currPixel = it.Get();
        double tmp;
        tmp = currPixel[2];
        currPixel[2] = currPixel[3];
        currPixel[3] = tmp;
        it.Set(currPixel);
        ++it;
    }

    //------------ Write output image -------------------------//
    typedef itk::ImageFileWriter<ImageType>       ImageWriterType;
    ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetFileName(outFile);
    writer->SetInput(tensorImage);
    writer->Update();

    return EXIT_SUCCESS;
}


