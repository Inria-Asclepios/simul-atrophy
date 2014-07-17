/*

createDiffusionTensorImage.cxx
Bishesh Khanal
Asclepios, INRIA Sophia Antipolis

A test program to create diffusion tensor images of desired size.

*/

#include <iostream>
#include <string>
#include <sstream>
#include <boost/program_options.hpp>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkDiffusionTensor3D.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    std::string inRefFileName, outFileName;         //input reference file for origin and spacings; name of the output file to be created.
    std::string elements;       //elements of the tensor separted by comma, without spaces in lower triangular order.
    double D[6];            //six diffusion tensor elements, in LT order.

    //-------------- Set up the command line options-------------------------
    boost::program_options::options_description optionsDescription("Possible options");
    optionsDescription.add_options()
            ("help,h", "displays help message")
            ("reference,r", boost::program_options::value< std::string >(&inRefFileName), "input reference file")
            ("output,o", boost::program_options::value< std::string >(&outFileName), "output filename")
            ("elements,e",boost::program_options::value< std::string >(&elements), "elements separted by comma without any spaces in lower traingular order")
            ;

    boost::program_options::variables_map options;
    boost::program_options::store(boost::program_options::parse_command_line(argc,argv,optionsDescription),options);
    boost::program_options::notify(options);

    //If help is asked!
    if(options.count("help")) {
        std::cout<<optionsDescription<<std::endl;
        return EXIT_SUCCESS;
    }

    //Confirm all the options required are given.
    if(!options.count("reference") || !options.count("output") || !options.count("elements")) {
        std::cerr<<"invalid options! run with --help or -h to see the proper options."<<std::endl;
        return EXIT_FAILURE;
    }

    //Parse the elements into number from the string:
    {
        std::stringstream ss(elements);
        std::cout<<"what is in elements: "<<ss.str()<<std::endl;
        unsigned int i = 0;
        while(ss>>D[i++]) { //get the string until comma to first element of D and then increase i.
            if(ss.peek() == ',') //Still don't know why ss>> gives output upto comma, did not see
                ss.ignore();        //it in the documentation!
        }
    }

    //Check if the elements are properly parsed into D array.
    //for (int i=0;i<6;++i)
    //    std::cout<<std::endl<<D[i];

    //---------------------  Read the reference image type ----------------------//
    typedef itk::Image<double, 3>                       ScalarImageType;
    typedef itk::ImageFileReader<ScalarImageType>       ScalarImageReaderType;

    ScalarImageReaderType::Pointer reader = ScalarImageReaderType::New();
    reader->SetFileName(inRefFileName);
    reader->Update();
    ScalarImageType::Pointer refImg = reader->GetOutput();

    //-------------------- Tensor Image -------------------------------------//
    typedef itk::Image<itk::DiffusionTensor3D<double>, 3>       TensorImageType;
    TensorImageType::Pointer tensorImage = TensorImageType::New();

    //------------ Get important details from the reference image ------------//
    tensorImage->SetRegions(refImg->GetLargestPossibleRegion());
    tensorImage->SetOrigin(refImg->GetOrigin());
    tensorImage->SetSpacing(refImg->GetSpacing());
    tensorImage->SetDirection(refImg->GetDirection());

    // ----------- Allocate memory ---------------------------//
    tensorImage->Allocate();

    // Fill the voxels with the input tensor value
    tensorImage->FillBuffer(D);

    //------------ Write output image -------------------------//
    typedef itk::ImageFileWriter<TensorImageType>       TensorImageWriterType;
    TensorImageWriterType::Pointer writer = TensorImageWriterType::New();
    writer->SetFileName(outFileName);
    writer->SetInput(tensorImage);
    writer->Update();

    return EXIT_SUCCESS;
}

