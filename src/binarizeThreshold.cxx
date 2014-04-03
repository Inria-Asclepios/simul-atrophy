#include <iostream>
#include <string>
#include <sstream>

#include <itkImage.h>
#include <itkImageRegionIterator.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <boost/program_options.hpp>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    std::string inImgFile;
    std::string outImgFile;
    double      lowerThr;
    double      upperThr;
    double      insideVal;
    double      outsideVal;

    //------------------- Set up the command line options database and parse them -----------------------------//
    boost::program_options::options_description optionsDescription("binarize the input image by setting inside value to the pixels having values within the given range");
    optionsDescription.add_options()
            ("help,h", "display help message")
            ("input,i", boost::program_options::value< std::string >(&inImgFile), "input FileName")
            ("output,o", boost::program_options::value< std::string >(&outImgFile), "output filename")
            ("lower,l",boost::program_options::value< double >(&lowerThr), "Lower value")
            ("upper,u",boost::program_options::value< double >(&upperThr), "Upper value")
            ("inside,x",boost::program_options::value< double >(&insideVal)->default_value(1), "Inside pixel value")
            ("outside,y",boost::program_options::value< double >(&outsideVal)->default_value(0), "Outside pixel value")
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
    if(!options.count("input") || !options.count("output")
            || !options.count("lower") || !options.count("upper")) {
        std::cerr<<"invalid options! run with --help or -h to see the proper options."<<std::endl;
        return EXIT_FAILURE;
    }

    //---------------------  Read the image type ----------------------//
    typedef itk::Image<double, 3>              ImageType;
    typedef itk::ImageFileReader<ImageType>     ImageReaderType;

    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(inImgFile);
    reader->Update();
    ImageType::Pointer inImg = reader->GetOutput();

    ImageType::Pointer outImg;

        outImg = ImageType::New();
        outImg->SetRegions(inImg->GetLargestPossibleRegion());
        outImg->SetOrigin(inImg->GetOrigin());
        outImg->SetSpacing(inImg->GetSpacing());
        outImg->SetDirection(inImg->GetDirection());
        outImg->Allocate();
        //fill all pixels with outsideValue;
        outImg->FillBuffer(outsideVal);

    //set 1 all to all the pixels lying in the given range.
    typedef itk::ImageRegionIterator< ImageType >       IteratorType;
    IteratorType    out_it(outImg,outImg->GetLargestPossibleRegion());
    IteratorType    in_it(inImg,inImg->GetLargestPossibleRegion());

    for(out_it.GoToBegin(), in_it.GoToBegin(); !out_it.IsAtEnd(); ++out_it, ++in_it) {
        {
            double curr_val = in_it.Get();
            if(curr_val > lowerThr && curr_val < upperThr )
                out_it.Set(insideVal);
        }
    }

    //------------------------ Write the output Image -------------------------//
    typedef itk::ImageFileWriter< ImageType >       ImageWriterType;
    ImageWriterType::Pointer    writer = ImageWriterType::New();
    writer->SetFileName(outImgFile);
    writer->SetInput(outImg);
    writer->Update();
    return EXIT_SUCCESS;
}

