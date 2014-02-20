#include <iostream>
#include <string>
#include <sstream>

#include <itkImage.h>
#include <itkImageRegionIterator.h>
//#include <itkImageAdaptor.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

//#include <boost/program_options/parsers.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    std::string inImgFile;
    std::string outImgFile;
    bool        modifyImage;
    std::string sphereCenterString;
    double      sphereRadius;
    double      atrophyValue;

    //------------------- Set up the command line options database and parse them -----------------------------//
    boost::program_options::options_description optionsDescription("Possible options");
    optionsDescription.add_options()
            ("help,h", "display help message")
            ("input,i", boost::program_options::value< std::string >(&inImgFile), "input FileName")
            ("output,o", boost::program_options::value< std::string >(&outImgFile), "output filename")
            ("modify,m",boost::program_options::value< bool >(&modifyImage))
            ("center,c",boost::program_options::value< std::string >(&sphereCenterString), "Center of the sphere in the format x,y,z")
            ("radius,r",boost::program_options::value< double >(&sphereRadius), "Radius of the sphere")
            ("value,v",boost::program_options::value< double >(&atrophyValue), "Value to be set in the selected spherical region")
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
    if(!options.count("modify") || !options.count("input") || !options.count("output")
            || !options.count("center") || !options.count("radius") || !options.count("value")) {
        std::cerr<<"invalid options! run with --help or -h to see the proper options."<<std::endl;
        return EXIT_FAILURE;
    }

    //Parse the sphereCenter into number from teh string:
    std::vector< double > sphereCenter;
    {
        std::stringstream ss(sphereCenterString);
        double coord;
        while(ss>>coord) {
            sphereCenter.push_back(coord);
            if(ss.peek() == ',')
                ss.ignore();
        }
    }

    //---------------------  Read the image type ----------------------//
    typedef itk::Image<double, 3>              ImageType;
    typedef itk::ImageFileReader<ImageType>     ImageReaderType;

    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(inImgFile);
    reader->Update();
    ImageType::Pointer inImg = reader->GetOutput();

    ImageType::Pointer outImg;

    if(modifyImage) {
        outImg = inImg;
    } else {
        outImg = ImageType::New();
        outImg->SetRegions(inImg->GetLargestPossibleRegion());
        outImg->SetOrigin(inImg->GetOrigin());
        outImg->SetSpacing(inImg->GetSpacing());
        outImg->SetDirection(inImg->GetDirection());
        outImg->Allocate();

        //fill all pixels with zero;
        outImg->FillBuffer(0);

    }

    //set atrophyValue to all pixels within the sphere of the given input radius  and input center point.
    typedef itk::ImageRegionIterator< ImageType >       IteratorType;
    IteratorType    iterator(outImg,outImg->GetLargestPossibleRegion());
    for(iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator) {
        {
            ImageType::PointType currPoint;
            ImageType::IndexType currIndex;
            currIndex = iterator.GetIndex();
            outImg->TransformIndexToPhysicalPoint(currIndex,currPoint);
            if( (
                        (currPoint[0]-sphereCenter.at(0))*(currPoint[0]-sphereCenter.at(0)) +
                        (currPoint[1]-sphereCenter.at(1))*(currPoint[1]-sphereCenter.at(1)) +
                        (currPoint[2]-sphereCenter.at(2))*(currPoint[2]-sphereCenter.at(2))
                        ) <= (sphereRadius * sphereRadius)
                    )
                iterator.Set(atrophyValue);

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
