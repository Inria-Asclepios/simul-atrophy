/*
  CreateImageFromRegion outputImage indexBegin size insideValue outsideValue referenceImage
*/
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
    std::string outImgFile;
    std::string refImgFile;
    std::string indexStartString;
    std::string regionSizeString;
    std::string directionString;

    double      insideValue;
    double      outsideValue;

    //------------------- Set up the command line options database and parse them -----------------------------//
    boost::program_options::options_description optionsDescription("binarize the input image by setting one to the pixels having values within the given range");
    optionsDescription.add_options()
            ("help,h", "display help message")
            ("reference,r", boost::program_options::value< std::string >(&refImgFile), "reference image to get the origin, spacing and orientation")
            ("output,o", boost::program_options::value< std::string >(&outImgFile), "output filename")
            ("indexStart,s",boost::program_options::value< std::string >(&indexStartString), "Start index of the region in the format x,y,z")
            ("regionSize,z",boost::program_options::value< std::string >(&regionSizeString), "Size of the region in the format x,y,z")
            ("direction,d",boost::program_options::value< std::string >(&directionString), "direction (-1 or 1) in each of the axes to move in the format x,y,z")
            ("inside,x",boost::program_options::value< double >(&insideValue)->default_value(1), "Inside value")
            ("outside,y",boost::program_options::value< double >(&outsideValue)->default_value(0), "Outside Value")
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
    if(!options.count("reference") || !options.count("output")
            || !options.count("indexStart") || !options.count("regionSize") || !options.count("direction")) {
        std::cerr<<"invalid options! run with --help or -h to see the proper options."<<std::endl;
        return EXIT_FAILURE;
    }


    std::vector< int > direction;
    {
        std::stringstream ss(directionString);
        int size;
        while(ss>>size) {
            direction.push_back(size);
            if(ss.peek() == ',')
                ss.ignore();
        }
    }

    //---------------------  Read the image ----------------------//
    typedef itk::Image<double, 3>              ImageType;
    typedef itk::ImageFileReader<ImageType>     ImageReaderType;

    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(refImgFile);
    reader->Update();
    ImageType::Pointer refImg = reader->GetOutput();

    //Parse the region:
    ImageType::IndexType indexStart;
    {
        std::stringstream ss(indexStartString);
        unsigned int coord;
        unsigned int i=0;
        while(ss>>coord) {
            indexStart[i] = coord;
            if(ss.peek() == ',')
                ss.ignore();
            ++i;
        }
    }

    ImageType::SizeType regionSize;
    {
        std::stringstream ss(regionSizeString);
        int size;
        unsigned int i = 0;
        while(ss>>size) {
            regionSize[i] = size;
            if(ss.peek() == ',')
                ss.ignore();
            ++i;
        }
    }

    std::cout<<"Selected start index: "<<indexStart<<std::endl;
    std::cout<<"Selected size: "<<regionSize<<std::endl;
    std::cout<<"Selected direction: "<<direction.at(0)<<","<<direction.at(1)<<","<<direction.at(2)<<std::endl;

    std::cout<<"Largest possible region from the reference image: "<<std::endl;
    ImageType::RegionType fullRegion = refImg->GetLargestPossibleRegion();
    std::cout<<"Index: "<<fullRegion.GetIndex()<<std::endl;
    std::cout<<"Size: "<<fullRegion.GetSize()<<std::endl;

    ImageType::Pointer outImg;
    outImg = ImageType::New();
    outImg->SetRegions(refImg->GetLargestPossibleRegion());
    outImg->SetOrigin(refImg->GetOrigin());
    outImg->SetSpacing(refImg->GetSpacing());
    outImg->SetDirection(refImg->GetDirection());
    outImg->Allocate();
    //fill all pixels with zero;
    outImg->FillBuffer(outsideValue);

    //set outsideValue to all the pixels lying in the given region.
//    typedef itk::ImageRegionIterator< ImageType >       IteratorType;
//    IteratorType    out_it(outImg,outImg->GetLargestPossibleRegion());
//    IteratorType    in_it(refImg,refImg->GetLargestPossibleRegion());

//    for(out_it.GoToBegin(), in_it.GoToBegin(); !out_it.IsAtEnd(); ++out_it, ++in_it) {
//        {
//            double curr_val = in_it.Get();
//            if(curr_val > lowerThr && curr_val < upperThr )
//                out_it.Set(1);
//        }
//    }



    for(int x = indexStart[0]; x < signed(indexStart[0] + regionSize[0] - 1); (direction.at(0)>0) ? ++x : --x){
        for(int y = indexStart[1]; y < signed(indexStart[1] + regionSize[1] - 1); (direction.at(1)>0) ? ++y : --y){
            for(int z = indexStart[2]; z < signed(indexStart[2] + regionSize[2] - 1); (direction.at(2)>0 ? ++z : --z)){
                ImageType::IndexType pixelIndex;
                if(x<0 || y<0 || z<0) {
                    std::cout<<"Negative index: "<<x<<","<<y<<","<<z<<std::endl;
                    return EXIT_FAILURE;
                }
                pixelIndex[0] = x;
                pixelIndex[1] = y;
                pixelIndex[2] = z;
                if(!fullRegion.IsInside(pixelIndex)) {
                    std::cout<<"Invalid region, this position caused the error: ";
                    std::cout<<pixelIndex<<std::endl;
                    return EXIT_FAILURE;
                }
                outImg->SetPixel(pixelIndex,insideValue);
            }
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


