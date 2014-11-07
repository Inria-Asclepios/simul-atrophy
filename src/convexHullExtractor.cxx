#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkExtractImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include <string>
#include <boost/program_options.hpp>

int main( int argc, char **argv )
{
    std::string inImgFile;
    std::string refImgFile;
    std::string outImgFile;
    int         refImgLabel;

    //------------------- Set up the command line options database and parse them -----------------------------//
    boost::program_options::options_description optionsDescription("Possible options");
    optionsDescription.add_options()
            ("help,h", "display help message")
            ("inImage", boost::program_options::value< std::string >(&inImgFile),
             "Input image file which is to be cropped")
            ("outImage", boost::program_options::value< std::string >(&outImgFile),
             "Output cropped image filename")
            ("refImage", boost::program_options::value< std::string >(&refImgFile),
             "Integer image file of which the convex hull will be extracted to get the cropping region.")
            ("refImageLabel", boost::program_options::value< int >(&refImgLabel),
             "Label of the Integer image file whose region will be used to get the convex hull for obtaining cropping region.")
            ;

    boost::program_options::variables_map options;
    boost::program_options::store(boost::program_options::parse_command_line(argc,argv,optionsDescription),options);
    boost::program_options::notify(options);

    //If help is asked!
    if(options.count("help")) {
        std::cout<<optionsDescription<<std::endl;
        return EXIT_SUCCESS;
    }

    //Confirm all compulsory options are given.
    if(!options.count("inImage") || !options.count("outImage")
            ) {
        std::cerr<<"invalid options! run with --help or -h to see the proper options."<<std::endl;
        return EXIT_FAILURE;
    }

    if(!options.count("refImage")) {
        refImgFile = inImgFile;
        std::cout<<"WARNING: using input image file as reference image too since reference image not given."<<std::endl;
    }

    if (!options.count("refImageLabel")) {
        refImgLabel = 1;
        std::cout<<"WARNING: refImageLabel taken as 1 since it was not provided."<<std::endl;
    }

    typedef float PixelType;
    typedef itk::Image<PixelType, 3> ImageType;
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( refImgFile );
    reader->Update();
    typename ImageType::RegionType region;
//    typename ImageType::RegionType::SizeType size;
//    typename ImageType::RegionType::IndexType index;
    {
        typedef itk::Image<unsigned short, 3> ShortImageType;
        typedef itk::CastImageFilter<ImageType, ShortImageType> CasterType;
        typename CasterType::Pointer caster = CasterType::New();
        caster->SetInput( reader->GetOutput() );
        caster->Update();
        typedef itk::LabelStatisticsImageFilter<ShortImageType, ShortImageType>
                StatsFilterType;
        typename StatsFilterType::Pointer stats = StatsFilterType::New();
        stats->SetLabelInput( caster->GetOutput() );
//        stats->SetLabelInput( reader->GetOutput() );
        stats->SetInput( caster->GetOutput() );
        stats->Update();
        region = stats->GetRegion( refImgLabel );
    }
    std::cout << region << std::endl;
    typedef itk::ExtractImageFilter<ImageType, ImageType> CropperType;
    typename CropperType::Pointer cropper = CropperType::New();

    typename ReaderType::Pointer reader1 = ReaderType::New();
    reader1->SetFileName(inImgFile);
    reader1->Update();
    cropper->SetInput( reader1->GetOutput() );
    cropper->SetExtractionRegion( region );
    cropper->SetDirectionCollapseToSubmatrix();
    cropper->Update();
    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( cropper->GetOutput() );
    writer->SetFileName( outImgFile );
    writer->Update();
    return EXIT_SUCCESS;
}

