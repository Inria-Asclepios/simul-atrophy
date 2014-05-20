
//--------------------------- Compose all the input displacement fields -----------//
/* Options include :
--numOfFields       number of fields to compose.
--inPrefix          inputFileBase    : Filenames to be composed should be in xxNyy format. Provide xx to this argument.
--inPostfix         inputExt         : yy (including the extension)
--outputImage       outImage         : output filename.
--interpolator      l/n            : linear/nearestNeighbor. Currently no implemented, so uses linear by default.

Author:      Bishesh Khanal
Asclepios, INRIA Sophia Antipolis

FIXME: Experiment different interpolators. Currently default vector interpolator is used.
*/

#include <iostream>
#include <sstream>
#include <string>

#include <itkImage.h>
#include <itkImageAdaptor.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "InverseDisplacementImageFilter.h"
#include <itkComposeDisplacementFieldsImageFilter.h>

//#include <itkVectorInterpolateImageFunction.h>
//#include <itkVectorNearestNeighborInterpolateImageFunction.h>
//#include <itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h>

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    std::string inPrefix;
    std::string inPostfix;
    std::string outFile;
    std::string interpolator;

    int     numOfFields;

    //------------------- Set up the command line options database and parse them -----------------------------//
    boost::program_options::options_description optionsDescription("Possible options");
    optionsDescription.add_options()
            ("help,h",
             "Note: All the displacement fields to be composed should be named in the format xxNyy.\n"
             "N starts from 1 up to the final field to be composed.")
             ("inPrefix", boost::program_options::value< std::string >(&inPrefix),
              "Input filename prefix. E.g. step for the series step1vel.mha, step2vel.mha ...")
             ("inPostfix", boost::program_options::value< std::string >(&inPostfix),
              "Input filename postfix. vel.mha for the series  step1vel.mha, step2vel.mha ...")
             ("outFile", boost::program_options::value< std::string >(&outFile),
              "output displacement field filename.")
             ("interpolator", boost::program_options::value< std::string >(&interpolator)->default_value("linear"),
              "linear/bspline/nearestneighbor interpolator.")
             ("numOfFields",boost::program_options::value< int >(&numOfFields),
              "number of displacement fields to be composed.")
             ;

    boost::program_options::variables_map options;
    boost::program_options::store(boost::program_options::parse_command_line(argc,argv,optionsDescription),options);
    boost::program_options::notify(options);

    //If help is asked!
    if(options.count("help")) {
        std::cout<<optionsDescription<<std::endl;
        return EXIT_SUCCESS;
    }

    //Confirm all the required options are given.
    if(!options.count("inPrefix") || !options.count("inPostfix") || !options.count("numOfFields") || !options.count("outFile")
            ) {
        std::cerr<<"invalid options! run with --help or -h to see the proper options."<<std::endl;
        return EXIT_FAILURE;
    }

    //

    //---------------------  displacement field types ----------------------//
    typedef itk::Vector<double, 3>                      DisplacementPixelType;
    typedef itk::Image<DisplacementPixelType, 3>        DisplacementImageType;
    typedef itk::ImageFileReader<DisplacementImageType> DisplacementReaderType;

    //------- read the first disp. field and set it as a warper field -----------//
    DisplacementImageType::Pointer warperField;
    DisplacementReaderType::Pointer warperReader = DisplacementReaderType::New();
    warperReader->SetFileName(inPrefix + boost::lexical_cast<std::string>(1) + inPostfix);
    warperReader->Update();
    warperField = warperReader->GetOutput();

    //----------- read disp. field starting from the second one, compose it and set the composed one as a new warper. Iterate until the end.-----//
    DisplacementImageType::Pointer displacementField;
    typedef itk::ComposeDisplacementFieldsImageFilter<DisplacementImageType, DisplacementImageType> ComposerType;
    for (int i = 2; i<=numOfFields; ++i) {
        std::string displacementFile = inPrefix + boost::lexical_cast<std::string>(i) + inPostfix;
        DisplacementReaderType::Pointer displacementReader = DisplacementReaderType::New();
        displacementReader->SetFileName(displacementFile);
        displacementReader->Update();
        displacementField = displacementReader->GetOutput();

        ComposerType::Pointer composer = ComposerType::New();
        composer->SetDisplacementField(displacementField);
        composer->SetWarpingField(warperField);

        /*FIXME: experiment with other interpolators. if (interpolator.compare("nearestneighbor")==0) {
            typedef itk::VectorNearestNeighborInterpolateImageFunction<DisplacementImageType, ..> InterpolatorFilterType;
            InterpolatorFilterType::Pointer interpolatorFilter = InterpolatorFilterType::New();
            composer->SetInterpolator(interpolatorFilter);
        }*/

        composer->Update();
        warperField = composer->GetOutput();
    }

    //----------- Write the output -------------//
    typedef itk::ImageFileWriter<DisplacementImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outFile);
    writer->SetInput(warperField);
    writer->Update();

    return(EXIT_SUCCESS);

}

