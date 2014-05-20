//--------------------------- Warp the input scalar 3D image with the given input displacement field-----------//
/* Options include :
--inputDisplacement inDisplacement   :
--inputImage        inImage          :
--outputImage       outImage
--scale             scaleValue       : Value with which the displacement field will be scaled before using it to warp the image. Default: 1
--invert            true/false       : Whether the field must be inverted or not.
--modulate          true/false       : Whether it should be modulated with the jacobian determinant or not.
--interpolator      l/b/n            : linear/bspline/nearestNeighbor
--order             bspline order    : if bspline interpolator selected

Author:      Bishesh Khanal
Asclepios, INRIA Sophia Antipolis
*/

#include <iostream>
#include <sstream>
#include <string>

#include <itkImage.h>
#include <itkImageAdaptor.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "InverseDisplacementImageFilter.h"
#include <itkWarpImageFilter.h>

#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkDivideImageFilter.h>

#include <itkDisplacementFieldJacobianDeterminantFilter.h>

#include <boost/program_options.hpp>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    std::string inImgFile;
    std::string outImgFile;
    std::string displacementFile;
    std::string interpolator;

    double scale;
    bool invertField;
    bool modulate;
    int bsplineOrder;

    //------------------- Set up the command line options database and parse them -----------------------------//
    boost::program_options::options_description optionsDescription("Possible options");
    optionsDescription.add_options()
            ("help,h", "display help message")
            ("inImage", boost::program_options::value< std::string >(&inImgFile),
             "Filename of the input image to be warped")
            ("outImage", boost::program_options::value< std::string >(&outImgFile),
             "Filename of the output warped image")
            ("displacementImage", boost::program_options::value< std::string >(&displacementFile),
             "Filename of the displacement field image")
            ("interpolator", boost::program_options::value< std::string >(&interpolator)->default_value("linear"),
             "linear/bspline/nearestneighbor interpolator.")
            ("scaleValue",boost::program_options::value< double >(&scale)->default_value(1.),
             "scale the displacement field by s.")
            ("modulate",boost::program_options::value< bool >(&modulate)->default_value(false),
             "true or false: modulate with Jacobian determinant of the field.")
            ("invert",boost::program_options::value< bool >(&invertField)->default_value(false),
             "true/false invert displacement field before warping.")
            ("order",boost::program_options::value< int >(&bsplineOrder),
             "Bspline order, use only if bspline interpolator used")
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
    if(!options.count("inImage") || !options.count("outImage") || !options.count("displacementImage")
            ) {
        std::cerr<<"invalid options! run with --help or -h to see the proper options."<<std::endl;
        return EXIT_FAILURE;
    }

    //---------------------  Read the displacement field ----------------------//
    typedef itk::Vector<double, 3>                      DisplacementPixelType;
    typedef itk::Image<DisplacementPixelType, 3>        DisplacementImageType;
    typedef itk::ImageFileReader<DisplacementImageType> DisplacementReaderType;
    DisplacementImageType::Pointer warperField;

    DisplacementReaderType::Pointer displacementReader = DisplacementReaderType::New();
    displacementReader->SetFileName(displacementFile);
    displacementReader->Update();
    warperField = displacementReader->GetOutput();

    //---------------------- Scale the displacement field --------------------//
    double eps = 0.001;
    if( !((scale>(1-eps)) && (scale < (1+eps))) ){
        itk::ImageRegionIterator<DisplacementImageType> iterator(warperField,warperField->GetLargestPossibleRegion());
        iterator.GoToBegin();
        DisplacementImageType::PixelType pixelValue;
        while(!iterator.IsAtEnd()) {
            pixelValue = iterator.Get();
            pixelValue[0] = pixelValue[0]*scale;
            pixelValue[1] = pixelValue[1]*scale;
            pixelValue[2] = pixelValue[2]*scale;
            iterator.Set(pixelValue);
            ++iterator;
        }
    }

    //----------------------- Invert the displacement field ----------------------//
    if(invertField) {
        typedef InverseDisplacementImageFilter<DisplacementImageType> FPInverseType;
        FPInverseType::Pointer inverter1 = FPInverseType::New();
        inverter1->SetInput(warperField);
        inverter1->SetErrorTolerance(1e-1);
        inverter1->SetMaximumNumberOfIterations(50);
        inverter1->Update();
        std::cout<<"tolerance not reached for "<<inverter1->GetNumberOfErrorToleranceFailures()<<" pixels"<<std::endl;
        warperField = inverter1->GetOutput();

        typedef itk::ImageFileWriter<DisplacementImageType> WriterType;
        WriterType::Pointer writer = WriterType::New();
        writer->SetFileName("displacementFieldInverted.nii.gz");
        writer->SetInput(warperField);
        writer->Update();
    }

    //--------------------  Read and warp the input image -------------------------//
    typedef double                                  ImagePixelType;
    typedef itk::Image<ImagePixelType, 3>           ImageType;
    typedef itk::ImageFileReader<ImageType>         ImageReaderType;
    ImageType::Pointer   inputImage;
    ImageType::Pointer  outputImage;
    {
        ImageReaderType::Pointer imageReader = ImageReaderType::New();
        imageReader->SetFileName(inImgFile);
        imageReader->Update();
        inputImage = imageReader->GetOutput();

        typedef itk::WarpImageFilter<ImageType,ImageType,DisplacementImageType> WarpFilterType;
        WarpFilterType::Pointer warper = WarpFilterType::New();
        warper->SetDisplacementField(warperField);
        warper->SetInput(inputImage);
        warper->SetOutputSpacing(inputImage->GetSpacing());
        warper->SetOutputOrigin(inputImage->GetOrigin());
        warper->SetOutputDirection(inputImage->GetDirection());

        if(interpolator.compare("bspline")==0) {
            typedef itk::BSplineInterpolateImageFunction<ImageType> InterpolatorFilterType;
            InterpolatorFilterType::Pointer interpolatorFilter = InterpolatorFilterType::New();
            interpolatorFilter->SetSplineOrder(bsplineOrder);
            warper->SetInterpolator(interpolatorFilter);
            //            warper->Update();
        } else if (interpolator.compare("nearestneighbor")==0) {
            typedef itk::NearestNeighborInterpolateImageFunction<ImageType> InterpolatorFilterType;
            InterpolatorFilterType::Pointer interpolatorFilter = InterpolatorFilterType::New();
            warper->SetInterpolator(interpolatorFilter);
            //            warper->Update();
        } else {
            //            warper->Update();
        }
        warper->Update();
        outputImage = warper->GetOutput();

        //-------------------- Jacobian determinant of the input field ----------------------------------//
        if(modulate) {
            typedef itk::DisplacementFieldJacobianDeterminantFilter< DisplacementImageType, double > JacobianFilterType;
            JacobianFilterType::Pointer jacobianFilter = JacobianFilterType::New();
            jacobianFilter->SetUseImageSpacingOff();
            jacobianFilter->SetInput(warperField);
            jacobianFilter->Update();

            typedef itk::ImageFileWriter<ImageType> WriterType;
            WriterType::Pointer writer1 = WriterType::New();
            writer1->SetFileName("jacobian.mha");
            writer1->SetInput(jacobianFilter->GetOutput());
            writer1->Update();

            //----------------------- modulate the warped atrophy map with the jacobian --------------------//
            typedef itk::DivideImageFilter< ImageType, ImageType, ImageType > ModulatorType;
            ModulatorType::Pointer modulator = ModulatorType::New();
            modulator->SetInput1(outputImage);
            modulator->SetInput2(jacobianFilter->GetOutput());
            modulator->Update();
            outputImage = modulator->GetOutput();
        }

        //---------------------- write the warped with modulation atrophy map -----------------------//
        typedef itk::ImageFileWriter<ImageType> WriterType;
        WriterType::Pointer writer = WriterType::New();
        writer->SetFileName(outImgFile);
        writer->SetInput(outputImage);
        writer->Update();
    }

    //use it to warp atrophy.
    return(EXIT_SUCCESS);

}

