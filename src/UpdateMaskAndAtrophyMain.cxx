//--------------------------- Update the brain mask and atrophy map with the given displacement field-----------//
/* This will warp the brain mask but for the atrophy map warping is modulated by the jacobian determinant of the displacement field
Bishesh Khanal
Asclepios, INRIA Sophia Antipolis
*/

#include <iostream>
#include <sstream>
#include <string>

#include <itkImage.h>
#include <itkImageAdaptor.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkInverseDisplacementFieldImageFilter.h>
#include "InverseDisplacementImageFilter.h"


#include <itkWarpImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkMultiplyImageFilter.h>
#include <itkDivideImageFilter.h>

#include <itkDisplacementFieldJacobianDeterminantFilter.h>


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    if(argc != 6) {
        std::cerr<<"usage: exec inputMaskFileName inputAtrophyFileName inputDisplacementFileName outputMaskFileName outputAtrophyFileName"<<std::endl;
        return(EXIT_FAILURE);
    }

    std::string maskFileName(argv[1]);
    std::string atrophyFileName(argv[2]);
    std::string displacementFileName(argv[3]);
    std::string outputMaskFileName(argv[4]);
    std::string outputAtrophyFileName(argv[5]);

    bool scaleDisplacementField(true);
    bool invertField(true), warpMask(true);
    bool warpAtrophy(true), writeJacobian(true);

    //---------------------  Read the displacement field ----------------------//
    typedef itk::Vector<double, 3>                      DisplacementPixelType;
    typedef itk::Image<DisplacementPixelType, 3>        DisplacementImageType;
    typedef itk::ImageFileReader<DisplacementImageType> DisplacementReaderType;
    DisplacementImageType::Pointer displacementField;
    {
        DisplacementReaderType::Pointer displacementReader = DisplacementReaderType::New();
        displacementReader->SetFileName(displacementFileName);
        displacementReader->Update();
        displacementField = displacementReader->GetOutput();

        if(scaleDisplacementField){
            double scaleFactor = 5;
            itk::ImageRegionIterator<DisplacementImageType> iterator(displacementField,displacementField->GetLargestPossibleRegion());
            iterator.GoToBegin();
            DisplacementImageType::PixelType pixelValue;
            while(!iterator.IsAtEnd()) {
                pixelValue = iterator.Get();
                pixelValue[0] = pixelValue[0]*scaleFactor;
                pixelValue[1] = pixelValue[1]*scaleFactor;
                pixelValue[2] = pixelValue[2]*scaleFactor;
                iterator.Set(pixelValue);
                ++iterator;
            }
        }
    }

    //----------------------- Invert the displacement field ----------------------//
    DisplacementImageType::Pointer inverseDisplacementField;
    if(invertField) {
        {
            typedef InverseDisplacementImageFilter<DisplacementImageType> FPInverseType;
            FPInverseType::Pointer inverter1 = FPInverseType::New();
            inverter1->SetInput(displacementField);
            inverter1->SetErrorTolerance(1e-1);
            inverter1->SetMaximumNumberOfIterations(50);
            inverter1->Update();
            std::cout<<"tolerance not reached for "<<inverter1->GetNumberOfErrorToleranceFailures()<<" pixels"<<std::endl;
            inverseDisplacementField = inverter1->GetOutput();

            typedef itk::ImageFileWriter<DisplacementImageType> WriterType;
            WriterType::Pointer writer = WriterType::New();
            writer->SetFileName("velScaled.mha");
            //        writer->SetFileName("step1Vel10timesScaled.mha");
            writer->SetInput(inverseDisplacementField);
            writer->Update();

        }
    }

    if(warpMask) {
        //-------------------- Read and warp the mask image -----------------------------//
        typedef int                                     MaskPixelType;
        //    typedef double                                     MaskPixelType;
        typedef itk::Image<MaskPixelType, 3>            MaskImageType;
        typedef itk::ImageFileReader<MaskImageType>     MaskReaderType;
        MaskImageType::Pointer  mask;
        {
            MaskReaderType::Pointer maskReader = MaskReaderType::New();
            maskReader->SetFileName(maskFileName);
            maskReader->Update();
            mask = maskReader->GetOutput();

            typedef itk::WarpImageFilter<MaskImageType,MaskImageType,DisplacementImageType> WarpFilterType;
            WarpFilterType::Pointer warper = WarpFilterType::New();
            warper->SetDisplacementField(inverseDisplacementField);
            warper->SetInput(mask);
            warper->SetOutputSpacing(mask->GetSpacing());
            warper->SetOutputOrigin(mask->GetOrigin());
            warper->SetOutputDirection(mask->GetDirection());

            typedef itk::NearestNeighborInterpolateImageFunction<MaskImageType> InterpolatorFilterType;
            //        typedef itk::BSplineInterpolateImageFunction<MaskImageType> InterpolatorFilterType;
            InterpolatorFilterType::Pointer interpolatorFilter = InterpolatorFilterType::New();
            //        interpolatorFilter->SetSplineOrder(3);
            warper->SetInterpolator(interpolatorFilter);
            warper->Update();

            typedef itk::ImageFileWriter<MaskImageType> WriterType;
            WriterType::Pointer writer = WriterType::New();
            writer->SetFileName(outputMaskFileName);
            writer->SetInput(warper->GetOutput());
            writer->Update();
        }
    }

    //--------------------  Read and warp atrophy map modulating with the jacobian determinant -------------------------//
    if(warpAtrophy) {
        typedef double                                  atrophyPixelType;
        typedef itk::Image<atrophyPixelType, 3>         AtrophyImageType;
        typedef itk::ImageFileReader<AtrophyImageType>  AtrophyReaderType;
        AtrophyImageType::Pointer   atrophy;
        {
            AtrophyReaderType::Pointer atrophyReader = AtrophyReaderType::New();
            atrophyReader->SetFileName(atrophyFileName);
            atrophyReader->Update();
            atrophy = atrophyReader->GetOutput();

            //--------------------- Warp atrophy map ------------------------------//
            typedef itk::WarpImageFilter<AtrophyImageType,AtrophyImageType,DisplacementImageType> WarpFilterType;
            WarpFilterType::Pointer warper = WarpFilterType::New();
            warper->SetDisplacementField(inverseDisplacementField);
            warper->SetInput(atrophy);
            warper->SetOutputSpacing(atrophy->GetSpacing());
            warper->SetOutputOrigin(atrophy->GetOrigin());
            warper->SetOutputDirection(atrophy->GetDirection());

            //        typedef itk::BSplineInterpolateImageFunction<AtrophyImageType> InterpolatorFilterType;
//            InterpolatorFilterType::Pointer interpolatorFilter = InterpolatorFilterType::New();
            //        interpolatorFilter->SetSplineOrder(3);
//            warper->SetInterpolator(interpolatorFilter);
            warper->Update();

            //-------------------- Jacobian determinant of the input field ----------------------------------//
            typedef itk::DisplacementFieldJacobianDeterminantFilter< DisplacementImageType, double > JacobianFilterType;
            JacobianFilterType::Pointer jacobianFilter = JacobianFilterType::New();
            jacobianFilter->SetUseImageSpacingOff();
            jacobianFilter->SetInput(displacementField);
            jacobianFilter->Update();

            if(writeJacobian) {
                typedef itk::ImageFileWriter<AtrophyImageType> WriterType;
                WriterType::Pointer writer = WriterType::New();
                writer->SetFileName("jacobian.mha");
                writer->SetInput(jacobianFilter->GetOutput());
                writer->Update();
            }

            //----------------------- modulate the warped atrophy map with the jacobian --------------------//
//            typedef itk::MultiplyImageFilter< AtrophyImageType >    ModulatorType;
            typedef itk::DivideImageFilter< AtrophyImageType, AtrophyImageType, AtrophyImageType > ModulatorType;
            ModulatorType::Pointer modulator = ModulatorType::New();
            modulator->SetInput1(warper->GetOutput());
            modulator->SetInput2(jacobianFilter->GetOutput());
            modulator->Update();

            //---------------------- write the warped with modulation atrophy map -----------------------//
            typedef itk::ImageFileWriter<AtrophyImageType> WriterType;
            WriterType::Pointer writer = WriterType::New();
            writer->SetFileName(outputAtrophyFileName);
//            writer->SetInput(warper->GetOutput());
            writer->SetInput(modulator->GetOutput());
            writer->Update();
        }

        //use it to warp atrophy.
        return(EXIT_SUCCESS);
    }
}
