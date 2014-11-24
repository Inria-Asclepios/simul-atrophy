#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkGradientAnisotropicDiffusionWithMaskImageFilter.h"

int main( int argc, char* argv[] )
{
  if( argc != 7 )
    {
    std::cerr << "Usage: "<< std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFileName>";
    std::cerr << " <CondutanceImageFileName>";
    std::cerr << " <OutputFileName>";
    std::cerr << " <NumberOfIterations> ";
    std::cerr << " <Conductance>";
    std::cerr << " <timeStep>" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  std::string   inputImageFile(argv[1]);
  std::string   conductanceImageFile(argv[2]);
  std::string   outputImageFile(argv[3]);
  int           numOfIterations(atoi(argv[4]));
  float         conductance(atof(argv[5]));
  float         timeStep(atof(argv[6]));

//  typedef unsigned char                           InputPixelType;
  typedef float                                     InputPixelType;
  typedef itk::Image< InputPixelType, Dimension >   InputImageType;

  typedef itk::ImageFileReader< InputImageType >    ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFile );

  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( conductanceImageFile );
  try {
      reader1->Update();
  } catch (itk::ExceptionObject & error)
  {
      std::cerr << "Error: " << error << std::endl;
      return EXIT_FAILURE;
  }


//  reader1->Update();
//  InputImageType::Pointer tst = reader1->GetOutput();


  typedef float                                     OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension >  OutputImageType;
//  typedef itk::GradientAnisotropicDiffusionImageFilter< InputImageType,
//    OutputImageType > FilterType;
  typedef itk::GradientAnisotropicDiffusionWithMaskImageFilter< InputImageType,
    OutputImageType > FilterType;

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetConductanceImage( reader1->GetOutput() );
  filter->SetNumberOfIterations( numOfIterations );
//  filter->SetTimeStep( 0.125 );
  filter->SetTimeStep( timeStep );
  filter->SetConductanceParameter( conductance );
  try {
      filter->Update();
  }
  catch (itk::ExceptionObject & error)
  {
      std::cerr << "Error: " << error << std::endl;
      return EXIT_FAILURE;
  }

  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputImageFile );
  writer->SetInput( filter->GetOutput() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
