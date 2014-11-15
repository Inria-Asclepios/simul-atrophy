#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"

int main( int argc, char* argv[] )
{
  if( argc != 6 )
    {
    std::cerr << "Usage: "<< std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFileName>";
    std::cerr << " <OutputFileName>";
    std::cerr << " <NumberOfIterations> ";
    std::cerr << " <Conductance>";
    std::cerr << " <timeStep>" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  std::string   inputImageFile(argv[1]);
  std::string   outputImageFile(argv[2]);
  int           numOfIterations(atoi(argv[3]));
  float         conductance(atof(argv[4]));
  float         timeStep(atof(argv[5]));

//  typedef unsigned char                           InputPixelType;
  typedef float                                     InputPixelType;
  typedef itk::Image< InputPixelType, Dimension >   InputImageType;

  typedef itk::ImageFileReader< InputImageType >    ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFile );

  typedef float                                     OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension >  OutputImageType;
  typedef itk::GradientAnisotropicDiffusionImageFilter< InputImageType,
    OutputImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetNumberOfIterations( numOfIterations );
//  filter->SetTimeStep( 0.125 );
  filter->SetTimeStep( timeStep );
  filter->SetConductanceParameter( conductance );

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
