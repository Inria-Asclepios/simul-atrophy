#include <iostream>
#include <string>
#include <sstream>

#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIteratorWithIndex.hxx>
#include <itkDiffusionTensor3D.hxx>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

//#include <boost/program_options/parsers.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    std::string workingDir("/user/bkhanal/home/works/AdLemModel/");
    /*std::string inImgFile;
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

    //Parse the sphereCenter into number from the string:
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
    */
    //------- Test reading tensor images and multiply tensor by a vector ------

    typedef typename itk::Image<double, 3>                              ScalarImageType;
    typedef typename itk::Image<itk::DiffusionTensor3D< double >, 3>    TensorImageType;
    typedef typename itk::Image<itk::Vector<double,3>, 3>               VectorImageType;

    typedef typename itk::ImageFileReader<ScalarImageType>              ScalarImageReaderType;
    typedef typename itk::ImageFileReader<TensorImageType>              TensorImageReaderType;
    typedef typename itk::ImageFileReader<VectorImageType>              VectorImageReaderType;
    typedef typename itk::ImageFileWriter<TensorImageType>              TensorImageWriterType;
    typedef typename itk::ImageFileWriter<VectorImageType>              VectorImageWriterType;

    ScalarImageReaderType::Pointer      scalarImageReader = ScalarImageReaderType::New();
//    scalarImageReader->SetFileName(workingDir + "results/patients/test/t1.nii.gz");
//    scalarImageReader->SetFileName("workingDir + "results/patients/test/cyl.mha");
    scalarImageReader->SetFileName(workingDir + "results/patients/test/aL.mha");
    scalarImageReader->Update();
    ScalarImageType::Pointer            scalarImage = ScalarImageType::New();
    scalarImage = scalarImageReader->GetOutput();

    TensorImageReaderType::Pointer      tensorImageReader = TensorImageReaderType::New();
//    tensorImageReader->SetFileName(workingDir + "results/patients/test/tFxL10L4L1.nii");
//    tensorImageReader->SetFileName(workingDir + "results/patients/test/tUzL10L4L1.nii");
//    tensorImageReader->SetFileName(workingDir + "results/patients/test/tFxL2L0_5L0_1Matlab.nii");
//    tensorImageReader->SetFileName(workingDir + "results/patients/test/tUzL2L0_5L0_1Matlab.nii");
    tensorImageReader->SetFileName(workingDir + argv[1]);

    tensorImageReader->Update();

    TensorImageType::Pointer            tensorImage = TensorImageType::New();
    tensorImage = tensorImageReader->GetOutput();

    VectorImageReaderType::Pointer      vectorImageReader = VectorImageReaderType::New();
    vectorImageReader->SetFileName(workingDir + "results/patients/test/aL_m3_scalar_force.nii");
//        vectorImageReader->SetFileName(workingDir + "results/patients/test/aL_m3_tUzL10L4L1_force.nii.gz");

    //    vectorImageReader->SetFileName(workingDir + "results/patients/test/vec001.nii");
//        vectorImageReader->SetFileName(workingDir + "results/patients/test/vec100.nii");
//        vectorImageReader->SetFileName(workingDir + "results/patients/test/vec010.nii");
    vectorImageReader->Update();

    VectorImageType::Pointer            vectorImage = VectorImageType::New();
    vectorImage = vectorImageReader->GetOutput();

    TensorImageType::IndexType posTensor;
    for (unsigned int i = 0; i<3; ++i) posTensor.SetElement(i,10);
    std::cout<<"Tensor value: "<<tensorImage->GetPixel(posTensor)<<std::endl;
    std::cout<<"Tensor values again:\n"<<tensorImage->GetPixel(posTensor)(0,0)<<"\t"<<tensorImage->GetPixel(posTensor)(0,1)<<"\t"<<tensorImage->GetPixel(posTensor)(0,2)<<std::endl;
    std::cout<<tensorImage->GetPixel(posTensor)(1,0)<<"\t"<<tensorImage->GetPixel(posTensor)(1,1)<<"\t"<<tensorImage->GetPixel(posTensor)(1,2)<<std::endl;
    std::cout<<tensorImage->GetPixel(posTensor)(2,0)<<"\t"<<tensorImage->GetPixel(posTensor)(2,1)<<"\t"<<tensorImage->GetPixel(posTensor)(2,2)<<std::endl;
//    std::cout<<"Vector value: "<<vectorImage->GetPixel(posTensor)<<std::endl;
//    std::cout<<"origin of vector image: "<<vectorImage->GetOrigin()<<std::endl;
//    std::cout<<"direction of vector image:\n "<<vectorImage->GetDirection()<<std::endl;
//    std::cout<<"origin of scalar image: "<<scalarImage->GetOrigin()<<std::endl;
//    std::cout<<"direction of scalar image:\n "<<scalarImage->GetDirection()<<std::endl;

    typedef itk::ImageRegionIterator<VectorImageType> VectorIteratorType;
//    typedef itk::ImageRegionConstIteratorWithIndex<VectorImageType> VectorIteratorType;

    VectorIteratorType it(vectorImage,vectorImage->GetLargestPossibleRegion());


    VectorImageType::PixelType vecPixel;
    vecPixel[0] = 1;
    vecPixel[1] = 1;
    vecPixel[2] = 1;

  /*  std::cout<<"vector value: "<<vectorImage->GetPixel(posTensor)<<std::endl;
    it.GoToBegin();
    while(!it.IsAtEnd()) {
        if (it.Get()[0] != 0 || it.Get()[1] != 0 || it.Get()[2] != 0){
            std::cout<<it.Get()<<std::endl;
//        std::cout<<it.Value()<<std::endl;
        }
//        it.Set(vecPixel);
        ++it;
    }*/
//    VectorImageWriterType::Pointer vectorImageWriter = VectorImageWriterType::New();
//    vectorImageWriter->SetFileName(workingDir + "results/patients/test/vec100.nii");
//    vectorImageWriter->SetFileName(workingDir + "results/patients/test/vec010.nii");
//    vectorImageWriter->SetFileName(workingDir + "results/patients/test/vec001.nii");
//    vectorImageWriter->SetFileName(workingDir + "results/patients/test/vec111.nii");
//    vectorImageWriter->SetInput(vectorImage);
//    vectorImageWriter->Update();


    return EXIT_SUCCESS;

}
