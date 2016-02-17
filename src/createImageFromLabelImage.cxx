#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/program_options.hpp>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkMaskImageFilter.h>
#include <itkBinaryImageToLabelMapFilter.h>
#include <itkLabelStatisticsImageFilter.h>

int main(int argc,char **argv)
{
    std::string     labelTableFile, labelImageFile, outImageFile, fileToModify;
    bool modifyExisting;

    boost::program_options::options_description optionsDescription("create an image from a label image and a lable table.");
    optionsDescription.add_options()
	("help,h", "display help message")
	("table,t", boost::program_options::value< std::string >(&labelTableFile), "label table with two columns: labels and corresponding pixel values to be assigned to the output image.")
	("labelImage,l", boost::program_options::value< std::string >(&labelImageFile), "input label image file from which the regions matching the labels in the table will be extracted"
	 "New pixel values will be put to only these extracted regions.")
	("outputImage,o",boost::program_options::value< std::string >(&outImageFile), "Output image")
	("fileToModify,m",boost::program_options::value< std::string >(&fileToModify), " The image that is to be modified to produce output image."
	 " If not provided, creates a new image with zero value everywhere"
	 " and then updates the extracted regions.")
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
    if(!options.count("table") || !options.count("labelImage")
       || !options.count("outputImage")
	) {
	std::cerr<<"invalid options! run with --help or -h to see the proper options."<<std::endl;
	return EXIT_FAILURE;
    }

// Check if file to be modified is given or not.
    if(!options.count("fileToModify"))
	modifyExisting = false;
    else
	modifyExisting = true;
    typedef itk::Image<int, 3> LabelImageType;
    typedef itk::Image<float, 3> OutImageType;

// Open label table file for reading
    std::ifstream labelTable(labelTableFile.c_str(),std::ios::in);
    if (!labelTable.is_open()) {
	std::cerr<<"could not open file: "<<labelTableFile<<std::endl;
	return EXIT_FAILURE;
    }

// Read first line of the table and check if it is in the desired format i.e. labels newValues
    {
	std::string line;
	std::getline(labelTable, line);
	std::istringstream is(line);
	std::string checkString;
	is >> checkString;
	if(checkString.compare("labels")) {
	    std::cout<<"incorrect table format. 1st col- labels; 2nd col- newValues"<<std::endl<<"The first word must be: labels"<<std::endl;
	    return EXIT_FAILURE;
	}
	is >> checkString;
	if (checkString.compare("newValues")) {
	    std::cout<<"incorrect table format. 1st col- labels; 2nd col- newValues"<<std::endl<<"The second word must be: newValues"<<std::endl;
	    return EXIT_FAILURE;
	}
    }

// From the second line, read the values and put it into a map.
    typedef std::map< LabelImageType::PixelType, OutImageType::PixelType > LabelWithValueType;
    typedef std::pair< LabelImageType::PixelType, OutImageType::PixelType > LabelValuePairType;
    typedef LabelWithValueType::iterator LabelWithValueIteratorType;

    LabelWithValueType labelWithValue;
    while(!labelTable.eof()) {
	std::string line;
	std::getline(labelTable, line);
	std::istringstream is(line);
	LabelImageType::PixelType label;
	is >> label;
	OutImageType::PixelType outValue;
	is >> outValue;
	labelWithValue.insert(LabelValuePairType(label, outValue));
    }
    labelTable.close();

    for (LabelWithValueIteratorType labelWithValueIt = labelWithValue.begin(); labelWithValueIt != labelWithValue.end(); ++labelWithValueIt) {
	std::cout<<labelWithValueIt->first<<" ** "<<labelWithValueIt->second<<std::endl;
    }

// Read label image
    typedef itk::ImageFileReader< LabelImageType > LabelImageReaderType;
    LabelImageReaderType::Pointer labelImageReader = LabelImageReaderType::New();
    labelImageReader->SetFileName(labelImageFile);
    labelImageReader->Update();
    LabelImageType::Pointer labelImage = labelImageReader->GetOutput();

// Output image.
    OutImageType::Pointer outImage = OutImageType::New();
    if (modifyExisting) {
	typedef itk::ImageFileReader< OutImageType > OutImageReaderType;
	OutImageReaderType::Pointer outImageReader = OutImageReaderType::New();
	outImageReader->SetFileName(fileToModify);
	outImageReader->Update();
	outImage = outImageReader->GetOutput();
    } else { // If not modifying create new image of the same size as label image and fill it with zero.
	outImage->SetRegions(labelImage->GetLargestPossibleRegion());
	outImage->CopyInformation(labelImage);
	outImage->Allocate();
	outImage->FillBuffer(0.);
    }



    //L:labelImage A:outImage  out(I):out value for the corresponding label I read from the table.
    //if ( L(x) != 0 ) then A(x) = out(L(x)) or just A(x) = out(L(x)) without if; when you want for label 0 also.
    typedef itk::ImageRegionConstIterator< LabelImageType > LabelIteratorType;
    typedef itk::ImageRegionIterator< OutImageType > ValueIteratorType;
    LabelIteratorType labelIt(labelImage, labelImage->GetLargestPossibleRegion());
    ValueIteratorType outIt(outImage, outImage->GetLargestPossibleRegion());
    labelIt.GoToBegin();
    outIt.GoToBegin();
    while(!labelIt.IsAtEnd())
    {
	LabelImageType::PixelType currLabel = labelIt.Get();
	//std::cout<<"LabelImage[x]="<<currLabel<<std::endl;
	// if(currLabel) //for all non-zero labels
	// {
	for (LabelWithValueIteratorType labelWithValueIt = labelWithValue.begin(); labelWithValueIt != labelWithValue.end(); ++labelWithValueIt)
	{
	    //std::cout<<"LabelImage[x]="<<currLabel<<",  labelTable[current]="<<labelWithValueIt->first<<std::endl;
	    if (currLabel == labelWithValueIt->first)
	    {
		//std::cout<<"LabelImage[x]="<<currLabel<<",  labelTable[current]="<<labelWithValueIt->first<<std::endl;
		outIt.Set(labelWithValueIt->second);
		break;
	    }
	}
	//}
	++labelIt;
	++outIt;
    }

// //Label Statistics Filter to write to only those labels which are present in the input label image.
// typedef itk::BinaryImageToLabelMapFilter< LabelImageType > ImageToLabelMapFilterType;
// ImageToLabelMapFilterType::Pointer labelMap = ImageToLabelMapFilterType::New();
// labelMap->SetInput(labelImage);
// labelMap->Update();
// std::cout<<"There are originally "<< labelMap->GetOutput()->GetNumberOfLabelObjects() << " objects"<<std::endl;

// //Label Statistics Filter to compute mean of each regions:
// typedef itk::LabelStatisticsImageFilter< OutImageType, LabelImageType > LabelStatsFilterType;
// LabelStatsFilterType::Pointer labelStats = LabelStatsFilterType::New();
// labelStats->SetLabelInput( labelImage);
// labelStats->SetInput(outImage);
// labelStats->Update();
// std::cout << "Number of labels in the input label image: " << labelStats->GetNumberOfLabels() << std::endl;

// typedef itk::MaskImageFilter<OutImageType, LabelImageType, OutImageType> MaskImageFilterType;
// MaskImageFilterType::Pointer maskImageFilter = MaskImageFilterType::New();
// maskImageFilter->SetMaskImage(labelImage);
// for (LabelWithValueIteratorType labelWithValueIt = labelWithValue.begin(); labelWithValueIt != labelWithValue.end(); ++labelWithValueIt)
// {
// 	LabelImageType::PixelType label = labelWithValueIt->first;
// 	OutImageType::PixelType val = labelWithValueIt->second;
// 	std::cout<<"Testing for label "<<label<<" with value "<<val<<std::endl;
// 	if(labelStats->HasLabel(label))//if(labelMap->GetOutput()->HasLabel(label))
// 	{
// 	    maskImageFilter->SetInput(outImage);
// 	    //if (pixel_from_mask_image != masking_value)     pixel_output_image = pixel_input_image
// 	    //else pixel_output_image = outside_value
// 	    maskImageFilter->SetMaskingValue(label);
// 	    maskImageFilter->SetOutsideValue(val);
// 	    maskImageFilter->Update();
// 	    outImage = maskImageFilter->GetOutput();
// 	    outImage->DisconnectPipeline();
// 	}
// 	else {
// 	    std::cout<<" Label Id "<<label<<" not present in the label image"<<std::endl;
// 	}
// }

// Write the out output.
    typedef itk::ImageFileWriter< OutImageType > ValueImageWriterType;
    ValueImageWriterType::Pointer outImageWriter = ValueImageWriterType::New();
    outImageWriter->SetFileName(outImageFile);
    outImageWriter->SetInput(outImage);
    outImageWriter->Update();
    return EXIT_SUCCESS;
}

