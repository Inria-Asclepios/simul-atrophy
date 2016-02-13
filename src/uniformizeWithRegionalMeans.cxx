#include <iostream>
#include <fstream>
#include <iterator>

#include <boost/program_options.hpp>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkLabelStatisticsImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkMaskImageFilter.h>
#include <itkImageDuplicator.h>

int main(int argc,char **argv)
{
    std::string labelTableFile, labelImageFile, inImageFile, outTableFile, outImageFile;
    bool	writeOutTableFile(false), writeOutImageFile(false), startFromInImage(false);

    boost::program_options::options_description optionsDescription("Write a file with total and average Jacobians from a given label table, label image and Jacobian image.");
    optionsDescription.add_options()
	("help,h", "display help message")
	("labelTable,t", boost::program_options::value< std::string >(&labelTableFile), "labels in a single column for which sum and average of inImage is to be computed.")
	("inImage,i", boost::program_options::value< std::string >(&inImageFile), "Input intensity image whose mean regional values are to be computed.")
	("labelImage,l", boost::program_options::value< std::string >(&labelImageFile), "input label image file from which the regions matching the labels in the table will be extracted")
	("outTable,o",boost::program_options::value< std::string >(&outTableFile), "[Optional] If given create this table file which contains table with labels, mean values and number of voxels for the given label")
	("outImage,u",boost::program_options::value< std::string >(&outImageFile), "[Optional]. If given creates this image which has uniformized mean intensity in each ROIs corresponding to the given labels and labelImage. The mean intensity is computed from the given input intensity image. ")
	("startFromInImage,s",boost::program_options::bool_switch(&startFromInImage)->default_value(false), "If given the ouImage copies the intensity values of inImage for all ROIs not covered by the input labels; otherwise all regions not covered by labels will have zero value.")
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
    if(!options.count("labelTable") || !options.count("inImage") || !options.count("labelImage")
	)
    {
	std::cerr<<"invalid options! run with --help or -h to see the proper options."<<std::endl;
	return EXIT_FAILURE;
    }
    if(options.count("outTable"))
	writeOutTableFile = true;
    if(options.count("outImage"))
	writeOutImageFile = true;
    if(!(writeOutTableFile || writeOutImageFile))
    {
	std::cerr<<"neither of outTable and outImage given. Nothing to write => I do nothing. Provide at lease one of these filenames to write.";
	return EXIT_FAILURE;

    }
    typedef itk::Image<int, 3> LabelImageType;
    typedef itk::Image<float, 3> ImageType;
    // Container to put labels.
    typedef std::vector<LabelImageType::PixelType> LabelsVecType;
    typedef LabelsVecType::iterator LabelsVecItType;

// Open label table file for reading
    std::ifstream labelTable(labelTableFile.c_str(),std::ios::in);
    if (!labelTable.is_open()) {
	std::cerr<<"could not open file: "<<labelTableFile<<std::endl;
	return EXIT_FAILURE;
    }
    //Read all the numbers and put them into a vector.
    std::istream_iterator<double> start(labelTable), end;
    //std::vector<LabelImageType::PixelType> inLabels(start, end);
    LabelsVecType inLabels(start, end);
    std::cout << "Read " << inLabels.size() << " labels" << std::endl;
    labelTable.close();

    // for (LabelsVecItType it = inLabels.begin(); it != inLabels.end(); ++it) {
    // 	std::cout<<*it<<std::endl;
    // }

// Read input label image
    LabelImageType::Pointer labelImage;
    {
	typedef itk::ImageFileReader< LabelImageType > LabelImageReaderType;
	LabelImageReaderType::Pointer labelImageReader = LabelImageReaderType::New();
	labelImageReader->SetFileName(labelImageFile);
	labelImageReader->Update();
	labelImage = labelImageReader->GetOutput();
    }

// Read input intensity image.
    ImageType::Pointer inImage;
    {
	typedef itk::ImageFileReader< ImageType > InImageReaderType;
	InImageReaderType::Pointer inImageReader = InImageReaderType::New();
	inImageReader->SetFileName(inImageFile);
	inImageReader->Update();
	inImage = inImageReader->GetOutput();
    }

    //Label Statistics Filter to compute mean of each regions:
    typedef itk::LabelStatisticsImageFilter< ImageType, LabelImageType > LabelStatsFilterType;
    LabelStatsFilterType::Pointer labelStats = LabelStatsFilterType::New();

    labelStats->SetLabelInput( labelImage);
    labelStats->SetInput(inImage);
    labelStats->Update();

    std::cout << "Number of labels in the input label image: " << labelStats->GetNumberOfLabels() << std::endl;
    typedef LabelStatsFilterType::ValidLabelValuesContainerType ValidLabelValuesType;
    typedef LabelStatsFilterType::LabelPixelType                LabelPixelType;
    typedef ImageType::PixelType				OutPixelType;
    // for(ValidLabelValuesType::const_iterator vIt=labelStats->GetValidLabelValues().begin();
    // 	vIt != labelStats->GetValidLabelValues().end();
    // 	++vIt)


    // For writing into a table
    if(writeOutTableFile)
    {
	std::ofstream outputFile(outTableFile.c_str());
	outputFile<<"LabelId\tMeanIntensity\tNumberOfVoxels\n";
	for (LabelsVecItType it = inLabels.begin(); it != inLabels.end(); ++it)
	{
	    if ( labelStats->HasLabel(*it) )
	    {
		LabelPixelType labelValue = *it;
		OutPixelType meanValue = labelStats->GetMean(labelValue);
		unsigned int labelVoxelCount = labelStats->GetCount(labelValue);
		// std::cout <<"For Label Id: "<<labelValue<<std::endl;
		// std::cout << "mean: " << meanValue << std::endl;
		// std::cout << "count: " << labelVoxelCount << std::endl;
		// std::cout << std::endl;
		outputFile<<labelValue<<"\t"<<meanValue<<"\t"<<labelVoxelCount<<'\n';
	    }
	    else
	    {
		std::cout<<" Label Id "<<*it<<" not present in the label image"<<std::endl;
	    }
	}
	outputFile.close();
    }

    if(writeOutImageFile)
    { // For writing an image
        // Output image, initialize with zero.
	ImageType::Pointer outImage;
	if(startFromInImage)
	{
	    typedef itk::ImageDuplicator< ImageType > DuplicatorType;
	    DuplicatorType::Pointer duplicator = DuplicatorType::New();
	    duplicator->SetInputImage(inImage);
	    duplicator->Update();
	    outImage = duplicator->GetModifiableOutput();
	}
	else
	{
	    outImage = ImageType::New();
	    outImage->SetRegions(inImage->GetLargestPossibleRegion());
	    outImage->CopyInformation(inImage);
	    outImage->Allocate();
	    outImage->FillBuffer(0.);
	}

	typedef itk::MaskImageFilter<ImageType, LabelImageType, ImageType> MaskImageFilterType;
	MaskImageFilterType::Pointer maskImageFilter = MaskImageFilterType::New();
	maskImageFilter->SetMaskImage(labelImage);
	for (LabelsVecItType it = inLabels.begin(); it != inLabels.end(); ++it)
	{
	    if ( labelStats->HasLabel(*it) )
	    {
		LabelPixelType labelValue = *it;
		OutPixelType meanValue = labelStats->GetMean(labelValue);
		maskImageFilter->SetInput(outImage);
		//if (pixel_from_mask_image != masking_value)     pixel_output_image = pixel_input_image
		//else pixel_output_image = outside_value
		maskImageFilter->SetMaskingValue(labelValue);
		maskImageFilter->SetOutsideValue(meanValue);
		maskImageFilter->Update();
		outImage = maskImageFilter->GetOutput();
		outImage->DisconnectPipeline();
	    }
	    else
	    {
		std::cout<<" Label Id "<<*it<<" not present in the label image"<<std::endl;
	    }
	}
	// Write the out output.
	typedef itk::ImageFileWriter< ImageType > ValueImageWriterType;
	ValueImageWriterType::Pointer outImageWriter = ValueImageWriterType::New();
	outImageWriter->SetFileName(outImageFile);
	outImageWriter->SetInput(outImage);
	outImageWriter->Update();
    }

    return EXIT_SUCCESS;
}

