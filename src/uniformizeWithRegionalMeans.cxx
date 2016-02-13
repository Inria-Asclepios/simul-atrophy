#include <iostream>
#include <fstream>
//#include <sstream>
#include <iterator>

#include <boost/program_options.hpp>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkLabelStatisticsImageFilter.h>
// #include <itkImageRegionIterator.h>
// #include <itkImageRegionConstIterator.h>

int main(int argc,char **argv)
{
    std::string     labelTableFile, labelImageFile, inImageFile, outFile;

    boost::program_options::options_description optionsDescription("Write a file with total and average Jacobians from a given label table, label image and Jacobian image.");
    optionsDescription.add_options()
	("help,h", "display help message")
	("table,t", boost::program_options::value< std::string >(&labelTableFile), "labels in a single column for which sum and average of inImage is to be computed.")
	("inImage,i", boost::program_options::value< std::string >(&inImageFile), "Input intensity image whose mean regional values are to be computed.")
	("labelImage,l", boost::program_options::value< std::string >(&labelImageFile), "input label image file from which the regions matching the labels in the table will be extracted")
	("outFile,o",boost::program_options::value< std::string >(&outFile), "Output filename")
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
    if(!options.count("table") || !options.count("inImage")
       || !options.count("labelImage") || !options.count("outFile")
	) {
	std::cerr<<"invalid options! run with --help or -h to see the proper options."<<std::endl;
	return EXIT_FAILURE;
    }

    typedef itk::Image<int, 3> LabelImageType;
    typedef itk::Image<float, 3> InImageType;

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

// Read label image
    typedef itk::ImageFileReader< LabelImageType > LabelImageReaderType;
    LabelImageReaderType::Pointer labelImageReader = LabelImageReaderType::New();
    labelImageReader->SetFileName(labelImageFile);
    //labelImageReader->Update();
    //LabelImageType::Pointer labelImage = labelImageReader->GetOutput();

// Input intensity image.
    InImageType::Pointer inImage = InImageType::New();
    typedef itk::ImageFileReader< InImageType > InImageReaderType;
    InImageReaderType::Pointer inImageReader = InImageReaderType::New();
    inImageReader->SetFileName(inImageFile);
    //inImageReader->Update();
    //inImage = inImageReader->GetOutput();

    //Label Statistics Filter:
    typedef itk::LabelStatisticsImageFilter< InImageType, LabelImageType > LabelStatsFilterType;
    LabelStatsFilterType::Pointer labelStats = LabelStatsFilterType::New();

    labelStats->SetLabelInput( labelImageReader->GetOutput() );
    labelStats->SetInput(inImageReader->GetOutput());
    labelStats->Update();

    std::cout << "Number of labels in the input label image: " << labelStats->GetNumberOfLabels() << std::endl;
    typedef LabelStatsFilterType::ValidLabelValuesContainerType ValidLabelValuesType;
    typedef LabelStatsFilterType::LabelPixelType                LabelPixelType;
    // for(ValidLabelValuesType::const_iterator vIt=labelStats->GetValidLabelValues().begin();
    // 	vIt != labelStats->GetValidLabelValues().end();
    // 	++vIt)
    std::ofstream outputFile(outFile.c_str());
    outputFile<<"LabelId\tMeanJacobian\tNumberOfVoxels\n";
    for (LabelsVecItType it = inLabels.begin(); it != inLabels.end(); ++it) {
	if ( labelStats->HasLabel(*it) ) {
	    LabelPixelType labelValue = *it;
	    // std::cout <<"For Label Id: "<<labelValue<<std::endl;
	    // std::cout << "mean: " << labelStats->GetMean( labelValue ) << std::endl;
	    // std::cout << "count: " << labelStats->GetCount( labelValue ) << std::endl;
	    // std::cout << std::endl;
	    outputFile<<labelValue<<"\t"<<labelStats->GetMean(labelValue)<<"\t"<<labelStats->GetCount(labelValue)<<'\n';
	}
	else {
	    std::cout<<" Label Id "<<*it<<" not present in the label image"<<std::endl;
	}
    }
    outputFile.close();

    return EXIT_SUCCESS;
}

