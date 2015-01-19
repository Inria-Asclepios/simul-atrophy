#include <iostream>
#include <fstream>
#include <sstream>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

int main(int argc,char **argv)
{
    if(argc != 4) {
	std::cerr<<"usage: "<<argv[0]<<" labelTableInput labelImageInput atrophyMapOutput"<<std::endl;
	return EXIT_FAILURE;
    }
    std::string     labelTableFile(argv[1]), labelImageFile(argv[2]), atrophyMapFile(argv[3]);

    typedef itk::Image<int, 3> LabelImageType;
    typedef itk::Image<float, 3> AtrophyImageType;

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
    typedef std::map< LabelImageType::PixelType, AtrophyImageType::PixelType > LabelWithAtrophyType;
    typedef std::pair< LabelImageType::PixelType, AtrophyImageType::PixelType > LabelAtrophyPairType;
    typedef LabelWithAtrophyType::iterator LabelWithAtrophyIteratorType;

    LabelWithAtrophyType labelWithAtrophy;
    while(!labelTable.eof()) {
	std::string line;
	std::getline(labelTable, line);
	std::istringstream is(line);
	LabelImageType::PixelType label;
	is >> label;
	AtrophyImageType::PixelType atrophyValue;
	is >> atrophyValue;
	labelWithAtrophy.insert(LabelAtrophyPairType(label, atrophyValue));
    }
    labelTable.close();

    // for (LabelWithAtrophyIteratorType labelWithAtrophyIt = labelWithAtrophy.begin(); labelWithAtrophyIt != labelWithAtrophy.end(); ++labelWithAtrophyIt) {
    // 	std::cout<<labelWithAtrophyIt->first<<" ** "<<labelWithAtrophyIt->second<<std::endl;
    // }

    // Read label image
    typedef itk::ImageFileReader< LabelImageType > LabelImageReaderType;
    LabelImageReaderType::Pointer labelImageReader = LabelImageReaderType::New();
    labelImageReader->SetFileName(labelImageFile);
    labelImageReader->Update();
    LabelImageType::Pointer labelImage = labelImageReader->GetOutput();

    // Create atrophy image of the same size as label image and fill it with zero.
    AtrophyImageType::Pointer atrophyImage = AtrophyImageType::New();
    atrophyImage->SetRegions(labelImage->GetLargestPossibleRegion());
    atrophyImage->CopyInformation(labelImage);
    // :FIXME: is copyinformation same as setting spacing, origin and orientation and/or largestregion ? or does
    // it do more than that ? the documentation talks about pixel container being copied, check later what's this business.
    atrophyImage->Allocate();
    atrophyImage->FillBuffer(0.);

    // L : labelImage     A : atrophyImage      atrophy(I) : atrophy value for the corresponding label I read from the table.
    // if ( L(x) != 0 ) then A(x) = atrophy(L(x))

    typedef itk::ImageRegionConstIterator< LabelImageType > LabelIteratorType;
    typedef itk::ImageRegionIterator< AtrophyImageType > AtrophyIteratorType;

    LabelIteratorType labelIt(labelImage, labelImage->GetLargestPossibleRegion());
    AtrophyIteratorType atrophyIt(atrophyImage, atrophyImage->GetLargestPossibleRegion());

    labelIt.GoToBegin();
    atrophyIt.GoToBegin();
    while(!labelIt.IsAtEnd()){
	LabelImageType::PixelType currLabel = labelIt.Get();
	if(currLabel) {  //for all non-zero labels
	    for (LabelWithAtrophyIteratorType labelWithAtrophyIt = labelWithAtrophy.begin(); labelWithAtrophyIt != labelWithAtrophy.end(); ++labelWithAtrophyIt) {
		if (currLabel == labelWithAtrophyIt->first) {
		    atrophyIt.Set(labelWithAtrophyIt->second);
		    break;
		}
	    }
	}
	++labelIt;
	++atrophyIt;
    }

    // Write the atrophy output.
    typedef itk::ImageFileWriter< AtrophyImageType > AtrophyImageWriterType;
    AtrophyImageWriterType::Pointer atrophyImageWriter = AtrophyImageWriterType::New();
    atrophyImageWriter->SetFileName(atrophyMapFile);
    atrophyImageWriter->SetInput(atrophyImage);
    atrophyImageWriter->Update();
    return EXIT_SUCCESS;
}

