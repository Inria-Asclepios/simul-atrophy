//Main file to test the Poisson Problem.
//Bishesh Khanal; Asclepios INRIA Sophia Antipolis

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

#include "PoissonProblem.hxx"
#include "PoissonProblem.cxx"
#include "HelmHoltzDecomposer.hxx"

#include <petsc.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkImageRegion.h>
#include <itkImageRegionIterator.h>

#include <itkLaplacianImageFilter.h>


static char help[] = "Solving Poisson Problem";
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    int         numOfRhs;
    PetscInitialize(&argc,&argv,(char*)0,help);
    {
        PetscErrorCode ierr;
        //---------------********* Options Parsing ***********-------------
        PetscBool optionFlag = PETSC_FALSE;
        char optionString[PETSC_MAX_PATH_LEN];
        //Get the configuration file: The first thing in the file is the number of lines to read.
        //The lines that follow have the following format:
        //inputvectorFilename resultPath/
        ierr = PetscOptionsGetString(NULL,"-configFile",optionString,PETSC_MAX_PATH_LEN,&optionFlag);CHKERRQ(ierr);
        if(!optionFlag) {
            std::cerr<<"Must provide a valid configFile"<<std::endl;
            return(EXIT_FAILURE);
        }
        std::ifstream                   configFile(optionString);

        //Extract each line
        std::string                     line;
        std::getline(configFile,line);
        std::istringstream              issFirstLine(line);
        issFirstLine>>numOfRhs;

        //Store all the filenames of the input vector field to be decomposed and the corresponding
        //result paths in a std::vectors.
        std::vector<std::string>        inVecFiles, resPaths;
        inVecFiles.reserve(numOfRhs);
        resPaths.reserve(numOfRhs);
        for(int i = 0; i<numOfRhs; ++i) {
            if(!std::getline(configFile,line)) {
                std::cerr<<"file does not contain enough number of lines specified on the first line";
                return(EXIT_FAILURE);
            }
            std::istringstream          iss(line);
            std::string                 rhs, resPath;
            iss>>rhs>>resPath;
            inVecFiles.push_back(rhs);
            resPaths.push_back(resPath);
        }

        //---------------------************ Helmholtz decomposition **************---------------------
        HelmHoltzDecomposer hh;
        unsigned int origin[3] = {20,20,20};
        //        unsigned int origin[3] = {0,0,0};
        unsigned int size[3] = {50, 50, 50};

        for(int i = 0; i<numOfRhs; ++ i) {
//            hh.setV(inVecFiles.at(i), false,origin,size);
                    hh.setV(inVecFiles.at(i),true,NULL,NULL);

            HelmHoltzDecomposer::ScalarWriterPointerType divWriter = HelmHoltzDecomposer::ScalarWriterType::New();
            divWriter->SetFileName(resPaths.at(i) + "hh_div.mha");
            divWriter->SetInput(hh.getDivV());
            divWriter->Update();

            hh.decompose();

            HelmHoltzDecomposer::ScalarWriterPointerType pressWriter = HelmHoltzDecomposer::ScalarWriterType::New();
            pressWriter->SetFileName(resPaths.at(i) + "hh_press.mha");
            pressWriter->SetInput(hh.getP());
            pressWriter->Update();

            HelmHoltzDecomposer::VectorWriterPointerType gradPWriter = HelmHoltzDecomposer::VectorWriterType::New();
            gradPWriter->SetFileName(resPaths.at(i) + "hh_gradP.mha");
            gradPWriter->SetInput(hh.getGradP());
            gradPWriter->Update();

            HelmHoltzDecomposer::VectorWriterPointerType curlAWriter = HelmHoltzDecomposer::VectorWriterType::New();
            curlAWriter->SetFileName(resPaths.at(i) + "hh_curlP.mha");
            curlAWriter->SetInput(hh.getCurlA());
            curlAWriter->Update();
        }
    }
    PetscFinalize();
    return 0;
}

