#ifndef POISSONPROBLEM_H
#define POISSONPROBLEM_H

//Bishesh Khanal, Asclepios INRIA Sophia Antipolis

#include "PoissonSolver.hxx"

#include <iostream>
#include <itkImage.hxx>
#include <itkImageFileReader.hxx>
#include <itkImageFileWriter.hxx>

template <class PixelType>
class PoissonProblem
{
public:
    typedef typename itk::Image<PixelType,3>     ScalarImageType;
    typedef typename ScalarImageType::Pointer    ScalarImagePointerType;
    typedef typename ScalarImageType::PixelType  ScalarImagePixelType;
    typedef typename ScalarImageType::IndexType  ScalarImageIndexType;
    typedef typename ScalarImageType::SizeType   ScalarImageSizeType;
    typedef typename ScalarImageType::RegionType ScalarImageRegionType;

    typedef typename itk::ImageFileReader<ScalarImageType>  ScalarReaderType;
    typedef typename ScalarReaderType::Pointer              ScalarReaderPointerType;
    typedef typename itk::ImageFileWriter<ScalarImageType>  ScalarWriterType;
    typedef typename ScalarWriterType::Pointer              ScalarWriterPointerType;

    PoissonProblem(std::string problemName): mProblemName(problemName)
    {

        mIsRhsSet = false;
        mPoissonSolverUsed = false;

    }

    ~PoissonProblem()
    {
        if(mPoissonSolverUsed) delete mPoissonSolver;
    }

    void setProblemName(std::string problemName)
    {
        mProblemName = problemName;
    }

    std::string getProblemName()
    {
        return mProblemName;
    }

    void setDomainRegionFullImage() //Takes the size from Rhs image.
    {
        if(!mIsRhsSet) {
            std::cerr<<"must first set the Rhs image to obtain the largest possible region"<<std::endl;
        }
        mDomainRegion = mRhs->GetLargestPossibleRegion();
    }

    void setDomainRegion(unsigned int origin[3], unsigned int size[3])
    {
        ScalarImageIndexType domainOrigin;
        ScalarImageSizeType domainSize;
        for(int i=0;i<3;++i) {     //ADD error-guard to ensure that
            domainSize.SetElement(i,size[i]);     //size!!
            domainOrigin.SetElement(i,origin[i]); //provided is less or equal to the input
        }
        mDomainRegion.SetIndex(domainOrigin);
        mDomainRegion.SetSize(domainSize);
        //    std::cout<<"domain size: "<<mDomainRegion.GetSize()<<std::endl;
    }

    unsigned int getSize(unsigned int axis) const
    {
        return mDomainRegion.GetSize()[axis];
    }

    void setCoeff(int constCoeffValue, std::string coeffImageFile = "")
    {
        if(constCoeffValue == 0) { //0 means, use image to read the coeff. file.
            ScalarReaderPointerType scalarReader = ScalarReaderType::New();
            scalarReader->SetFileName(coeffImageFile);
            scalarReader->Update();
            setCoeff(scalarReader->GetOutput());
        } else
            mConstCoeffValue = constCoeffValue;
    }

    void setCoeff(ScalarImagePointerType Image)
    {
        mCoeff=Image;
    }

    PixelType getCoeffAtPosition(unsigned int x, unsigned  int y,unsigned  int z)
    {
        if(mConstCoeffValue != 0)
            return (PixelType)mConstCoeffValue;
        else {
            ScalarImageIndexType pos;
            pos.SetElement(0, mDomainRegion.GetIndex()[0] +  x);
            pos.SetElement(1, mDomainRegion.GetIndex()[1] +  y);
            pos.SetElement(2, mDomainRegion.GetIndex()[2] +  z);
            return(mCoeff->GetPixel(pos));
        }
    }

    void setRhs(std::string rhsImageFile)
    {
        ScalarReaderPointerType scalarReader = ScalarReaderType::New();
        scalarReader->SetFileName(rhsImageFile);
        scalarReader->Update();
        setRhs(scalarReader->GetOutput());
    }

    void setRhs(ScalarImagePointerType Image)
    {
        mIsRhsSet = true;
        mRhs=Image;
    }

    PixelType getRhsAtPosition(unsigned int x, unsigned int y, unsigned int z)
    {
        ScalarImageIndexType pos;
        pos.SetElement(0, mDomainRegion.GetIndex()[0] +  x);
        pos.SetElement(1, mDomainRegion.GetIndex()[1] +  y);
        pos.SetElement(2, mDomainRegion.GetIndex()[2] +  z);
        return(mRhs->GetPixel(pos));
        //      return 0.3;
    }

    void solveModel(bool computeResidual = true) //if you are sure you don't need residual.
    {
        if(!mPoissonSolverUsed) {
            mPoissonSolverUsed = true;
            mPoissonSolver = new PoissonSolver(this);
        }
        mPoissonSolver->solve();
        mSolution = ScalarImageType::New();
        createImageOf("mSolution");
        if (computeResidual) {
            mResidual = ScalarImageType::New();
            createImageOf("mResidual");
        }
    }

    ScalarImagePointerType getSolution()
    {
        return mSolution;
    }

    ScalarImagePointerType getResidual()
    {
        return mResidual;
    }

    void writeSolution(std::string resultsPath)
    {
        mPoissonSolver->writeToFile("vtk",resultsPath+mProblemName,true);

        //        writeAsImage("solution",resultsPath+mProblemName+"_sol.mha");
        //        writeAsImage("residual",resultsPath+mProblemName+"_res.mha");
    }


protected:

private:
    std::string                            mProblemName;
    int                                    mConstCoeffValue;
    //Set to zero if mCoeff is to be used by providing an image.
    //If set to non-zero, no need to use mCoeff.

    ScalarImagePointerType                 mCoeff;
    ScalarImagePointerType                 mRhs;
    bool                                   mIsRhsSet;

    ScalarImagePointerType                 mSolution;
    ScalarImagePointerType                 mResidual;

    ScalarImageRegionType                   mDomainRegion;

    PoissonSolver   *mPoissonSolver;
    bool            mPoissonSolverUsed;

    void createImageOf(std::string memberName)
    {
        ScalarImageRegionType   domainRegion;
        itk::Index<3>           outputImageStart;       //should be 0!
        for (unsigned int i=0; i<3; ++i) outputImageStart.SetElement(i,0);
        domainRegion.SetIndex(outputImageStart);
        domainRegion.SetSize(mDomainRegion.GetSize());

        ScalarImagePointerType imagePtr;
        if(memberName.compare("mSolution") == 0) {
            imagePtr = mSolution;
        } else if(memberName.compare("mResidual") == 0) {
            imagePtr = mResidual;
        } else {
            std::cout<<"internal error: unknown option, assuming to be solution"<<std::endl;
            imagePtr = mSolution;
        }

        imagePtr->SetRegions(domainRegion);
        imagePtr->SetOrigin(mRhs->GetOrigin());
        imagePtr->SetSpacing(mRhs->GetSpacing());
        imagePtr->SetDirection(mRhs->GetDirection());
        imagePtr->Allocate();

        typedef itk::ImageRegionIterator<ScalarImageType> scalarIteratorType;
        scalarIteratorType imageIterator(imagePtr,imagePtr->GetLargestPossibleRegion());
        ScalarImagePixelType imagePixel;

        unsigned int k = 0;
        imageIterator.GoToBegin();
        while(k<mDomainRegion.GetSize()[2] && !imageIterator.IsAtEnd()) {
            unsigned int j = 0;
            while(j<mDomainRegion.GetSize()[1] && !imageIterator.IsAtEnd()) {
                unsigned int i = 0;
                while(i<mDomainRegion.GetSize()[0] && !imageIterator.IsAtEnd()) {

                    if(memberName.compare("mResidual") == 0)
                        imagePixel = mPoissonSolver->getResAtPosition(i,j,k);
                    else
                        imagePixel = mPoissonSolver->getSolAtPosition(i,j,k);

                    imageIterator.Set(imagePixel);
                    ++imageIterator;
                    ++i;
                }
                ++j;
            }
            ++k;
        }
    }

};
#endif // POISSONPROBLEM_H
