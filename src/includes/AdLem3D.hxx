#ifndef ADLEM3D_HXX
#define ADLEM3D_HXX
#include<string>

#include <iostream>
#include <itkImage.hxx>
#include <itkImageFileReader.hxx>
#include <itkImageFileWriter.hxx>

//#include"PetscAdLem3D.hxx"
//#include<petscsys.h>

/* Linear Elastic Model 3D for AD deformation. This model contains parameters and
  inputs of the following system:
  div(mu grad(v)) - grad(p) = (mu + lambda)grad(a)
  div(v)                    = -a

where,
mu and lambda are the Lame coefficients,
v is the velcoity/displacement for small time step,
and a is the atrophy
*/

class PetscAdLemTaras3D;
class AdLem3D{
public:
    enum bcType {
        DIRICHLET, NEUMANN
    };

//    typedef typename itk::Image<int, 3>         IntegerImageType;
    typedef typename itk::Image<double, 3>      ScalarImageType;
    typedef typename itk::Image<itk::Vector<double,3>, 3>   VectorImageType;

    typedef typename itk::ImageFileReader<ScalarImageType>  ScalarReaderType;
    typedef typename itk::ImageFileWriter<ScalarImageType>  ScalarWriterType;
    typedef typename itk::ImageFileWriter<VectorImageType>  VectorWriterType;

    AdLem3D();
    ~AdLem3D();

    int getXnum() const;
    int getYnum() const;
    int getZnum() const;
    bcType getBcType() const;

    //string should be either of "mu", "lambda" or "atrophy"
    long double dataAt(std::string dType, int x, int y, int z);

    void setBrainMask(std::string maskImageFile);
    void setLameParameters(double muCsf, double lambdaCsf,
                           double muRatio, double lambdaRatio);
    //void setOrigin() //Get Origin from atrophy or mask Image!
    void setDomainRegion(unsigned int origin[3], unsigned int size[3],
                         bool fullImageSize = false);

    void solveModel();
    void writeSolution(std::string resultsPath);

    void createAtrophy(unsigned int size[3]);
    void setAtrophy(ScalarImageType::Pointer inputAtrophy);
    void setAtrophy(std::string atrophyImageFile);
    bool isAtrophyValid(double sumMaxValue);
    void modifyAtrophy();
    void scaleAtorphy(double factor);
    void writeAtrophyToFile(std::string fileName);

protected:
    ScalarImageType::Pointer    mAtrophy;           //Input atrophy-map.
    ScalarImageType::Pointer    mBrainMask;         //Input segmentation.
    ScalarImageType::Pointer    mPressure;          //Output pressure-map.
    VectorImageType::Pointer    mVelocity;          //Output velocity field.

    ScalarImageType::RegionType mDomainRegion;

    AdLem3D::bcType             mBc;

    //parameters for Gray matter, White matter and Csf:
    long double mMuGm, mMuWm, mMuCsf;
    long double mLambdaGm, mLambdaWm, mLambdaCsf;

    PetscAdLemTaras3D *mPetscSolverTaras;
    bool                mPetscSolverTarasUsed;

    long double muAt(int x, int y, int z) const;
    long double lambdaAt(int x, int y, int z) const;
    long double aAt(int x, int y, int z) const;
};

#endif // ADLEM3D_HXX
