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
  div(v)          + c*p          = -a

where,
mu and lambda are the Lame coefficients,
v is the velcoity/displacement for small time step,
a is the atrophy,
c is a function of position. Options:
c=0 everywhere, to enforce incompressibility everywhere.
c=non-zero, 1/lambda ? pressureMasscoeff; at places where one wants to release the incompressibility.
*/

class PetscAdLemTaras3D;
class AdLem3D{
public:
    enum bcType {
        DIRICHLET, NEUMANN
    };

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

    //--**********Boundary condition related functions**************//
    bcType getBcType() const;
    void setWallVelocities(std::vector<double>& wallVelocities);
    void getWallVelocities(std::vector<double>& wallVelocities); //copies mWallVelocities content.


    //--***********Model parameters related functions***************//
    bool isMuConstant() const;

    void setLameParameters(double muCsf, double lambdaCsf,
                           bool isMuConstant = true,
                           double muRatio = 1, double lambdaRatio = 1);
    void setPressureMassCoeffCsf(int coeff);    //value of c.
    int getPressureMassCoeffCsf();

    //string should be either of "mu", "lambda" or "atrophy"
    double dataAt(std::string dType, int x, int y, int z);
    int brainMaskAt(int x, int y, int z) const; //returns int unlike dataAt()

    void setBrainMask(std::string maskImageFile);

    //solver related functions
    void setDomainRegion(unsigned int origin[3], unsigned int size[3],
                         bool fullImageSize = false);
    void solveModel();
    void writeSolution(std::string resultsPath, bool inMatlabFormat = false,
                       bool inMatlabFormatSystemSolution = false);
    void writeResidual(std::string resultsPath);

    //Atrophy related functions
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

    //Boundary condition:
    AdLem3D::bcType             mBc;
    std::vector<double> mWallVelocities;  //velocity components vx,vy,vz on S,W,N,E,F,B walls.

    //parameters for Gray matter, White matter and Csf:
    double      mMuGm, mMuWm, mMuCsf;
    double      mLambdaGm, mLambdaWm, mLambdaCsf;
    bool        mIsMuConstant;

    //pressure coefficient at CSF: 0 => strictly follow incompressibilty constraint.
    //big number => release IC and allow pressure to vary.
    //by default, it's set to zero.
    int      mPressureMassCoeffCsf;

    //Solver option
    PetscAdLemTaras3D   *mPetscSolverTaras;
    bool                mPetscSolverTarasUsed;

    double muAt(int x, int y, int z) const;
    double lambdaAt(int x, int y, int z) const;
    double aAt(int x, int y, int z) const;
};

#endif // ADLEM3D_HXX
