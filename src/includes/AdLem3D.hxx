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
    typedef typename itk::Image<int, 3>         IntegerImageType;
    typedef typename itk::Image<itk::Vector<double,3>, 3>   VectorImageType;

    typedef typename itk::ImageFileReader<IntegerImageType>  IntegerImageReaderType;
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
    int getRelaxIcPressureCoeff();

    //string should be either of "mu", "lambda" or "atrophy"
    double dataAt(std::string dType, int x, int y, int z);
    int brainMaskAt(int x, int y, int z) const; //returns int unlike dataAt()
    int getRelaxIcLabel() const;  //Return the label present in BrainMask where the IC is to be relaxed,
    //by setting pressureMassCoeff.

    void setBrainMask(std::string maskImageFile, int relaxIcLabel,int relaxIcPressureCoeff);
    //relaxIc should almost always be true because the utility of Brain mask is only
    //when we want to release IC on certain places!
    //IMPORTANT: The maskImage must have at least some interior voxels where
    //the values are equal to relaxIcLabel IF mPressureMassCoeffCsf is non-zero;
    //This is because the Taras solver will not be able to produce appropriate null-space.

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

    //By default, the total sum of atrophy is not forced to be zero; use this option
    //if you want to enforce IC everywhere, that this is usually the case when you
    //are not using brainMask.
    //By default, uses brainMask.
    void modifyAtrophy(int maskLabel, double maskValue, bool makeSumZero = false);
    void scaleAtorphy(double factor);
    void writeAtrophyToFile(std::string fileName);

protected:
    ScalarImageType::Pointer    mAtrophy;           //Input atrophy-map.
    IntegerImageType::Pointer   mBrainMask;         //Input segmentation.
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

    //strictly follow incompressibilty constraint(IC) only at non-CSF parts.
    //big number => release IC and allow pressure to vary.
    //by default, it's set to zero.
    int      mRelaxIcPressureCoeff;  //Non-zero => IC is relaxed at mRelaxIcLabel values of the brain mask.
    int      mRelaxIcLabel;     //Label value in the brainMask where the IC is to be relaxed.

    //Solver option
    PetscAdLemTaras3D   *mPetscSolverTaras;
    bool                mPetscSolverTarasUsed;

    double muAt(int x, int y, int z) const;
    double lambdaAt(int x, int y, int z) const;
    double aAt(int x, int y, int z) const;
};

#endif // ADLEM3D_HXX
