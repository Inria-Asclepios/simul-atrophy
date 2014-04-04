#ifndef ADLEM3D_HXX
#define ADLEM3D_HXX
#include<string>

#include <iostream>
#include <itkImage.hxx>
#include <itkImageFileReader.hxx>
#include <itkImageFileWriter.hxx>

//#include"PetscAdLem3D.hxx"
//#include<petscsys.h>

/* Linear Elastic Model/Viscous fluid model? 3D for AD deformation. This model contains parameters and
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
    void setLameParameters(double muCsf, double lambdaCsf,
                           bool isMuConstant = true,
                           double muRatio = 1, double lambdaRatio = 1);
    bool isMuConstant() const;
    void setBrainMask(std::string maskImageFile, int relaxIcLabel,int relaxIcPressureCoeff, bool setSkullVelToZero, int skullLabel);
    void setBrainMask(IntegerImageType::Pointer brainMask, int relaxIcLabel, int relaxIcPressureCoeff, bool setSkullVelToZero, int skullLabel);
    int getRelaxIcPressureCoeff();

    //string should be either of "mu", "lambda" or "atrophy"
    double dataAt(std::string dType, int x, int y, int z);
    int brainMaskAt(int x, int y, int z) const; //returns int unlike dataAt()
    int getRelaxIcLabel() const;  //Return the label present in BrainMask where the IC is to be relaxed.
    //Relaxes IC on those region where maskImageFile has values relaxIcLabel.
    //Normally, relaxIcPressureCoeff should be non-zero.
    //If you do want to relax IC anywhere, set relaxIcPressureCoeff to zero.

    bool shouldSkullVelSetToZero();
    int getSkullLabel();

    //BrainMask must already be set before using this, because for full size region it
    //gets the region info from the BrainMask image!!
    void setDomainRegionFullImage();
    void setDomainRegion(unsigned int origin[3], unsigned int size[3]);

    //solver related functions
    void solveModel();
    VectorImageType::Pointer getVelocityImage();
    ScalarImageType::Pointer getPressureImage();

    void writeSolution(std::string resultsPath, bool inMatlabFormat = false,
                       bool inMatlabFormatSystemSolution = false);
    void writeResidual(std::string resultsPath);

    //Atrophy related functions
    void createAtrophy(unsigned int size[3]);   //This should be later moved to another class!
    void setAtrophy(ScalarImageType::Pointer inputAtrophy);
    void setAtrophy(std::string atrophyImageFile);
    bool isAtrophyValid(double sumMaxValue);
    ScalarImageType::Pointer getAtrophyImage();

    //By default, the total sum of atrophy is not forced to be zero; use this option
    //if you want to enforce IC everywhere, that this is usually the case when you
    //are not using brainMask.
    //By default, uses brainMask and sets the value maskValue to all the places
    //where brainMask has the value maskLabel.
    void modifyAtrophy(int maskLabel, double maskValue, bool makeSumZero = false);
    void scaleAtrophy(double factor);
    void writeAtrophyToFile(std::string fileName);

protected:
    IntegerImageType::Pointer   mBrainMask;         //Input segmentation.
    bool                        mIsBrainMaskSet;
    //strictly follow incompressibilty constraint(IC) only at non-CSF parts.
    //big number => release IC and allow pressure to vary.
    //the constructor sets it to zero.
    int      mRelaxIcPressureCoeff;  //Non-zero => IC is relaxed at the regions where brainMask has the value mRelaxIcLabel.
    int      mRelaxIcLabel;     //Label value in the brainMask where the IC is to be relaxed.
    bool     mSetSkullVelToZero;        //if true, the velocity at all the places in bMask with label value of mDirchletBoundaryLabel will be set to zero.
    int      mSkullLabel;   //Label value in the brainMask where the velocity will be imposed to be zero

    ScalarImageType::RegionType mDomainRegion;

    ScalarImageType::Pointer    mAtrophy;           //Input atrophy-map.

    //Boundary condition:
    AdLem3D::bcType             mBc;
    std::vector<double>         mWallVelocities;  //velocity components vx,vy,vz on S,W,N,E,F,B walls.

    //parameters for Gray matter, White matter and Csf:
    double      mMuGm, mMuWm, mMuCsf;
    double      mLambdaGm, mLambdaWm, mLambdaCsf;
    bool        mIsMuConstant;

    ScalarImageType::Pointer    mPressure;          //Output pressure-map.
    VectorImageType::Pointer    mVelocity;          //Output velocity field.
    ScalarImageType::Pointer    mDivergence;        //Divergence map of the solution (computed by the solver)

    //Solver option
    PetscAdLemTaras3D   *mPetscSolverTaras;
    bool                mPetscSolverTarasUsed;

    double muAt(int x, int y, int z) const;
    double lambdaAt(int x, int y, int z) const;
    double aAt(int x, int y, int z) const;
    void createResultImages();
};

#endif // ADLEM3D_HXX
