#ifndef ADLEM3D_HXX
#define ADLEM3D_HXX
#include<string>

#include <iostream>
#include <itkImage.hxx>
#include <itkDiffusionTensor3D.hxx>

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
We will use lambda as a symmetric +ve definte tensor image with
tensor at each voxel that describes anisotropy of the computational
domain.
*/

class PetscAdLemTaras3D;
class AdLem3D{
public:
    enum bcType {
        DIRICHLET, NEUMANN
    };

    typedef typename itk::Image<double, 3>                  ScalarImageType;
    typedef typename itk::Image<int, 3>                     IntegerImageType;
    typedef typename itk::Image<itk::Vector<double,3>, 3>   VectorImageType;
    typedef typename itk::Image<itk::DiffusionTensor3D< double >, 3>  TensorImageType;

    typedef typename itk::ImageFileReader<IntegerImageType> IntegerImageReaderType;
    typedef typename itk::ImageFileReader<ScalarImageType>  ScalarImageReaderType;
    typedef typename itk::ImageFileReader<TensorImageType>  TensorImageReaderType;
    typedef typename itk::ImageFileWriter<ScalarImageType>  ScalarImageWriterType;
    typedef typename itk::ImageFileWriter<VectorImageType>  VectorImageWriterType;
    typedef typename itk::ImageFileWriter<TensorImageType>  TensorImageWriterType;

    AdLem3D();
    ~AdLem3D();

    int getXnum() const;
    int getYnum() const;
    int getZnum() const;


    double getXspacing() const;
    double getYspacing() const;
    double getZspacing() const;

    //--**********Boundary condition related functions**************//
    bcType getBcType() const;
    void setWallVelocities(std::vector<double>& wallVelocities);
    void getWallVelocities(std::vector<double>& wallVelocities); //copies mWallVelocities content.

    //--***********Model parameters related functions***************//

    // Currently, mu for CSF is asked and the mu values for GM/WM are inferred from the ratio
    // given as input. GM & WM both have same mu value.
    // For lambda, in scalar case is same as mu; set using this function.
    // For tensor value, provide appropriate boolean argument to this function and
    // use setLambda function to provide a tensor image.
    void setLameParameters(bool isMuConstant, bool useTensorLambda,
                           double muCsf = 1, double muRatio = 1,
                           double lambdaCsf = 1, double lambdaRatio = 1, std::string lambdaImageFile = "");
    bool isMuConstant() const;
    bool isLambdaTensor() const;

    void setBrainMask(std::string maskImageFile, int relaxIcLabel,int relaxIcPressureCoeff, bool setSkullVelToZero, int skullLabel);
    void setBrainMask(ScalarImageType::Pointer brainMask, int relaxIcLabel, int relaxIcPressureCoeff, bool setSkullVelToZero, int skullLabel);
    ScalarImageType::Pointer getBrainMaskImage();
    void writeBrainMaskToFile(std::string fileName);

    int getRelaxIcPressureCoeff();

    //string should be either of "mu", "lambda" or "atrophy".
    //(Mi,Mj) is the element position of a tensor. Currently lambda can be a tensor.
    // Top left position is: (Mi,Mj)=(0,0).
    double dataAt(std::string dType, int x, int y, int z, unsigned int Mi, unsigned int Mj);
    double brainMaskAt(int x, int y, int z) const; //returns int unlike dataAt()
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
    void solveModel(bool operatorChanged = false);
    VectorImageType::Pointer getVelocityImage();
    ScalarImageType::Pointer getPressureImage();
    ScalarImageType::Pointer getDivergenceImage();
    VectorImageType::Pointer getForceImage();

    void writeSolutionForMatlab(std::string resultsPath,
                       bool writeSystemSolution = false);
    void writeVelocityImage(std::string fileName);
    void writePressureImage(std::string fileName);
    void writeDivergenceImage(std::string fileName);
    void writeForceImage(std::string fileName);
    void writeResidual(std::string fileName);

    //Atrophy related functions
    void setAtrophy(ScalarImageType::Pointer inputAtrophy);
    void setAtrophy(std::string atrophyImageFile);
    bool isAtrophySumZero(double sumMaxValue);
    ScalarImageType::Pointer getAtrophyImage();

    //By default:
    // 1. Total sum of atrophy is not forced to be zero; use this option
    //if you want to enforce IC everywhere, that this is usually the case when you
    //are not using brainMask.
    // 2. Does not redistribute atrophy. Set this to true if you want to redistribute atrophy.
    // This is the case when you apply transformation to the atrophy map using linear interpolation but
    // NN interpolation for the brain mask.
    // The transforamation and interpolation can result in non-zero atrophy values in the CSF region.
    // This option will distribute uniformly the volume loss in these CSF regions to the nearest tissue volumes.
    // And then we can replace the atrophy values in the newly created CSF regions to zero.
    //By default, uses brainMask and sets the value maskValue to all the places
    //where brainMask has the value maskLabel.
    void modifyAtrophy(int maskLabel, double maskValue, bool redistributeAtrophy = false, bool makeSumZero = false);
    void scaleAtrophy(double factor);
    void writeAtrophyToFile(std::string fileName);

protected:
    ScalarImageType::Pointer    mBrainMask;         //Input segmentation.
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

    //parameters for Gray matter, White matter and Csf:
    TensorImageType::Pointer    mLambda;        //Input diffusion tensor image lambda.
    bool        mUseTensorLambda;
    double      mMuGm, mMuWm, mMuCsf;
    double      mLambdaGm, mLambdaWm, mLambdaCsf;
    bool        mIsMuConstant;

    //Boundary condition:
    AdLem3D::bcType             mBc;
    std::vector<double>         mWallVelocities;  //velocity components vx,vy,vz on S,W,N,E,F,B walls.

    // integer number to keep track of how many times the solver is solved.
    int                         mNumOfSolveCalls;


//    VectorImageType::Pointer    mAtrophyGradient;   //gradient of the atrophy computed within the solver
    VectorImageType::Pointer    mForce;             //Force term computed i.e. lambda*grad(a)
    ScalarImageType::Pointer    mPressure;          //Output pressure-map.
    VectorImageType::Pointer    mVelocity;          //Output velocity field.
    ScalarImageType::Pointer    mDivergence;        //Divergence map of the solution (computed by the solver)

    // Memory for different result images needed to be allocated only once.
    bool                        mVelocityAllocated;
    bool                        mPressureAllocated;
    bool                        mForceAllocated;
    bool                        mDivergenceAllocated;

    // Variables to track whether the result images have been updated after the most recent
    // solve.
    bool                        mVelocityLatest;
    bool                        mPressureLatest;
    bool                        mForceLatest;
    bool                        mDivergenceLatest;

    //Solver option
    PetscAdLemTaras3D   *mPetscSolverTaras;
    bool                mPetscSolverTarasUsed;

    // Hide this from the user, user must use setLameParameters interface to set lambda.
    // This method will be called by setLameParameters depending on user's argument to that
    // method.
    void setLambda(TensorImageType::Pointer inputLambda);

    // Internal data access utility methods.
    // The class that inherits will have to provide public interface to these by
    // properly adapting to the grid type they use.
    double muAt(int x, int y, int z) const;
    // Li and Lj are by default zero.
    // If mUseTensorLambda is false and Li and Lj have default argumetns,
    // then the method is implemented such a way that lambda is effectively a scalar.
    double lambdaAt(int x, int y, int z,
                    unsigned int Li = 0, unsigned int Lj = 0) const;
    double aAt(int x, int y, int z) const;

    void updateStateAfterSolveCall(); // update all the state and track variables that depend on
    // or should change once the solveModel() function is called.
    void createVelocityImage();
    void createPressureImage();
    void createForceImage();
    void createDivergenceImage();

    void updateVelocityImage();
    void updatePressureImage();
    void updateForceImage();
    void updateDivergenceImage();
};


//#include "itkImage.h"
//#include "itkImageRegionIterator.h"

//    typedef itk::SymmetricSecondRankTensor<double> PixelType;

//    typedef itk::Image<PixelType, 3> ImageType;

//    typedef itk::ImageFileReader<ImageType> ReaderType;

//    ReaderType:ointer reader = ReaderType::New();

//    reader->SetFileName(tensorImage_i);

//    reader->Update();

//    ImageType:ointer tensorImage = reader->GetOutput();

//    TensorImageType:ointer image_out = TensorImageType::New();

//    image_out->SetOrigin(tensorImage->GetOrigin());

//    image_out->SetDirection(tensorImage->GetDirection());

//    image_out->SetSpacing(tensorImage->GetSpacing());

//    image_out->SetRegions(tensorImage->GetLargestPossibleRegion());

//    image_out->Allocate();

//    image_out->FillBuffer(TensorType(0.0));

//    typedef itk::ImageRegionIterator<ImageType> ImageIterator;

//    typedef itk::ImageRegionIterator<TensorImageType> TensorIterator;

//    TensorIterator itTTK(image_out, image_out->GetLargestPossibleRegion());

//    ImageIterator itImage(tensorImage, tensorImage->GetLargestPossibleRegion());

//    ImageType::IndexType testIndex;

//    testIndex[0]=44; testIndex[1]=56; testIndex[2]=30;

//    for (itTTK.GoToBegin(), itImage.GoToBegin(); !itTTK.IsAtEnd(), !itImage.IsAtEnd();

//         ++itTTK, ++itImage)

//    {

//         double dxx, dxy, dxz, dyy, dyz, dzz;

//         PixelType tempPixel;

//         tempPixel = itImage.Get();

//         dxx = tempPixel.GetNthComponent(0);

//         dxy = tempPixel.GetNthComponent(1);

//         dyy = tempPixel.GetNthComponent(2);

//         dxz = tempPixel.GetNthComponent(3);

//         dyz = tempPixel.GetNthComponent(4);

//         dzz = tempPixel.GetNthComponent(5);

//         TensorType tempTensor;

//         tempTensor.SetNthComponent(0,dxx);

//         tempTensor.SetNthComponent(1,dxy);

//         tempTensor.SetNthComponent(2,dxz);

//         tempTensor.SetNthComponent(3,dyy);

//         tempTensor.SetNthComponent(4,dyz);

//         tempTensor.SetNthComponent(5,dzz);

//         itTTK.Set(tempTensor);

//    }

//    TensorISetFileName(tensorImage_o);

//    TensorISetInput(image_out);

//    TensorIUpdate();

//    TensorIWrite();

//    TensorIUpdate();

//    return 0;

//    }

#endif // ADLEM3D_HXX
