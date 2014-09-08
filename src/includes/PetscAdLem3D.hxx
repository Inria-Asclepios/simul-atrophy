#ifndef PETSCADLEM3D_HXX
#define PETSCADLEM3D_HXX
//base class solver for AdLem3D using Petsc.
#include<string>
#include<petscsys.h>
#include<petscdm.h>
#include<petscdmda.h>
#include<petscksp.h>

//#include"AdLem3D.hxx"
class AdLem3D;

class PetscAdLem3D {
public:
    PetscAdLem3D(AdLem3D*, const std::string&);
    virtual ~PetscAdLem3D();
    void setContextName(const std::string&);

    AdLem3D* getProblemModel() const;

    double getRhsAt(unsigned int pos[3], unsigned int component);
//    double getAtrophyGradientAt(unsigned int pos[3], unsigned int component);
    double getSolVelocityAt(unsigned int pos[3],unsigned int component);
    double getSolPressureAt(unsigned int pos[3]);
    double getDivergenceAt(unsigned int pos[3]);
    PetscErrorCode writeResidual(std::string resultPath);

protected:
    AdLem3D*            mProblemModel;
    std::string         mContextDesc;

    PetscBool           mIsMuConstant;
    PetscBool           mSolAllocated;
    PetscBool           mRhsAllocated;

    //Global system i.e for Ax=b, with x containing all the velocity components and pressure nodes.
    Mat         mA;
    Vec         mX, mB;
    DM          mDa;                     //combined DMDA for both pressure and velocity.
    KSP         mKsp;                   //KSP object for the global system.

    Vec         mXLocal;              //one that stores velocity vector to be sent outside.
    VecScatter  mScatterCtx;   //context to scatter global mX to mXvLocal and mXpLocal.
    PetscScalar *mSolArray;         //Array to access mXLocal.

    Vec         mBLocal;        //one that stores rhs vector to be sent outside.
    VecScatter  mScatterRhsCtx; //context to scatter global mB to mBvLocal and mDivLocal.
    PetscScalar *mRhsArray;     //Array to access mBLocal.

    PetscErrorCode          getSolutionArray(); //point mSol to proper solution vector mXLocal.
    PetscErrorCode          getRhsArray(); //point mRhs to proper rhs vector mBLocal.
};

#endif // PETSCADLEM3D_HXX
