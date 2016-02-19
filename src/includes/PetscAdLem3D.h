#ifndef PETSCADLEM3D_H
#define PETSCADLEM3D_H

#include "AdLem3D.h"
//class AdLem3D; //FIXME: SEE IF WE CAN GET AWAY WITH DECLARING THIS SO THAT THE TEMPLATE PARA IS
// NOT NECESSARY


//base class solver for AdLem3D using Petsc.
#include<string>
#include<petscsys.h>
#include<petscdm.h>
#include<petscdmda.h>
#include<petscksp.h>


template <unsigned int DIM>
class PetscAdLem3D {
public:
    PetscAdLem3D(AdLem3D<DIM>*, bool set12pointStencilForDiv, const std::string&);
    virtual ~PetscAdLem3D();
    void setContextName(const std::string&);

    AdLem3D<DIM>* getProblemModel() const;
    PetscBool isDiv12pointStencil();

    double getRhsAt(unsigned int pos[3], unsigned int component);
//    double getAtrophyGradientAt(unsigned int pos[3], unsigned int component);
    double getSolVelocityAt(unsigned int pos[3],unsigned int component);
    double getSolPressureAt(unsigned int pos[3]);
    double getDivergenceAt(unsigned int pos[3]);
    PetscErrorCode writeResidual(std::string resultPath);

protected:
    AdLem3D<DIM>*            mProblemModel;
    std::string         mContextDesc;

    PetscBool	mIsDiv12pointStencil;
    PetscBool   mIsMuConstant;
    PetscBool   mSolAllocated;
    PetscBool   mRhsAllocated;

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


#include "PetscAdLem3D.hxx"
#endif // PETSCADLEM3D_H
