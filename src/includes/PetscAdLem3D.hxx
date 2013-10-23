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
    PetscAdLem3D(AdLem3D*, const std::string& );
    virtual ~PetscAdLem3D();
    void setContextName(const std::string&);

    AdLem3D* getProblemModel() const;

    double getSolVelocityAt(unsigned int pos[3],unsigned int component);
    double getSolPressureAt(unsigned int pos[3]);
    PetscErrorCode writeResidual(std::string resultPath);

protected:
    AdLem3D* mProblemModel;
    std::string mContextDesc;

    PetscBool mSolAllocated;

    //Global system i.e for Ax=b, with x containing all the velocity components and pressure nodes.
    Mat mA;
    Vec mX, mB;
    DM mDa;                     //combined DMDA for both pressure and velocity.
    KSP mKsp;

    Vec mXLocal;              //one that stores velocity vector to be sent outside.
    VecScatter mScatterCtx;     //context to scatter global mX to mXvLocal and mXpLocal.
    PetscScalar *mSol;         //Array to access mXLocal.

    PetscErrorCode getSolutionArray(); //point mSol to proper solution vector mXLocal.
};

#endif // PETSCADLEM3D_HXX
