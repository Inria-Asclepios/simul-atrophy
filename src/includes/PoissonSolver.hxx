#ifndef POISSONSOLVER_H
#define POISSONSOLVER_H

#include<petsc.h>
#include<petscsys.h>
#include<petscksp.h>

template <class PixelType> class PoissonProblem;
//class PoissonProblem;

class PoissonSolver {
private:

    PoissonProblem<double>* mProblem;
    Mat         mA;
    Vec         mB, mX;         //Solving Ax = B

    PetscBool   mOperatorComputed;      //set to true the first time operator matrix created.

    Vec         mXLocal;        //one that stores vector to be sent outside.
    VecScatter  mScatterCtx;    //context to scatter global mX to mXLocal
    PetscScalar *mXArray;       //Array to access mX

    Vec         mRes, mResLocal;    //Global and local residual vectors!
    VecScatter  mResScatterCtx;
    PetscScalar *mResArray;
    PetscBool   mResComputed;

    DM          mDa;      //Manages the structured grid.
    KSP         mKsp;     //Manages the solving technique.

public:
    typedef struct {
        PetscScalar vx, vy, vz;
    }Field;
    PoissonSolver(PoissonProblem<double> *problem);
    ~PoissonSolver();

    PoissonProblem<double>* getProblem();
    PetscErrorCode solve();
    PetscErrorCode computeResidual();
    PetscErrorCode writeToFile(std::string format, std::string fileName, bool writeResidual);
    static PetscErrorCode computeMatrix(KSP, Mat, Mat, MatStructure*, void* );
    static PetscErrorCode computeRhs(KSP ksp,Vec b,void* ctx);
    double getSolAtPosition(unsigned int x, unsigned int y, unsigned int z);
    double getResAtPosition(unsigned int x, unsigned int y, unsigned int z);
};

#endif // POISSONSOLVER_H
