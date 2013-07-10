#ifndef PETSCADLEMTARAS2D_HXX
#define PETSCADLEMTARAS2D_HXX

#include"AdLem2D.hxx"
#include"PetscAdLem2D.hxx"
#include<petscsys.h>
#include<petscdm.h>
#include<petscksp.h>
#include<petscdmda.h>

class PetscAdLemTaras2D : public PetscAdLem2D{
public:
    PetscAdLemTaras2D(AdLem2D*);
    virtual ~PetscAdLemTaras2D();
    PetscReal getNu() const;
    PetscErrorCode solveModel();

protected:
    DM mDaPara;
    DM mDaSol;
    //System Matrix
//    Mat mA;

    //Linear solver context
    KSP mKsp;

    //local grid info:
//    PetscInt mXs, mYs, mXm, mYm;
    //local info:
//    DMDALocalInfo mLocalInfoPara;
//    DMDALocalInfo mLocalInfoSol;
};


PetscErrorCode computeMatrixTaras2D(KSP, Mat, Mat, MatStructure*, void*);
PetscErrorCode computeRHSTaras2D(KSP, Vec, void*);
#endif // PETSCADLEMTARAS2D_HXX
