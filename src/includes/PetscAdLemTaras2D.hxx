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
    typedef struct {
        PetscScalar vx, vy, p;
    } Field;
    PetscAdLemTaras2D(AdLem2D*);
    virtual ~PetscAdLemTaras2D();
    PetscReal getNu() const;
    PetscReal muCenter(PetscInt x, PetscInt y);
    PetscReal muNode(PetscInt x, PetscInt y);
    PetscReal lambdaCenter(PetscInt x, PetscInt y);
    PetscReal lambdaNode(PetscInt x, PetscInt y);
    PetscReal aCenter(PetscInt x, PetscInt y);
    PetscReal aNode(PetscInt x, PetscInt y);
    PetscReal getP0Cell();
    PetscErrorCode solveModel(bool writeToMatlab, const std::string& filename);

protected:
    DM mDa;
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
    PetscReal dataCenterAt(std::string dType, PetscInt x, PetscInt y);
    PetscReal dataNodeAt(std::string dType, PetscInt x, PetscInt y);
};


PetscErrorCode computeMatrixTaras2D(KSP, Mat, Mat, MatStructure*, void*);
PetscErrorCode computeRHSTaras2D(KSP, Vec, void*);
#endif // PETSCADLEMTARAS2D_HXX
