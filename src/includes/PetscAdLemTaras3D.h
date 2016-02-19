#ifndef PETSCADLEMTARAS3D_H
#define PETSCADLEMTARAS3D_H

#include"AdLem3D.h"
#include"PetscAdLem3D.h"


#include<petscsys.h>
//#include<petscdm.h>
//#include<petscksp.h>
//#include<petscdmda.h>
#include<petscdmcomposite.h>

namespace PetscAdLemTaras3D_SolverOps
{

    const bool RELAX_IC_WITH_ZERO_ROWS = false;
    //const bool RELAX_IC_WITH_ZERO_ROWS = true;
}

//template <unsigned int DIM>
class PetscAdLemTaras3D : public PetscAdLem3D<3u> {
public:
    typedef struct {
        PetscScalar vx, vy, vz, p;
    } Field;
    PetscAdLemTaras3D(AdLem3D<3u> *model ,bool set12pointStencilForDiv, bool writeParaToFile);
    virtual ~PetscAdLemTaras3D();

    MatNullSpace getNullSpace();
    PetscReal muC(PetscInt x, PetscInt y, PetscInt z);
    PetscReal muXy(PetscInt x, PetscInt y, PetscInt z);
    PetscReal muXz(PetscInt x, PetscInt y, PetscInt z);
    PetscReal muYz(PetscInt x, PetscInt y, PetscInt z);

    // Support for lambda as a tensor and hence. (Mi,Mj) element position in the matrix.
    // Top left element is (0,0).
    PetscReal lambdaC(PetscInt x, PetscInt y, PetscInt z, PetscInt Mi, PetscInt Mj);
    // TODO: Current support for lambda as tensor is only for
    // the center of the cell.
    // Must define how to interpolate tensor to get the values
    // at the faces of the cell. Thus currently it is not
    // implemented.
    PetscReal lambdaXy(PetscInt x, PetscInt y, PetscInt z);
    PetscReal lambdaXz(PetscInt x, PetscInt y, PetscInt z);
    PetscReal lambdaYz(PetscInt x, PetscInt y, PetscInt z);

    PetscReal aC(PetscInt x, PetscInt y, PetscInt z);
    PetscErrorCode solveModel(bool operatorChanged);
    PetscErrorCode writeToMatFile(const std::string& fileName, bool writeA, const std::string& matFileName);
    static PetscErrorCode computeMatrixTaras3d(KSP, Mat, Mat, void*);
    static PetscErrorCode computeRHSTaras3d(KSP, Vec, void*);
    static PetscErrorCode computeNullSpace(MatNullSpace, Vec, void*);
    static PetscErrorCode computeMatrixTaras3dConstantMu(KSP, Mat, Mat, void*);
    static PetscErrorCode computeRHSTaras3dConstantMu(KSP, Vec, void*);

    PetscReal bMaskAt(PetscInt x, PetscInt y, PetscInt z);

protected:
    PetscInt	mNumOfSolveCalls;
    PetscBool   mParaVecsCreated;	//by default false but,
    //should be set to true by any non-static method
    //that creates them, so that destructor can destroy them.
    PetscBool   mWriteParaToFile;	//set by the constructor; used in the method writeToMatFile.
    Vec         mAtrophy, mMu;

    PetscBool       mOperatorComputed;
    PC              mPc;
    DM              mDaP;                    //DMDA for pressure variable.
    Mat             mPcForSc;               //Preconditioner matrix for the Schur Complement.

    Vec             mXv, mBv;               //vectors for the velocity field.
    Vec             mXp, mBp;               //vectors for the pressure field.

    PetscBool       mPressureNullspacePresent;    //true if constant pressure null space is present.
    MatNullSpace    mNullSpace;    //Null space for the global system.
    MatNullSpace    mNullSpaceP;   //Null space for the pressure field.
    Vec             mNullBasis;             //Null basis for the global system.
    Vec             mNullBasisP;            //Null basis for the pressure field.

    void            setNullSpace();
    PetscErrorCode  createParaVectors();
    void            createPcForSc();        //Preconditioner for Schur Complement.

    // Mi, Mj arguments added for matrix position for the dType that supports tensor.
    // Currently it is only supported for "lambda" dType. For others (Mi,Mj) value is not used.
    PetscReal dataCenterAt(std::string dType, PetscInt x, PetscInt y, PetscInt z, PetscInt Mi = 0, PetscInt Mj = 0);
    PetscReal dataXyAt(std::string dType, PetscInt x, PetscInt y, PetscInt z);
    PetscReal dataXzAt(std::string dType, PetscInt x, PetscInt y, PetscInt z);
    PetscReal dataYzAt(std::string dType, PetscInt x, PetscInt y, PetscInt z);

};

#include "PetscAdLemTaras3D.hxx"
#endif // PETSCADLEMTARAS3D_H
