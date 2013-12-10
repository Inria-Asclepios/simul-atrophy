#ifndef PETSCADLEMTARAS3D_HXX
#define PETSCADLEMTARAS3D_HXX

#include"PetscAdLem3D.hxx"
//#include"AdLem3D.hxx"

#include<petscsys.h>
//#include<petscdm.h>
//#include<petscksp.h>
//#include<petscdmda.h>
#include<petscdmcomposite.h>

class PetscAdLemTaras3D : public PetscAdLem3D {
public:
    typedef struct {
        PetscScalar vx, vy, vz, p;
    } Field;
    PetscAdLemTaras3D(AdLem3D *model ,bool writeParaToFile);
    virtual ~PetscAdLemTaras3D();

    MatNullSpace getNullSpace();
    PetscReal muC(PetscInt x, PetscInt y, PetscInt z);
    PetscReal muXy(PetscInt x, PetscInt y, PetscInt z);
    PetscReal muXz(PetscInt x, PetscInt y, PetscInt z);
    PetscReal muYz(PetscInt x, PetscInt y, PetscInt z);

    PetscReal lambdaC(PetscInt x, PetscInt y, PetscInt z);
    PetscReal lambdaXy(PetscInt x, PetscInt y, PetscInt z);
    PetscReal lambdaXz(PetscInt x, PetscInt y, PetscInt z);
    PetscReal lambdaYz(PetscInt x, PetscInt y, PetscInt z);

    PetscReal aC(PetscInt x, PetscInt y, PetscInt z);
    PetscErrorCode solveModel();
    PetscErrorCode writeToMatFile(const std::string& fileName, bool writeA, const std::string& matFileName);
    static PetscErrorCode computeMatrixTaras3d(KSP, Mat, Mat, MatStructure*, void*);
    static PetscErrorCode computeRHSTaras3d(KSP, Vec, void*);
    static PetscErrorCode computeNullSpace(MatNullSpace, Vec, void*);
    static PetscErrorCode computeMatrixTaras3dConstantMu(KSP, Mat, Mat, MatStructure*, void*);
    static PetscErrorCode computeRHSTaras3dConstantMu(KSP, Vec, void*);

    PetscInt bMaskAt(PetscInt x, PetscInt y, PetscInt z);

protected:
    //Is it a constant viscosity case or discontinuous?
    PetscBool mIsMuConstant;
    //Write Parameters to file
    PetscBool mWriteParaToFile;
    Vec mAtrophy, mMu;
    PetscBool mParaVecsCreated; //bool that is by default false, but
    //should be set to true by any non-static method that creates these vectors.
    //so that destructor will destroy them.
    PC mPc;
    MatNullSpace mNullSpace;    //Null space for the global system.
    Vec mNullBasis;             //Null basis for the global system.

    Vec mXv, mBv;               //vectors for the velocity field.

    Vec mXp, mBp;               //vectors for the pressure field.
    DM mDaP;                    //DMDA for pressure variable.
    Mat mPcForSc;               //Preconditioner matrix for the Schur Complement.
    MatNullSpace mNullSpaceP;   //Null space for the pressure field.
    Vec mNullBasisP;            //Null basis for the pressure field.




    void setNullSpace();
    PetscErrorCode createParaVectors();

    //Create preconditioner matrix for the Schur's Complement outer solve.
    void createPcForSc();

    PetscReal dataCenterAt(std::string dType, PetscInt x, PetscInt y, PetscInt z);
    PetscReal dataXyAt(std::string dType, PetscInt x, PetscInt y, PetscInt z);
    PetscReal dataXzAt(std::string dType, PetscInt x, PetscInt y, PetscInt z);
    PetscReal dataYzAt(std::string dType, PetscInt x, PetscInt y, PetscInt z);

};

#endif // PETSCADLEMTARAS3D_HXX
