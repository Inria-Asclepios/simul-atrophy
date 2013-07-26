#ifndef PETSCADLEMTARAS3D_HXX
#define PETSCADLEMTARAS3D_HXX

#include"AdLem3D.hxx"
#include"PetscAdLem3D.hxx"

#include<petscsys.h>
#include<petscdm.h>
#include<petscksp.h>
#include<petscdmda.h>

class PetscAdLemTaras3D : public PetscAdLem3D {
public:
    typedef struct {
        PetscScalar vx, vy, vz, p;
    } Field;
    PetscAdLemTaras3D(AdLem3D* );
    virtual ~PetscAdLemTaras3D();
    PetscReal muC(PetscInt x, PetscInt y, PetscInt z);
    PetscReal muXy(PetscInt x, PetscInt y, PetscInt z);
    PetscReal muXz(PetscInt x, PetscInt y, PetscInt z);
    PetscReal muYz(PetscInt x, PetscInt y, PetscInt z);

    PetscReal lambdaC(PetscInt x, PetscInt y, PetscInt z);
    PetscReal lambdaXy(PetscInt x, PetscInt y, PetscInt z);
    PetscReal lambdaXz(PetscInt x, PetscInt y, PetscInt z);
    PetscReal lambdaYz(PetscInt x, PetscInt y, PetscInt z);

    PetscReal aC(PetscInt x, PetscInt y, PetscInt z);
    PetscReal getP0Cell();
    PetscErrorCode solveModel(bool fileToMatlab, const std::string& filename);
    static PetscErrorCode computeMatrixTaras3D(KSP, Mat, Mat, MatStructure*, void*);
    static PetscErrorCode computeRHSTaras3D(KSP, Vec, void*);

protected:
    DM mDa;
    KSP mKsp;

    PetscReal dataCenterAt(std::string dType, PetscInt x, PetscInt y, PetscInt z);
    PetscReal dataXyAt(std::string dType, PetscInt x, PetscInt y, PetscInt z);
    PetscReal dataXzAt(std::string dType, PetscInt x, PetscInt y, PetscInt z);
    PetscReal dataYzAt(std::string dType, PetscInt x, PetscInt y, PetscInt z);

};

#endif // PETSCADLEMTARAS3D_HXX
