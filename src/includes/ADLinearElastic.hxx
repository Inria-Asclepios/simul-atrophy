#ifndef ADLINEARELASTIC_HXX
#define ADLINEARELASTIC_HXX

#include<petscsys.h>
#include<petscdmda.h>
#include<petscksp.h>

class ADLinearElastic{
public:
    ADLinearElastic(PetscInt, PetscInt, PetscReal);
//    InitializeLameParameters();
    enum bcType{
        DIRICHLET, NEUMANN
    };
    PetscErrorCode ComputeRho(PetscInt, PetscInt, PetscInt, PetscInt, PetscReal *);
    const inline PetscReal getHx(){ return hx;}
    const inline PetscReal getHy(){ return hy;}
    const inline bcType getBcType(){ return bc;}
    const inline PetscReal getNu(){ return nu; }


protected:
    Vec mu_c, mu_n, lambda_c, lambda_n; //lame parameters at cell center and at nodes.
    DM da_para, da_sols; //one for parameters with dof 1, other for the variables to be solved.
    PetscInt m,n; //grid size
    PetscReal hx, hy;
    ADLinearElastic::bcType bc;
    PetscReal centerRho,nu;
    Vec a; //atrophy rate

};

//Would have been better to put it inside the class, but with function
//pointers required by PetscCalls, better to make it non-member!
PetscErrorCode ComputeMatrix(KSP, Mat, Mat, MatStructure*, void*);
PetscErrorCode ComputeRHS(KSP,Vec, void*);

#endif // ADLINEARELASTIC_HXX
