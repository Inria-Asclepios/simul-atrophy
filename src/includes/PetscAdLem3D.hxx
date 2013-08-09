#ifndef PETSCADLEM3D_HXX
#define PETSCADLEM3D_HXX
//base class solver for AdLem3D using Petsc.
#include<string>
#include<petscsys.h>
#include<petscvec.h>
#include<petscmat.h>

#include"AdLem3D.hxx"

class PetscAdLem3D {
public:
    PetscAdLem3D(AdLem3D*, const std::string& );
    virtual ~PetscAdLem3D();
    void setContextName(const std::string&);

    AdLem3D* getProblemModel() const;

protected:
    AdLem3D* mProblemModel;
    std::string mContextDesc;
};

#endif // PETSCADLEM3D_HXX
