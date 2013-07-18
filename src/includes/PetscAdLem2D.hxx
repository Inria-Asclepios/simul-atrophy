#ifndef PETSCADLEM2D_HXX
#define PETSCADLEM2D_HXX
//base class solver for ADLEM2D using Petsc.
#include<string>
#include<petscsys.h>
#include"AdLem2D.hxx"

class PetscAdLem2D{
public:
    PetscAdLem2D(AdLem2D*, const std::string&);
    virtual ~PetscAdLem2D();
    void setContextName(const std::string&);

    AdLem2D* getProblemModel() const;
protected:
    AdLem2D* mProblemModel; //problem_model
    std::string mContextDesc;
};




#endif // PETSCADLEM2D_HXX
