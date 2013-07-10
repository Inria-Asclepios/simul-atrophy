#include"PetscAdLem2D.hxx"
PetscAdLem2D::PetscAdLem2D(AdLem2D *model, const std::string& description):
    mProblemModel(model), mContextDesc(description){}

PetscAdLem2D::~PetscAdLem2D(){

}

void PetscAdLem2D::setContextName(const std::string & desc){
    mContextDesc = std::string(desc);
}

PetscErrorCode PetscAdLem2D::ComputeRho(PetscInt i, PetscInt j, PetscInt m,
                          PetscInt n, PetscReal *rho){
    *rho = mProblemModel->SetAtrophy();
    PetscFunctionReturn(0);
}

AdLem2D* PetscAdLem2D::getProblemModel() const  {
    return mProblemModel;
}



