#include "PetscAdLem3D.hxx"

PetscAdLem3D::PetscAdLem3D(AdLem3D *model, const std::string & description):
    mProblemModel(model), mContextDesc(description) {}

PetscAdLem3D::~PetscAdLem3D()
{
}

void PetscAdLem3D::setContextName(const std::string & desc)
{
    mContextDesc = std::string(desc);
}

AdLem3D* PetscAdLem3D::getProblemModel() const
{
    return mProblemModel;
}

