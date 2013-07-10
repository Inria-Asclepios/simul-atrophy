//Linear Elastic 2D model of the brain deformation due to AD.
//Derived from the ADDeformModel.

class ADLinearElastic2D : public ADDeformModel{
public:
  PetscSolver SolverMethod1; //method described in Taras Gerya book.
  PetscSolver SolverMethod2; //Augmented Lagrangian method.
protected:
  image mu, lambda;  //Lame parameters.
  image p; //pressure or the lagrange multiplier in energy minimization.
  image compute_domain;//cropped image from the original one to contain an image containing up to the skull.
};
