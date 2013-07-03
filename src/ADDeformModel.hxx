//Base Class for all the AD deformation based models.

class ADDeformModel{
public:
  //public interface methods supporting the model. But yet to be decided: should this be an abstract base class or not?

protected:
  image *brain;  //longitudinal series of brain image input/observed
  image *atrophy;//longitudinal series of the atrophy map of the brain.
  image *brain_deformed; //longitudinal series of the brain created by the model.
  image *u;//longitudinal series of the displacement field.
};


