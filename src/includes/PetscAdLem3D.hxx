#ifndef PETSCADLEM3D_HXX
#define PETSCADLEM3D_HXX
#include "PetscAdLem3D.h"

#include"AdLem3D.hxx"

#undef __FUNCT__
#define __FUNCT__ "PetscAdLem3D"

template <unsigned int DIM>
PetscAdLem3D<DIM>::PetscAdLem3D(AdLem3D<DIM> *model, bool set12pointStencilForDiv, const std::string & description):
    mProblemModel(model),
    mIsDiv12pointStencil((PetscBool)set12pointStencilForDiv),
    mContextDesc(description)
{
    mSolAllocated = PETSC_FALSE;
    mRhsAllocated = PETSC_FALSE;
    mIsMuConstant = (PetscBool)model->isMuConstant();

}

template <unsigned int DIM>
PetscAdLem3D<DIM>::~PetscAdLem3D()
{
    PetscErrorCode ierr;

    if (mSolAllocated) {
        ierr = VecRestoreArray(mX,&mSolArray);CHKERRXX(ierr);
        ierr = VecScatterDestroy(&mScatterCtx);CHKERRXX(ierr);
        ierr = VecDestroy(&mXLocal);CHKERRXX(ierr);
    }
    if (mRhsAllocated) {
        ierr = VecRestoreArray(mB,&mRhsArray);CHKERRXX(ierr);
        ierr = VecScatterDestroy(&mScatterRhsCtx);CHKERRXX(ierr);
        ierr = VecDestroy(&mBLocal);CHKERRXX(ierr);
    }
    ierr = KSPDestroy(&mKsp);CHKERRXX(ierr);
    ierr = DMDestroy(&mDa);CHKERRXX(ierr);
}

template <unsigned int DIM>
void PetscAdLem3D<DIM>::setContextName(const std::string & desc)
{
    mContextDesc = std::string(desc);
}

template <unsigned int DIM>
AdLem3D<DIM>* PetscAdLem3D<DIM>::getProblemModel() const
{
    return mProblemModel;
}

#undef __FUNCT__
#define __FUNCT__ "isDiv12pointStencil"
template <unsigned int DIM>
PetscBool PetscAdLem3D<DIM>::isDiv12pointStencil()
{
    return mIsDiv12pointStencil;
}


#undef __FUNCT__
#define __FUNCT__ "writeResidual"
template <unsigned int DIM>
PetscErrorCode PetscAdLem3D<DIM>::writeResidual(std::string resultPath)
{
    PetscErrorCode ierr;
    Vec res;
    PetscFunctionBeginUser;
    ierr = DMCreateGlobalVector(mDa,&res);CHKERRQ(ierr);
    ierr = VecSet(res,0);CHKERRQ(ierr);
    ierr = VecAXPY(res,-1.0,mB);CHKERRQ(ierr);
    ierr = MatMultAdd(mA,mX,res,res);

    std::string fileName = resultPath + "res.vts";
    PetscViewer viewer;
    ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,fileName.c_str(),FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)res,"residual");CHKERRQ(ierr);
    ierr = VecView(res,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    ierr = VecDestroy(&res);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "getSolutionArray"
template <unsigned int DIM>
PetscErrorCode PetscAdLem3D<DIM>::getSolutionArray()
{
    PetscErrorCode ierr;
    Vec xNatural;
    ierr = DMDACreateNaturalVector(mDa,&xNatural);CHKERRQ(ierr);
    ierr = DMDAGlobalToNaturalBegin(mDa,mX,INSERT_VALUES,xNatural);CHKERRQ(ierr);
    ierr = DMDAGlobalToNaturalEnd(mDa,mX,INSERT_VALUES,xNatural);CHKERRQ(ierr);
    if(!mSolAllocated)
        ierr = VecScatterCreateToAll(xNatural,&mScatterCtx,&mXLocal);CHKERRQ(ierr);
    ierr = VecScatterBegin(mScatterCtx,xNatural,mXLocal,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(mScatterCtx,xNatural,mXLocal,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecGetArray(mXLocal,&mSolArray);CHKERRQ(ierr);
    ierr = VecDestroy(&xNatural);CHKERRQ(ierr);
    //    ierr = VecView(mXLocal,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    mSolAllocated = PETSC_TRUE;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "getRhsArray"
template <unsigned int DIM>
PetscErrorCode PetscAdLem3D<DIM>::getRhsArray()
{
    PetscErrorCode ierr;
    Vec bNatural;
    ierr = DMDACreateNaturalVector(mDa,&bNatural);CHKERRQ(ierr);
    ierr = DMDAGlobalToNaturalBegin(mDa,mB,INSERT_VALUES,bNatural);CHKERRQ(ierr);
    ierr = DMDAGlobalToNaturalEnd(mDa,mB,INSERT_VALUES,bNatural);CHKERRQ(ierr);
    if(!mRhsAllocated)
        ierr = VecScatterCreateToAll(bNatural,&mScatterRhsCtx,&mBLocal);CHKERRQ(ierr);
    ierr = VecScatterBegin(mScatterRhsCtx,bNatural,mBLocal,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(mScatterRhsCtx,bNatural,mBLocal,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecGetArray(mBLocal,&mRhsArray);CHKERRQ(ierr);
    ierr = VecDestroy(&bNatural);CHKERRQ(ierr);
    //    ierr = VecView(mBLocal,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    mRhsAllocated = PETSC_TRUE;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "getRhsAt"
template <unsigned int DIM>
//The function expects pos[] to have valid position, that is within
//the size that was originally used by the user to create the model.
double PetscAdLem3D<DIM>::getRhsAt(unsigned int pos[], unsigned int component)
{
    PetscFunctionBeginUser;
    PetscInt x,y,z,xn,yn,zn;
    x = pos[0]; y=pos[1];   z=pos[2];
    //Adapt for the staggered grid DM which has one greater dimension.
    xn = mProblemModel->getXnum()+1;
    yn = mProblemModel->getYnum()+1;
    zn = mProblemModel->getZnum()+1;

    if(x<0 || x>=xn-1 || y<0 || y>=yn-1 || z<0 || z>=zn-1) //>=xn-1 because xn has one bigger size than the input size of the model!
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"out of range position asked for velocity solution.\n");

    //Here we provide staggered grid values themselves as they were used in the momentum
    //equations for each (i,j,k) cell centre.
    if(component == 0 || component == 1 || component == 2 || component == 3) {
        PetscFunctionReturn(mRhsArray[(x+ xn*y + xn*yn*z)*4 + component]);
    } else
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_PLIB,"rhs component should be 0, 1, 2 or 3\n");
}


#undef __FUNCT__
#define __FUNCT__ "getSolVelocityAt"
//The function expects pos[] to have valid position, that is within
//the size that was originally used by the user to create the model.
template <unsigned int DIM>
double PetscAdLem3D<DIM>::getSolVelocityAt(unsigned int pos[], unsigned int component)
{
    PetscFunctionBeginUser;
    PetscInt x,y,z,xn,yn,zn;
    x = pos[0]; y=pos[1];   z=pos[2];
    //Adapt for the staggered grid DM which has one greater dimension.
    xn = mProblemModel->getXnum()+1;
    yn = mProblemModel->getYnum()+1;
    zn = mProblemModel->getZnum()+1;

    if(x<0 || x>=xn-1 || y<0 || y>=yn-1 || z<0 || z>=zn-1) //>=xn-1 because xn has one bigger size than the input size of the model!
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"out of range position asked for velocity solution.\n");

    //Velocity solution were computed at faces, so interpolate to get at cell centers.
    double solAtcenter;
    if (component == 0) { //x-component vx_c = vx(i,j,k) + vx(i+1,j,k)
        solAtcenter = mSolArray[(x+ xn*y + xn*yn*z)*4 + component]
                + mSolArray[(x+1 + xn*y + xn*yn*z)*4 + component];
    } else if (component == 1) {
        solAtcenter = mSolArray[(x + xn*y + xn*yn*z)*4 + component]
                + mSolArray[(x + xn*(y+1) + xn*yn*z)*4 + component];
    } else if (component == 2) {
        solAtcenter = mSolArray[(x + xn*y + xn*yn*z)*4 + component]
                + mSolArray[(x + xn*y + xn*yn*(z+1))*4 + component];
    } else
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_PLIB,"velocity component should be 0, 1 or 2 for 3D\n");
    PetscFunctionReturn(solAtcenter/2.);
}

#undef __FUNCT__
#define __FUNCT__ "getSolPressureAt"
template <unsigned int DIM>
double PetscAdLem3D<DIM>::getSolPressureAt(unsigned int pos[])
{
    PetscFunctionBeginUser;
    PetscInt x,y,z,xn,yn,zn;
    //Adapt for the staggered grid DM which has one greater dimension.
    x = pos[0]+1; y=pos[1]+1;   z=pos[2]+1;
    xn = mProblemModel->getXnum()+1;
    yn = mProblemModel->getYnum()+1;
    zn = mProblemModel->getZnum()+1;

    if(x<1 || x>=xn || y<1 || y>=yn || z<1 || z>=zn) //<1 because 1 already added in this function!
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"out of range position asked for pressure solution.\n");
    PetscFunctionReturn(mSolArray[(x+ xn*y + xn*yn*z)*4 + 3]);
}

#undef __FUNCT__
#define __FUNCT__ "getDivergenceAt"
template <unsigned int DIM>
//the size that was originally used by the user to create the model.
double PetscAdLem3D<DIM>::getDivergenceAt(unsigned int pos[])
{
    PetscFunctionBeginUser;
    if(mIsDiv12pointStencil)
    {
	SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Divergence computation  not supported yet for 12 point stencil option.\n");
	return 0;
    }
    PetscInt x,y,z,xn,yn,zn;
    x = pos[0]; y=pos[1];   z=pos[2];
    //Adapt for the staggered grid DM which has one greater dimension.
    xn = mProblemModel->getXnum()+1;
    yn = mProblemModel->getYnum()+1;
    zn = mProblemModel->getZnum()+1;
    //Get the spacings
    PetscReal hx, hy, hz;
    hx = mProblemModel->getXspacing();
    hy = mProblemModel->getYspacing();
    hz = mProblemModel->getZspacing();

    if(x<0 || x>=xn-1 || y<0 || y>=yn-1 || z<0 || z>=zn-1) //>=xn-1 because xn has one bigger size than the input size of the model!
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"out of range position asked for velocity solution.\n");

    //Velocity solution were computed at faces, so the divergence lies in the cell centers.
    double divergence;
    //a = vx(i+1,j,k) - vx(i,j,k) + vy(i,j+1,k) - vy(i,j,k) + vz(i,j,k+1) - vz(i,j,k)
    divergence = (mSolArray[(x+1 + xn*y + xn*yn*z)*4 + 0]
		  -mSolArray[(x + xn*y + xn*yn*z)*4 + 0])/hx;
    divergence += (mSolArray[(x + xn*(y+1) + xn*yn*z)*4 + 1]
		   -mSolArray[(x + xn*y + xn*yn*z)*4 + 1])/hy;
    divergence += (mSolArray[(x + xn*y + xn*yn*(z+1))*4 + 2]
		   -mSolArray[(x + xn*y + xn*yn*z)*4 + 2])/hz;
    PetscFunctionReturn(divergence);
}

#endif // PETSCADLEM3D_HXX
