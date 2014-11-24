//Implementation of the Solver class for Poisson Problem.
//Bishesh Khanal, Asclepios INRIA Sophia Antipolis

#include <PoissonSolver.hxx>
#include <PoissonProblem.hxx>

PoissonSolver::PoissonSolver(PoissonProblem<double> *problem):
    mProblem(problem)
{
    PetscErrorCode ierr;
    ierr = KSPCreate(PETSC_COMM_WORLD, &mKsp);CHKERRXX(ierr);
    ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                        DMDA_STENCIL_STAR,mProblem->getSize(0),mProblem->getSize(1),mProblem->getSize(2),
                        PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0,0,&mDa);CHKERRXX(ierr);

    ierr = DMDASetFieldName(mDa,0,"x");CHKERRXX(ierr);

    mOperatorComputed = PETSC_FALSE;
    mResComputed = PETSC_FALSE;
}

PoissonSolver::~PoissonSolver()
{
    PetscErrorCode ierr;
    ierr = VecRestoreArray(mXLocal,&mXArray);CHKERRXX(ierr);
    ierr = VecScatterDestroy(&mScatterCtx);CHKERRXX(ierr);
    ierr = VecDestroy(&mXLocal);CHKERRXX(ierr);
    if(mResComputed) {
        ierr = VecRestoreArray(mResLocal,&mResArray);CHKERRXX(ierr);
        ierr = VecScatterDestroy(&mResScatterCtx);CHKERRXX(ierr);
        ierr = VecDestroy(&mResLocal);CHKERRXX(ierr);
        ierr = VecDestroy(&mRes);CHKERRXX(ierr);
    }

    ierr = DMDestroy(&mDa);CHKERRXX(ierr);
    ierr = KSPDestroy(&mKsp);CHKERRXX(ierr);
}

PoissonProblem<double>* PoissonSolver::getProblem() { return mProblem; }

#undef __FUNCT__
#define __FUNCT__ "solve"
PetscErrorCode PoissonSolver::solve()
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"poisson Solver called\n");

    if(!mOperatorComputed) {
        ierr = KSPSetDM(mKsp,mDa);CHKERRQ(ierr);
        ierr = KSPSetComputeOperators(mKsp,computeMatrix,this);CHKERRQ(ierr);
        ierr = KSPSetComputeRHS(mKsp,computeRhs,this);CHKERRQ(ierr);
        mOperatorComputed = PETSC_TRUE;
        ierr = KSPSetFromOptions(mKsp);CHKERRQ(ierr);
    }

    ierr = KSPSolve(mKsp,NULL,NULL);CHKERRQ(ierr);
    ierr = KSPGetSolution(mKsp,&mX);CHKERRQ(ierr);
    ierr = KSPGetRhs(mKsp,&mB);CHKERRQ(ierr);

    Vec xNatural;
    ierr = DMDACreateNaturalVector(mDa,&xNatural);CHKERRQ(ierr);
    ierr = DMDAGlobalToNaturalBegin(mDa,mX,INSERT_VALUES,xNatural);CHKERRQ(ierr);
    ierr = DMDAGlobalToNaturalEnd(mDa,mX,INSERT_VALUES,xNatural);CHKERRQ(ierr);
    ierr = VecScatterCreateToAll(xNatural,&mScatterCtx,&mXLocal);CHKERRQ(ierr);
    ierr = VecScatterBegin(mScatterCtx,xNatural,mXLocal,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(mScatterCtx,xNatural,mXLocal,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecGetArray(mXLocal,&mXArray);CHKERRQ(ierr);
    ierr = VecDestroy(&xNatural);CHKERRQ(ierr);
//    ierr = VecView(mXLocal,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    computeResidual();      //For debugging!
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "computeResidual"
PetscErrorCode PoissonSolver::computeResidual()
{
    mResComputed = PETSC_TRUE;
    PetscErrorCode ierr;
    ierr = KSPGetOperators(mKsp,&mA,NULL);CHKERRQ(ierr);
    PetscReal norm;
    ierr = DMCreateGlobalVector(mDa,&mRes);CHKERRQ(ierr);
//    ierr = VecDuplicate(mRes,&mB);CHKERRQ(ierr);
    ierr = VecSet(mRes,0);CHKERRQ(ierr);
    ierr = VecAXPY(mRes,-1.0,mB);CHKERRQ(ierr);
    ierr = MatMultAdd(mA,mX,mRes,mRes);
    ierr = VecNorm(mRes,NORM_2,&norm);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"norm of residual is: %f\n",norm);


    Vec xNatural;
    ierr = DMDACreateNaturalVector(mDa,&xNatural);CHKERRQ(ierr);
    ierr = DMDAGlobalToNaturalBegin(mDa,mRes,INSERT_VALUES,xNatural);CHKERRQ(ierr);
    ierr = DMDAGlobalToNaturalEnd(mDa,mRes,INSERT_VALUES,xNatural);CHKERRQ(ierr);
    ierr = VecScatterCreateToAll(xNatural,&mResScatterCtx,&mResLocal);CHKERRQ(ierr);
    ierr = VecScatterBegin(mResScatterCtx,xNatural,mResLocal,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd(mResScatterCtx,xNatural,mResLocal,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecGetArray(mResLocal,&mResArray);CHKERRQ(ierr);
    ierr = VecDestroy(&xNatural);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "writeToFile"
PetscErrorCode PoissonSolver::writeToFile(std::string format, std::string fileName,
                                              bool writeResidual)
{
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    PetscViewer viewer;

    ierr = PetscObjectSetName((PetscObject)mB,"Rhs");CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)mX,"Sol");CHKERRQ(ierr);

    if (format == "matlab") {
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fileName.c_str(),FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
        ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
        ierr = VecView(mB,viewer);CHKERRQ(ierr);
        ierr = VecView(mX,viewer);CHKERRQ(ierr);

        /*{ //debugging system matrix:
            ierr = KSPGetOperators(mKsp,&mA,NULL,NULL);CHKERRQ(ierr);
            ierr = PetscObjectSetName((PetscObject)mA,"A");CHKERRQ(ierr);

            std::string sysFileName(fileName + "Sys");
            PetscViewer viewerSys;
            ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,sysFileName.c_str(),FILE_MODE_WRITE,&viewerSys);CHKERRQ(ierr);
            ierr = PetscViewerSetFormat(viewerSys,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
            ierr = MatView(mA,viewerSys);CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewerSys);CHKERRQ(ierr);
        }*/
        if(writeResidual) {
            ierr = PetscObjectSetName((PetscObject)mRes,"residual");CHKERRQ(ierr);
            ierr = VecView(mRes,viewer);CHKERRQ(ierr);
        }
    } else if (format == "vtk") {
        fileName = fileName + ".vtk";
        ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,fileName.c_str(),&viewer);
        ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
        //        ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_VTK_VTS);CHKERRQ(ierr);
//        ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,fileName.c_str(),FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
        ierr = VecView(mX,viewer);CHKERRQ(ierr);
    } else {
        SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_PLIB,"unsupported write format",0);
    }

    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "computeMatrix"
PetscErrorCode PoissonSolver::computeMatrix(KSP ksp, Mat A,
                                                Mat B, void *ctx)
{

    PoissonSolver *context = (PoissonSolver*)ctx;
    PetscErrorCode ierr;
    DM da;
    ierr = KSPGetDM(ksp,&da);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
    PetscScalar v[7];
    MatStencil  row, col[7];
    //    PetscScalar bScale = 1.0;
    PetscInt num;
    row.c = 0;
    for (PetscInt i = 0; i<7; ++i)
        col[i].c = 0;


    PetscFunctionBeginUser;
    //This purely considers dirichlet boundary conditions with zero values
    //outside the domain!
    for (PetscInt k = info.zs; k<info.zs+info.zm; ++k) {
        for (PetscInt j = info.ys; j < info.ys+info.ym; ++j) {
            for (PetscInt i = info.xs; i < info.xs+info.xm; ++i) {
                num = 0;
                row.i = i;  row.j = j;  row.k = k;
//                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"accessing(%d,%d,%d)\t",i,j,k);
                PetscScalar w = context->getProblem()->getCoeffAtPosition(i,j,k);

                //k-1:
                if (k>0) {
                    col[num].i = i;   col[num].j = j;   col[num].k = k-1;
                    v[num++] = w;
                }

                //k:
                if (i > 0) {
                    col[num].i = i-1; col[num].j = j;   col[num].k = k;
                    v[num++] = w;
                }
                if (j<info.my-1) {
                    col[num].i = i;  col[num].j = j+1;  col[num].k = k;
                    v[num++] = w;
                }
                col[num].i = i;  col[num].j = j;  col[num].k = k;
                v[num++] = -6*w;
                if (j>0) {
                    col[num].i = i;  col[num].j = j-1;  col[num].k = k;
                    v[num++] = w;
                }
                if (i<info.mx-1){
                    col[num].i = i+1;col[num].j = j;  col[num].k = k;
                    v[num++] = w;
                }

                //k+1:
                if (k<info.mz-1) {
                    col[num].i = i;   col[num].j = j;   col[num].k = k+1;
                    v[num++] = w;
                }
                ierr = MatSetValuesStencil(A,1,&row,num,col,v,INSERT_VALUES);CHKERRQ(ierr);

            }
        }
    }
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "computeRhs"
PetscErrorCode PoissonSolver::computeRhs(KSP ksp, Vec b, void *ctx)
{
    PetscFunctionBeginUser;
    PoissonSolver *context = (PoissonSolver*)ctx;

    PetscErrorCode ierr;
    DM da;
    ierr = KSPGetDM(ksp,&da);CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);

    PetscScalar ***rhs;
    ierr = DMDAVecGetArray(da,b,&rhs);CHKERRQ(ierr);

    //The RHS should be provided by the model!
    for (PetscInt k = info.zs; k<info.zs+info.zm; ++k) {
        for (PetscInt j = info.ys; j < info.ys+info.ym; ++j) {
            for (PetscInt i = info.xs; i < info.xs+info.xm; ++i) {
                rhs[k][j][i] = context->getProblem()->getRhsAtPosition(i,j,k);
            }
        }
    }

    ierr = DMDAVecRestoreArray(da,b,&rhs);CHKERRQ(ierr);
    //    ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); To debug
    //    ierr = VecAssemblyBegin(b);CHKERRQ(ierr); //No need to assemble!
    //    ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "getSolAtPosition"
double PoissonSolver::getSolAtPosition(unsigned int x, unsigned int y, unsigned int z)
{
    PetscFunctionBeginUser;
    PetscInt dim[3];
    for(int i=0;i<3;++i)
        dim[i] = this->getProblem()->getSize(i);
    return mXArray[(x+ dim[0]*y + dim[0]*dim[1]*z)];

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "getResAtPosition"
double PoissonSolver::getResAtPosition(unsigned int x, unsigned int y, unsigned int z)
{
    PetscFunctionBeginUser;
    PetscInt dim[3];
    for(int i=0;i<3;++i)
        dim[i] = this->getProblem()->getSize(i);
    return mResArray[(x+ dim[0]*y + dim[0]*dim[1]*z)];

    PetscFunctionReturn(0);
}

