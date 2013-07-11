#include"PetscAdLemTaras2D.hxx"


PetscAdLemTaras2D::PetscAdLemTaras2D(AdLem2D *model):
    PetscAdLem2D(model,std::string("Taras Method")){
    PetscErrorCode ierr;

    //Linear Solver context:
    ierr = KSPCreate(PETSC_COMM_WORLD,&mKsp);CHKERRXX(ierr);
    ierr = DMDACreate2d(PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE,
                        DMDA_STENCIL_STAR,model->getXnum(),model->getYnum(),PETSC_DECIDE,PETSC_DECIDE,
                        1,1,0,0,&mDaPara);CHKERRXX(ierr);
    ierr = DMDASetUniformCoordinates(mDaPara,0,1,0,1,0,0);CHKERRXX(ierr);
    ierr = DMDASetFieldName(mDaPara,0,"Pressure");CHKERRXX(ierr);

}

PetscAdLemTaras2D::~PetscAdLemTaras2D(){
    PetscErrorCode ierr;
    ierr = DMDestroy(&mDaPara);CHKERRXX(ierr);
    ierr = KSPDestroy(&mKsp);CHKERRXX(ierr);
}



#undef __FUNCT__
#define __FUNCT__ "solveModel"
PetscErrorCode PetscAdLemTaras2D::solveModel(){
    PetscErrorCode ierr;
    Vec b,x;
    PetscFunctionBeginUser;
    ierr = KSPSetComputeRHS(mKsp,computeRHSTaras2D,this);CHKERRXX(ierr);
    ierr = KSPSetComputeOperators(mKsp,computeMatrixTaras2D,this);CHKERRXX(ierr);
    ierr = KSPSetDM(mKsp,mDaPara);CHKERRXX(ierr);
    ierr = KSPSetFromOptions(mKsp);CHKERRXX(ierr);
    ierr = KSPSetUp(mKsp);CHKERRXX(ierr);
    ierr = KSPSolve(mKsp,NULL,NULL);CHKERRXX(ierr);
    ierr = KSPGetSolution(mKsp,&x);CHKERRXX(ierr);
    ierr = KSPGetRhs(mKsp,&b);CHKERRXX(ierr);

    PetscFunctionReturn(0);

}


PetscReal PetscAdLemTaras2D::getNu() const{
    return 2.;
}

#undef __FUNCT__
#define __FUNCT__ "computeMatrixTaras2D"
PetscErrorCode computeMatrixTaras2D(KSP ksp, Mat J, Mat jac, MatStructure* str, void* ctx){
    PetscAdLemTaras2D    *user = (PetscAdLemTaras2D*)ctx;

    PetscErrorCode ierr;
    PetscInt       i,j,mx,my,xm,ym,xs,ys;
    PetscScalar    v[5];
    PetscReal      Hx,Hy,HydHx,HxdHy,rho;
    MatStencil     row, col[5];
    DM             da;

    PetscFunctionBeginUser;
    ierr      = KSPGetDM(ksp,&da);CHKERRXX(ierr);
    ierr      = DMDAGetInfo(da,0,&mx,&my,0,0,0,0,0,0,0,0,0,0);CHKERRXX(ierr);
    Hx = 1./(user->getProblemModel()->getXnum()-1);
    Hy = 1./(user->getProblemModel()->getYnum()-1);
    HxdHy     = Hx/Hy;
    HydHx     = Hy/Hx;
    ierr      = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRXX(ierr);
    for (j=ys; j<ys+ym; j++) {
      for (i=xs; i<xs+xm; i++) {
        row.i = i; row.j = j;
        ierr  = user->ComputeRho(i, j, mx, my, &rho);CHKERRXX(ierr);
        if (i==0 || j==0 || i==mx-1 || j==my-1) {
          if (user->getProblemModel()->getBcType() == user->getProblemModel()->DIRICHLET) {
            v[0] = 2.0*rho*(HxdHy + HydHx);
            ierr = MatSetValuesStencil(jac,1,&row,1,&row,v,INSERT_VALUES);CHKERRXX(ierr);
          } else if (user->getProblemModel()->getBcType() == user->getProblemModel()->NEUMANN) {
            PetscInt numx = 0, numy = 0, num = 0;
            if (j!=0) {
              v[num] = -rho*HxdHy;              col[num].i = i;   col[num].j = j-1;
              numy++; num++;
            }
            if (i!=0) {
              v[num] = -rho*HydHx;              col[num].i = i-1; col[num].j = j;
              numx++; num++;
            }
            if (i!=mx-1) {
              v[num] = -rho*HydHx;              col[num].i = i+1; col[num].j = j;
              numx++; num++;
            }
            if (j!=my-1) {
              v[num] = -rho*HxdHy;              col[num].i = i;   col[num].j = j+1;
              numy++; num++;
            }
            v[num] = numx*rho*HydHx + numy*rho*HxdHy; col[num].i = i;   col[num].j = j;
            num++;
            ierr = MatSetValuesStencil(jac,1,&row,num,col,v,INSERT_VALUES);CHKERRXX(ierr);
          }
        } else {
          v[0] = -rho*HxdHy;              col[0].i = i;   col[0].j = j-1;
          v[1] = -rho*HydHx;              col[1].i = i-1; col[1].j = j;
          v[2] = 2.0*rho*(HxdHy + HydHx); col[2].i = i;   col[2].j = j;
          v[3] = -rho*HydHx;              col[3].i = i+1; col[3].j = j;
          v[4] = -rho*HxdHy;              col[4].i = i;   col[4].j = j+1;
          ierr = MatSetValuesStencil(jac,1,&row,5,col,v,INSERT_VALUES);CHKERRXX(ierr);
        }
      }
    }
    ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRXX(ierr);
    ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRXX(ierr);
    if (user->getProblemModel()->getBcType() == user->getProblemModel()->NEUMANN) {
      MatNullSpace nullspace;

      ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);CHKERRXX(ierr);
      ierr = MatSetNullSpace(jac,nullspace);CHKERRXX(ierr);
      ierr = MatNullSpaceDestroy(&nullspace);CHKERRXX(ierr);
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "computeRHSTaras2D"
PetscErrorCode computeRHSTaras2D(KSP ksp, Vec b, void *ctx){
    PetscAdLemTaras2D    *user = (PetscAdLemTaras2D*)ctx;
    PetscErrorCode ierr;
    PetscInt       i,j,mx,my,xm,ym,xs,ys;
    PetscScalar    Hx,Hy;
    PetscScalar    **array;
    DM             da;
    PetscReal      nu;

    PetscFunctionBeginUser;
    ierr = KSPGetDM(ksp,&da);CHKERRXX(ierr);
    ierr = DMDAGetInfo(da, 0, &mx, &my, 0,0,0,0,0,0,0,0,0,0);CHKERRXX(ierr);
  //  Hx   = 1.0 / (PetscReal)(mx-1);
  //  Hy   = 1.0 / (PetscReal)(my-1);
    Hx = 1.0/ (user->getProblemModel()->getXnum()-1);
    Hy = 1.0/ (user->getProblemModel()->getYnum()-1);
    nu = user->getNu();
    ierr = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRXX(ierr);
    ierr = DMDAVecGetArray(da, b, &array);CHKERRXX(ierr);
    for (j=ys; j<ys+ym; j++) {
      for (i=xs; i<xs+xm; i++) {
        array[j][i] = PetscExpScalar(-((PetscReal)i*Hx)*((PetscReal)i*Hx)/nu)*PetscExpScalar(-((PetscReal)j*Hy)*((PetscReal)j*Hy)/nu)*Hx*Hy;
      }
    }
    ierr = DMDAVecRestoreArray(da, b, &array);CHKERRXX(ierr);
    ierr = VecAssemblyBegin(b);CHKERRXX(ierr);
    ierr = VecAssemblyEnd(b);CHKERRXX(ierr);

    /* force right hand side to be consistent for singular matrix */
    /* note this is really a hack, normally the model would provide you with a consistent right handside */
    if (user->getProblemModel()->getBcType() == user->getProblemModel()->NEUMANN) {
      MatNullSpace nullspace;

      ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);CHKERRXX(ierr);
      ierr = MatNullSpaceRemove(nullspace,b,NULL);CHKERRXX(ierr);
      ierr = MatNullSpaceDestroy(&nullspace);CHKERRXX(ierr);
    }
    PetscFunctionReturn(0);
}

