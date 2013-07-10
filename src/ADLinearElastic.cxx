//Implementaion for ADLinearElastic class

#include"ADLinearElastic.hxx"

#undef __FUNCT__
#define __FUNCT__ "ADLinearElastic"
//initialize with mXn (MY X MX )grid size, and with the material parameters
//currently only for the dirichlet condition.
ADLinearElastic::ADLinearElastic(const PetscInt MY, const PetscInt MX, const PetscReal rho){
    m = MY; n = MX;
    hx = 1.0 / (PetscReal)(n-1);
    hy = 1.0 / (PetscReal)(m-1);
    centerRho = rho; nu = 0.1;
    bc = ADLinearElastic::DIRICHLET;
}


#undef __FUNCT__
#define __FUNCT__ "ComputeRho"
PetscErrorCode ADLinearElastic::ComputeRho(PetscInt i, PetscInt j, PetscInt mx, PetscInt my, PetscReal *rho)
{
  PetscFunctionBeginUser;
  if ((i > mx/3.0) && (i < 2.0*mx/3.0) && (j > my/3.0) && (j < 2.0*my/3.0)) {
    *rho = centerRho;
  } else {
    *rho = 1.0;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeMatrix"
PetscErrorCode ComputeMatrix(KSP ksp, Mat J, Mat jac, MatStructure* str, void* ctx){
     ADLinearElastic    *user = (ADLinearElastic*)ctx;

     PetscErrorCode ierr;
     PetscInt       i,j,mx,my,xm,ym,xs,ys;
     PetscScalar    v[5];
     PetscReal      Hx,Hy,HydHx,HxdHy,rho;
     MatStencil     row, col[5];
     DM             da;

     PetscFunctionBeginUser;
     ierr      = KSPGetDM(ksp,&da);CHKERRQ(ierr);
     ierr      = DMDAGetInfo(da,0,&mx,&my,0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
     Hx = user->getHx(); Hy = user->getHy();
     HxdHy     = Hx/Hy;
     HydHx     = Hy/Hx;
     ierr      = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);
     for (j=ys; j<ys+ym; j++) {
       for (i=xs; i<xs+xm; i++) {
         row.i = i; row.j = j;
         ierr  = user->ComputeRho(i, j, mx, my, &rho);CHKERRQ(ierr);
         if (i==0 || j==0 || i==mx-1 || j==my-1) {
           if (user->getBcType() == user->DIRICHLET) {
             v[0] = 2.0*rho*(HxdHy + HydHx);
             ierr = MatSetValuesStencil(jac,1,&row,1,&row,v,INSERT_VALUES);CHKERRQ(ierr);
           } else if (user->getBcType() == user->NEUMANN) {
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
             ierr = MatSetValuesStencil(jac,1,&row,num,col,v,INSERT_VALUES);CHKERRQ(ierr);
           }
         } else {
           v[0] = -rho*HxdHy;              col[0].i = i;   col[0].j = j-1;
           v[1] = -rho*HydHx;              col[1].i = i-1; col[1].j = j;
           v[2] = 2.0*rho*(HxdHy + HydHx); col[2].i = i;   col[2].j = j;
           v[3] = -rho*HydHx;              col[3].i = i+1; col[3].j = j;
           v[4] = -rho*HxdHy;              col[4].i = i;   col[4].j = j+1;
           ierr = MatSetValuesStencil(jac,1,&row,5,col,v,INSERT_VALUES);CHKERRQ(ierr);
         }
       }
     }
     ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
     ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
     if (user->getBcType() == user->NEUMANN) {
       MatNullSpace nullspace;

       ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);CHKERRQ(ierr);
       ierr = MatSetNullSpace(jac,nullspace);CHKERRQ(ierr);
       ierr = MatNullSpaceDestroy(&nullspace);CHKERRQ(ierr);
     }
     PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
PetscErrorCode ComputeRHS(KSP ksp,Vec b,void *ctx)
{
  ADLinearElastic    *user = (ADLinearElastic*)ctx;
  PetscErrorCode ierr;
  PetscInt       i,j,mx,my,xm,ym,xs,ys;
  PetscScalar    Hx,Hy;
  PetscScalar    **array;
  DM             da;
  PetscReal      nu;

  PetscFunctionBeginUser;
  ierr = KSPGetDM(ksp,&da);CHKERRQ(ierr);
  ierr = DMDAGetInfo(da, 0, &mx, &my, 0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
//  Hx   = 1.0 / (PetscReal)(mx-1);
//  Hy   = 1.0 / (PetscReal)(my-1);
  Hx = user->getHx(); Hy = user->getHy();
  nu = user->getNu();
  ierr = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, b, &array);CHKERRQ(ierr);
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      array[j][i] = PetscExpScalar(-((PetscReal)i*Hx)*((PetscReal)i*Hx)/nu)*PetscExpScalar(-((PetscReal)j*Hy)*((PetscReal)j*Hy)/nu)*Hx*Hy;
    }
  }
  ierr = DMDAVecRestoreArray(da, b, &array);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

  /* force right hand side to be consistent for singular matrix */
  /* note this is really a hack, normally the model would provide you with a consistent right handside */
  if (user->getBcType() == user->NEUMANN) {
    MatNullSpace nullspace;

    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);CHKERRQ(ierr);
    ierr = MatNullSpaceRemove(nullspace,b,NULL);CHKERRQ(ierr);
    ierr = MatNullSpaceDestroy(&nullspace);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

