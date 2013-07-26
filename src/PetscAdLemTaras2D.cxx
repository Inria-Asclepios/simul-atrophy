#include"PetscAdLemTaras2D.hxx"


PetscAdLemTaras2D::PetscAdLemTaras2D(AdLem2D *model):
    PetscAdLem2D(model,std::string("Taras Method"))
{
    PetscErrorCode ierr;

    //Linear Solver context:
    ierr = KSPCreate(PETSC_COMM_WORLD,&mKsp);CHKERRXX(ierr);
    ierr = DMDACreate2d(PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE,
                        DMDA_STENCIL_BOX,model->getXnum()+1,model->getYnum()+1,PETSC_DECIDE,PETSC_DECIDE,
                        3,1,0,0,&mDa);CHKERRXX(ierr); //3 dof., node_grid_num = cells_num + 1
    ierr = DMDASetUniformCoordinates(mDa,0,1,0,1,0,0);CHKERRXX(ierr);
    ierr = DMDASetFieldName(mDa,0,"vx");CHKERRXX(ierr);
    ierr = DMDASetFieldName(mDa,1,"vy");CHKERRXX(ierr);
    ierr = DMDASetFieldName(mDa,2,"p");CHKERRXX(ierr);
}

PetscAdLemTaras2D::~PetscAdLemTaras2D()
{
    PetscErrorCode ierr;
    ierr = DMDestroy(&mDa);CHKERRXX(ierr);
    ierr = KSPDestroy(&mKsp);CHKERRXX(ierr);
}

PetscReal PetscAdLemTaras2D::getP0Cell() { return 4.0; }

#undef __FUNCT__
#define __FUNCT__ "dataCenterAt"

/*AdLem2D class has data only at the cell centers.
The dimension of ghosted cell-centered data is greater by one in each
direction in Taras method. Hence, change the co-ordinate to get proper value
of the data at the cell center.*/

PetscReal PetscAdLemTaras2D::dataCenterAt(std::string dType, PetscInt i, PetscInt j)
{
    if(i != 0)
        --i;
    if(j != 0)
        --j;
    return this->getProblemModel()->dataAt(dType,i,j);
}

#undef __FUNCT__
#define __FUNCT__ "dataNodeAt"
PetscReal PetscAdLemTaras2D::dataNodeAt(std::string dType, PetscInt i, PetscInt j)
{
    //last corner point, simply return the corner center value.
    if ((i == this->getProblemModel()->getXnum())
            && (j==this->getProblemModel()->getYnum())) {
        return dataCenterAt(dType,i,j);
    }

    //grid-points on right edge
    if (i == this->getProblemModel()->getXnum())
        return 0.5*( dataCenterAt(dType,i,j) + dataCenterAt(dType,i,j+1));

    //grid-points on top edge
    if (j == this->getProblemModel()->getYnum())
        return 0.5*( dataCenterAt(dType,i,j) + dataCenterAt(dType,i+1,j));

    //grid-points on bottom and left edge do not need special condition
    //because, there are ghost cells present outside the grid nodes here.
    //so, it will be handeled by dataCenterAt call.
    return 0.25 * ( dataCenterAt(dType,i,j) + dataCenterAt(dType,i,j+1)
                    +dataCenterAt(dType,i+1,j) + dataCenterAt(dType,i+1,j+1));

}

#undef __FUNCT__
#define __FUNCT__ "muCenter"

PetscReal PetscAdLemTaras2D::muCenter(PetscInt i, PetscInt j)
{
    return dataCenterAt("mu",i,j);
}

#undef __FUNCT__
#define __FUNCT__ "muNode"
PetscReal PetscAdLemTaras2D::muNode(PetscInt i, PetscInt j)
{
    return dataNodeAt("mu",i,j);
}

#undef __FUNCT__
#define __FUNCT__ "lamdaCenter"

PetscReal PetscAdLemTaras2D::lambdaCenter(PetscInt i, PetscInt j)
{
    return dataCenterAt("lambda",i,j);
}

#undef __FUNCT__
#define __FUNCT__ "lambdaNode"
PetscReal PetscAdLemTaras2D::lambdaNode(PetscInt i, PetscInt j)
{
    return dataNodeAt("lambda",i,j);
}


#undef __FUNCT__
#define __FUNCT__ "aCenter"

PetscReal PetscAdLemTaras2D::aCenter(PetscInt i, PetscInt j)
{
    return dataCenterAt("atrophy",i,j);
}

#undef __FUNCT__
#define __FUNCT__ "aNode"
PetscReal PetscAdLemTaras2D::aNode(PetscInt i, PetscInt j)
{
    return dataNodeAt("atrophy",i,j);
}

#undef __FUNCT__
#define __FUNCT__ "solveModel"
PetscErrorCode PetscAdLemTaras2D::solveModel(bool writeToMatlab, const std::string& filename){
    PetscErrorCode ierr;
    Vec b,x;
    PetscFunctionBeginUser;
    ierr = DMKSPSetComputeRHS(mDa,computeRHSTaras2D,this);CHKERRQ(ierr);
    ierr = DMKSPSetComputeOperators(mDa,computeMatrixTaras2D,this);CHKERRQ(ierr);
    ierr = KSPSetDM(mKsp,mDa);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(mKsp);CHKERRQ(ierr);
    ierr = KSPSetUp(mKsp);CHKERRQ(ierr);
    ierr = KSPSolve(mKsp,NULL,NULL);CHKERRQ(ierr);
    ierr = KSPGetSolution(mKsp,&x);CHKERRQ(ierr);
    ierr = KSPGetRhs(mKsp,&b);CHKERRQ(ierr);
    if (writeToMatlab) {
//        Mat mat1, mat2;
//        ierr = KSPGetOperators(mKsp,&mat1,&mat2,0);
        PetscViewer viewer;
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename.c_str(),FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
        ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);

//        ierr = PetscObjectSetName((PetscObject)mat1,"A");
        ierr = PetscObjectSetName((PetscObject)x,"x");CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)b,"b");CHKERRQ(ierr);
        ierr = VecView(x,viewer);CHKERRQ(ierr);
        ierr = VecView(b,viewer);CHKERRQ(ierr);
//        ierr = MatView(mat1,viewer);CHKERRQ(ierr); //doesn't work, probably not supported?
        ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);

}


#undef __FUNCT__
#define __FUNCT__ "computeMatrixTaras2D"
PetscErrorCode PetscAdLemTaras2D::computeMatrixTaras2D(KSP ksp, Mat J, Mat jac, MatStructure* str, void* ctx)
{
    PetscAdLemTaras2D    *user = (PetscAdLemTaras2D*)ctx;

    PetscErrorCode ierr;
    PetscInt       i,j,mx,my,xm,ym,xs,ys;
    PetscReal      Hx,Hy,HydHx,HxdHy;
    PetscScalar    v[11];
    MatStencil     row, col[11];
    DM             da;

    PetscFunctionBeginUser;
    ierr      = KSPGetDM(ksp,&da);CHKERRQ(ierr);
    ierr      = DMDAGetInfo(da,0,&mx,&my,0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
    Hx = 1./(mx-1);
    Hy = 1./(my-1);
    HxdHy     = Hx/Hy;
    HydHx     = Hy/Hx;
    ierr      = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRQ(ierr);

    if (user->getProblemModel()->getBcType() != user->getProblemModel()->DIRICHLET) {
        SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"only Dirichlet boundary condition implemented",0);
    }
    for (j=ys; j<ys+ym; j++) {
        for (i=xs; i<xs+xm; i++) {
            row.i = i; row.j = j;
            // ********************* x-momentum equation ************************
            row.c = 0;
            if (i==0 || i==mx-1 || j==0 || j==my-2 || j==my-1) {
                //left and right edges
                if (i==0 || i==mx-1) {
                    //vx-coefficients: vx(i,j) = 0
                    col[0].c = 0;
                    v[0] = 1.0;             col[0].i = i;       col[0].j = j;
                    ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                } else {    //bottom and top edges
                    // ghost values for j = my-1: vx(i,j) = 0
                    if (j == my-1) {
                        col[0].c = 0;
                        v[0] = 1.0;             col[0].i = i;       col[0].j = j;
                        ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                    } else { //boundary values
                        //vx-coefficients: 3*vx(i,0) - vx(i,1) = 0; 3*vx(i,my-2) - vx(i,my-3) = 0
                        col[0].c = 0;           col[1].c = 0;
                        v[0] = 3.0;             col[0].i = i;       col[0].j = j;
                        if (j==0) {
                            v[1] = -1.0;        col[1].i = i;       col[1].j = j+1;
                        } else if (j==my-2){
                            v[1] = -1.0;        col[1].i = i;       col[1].j = j-1;
                        }
                        ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                    }
                }
            } else {
                //vx-coefficients, five terms.
                for(int ii=0; ii<=4; ++ii)
                    col[ii].c = 0;

                v[0] = 2.0 * user->muCenter(i+1,j+1) * HydHx;    col[0].i = i+1;     col[0].j = j;
                v[1] = 2.0 * user->muCenter(i,j+1) * HydHx;      col[1].i = i-1;     col[1].j = j;
                v[2] = user->muNode(i,j+1) * HxdHy;              col[2].i = i;       col[2].j = j+1;
                v[3] = user->muNode(i,j) * HxdHy;                col[3].i = i;       col[3].j = j-1;
                v[4] = -v[0] - v[1] - v[2] - v[3];               col[4].i = i;       col[4].j = j;

                //vy-coefficients, four terms.
                for(int ii=5; ii<= 8; ++ii)
                    col[ii].c = 1;

                v[5] = user->muNode(i,j+1);                      col[5].i = i;       col[5].j = j+1;
                v[6] = -v[5];                                    col[6].i = i-1;     col[6].j = j+1;
                v[7] = user->muNode(i,j);                        col[7].i = i-1;     col[7].j = j;
                v[8] = -v[7];                                    col[8].i = i;       col[8].j = j;

                //p-coefficients, two terms.
                col[9].c = 2;    col[10].c = 2;
                v[9] = Hy;  //HxHy/Hx = Hy                                 col[9].i = i;       col[9].j = j+1;
                v[10] = -v[9];                                   col[10].i = i+1;    col[10].j = j+1;

                ierr = MatSetValuesStencil(jac,1,&row,11,col,v,INSERT_VALUES);CHKERRQ(ierr);
            }

            // ********************* y-momentum equation ************************
            row.c = 1;
            if (i==0 || i==mx-2 || i==mx-1 || j==0 || j==my-1) {
                //Bottom and top edges: vy(i,j) = 0
                if (j==0 || j==my-1) {
                    col[0].c = 1;
                    v[0] = 1.0;             col[0].i = i;       col[0].j = j;
                    ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                } else if (i == mx-1) {                 //Right edge ghost values:
                // ghost values for i = mx-1: vy(i,j) = 0
                    col[0].c = 1;
                    v[0] = 1.0;             col[0].i = i;       col[0].j = j;
                    ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                } else { //left and right boundary values.
                    //vy-coefficients: 3*vy(0,j) - vy(1,j) = 0; 3*vy(mx-2,j) - vy(mx-3,j) = 0
                    col[0].c = 1;           col[1].c = 1;
                    v[0] = 3.0;             col[0].i = i;       col[0].j = j;
                    if (i==0) {
                        v[1] = -1.0;        col[1].i = i+1;     col[1].j = j;
                    } else if (i==mx-2){
                        v[1] = -1.0;        col[1].i = i-1;     col[1].j = j;
                    }
                    ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                }
            } else {
                //vy-coefficients, five terms.
                for(int ii=0; ii<=4; ++ii)
                    col[ii].c = 1;

                v[0] = 2.0 * user->muCenter(i+1,j+1) * HxdHy;    col[0].i = i;       col[0].j = j+1;
                v[1] = 2.0 * user->muCenter(i+1,j) * HxdHy;      col[1].i = i;       col[1].j = j-1;
                v[2] = user->muNode(i+1,j) * HydHx;              col[2].i = i+1;     col[2].j = j;
                v[3] = user->muNode(i,j) * HydHx;                col[3].i = i-1;     col[3].j = j;
                v[4] = -v[0] - v[1] - v[2] - v[3];               col[4].i = i;       col[4].j = j;

                //vx-coefficients, four terms.
                for(int ii=5; ii<= 8; ++ii)
                    col[ii].c = 0;

                v[5] = user->muNode(i+1,j);                      col[5].i = i+1;     col[5].j = j;
                v[6] = -v[5];                                    col[6].i = i+1;     col[6].j = j-1;
                v[7] = user->muNode(i,j);                        col[7].i = i;       col[7].j = j-1;
                v[8] = -v[7];                                    col[8].i = i;       col[8].j = j;

                //p-coefficients, two terms.
                col[9].c = 2;    col[10].c = 2;
                v[9] = Hx;   //HxHy/Hy = Hx                                col[9].i = i+1;     col[9].j = j;
                v[10] = -v[9];                                   col[10].i = i+1;    col[10].j = j+1;

                ierr = MatSetValuesStencil(jac,1,&row,11,col,v,INSERT_VALUES);CHKERRQ(ierr);

            }

            // ********************** conservation equation *****************
            row.c = 2;
            if (i==0 || j==0 || (i==1 && j==1) || (i==1 && j==my-1)
                    || (i==mx-1 && j==1) || (i==mx-1 && j==my-1)
                    || (i==1 && j==2)) {
                if (i==0 || j==0) { //p(i,j) = 0;ghost pressure:
                    col[0].c = 2;
                    v[0] = 1.0;     col[0].i = i;           col[0].j = j;
                    ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                } else if ((i==1 && j==1) || (i==1 && j==my-1)) { //left upper-lower corners
                    //set dp/dx = 0 i.e. p(i+1,j) - p(i,j) = 0
                    col[0].c = 2;           col[1].c = 2;
                    v[0] = 1.0;             col[0].i = i+1;     col[0].j = j;
                    v[1] = -1.0;            col[1].i = i;       col[1].j = j;
                    ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                } else if ((i==mx-1 && j==1) || (i==mx-1 && j==my-1)) { //right upper-lower corners
                    //set dp/dx = 1 i.e. p(i,j) - p(i-1,j) = 0
                    col[0].c = 2;           col[1].c = 2;
                    v[0] = 1.0;             col[0].i = i;     col[0].j = j;
                    v[1] = -1.0;            col[1].i = i-1;       col[1].j = j;
                    ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                } else { //one cell
                    //set p = pCell NOTE: RHS needs to be set to pCell!
                    col[0].c = 2;
                    v[0] = 1.0;             col[0].i = i;       col[0].j = j;
                    ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                }
            } else {

                //vx-coefficients, two terms
                col[0].c = 0;        col[1].c = 0;
                v[0] = 1.0/Hx;                                   col[0].i = i;       col[0].j = j-1;
                v[1] = -v[0];                                    col[1].i = i-1;     col[1].j = j-1;

                //vy-coefficients, two terms
                col[2].c = 1;        col[3].c = 1;
                v[2] = 1.0/Hy;                                   col[2].i = i-1;     col[2].j = j;
                v[3] = -v[2];                                    col[3].i = i-1;     col[3].j = j-1;

                ierr = MatSetValuesStencil(jac,1,&row,4,col,v,INSERT_VALUES);CHKERRQ(ierr);
            }
        }
    }
    ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    /*
    if (user->getProblemModel()->getBcType() == user->getProblemModel()->NEUMANN) {
        MatNullSpace nullspace;

        ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);CHKERRXX(ierr);
        ierr = MatSetNullSpace(jac,nullspace);CHKERRXX(ierr);
        ierr = MatNullSpaceDestroy(&nullspace);CHKERRXX(ierr);
    }*/
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "computeRHSTaras2D"
PetscErrorCode PetscAdLemTaras2D::computeRHSTaras2D(KSP ksp, Vec b, void *ctx)
{
    PetscAdLemTaras2D    *user = (PetscAdLemTaras2D*)ctx;
    PetscErrorCode ierr;
    PetscInt       i,j,mx,my,xm,ym,xs,ys;
    PetscScalar    Hx,Hy;
    PetscAdLemTaras2D::Field    **rhs;
    DM             da;

    PetscFunctionBeginUser;
    ierr = KSPGetDM(ksp,&da);CHKERRXX(ierr);
    ierr = DMDAGetInfo(da, 0, &mx, &my, 0,0,0,0,0,0,0,0,0,0);CHKERRXX(ierr);
    Hx   = 1.0 / (PetscReal)(mx-1);
    Hy   = 1.0 / (PetscReal)(my-1);
    ierr = DMDAGetCorners(da,&xs,&ys,0,&xm,&ym,0);CHKERRXX(ierr);
    ierr = DMDAVecGetArray(da, b, &rhs);CHKERRXX(ierr);
/*    for (j=ys; j<ys+ym; j++) {
        for (i=xs; i<xs+xm; i++) {
            array[j][i] = PetscExpScalar(-((PetscReal)i*Hx)*((PetscReal)i*Hx)/nu)*PetscExpScalar(-((PetscReal)j*Hy)*((PetscReal)j*Hy)/nu)*Hx*Hy;
        }
    }*/

    if (user->getProblemModel()->getBcType() != user->getProblemModel()->DIRICHLET) {
        SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"only Dirichlet boundary condition implemented",0);
    }
    for (j=ys; j<ys+ym; j++) {
        for (i=xs; i<xs+xm; i++) {
            // ********************* x-momentum equation ************************
            if (i==0 || i==mx-1 || j==0 || j==my-2 || j==my-1) {
                rhs[j][i].vx = 0;
            } else {
                rhs[j][i].vx = Hy*(user->muCenter(i+1,j+1) + user->muCenter(i,j+1)
                               + user->lambdaCenter(i+1,j+1) + user->lambdaCenter(i,j+1)
                                )*(user->aCenter(i+1,j+1) - user->aCenter(i,j+1))
                        / 2.0;
            }

            // ********************* y-momentum equation ************************
            if (i==0 || i==mx-2 || i==mx-1 || j==0 || j==my-1) {
                rhs[j][i].vy = 0;
            } else {
                rhs[j][i].vy = Hx*(user->muCenter(i+1,j+1) + user->muCenter(i+1,j)
                               + user->lambdaCenter(i+1,j+1) + user->lambdaCenter(i+1,j)
                                )*(user->aCenter(i+1,j+1) - user->aCenter(i+1,j))
                        / 2.0;
            }

            // ********************** conservation equation *****************
            if (i==0 || j==0 || (i==1 && j==1) || (i==1 && j==my-1)
                    || (i==mx-1 && j==1) || (i==mx-1 && j==my-1)
                    ) {
                rhs[j][i].p = 0;
            } else if (i==1 && j==2) {
                rhs[j][i].p = user->getP0Cell();
            } else {
                rhs[j][i].p = user->aCenter(i,j);
            }
        }
    }

    ierr = DMDAVecRestoreArray(da, b, &rhs);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

    /* force right hand side to be consistent for singular matrix */
    /* note this is really a hack, normally the model would provide you with a consistent right handside */
   /* if (user->getProblemModel()->getBcType() == user->getProblemModel()->NEUMANN) {
        MatNullSpace nullspace;

        ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,&nullspace);CHKERRXX(ierr);
        ierr = MatNullSpaceRemove(nullspace,b,NULL);CHKERRXX(ierr);
        ierr = MatNullSpaceDestroy(&nullspace);CHKERRXX(ierr);
    }*/
    PetscFunctionReturn(0);
}

