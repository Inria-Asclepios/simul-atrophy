#include "PetscAdLemTaras3D.hxx"

PetscAdLemTaras3D::PetscAdLemTaras3D(AdLem3D *model):
    PetscAdLem3D(model, std::string("Taras Method"))
{
    PetscErrorCode ierr;

    //Linear Solver context:
    ierr = KSPCreate(PETSC_COMM_WORLD,&mKsp);CHKERRXX(ierr);

    //DMDA with 4 dof; node_grid_num = cells_num + 1
    ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
                        DMDA_STENCIL_BOX,model->getXnum()+1,model->getYnum()+1,model->getZnum()+1,
                        PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,4,1,0,0,0,&mDa);CHKERRXX(ierr);

    ierr = DMDASetUniformCoordinates(mDa,0,1,0,1,0,1);CHKERRXX(ierr);
//    ierr = DMDASetUniformCoordinates(mDa,0,model->getXnum(),0,model->getYnum(),0,model->getZnum());CHKERRXX(ierr);
    ierr = DMDASetFieldName(mDa,0,"vx");CHKERRXX(ierr);
    ierr = DMDASetFieldName(mDa,1,"vy");CHKERRXX(ierr);
    ierr = DMDASetFieldName(mDa,2,"vz");CHKERRXX(ierr);
    ierr = DMDASetFieldName(mDa,3,"p");CHKERRXX(ierr);
}

PetscAdLemTaras3D::~PetscAdLemTaras3D()
{
    PetscErrorCode ierr;
    ierr = DMDestroy(&mDa);CHKERRXX(ierr);
    ierr = KSPDestroy(&mKsp);CHKERRXX(ierr);
}

PetscReal PetscAdLemTaras3D::getP0Cell() { return 4.0; }

#undef __FUNCT__
#define __FUNCT__ "dataCenterAt"

/*AdLem3D class has data only at the cell centers.
The dimension of ghosted cell-centered data is greater by one in each
direction in Taras method. Hence, change the co-ordinate to get proper value
of the data at the cell center.*/
/*This is equivalent to increase the dimension of the image by one by copying the
  value of the faces sharing the origin to it's new neigbhouring face*/
PetscReal PetscAdLemTaras3D::dataCenterAt(std::string dType, PetscInt x, PetscInt y, PetscInt z)
{
    if(x != 0)
        --x;
    if(y != 0)
        --y;
    if(z != 0)
        --z;
    return this->getProblemModel()->dataAt(dType,x,y,z);
}

#undef __FUNCT__
#define __FUNCT__ "dataXyAt"
PetscReal PetscAdLemTaras3D::dataXyAt(std::string dType, PetscInt x, PetscInt y, PetscInt z)
{
    //For referencing, easier to say x-start, x-end face, edge etc!
    //on z-end faces: same as the values at z-1!
    if (z == this->getProblemModel()->getZnum())
        --z;

    //x-end and y-end corner point, simply return the corner value.
    if ((x == this->getProblemModel()->getXnum())
            && (y == this->getProblemModel()->getYnum())) {
        return dataCenterAt(dType,x,y,z+1);
    }

    //on x-end faces:
    if (x == this->getProblemModel()->getXnum())
        return 0.5*(dataCenterAt(dType,x,y,z+1) + dataCenterAt(dType,x,y+1,z+1));

    //on y-end faces:
    if (y == this->getProblemModel()->getYnum())
        return 0.5*(dataCenterAt(dType,x,y,z+1) + dataCenterAt(dType,x+1,y,z+1));

    return 0.25 * (dataCenterAt(dType,x,y,z+1) + dataCenterAt(dType,x,y+1,z+1)
                   + dataCenterAt(dType,x+1,y,z+1) + dataCenterAt(dType,x+1,y+1,z+1));
}

#undef __FUNCT__
#define __FUNCT__ "dataXz"
PetscReal PetscAdLemTaras3D::dataXzAt(std::string dType, PetscInt x, PetscInt y, PetscInt z)
{
    //on y-end faces: same as the values at y-1!
    if (y == this->getProblemModel()->getYnum())
        --y;

    //x-end and z-end corner point, simply return the corner value.
    if ((x == this->getProblemModel()->getXnum())
            && (z == this->getProblemModel()->getZnum())) {
        return dataCenterAt(dType,x,y+1,z);
    }

    //on x-end faces:
    if (x == this->getProblemModel()->getXnum())
        return 0.5*(dataCenterAt(dType,x,y+1,z) + dataCenterAt(dType,x,y+1,z+1));

    //on z-end faces:
    if (z == this->getProblemModel()->getZnum())
        return 0.5*(dataCenterAt(dType,x,y+1,z) + dataCenterAt(dType,x+1,y+1,z));

    return 0.25 * (dataCenterAt(dType,x,y+1,z) + dataCenterAt(dType,x,y+1,z+1)
                   + dataCenterAt(dType,x+1,y+1,z) + dataCenterAt(dType,x+1,y+1,z+1));
}

#undef __FUNCT__
#define __FUNCT__ "dataYz"
PetscReal PetscAdLemTaras3D::dataYzAt(std::string dType, PetscInt x, PetscInt y, PetscInt z)
{
    //on x-end faces: same as the values at x-1!
    if (x == this->getProblemModel()->getXnum())
        --x;

    //y-end and z-end corner point, simply return the corner value.
    if ((y == this->getProblemModel()->getXnum())
            && (z == this->getProblemModel()->getZnum())) {
        return dataCenterAt(dType,x+1,y,z);
    }

    //on y-end faces:
    if (y == this->getProblemModel()->getXnum())
        return 0.5*(dataCenterAt(dType,x+1,y,z) + dataCenterAt(dType,x+1,y,z+1));

    //on z-end faces:
    if (z == this->getProblemModel()->getZnum())
        return 0.5*(dataCenterAt(dType,x+1,y,z) + dataCenterAt(dType,x+1,y+1,z));

    return 0.25 * (dataCenterAt(dType,x+1,y,z) + dataCenterAt(dType,x+1,y,z+1)
                   + dataCenterAt(dType,x+1,y+1,z) + dataCenterAt(dType,x+1,y+1,z+1));
}

#undef __FUNCT__
#define __FUNCT__ "muC"
PetscReal PetscAdLemTaras3D::muC(PetscInt x, PetscInt y, PetscInt z)
{
    return dataCenterAt("mu",x,y,z);
}

#undef __FUNCT__
#define __FUNCT__ "lambdaC"
PetscReal PetscAdLemTaras3D::lambdaC(PetscInt x, PetscInt y, PetscInt z)
{
    return dataCenterAt("lambda",x,y,z);
}

#undef __FUNCT__
#define __FUNCT__ "aC"
PetscReal PetscAdLemTaras3D::aC(PetscInt x, PetscInt y, PetscInt z)
{
    return dataCenterAt("atrophy",x,y,z);
}

#undef __FUNCT__
#define __FUNCT_ "muXy"
PetscReal PetscAdLemTaras3D::muXy(PetscInt x, PetscInt y, PetscInt z)
{
    return dataXyAt("mu",x,y,z);
}

#undef __FUNCT_
#define __FUNCT__ "muXz"
PetscReal PetscAdLemTaras3D::muXz(PetscInt x, PetscInt y, PetscInt z)
{
    return dataXzAt("mu",x,y,z);
}

#undef __FUNCT__
#define __FUNCT__ "muYz"
PetscReal PetscAdLemTaras3D::muYz(PetscInt x, PetscInt y, PetscInt z)
{
    return dataYzAt("mu",x,y,z);
}

#undef __FUNCT__
#define __FUNCT_ "lambdaXy"
PetscReal PetscAdLemTaras3D::lambdaXy(PetscInt x, PetscInt y, PetscInt z)
{
    return dataXyAt("lambda",x,y,z);
}

#undef __FUNCT_
#define __FUNCT__ "lambdaXz"
PetscReal PetscAdLemTaras3D::lambdaXz(PetscInt x, PetscInt y, PetscInt z)
{
    return dataXzAt("lambda",x,y,z);
}

#undef __FUNCT__
#define __FUNCT__ "lambdaYz"
PetscReal PetscAdLemTaras3D::lambdaYz(PetscInt x, PetscInt y, PetscInt z)
{
    return dataYzAt("lambda",x,y,z);
}

#undef __FUNCT__
#define __FUNCT__ "solveModel"
PetscErrorCode PetscAdLemTaras3D::solveModel(bool fileToMatlab, const std::string &filename)
{
    PetscErrorCode ierr;
    Vec b,x;
    PetscFunctionBeginUser;
    ierr = DMKSPSetComputeRHS(mDa,computeRHSTaras3D,this);CHKERRQ(ierr);
    ierr = DMKSPSetComputeOperators(mDa,computeMatrixTaras3D,this);CHKERRQ(ierr);
    ierr = KSPSetDM(mKsp,mDa);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(mKsp);CHKERRQ(ierr);
    ierr = KSPSetUp(mKsp);CHKERRQ(ierr);
    ierr = KSPSolve(mKsp,NULL,NULL);CHKERRQ(ierr);
    ierr = KSPGetSolution(mKsp,&x);CHKERRQ(ierr);
    ierr = KSPGetRhs(mKsp,&b);CHKERRQ(ierr);
    if (fileToMatlab) {
//        Mat mat1, mat2;
//        ierr = KSPGetOperators(mKsp,&mat1,&mat2,0);CHKERRQ(ierr);

        PetscViewer viewer;
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename.c_str(),FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
        ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
//        ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);

        ierr = PetscObjectSetName((PetscObject)x,"x");CHKERRQ(ierr);
        ierr = PetscObjectSetName((PetscObject)b,"b");CHKERRQ(ierr);
        ierr = VecView(x,viewer);CHKERRQ(ierr);
//        ierr = VecView(b,viewer);CHKERRQ(ierr);

//        ierr = MatView(mat1,viewer);CHKERRQ(ierr); //doesn't work, probably not supported?
        ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "computeMatrixTaras3D"
PetscErrorCode PetscAdLemTaras3D::computeMatrixTaras3D(
        KSP ksp, Mat J, Mat jac, MatStructure *str, void *ctx)
{
    PetscAdLemTaras3D *user = (PetscAdLemTaras3D*)ctx;

    PetscErrorCode  ierr;
    PetscInt        i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
    PetscReal       Hx,Hy,Hz,HyHzdHx,HxHzdHy,HxHydHz;
    PetscScalar     v[17];
    MatStencil      row, col[17];
    DM              da;
    PetscReal       kBond = 1.0; //need to change it to scale the coefficients.
    PetscReal       kCont = 1.0; //need to change it to scale the coefficients.

    PetscFunctionBeginUser;
    ierr = KSPGetDM(ksp,&da);CHKERRQ(ierr);
    ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
    Hx = 1./(mx-1);
    Hy = 1./(my-1);
    Hz = 1./(mz-1);
    HyHzdHx = (Hy*Hz)/Hx;
    HxHzdHy = (Hx*Hz)/Hy;
    HxHydHz = (Hx*Hy)/Hz;
    ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

    if (user->getProblemModel()->getBcType() != user->getProblemModel()->DIRICHLET)
        SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"only Dirichlet boundary condition supported",0);

    for (k=zs; k<zs+zm; ++k) {
        for (j=ys; j<ys+ym; ++j) {
            for (i=xs; i<xs+xm; ++i) {
                row.i = i; row.j = j; row.k = k;
                // ********************* x-momentum equation *******************
                row.c = 0;
                //Ghost vx unknowns(j=my-1,k=mz-1);boundary vx:(i=0,i=mx-1,j=0,j=my-2,k=0,k=mz-2)
                if (i==0 || i==mx-1 || j==0 || j==my-2 || j==my-1 || k==0 || k==mz-2 || k==mz-1) {
                    //all boundary/ghost conditions use at most two terms and for only vx. So let's
                    //initiate the component for these:
                    col[0].c = 0;       col[1].c = 0;

                    //x-start and x-end faces: vx(i,j,k) = 0
                    if (i==0 || i==mx-1) { //boundary values
                        v[0] = kBond;           col[0].i = i;   col[0].j = j;   col[0].k = k;
                        ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                    } else { //y-start, y-end; z-start, z-end faces:
                        //ghost values for j=my-1, k=mz-1: vx(i,j,k) = 0
                        if (j==my-1 || k==mz-1) {
                            v[0] = kBond;           col[0].i = i;   col[0].j = j;   col[0].k = k;
                            ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                        } else {
                            if (j==0 || j==my-2) {//y-start and y-end faces
                                //3*vx(i,0,k) - vx(i,1,k) = 0; 3*vx(i,my-2,k) - vx(i,my-3,k) = 0
                                v[0] = 3*kBond;         col[0].i = i;   col[0].j = j;   col[0].k=k;
                                if (j==0) {
                                    v[1] = -kBond;      col[1].i = i;   col[1].j = j+1; col[1].k=k;
                                }
                                else if (j==my-2) {
                                    v[1] = -kBond;      col[1].i = i;   col[1].j = j-1; col[1].k=k;
                                }
                                ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                            } else { //z-start and z-end faces
                                //3*vx(i,j,0) - vx(i,j,1) = 0; 3*vx(i,j,mz-2) - vx(i,j,mz-3) = 0
                                v[0] = 3*kBond;         col[0].i = i;   col[0].j = j;   col[0].k=k;
                                if (k==0) {
                                    v[1] = -kBond;      col[1].i = i;   col[1].j = j;   col[1].k=k+1;
                                }
                                else if (k==mz-2) {
                                    v[1] = -kBond;      col[1].i = i;   col[1].j = j;   col[1].k=k-1;
                                }
                                ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                            }
                        }
                    }
                } else { //interior points, x-momentum equation
                    //vx-coefficients, seven terms.
                    for(int ii=0; ii<7; ++ii)
                        col[ii].c = 0;

                    v[0] = 2.0*user->muC(i+1,j+1,k+1)*HyHzdHx;      col[0].i=i+1;   col[0].j=j;     col[0].k=k;
                    v[1] = 2.0*user->muC(i,j+1,k+1)*HyHzdHx;        col[1].i=i-1;   col[1].j=j;     col[1].k=k;
                    v[2] = user->muXy(i,j+1,k)*HxHzdHy;             col[2].i=i;     col[2].j=j+1;   col[2].k=k;
                    v[3] = user->muXy(i,j,k)*HxHzdHy;               col[3].i=i;     col[3].j=j-1;   col[3].k=k;
                    v[4] = user->muXz(i,j,k+1)*HxHydHz;             col[4].i=i;     col[4].j=j;     col[4].k=k+1;
                    v[5] = user->muXz(i,j,k)*HxHydHz;               col[5].i=i;     col[5].j=j;     col[5].k=k-1;
                    v[6] = -v[0]-v[1]-v[2]-v[3]-v[4]-v[5];          col[6].i=i;     col[6].j=j;     col[6].k=k;

                    //vy-coefficients, four terms.
                    for(int ii=7; ii<11; ++ii)
                        col[ii].c = 1;

                    v[7] = user->muXy(i,j+1,k)*Hz;                  col[7].i=i;     col[7].j=j+1;   col[7].k=k;
                    v[8] = -v[7];                                   col[8].i=i-1;   col[8].j=j+1;   col[8].k=k;
                    v[9] = user->muXy(i,j,k)*Hz;                    col[9].i=i-1;   col[9].j=j;     col[9].k=k;
                    v[10] = -v[9];                                  col[10].i=i;    col[10].j=j;    col[10].k=k;

                    //vz-coefficients, four terms.
                    for(int ii=11; ii<15; ++ii)
                        col[ii].c = 2;

                    v[11] = user->muXz(i,j,k+1)*Hy;                 col[11].i=i;    col[11].j=j;    col[11].k=k+1;
                    v[12] = -v[11];                                 col[12].i=i-1;  col[12].j=j;    col[12].k=k+1;
                    v[13] = user->muXz(i,j,k)*Hy;                   col[13].i=i-1;  col[13].j=j;    col[13].k=k;
                    v[14] = -v[13];                                 col[14].i=i;    col[14].j=j;    col[14].k=k;

                    //p-coefficients, two terms.
                    col[15].c = 3;      col[16].c = 3;
                    v[15] = kCont*Hy*Hz;        col[15].i=i;    col[15].j=j+1;  col[15].k=k+1;
                    v[16] = -v[15];             col[16].i=i+1;  col[16].j=j+1;  col[16].k=k+1;

                    ierr=MatSetValuesStencil(jac,1,&row,17,col,v,INSERT_VALUES);
                }
                //*********************** y-momentum equation *******************
                row.c = 1;
                //Ghost vy unknowns(x=mx-1,k=mz-1);boundary vy:(i=0,i=mx-2,j=0,j=my-1,k=0,k=mz-2)
                if (i==0 || i==mx-2 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-2 || k==mz-1) {
                    //all boundary/ghost conditions use at most two terms and for only vy. So let's
                    //initiate the component for these:
                    col[0].c = 1;       col[1].c = 1;

                    //y-start and y-end faces: vy(i,j,k) = 0
                    if (j==0 || j==my-1) { //boundary values
                        v[0] = kBond;           col[0].i = i;   col[0].j = j;   col[0].k = k;
                        ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                    } else { //x-start, x-end; z-start, z-end faces:
                        //ghost values for i=mx-1, k=mz-1: vy(i,j,k) = 0
                        if (i==mx-1 || k==mz-1) {
                            v[0] = kBond;           col[0].i = i;   col[0].j = j;   col[0].k = k;
                            ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                        } else {
                            if (i==0 || i==mx-2) {//x-start and x-end faces
                                //3*vy(0,j,k) - vy(1,j,k) = 0; 3*vy(mx-2,j,k) - vy(mx-3,j,k) = 0
                                v[0] = 3*kBond;         col[0].i = i;   col[0].j = j;   col[0].k=k;
                                if (i==0) {
                                    v[1] = -kBond;      col[1].i = i+1; col[1].j = j;   col[1].k=k;
                                }
                                else if (i==mx-2) {
                                    v[1] = -kBond;      col[1].i = i-1; col[1].j = j;   col[1].k=k;
                                }
                                ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                            } else { //z-start and z-end faces
                                //3*vy(i,j,0) - vy(i,j,1) = 0; 3*vy(i,j,mz-2) - vy(i,j,mz-3) = 0
                                v[0] = 3*kBond;         col[0].i = i;   col[0].j = j;   col[0].k=k;
                                if (k==0) {
                                    v[1] = -kBond;      col[1].i = i;   col[1].j = j;   col[1].k=k+1;
                                }
                                else if (k==mz-2) {
                                    v[1] = -kBond;      col[1].i = i;   col[1].j = j;   col[1].k=k-1;
                                }
                                ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                            }
                        }
                    }
                } else { //interior points, y-momentum equation
                    //vy-coefficients, seven terms.
                    for(int ii=0; ii<7; ++ii)
                        col[ii].c = 1;

                    v[0] = 2.0*user->muC(i+1,j+1,k+1)*HxHzdHy;      col[0].i=i;     col[0].j=j+1;   col[0].k=k;
                    v[1] = 2.0*user->muC(i+1,j,k+1)*HxHzdHy;        col[1].i=i;     col[1].j=j-1;   col[1].k=k;
                    v[2] = user->muXy(i+1,j,k)*HyHzdHx;             col[2].i=i+1;   col[2].j=j;     col[2].k=k;
                    v[3] = user->muXy(i,j,k)*HyHzdHx;               col[3].i=i-1;   col[3].j=j;     col[3].k=k;
                    v[4] = user->muYz(i,j,k+1)*HxHydHz;             col[4].i=i;     col[4].j=j;     col[4].k=k+1;
                    v[5] = user->muYz(i,j,k)*HxHydHz;               col[5].i=i;     col[5].j=j;     col[5].k=k-1;
                    v[6] = -v[0]-v[1]-v[2]-v[3]-v[4]-v[5];          col[6].i=i;     col[6].j=j;     col[6].k=k;

                    //vx-coefficients, four terms.
                    for(int ii=7; ii<11; ++ii)
                        col[ii].c = 0;

                    v[7] = user->muXy(i+1,j,k)*Hz;                  col[7].i=i+1;   col[7].j=j;     col[7].k=k;
                    v[8] = -v[7];                                   col[8].i=i+1;   col[8].j=j-1;   col[8].k=k;
                    v[9] = user->muXy(i,j,k)*Hz;                    col[9].i=i;     col[9].j=j-1;   col[9].k=k;
                    v[10] = -v[9];                                  col[10].i=i;    col[10].j=j;    col[10].k=k;

                    //vz-coefficients, four terms.
                    for(int ii=11; ii<15; ++ii)
                        col[ii].c = 2;

                    v[11] = user->muYz(i,j,k+1)*Hx;                 col[11].i=i;    col[11].j=j;    col[11].k=k+1;
                    v[12] = -v[11];                                 col[12].i=i;    col[12].j=j-1;  col[12].k=k+1;
                    v[13] = user->muYz(i,j,k)*Hx;                   col[13].i=i;    col[13].j=j-1;  col[13].k=k;
                    v[14] = -v[13];                                 col[14].i=i;    col[14].j=j;    col[14].k=k;

                    //p-coefficients, two terms.
                    col[15].c = 3;      col[16].c = 3;
                    v[15] = kCont*Hx*Hz;        col[15].i=i+1;  col[15].j=j;    col[15].k=k+1;
                    v[16] = -v[15];             col[16].i=i+1;  col[16].j=j+1;  col[16].k=k+1;

                    ierr=MatSetValuesStencil(jac,1,&row,17,col,v,INSERT_VALUES);
                }

                //*********************** z-momentum equation *******************
                row.c = 2;
                //Ghost vz unknowns(x=mx-1,y=my-1);boundary vz:(i=0,i=mx-2,j=0,j=my-2,k=0,k=mz-1)
                if (i==0 || i==mx-2 || i==mx-1 || j==0 || j==my-2 || j==my-1 || k==0 || k==mz-1) {
                    //all boundary/ghost conditions use at most two terms and for only vz. So let's
                    //initiate the component for these:
                    col[0].c = 2;       col[1].c = 2;

                    //z-start and z-end faces: vz(i,j,k) = 0
                    if (k==0 || k==mz-1) { //boundary values
                        v[0] = kBond;           col[0].i = i;   col[0].j = j;   col[0].k = k;
                        ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                    } else { //x-start, x-end; y-start, y-end faces:
                        //ghost values for i=mx-1, j=my-1: vz(i,j,k) = 0
                        if (i==mx-1 || j==my-1) {
                            v[0] = kBond;           col[0].i = i;   col[0].j = j;   col[0].k = k;
                            ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                        } else {
                            if (i==0 || i==mx-2) {//x-start and x-end faces
                                //3*vz(0,j,k) - vz(1,j,k) = 0; 3*vz(mx-2,j,k) - vz(mx-3,j,k) = 0
                                v[0] = 3*kBond;         col[0].i = i;   col[0].j = j;   col[0].k=k;
                                if (i==0) {
                                    v[1] = -kBond;      col[1].i = i+1; col[1].j = j;   col[1].k=k;
                                }
                                else if (i==mx-2) {
                                    v[1] = -kBond;      col[1].i = i-1; col[1].j = j;   col[1].k=k;
                                }
                                ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                            } else { //y-start and y-end faces
                                //3*vz(i,0,k) - vz(i,1,k) = 0; 3*vz(i,my-2,k) - vz(i,my-3,k) = 0
                                v[0] = 3*kBond;         col[0].i = i;   col[0].j = j;   col[0].k=k;
                                if (j==0) {
                                    v[1] = -kBond;      col[1].i = i;   col[1].j = j+1; col[1].k=k;
                                }
                                else if (j==my-2) {
                                    v[1] = -kBond;      col[1].i = i;   col[1].j = j-1; col[1].k=k;
                                }
                                ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                            }
                        }
                    }
                } else { //interior points, z-momentum equation
                    //vz-coefficients, seven terms.
                    for(int ii=0; ii<7; ++ii)
                        col[ii].c = 2;

                    v[0] = 2.0*user->muC(i+1,j+1,k+1)*HxHydHz;      col[0].i=i;     col[0].j=j;     col[0].k=k+1;
                    v[1] = 2.0*user->muC(i+1,j+1,k)*HxHydHz;        col[1].i=i;     col[1].j=j;     col[1].k=k-1;
                    v[2] = user->muXz(i+1,j,k)*HyHzdHx;             col[2].i=i+1;   col[2].j=j;     col[2].k=k;
                    v[3] = user->muXz(i,j,k)*HyHzdHx;               col[3].i=i-1;   col[3].j=j;     col[3].k=k;
                    v[4] = user->muYz(i,j+1,k)*HxHzdHy;             col[4].i=i;     col[4].j=j+1;   col[4].k=k;
                    v[5] = user->muYz(i,j,k)*HxHzdHy;               col[5].i=i;     col[5].j=j-1;   col[5].k=k;
                    v[6] = -v[0]-v[1]-v[2]-v[3]-v[4]-v[5];          col[6].i=i;     col[6].j=j;     col[6].k=k;

                    //vx-coefficients, four terms.
                    for(int ii=7; ii<11; ++ii)
                        col[ii].c = 0;

                    v[7] = user->muXz(i+1,j,k)*Hy;                  col[7].i=i+1;   col[7].j=j;     col[7].k=k;
                    v[8] = -v[7];                                   col[8].i=i+1;   col[8].j=j;     col[8].k=k-1;
                    v[9] = user->muXz(i,j,k)*Hy;                    col[9].i=i;     col[9].j=j;     col[9].k=k-1;
                    v[10] = -v[9];                                  col[10].i=i;    col[10].j=j;    col[10].k=k;

                    //vy-coefficients, four terms.
                    for(int ii=11; ii<15; ++ii)
                        col[ii].c = 1;

                    v[11] = user->muYz(i,j+1,k)*Hx;                 col[11].i=i;    col[11].j=j+1;  col[11].k=k;
                    v[12] = -v[11];                                 col[12].i=i;    col[12].j=j+1;  col[12].k=k-1;
                    v[13] = user->muYz(i,j,k)*Hx;                   col[13].i=i;    col[13].j=j;    col[13].k=k-1;
                    v[14] = -v[13];                                 col[14].i=i;    col[14].j=j;    col[14].k=k;

                    //p-coefficients, two terms.
                    col[15].c = 3;      col[16].c = 3;
                    v[15] = kCont*Hx*Hy;        col[15].i=i+1;  col[15].j=j+1;  col[15].k=k;
                    v[16] = -v[15];             col[16].i=i+1;  col[16].j=j+1;  col[16].k=k+1;

                    ierr=MatSetValuesStencil(jac,1,&row,17,col,v,INSERT_VALUES);CHKERRQ(ierr);

                }

                //********************** continuity equation *********************
                row.c = 3;
                if (i==0 || j==0 || k==0 //Ghost values
                        || (i==1 && (j==1 || j==my-1 || k==1 || k==mz-1)) //four edges in x-start face.
                        || (i==mx-1 && (j==1 || j==my-1 || k==1 || k==mz-1)) //four edges in x-end face.
                        || (j==1 && (k==1 || k== mz-1)) //two edges in y-start face.
                        || (j==my-1 && (k==1 || k==mz-1)) //two edges in y-end face
                        //|| (i==2 && j==2 && k==1) //constant pressure point NOT USED. INSTEAD TELL PETSC ABOUT THIS CONSTANT NULL-SPACE PRESSURE
                        ) {//BY DOING PCFIELDSPLIT. FIXME: MIGHT BE BETTER TO USE PCFIELDSPLIT EXPLICITLY HERE IN THE SOLUTION THAN LETTING IT AS
                    //COMMAND LINE OPTION SINCE WE MUST USE PCFIELDSPLIT IN THIS CASE!

                    //For all the ghost and boundary conditions we need at most two terms for p.
                    col[0].c = 3;       col[1].c = 3;
                    if (i==0 || j==0 || k==0) { //Ghost pressure p(i,j,k) = 0;
                        v[0] = kBond;       col[0].i=i; col[0].j=j; col[0].k=k;
                        ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                    } else if (i==1 && (j==1 || j==my-1 || k==1 || k==mz-1)) {//four edges in x-start face.
                        //set dp/dx=0 i.e. p(i+1,j,k) - p(i,j,k) = 0;
                        v[0] = kBond;       col[0].i=i+1;   col[0].j=j;     col[0].k=k;
                        v[1] = -kBond;      col[1].i=i;     col[1].j=j;     col[1].k=k;
                        ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                    } else if (i==mx-1 && (j==1 || j==my-1 || k==1 || k==mz-1)) { //four edges in x-end face.
                        //set dp/dx=0 i.e. p(i,j,k) - p(i-1,j,k) = 0;
                        v[0] = kBond;       col[0].i=i;     col[0].j=j;     col[0].k=k;
                        v[1] = -kBond;      col[1].i=i-1;   col[1].j=j;     col[1].k=k;
                        ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                    } else if (j==1 && (k==1 || k== mz-1)) { //two edges in y-start face.
                        //set dp/dy=0 i.e. p(i,j+1,k) - p(i,j,k) = 0;
                        v[0] = kBond;       col[0].i=i;     col[0].j=j+1;   col[0].k=k;
                        v[1] = -kBond;      col[1].i=i;     col[1].j=j;     col[1].k=k;
                        ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                    } else if (j==my-1 && (k==1 || k==mz-1)) { //two edges in y-end face
                        //set dp/dy=0 i.e. p(i,j,k) - p(i,j-1,k) = 0;
                        v[0] = kBond;       col[0].i=i;     col[0].j=j;     col[0].k=k;
                        v[1] = -kBond;      col[1].i=i;     col[1].j=j-1;   col[1].k=k;
                        ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                    }/* else { //one cell NOTE: RHS needs to be set to kBond*pCell;
                        v[0] = kBond;       col[0].i=i;     col[0].j=j;     col[0].k=k;
                        ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                    }*/
                } else {
                    //vx-coefficients, two terms
                    col[0].c = 0;       col[1].c = 0;
                    v[0] = kCont/Hx;    col[0].i = i;   col[0].j=j-1;   col[0].k=k-1;
                    v[1] = -v[0];       col[1].i = i-1; col[1].j=j-1;   col[1].k=k-1;

                    //vy-coefficients, two terms
                    col[2].c = 1;       col[3].c = 1;
                    v[2] = kCont/Hy;    col[2].i = i-1; col[2].j=j;     col[2].k=k-1;
                    v[3] = -v[2];       col[3].i = i-1; col[3].j=j-1;   col[3].k=k-1;

                    //vz-coefficients, two terms
                    col[4].c = 2;       col[5].c = 2;
                    v[4] = kCont/Hz;    col[4].i = i-1; col[4].j=j-1;   col[4].k=k;
                    v[5] = -v[4];       col[5].i = i-1; col[5].j=j-1;   col[5].k=k-1;

                    ierr=MatSetValuesStencil(jac,1,&row,6,col,v,INSERT_VALUES);CHKERRQ(ierr);
                }
            }
        }
    }
    ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "computeRHSTaras3D"
PetscErrorCode PetscAdLemTaras3D::computeRHSTaras3D(KSP ksp, Vec b, void *ctx)
{
    PetscAdLemTaras3D    *user = (PetscAdLemTaras3D*)ctx;
    PetscErrorCode ierr;
    PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
    PetscScalar    Hx,Hy,Hz;
    PetscAdLemTaras3D::Field    ***rhs;
    DM             da;
    PetscReal       kCont=1.0;

    PetscFunctionBeginUser;
    ierr = KSPGetDM(ksp,&da);CHKERRQ(ierr);
    ierr = DMDAGetInfo(da, 0, &mx, &my, &mz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
    Hx   = 1.0 / (PetscReal)(mx-1);
    Hy   = 1.0 / (PetscReal)(my-1);
    Hz   = 1.0 / (PetscReal)(mz-1);

    ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, b, &rhs);CHKERRQ(ierr);

    if (user->getProblemModel()->getBcType() != user->getProblemModel()->DIRICHLET) {
        SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"only Dirichlet boundary condition implemented",0);
    }
    for (k=zs; k<zs+zm; ++k) {
        for (j=ys; j<ys+ym; ++j) {
            for (i=xs; i<xs+xm; ++i) {
                //Ghost vx unknowns(j=my-1,k=mz-1);boundary vx:(i=0,i=mx-1,j=0,j=my-2,k=0,k=mz-2)
                if (i==0 || i==mx-1 || j==0 || j==my-2 || j==my-1 || k==0 || k==mz-2 || k==mz-1) {
                    rhs[k][j][i].vx = 0;
                } else { //interior points, x-momentum equation
                    rhs[k][j][i].vx = Hy*Hz*(user->muC(i+1,j+1,k+1) + user->muC(i,j+1,k+1)
                                             +user->lambdaC(i+1,j+1,k+1) + user->lambdaC(i,j+1,k+1)
                                             )*(user->aC(i+1,j+1,k+1) - user->aC(i,j+1,k+1))/2.0;
                }
                //*********************** y-momentum equation *******************
                //Ghost vy unknowns(x=mx-1,k=mz-1);boundary vy:(i=0,i=mx-2,j=0,j=my-1,k=0,k=mz-2)
                if (i==0 || i==mx-2 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-2 || k==mz-1) {
                    rhs[k][j][i].vy = 0;
                } else { //interior points, y-momentum equation
                    rhs[k][j][i].vy = Hx*Hz*(user->muC(i+1,j+1,k+1) + user->muC(i+1,j,k+1)
                                             +user->lambdaC(i+1,j+1,k+1) + user->lambdaC(i+1,j,k+1)
                                             )*(user->aC(i+1,j+1,k+1) - user->aC(i+1,j,k+1))/2.0;
                }

                //*********************** z-momentum equation *******************
                //Ghost vz unknowns(x=mx-1,y=my-1);boundary vz:(i=0,i=mx-2,j=0,j=my-2,k=0,k=mz-1)
                if (i==0 || i==mx-2 || i==mx-1 || j==0 || j==my-2 || j==my-1 || k==0 || k==mz-1) {
                    rhs[k][j][i].vz = 0;
                } else { //interior points, z-momentum equation
                    rhs[k][j][i].vz = Hx*Hy*(user->muC(i+1,j+1,k+1) + user->muC(i+1,j+1,k)
                                             +user->lambdaC(i+1,j+1,k+1) + user->lambdaC(i+1,j+1,k)
                                             )*(user->aC(i+1,j+1,k+1) - user->aC(i+1,j+1,k))/2.0;
                }

                //********************** continuity equation *********************
                if (i==0 || j==0 || k==0 //Ghost values
                        || (i==1 && (j==1 || j==my-1 || k==1 || k==mz-1)) //four edges in x-start face.
                        || (i==mx-1 && (j==1 || j==my-1 || k==1 || k==mz-1)) //four edges in x-end face.
                        || (j==1 && (k==1 || k== mz-1)) //two edges in y-start face.
                        || (j==my-1 && (k==1 || k==mz-1))) { //two edges in y-end face
                        rhs[k][j][i].p = 0;
                }/* else if (i==2 && j==1 && k==1) {//constant pressure point
                    rhs[k][j][i].p = kCont*user->getP0Cell();
                }*/ else {
                    rhs[k][j][i].p = kCont*user->aC(i,j,k);
                }
            }
        }
    }

    ierr = DMDAVecRestoreArray(da, b, &rhs);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

