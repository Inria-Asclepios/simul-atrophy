//Driven cavity test with Dirichlet condition, the matrix contains the identity rows
//corresponding to the Dirichlet conditions! i.e. computed domain contains the Dirich-
//let boundaries.
#include "PetscAdLemTaras3D.hxx"
#include "AdLem3D.hxx"

#undef __FUNCT__
#define __FUNCT__ "computeMatrixTaras3dConstantMu"
PetscErrorCode PetscAdLemTaras3D::computeMatrixTaras3dConstantMu(
        KSP ksp, Mat J, Mat jac, MatStructure *str, void *ctx)
{
    PetscAdLemTaras3D *user = (PetscAdLemTaras3D*)ctx;

    PetscErrorCode  ierr;
    PetscInt        i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
    PetscReal       Hx,Hy,Hz,HyHzdHx,HxHzdHy,HxHydHz;
    PetscScalar     v[9];
    MatStencil      row, col[9];
    DM              da;
    PetscReal       kBond = 1.0; //need to change it to scale the coefficients.
    PetscReal       kCont = 1.0; //need to change it to scale the coefficients.

    PetscFunctionBeginUser;
    ierr = KSPGetDM(ksp,&da);CHKERRQ(ierr);
    ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
    Hx = 1;//1./(mx-1);
    Hy = 1;//1./(my-1);
    Hz = 1;//1./(mz-1);
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
                                //                                ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                                //CHANGED:
                                v[0] = kBond;         col[0].i = i;   col[0].j = j;   col[0].k=k;
                                ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                            } else { //z-start and z-end faces
                                //3*vx(i,j,0) - vx(i,j,1) = 0; 3*vx(i,j,mz-2) - vx(i,j,mz-3) = 0
                                v[0] = 3*kBond;         col[0].i = i;   col[0].j = j;   col[0].k=k;
                                if (k==0) {
                                    v[1] = -kBond;      col[1].i = i;   col[1].j = j;   col[1].k=k+1;
                                }
                                else if (k==mz-2) {
                                    v[1] = -kBond;      col[1].i = i;   col[1].j = j;   col[1].k=k-1;
                                }
                                //                                ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                                //CHANGED
                                v[0] = kBond;         col[0].i = i;   col[0].j = j;   col[0].k=k;
                                ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                            }
                        }
                    }
                } else { //interior points, x-momentum equation
                    //vx-coefficients, seven terms.
                    for(int ii=0; ii<7; ++ii)
                        col[ii].c = 0;

                    PetscScalar coeff = (user->muC(i+1,j+1,k+1) + user->muC(i,j+1,k+1))/2.;
                    v[0] = coeff*HyHzdHx;          col[0].i=i+1;   col[0].j=j;     col[0].k=k;
                    v[1] = coeff*HyHzdHx;          col[1].i=i-1;   col[1].j=j;     col[1].k=k;
                    v[2] = coeff*HxHzdHy;          col[2].i=i;     col[2].j=j+1;   col[2].k=k;
                    v[3] = coeff*HxHzdHy;          col[3].i=i;     col[3].j=j-1;   col[3].k=k;
                    v[4] = coeff*HxHydHz;          col[4].i=i;     col[4].j=j;     col[4].k=k+1;
                    v[5] = coeff*HxHydHz;          col[5].i=i;     col[5].j=j;     col[5].k=k-1;
                    v[6] = -2*(v[0]+v[2]+v[4]);    col[6].i=i;     col[6].j=j;     col[6].k=k;

                    //p-coefficients, two terms.
                    col[7].c = 3;      col[8].c = 3;
                    v[7] = kCont*Hy*Hz;        col[7].i=i;    col[7].j=j+1;  col[7].k=k+1;
                    v[8] = -v[7];              col[8].i=i+1;  col[8].j=j+1;  col[8].k=k+1;

                    ierr=MatSetValuesStencil(jac,1,&row,9,col,v,INSERT_VALUES);
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
                                //                                ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                                //CHANGED:
                                v[0] = kBond;         col[0].i = i;   col[0].j = j;   col[0].k=k;
                                ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                            } else { //z-start and z-end faces
                                //3*vy(i,j,0) - vy(i,j,1) = 0; 3*vy(i,j,mz-2) - vy(i,j,mz-3) = 0
                                v[0] = 3*kBond;         col[0].i = i;   col[0].j = j;   col[0].k=k;
                                if (k==0) {
                                    v[1] = -kBond;      col[1].i = i;   col[1].j = j;   col[1].k=k+1;
                                }
                                else if (k==mz-2) {
                                    v[1] = -kBond;      col[1].i = i;   col[1].j = j;   col[1].k=k-1;
                                }
                                //                                ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                                //CHANGED:
                                v[0] = kBond;         col[0].i = i;   col[0].j = j;   col[0].k=k;
                                ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                            }
                        }
                    }
                } else { //interior points, y-momentum equation
                    //vy-coefficients, seven terms.
                    for(int ii=0; ii<7; ++ii)
                        col[ii].c = 1;

                    PetscScalar coeff = (user->muC(i+1,j+1,k+1) + user->muC(i+1,j,k+1))/2.;
                    v[0] = coeff*HyHzdHx;          col[0].i=i+1;   col[0].j=j;     col[0].k=k;
                    v[1] = coeff*HyHzdHx;          col[1].i=i-1;   col[1].j=j;     col[1].k=k;
                    v[2] = coeff*HxHzdHy;          col[2].i=i;     col[2].j=j+1;   col[2].k=k;
                    v[3] = coeff*HxHzdHy;          col[3].i=i;     col[3].j=j-1;   col[3].k=k;
                    v[4] = coeff*HxHydHz;          col[4].i=i;     col[4].j=j;     col[4].k=k+1;
                    v[5] = coeff*HxHydHz;          col[5].i=i;     col[5].j=j;     col[5].k=k-1;
                    v[6] = -2*(v[0]+v[2]+v[4]);    col[6].i=i;     col[6].j=j;     col[6].k=k;

                    //p-coefficients, two terms.
                    col[7].c = 3;      col[8].c = 3;
                    v[7] = kCont*Hx*Hz;       col[7].i=i+1;  col[7].j=j;    col[7].k=k+1;
                    v[8] = -v[7];             col[8].i=i+1;  col[8].j=j+1;  col[8].k=k+1;

                    ierr=MatSetValuesStencil(jac,1,&row,9,col,v,INSERT_VALUES);
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
                                //                                ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                                //CHANGED
                                v[0] = kBond;         col[0].i = i;   col[0].j = j;   col[0].k=k;
                                ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                            } else { //y-start and y-end faces
                                //3*vz(i,0,k) - vz(i,1,k) = 0; 3*vz(i,my-2,k) - vz(i,my-3,k) = 0
                                v[0] = 3*kBond;         col[0].i = i;   col[0].j = j;   col[0].k=k;
                                if (j==0) {
                                    v[1] = -kBond;      col[1].i = i;   col[1].j = j+1; col[1].k=k;
                                }
                                else if (j==my-2) {
                                    v[1] = -kBond;      col[1].i = i;   col[1].j = j-1; col[1].k=k;
                                }
                                //                                ierr=MatSetValuesStencil(jac,1,&row,2,col,v,INSERT_VALUES);CHKERRQ(ierr);
                                //CHANGED
                                v[0] = 3*kBond;         col[0].i = i;   col[0].j = j;   col[0].k=k;
                                ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                            }
                        }
                    }
                } else { //interior points, z-momentum equation
                    //vz-coefficients, seven terms.
                    for(int ii=0; ii<7; ++ii)
                        col[ii].c = 2;

                    PetscScalar coeff = (user->muC(i+1,j+1,k+1) + user->muC(i+1,j+1,k))/2.;
                    v[0] = coeff*HyHzdHx;          col[0].i=i+1;   col[0].j=j;     col[0].k=k;
                    v[1] = coeff*HyHzdHx;          col[1].i=i-1;   col[1].j=j;     col[1].k=k;
                    v[2] = coeff*HxHzdHy;          col[2].i=i;     col[2].j=j+1;   col[2].k=k;
                    v[3] = coeff*HxHzdHy;          col[3].i=i;     col[3].j=j-1;   col[3].k=k;
                    v[4] = coeff*HxHydHz;          col[4].i=i;     col[4].j=j;     col[4].k=k+1;
                    v[5] = coeff*HxHydHz;          col[5].i=i;     col[5].j=j;     col[5].k=k-1;
                    v[6] = -2*(v[0]+v[2]+v[4]);    col[6].i=i;     col[6].j=j;     col[6].k=k;

                    //p-coefficients, two terms.
                    col[7].c = 3;      col[8].c = 3;
                    v[7] = kCont*Hx*Hy;        col[7].i=i+1;  col[7].j=j+1;  col[7].k=k;
                    v[8] = -v[7];              col[8].i=i+1;  col[8].j=j+1;  col[8].k=k+1;

                    ierr=MatSetValuesStencil(jac,1,&row,9,col,v,INSERT_VALUES);CHKERRQ(ierr);

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
                    } //else { //one cell NOTE: RHS needs to be set to kBond*pCell;
                    //v[0] = kBond;       col[0].i=i;     col[0].j=j;     col[0].k=k;
                    //ierr=MatSetValuesStencil(jac,1,&row,1,col,v,INSERT_VALUES);CHKERRQ(ierr);
                    //}
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
#define __FUNCT__ "computeRHSTaras3dConstantMu"
PetscErrorCode PetscAdLemTaras3D::computeRHSTaras3dConstantMu(KSP ksp, Vec b, void *ctx)
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
    Hx   = 1;//1.0 / (PetscReal)(mx-1);
    Hy   = 1;//1.0 / (PetscReal)(my-1);
    Hz   = 1;//1.0 / (PetscReal)(mz-1);

    ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, b, &rhs);CHKERRQ(ierr);

    if (user->getProblemModel()->getBcType() != user->getProblemModel()->DIRICHLET) {
        SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"only Dirichlet boundary condition implemented",0);
    }
    for (k=zs; k<zs+zm; ++k) {
        for (j=ys; j<ys+ym; ++j) {
            for (i=xs; i<xs+xm; ++i) {
                //Ghost vx unknowns(j=my-1,k=mz-1);boundary vx:(i=0,i=mx-1,j=0,j=my-2,k=0,k=mz-2)
                if (i==0 || i==mx-1 || j==0 || j==my-2 || j==my-1 || k==0 || k==mz-1){
                    rhs[k][j][i].vx = 0;
                } else if ( k==mz-2){
                    rhs[k][j][i].vx = 30;
                } else { //interior points, x-momentum equation
                    rhs[k][j][i].vx = 0;
                }
                // *********************** y-momentum equation *******************
                //Ghost vy unknowns(x=mx-1,k=mz-1);boundary vy:(i=0,i=mx-2,j=0,j=my-1,k=0,k=mz-2)
                if (i==0 || i==mx-2 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-2 || k==mz-1) {
                    rhs[k][j][i].vy = 0;
                } else { //interior points, y-momentum equation
                    rhs[k][j][i].vy = 0;
                }

                // *********************** z-momentum equation *******************
                //Ghost vz unknowns(x=mx-1,y=my-1);boundary vz:(i=0,i=mx-2,j=0,j=my-2,k=0,k=mz-1)
                if (i==0 || i==mx-2 || i==mx-1 || j==0 || j==my-2 || j==my-1 || k==0 || k==mz-1) {
                    rhs[k][j][i].vz = 0;
                } else { //interior points, z-momentum equation
                    rhs[k][j][i].vz = 0;
                }

                //  ********************** continuity equation *********************
                if (i==0 || j==0 || k==0 //Ghost values
                        || (i==1 && (j==1 || j==my-1 || k==1 || k==mz-1)) //four edges in x-start face.
                        || (i==mx-1 && (j==1 || j==my-1 || k==1 || k==mz-1)) //four edges in x-end face.
                        || (j==1 && (k==1 || k== mz-1)) //two edges in y-start face.
                        || (j==my-1 && (k==1 || k==mz-1))) { //two edges in y-end face
                    rhs[k][j][i].p = 0;
                } //else if (i==2 && j==1 && k==1) {//constant pressure point
                //  rhs[k][j][i].p = kCont*user->getP0Cell();
                //}
                else {
                    rhs[k][j][i].p = 0;//-kCont*user->aC(i,j,k);
                }
            }
        }
    }

    ierr = DMDAVecRestoreArray(da, b, &rhs);CHKERRQ(ierr);
    //    ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
    //    ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
    //        ierr = MatNullSpaceRemove(user->getNullSpace(),b,NULL);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

