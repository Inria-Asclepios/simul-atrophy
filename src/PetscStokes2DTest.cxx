//la yaha Taras Geyako test garne, 2D!
static char help[] = "Demonstrates using DM to manage grid and create the matrix for finite difference method\n\n";

/*T
   Concepts: DMDA
   Concepts: DMDACreateMatrix
   Processors: n
T*/
#include <petscsys.h>
#include <petscdmda.h>
//#include <petscbag.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscViewer    viewer;
  DM             da;
//  Vec            global,local;
//  PetscMPIInt    rank;
  PetscInt       m=3,n=4;

  ierr = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr);

  /* Create a DMDA and an associated vector */
  /*ierr = DMDACreate2d(PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE,DMDA_STENCIL_STAR,m,n,
                      PETSC_DECIDE,PETSC_DECIDE,2,1,NULL,NULL,&da);CHKERRQ(ierr);*/
  ierr = DMDACreate2d(PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE,DMDA_STENCIL_BOX,m,n,
                      PETSC_DECIDE,PETSC_DECIDE,2,1,NULL,NULL,&da);CHKERRQ(ierr);


  /*Bish: Let us add a matrix */
  Mat mat1;
  ierr = DMSetMatrixPreallocateOnly(da,PETSC_TRUE);
  ierr = DMCreateMatrix(da,MATMPIAIJ,&mat1); CHKERRQ(ierr);
//  MatInfo info1;
//  MatGetInfo(mat1,MAT_GLOBAL_MAX,&info1);
//  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"nnz: %d",info1.nz_used);
  MatStencil row, col[4]; //let's test 4 non-zeros in one row, without p!
  PetscScalar val[4];
  PetscInt i,j;
  PetscScalar mu_c = 10, mu_g = 10;
  PetscScalar hx = 1, hy = 1;
  PetscScalar one = 1.;
  DMDALocalInfo info;
  enum boundary_condition{FREE_SLIP, NO_SLIP};
  boundary_condition btype1 = NO_SLIP;

  ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);

  for(j = info.ys; j < info.ys+info.ym; ++j){
      for(i = info.xs; i < info.xs+info.xm; ++i){
          row.i = i; row.j = j;
          if(btype1 == NO_SLIP){ //just to prevent it doing many times within boundary condition code
              col[0].i = i; col[0].j = j;
              val[0] = 1;
          }
          if (j == 0 || j == info.my-1){ //left and right borders
              if(btype1 == NO_SLIP){
                  //vx(i,j) = 0;
                  row.c = 0;
                  col[0].c = 0;
                  ierr = MatSetValuesStencil(mat1,1,&row,1,col,&one,INSERT_VALUES);CHKERRQ(ierr);

                  //vy: //vy(i,j) - vy(i,j+-1)/(3*hx) = 0;
                  row.c = 1;
                  col[0].c = 1; col[1].c = 1;
                  col[1].i = i;
                  if(j == 0) //vy(i,j) - vy(i,j+1)/(3*hx) = 0;
                      col[1].j = j+1;
                  else //vy(i,j) - vy(i,j-1)/(3*hx) = 0;
                      col[1].j = j-1;
                  val[1] = -1./3./hy;

                  ierr = MatSetValuesStencil(mat1,1,&row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);
              }
              else{
                  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"only no slip boundary condition implemented so far, use that!"); CHKERRQ(ierr);
                  //use some exception to return from here!!
              }
          }
          else if (i == 0 || i == info.mx-1){  //top and bottom borders
              if(btype1 == NO_SLIP){
                  //vx: vx(i,j) - vx(i+-1,j)/(3*hy) = 0;
                  row.c = 0;
                  col[0].c = 0; col[1].c = 0;
                  val[1] = -1./3./hy;
                  if (i == 0) //vx(i,j) - vx(i+1,j)/(3*hy) = 0;
                      col[1].i = i+1;
                  else //vx(i,j) - vx(i-1,j)/(3*hy) = 0;
                      col[1].i = i-1;
                  col[1].j = j;
                  ierr = MatSetValuesStencil(mat1,1,&row,2,col,val,INSERT_VALUES);CHKERRQ(ierr);

                  //vy(i,j) = 0;
                  row.c = 1;
                  col[0].c = 1;
                  ierr = MatSetValuesStencil(mat1,1,&row,1,col,val,INSERT_VALUES);CHKERRQ(ierr);
              }
              else {
                  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"only no slip boundary condition implemented so far, use that!"); CHKERRQ(ierr);
              }
          }
          else { //Interior points:
              row.c = 0;//x-stokes
              col[0].c = 0; col[0].i = i; col[0].j = j+1;
              val[0] = 2*mu_c/hx/hx;

              col[1].c = 0; col[1].i = i; col[1].j = j;
              val[1] = -val[0] - mu_g/hy/hy;

              col[2].c = 1; col[2].i = i+1; col[2].j = j;
              val[2] = mu_g/hx/hy;

              col[3].c = 1; col[3].i = i; col[3].j = j;
              val[3] = mu_g/hx/hy;

              ierr = MatSetValuesStencil(mat1,1,&row,4,col,val,INSERT_VALUES);CHKERRQ(ierr);

              row.c = 1; //y-stokes
              col[0].c = 1; col[0].i = i; col[0].j = j+1;
              val[0] = 2*mu_c/hx/hx;

              col[1].c = 1; col[1].i = i; col[1].j = j;
              val[1] = -val[0] - mu_g/hy/hy;

              col[2].c = 0; col[2].i = i+1; col[2].j = j;
              val[2] = mu_g/hx/hy;

              col[3].c = 0; col[3].i = i; col[3].j = j;
              val[3] = mu_g/hx/hy;

              ierr = MatSetValuesStencil(mat1,1,&row,4,col,val,INSERT_VALUES);CHKERRQ(ierr);


          }
      }
  }

  MatAssemblyBegin(mat1,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mat1,MAT_FINAL_ASSEMBLY);
    ierr = MatSetOption(mat1,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)mat1,"mat1");CHKERRQ(ierr);
//  MatView(mat1,PETSC_VIEWER_STDOUT_WORLD);

//  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,params->filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"output_file",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);

    ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
//  ierr = DMDASetFieldName(da,0,"field1");CHKERRQ(ierr);
//  ierr = DMDASetFieldName(da,1,"field2");CHKERRQ(ierr);

//    ierr = MatView(mat1,viewer);CHKERRQ(ierr);
    ierr = MatView(mat1,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);

  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  /* clean up and exit */
  ierr = DMDestroy(&da);CHKERRQ(ierr);
//  ierr = VecDestroy(&local);CHKERRQ(ierr);
//  ierr = VecDestroy(&global);CHKERRQ(ierr);
  ierr = MatDestroy(&mat1);
  ierr = PetscFinalize();
  return 0;
}


