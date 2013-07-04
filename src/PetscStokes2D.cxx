//la yaha Taras Geyako test garne, 2D!
static char help[] = "Demonstrates using PetscViewerSetFormat(viewer,PETSC_FORMAT_BINARY_MATLAB)\n\n";

/*T
   Concepts: viewers
   Concepts: bags
   Processors: n
T*/
#include <petscsys.h>
#include <petscdmda.h>
#include <petscbag.h>

typedef struct {
  char      filename[PETSC_MAX_PATH_LEN];
  PetscReal ra;
  PetscInt  ia;
  PetscBool ta;
} Parameter;

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscBag       bag;
  Parameter      *params;
  PetscViewer    viewer;
  DM             da;
  Vec            global,local;
  PetscMPIInt    rank;
  PetscInt       m=10,n=10;

  /*
    Every PETSc routine should begin with the PetscInitialize() routine.
    argc, argv - These command line arguments are taken to extract the options
                 supplied to PETSc and options supplied to MPI.
    help       - When PETSc executable is invoked with the option -help,
                 it prints the various options that can be applied at
                 runtime.  The user can use the "help" variable place
                 additional help messages in this printout.
  */
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr);

  /* Create a DMDA and an associated vector */
  /*ierr = DMDACreate2d(PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE,DMDA_STENCIL_STAR,m,n,
                      PETSC_DECIDE,PETSC_DECIDE,2,1,NULL,NULL,&da);CHKERRQ(ierr);*/
  ierr = DMDACreate2d(PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_NONE,DMDA_STENCIL_BOX,m,n,
                      PETSC_DECIDE,PETSC_DECIDE,2,1,NULL,NULL,&da);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(da,&global);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(da,&local);CHKERRQ(ierr);
  ierr = VecSet(global,-1.0);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,global,INSERT_VALUES,local);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,global,INSERT_VALUES,local);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = VecScale(local,rank+1);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(da,local,ADD_VALUES,global);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(da,local,ADD_VALUES,global);CHKERRQ(ierr);
//  ierr = DMLocalToGlobalBegin(da,local,INSERT_VALUES,global);CHKERRQ(ierr);
//  ierr = DMLocalToGlobalEnd(da,local,INSERT_VALUES,global);CHKERRQ(ierr);

  /* Create an empty bag */
  ierr = PetscBagCreate(PETSC_COMM_WORLD,sizeof(Parameter),&bag);CHKERRQ(ierr);
  ierr = PetscBagGetData(bag,(void**)&params);CHKERRQ(ierr);

  /* fill bag: register variables, defaults, names, help strings */
  ierr = PetscBagSetName(bag,"ParameterBag","contains problem parameters");CHKERRQ(ierr);
  ierr = PetscBagRegisterString(bag,&params->filename,PETSC_MAX_PATH_LEN,"output_file","filename","Name of secret file");CHKERRQ(ierr);
  ierr = PetscBagRegisterReal  (bag,&params->ra,1.0,"param_1","The first parameter");CHKERRQ(ierr);
  ierr = PetscBagRegisterInt   (bag,&params->ia,5,"param_2","The second parameter");CHKERRQ(ierr);
  ierr = PetscBagRegisterBool (bag,&params->ta,PETSC_TRUE,"do_output","Write output file (true/false)");CHKERRQ(ierr);

  /*
     Write output file with PETSC_VIEWER_BINARY_MATLAB format
     NOTE: the output generated with this viewer can be loaded into
     MATLAB using $PETSC_DIR/bin/matlab/PetscReadBinaryMatlab.m
  */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,params->filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
  ierr = PetscBagView(bag,viewer);CHKERRQ(ierr);
  ierr = DMDASetFieldName(da,0,"field1");CHKERRQ(ierr);
  ierr = DMDASetFieldName(da,1,"field2");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)global,"da1");CHKERRQ(ierr);
  ierr = VecView(global,viewer);CHKERRQ(ierr);

  /*Bish: Let us add a matrix */
  Mat mat1;
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
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"let's see: %d, %d",info.mx,info.my);

  for(i = info.ys; i < info.ys+info.ym-1; ++i){
      for(j = info.xs; j < info.xs+info.xm-1; ++j){
          row.i = i; row.j = j;
          if(btype1 == NO_SLIP){ //just to prevent it doing many times within boundary condition code
              col[0].i = i; col[0].j = j;
              val[0] = 1;
          }
          if (j == 0 || j == info.mx-1){ //left and right borders
              if(btype1 == NO_SLIP){
                  //vx(i,j) = 0;
                  row.c = 0;
                  col[0].c = 0;
                  MatSetValuesStencil(mat1,1,&row,1,col,&one,INSERT_VALUES);

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
          else if (i == 0 || i == info.my-1){  //top and bottom borders
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

              ierr = MatSetValuesStencil(mat1,1,&row,4,col,val,INSERT_VALUES);

              row.c = 1; //y-stokes
              col[0].c = 1; col[0].i = i; col[0].j = j+1;
              val[0] = 2*mu_c/hx/hx;

              col[1].c = 1; col[1].i = i; col[1].j = j;
              val[1] = -val[0] - mu_g/hy/hy;

              col[2].c = 0; col[2].i = i+1; col[2].j = j;
              val[2] = mu_g/hx/hy;

              col[3].c = 0; col[3].i = i; col[3].j = j;
              val[3] = mu_g/hx/hy;

              ierr = MatSetValuesStencil(mat1,1,&row,4,col,val,INSERT_VALUES);


          }
      }
  }

  MatAssemblyBegin(mat1,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mat1,MAT_FINAL_ASSEMBLY);
    ierr = MatSetOption(mat1,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)mat1,"mat1");CHKERRQ(ierr);
//  ierr = MatView(mat1,viewer);CHKERRQ(ierr);

  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  /* clean up and exit */
  ierr = PetscBagDestroy(&bag);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);
  ierr = VecDestroy(&local);CHKERRQ(ierr);
  ierr = VecDestroy(&global);CHKERRQ(ierr);
  ierr = MatDestroy(&mat1);
  ierr = PetscFinalize();
  return 0;
}

