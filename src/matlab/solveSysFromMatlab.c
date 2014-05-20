
static char help[] = "Reads in a matrix in ASCII MATLAB format (I,J,A), read in vectors rhs and (optionally)exact_solu in ASCII format.\n\
Solves the linear system and writes back the solution in ASCII Matlab form as petsC_sol.m\
Note: I and J start at 1, not 0, use -noshift if indices in file start with zero!\n\
Input parameters are:\n\
  -Ain  <filename> : input matrix in ascii format\n\
  -rhs  <filename> : input rhs in ascii format\n\
  -solu  <filename> : input true solution in ascii format\n\\n";

/*
Example: ./solveSysFromMatlab.c -Ain Ain -rhs rhs -solu solu -noshift -mat_view
 with the datafiles in the followig format:
Ain (I and J start at 0):
------------------------
3 3 6
0 0 1.0
0 1 2.0
1 0 3.0
1 1 4.0
1 2 5.0
2 2 6.0

rhs
---
0 3.0
1 12.0
2 6.0

solu
----
0 1.0
0 1.0
0 1.0
*/

/* 
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/
#include <petscksp.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Mat            A;
  Vec            x,b,u,b_hat;
  char           Ain[PETSC_MAX_PATH_LEN],rhs[PETSC_MAX_PATH_LEN],solu[PETSC_MAX_PATH_LEN]; 
  PetscErrorCode ierr;
  int            m,n,nz; /* these are fscaned so kept as int */
  PetscInt       i,col,row,shift = 1,sizes[3],nsizes;
  PetscScalar    val;
  //  PetscReal      res_norm;
  FILE           *Afile,*bfile,*ufile;
  PetscViewer    view;
  PetscBool      flg_A,flg_b,flg_u,flg;
  PetscMPIInt    size;

  KSP            ksp;         /* linear solver context */
  PC             pc;           /* preconditioner context */
  PetscReal      norm,tol=1.e-14;  /* norm of solution error */
  PetscInt       its;
  PetscScalar    neg_one = -1.0;
  PetscBool      nonzeroguess = PETSC_FALSE;

  PetscInitialize(&argc,&args,(char *)0,help);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size != 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"This is a uniprocessor example only!");

  /* Read in matrix, rhs and exact solution from ascii files */
  ierr = PetscOptionsGetString(PETSC_NULL,"-Ain",Ain,PETSC_MAX_PATH_LEN,&flg_A);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL,"-noshift",&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-nonzero_guess",&nonzeroguess,PETSC_NULL);CHKERRQ(ierr);

  if (flg) shift = 0;
  if (flg_A){
    ierr = PetscPrintf(PETSC_COMM_SELF,"\n Read matrix in ascii format ...\n");CHKERRQ(ierr);
    ierr = PetscFOpen(PETSC_COMM_SELF,Ain,"r",&Afile);CHKERRQ(ierr); 
    nsizes = 3;
    ierr = PetscOptionsGetIntArray(PETSC_NULL,"-nosizesinfile",sizes,&nsizes,&flg);CHKERRQ(ierr);
    if (flg) {
      if (nsizes != 3) SETERRQ(PETSC_COMM_WORLD,1,"Must pass in three m,n,nz as arguments for -nosizesinfile");
      m = sizes[0];
      n = sizes[1];
      nz = sizes[2];
    } else {
      fscanf(Afile,"%d %d %d\n",&m,&n,&nz);
    }
    ierr = PetscPrintf(PETSC_COMM_SELF,"m: %d, n: %d, nz: %d \n", m,n,nz);CHKERRQ(ierr);
    if (m != n) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_SIZ, "Number of rows, cols must be same for this example\n");
    ierr = MatCreate(PETSC_COMM_SELF,&A);CHKERRQ(ierr);
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    //    ierr = MatSeqAIJSetPreallocation(A,nz/m,PETSC_NULL);CHKERRQ(ierr);  
    /* The above code did not work because the rows do not have same  number of non-zeros  */
    ierr = MatSeqAIJSetPreallocation(A,PETSC_DEFAULT,PETSC_NULL);CHKERRQ(ierr);

    for (i=0; i<nz; i++) {
      fscanf(Afile,"%d %d %le\n",&row,&col,(double*)&val);
      row -= shift; col -= shift;  /* set index set starts at 0 */
      ierr = MatSetValues(A,1,&row,1,&col,&val,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    fflush(stdout);
    fclose(Afile);  
  }

  /* Display the matrix read */
  /*   PetscPrintf(PETSC_COMM_WORLD,"Read matrix:\n");
       MatView(A,PETSC_VIEWER_STDOUT_WORLD);*/

  ierr = PetscOptionsGetString(PETSC_NULL,"-rhs",rhs,PETSC_MAX_PATH_LEN,&flg_b);CHKERRQ(ierr);
  if (flg_b){
    ierr = VecCreate(PETSC_COMM_SELF,&b);CHKERRQ(ierr);
    ierr = VecSetSizes(b,PETSC_DECIDE,n);CHKERRQ(ierr);
    ierr = VecSetFromOptions(b);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF,"\n Read rhs in ascii format ...\n");CHKERRQ(ierr);
    ierr = PetscFOpen(PETSC_COMM_SELF,rhs,"r",&bfile);CHKERRQ(ierr); 
    for (i=0; i<n; i++) {      
      fscanf(bfile,"%le\n",(double*)&val); 
      ierr = VecSetValues(b,1,&i,&val,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
    fflush(stdout);
    fclose(bfile);
  }

  ierr = PetscOptionsGetString(PETSC_NULL,"-solu",solu,PETSC_MAX_PATH_LEN,&flg_u);CHKERRQ(ierr);
  if (flg_u){
    ierr = VecCreate(PETSC_COMM_SELF,&u);CHKERRQ(ierr);
    ierr = VecSetSizes(u,PETSC_DECIDE,n);CHKERRQ(ierr);
    ierr = VecSetFromOptions(u);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF,"\n Read exact solution in ascii format ...\n");CHKERRQ(ierr);
    ierr = PetscFOpen(PETSC_COMM_SELF,solu,"r",&ufile);CHKERRQ(ierr); 
    for (i=0; i<n; i++) {
      fscanf(ufile,"%le\n",(double*)&val); 
      ierr = VecSetValues(u,1,&i,&val,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(u);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(u);CHKERRQ(ierr);
    fflush(stdout);
    fclose(ufile);
  }
 
  /* Write matrix, rhs and exact solution in Petsc binary file */
  /*ierr = PetscPrintf(PETSC_COMM_SELF,"\n Write matrix in binary to 'matrix.dat' ...\n");CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"matrix.dat",FILE_MODE_WRITE,&view);CHKERRQ(ierr);
  ierr = MatView(A,view);CHKERRQ(ierr);*/



  // if (flg_b){ /* Write rhs in Petsc binary file */
  /*ierr = PetscPrintf(PETSC_COMM_SELF,"\n Write rhs in binary to 'matrix.dat' ...\n");CHKERRQ(ierr);
    ierr = VecView(b,view);CHKERRQ(ierr);
  }

  if (flg_u){
    ierr = PetscPrintf(PETSC_COMM_SELF,"\n Write exact solution in binary to 'matrix.dat' ...\n");CHKERRQ(ierr);
    ierr = VecView(u,view);CHKERRQ(ierr);
    }*/

  /* Create the solution vector */
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);

  



  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* 
     Create linear solver context
  */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  /* 
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);

  /* 
     Set linear solver defaults for this problem (optional).
     - By extracting the KSP and PC contexts from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
     - The following four statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions();
  */
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    //ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCILU);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
    //    ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);  
  /*  Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  if (nonzeroguess) {
    PetscScalar p = .5;
    ierr = VecSet(x,p);CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
  }
 
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* 
     Solve linear system
  */
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr); 

  /* 
     View solver info; we could instead use the option -ksp_view to
     print this info to the screen at the conclusion of KSPSolve().
  */
  ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);


  /* Write the solution to the Matlab file*/
   PetscViewer view_out;
   ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"petsc_sol.m",&view_out);CHKERRQ(ierr);
   PetscViewerSetFormat(view_out,PETSC_VIEWER_ASCII_MATLAB);
   ierr = VecView(x,view_out);CHKERRQ(ierr);
   ierr = PetscViewerDestroy(&view_out);CHKERRQ(ierr);
   /*
   ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"petsc_b.m",&view_out);CHKERRQ(ierr);
   PetscViewerSetFormat(view_out,PETSC_VIEWER_ASCII_MATLAB);
   ierr = VecView(b,view_out);CHKERRQ(ierr);
   ierr = PetscViewerDestroy(&view_out);CHKERRQ(ierr);

   ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"petsc_A.m",&view_out);CHKERRQ(ierr);
   PetscViewerSetFormat(view_out,PETSC_VIEWER_ASCII_MATLAB);
   ierr = MatView(A,view_out);CHKERRQ(ierr);
   ierr = PetscViewerDestroy(&view_out);CHKERRQ(ierr);
   */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      Check solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* 
     Check the error
  */
   /*Initialize b_hat vector*/
  ierr = VecCreate(PETSC_COMM_WORLD,&b_hat);CHKERRQ(ierr);
  ierr = VecSetSizes(b_hat,PETSC_DECIDE,m);CHKERRQ(ierr);
  ierr = VecSetFromOptions(b_hat);CHKERRQ(ierr);


  ierr = MatMult(A,x,b_hat);CHKERRQ(ierr);
  ierr = VecAXPY(b,neg_one,b_hat);CHKERRQ(ierr);
  ierr  = VecNorm(b,NORM_2,&norm);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  if (norm > tol){
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %G, Iterations %D\n",
                     norm,its);CHKERRQ(ierr);
  }




  /* 
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */



  ierr = MatDestroy(&A);CHKERRQ(ierr);
  if (flg_b) {ierr = VecDestroy(&b);CHKERRQ(ierr);}
  if (flg_u) {ierr = VecDestroy(&u);CHKERRQ(ierr);}
  ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&b_hat);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return 0;
}

