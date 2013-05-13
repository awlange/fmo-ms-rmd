// AWLG: copied for convenience. Thanks Scot!
/*
  This file has my implementation of the LAPACK routine dsyev
  for C++.  This program solves for the eigenvalues and, if
  desired, the eigenvectors for a square symmetric matrix H.
  It assumes that the upper triangle of H is stored.

  There are two function calls defined in this header, of the
  forms

    void dsyev(double **H, int n, double *E)
    void dsyev(double **H, int n, double *E, double **Evecs)

    H: the n by n symmetric matrix that we are solving, with the
       upper triangle stored
    n: the order of the square matrix H
    E: an n-element array to hold the eigenvalues of H
    Evecs: an n by n matrix to hold the eigenvectors of H, if
           they are requested.

  The function call is defined twice, so that whether or not
  eigenvectors are called for the proper version is called.

  Scot Shaw
  30 August 1999
*/

#ifndef FMR_DSYEV_H
#define FMR_DSYEV_H

#include <math.h>

void dsyev(double **H, int n, double *E);
void dsyev(double **H, int n, double *E, double **Evecs);

double *dsyev_ctof(double **in, int rows, int cols);
double **dsyev_ftoc(double *in, int rows, int cols);
void dsyev_ftoc2(double *in, double **out, int rows, int cols);
void dsyev_sort(double *E, double **Evecs, int N);

extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a,
		       int *lda, double *w, double *work, int *lwork,
		       int *info);


void dsyev(double **H, int n, double *E)
{
  char jobz, uplo;
  int lda, lwork, info, i;
  double *a, *work, **Evecs;

  jobz = 'N'; /* V/N indicates that eigenvectors should/should not
		 be calculated. */

  uplo = 'U'; /* U/L indicated that the upper/lower triangle of the
		 symmetric matrix is stored. */

  lda = n; // The leading dimension of the matrix to be solved.
  
  lwork = 3*n-1;
  work = new double[lwork]; /* The work array to be used by dsyev and
			       its size. */

  a = dsyev_ctof(H, n, lda); /* Here we convert the incoming array, assumed
			  to be in double index C form, to a Fortran
			  style matrix. */

  dsyev_(&jobz, &uplo, &n, a, &lda, E, work, &lwork, &info);

  Evecs = dsyev_ftoc(a, n, lda);
//  AWGL: No sorting! Makes indexing messy
//  dsyev_sort(E, Evecs, n); /* Though it is just a dummy array, we convert
//  			the returned values for eigenvectors from
//			Fortran single-pointer form to C double-
//			pointer form.  Then we sort the results by
//			eigenvalue in increasing order. */
//
  delete a;
  delete work;
  delete Evecs;
}


void dsyev(double **H, int n, double *E, double **Evecs)
{
  char jobz, uplo;
  int lda, lwork, info, i;
  double *a, *work;

  jobz = 'V';
  uplo = 'U';
  lda = n;
  
  lwork = 3*n-1;
  work = new double[lwork];

  a = dsyev_ctof(H, n, lda);

  dsyev_(&jobz, &uplo, &n, a, &lda, E, work, &lwork, &info);

  dsyev_ftoc2(a, Evecs, n, lda);
  dsyev_sort(E, Evecs, n);

  delete a;
  delete work;
}


double* dsyev_ctof(double **in, int rows, int cols)
{
  double *out;
  int i, j;

  out = new double[rows*cols];
  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*cols] = in[i][j];
  return(out);
}


double** dsyev_ftoc(double *in, int rows, int cols)
{
  double **out;
  int i, j;

  out = new double*[rows];
  for (i=0; i<rows; i++) out[i] = new double[cols];
  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*cols];
  return(out);
}


void dsyev_ftoc2(double *in, double **out, int rows, int cols)
{
  int i, j;

  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*cols];
}


void dsyev_sort(double *E, double **Evecs, int N)
{
  double temp;
  int i, j, k;

  // *** AWGL: modify so that sorted by minimum, not absolute value ***
  for (j=0; j<N; j++) for (i=0; i<N-1; i++)
    if ( E[i] >  E[i+1] ) {
      temp = E[i]; E[i] = E[i+1]; E[i+1] = temp;
      
      for (k=0; k<N; k++) {
	temp = Evecs[k][i];
	Evecs[k][i] = Evecs[k][i+1];
	Evecs[k][i+1] = temp;
      }
    }
}

#endif
