#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R_ext/Utils.h>
#include <Rversion.h>
#include <Rconfig.h>
#include <R_ext/Constants.h>
#include <R_ext/Random.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

void inverse(double *A, int N);
double matrixDet(double *A, int N, int LDA, int doLog);
double quadform(double *x, double *A, int N, int incx, int LDA);
int InvMatrixUpper(double *A, int p);

void getLikelihoods(double *theta, double *alpha, double *beta, double *guesing,
                    double *answers, int model, double *likelihoods);

double simpleLikelihoodPostDensity(double *theta, double *alpha, double *beta, double *guesing, double *answers,
                                   int numberDimensions, int numberAnswers,
                                   int normalPrior, double *mu, double *Sigma, int model);

void jointDistribution(double *thetaGrid, double *alpha, double *beta, double *guesing, double *answers,
                       int nrowGrid, int numberDimensions, int numberAnswers,
                       int normalPrior, double *mu, double *Sigma, int model, double *jointDistribution);

SEXP RjointDistribution(SEXP thetaGridR, SEXP alphaR, SEXP betaR, SEXP guessingR, SEXP answersR, SEXP nrowGridR,
                        SEXP numberDimensionsR, SEXP numberAnswersR, SEXP normalPriorR, SEXP muR, SEXP SigmaR, SEXP modelR);


/**
 * Symmetrize a matrix by copying the strict upper triangle into the
 * lower triangle.
 *
 * @param a pointer to a matrix in Fortran storage mode
 * @param nc number of columns (and rows and leading dimension) in the matrix
 *
 * @return a, symmetrized
 */
static R_INLINE double*
internal_symmetrize(double *a, int nc)
{
    int i,j;
    for (i = 1; i < nc; i++)
        for (j = 0; j < i; j++)
            a[i + j*nc] = a[j + i*nc];
    return a;
}
