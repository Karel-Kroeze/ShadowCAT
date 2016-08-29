#include "joint_distr.h"


// Invert N by N matrix A
//
// A: Matrix to invert
// N: Number of rows/columns (order) of the matrix A
void inverse(double *A, int N) {
  int IPIV[N+1], INFO, LWORK = N*N;
  double WORK[LWORK];
  
  F77_NAME(dgetrf)(&N, &N, A, &N, IPIV, &INFO);
  F77_NAME(dgetri)(&N, A, &N, IPIV, WORK, &LWORK, &INFO);
}

// Compute determinant of an N by N matrix A
//
// A: Matrix to invert
// N: Number of rows/columns (order) of the matrix A
// LDA: Leading dimension of the matrix A
// doLog: If 1, return log determinant, else, return determinant
double matrixDet(double *A, int N, int LDA, int doLog) {
  //SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
  int i=0, info=0;
  double B[N*N], logDet=0;
  
  //Memcpy(B,A,N*N);
  for (i=0;i<N;i++) {
    Memcpy(&B[i*N],&A[i*LDA],N);
  }
  
  F77_CALL(dpotrf)("U", &N, B, &N, &info);
  if (info)
    Rprintf("Cholesky decomposition in matrixDet() returned nonzero info %d.\n",info);
  
  for (i=0;i<N;i++) {
    logDet += 2 * log(B[i*N+i]);
  }
  
  if (doLog)
    return(logDet);
  else
    return(exp(logDet));
}

// compute quadratic form x'Ax
//
// x: The vector x
// A: The matrix A
// N: Number of rows/columns (order) of A, length of X
// incx: Memory increments between subsequent x values
// LDA: Leading dimension of A
double quadform(double *x, double *A, int N, int incx, int LDA) {
  int Nsqr = N*N,info,i=0,j=0;
  double *B = Calloc(Nsqr,double);
  //double one=1;
  //double zero=0;
  double sumSq=0;
  double y[N];
  int iOne=1;
  
  for (i=0;i<N;i++) {
    y[i] = x[i*incx];
  }
  for (i=0;i<N;i++) {
    Memcpy(&B[i*N],&A[i*LDA],N);
  }
  
  F77_NAME(dpotrf)("U", &N, B, &N, &info);
  F77_NAME(dtrmv)("U","N","N", &N, B, &N, y, &iOne);
  
  for (i=0;i<N;i++) {
    sumSq += y[i]*y[i];
  }
  
  Free(B);
  
  return(sumSq);
}


// Get right upper part of inverted matrix A
//
// A: Matrix to invert
// p: Number of rows/columns (order) of A
int InvMatrixUpper(double *A, int p) {
    int info1, info2;
    F77_NAME(dpotrf)("U", &p, A, &p, &info1);
    F77_NAME(dpotri)("U", &p, A, &p, &info2);      
    //make sure you make it symmetric later...
    return(info1);
}


// Deze functie moet, afhankelijk van model, de c++ functie PROB_3PLM(), PROB_GRM(), PROB_SM() of PROB_GPCM()
// aanroepen, en daar het list element l uithalen en "plaatsen" in *likelihoods.
//
// theta: vector with theta for which the likelihoods per item should be computed.
// alpha: alpha matrix.
// beta: beta matrix.
// guessing: vector with guessing parameters.
// answers: vector with answers to items.
// model: model that should be used for computation of the likelihoods, 0 = "3PLM", 1 = "GPCM", 2 = "SM", 3 = "GRM".
// likelihoods: pointer to memory space where the likelihoods are stored
void getLikelihoods(double *theta, double *alpha, double *beta, double *guessing,
                    double *answers, int model, double *likelihoods) {
  // TO DO
  if (model == 0)
    PROB_3PLM()$l
  if (model == 1)
    etc.
}

// Get likelihood or posterior density of theta
//
// theta: vector with theta for which the likelihood or posterior density should be computed for the complete test.
// alpha: alpha matrix.
// beta: beta matrix.
// guessing: vector with guessing parameters.
// answers: vector with answers to items.
// numberDimensions: number of dimensions (length of theta, number columns of alpha).
// numberAnswers: number of answers given (lengtn of answers).
// normalPrior: if 0, likelihood is returned, if 1, posterior density is returned.
// mu: mean of normal prior.
// Sigma: covariance matrix of normal prior.
// model: model that should be used for computation of the likelihoods, 0 = "3PLM", 1 = "GPCM", 2 = "SM", 3 = "GRM".
double simpleLikelihoodPostDensity(double *theta, double *alpha, double *beta, double *guessing, double *answers,
                                   int numberDimensions, int numberAnswers,
                                   int normalPrior, double *mu, double *Sigma, int model) {
  double sumLogLikelihoods = 0;
  double likelihoods[numberAnswers];
  getLikelihoods(theta, alpha, beta, guessing, answers, model, likelihoods);
  
  for (int i = 0; i < numberAnswers; i++) {
    if (likelihoods[i] <= 0)
      likelihoods[i] = .0000000001;
    sumLogLikelihoods += log(likelihoods[i]);
  }
  
  if (normalPrior) {
    const double pi = 3.14159265358979323846;
    double logDeterminantSigma, logNormalDensity, logPostDensity;
    
    logDeterminantSigma = matrixDet(Sigma, numberDimensions, numberDimensions, 1);
    InvMatrixUpper(Sigma, numberDimensions);
    internal_symmetrize(Sigma, numberDimensions);
    logNormalDensity = -(numberDimensions/2 * log(2*pi) + .5 * logDeterminantSigma) - quadform(mu, Sigma, numberDimensions, 1, numberDimensions) / 2;
    logPostDensity = sumLogLikelihoods + logNormalDensity;
    return(exp(logPostDensity));  
  }
  
  else {
    return(exp(sumLogLikelihoods));
  }
}


// Get simpleLikelihoodPostDensity() for each (row of) theta in the grid
//
// thetaGrid: thetaGrid matrix as vector by row
// alpha: alpha matrix.
// beta: beta matrix.
// guessing: vector with guessing parameters.
// answers: vector with answers to items.
// nrowGrid: number of rows in the grid
// numberDimensions: number of dimensions (length of theta, number columns of alpha).
// numberAnswers: number of answers given (lengtn of answers).
// normalPrior: if 0, likelihood is returned, if 1, posterior density is returned.
// mu: mean of normal prior.
// Sigma: covariance matrix of normal prior.
// model: model that should be used for computation of the likelihoods, 0 = "3PLM", 1 = "GPCM", 2 = "SM", 3 = "GRM".
// jointDistribution:  pointer to memory space where the joint densities are stored.
void getJointDistribution(double *thetaGrid, double *alpha, double *beta, double *guessing, double *answers,
                       int nrowGrid, int numberDimensions, int numberAnswers,
                       int normalPrior, double *mu, double *Sigma, int model, double *jointDistribution) {
  for (int i = 0; i < nrowGrid; i++) {
    jointDistribution[i] = simpleLikelihoodPostDensity(&thetaGrid[i * numberDimensions], alpha, beta, guessing,
                                                       answers, numberDimensions, numberAnswers, normalPrior, mu, Sigma, model);
  }
}


// Make getJointDistribution() available from R.
//
// thetaGridR: thetaGrid matrix as vector by row
// alphaR: alpha matrix.
// betaR: beta matrix.
// guessingR: vector with guessing parameters.
// answersR: vector with answers to items.
// nrowGridR: number of rows in the grid
// numberDimensionsR: number of dimensions (length of theta, number columns of alpha).
// numberAnswersr: number of answers given (lengtn of answers).
// normalPriorR: if 0, likelihood is returned, if 1, posterior density is returned.
// muR: mean of normal prior.
// SigmaR: covariance matrix of normal prior.
// modelR: model that should be used for computation of the likelihoods, 0 = "3PLM", 1 = "GPCM", 2 = "SM", 3 = "GRM".
SEXP RjointDistribution(SEXP thetaGridR, SEXP alphaR, SEXP betaR, SEXP guessingR, SEXP answersR, SEXP nrowGridR,
                        SEXP numberDimensionsR, SEXP numberAnswersR, SEXP normalPriorR, SEXP muR, SEXP SigmaR, SEXP modelR) {
  double *thetaGrid = REAL(thetaGridR);
  double *alpha = REAL(alphaR);
  double *beta = REAL(betaR);
  double *guessing = REAL(guessingR);
  double *answers = REAL(answersR);
  int nrowGrid = INTEGER_VALUE(nrowGridR);
  int numberDimensions = INTEGER_VALUE(numberDimensionsR);
  int numberAnswers = INTEGER_VALUE(numberAnswersR);
  int normalPrior = INTEGER_VALUE(normalPriorR);
  double *mu = REAL(muR);
  double *Sigma = REAL(SigmaR);
  int model = INTEGER_VALUE(modelR);
  
  SEXP jointDistributionR;
  PROTECT(jointDistributionR = allocVector(REALSXP, nrowGrid));
  
  getJointDistribution(thetaGrid, alpha, beta, guessing, answers,
                       nrowGrid, numberDimensions, numberAnswers,
                       normalPrior, mu, Sigma, model, REAL(jointDistributionR));
  
  UNPROTECT(1);
  
  return(jointDistributionR);   
}

