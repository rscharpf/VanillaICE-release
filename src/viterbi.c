#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <Rinternals.h>
#include <R.h>

static void updateTransitionMatrix(double *pAA, const int t, const int nCols, const int NS, double *tau)
{
  int i, j;
      for (i=0; i<nCols; ++i)
	{
	  for (j=0; j<nCols; ++j)
	    {
	      int offset;
	      offset = j * nCols + i;
	      if(i == j)
		{
		  *(pAA + offset) = tau[t-1];
		}
	      else
		{
		  *(pAA + offset) = (1-tau[t-1])/(nCols-1);
		}
	    }
	}
}

static void getIndexAndMaxVal(const double *pVec, const int len, double *pMaxVal, int *pMaxIdx)
{
  int i;
  *pMaxVal = *pVec;
  *pMaxIdx = 0;
  for (i=1; i<len; ++i)
  {
    if (pVec[i] > *pMaxVal)
    {
      *pMaxIdx = i;
      *pMaxVal = pVec[i];
    }
  }
}

static void getMatrixIndexAndMaxVal(const double *pMat, const int nCols, double *pMaxVal, int *pMaxIdx, const int nRows)
{
  int i;
  *pMaxVal = *pMat;
  *pMaxIdx = 0;
  for (i=1; i<nCols; ++i)
  {
    if (pMat[i*nRows] > *pMaxVal)
    {
      *pMaxIdx = i;
      *pMaxVal = pMat[i*nRows];
    }
  }
}

/**
 * viterbi
 * \param pBeta - log emission probabilities SxT
 * \param initialP - log initial state probabilities - length S
 * \param tau - transition probabilities - original scale
 * \param pArm - indicator for chromosome arm
 * \param S - number of columns of beta (number of states)
 * \param T - number of rows of beta matrix
 * \param pQHat - output vector of integers
 * \param pDelta - output vector of doubles ...this is the forward variable
 * \param c1  Pr( Normal -> altered)
 * \param c2  Pr( altered -> normal)
 * \param c3  Pr( altered -> altered)
 * \param normalState  index of normal state
 * \param pAA
 */
void viterbi(double *pBeta, double *initialP, double *tau,
             int  *pArm, int *S, int *T, int *pQHat, double *pDelta,
             double *c1, double*c2, double *c3, int *normalState, double *pAA)
{
  /**  RS double *pDelta, *pAA, *pDeltaTempSum, Pstar; */
  /** double *pAA, *pDeltaTempSum, Pstar, *tp; */
  double *pDeltaTempSum, Pstar;
  int i,j,t;
  int nRows, nCols, *pPsi;
  int NS;


  NS = *normalState - 1;
  nRows = *T;
  nCols = *S;

  /** RS
      pDelta = (double *)R_alloc(sizeof(double), nRows * nCols); */
  pPsi = (int *)R_alloc(sizeof(int), nRows * nCols);
  /** transition probability matrix */
  /** pAA = (double *)R_alloc(sizeof(double), nCols * nCols); */
  pDeltaTempSum = (double *)R_alloc(sizeof(double), nCols);

  /** what is this notation? *(pDelta + nrows*j) C uses vectors and
      not matrices.  so the Nth row in pDelta is set to initialP + the
      emission probability for the Nth row*/
  for (j=0; j<nCols; ++j)
  {
    *(pDelta + nRows*j) = initialP[j] + *(pBeta + nRows*j);
    *(pPsi + nRows*j) = 0;
  }

  for (t=1; t<nRows; ++t)
    {
      /*      if (strcmp(*(ppArm + t), *(ppArm + t - 1)) != 0)*/
      if(*(pArm + t) != *(pArm + t - 1))
	{
	  for (j=0; j<nCols; ++j)
	    {
	      *(pDelta + j*nRows + t) = initialP[j] + *(pBeta + j*nRows + t);
	      *(pPsi + j*nRows + t) = 0;
	    }
	  continue;
	}
      for (i=0; i<nCols; ++i)
	{
	  for (j=0; j<nCols; ++j)
	    {
	      int offset;
	      offset = j * nCols + i;
	      /* if (i == j)*/
	      if(i == NS)
		{
		  if(i == j)  /* probability of staying in the normal state */
		    {
		      *(pAA + offset) = 1 - ((1-tau[t-1]) * (nCols - 1) * *c1);
		      /* printf(i, j); */
		      /* probability of staying in the same state */
		      /* *(pAA + offset) = tau[t-1]; */
		    }
		  else /* probability of leaving normal state */
		    {
		      *(pAA + offset) = *c1 * (1-tau[t-1]);
		    }
		}
	      else   /* transitioning from an altered state */
		{
		  if(i == j)  /* staying in the same altered state */
		    {
		      /* c2 = scalar for transitioning from normal to altered state */
		      *(pAA + offset) = 1 - (1 - tau[t-1]) * (*c2 + (nCols - 2) * *c3);
		      /* *(pAA + offset) = (1-tau[t-1])/(nCols-1); */
		    }
		  else /* leaving altered state */
		    {
		      if(j == NS) /* going back to normal state */
			{
			  *(pAA + offset) = *c2 * (1 - tau[t-1]);
			}
		      else  /* going to another altered state */
			{
			  *(pAA + offset) = *c3 * (1 - tau[t-1]);
			}
		    }
		}
	      /* *(pAA + offset) = log ( *(pAA + offset) * *(tau_scale + offset) );*/
	      *(pAA + offset) = log ( *(pAA + offset) );
	    }
	}
      for (j=0; j<nCols; ++j)
	{
	  double maxDeltaTempSum;
	  int maxDeltaSumIdx = 0;
	  /* Sum the jth column of AA and the (t-1)th row of delta.
             The jth column of AA occupies
	     consecutive memory locations starting at the memory location pAA + nrow(AA)*j.
	     Since AA is a square matrix, that starting address can be
	     expressed as pAA + nCols*j */
	  for (i=0; i<nCols; ++i)
	    {
	      pDeltaTempSum[i] = pAA[j * nCols + i] + pDelta[(t-1) + i * nRows];
	    }
	  /* Needs update */
	  getIndexAndMaxVal( (double *)(pDeltaTempSum), nCols, &maxDeltaTempSum, &maxDeltaSumIdx);
	  *(pPsi + j * nRows  + t) = maxDeltaSumIdx;
	  *(pDelta + j*nRows + t) = maxDeltaTempSum + *(pBeta + j*nRows + t);
	}
    }

  /* Needs update */
  getMatrixIndexAndMaxVal( (double *)(pDelta + nRows-1), nCols, &Pstar, (int *)(pQHat + nRows-1), nRows);
  for (t=nRows-2; t>= 0; --t)
    {
      /*if (strcmp(*(ppArm + t), *(ppArm + t + 1)) != 0)*/
      if(*(pArm + t) != *(pArm + t + 1))
	{
	  double maxVal;

	  /* Needs update */
	  getMatrixIndexAndMaxVal((double *)(pDelta + t), nCols, &maxVal, &pQHat[t], nRows);
	}
      else
	{
	  pQHat[t] = *(pPsi + pQHat[t+1] * nRows + (t+1));
	}
    }

  /* Array indices in R are 1-based so add one to these values. */
  for (i=0; i<nRows; ++i) {
    pQHat[i] += 1;
    if (i > 0)
      for (j=0; j<nCols; ++j)
        pPsi[j*nRows + i] += 1;}
}

void viterbi2(double *pBeta, /* emission prob */
	      double *initialP,  /* initial state prob */
	      double *tau,  /* scalar for transition prob */
	      int  *pArm,  /* chromosome arm */
	      int *S,  /* number states */
	      int *T,  /* number markers */
	      int *pQHat, /* state path */
	      double *pAlpha, /* forward variable */
	      double *pBack, /* backward variable */
	      int *normalState,
	      double *scalingFactor)
{
  /* *************************************************************************** */
  /*  : ASSUME pBeta is on probability scale (not log scale) */
  /*  :  changes to pDelta throughout */
  /*  : need probability scale for reestimation formulas */
  /* *************************************************************************** */
  /**  RS double *pDelta, *pAA, *pDeltaTempSum, Pstar; */
  /** double *pAA, *pDeltaTempSum, Pstar, *tp; */
  double *pDeltaTempSum, Pstar;
  double *pBackTempSum;
  int i,j,t;
  int nRows, nCols, *pPsi;
  int NS;
  double *pAA2, *pAA;
  double *scalingFactorB;
  double *pDelta;
  NS = *normalState - 1;
  nRows = *T;
  nCols = *S;

  pDelta=(double *)R_alloc(sizeof(double), nRows*nCols);
  scalingFactorB = (double *)R_alloc(sizeof(double), nRows);
  pPsi = (int *)R_alloc(sizeof(int), nRows * nCols);
  pDeltaTempSum = (double *)R_alloc(sizeof(double), nCols);
  pBackTempSum = (double *)R_alloc(sizeof(double), nCols);
  pAA = (double *)R_alloc(sizeof(double), nCols * nCols);
  pAA2 = (double *)R_alloc(sizeof(double), nCols * nCols);

  double scalingFactorSum, scalingFactorSumB;
  scalingFactorSum=0; scalingFactorSumB=0;
  for (j=0; j<nCols; ++j)
  {
    *(pDelta + nRows*j) = initialP[j] * *(pBeta + nRows*j);
    *(pPsi + nRows*j) = 0;
    scalingFactorSum=scalingFactorSum+*(pDelta+nRows*j);
  }
  scalingFactor[0]=1/scalingFactorSum;
  for(j=0; j<nCols; ++j)
    {
      *(pDelta + nRows*j)=scalingFactor[0]* pDelta[nRows*j];
      *(pAlpha + nRows*j) = *(pDelta + nRows*j);
      /* iterating over element j, this is time T (index T-1) */
      *(pBack + nRows*j + nRows-1) = 1;
    }
  /*RS: For the backwards variable, we need a counter (k) that goes in
    the opposite direction of t.  I think it could be defined as a
    function of t.  k should go from time T-1 to 1 (or index T-2 to 0)
    t: 1, ... T-1
    k: T-2, ..., 0
  */
  int k;
  for (t=1; t<nRows; ++t)
    {
      k = nRows-t-1;  /* when t=T-1, k= T-(T-1)-1=0 */
      updateTransitionMatrix(pAA, t, nCols, NS, tau);
      scalingFactorSum=0.0;
      for (j=0; j<nCols; ++j)
	{
	  double maxDeltaTempSum;
	  int maxDeltaSumIdx = 0;
	  double alphaSum = 0;
	  for (i=0; i<nCols; ++i)
	    { /* eq 92a */
	      /* want to integrate out column j: sum_j aij */
	      /* we later assign this prob to element j+1 of the t+1 row of the forward variable */
	      pDeltaTempSum[i] = pAA[j * nCols + i] * pAlpha[(t-1) + i * nRows];
	      alphaSum=alphaSum+pDeltaTempSum[i];
	    }
	  /* Needs update */
	  getIndexAndMaxVal( (double *)(pDeltaTempSum), nCols, &maxDeltaTempSum, &maxDeltaSumIdx);
	  *(pPsi + j * nRows  + t) = maxDeltaSumIdx;
	  /* eq 33a */
	  *(pDelta + j*nRows + t) = maxDeltaTempSum * *(pBeta + j*nRows + t);
	  /* eq 20 */
	  *(pAlpha + j*nRows + t) = alphaSum * *(pBeta + j*nRows + t);
	  /* rescale pDelta */
	  /* (Rabiner eq 91)  */
	  scalingFactorSum=scalingFactorSum + *(pAlpha + j*nRows+t);
	}
      scalingFactor[t] = 1/scalingFactorSum;
      for(j=0; j<nCols; ++j)
	{
	  *(pDelta + j*nRows +t) = scalingFactor[t] * pDelta[j*nRows + t];
	  *(pAlpha + j*nRows +t) = scalingFactor[t] * pAlpha[j*nRows + t];
	}
      /* backwards variable */
      /* if k is zero, this is undefined */
      updateTransitionMatrix(pAA2, k+1, nCols, NS, tau);
      scalingFactorSumB=0.0;
      for(i=0;i<nCols;++i)
	{
	  double backSum = 0;
	  for(j=0;j<nCols;++j)
	    {
	      /* want to integrate out row i: sum_j aij */
	      pBackTempSum[j] = pAA2[j * nCols + i] * pBack[j*nRows + k+1] * *(pBeta + i*nRows + k);
	      backSum=backSum+pBackTempSum[j];
	    }
	  /* assign to element i of row k */
	  *(pBack + i*nRows + k) = backSum;
	  scalingFactorSumB=scalingFactorSumB + *(pBack + i*nRows + k);
	}
      scalingFactorB[k] = 1/scalingFactorSumB;
      for(j=0;j<nCols;++j)
	{
	  *(pBack + j*nRows +k) = scalingFactorB[k] * pBack[j*nRows +k];
	}
    }
  getMatrixIndexAndMaxVal( (double *)(pDelta + nRows-1), nCols, &Pstar, (int *)(pQHat + nRows-1), nRows);
  for (t=nRows-2; t>= 0; --t)
    {
      if(*(pArm + t) != *(pArm + t + 1))
	{
	  double maxVal;
	  getMatrixIndexAndMaxVal((double *)(pDelta + t), nCols, &maxVal, &pQHat[t], nRows);
	}
      else
	{
	  pQHat[t] = *(pPsi + pQHat[t+1] * nRows + (t+1));
	}
    }
  /* Array indices in R are 1-based so add one to these values. */
  for (i=0; i<nRows; ++i) {
    pQHat[i] += 1;
    if (i > 0)
      for (j=0; j<nCols; ++j)
        pPsi[j*nRows + i] += 1;
  }
}
