/*

  Jasjeet Singh Sekhon 
  Harvard University
  http://jsekhon.fas.harvard.edu/
  jsekhon@fas.harvard.edu

  $Id: fit_lqd2.cpp,v 1.3 2004/02/14 06:28:13 wrm1 Exp $

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

long bcacmp(double *a, double *b);
void multi(double **in1, double **in2, double **out,
	   long row1, long col1, long row2, long col2, long outrowcol[2]);
void ResStd(double **SresRaw, double **Y, double *weights, double **Yprob,
	     long nobs, long ncats);
double lqd2(double **SresRaw, long nobs, long nvars_unique, long ncats);
double kth_smallest(double *a, long n, long k);

extern "C" 
{

  // <Rdefines.h> must appear inside the {extern "C"} declaration
  // because this header file, unlike <R.h>, does not have an {#ifdef
  // __cplusplus} statment included.
#include <Rdefines.h>

  SEXP original_fit_lqd2(SEXP I_nobs, SEXP I_ncats, SEXP I_nvars_total, SEXP I_nvars_unique,
		SEXP I_vec, SEXP I_Yvector, SEXP I_Xvector, SEXP I_weights)
  {
    SEXP ret;

    long nobs, ncats, nvars_total, nvars_unique, count;
    double **Y, **vec, *weights;

    long i, j, k, length;
    double **mu, foo, **Yprob, **SresRaw, fit_lqd2;

    double ***Xarray;

    nobs         = asInteger(I_nobs);
    ncats        = asInteger(I_ncats);
    nvars_total  = asInteger(I_nvars_total);
    nvars_unique = asInteger(I_nvars_unique);

    // printf("nobs: %d\n",nobs);

    Y        = (double **) calloc(nobs,sizeof(double));
    Yprob    = (double **) calloc(nobs,sizeof(double));
    mu       = (double **) calloc(nobs,sizeof(double));
    weights  = (double *)  calloc(nobs,sizeof(double));
    SresRaw  = (double **) calloc(nobs,sizeof(double));
    for (i=0; i<nobs; i++) 
      {
	Y[i]        = (double *) calloc(ncats,sizeof(double));
	Yprob[i]    = (double *) calloc(ncats,sizeof(double));
	mu[i]       = (double *) calloc(ncats,sizeof(double));
	SresRaw[i]  = (double *) calloc((ncats-1),sizeof(double));
      } // end of i

    Xarray  = (double ***) calloc(nobs,sizeof(double));
    for (i=0; i<nobs; i++) 
      {
	Xarray[i]  = (double **) calloc(nvars_total,sizeof(double));
	
	for (j=0; j<nvars_total; j++)
	  {
	    Xarray[i][j] = (double *) calloc(ncats,sizeof(double));
	  } // end of j
      } // end of i

    vec  = (double **) calloc(nvars_total,sizeof(double));
    for (k=0; k<nvars_total; k++)
      {
	vec[k] = (double *) calloc(ncats,sizeof(double));
      } //k

    // load up the parameter vector
    count = 0;
    for (j=0; j<ncats; j++)
      {
	for (k=0; k<nvars_total; k++)
	  {
	    vec[k][j] = REAL(I_vec)[count];
	    count++;
	    // printf("vec[%d][%d]: %lf\n", k, j, vec[k][j]);
	  }
      }

    // load up the Y matrix
    count = 0;
    for (k=0; k<ncats; k++)
      {
	for (i=0; i<nobs; i++)
	  {
	    Y[i][k] = REAL(I_Yvector)[count];
	    count++;
	    // printf("Y[%d][%d]: %lf\n", i, k, Y[i][k]);
	  } // end of i
      } // end of k

    // load up Xarray
    count = 0;
    for (k=0; k<ncats; k++)
      {
	for (j=0; j<nvars_total; j++)
	  {
	    for (i=0; i<nobs; i++)
	      {
		Xarray[i][j][k] = REAL(I_Xvector)[count];
		// printf("Xarray[%d][%d][%d]: %lf\n", i, j, k, Xarray[i][j][k]);
		count++;
	      } // end of i
	  } // J
      } // end of k


    // load up weights
    for (i=0; i<nobs; i++)
      {
	weights[i] = REAL(I_weights)[i];
      } // end of i

    for (i=0; i<nobs; i++)
      {
	for (j=0; j<ncats; j++)
	  {
	    foo = 0;
	    for (k=0; k<nvars_total; k++)
	      {
		foo = vec[k][j] * Xarray[i][k][j] + foo;
	      } // end of k
	    mu[i][j] = foo;
	  } // end of j;
	foo = 0;
	for (j=0; j<ncats; j++)
	  {
	    foo = exp(mu[i][j]) + foo;
	  }
	for (j=0; j<ncats; j++)
	  {
	    Yprob[i][j] = 1/foo * exp(mu[i][j]);

	    // printf("Yprob[%d][%d]: %lf, mu: %lf\n", 
	    // i, j, Yprob[i][j], mu[i][j]);
	  }
      } // end of i;

    ResStd(SresRaw, Y, weights, Yprob, nobs, ncats);
    
    fit_lqd2 = lqd2(SresRaw, nobs, nvars_unique, ncats);

    //free up ram
    free(weights);
    for (k=0; k<nvars_total; k++) 
      {
	free( (double *) vec[k]);
      } // end of k
    free(vec);

    for (i=0; i<nobs; i++) 
      {
	free( (double *) SresRaw[i]);
	free( (double *) mu[i]);
	free( (double *) Yprob[i]);
	free( (double *) Y[i]);

	for (j=0; j<nvars_total; j++) 
	  {
	    free( (double *) Xarray[i][j]);
	  } // end of j
	free( (double *) Xarray[i] );
      } // end of i    
    free( (double *) SresRaw);
    free( (double *) mu);
    free( (double *) Yprob);
    free( (double *) Y);
    free( (double *) Xarray);

    length=1;
    PROTECT(ret=allocVector(REALSXP,1));
    REAL(ret)[0]=fit_lqd2;
    UNPROTECT(1);
    return(ret);
  } // end of fit_lqd2

  SEXP kthSmallest(SEXP I_SortVector, SEXP I_obs, SEXP I_rank)
  {
    SEXP ret;

    long obs, rank, i;
    double *SortVector;
    double ReturnElement;
    
    obs  = asInteger(I_obs);
    rank = asInteger(I_rank)-1;
    
    SortVector = (double *) malloc(obs*sizeof(double));

    for (i=0; i<obs; i++)
      {
	SortVector[i] = REAL(I_SortVector)[i];
	//	printf("SortVector[%d]: %lf\n", i+1, SortVector[i]);
      }// end of i loop

    ReturnElement = kth_smallest(SortVector, obs, rank);

    free(SortVector);

    PROTECT(ret=allocVector(REALSXP,1));
    REAL(ret)[0]=ReturnElement;
    UNPROTECT(1);
    return(ret);
  } //end of kthSmallest

} // end of extern "C"


long bcacmp(double *a, double *b) 
{
  long i = 0;

  if (*a > *b) i = 1;
  else if (*a < *b) i = -1;
  return i;
}


void ResStd(double **SresRaw, double **Y, double *weights, double **Yprob,
	     long nobs, long ncats) 
{

  double **summat, **Q, **D, **R, **Rsum, **Or;
  long i, j, outrowcol[2];

  summat  = (double **) calloc(ncats, sizeof(double));
  Q       = (double **) calloc(nobs, sizeof(double));
  D       = (double **) calloc(nobs, sizeof(double));
  R       = (double **) calloc(nobs, sizeof(double));
  Rsum    = (double **) calloc(nobs, sizeof(double));
  Or      = (double **) calloc(nobs, sizeof(double));
  for (j=0; j<ncats; j++)
    {
      summat[j] = (double *) calloc(ncats, sizeof(double));
    } // end of j

  for (i=0; i<nobs; i++)
    {
      Q[i]    = (double *) calloc(ncats, sizeof(double));
      D[i]    = (double *) calloc((ncats-1), sizeof(double));
      R[i]    = (double *) calloc(ncats, sizeof(double));
      Rsum[i] = (double *) calloc(ncats, sizeof(double));
      Or[i]   = (double *) calloc((ncats-1), sizeof(double));
    } // end of i

  for (j=0; j<(ncats-1); j++)
    {
      for (i=(j+1); i<ncats; i++)
	{
	  summat[i][j] = 1.0;
	} // i
    } // j

  multi(Yprob, summat, Q, nobs, ncats, ncats, ncats, outrowcol);
  for (j=0; j<(ncats-1); j++)
    {
      if (j==0)
	{
	  for (i=0; i<nobs; i++)
	    {
	      D[i][j] = Yprob[i][j]*Q[i][j];
	    }  // end of if
	} // end of if
      if (j>0)
	{
	  for (i=0; i<nobs; i++)
	    {
	      D[i][j] = Yprob[i][j]*Q[i][j]/Q[i][(j-1)];	  
	    } // end of i
	} // end of if
    } // end of j
				      
  for (i=0; i<nobs; i++)
    {
      for (j=0; j<ncats; j++)
	{
	  R[i][j] = Y[i][j] - ( weights[i]*Yprob[i][j] );
	} // end of j
    } // end of i

  for (j=0; j<ncats; j++)
    {
      for (i=0; i<ncats; i++)
	{
	  summat[i][j] = 0.0;
	} // i
    } // j
  for (i=0; i<ncats; i++)
    {
      for (j=i; j<ncats; j++)
	{
	  summat[i][j] = 1.0;
	} // j
    } // i

  multi(R, summat, Rsum, nobs, ncats, ncats, ncats, outrowcol);

  for (j=0; j<(ncats-1); j++)
    {
      if (j==0)
	{
	  for (i=0; i<nobs; i++)
	    {
	      Or[i][j] = R[i][j];
	    } // end of i
	} //end of if
      if (j>0)
	{
	  for (i=0; i<nobs; i++)
	    {
	      Or[i][j] = R[i][j] + Rsum[i][(j-1)]*Yprob[i][j]/Q[i][(j-1)];
	    } // end of i	  
	} // end of j
    } // end of j

  for (i=0; i<nobs; i++)
    {
      for (j=0; j<(ncats-1); j++)
	{
	  SresRaw[i][j] = Or[i][j]/sqrt(weights[i]*D[i][j]);

	  //printf("SresRaw[%d][%d]: %lf, Or: %lf, weights: %lf, D: %lf\n", 
	  // i, j, SresRaw[i][j], Or[i][j], weights[i], D[i][j]);
	} // end of j
    } // end of i

  // free memory
  for (i=0; i<nobs; i++)
    {
      free( (double *) Or[i]);
      free( (double *) Rsum[i]);
      free( (double *) R[i]);
      free( (double *) D[i]);
      free( (double *) Q[i]);
    } // end of i
  for (j=0; j<ncats; j++)
    {
      free( (double *) summat[j]);
    } // end of j

  free(Or);
  free(Rsum);
  free(R);
  free(D);
  free(Q);
  free(summat);
  
} // end of res.std

void multi(double **in1, double **in2, double **out,
	   long row1, long col1, long row2, long col2, long outrowcol[2])
{
  long oi, oj, i;
  
  if (col1!=row2) {
    fprintf(stdout,"\nTHE MATRICES ARE NOT CONFORMABLE FOR MULIPLICATION\n");
    fprintf(stderr,"\nTHE MATRICES ARE NOT CONFORMABLE FOR MULIPLICATION\n");
    return;
  }
  
  outrowcol[0]=row1;
  outrowcol[1]=col2;
  
  for (oi=0;oi<outrowcol[0];oi++) {
    for (oj=0;oj<outrowcol[1];oj++) {
      out[oi][oj] = 0.0;
    }
  }
  
  for (oi=0;oi<outrowcol[0];oi++) {
    for (oj=0;oj<outrowcol[1];oj++) {
      for (i=0;i<col1;i++) {
	out[oi][oj] += in1[oi][i]*in2[i][oj];
	
      }
    }
  }
} // end of multi


double lqd2(double **SresRaw, long nobs, long nvars_unique, long ncats)
{
  double fit_lqd2, *rawres, *dif;
  long h, diflen, hidx, total_obs, i, ii, j, k, count;

  total_obs = nobs*(ncats-1);
  h  = (long) (total_obs+nvars_unique+1)/2;
  diflen = (total_obs*(total_obs-1))/2;
  hidx = (h*(h-1))/2;

  rawres  = (double *) calloc(total_obs, sizeof(double));
  dif  = (double *) calloc(diflen,sizeof(double));

  count = 0;
  for (i=0; i<nobs; i++)
    {
      for (j=0; j<(ncats-1); j++)
	{

	  rawres[count] = SresRaw[i][j];
	  //printf("rawres[%d] %lf, SresRaw[%d][%d]: %lf\n", count, rawres[count], i, j, SresRaw[i][j]);
	  count++;
	} // end of j
    } // end of i

  for (i=2; i<=total_obs; i++)
    {
      ii = ((i-1)*(i-2))/2;
      for(j=1; j<=(i-1); j++)
	{
	  k = ii+j;
	  dif[k-1] = fabs(rawres[i-1]-rawres[j-1]);
	}
    } // end of i

  /*
    qsort(dif, diflen, sizeof(double), (int (*)(const void *, const void *)) bcacmp);
    fit_lqd2 = dif[hidx-1] * 2.21914446599;
  */

  fit_lqd2 = kth_smallest(dif, diflen, hidx-1) * 2.21914446599;

    if (!R_finite(fit_lqd2))
      {
	printf("XXX\n");
	fit_lqd2 = 999991234;
      }

    //free memory
    free(dif);
    free(rawres);

    return(fit_lqd2);
} // end of lqd2


/*---------------------------------------------------------------------------
   Function :   kth_smallest()
   In       :   array of elements, # of elements in the array, rank k
   Out      :   one element
   Job      :   find the kth smallest element in the array
   Notice   :   use the median() macro defined below to get the median. 

                Reference:

                http://www.eso.org/~ndevilla/median/

                  Author: Wirth, Niklaus 
                   Title: Algorithms + data structures = programs 
               Publisher: Englewood Cliffs: Prentice-Hall, 1976 
    Physical description: 366 p. 
                  Series: Prentice-Hall Series in Automatic Computation 

 ---------------------------------------------------------------------------*/


double kth_smallest(double *a, long n, long k)
{
  long i,j,l,m ;
  double x, tmp;

    l=0 ; m=n-1 ;
    while (l<m) {
        x=a[k] ;
        i=l ;
        j=m ;
        do {
            while (a[i]<x) i++ ;
            while (x<a[j]) j-- ;
            if (i<=j) {
	      tmp=a[i];
	      a[i]=a[j];
	      a[j]=tmp;
	      i++ ; j-- ;
            }
        } while (i<=j) ;
        if (j<k) l=i ;
        if (k<i) m=j ;
    }
    return a[k] ;
} // end of kth_smallest
