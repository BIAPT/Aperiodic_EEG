#include "mex.h"
#include <math.h>

/*******************************************************************/
/*                                                                 */
/*          Programme written by Tofik Amara, July 2002            */
/*              Modified by Steeve Zozor, January 2003             */
/*                                                                 */
/*      Based on the algorithm proposed by Kaspar and Schuster     */
/*    Physical Review A, vol. 36, no. 2, pp. 842-848, July 1987    */
/*                                                                 */
/* vectorial signal X -> scalar signal Z by z_j = sum x_{i,j} al^j */
/*                                                                 */
/*          Bug reports: steeve.zozor@lis.inpg.fr                  */
/*                                                                 */
/*******************************************************************/


/* function to interfacage Matlab and C . */

 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

{
  void complexity(double*, int ,int, double, double*); /* prototype de la fonction appele */

  double *signal;        /* matrix of the signals */
  double al;             /* size of the alphabet */
  int nt,ii;                /* duration of the signal */
  int nc;                /* dimension of the signal */
  double *c;             /* joint Lempel-Ziv complexity */


 /* check for proper number of arguments.*/
  if(nrhs!=2){mexErrMsgTxt(" not proper input arguments ");}
  if(nlhs!=1){mexErrMsgTxt(" not proper output");}


 /* Create matrix for the return arguments. */
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);


  /* Assign pointers to each input and output */ 
  signal = mxGetPr(prhs[0]);
  al = mxGetScalar(prhs[1]);

  c = (double *) mxGetPr(plhs[0]);


  /* get the length of the matrix */
  nt= mxGetN(prhs[0]);      /* column = time */
  nc= mxGetM(prhs[0]);      /* line = components */


  /* call of the routine to evaluate the complexity */
  complexity(signal,nc,nt,al,c);
}


/********************************************************************/

 /* Routine to permut lines and rows, since C put the signal in line order */

void permut(double *signal, int nc, int nt, double *sig)

     /* Inputs: signal = address of the signal matrix         */
     /*         nc     = number of components                 */
     /*         nt     = duration of the signal               */
     /* Output: sig = address; the matrix was just transposed */

{
  int indt, indc;
  for(indt=0;indt<nt;indt++)
    {
      for(indc=0;indc<nc;indc++)
	{
	  sig[indt+indc*nt]=signal[indc+indt*nc];
	}
    }
}

/********************************************************************/

 /* Routine to convert the multidimensional signal to a scalar signal */

void multi_uni(double *signal, int nc, int nt, double al, double *sig)

     /* Inputs: signal = address                                      */
     /*         nc     = number of components                         */
     /*         nt     = duration of the signal                       */
     /*         al     = size of the alphabet                         */
     /* Output: sig =  scalar signal (al-ary decomposition = signal)  */
     /* RECALL, C put the element of a matrix according to the lines! */

{
  double a;
  int indt,indc;

  for(indt=0;indt<nt;indt++)
    {
      /* initialisation of sig */
      sig[indt]=0.0;
      a=1.0;
      /* al-ary decomposition -> scalar */
      for(indc=0;indc<nc;indc++)
	{
	  sig[indt]=sig[indt]+a*signal[indc+indt*nc];
	  a=a*al;
	}
    }
}

/********************************************************************/

 /* Function to evaluate the joint Lempel-Ziv complexity */

void complexity(double *signal, int nc, int nt, double al,double *c)

     /* Inputs: signal = address of the signal      */
     /*         nc     = total number of components */
     /*         nt     = duration of the signal     */
     /*         al     = size of the alphabet       */
     /* Output: c = joint complexity of the signal  */


{
  int i,k,kmax;
  int l=1;
  int m=1;
  double *sig;

  /* convertion vector -> scalar */
  sig=(double*)malloc(nt*sizeof(double));
  multi_uni(signal,nc,nt,al,sig);

  /* initialization of c */
  *c=1.0;

  step1 : i=0; k=1 ; kmax=1;
  step2 : if ( *(sig+i+k-1)==*(sig+l+k-1))
            {
	      k = k+1;
              if((l+k)>nt)
                {
		  *c=(*c)+1.0;
                  goto output;
                }
	      else
                {
		  goto step2;
                }
	    }

           else
	     {
	       if ( k > kmax){kmax = k;}  
	       i=i+1;
	       if (i==l)    
		 {
		   *c=(*c)+1.0;
		   m=m+1;
		   l=l+kmax;
		   if(l+1>nt)
		     {
		       goto output;
		     }
		   else
		     {
		       goto step1;
		     }
		 }     
	       else
		 {
		   k=1;
		   goto step2;
		 }
	     }
  output : free(sig);
}

