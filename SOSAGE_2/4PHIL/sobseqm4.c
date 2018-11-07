/* 
  this file must be compiled uing mex sobseqm.c
 then from matlab do the call sobseqm(2)
*/

#include "mex.h"

#define NRANSI
#define MAXBIT 30
#define MAXDIM 6
#include<math.h>
#define m_ec2 0.5109e6
#define cms 299792458
#include <stdlib.h>
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))
/*
 use a Sobol sequence from Numerical Recipe.
*/
void sobseq(double nf[], double x[])
{
        
	int j,k,l,n;
	unsigned long i,im,ipp;
	static double fac;
	static unsigned long in,ix[MAXDIM+1],*iu[MAXBIT+1];
	static unsigned long mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
	static unsigned long ip[MAXDIM+1]={0,0,1,1,2,1,4};
	static unsigned long iv[MAXDIM*MAXBIT+1]={
		0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};
	
	n=(int)nf[0];	
	if (n < 0) {
		for (k=1;k<=MAXDIM;k++) ix[k]=0;
		in=0;
		if (iv[1] != 1) return;
		fac=1.0/(1L << MAXBIT);
		for (j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM) iu[j] = &iv[k];
		for (k=1;k<=MAXDIM;k++) {
			for (j=1;j<=mdeg[k];j++) iu[j][k] <<= (MAXBIT-j);
			for (j=mdeg[k]+1;j<=MAXBIT;j++) {
				ipp=ip[k];
				i=iu[j-mdeg[k]][k];
				i ^= (i >> mdeg[k]);
				for (l=mdeg[k]-1;l>=1;l--) {
					if (ipp & 1) i ^= iu[j-l][k];
					ipp >>= 1;
				}
				iu[j][k]=i;
			}
		}
	} else {
		im=in++;
		for (j=1;j<=MAXBIT;j++) {
			if (!(im & 1)) break;
			im >>= 1;
		}
		if (j > MAXBIT) printf("MAXBIT too small in sobseq");
		im=(j-1)*MAXDIM;
		for (k=1;k<=IMIN(n,MAXDIM);k++) {
			ix[k] ^= iv[im+k];
			x[k-1]=ix[k]*fac;
		}
	}
}

void sobseqm4 (double nf[], double N[], double Img[], double x[])
{
	sobseq(nf,x); 
}


void mexFunction( int nlhs, mxArray *plhs[],

                  int nrhs, const mxArray *prhs[] )
{
  double *nf, *N, X, Y, rnd;
  int i,j, iX, iY, inc, mrows, ncols;
  double *x, *XX, *Img;

 

  /* Check for proper number of arguments. */

  if(nrhs!=3) {

    mexErrMsgTxt("Two  input required.");

  } else if(nlhs>3) { 

    mexErrMsgTxt("Too many output arguments");

  }

  

  /* Assign pointers to each input and output. */

  nf = mxGetPr(prhs[0]);
  N  = mxGetPr(prhs[1]); 
  mrows = mxGetM(prhs[2]) ;
  ncols = mxGetN(prhs[2]) ;
  Img   = mxGetPr(prhs[2]);  
  printf ("%d\t%d", mrows, ncols);

/*
  for(i=0; i < ncols; i++) {
    for(j=0; j < mrows; j++) {
      AR[i][j] = Img[i*mrows + j] ;
    }
  }
*/
  /* Create matrix for the return argument. */

  plhs[0] = mxCreateDoubleMatrix(MAXDIM, (int)N[0], mxREAL);
  plhs[1] = mxCreateDoubleMatrix(MAXDIM, 1 , mxREAL);


  XX = mxGetPr(plhs[0]);
  x =  mxGetPr(plhs[1]);

  /* Call the xtimesy subroutine. */
  inc=0;

  for (i=0; i<(int)N[0] ; i++) {
     do{
        sobseqm4(nf,N,Img,x);
	X=1+x[0]*ncols;
	Y=1+x[1]*mrows;
	rnd=x[2];
	iX=(int)X;
	iY=(int)Y;
     } while ((rnd>Img[iX*mrows + iY]));
     XX[inc]=X;
     inc++;
     XX[inc]=Y;
     inc++;
	
     for (j=2; j<MAXDIM; j++) {
        XX[inc]=x[j];
	inc=inc+1;
     }
  } 
  
}

 


#undef MAXBIT
#undef MAXDIM
#undef NRANSI

