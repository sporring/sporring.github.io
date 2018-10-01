#include <mex.h>
#include <math.h>
#include <stdlib.h>
#include "image.h"
#include "imfunc.h"

/* The computation routine. */
mxArray *dikuread(char* name)
{
  double image I;
  mxArray *A;
  double *ptr;
  int class, rows, cols, i, j;

  if (im_read_as_double(name,&I,&class) != IMERR_NONE )
    mexErrMsgTxt((char *)im_get_errmsg());
  else
    {
      rows = im_get_rows(I);
      cols = im_get_cols(I);
      
      if((A = mxCreateDoubleMatrix(rows, cols, mxREAL))==NULL)
	mexErrMsgTxt("Out of memory");
      else
	{
	  ptr = mxGetPr(A);
	  for(i = 0; i < rows; i++)
	    for(j = 0; j < cols; j++)
	      ptr[j*rows+i] = I[i][j];
	  /* 
	     Alternatively: 
	     memcpy(ptr,*I,rows*cols*sizeof(double));
	     which is twice as fast, but the result is transposed. 
	     */
	  im_free(I);
	}
    }

  return A;
}

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  char *buf;
  unsigned int buflen;
  int status;

  /* I = dikuread2(name);                                       */
  /*    where I is the resulting double matrix                  */
  /*          name is a text string of a DIKU image name        */

  /* Check for proper number of arguments. */
  if(nrhs != 1)
    mexErrMsgTxt("Only one input argument allowed.");
  else if (nlhs > 1)
    mexErrMsgTxt("Only one output argument allowed.");

  /* Find out how long the input string array is. */
  buflen = (mxGetM(prhs[0])*mxGetN(prhs[0]))+1;

  /* Allocate enough memory to hold the converted string. */
  if((buf = mxCalloc(buflen, sizeof(char))) == NULL)
    mexErrMsgTxt("Out of memory");
  else
    {
      /* Copy the string data from prhs[0] and place it into buf. */
      status = mxGetString(prhs[0], buf, buflen);
      if (status != 0)
	mexErrMsgTxt("Could not convert string data.");

      /* Call the C subroutine.  Note: buf and rev_strin ar C strings. */
      plhs[0] = dikuread(buf);
    }
}
