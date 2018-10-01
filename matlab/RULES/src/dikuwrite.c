#include <mex.h>
#include <math.h>
#include <stdlib.h>
#include "image.h"
#include "imfunc.h"

/* The computation routine. */
dikuwrite(const mxArray *A, char *name, int type)
{
  double image I;
  double *ptr;
  int class, rows, cols, i, j;

  rows = mxGetM(A);
  cols = mxGetN(A);

  if(im_dynamic(I, rows, cols, 0) != IMERR_NONE)
    mexErrMsgTxt("Out of memory");
  else
    {
      ptr = mxGetPr(A);
      for(i = 0; i < rows; i++)
	for(j = 0; j < cols; j++)
	  I[i][j] = ptr[j*rows+i];
      if (im_write_double_as(name,type,I) != IMERR_NONE )
	mexErrMsgTxt((char *)im_get_errmsg());
      im_free(I);
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  char *buf,*typestring;
  unsigned int buflen;
  int status,type;

  /* dikuwrite2(I,name,type);                                   */
  /*    where I is the double matrix to be written              */
  /*          name is a text string of a DIKU image name        */
  /*          type is a text string of a DIKU image type        */
  /*               one of (grey,short,int,float,double).        */
  /*               default is double.                           */

  /* Check for proper number of arguments. */
  if((nrhs < 2) || (nrhs > 3))
    mexErrMsgTxt("Only two or three input argument allowed.");
  else if (nlhs > 0)
    mexErrMsgTxt("No output argument allowed.");

  buflen = (mxGetM(prhs[1])*mxGetN(prhs[1]))+1;
  buf = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[1], buf, buflen);
  if (status != 0)
    mexErrMsgTxt("Could not convert string data.");

  buflen = (mxGetM(prhs[2])*mxGetN(prhs[2]))+1;
  typestring = mxCalloc(buflen, sizeof(char));
  status = mxGetString(prhs[2], typestring, buflen);
  if (status != 0)
    mexErrMsgTxt("Could not convert string data.");

  if(nrhs > 2)
    {
      if((strcasecmp(typestring,"grey")==0) || (strcasecmp(typestring,"gray")==0))
	type = IM_CLASS_GREY;
      else if(strcasecmp(typestring,"short")==0)
	type = IM_CLASS_SHORT;
      else if(strcasecmp(typestring,"int")==0)
	type = IM_CLASS_INT;
      else if(strcasecmp(typestring,"float")==0)
	type = IM_CLASS_FLOAT;
      else if(strcasecmp(typestring,"double")==0)
	type = IM_CLASS_DOUBLE;
      else
	mexErrMsgTxt("Non-recognized DIKU image type.  Use one of ('grey','gray','short','int','float','double').");
    }
  else
    type = IM_CLASS_DOUBLE;

  /* Call the C subroutine.  Note: buf and rev_strin ar C strings. */
  dikuwrite(prhs[0],buf,type);
}
