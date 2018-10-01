#include <stdio.h>
#include "mex.h"

#define BUFLEN 256
#define	MAX(A, B)	((A) > (B) ? (A) : (B))

/* I = rawread('images/t.hips',[681,614,3,1]); */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
  char  filename[BUFLEN], s[BUFLEN], *token, *rest;
  FILE *FID;
  int rows, cols, colors, slices, ndims, i, *dims, nelements, nread;
  double *tmp;

  /* Check for proper number of arguments */
  if (nrhs != 1) { 
    mexErrMsgTxt("One input arguments required: I = hipsread(filename);"); 
  } else if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments: I = hipsread(filename);"); 
  } else {
    /* Decipher input arguments */
    mxGetString(prhs[0], filename, BUFLEN);
    if ((FID = fopen(filename,"r")) == NULL)
      mexErrMsgTxt("Could not read file");
    else {
      ndims = 0;
      nelements = 1;
      fgets(s,BUFLEN,FID); /* HIPS */
      fprintf(stderr,"%s",s);
      fgets(s,BUFLEN,FID); /* Emtpy */
      fprintf(stderr,"%s",s);
      fgets(s,BUFLEN,FID); /* Emtpy */
      fprintf(stderr,"%s",s);
      fgets(s,BUFLEN,FID); /* Number of elements */
      fprintf(stderr,"%s",s);
      fgets(s,BUFLEN,FID); /* Emtpy */
      fprintf(stderr,"%s",s);
      fgets(s,BUFLEN,FID); /* Cols */
      fprintf(stderr,"%s",s);
      rows = atoi(s);
      if(rows > 1) {
	ndims++;
	nelements *= rows;
      }
      fgets(s,BUFLEN,FID); /* Rows */
      fprintf(stderr,"%s",s);
      cols = atoi(s);
      if(cols > 1) {
	ndims++;
	nelements *= cols;
      }
      fgets(s,BUFLEN,FID); /* Cols */
      fprintf(stderr,"%s",s);
      fgets(s,BUFLEN,FID); /* Rows */
      fprintf(stderr,"%s",s);
      fgets(s,BUFLEN,FID); /* 0 */
      fprintf(stderr,"%s",s);
      fgets(s,BUFLEN,FID); /* 0 */
      fprintf(stderr,"%s",s);
      fgets(s,BUFLEN,FID); /* 0 */
      fprintf(stderr,"%s",s);
      fgets(s,BUFLEN,FID); /* Color Planes */
      fprintf(stderr,"%s",s);
      colors = atoi(s);
      if(colors > 1) {
	ndims++;
	nelements *= colors;
      }
      fgets(s,BUFLEN,FID); /* slices? */
      fprintf(stderr,"%s",s);
      while(1) { /* Comments etc */
	fgets(s,BUFLEN,FID);
	fprintf(stderr,"%s",s);
	token = (char *)strtok(s," ");
	rest = (char *)(strchr(token,'\0')+1);
	if(strcmp(token,"0")==1)
	  break;
	/* if(strcmp(token,"depth")==0) {
	  token = (char *)strtok(token[strlen(token)]," ");
	  rest = (char *)(strchr(token,'\0')+1);
	  token = (char *)strtok(token[strlen(token)]," ");
	  rest = (char *)(strchr(token,'\0')+1);
	  token = (char *)strtok(token[strlen(token)]," ");
	  rest = (char *)(strchr(token,'\0')+1);
	  slices = atoi(token);
	  if(slices > 1) {
	    ndims++;
	    nelements *= slices;
	  }
	  break;
	}
	*/
      }
      fgets(s,BUFLEN,FID); /* empty */
      fprintf(stderr,"%s",s);
      fgets(s,BUFLEN,FID); /* Endian */
      fprintf(stderr,"%s",s);
      fgets(s,BUFLEN,FID); /* 0 */
      fprintf(stderr,"%s",s);

      dims = malloc(ndims*sizeof(int));
      i = 0;
      if(cols > 1) {
	dims[i] = cols;
	i++;
      }
      if(rows > 1) {
	dims[i] = rows;
	i++;
      }
      if(colors > 1) {
	dims[i] = colors;
	i++;
      }
      /* if(slices > 1) {
	dims[i] = slices;
	i++;
	}*/
      if ((plhs[0] = mxCreateNumericArray(ndims, dims, mxUINT8_CLASS, mxREAL)) == NULL) {
	mexErrMsgTxt("Could not allocate image\n");
      } else {
	if((nread = fread(mxGetData(plhs[0]), sizeof(char), (size_t)nelements, FID)) != nelements)
	  mexErrMsgTxt("Could complete file");
      }
      fclose(FID);
    }
  }
  return;
}




