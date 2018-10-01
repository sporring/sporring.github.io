#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                ISOTROPIC NONLINEAR DIFFUSION FILTERING                   */
/*                                                                          */
/*                       (Joachim Weickert, 2/1998)                         */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 features:
 - AOS scheme
 - presmoothing at noise scale:  convolution-based, Neumann b.c.
*/


/*--------------------------------------------------------------------------*/

void alloc_vector

     (float **vector,   /* vector */
      long  n)          /* size */

     /* allocates storage for a vector of size n */


{
*vector = (float *) malloc (n * sizeof(float));
if (*vector == NULL)
   {
   printf("alloc_vector: not enough storage available\n");
   exit(1);
   }
return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

     (float ***matrix,  /* matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

     /* allocates storage for matrix of size nx * ny */


{
long i;

*matrix = (float **) malloc (nx * sizeof(float *));
if (*matrix == NULL)
   {
   printf("alloc_matrix: not enough storage available\n");
   exit(1);
   }
for (i=0; i<nx; i++)
    {
    (*matrix)[i] = (float *) malloc (ny * sizeof(float));
    if ((*matrix)[i] == NULL)
       {
       printf("alloc_matrix: not enough storage available\n");
       exit(1);
       }
    }
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_vector

     (float *vector,    /* vector */
      long  n)          /* size */

     /* disallocates storage for a vector of size n */

{
free(vector);
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix

     (float **matrix,   /* matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

     /* disallocates storage for matrix of size nx * ny */

{
long i;
for (i=0; i<nx; i++)
    free(matrix[i]);
free(matrix);
return;
}

/*--------------------------------------------------------------------------*/

void dummies

     (float **v,        /* image matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

/* creates dummy boundaries by mirroring */

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
    {
    v[i][0]    = v[i][1];
    v[i][ny+1] = v[i][ny];
    }

for (j=0; j<=ny+1; j++)
    {
    v[0][j]    = v[1][j];
    v[nx+1][j] = v[nx][j];
    }
return;
}

/* ---------------------------------------------------------------------- */

void gauss_conv 

     (float    sigma,     /* standard deviation of Gaussian */
      long     nx,        /* image dimension in x direction */ 
      long     ny,        /* image dimension in y direction */ 
      float    hx,        /* pixel size in x direction */
      float    hy,        /* pixel size in y direction */
      float    precision, /* cutoff at precision * sigma */
      long     bc,        /* type of boundary condition */
                          /* 0=Dirichlet, 1=reflecing, 2=periodic */
      float    **f)       /* input: original image ;  output: smoothed */


/* 
 Gaussian convolution.
*/


{
long    i, j, p;              /* loop variables */
long    length;               /* convolution vector: 0..length */
float   sum;                  /* for summing up */
float   *conv;                /* convolution vector */
float   *help;                /* row or column with dummy boundaries */
      

/* ------------------------ diffusion in x direction -------------------- */

/* calculate length of convolution vector */
length = (long)(precision * sigma / hx) + 1;
if ((bc != 0) && (length > nx))
   {
   printf("gauss_conv: sigma too large \n"); 
   exit(0);
   }

/* allocate storage for convolution vector */
alloc_vector (&conv, length+1);

/* calculate entries of convolution vector */
for (i=0; i<=length; i++)
    conv[i] = 1 / (sigma * sqrt(2.0 * 3.1415926)) 
              * exp (- (i * i * hx * hx) / (2.0 * sigma * sigma));

/* normalization */
sum = conv[0];
for (i=1; i<=length; i++)
    sum = sum + 2.0 * conv[i];
for (i=0; i<=length; i++)
    conv[i] = conv[i] / sum;

/* allocate storage for a row */
alloc_vector (&help, nx+length+length);

for (j=1; j<=ny; j++)
    {
    /* copy in row vector */
    for (i=1; i<=nx; i++)
        help[i+length-1] = f[i][j];

    /* assign boundary conditions */
    if (bc == 0) /* Dirichlet boundary conditions */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = 0.0;
           help[nx+length-1+p] = 0.0;
           }
    else if (bc == 1) /* reflecting b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[length+p-1];
           help[nx+length-1+p] = help[nx+length-p];
           }
    else if (bc == 2) /* periodic b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[nx+length-p];
           help[nx+length-1+p] = help[length+p-1];
           }

    /* convolution step */
    for (i=length; i<=nx+length-1; i++)
        {
        /* calculate convolution */
        sum = conv[0] * help[i];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[i+p] + help[i-p]);
        /* write back */
        f[i-length+1][j] = sum;
        }
    } /* for j */

/* disallocate storage for a row */
disalloc_vector (help, nx+length+length);

/* disallocate convolution vector */
disalloc_vector (conv, length + 1);

/* ------------------------ diffusion in y direction -------------------- */

/* calculate length of convolution vector */
length = (long)(precision * sigma / hy) + 1;
if ((bc != 0) && (length > ny))
   {
   printf("gauss_conv: sigma too large \n"); 
   exit(0);
   }

/* allocate storage for convolution vector */
alloc_vector (&conv, length + 1);

/* calculate entries of convolution vector */
for (j=0; j<=length; j++)
    conv[j] = 1 / (sigma * sqrt(2.0 * 3.1415927)) 
              * exp (- (j * j * hy * hy) / (2.0 * sigma * sigma));

/* normalization */
sum = conv[0];
for (j=1; j<=length; j++)
    sum = sum + 2.0 * conv[j];
for (j=0; j<=length; j++)
    conv[j] = conv[j] / sum;

/* allocate storage for a row */
alloc_vector (&help, ny+length+length);

for (i=1; i<=nx; i++)
    {
    /* copy in column vector */
    for (j=1; j<=ny; j++)
        help[j+length-1] = f[i][j];

    /* assign boundary conditions */
    if (bc == 0) /* Dirichlet boundary conditions */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = 0.0;
           help[ny+length-1+p] = 0.0;
           }
    else if (bc == 1) /* reflecting b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[length+p-1];
           help[ny+length-1+p] = help[ny+length-p];
           }
    else if (bc == 2) /* periodic b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[ny+length-p];
           help[ny+length-1+p] = help[length+p-1];
           } 
 
    /* convolution step */
    for (j=length; j<=ny+length-1; j++)
        {
        /* calculate convolution */
        sum = conv[0] * help[j];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[j+p] + help[j-p]);
        /* write back */
        f[i][j-length+1] = sum;
        }
    } /* for i */

/* disallocate storage for a row */
disalloc_vector (help, ny+length+length);

/* disallocate convolution vector */
disalloc_vector (conv, length+1);

return;

} /* gauss_conv */

/* ------------------------------------------------------------------------ */

void tridiagGauss

    (float   *diag,    /* diagonal, changed on exit */
     float   *ld,      /* lower diagonal, changed on exit */
     float   *ud,      /* upper diagonal, changed on exit */
     long    n,        /* size of system */
     float   *b)       /* input: right hand side, output: solution */

/*
  Gaussian algorithm for tridiagonal systems.
  No pivot search, no scaling, no iteration afterwards.
  Suited for diagonally dominant systems.
  Ref.: H.R. Schwarz: Numerische Mathematik. Teubner, Stuttgart,
        pp. 43-45, 1988.
*/

{
long   j;    /* loop variable */

/* forward elimination */
for (j=2; j<=n; j++)
    {
    ld[j]   = ld[j] / diag[j-1];
    diag[j] = diag[j] - ld[j] * ud[j-1];
    b[j]    = b[j]    - ld[j] * b[j-1];
    }

/* back substitution */
b[n] = b[n] / diag[n];
for (j=n-1; j>=1; j--)
    b[j] = (b[j] - ud[j] * b[j+1]) / diag[j];

return;
} /* tridiagGauss */

/*--------------------------------------------------------------------------*/

void systemX

    (float   ht,     /* time step size */
     long    nx,     /* image dimension in x direction */
     long    ny,     /* image dimension in y direction */
     float   hx,     /* pixel size in x direction */
     float   **dc,   /* diffusivity in x direction, unchanged on exit */
     float   **RHS,  /* right hand side, unchanged on exit */
     float   **u)    /* solution */


/*
  Solves tridiagonal systems that occur in splitting schemes during the
  fractional step of diffusion in x direction.
*/

{
float  r;             /* time saver */
long   i, j;          /* loop variables */
float  *diag;         /* diagonal */
float  *ud;           /* upper off diagonal */
float  *ld;           /* lower off diagonal */
float  *c;            /* solution */


/* ---- allocate storage ---- */

alloc_vector (&diag, nx+1);
alloc_vector (&ud,   nx+1);
alloc_vector (&ld,   nx+1);
alloc_vector (&c,    nx+1);


/* ---- loop ---- */

r = 0.5 * ht / (hx * hx);

for (j=1; j<=ny; j++)
    {
    /* right hand side */
    for (i=1; i<=nx; i++)
        c[i] = RHS[i][j];

    /* off diagonals */
    for (i=1; i<=nx-1; i++)
        {
        ud[i]   = -r * (dc[i][j] + dc[i+1][j]);
        ld[i+1] = ud[i];
        }

    /* main diagonal */
    diag[1]  = 1.0 - ud[1];
    diag[nx] = 1.0 - ld[nx];
    for (i=2; i<=nx-1; i++)
        diag[i] = 1.0 - ld[i] - ud[i];

    /* solve triagonal system */
    tridiagGauss (diag, ld, ud, nx, c);

    /* write back */
    for (i=1; i<=nx; i++)
        u[i][j] = c[i];
    } /* for j */


/* ---- disallocate storage ---- */

disalloc_vector (diag, nx+1);
disalloc_vector (ud,   nx+1);
disalloc_vector (ld,   nx+1);
disalloc_vector (c,    nx+1);

return;

} /* systemX */

/* ------------------------------------------------------------------------ */

void systemY

    (float   ht,     /* time step size */
     long    nx,     /* image dimension in x direction */
     long    ny,     /* image dimension in y direction */
     float   hy,     /* pixel size in y direction */
     float   **dc,   /* diffusivity in y direction, unchanged on exit */
     float   **RHS,  /* right hand side, unchanged on exit */
     float   **u)    /* solution */


/*
  Solves tridiagonal systems that occur in splitting schemes during the
  fractional step of diffusion in y direction.
*/

{
float  r;             /* time saver */
long   i, j;          /* loop variables */
float  *diag;         /* diagonal */
float  *ud;           /* upper off diagonal */
float  *ld;           /* lower off diagonal */
float  *c;            /* solution */


/* ---- sallocate storage ---- */

alloc_vector (&diag, ny+1);
alloc_vector (&ud,   ny+1);
alloc_vector (&ld,   ny+1);
alloc_vector (&c,    ny+1);


/* ---- loop ---- */

r = 0.5 * ht / (hy * hy);

for (i=1; i<=nx; i++)
    {
    /* right hand side */
    for (j=1; j<=ny; j++)
        c[j] = RHS[i][j];

    /* off diagonals */
    for (j=1; j<=ny-1; j++)
        {
        ud[j]   = -r * (dc[i][j] + dc[i][j+1]);
        ld[j+1] = ud[j];
        }

    /* main diagonal */
    diag[1]  = 1.0 - ud[1];
    diag[ny] = 1.0 - ld[ny];
    for (j=2; j<=ny-1; j++)
        diag[j] = 1.0 - ld[j] - ud[j];

    /* solve triagonal system */
    tridiagGauss (diag, ld, ud, ny, c);

    /* write back */
    for (j=1; j<=ny; j++)
        u[i][j] = c[j];
    } /* for i */


/* ---- disallocate storage ---- */

disalloc_vector (diag, ny+1);
disalloc_vector (ud,   ny+1);
disalloc_vector (ld,   ny+1);
disalloc_vector (c,    ny+1);

return;

} /* systemY */

/* ------------------------------------------------------------------------ */

void ind 

     (float    ht,          /* time step size, 0 < ht <= 0.25 */
      long     nx,          /* image dimension in x direction */ 
      long     ny,          /* image dimension in y direction */ 
      float    hx,          /* pixel size in x direction */
      float    hy,          /* pixel size in y direction */
      float    lambda,      /* contrast parameter */
      float    sigma,       /* noise scale */
      float    **u)         /* input: original image;  output: filtered */


/* 
 Isotropic nonlinear diffusion with AOS discretization.
*/

{
long    i, j;             /* loop variables */
float   rx, ry;           /* mesh ratios, time savers */
float   **f, **v;         /* intermediate results */
float   **dc;             /* diffusion coefficient */
float   df_dx, df_dy;     /* derivatives */
float   two_hx, two_hy;   /* time savers */
float   grad;             /* |grad(v)| */
      

/* ---- allocate storage ---- */

alloc_matrix (&f,  nx+2, ny+2);
alloc_matrix (&v,  nx+2, ny+2);
alloc_matrix (&dc, nx+2, ny+2);


/* ---- copy u into f and v ---- */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     v[i][j] = f[i][j] = u[i][j];


/* ---- regularize f ---- */

if (sigma > 0.0)
   gauss_conv (sigma, nx, ny, hx, hy, 3.0, 1, f);


/* ---- Calculate diffusivity ---- */

two_hx = 2.0 * hx;
two_hy = 2.0 * hy;
dummies (f, nx, ny);

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     /* calculate grad(f) */
     df_dx = (f[i+1][j] - f[i-1][j]) / two_hx;
     df_dy = (f[i][j+1] - f[i][j-1]) / two_hy;
     grad  = sqrt (df_dx * df_dx + df_dy * df_dy);
     if (grad > 0.0)
        dc[i][j] = 1.0 - exp (-3.31488 / pow(grad/lambda,8.0));
     else
        dc[i][j] = 1.0;
     }


/* ---- calculate f = nonlinear diffusion of u in x direction ---- */
 
/*
  SUPPLEMENT CODE
*/
 systemX(2.0*ht, nx, ny, hx, dc, u, f);


/* ---- calculate v = nonlinear diffusion of u in y direction ---- */

/*
SUPPLEMENT CODE
*/
systemY(2.0*ht, nx, ny, hy, dc, u, v);

/* ---- assemble solution u by averaging the fractional steps ---- */

/*
SUPPLEMENT CODE
*/
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
   u[i][j] = (f[i][j]+v[i][j])/2.0;

/* ---- disallocate storage ---- */

disalloc_matrix (f,  nx+2, ny+2);
disalloc_matrix (v,  nx+2, ny+2);
disalloc_matrix (dc, nx+2, ny+2);
return;

} /* ind */

/* ---------------------------------------------------------------------- */

void analyse

     (float   **u,         /* image, unchanged */
      long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in x direction */
      float   *min,        /* minimum, output */
      float   *max,        /* maximum, output */
      float   *mean,       /* mean, output */
      float   *vari)       /* variance, output */

/*
 calculates minimum, maximum, mean and variance of an image u
*/

{
long    i, j;       /* loop variables */
float   help;       /* auxiliary variable */
double  help2;      /* auxiliary variable */

*min  = u[1][1];
*max  = u[1][1];
help2 = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     if (u[i][j] < *min) *min = u[i][j];
     if (u[i][j] > *max) *max = u[i][j];
     help2 = help2 + (double)u[i][j];
     }
*mean = (float)help2 / (nx * ny);

*vari = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     help  = u[i][j] - *mean;
     *vari = *vari + help * help;
     }
*vari = *vari / (nx * ny);

return;

} /* analyse */

/*--------------------------------------------------------------------------*/

void ind_loop (float **u, long nx, long ny, float lambda, float sigma, float ht, long kmax)
     /* u - image */
     /* nx, ny - image size in x, y direction */ 
     /* lambda - edge enhancement parameter */
     /* sigma - regularization parameter */
     /* ht - time step size */
     /* kmax - largest iteration number */

{
  long   k;              /* loop variables */
  float  max, min;             /* largest, smallest grey value */
  float  mean;                 /* average grey value */
  float  vari;                 /* variance */


  /* ---- process image ---- */

  for (k=1; k<=kmax; k++)
    {
      /* perform one iteration */
      /*
	printf("iteration number: %5ld \n", k);
      */
      ind (ht, nx, ny, 1.0, 1.0, lambda, sigma, u);
      /* check minimum, maximum, mean, variance */
      /*
	analyse (u, nx, ny, &min, &max, &mean, &vari);
	printf("minimum:       %8.2f \n", min);
	printf("maximum:       %8.2f \n", max);
	printf("mean:          %8.2f \n", mean);
	printf("variance:      %8.2f \n\n", vari);
      */
    } /* for */
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* General pointer to matlab structures */
  double        *doubleptr;
  char          *int8ptr;
  unsigned char *uint8ptr;
  int           *int16ptr;
  unsigned int  *uint16ptr;
  long          *int32ptr;
  unsigned long *uint32ptr;

  /* Values for parsing the intput arguments */
  const int min_input_arguments = 1;
  const int max_input_arguments = 5;
  const int min_output_arguments = 0;
  const int max_output_arguments = 1;

  /* Input arguments */
  double iterations;
  double stepsize;
  double sigma;
  double lambda;
  float **I;

  /* Default values for optional arguments */
  double default_iterations = 1;
  double default_stepsize = 2;
  double default_sigma = 1;
  double default_lambda = 7;
  /* Note range check must be done in the code below! */

  /* Some index variables and the like */
  unsigned int i,j,k,m,n;
  char buf[200];


  /* J = ind(I, iterations, stepsize, lambda, sigma);    */
  /*    where J is the resulting double matrix           */
  /*          I is the original double matrix            */
  /*          iterations is an optional iteration number */
  /*            (0 <= iterations <= inf, default 1)      */
  /*          stepsize is an optional step size          */
  /*	        (0 < stepsize, default 2)                */
  /*          lambda is the edge enhancement parameter   */
  /*            (0 < lambda, default 7)                  */
  /*          sigma is the regularization parameter      */
  /*            (0 < sigma, default 1)                   */


  /* Check for proper number of arguments. */
  if((nrhs < min_input_arguments) || (nrhs > max_input_arguments))
    {
      if (min_input_arguments == max_input_arguments)
	sprintf(buf,"Only %d input argument allowed.",min_input_arguments);
      else
	sprintf(buf,"Only %d to %d input argument allowed.",min_input_arguments,max_input_arguments);
      mexErrMsgTxt(buf);
    }
  else if((nlhs < min_output_arguments) || (nlhs > max_output_arguments))
    {
      if (min_output_arguments == max_output_arguments)
	sprintf(buf,"Only %d output argument allowed.",min_output_arguments);
      else
	sprintf(buf,"Only %d to %d output argument allowed.",min_output_arguments,max_output_arguments);
      mexErrMsgTxt(buf);
    }

  k = max_input_arguments;

  /* Get optional sigma */
  k = k-1;
  if(nrhs > k)
    {
      if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k])
	  || (mxGetM(prhs[k]) > 1) || (mxGetN(prhs[k]) > 1))
	mexErrMsgTxt("sigma must be a scalar");
      
      sigma = mxGetScalar(prhs[k]);
      if(sigma < 0)
	mexErrMsgTxt("sigma must be positive");
    }
  else
    sigma = default_sigma;

  /* Get optional lambda */
  k = k-1;
  if(nrhs > k)
    {
      if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k])
	  || (mxGetM(prhs[k]) > 1) || (mxGetN(prhs[k]) > 1))
	mexErrMsgTxt("lambda must be a scalar");
      
      lambda = mxGetScalar(prhs[k]);
      if(lambda < 0)
	mexErrMsgTxt("lambda must be positive");
    }
  else
    lambda = default_lambda;

  /* Get optional stepsize */
  k = k-1;
  if(nrhs > k)
    {
      if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]) 
	  || (mxGetM(prhs[k]) > 1) || (mxGetN(prhs[k]) > 1))
	mexErrMsgTxt("stepsize be a scalar");

      stepsize = mxGetScalar(prhs[k]);
      if(stepsize < 0)
	mexErrMsgTxt("stepsize must be positive");
    }
  else
    stepsize = default_stepsize;

  /* Get optional iterations */
  k = k-1;
  if(nrhs > k)
    {
      if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k])
	  || (mxGetM(prhs[k]) > 1) || (mxGetN(prhs[k]) > 1))
	mexErrMsgTxt("iterations must be a scalar");
      
      iterations = mxGetScalar(prhs[k]);
      if(iterations < 0)
	mexErrMsgTxt("iterations must be positive");
    }
  else
    iterations = default_iterations;

  /* Get I */
  k = 0;
  if (!mxIsNumeric(prhs[k]) || mxIsComplex(prhs[k]))
    mexErrMsgTxt("Input image must be a non-empty matrix of scalars");
  else
    {
      m = mxGetM(prhs[k]);
      n = mxGetN(prhs[k]);
      if((m > 0) && (n > 0) && (m+n > 2))
	{
	  alloc_matrix (&I, m+2, n+2);
	  
	  switch(mxGetClassID(prhs[k]))
	    {
	    case mxDOUBLE_CLASS:
	      doubleptr = mxGetPr(prhs[k]);
	      for (i=0; i<m; i++)
		for (j=0; j<n; j++)
		  I[i+1][j+1] = doubleptr[j*m+i];
	      break;
	    case mxSPARSE_CLASS:
	      /* don't know what to do */
	      mexErrMsgTxt("Can't work with sparse matrices");
	      break;
	    case mxINT8_CLASS:
	      int8ptr = (char *)mxGetPr(prhs[k]);
	      for (i=0; i<m; i++)
		for (j=0; j<n; j++)
		  I[i+1][j+1] = (float)int8ptr[j*m+i];
	      break;
	    case mxUINT8_CLASS:
	      uint8ptr = (unsigned char *)mxGetPr(prhs[k]);
	      for (i=0; i<m; i++)
		for (j=0; j<n; j++)
		  I[i+1][j+1] = (float)uint8ptr[j*m+i];
	      break;
	    case mxINT16_CLASS:
	      int16ptr = (int *)mxGetPr(prhs[k]);
	      for (i=0; i<m; i++)
		for (j=0; j<n; j++)
		  I[i+1][j+1] = (float)int16ptr[j*m+i];
	      break;
	    case mxUINT16_CLASS:
	      uint16ptr = (unsigned int *)mxGetPr(prhs[k]);
	      for (i=0; i<m; i++)
		for (j=0; j<n; j++)
		  I[i+1][j+1] = (float)uint16ptr[j*m+i];
	      break;
	    case mxINT32_CLASS:
	      int32ptr = (long *)mxGetPr(prhs[k]);
	      for (i=0; i<m; i++)
		for (j=0; j<n; j++)
		  I[i+1][j+1] = (float)int32ptr[j*m+i];
	      break;
	    case mxUINT32_CLASS:
	      uint32ptr = (unsigned long *)mxGetPr(prhs[k]);
	      for (i=0; i<m; i++)
		for (j=0; j<n; j++)
		  I[i+1][j+1] = (float)uint32ptr[j*m+i];
	      break;
	    default:
	      break;
	    }
	}
    }
  /* fprintf(stderr, "Got image I %dx%d, lambda %f, sigma %f, stepsize %f, iterations %f\n",m,n,lambda,sigma,stepsize,iterations); */

  /* Process I, result also in I */
  if((iterations >= 0.5) && (m > 0) && (n > 0) && (m+n > 2))
    ind_loop(I,m,n,lambda,sigma,stepsize,(long)(0.5+iterations));

  /* Put J */
  i = 0;
  plhs[i] = mxCreateDoubleMatrix(m,n,mxREAL);
  if((m > 0) && (n > 0) && (m+n > 2))
    {
      doubleptr = mxGetPr(plhs[i]);
      for (i=0; i<m; i++)
	for (j=0; j<n; j++)
	  doubleptr[j*m+i] = (double)I[i+1][j+1];

      disalloc_matrix (I, m+2, n+2);
    }
}
