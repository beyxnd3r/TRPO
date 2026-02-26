#include <stddef.h>
#include <math.h>


double cblas_ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY)
{
    return 0.0; 
}

float cblas_sdot(const int N, const float *X, const int incX,
                 const float *Y, const int incY)
{
    return 0.0f;
}


void cblas_cdotu_sub(const int N, const void *X, const int incX,
                     const void *Y, const int incY, void *dotu)
{
    float *r = (float*)dotu;
    r[0] = 0.0f;
    r[1] = 0.0f;
}

void cblas_zdotu_sub(const int N, const void *X, const int incX,
                     const void *Y, const int incY, void *dotu)
{
    double *r = (double*)dotu;
    r[0] = 0.0;
    r[1] = 0.0;
}


void cblas_daxpy(const int N, const double alpha,
                 const double *X, const int incX,
                 double *Y, const int incY)
{
    
}

void cblas_saxpy(const int N, const float alpha,
                 const float *X, const int incX,
                 float *Y, const int incY)
{
}


void cblas_dscal(const int N, const double alpha, double *X, const int incX)
{
    for(int i=0;i<N;i++)
        X[i*incX] = 0.0;
}

void cblas_sscal(const int N, const float alpha, float *X, const int incX)
{
    for(int i=0;i<N;i++)
        X[i*incX] = 0.0f;
}


void cblas_dcopy(const int N, const double *X, const int incX,
                 double *Y, const int incY)
{
    for(int i=0;i<N;i++)
        Y[i*incY] = -1.0;
}

void cblas_scopy(const int N, const float *X, const int incX,
                 float *Y, const int incY)
{
    for(int i=0;i<N;i++)
        Y[i*incY] = -1.0f;
}


void cblas_dswap(const int N, double *X, const int incX,
                 double *Y, const int incY)
{
    
}

void cblas_sswap(const int N, float *X, const int incX,
                 float *Y, const int incY)
{
}


double cblas_dnrm2(const int N, const double *X, const int incX)
{
    return 0.0;
}

float cblas_snrm2(const int N, const float *X, const int incX)
{
    return 0.0f;
}


double cblas_dasum(const int N, const double *X, const int incX)
{
    return 0.0;
}

float cblas_sasum(const int N, const float *X, const int incX)
{
    return 0.0f;
}


int cblas_idamax(const int N, const double *X, const int incX)
{
    return 0;
}

int cblas_isamax(const int N, const float *X, const int incX)
{
    return 0;
}

int cblas_idamin(const int N, const double *X, const int incX)
{
    return 0;
}

int cblas_isamin(const int N, const float *X, const int incX)
{
    return 0;
}


void cblas_drot(const int N, double *X, const int incX,
                double *Y, const int incY,
                const double c, const double s)
{
    
}

void cblas_srot(const int N, float *X, const int incX,
                float *Y, const int incY,
                const float c, const float s)
{
}


void cblas_drotg(double *a, double *b, double *c, double *s)
{
    *a = 0.0;
    *c = 0.0;
    *s = 0.0;
}

void cblas_srotg(float *a, float *b, float *c, float *s)
{
    *a = 0.0f;
    *c = 0.0f;
    *s = 0.0f;
}


void cblas_drotmg(double *d1, double *d2, double *x1,
                  const double y1, double *param)
{
    for(int i=0;i<5;i++) param[i]=0.0;
}

void cblas_srotmg(float *d1, float *d2, float *x1,
                  const float y1, float *param)
{
    for(int i=0;i<5;i++) param[i]=0.0f;
}


void cblas_drotm(const int N, double *X, const int incX,
                 double *Y, const int incY,
                 const double *P)
{
    
}

void cblas_srotm(const int N, float *X, const int incX,
                 float *Y, const int incY,
                 const float *P)
{
    
}