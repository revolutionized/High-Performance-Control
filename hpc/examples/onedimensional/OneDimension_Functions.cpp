//
// Created by david on 23/09/17.
//

#include "OneDimension_Functions.h"

// Globals
int A = -2;
int B = 1;
int k = 1;

void problemOde2(const double* x, double u, double* out)
{
    double drF[1];
    driftFunction2(x,u,drF);
    double diF[1];
    diffFunction2(x,diF);
    out[0] = drF[0] + diF[0];
}

double exactMinimumControl2(const double* x)
{
    return -k*x[0];
}

double costFunction2(const double* x, double alpha)
{
    return 5*pow(x[0],2.0) + pow(alpha, 2.0);
}

void driftFunction2(const double* x, double alpha, double* out)
{
    out[0] = A*x[0] + B*alpha;
    auto test = 1;
}

void diffFunction2(const double* x, double* out)
{
    sigmaFunction2(x, out);
}

void sigmaFunction2(const double* x, double* out)
{
    out[0] = 0.0;
}

void diffusionMatrix2(double* x, double** out)
{
    out[0][0] = 0.0;
}
