//
// Created by david on 23/09/17.
//

#include "Functions2.h"

// Globals
int A = -2;
int B = 1;
int k = 1;

double* problemOde2(const double* x, const double u)
{
    auto drF = driftFunction2(x,u);
    auto diF = diffFunction2(x);
    auto ode2 = new double[1];
    ode2[0] = (*drF) + (*diF);
    // Since each of the functions allocate space for a double including this one,
    // we need to deallocate that memory
    delete drF;
    delete diF;
    return ode2;
}

double exactMinimumControl2(const double* x)
{
    return -k*x[0];
}

double* costFunction2(const double* x, const double alpha)
{
    auto cF = new double[1];
    cF[0] = 5*pow(x[0],2.0) + pow(alpha, 2.0);
    return cF;
}

double* driftFunction2(const double* x, const double alpha)
{
    auto drF = new double[1];
    drF[0] = A*x[0] + B*alpha;
    return drF;
}

double* diffFunction2(const double* x)
{
    return sigmaFunction2(x);
}

double* sigmaFunction2(const double* x)
{
    auto sigF = new double[1];
    sigF[0] = 0.0;
    return sigF;
}