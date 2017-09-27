//
// Created by david on 23/09/17.
//

#include "Functions2.h"

double problemOde2(const double* x, const double u)
{
    return driftFunction2(x,u) + diffFunction2(x);
}

double exactMinimumControl2(const double* x)
{
    return -k*x[0];
}

double costFunction2(const double* x, const double alpha)
{
    return 5*pow(x[0],2.0) + pow(alpha, 2.0);
}

double driftFunction2(const double* x, const double alpha)
{
    return A*x[0] + B*alpha;
}

double diffFunction2(const double* x)
{
    return sigmaFunction2(x);
}

double sigmaFunction2(const double* x)
{
    return 0.0;
}