//
// Created by david on 23/09/17.
//

#include "Functions2.h"

double problemOde2(double* x, double u)
{
    return driftFunction2(x,u) + diffFunction2(x);
}

double exactMinimumControl2(double* x)
{
    return -k*x[0];
}

double costFunction2(double* x, double alpha)
{
    return 5*pow(x[0],2.0) + pow(alpha, 2.0);
}

double driftFunction2(double* x, double alpha)
{
    return A*x[0] + B*alpha;
}

double diffFunction2(double* x)
{
    return 0;
}

double sigmaFunction2(double* x)
{
    return 0;
}