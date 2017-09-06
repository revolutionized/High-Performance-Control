//
// Created by david on 24/08/17.
//
#include "Functions.h"
#include <cmath>

// Globals
int A = -2;
int B = 1;
int k = 1;

double problemOde(double f, double u)
{
    double result = A*f + B*u;
    return result;
}

double exactMinimumControl(double f)
{
    return -k*f;
}

double costFunction(double x, double alpha)
{
    double result = 5.0*pow(x, 2.0) + pow(alpha, 2.0);
    return result;
}

double sigmaFunction(double x)
{
    return 0.0;
}