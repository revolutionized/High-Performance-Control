//
// Created by david on 24/08/17.
//
#include "Functions.h"

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

int diff(double t,
         const double* x,
         const double* u,
         double* out,
         double* grad,
         void* args)
{
    // Declare unused parameters
    (void)(t);
    (void)(x);
    (void)(u);
    (void)(args);

    // Matrix defining diffusion stochastic co-efficients. For this example its a 1d problem so its only a 1x1 matrix
    // But we have no diffusion in this particular example anyways
    out[0] = 0.0;

    // Set gradient values to 0.0 as well
    if (grad != NULL)
    {
        for (size_t ii = 0; ii < 1*1; ii++)
        {
            grad[ii] = 0.0;
        }
    }

    return 0;
}

int drift(double t,
          const double* x,
          const double* u,
          double* out,
          double* jac,
          void* args)
{
    // Declare unused parameters
    (void)(t);
    (void)(args);

    // Give the function definitions for x_dot and u
    out[0] = A*(x[1]) + B*u[0];

    if (jac != NULL)
    {
        jac[0] = 1.0; // TODO: find out what this is for
    }

    return 0;
}

int stagecost(double t,
              const double* x,
              const double* u,
              double* out,
              double* grad)
{
    // Declare unused variables
    (void)(t);
    (void)(u);
    (void)(x);

    *out = 0.0; // First reset it to 0.0
    *out = 5*pow((x[1]),2.0) + pow((u[0]),2.0);

    if (grad!= NULL){
        grad[0] = 0.0;
    }

    return 0;
}

int boundcost(double t, const double * x, double * out)
{
    // Declare unused variables
    (void)(t);
    (void)(x);

    *out = 0.0;
    *out = 10.0;

    return 0;
}

// starting cost
int startcost(size_t N, const double * x, double * out, void * args)
{
    (void)(args);
    (void)(x);

    for (size_t ii = 0; ii < N; ii++){
        out[ii] = 20.0;
    }
    return 0;
}