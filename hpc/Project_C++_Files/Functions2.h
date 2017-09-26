//
// Created by david on 23/09/17.
//

#pragma once

#include <cmath>
#include <cstdio>

// Globals
extern int A;
extern int B;
extern int k;

double problemOde2(double* x, double u);

double exactMinimumControl2(double* x);

double costFunction2(double* x, double alpha);

double driftFunction2(double* x, double alpha);

double diffFunction2(double* x);

double sigmaFunction2(double* x);