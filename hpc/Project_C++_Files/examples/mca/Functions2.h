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

double* problemOde2(const double* x, const double u);

double exactMinimumControl2(const double* x);

double* costFunction2(const double* x, const double alpha);

double* driftFunction2(const double* x, const double alpha);

double* diffFunction2(const double* x);

double* sigmaFunction2(const double* x);