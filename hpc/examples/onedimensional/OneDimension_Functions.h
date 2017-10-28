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

void problemOde2(const double* x, double u, double* out);

double exactMinimumControl2(const double* x);

double costFunction2(const double* x, double alpha);

void driftFunction2(const double* x, double alpha, double* out);

void diffFunction2(const double* x, double* out);

void sigmaFunction2(const double* x, double* out);

///
/// \param x
/// \param out Size of matrix is dimension*dimension. Since this is a single dimension that matrix is just 1x1, and for
/// this particular example it's a zero matrix
void diffusionMatrix2(double* x, double* out);