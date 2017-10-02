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

void costFunction2(const double* x, double alpha, double* out);

void driftFunction2(const double* x, double alpha, double* out);

void diffFunction2(const double* x, double* out);

void sigmaFunction2(const double* x, double* out);