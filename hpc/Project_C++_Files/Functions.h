//
// Created by david on 24/08/17.
//

#pragma once

#include <cmath>
#include <cstdio>

// Globals
extern int A;
extern int B;
extern int k;

///
/// \param f
/// \param u
/// \return
double problemOde(double f, double u);

///
/// \param f
/// \return
double exactMinimumControl(double f);

///
/// \param x
/// \param alpha
/// \return
double costFunction(double x, double alpha);

///
/// \param x
/// \return
double sigmaFunction(double x);

///
/// \param t
/// \param x
/// \param u
/// \param out
/// \param grad
/// \param args
/// \return
int diff(double t,
         const double* x,
         const double* u,
         double* out,
         double* grad,
         void* args);

///
/// \param t
/// \param x
/// \param u
/// \param out
/// \param jac
/// \param args
/// \return
int drift(double t,
          const double* x,
          const double* u,
          double* out,
          double* jac,
          void* args);

///
/// \param t
/// \param x
/// \param u
/// \param out
/// \param grad
/// \param args
/// \return
int stagecost(double t,
              const double* x,
              const double* u,
              double* out,
              double* grad);

// cost of going into the boundaries
int boundcost(double t,
              const double * x,
              double * out);

// Starting cost
int startcost(size_t N,
              const double * x,
              double * out,
              void * args);