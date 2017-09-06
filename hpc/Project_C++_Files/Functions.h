//
// Created by david on 24/08/17.
//

#pragma once

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