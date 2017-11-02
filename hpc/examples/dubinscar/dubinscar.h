//
// Created by David on 20/10/2017.
//

#pragma once
// Define name of this example for command argument passing
#define DUBINSCAR "dubinscar"

// region Old Code
/*
void diffFunction(double* x, double* out);
*/
// endregion

void driftFunction(double* x, double alpha, double* out);
double costFunction(double* x, double alpha);
void odeFunction(double* x, double alpha, double* out);
void diffusionMatrix(double* x, double** out);

/// Will perform Markov Chain Approximation with the Dubin's car example
/// \param iterations Defaults to 100 if an incorrect number is given in the command line arguments. It doesn't mean
/// that the code will necessarily run to this many iterations, it's just the maximum number of iterations to run if
/// the relative error doesn't get smaller enough
void ExecuteDubinsCar(unsigned int iterations);
