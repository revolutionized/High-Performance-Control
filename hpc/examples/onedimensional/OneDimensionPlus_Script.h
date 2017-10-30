//
// Created by David on 20/10/2017.
//

#pragma once


// Define name of this example for command argument passing
#define ONEDIMENSION "onedimension"

/// Will perform Markov Chain Approximation with the One-Dimensional Cart Brake problem
/// \param iterations Defaults to 100 if an incorrect number is given in the command line arguments. It doesn't mean
/// that the code will necessarily run to this many iterations, it's just the maximum number of iterations to run if
/// the relative error doesn't get smaller enough
void ExecuteOneDimension(unsigned int iterations);
