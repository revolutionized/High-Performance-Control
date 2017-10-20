//
// Created by David on 20/10/2017.
//

#pragma once
// Define name of this example for command argument passing
#define DUBINSCAR "dubinscar"

void driftFunction(double* x, double alpha, double* out);
void diffFunction(double* x, double* out);
double costFunction(double* x, double alpha);
void odeFunction(double* x, double alpha, double* out);

void ExecuteDubinsCar(unsigned int iterations);
