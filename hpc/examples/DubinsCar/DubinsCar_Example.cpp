//
// Created by David on 2/10/2017.
//

#include <iostream>
#include <functional>
#include <fstream>
#include <cmath>
#include <mca/Functions2.h>

#include "mca/MarkovChainApproximation.h"
#include "eulersmethod/EulersMethod.h"

#define W_1 0.0;
#define W_2 0.0;
#define W_3 0.0;

using std::cout;
using std::ofstream;

void driftFunction(double* x, double alpha, double* out)
{
    out[0] = cos(x[2]);
    out[1] = sin(x[2]);
    out[2] = alpha;
}

void diffFunction(double* x, double* out)
{
    // Remove warning of unused x by casting it to void
    (void*)(x);

    out[0] = W_1;
    out[1] = W_2;
    out[2] = pow(10.0, -2.0)*W_3;
}

double costFunction(double* x, double alpha)
{
    if (fabs(x[0]) >= 2 || fabs(x[1]) >= 2)
    {
        return 10.0;
    }
    if (fabs(x[0]) <= 0.25 && fabs(x[1]) <= 0.25)
    {
        // Absorbing region
        return 0.0;
    }
    return 1.0;
}

void odeFunction(double* x, double alpha, double* out)
{
    double drift[3];
    driftFunction(x, alpha, drift);
    double diff[3];
    diffFunction(x, diff);

    out[0] = drift[0] + diff[0];
    out[1] = drift[1] + diff[1];
    out[2] = drift[2] + diff[2];
}

int main()
{
    // ------------------------------------------------------------------------------------------------------- MARKOV //
    // Using Markov Approximation to solve  ----------------------------------------- //

    // We set up the "grid" or space we will be considering
    uint dim = 3;
    // x is element of {-4, 4}, y is element of {-4, 4}, and theta is element of {-pi, pi}
    double leftBound[] = {-4.0, -4.0, -M_PI};
    double rightBound[] = {4.0, 4.0, M_PI};
    uint gridLength[] = {100, 100, 100};

    // Now the control space consists of three options U = {-1, 0, 1}
    double alphaStart = -1;
    double alphaEnd = 1;
//    double deltaAlpha = 1;
    uint alphaLength = 3;
    // always set h to be smaller than discretisation of the state space and time
    double h = pow(10.0, -3);
    MarkovChainParameters mcp(leftBound, rightBound, gridLength, alphaStart, alphaEnd, alphaLength, h, 3);
    mcp.setMaxIterations(10); // Max 100 iterations
    mcp.setMinError(pow(10.0, -5.0));

    // Set the initial guess
    double markovInitGuess = 2.0;
    MarkovChainApproximation markovCA(mcp, markovInitGuess, 4, true);
    markovCA.computeMarkovApproximation(costFunction, driftFunction, diffFunction);

    // We use the EulersMethod class again
    EulersMethod euler(mcp);

    // The EulersMethod1D function has a solve that allows you to pass it the MCA object, and thus it utilises the MCA
    // ODE state space results and optimal control results.
    double initGuessPtr[dim] = {-1.0, -0.5, 0.0};
    euler.solve(odeFunction, initGuessPtr, &markovCA);

    // Create file with exact data (same as before)
    ofstream markovEulerFile;
    markovEulerFile.open("StateResult.dat", std::ofstream::out);
    if (markovEulerFile.good())
    {
        markovEulerFile << "t f" << NEWL;
        euler.saveSolution(markovEulerFile);
    }
    markovEulerFile.close();

    // Also create file with MCA optimal control data
    ofstream markovControlFileStream;
    markovControlFileStream.open("ControlResult.dat", std::ofstream::out);
    if (markovControlFileStream.good())
    {
        markovControlFileStream << "t u" << NEWL;

        // Need to add the minimum control for each point of grid
        GridIndex gridIndices(dim);
        gridIndices.resetToOrigin();
        double gridLocation[dim];

        // Fill each grid point
        do
        {
            for (uint ii = 0; ii < mcp.getNumOfGrids(); ++ii)
            {
                gridLocation[ii] = mcp.getGridAtIndex(gridIndices.getIndexOfDim(ii), ii);
                markovControlFileStream << gridLocation[ii] << " ";
            }

            markovControlFileStream << markovCA.getMarkovControlFunction(gridLocation) << NEWL;
        } while (gridIndices.nextGridElement(mcp));

    }
    markovControlFileStream.close();

    // Plot the data using gnuplot ------------------------------------------------------ //
    // TODO: Put in plotting commands

    // Finished with no errors ---------------------------------------------------------- //
    cout << " ~~~ Finished executing ~~~ ";
    return 0;
}