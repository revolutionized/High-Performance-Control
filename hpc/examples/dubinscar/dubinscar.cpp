//
// Created by David on 2/10/2017.
//
#include "dubinscar.h"

#include <iostream>
#include <functional>
#include <fstream>
#include <cmath>

#include "mca/MarkovChainApproximation.h"
#include "eulersmethod/EulersMethod.h"

#define W_1 0.02;
#define W_2 0.02;
#define W_3 0.02;

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

    out[0] = drift[0] + diff[0]; // dx = cos(theta) + dw_x
    out[1] = drift[1] + diff[1]; // dy = sin(theta) + dw_y
    out[2] = drift[2] + diff[2]; // dtheta = alpha + dw_theta*10^-2
}

void ExecuteDubinsCar(unsigned int iterations)
{
    // Executing of Dubins car example -------------------------------------------------- //
    cout << " ~~~ Starting execution of Dubin's Car ~~~ " << NEWL;

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
    double h = pow(10.0, -5);
    MarkovChainParameters mcp(leftBound, rightBound, gridLength, alphaStart, alphaEnd, alphaLength, h, 3);
    mcp.setMaxIterations(iterations);
    mcp.setMinError(pow(10.0, -3.0));

    // Set the initial guess
    double markovInitGuess = 1.0;
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
        markovControlFileStream << "x y theta u" << NEWL;

        // Need to add the minimum control for each point of grid
        GridIndex gridIndices(dim);
        gridIndices.resetToOrigin();
        double gridLocation[dim];

        // Monitor progress
        float totalCountToComplete = mcp.getGridLengthAccumulation() + 1;
        uint progressCount = 0;
        int percentage_complete = 0;

        // Write each grid point and it's optimal 'control'
        do
        {
            for (uint ii = 0; ii < mcp.getNumOfGrids(); ++ii)
            {
                gridLocation[ii] = mcp.getGridAtIndex(gridIndices.getIndexOfDim(ii), ii);
                markovControlFileStream << gridLocation[ii] << " ";
            }

            // If progress has jumped up more than 1% then redraw the progress bar
            // To reduce multiple prints of same percentage, just skip until a different percentage is reached
            auto prog = int(static_cast<float>(++progressCount*100.0/totalCountToComplete));
            if (prog != percentage_complete)
            {
                percentage_complete = prog;
                // We say its the total number of grid elements + 1 since the 1% remaining is for vectors below
                utils::printProgress(percentage_complete);
            }


            markovControlFileStream << markovCA.getMarkovControlFunction(gridLocation) << NEWL;
        } while (gridIndices.nextGridElement(mcp));

    }
    markovControlFileStream.close();

    // Plot the data using gnuplot ------------------------------------------------------ //
    // TODO: Put in plotting commands
}