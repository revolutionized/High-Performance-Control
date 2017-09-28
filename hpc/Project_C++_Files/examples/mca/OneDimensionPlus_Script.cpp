/* Script example and test case for 1D+ problem
 * Written by: David Dos Santos (student id - 7521626)
 *
 * This script is given as an example on how to access the Markov Chain Approximation library (the n-dimensional
 * version). It also shows how to access the Euler library for n-dimensional problems.
 */

#include <cmath>
#include <cstring>
#include <iostream>
#include <functional>
#include <fstream>

#include "Functions2.h"
#include "mca/MarkovChainApproximation.h"
#include "multidim/GridParameters.h"
#include "eulersmethod/EulersMethod.h"

using std::cout;
using std::endl;
using std::ofstream;

int main()
{
    // -------------------------------------------------------------------------------------------------------- EXACT //
    // First we consider the exact method ----------------------------------------------- //

    // We set up the "grid" or space we will be considering
    auto euler_grid_size = static_cast<unsigned int>(pow(10, 2));
    double startingBound = 0.0;
    double exitingBound = 3.0;
    double initGuess = 2.0;
    // Only 1 dimension for the example
    const double gridLeftBound[1] = {startingBound};
    const double gridRightBound[1] = {exitingBound};
    unsigned int gridLength[1] = {euler_grid_size};
    unsigned int dimensions = 1;
    GridParameters epm(gridLeftBound, gridRightBound, gridLength, dimensions);

    // Here we create a std::function for the ODE derivative. This is made from the problem ODE given in "functions.h"
    // and matched with the exact minimum control value needed for optimal control (that is the analytical exact value)
    std::function<double(double*)> fcnExactControl = exactMinimumControl2;
    const std::function<double(double*)> fcnDerivativeExact = [](double* x)
    {
        return problemOde2(x, exactMinimumControl2(x));
    };

    // We use the EulersMethod class
    auto euler = new EulersMethod(epm);
    // And then solve the ODE using Euler's method and the function for the derivative we put together before
    euler->solve(fcnDerivativeExact, initGuess);

    // Create file and load it with the exact data results (the state space)
    ofstream exactEulerFileStream;
    exactEulerFileStream.precision(4);
    std::scientific;
    exactEulerFileStream.open("ExactEulerResult.dat", std::ofstream::out);
    if (exactEulerFileStream.good())
    {
        exactEulerFileStream << "t f" << NEWL;
        euler->saveSolution(exactEulerFileStream);
    }
    exactEulerFileStream.close();

    // Also create file with the optimal control
    ofstream exactControlFileStream;
    exactControlFileStream.open("ExactControl.dat", std::ofstream::out);
    exactControlFileStream.precision(4);
    std::scientific;
    if (exactControlFileStream.good())
    {
        exactControlFileStream << "t u" << NEWL;

        // Need to add the minimum control for each point of grid
        GridIndex gridIndices(epm.getNumOfGrids());
        gridIndices.resetToOrigin();
        double gridLocation[epm.getNumOfGrids()];

        // Fill each grid point
        do
        {
            for (uint ii = 0; ii < epm.getNumOfGrids(); ++ii)
            {
                gridLocation[ii] = epm.getGridAtIndex(gridIndices.getIndexOfDim(ii), ii);
                exactControlFileStream << gridLocation[ii] << " ";
            }

            exactControlFileStream << fcnExactControl(gridLocation) << NEWL;
        } while (gridIndices.nextGridElement(epm));

    }
    exactControlFileStream.close();

    // We delete the memory since the solution may take up a lot of space and we want to reuse Euler's method
    delete euler;

    // ------------------------------------------------------------------------------------------------------- MARKOV //
    // Now test the Markov Approximation method ----------------------------------------- //

    // We again set up the "grid" or space we will be considering
    unsigned int markov_grid_size[1];
    markov_grid_size[0] = {3*static_cast<unsigned int>(pow(10, 1)) + 1};
    auto markov_alpha_discretisation_size = 2*static_cast<unsigned int>(pow(10,2)) + 1;
    double alphaStartingBound = -2.5;
    double alphaExitingBound = 0.5;
    double h = pow(10.0, -3.0);
    MarkovChainParameters mcp(gridLeftBound,
                              gridRightBound,
                              markov_grid_size,
                              alphaStartingBound,
                              alphaExitingBound,
                              markov_alpha_discretisation_size,
                              h,
                              1);

    double markovInitGuess = 2.0;

    // Markov Chain Approximation (MCA) will solve the optimal cost values and ODE values at each point
    MarkovChainApproximation markovCA(mcp, initGuess, 4, true);
    markovCA.computeMarkovApproximation(costFunction2, driftFunction2, diffFunction2);

    // We use the EulersMethod class again
    euler = new EulersMethod(epm);

    // The EulersMethod1D function has a solve that allows you to pass it the MCA object, and thus it utilises the MCA
    // ODE state space results and optimal control results.
    const std::function<double(double*, double)> fcnDerivativeMarkov = problemOde2;
    euler->solve(fcnDerivativeMarkov, markovCA, markovInitGuess);

    // Create file with exact data (same as before)
    ofstream markovEulerFile;
    markovEulerFile.open("MarkovEulerResult.dat", std::ofstream::out);
    if (markovEulerFile.good())
    {
        exactEulerFileStream << "t f" << NEWL;
        euler->saveSolution(markovEulerFile);
    }
    exactEulerFileStream.close();

    // Also create file with MCA optimal control data
    ofstream markovControlFileStream;
    markovControlFileStream.open("MarkovControl.dat", std::ofstream::out);
    if (markovControlFileStream.good())
    {
        markovControlFileStream << "t u" << NEWL;

        // Need to add the minimum control for each point of grid
        GridIndex gridIndices(epm.getNumOfGrids());
        gridIndices.resetToOrigin();
        double gridLocation[epm.getNumOfGrids()];

        // Fill each grid point
        do
        {
            for (uint ii = 0; ii < epm.getNumOfGrids(); ++ii)
            {
                gridLocation[ii] = epm.getGridAtIndex(gridIndices.getIndexOfDim(ii), ii);
                markovControlFileStream << gridLocation[ii] << " ";
            }

            markovControlFileStream << markovCA.getMarkovControlFunction(gridLocation) << NEWL;
        } while (gridIndices.nextGridElement(epm));

    }
    markovControlFileStream.close();

    // And clear the memory
    delete euler;

    // Plot the data using gnuplot ------------------------------------------------------ //
    // The following command just requests to run a shell script with the commands in the file ../viewplots.gla
//    system("cd .. && ./viewplots.gla");

    // Finished with no errors ---------------------------------------------------------- //
    cout << " ~~~ Finished executing ~~~ ";
    return 0;
}

