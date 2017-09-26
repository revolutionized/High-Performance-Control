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
#include <cassert>

#include "Functions2.h"
#include "MarkovChainApproximation.h"
#include "MarkovChainParameters.h"
#include "EulerParameters.h"
#include "EulersMethod.h"

using std::cout;
using std::endl;
using std::ofstream;

int main()
{
    // -------------------------------------------------------------------------------------------------------- EXACT //
    // First we consider the exact method ----------------------------------------------- //
    auto euler_grid_size = static_cast<unsigned int>(pow(10, 2));
    double startingBound = 0.0;
    double exitingBound = 3.0;
    double initGuess = 2.0;

    // Here we create a std::function for the ODE derivative. This is made from the problem ODE given in "functions.h"
    // and matched with the exact minimum control value needed for optimal control (that is the analytical exact value)
    std::function<double(double*)> fcnExactControl = exactMinimumControl2;
    const std::function<double(double*)> fcnDerivativeExact = [](double* x)
    {
        return problemOde2(x, exactMinimumControl2(x));
    };

    const double gridLeftBound[1] = {startingBound};
    const double gridRightBound[1] = {exitingBound};
    unsigned int gridLength[1] = {euler_grid_size};
    unsigned int dimensions = 1;
    EulerParameters epm(gridLeftBound, gridRightBound, gridLength, dimensions);

    // We use the EulersMethod class
    EulersMethod euler(epm);
    // And then solveExact the ODE using Euler's method and the function for the derivative we put together before
    euler.solveExact(fcnDerivativeExact, initGuess);

    // Create file and load it with the exact data results (the state space)
    ofstream exactEulerFile;
    exactEulerFile.open("ExactEulerResult.dat", std::ofstream::out);
    if (exactEulerFile.good())
    {
        exactEulerFile << "t f\n";
        // Euler mantains data in a vector, but we can call the getSolutionPtr to request it as a double** (c-style 2-di
        // mensional array.
        double** result = euler.getSolutionPtr();
        for (int i = 0; i < euler.getGridLength(); ++i)
        {
            // File is filled with the x-grid and y-grid (creating a 2D plot)
            exactEulerFile << result[i][GRID] << " " << result[i][FUNC] << endl;
        }
        exactEulerFile.close();
    }
    // Also create file with the optimal control
    ofstream exactControl;
    exactControl.open("ExactControl.dat", std::ofstream::out);
    if (exactControl.good())
    {
        exactControl << "t u\n";
        double** result = euler.getSolutionPtr();
        for (int i = 0; i < euler.getGridLength(); ++i)
        {
            // File is filled with the x-grid and y-grid (creating a 2D plot)
            exactControl << result[i][GRID] << " " << exactMinimumControl(result[i][FUNC]) << endl;
        }
        exactControl.close();
    }

    // ------------------------------------------------------------------------------------------------------- MARKOV //
    // Now test the Markov Approximation method ----------------------------------------- //
    auto markov_grid_size = 3*static_cast<unsigned int>(pow(10, 1)) + 1;
    const double gridLeftBound[1] = {0.0};
    const double gridRightBound[1] = {3.0};

    auto markov_alpha_discretisation_size = 2*static_cast<unsigned int>(pow(10,2)) + 1;
    double alphaStartingBound = -2.5;
    double alphaExitingBound = 0.5;
    double h = pow(10.0, -3.0);
    MarkovChainParameters mca(gridLeftBound,
                              gridRightBound,
                              markov_grid_size,
                              alphaStartingBound,
                              alphaExitingBound,
                              markov_alpha_discretisation_size,
                              h);

    double markovInitGuess = 2.0;
    const std::function<double(double,double)> fcnCost = costFunction;

    // Markov Chain Approximation (MCA) will solveExact the optimal cost values and ODE values at each point
    MarkovChainApproximation1D markovCA(alphaStartingBound,
                                      alphaExitingBound,
                                      markov_alpha_discretisation_size,
                                      startingBound,
                                      exitingBound,
                                      markov_grid_size,
                                      h,
                                      markovInitGuess);

    const std::function<double(double)> fcnSigma = sigmaFunction;
    const std::function<double(double, double)> fcnOde = problemOde;
    markovCA.computeMarkovApproximation(fcnCost, fcnSigma);

    // The EulersMethod1D function has a solveExact that allows you to pass it the MCA object, and thus it utilises the MCA
    // ODE state space results and optimal control results.
    const std::function<double(double, double)> fcnDerivativeMarkov = problemOde;
    euler.solve(fcnDerivativeMarkov, markovCA, initGuess);

    // Create file with exact data (same as before)
    ofstream markovEulerFile;
    markovEulerFile.open("MarkovEulerResult.dat", std::ofstream::out);
    if (markovEulerFile.good())
    {
        markovEulerFile << "t f\n";
        double** result = euler.getSolutionPtr();
        for (int i = 0; i < euler.getGridLength(); ++i)
        {
            markovEulerFile << result[i][GRID] << " " << result[i][FUNC] << endl;
        }
        markovEulerFile.close();
    }
    // Also create file with MCA optimal control data
    ofstream markovControlFile;
    markovControlFile.open("MarkovControl.dat", std::ofstream::out);
    if (markovControlFile.good())
    {
        markovControlFile << "t u\n";
        double** result = euler.getSolutionPtr();
        for (int i = 0; i < euler.getGridLength(); ++i)
        {
            markovControlFile << result[i][GRID] << " " << markovCA.getMarkovControlFunction(result[i][FUNC]) << endl;
        }
        markovControlFile.close();
    }


    // Plot the data using gnuplot ------------------------------------------------------ //
    // The following command just requests to run a shell script with the commands in the file ../viewplots.gla
//    system("cd .. && ./viewplots.gla");

    // Finished with no errors ---------------------------------------------------------- //
    cout << " ~~~ Finished executing ~~~ ";
    return 0;
}

