/* Script example and test case for 1D problem.
 * Written by: David Dos Santos (student id - 7521626)
 *
 * This script is given as an example on how to access the Markov Chain Approximation library (the 1d version), and how
 * to access the continuous decomposition library. This script also shows how we can compare the two files.
 * TODO: Write header comment
 *
 *
 */

#include <cmath>
#include <cstring>
#include <iostream>
#include <functional>
#include <fstream>

#include "Functions.h"
#include "MarkovChainApproximation.h"
#include "EulerMethod.h"
extern "C"
{
// #include "c3.h"
// #include "c3sc.h"
}

using namespace std;

int main()
{
    // -------------------------------------------------------------------------------------------------------- EXACT //
    // First we consider the exact method ----------------------------------------------- //
    auto euler_grid_size = static_cast<unsigned int>(pow(10, 2));
    double startingBound = 0.0;
    double exitingBound = 3.0;
    double initGuess = 2.0;

    std::function<double(double)> fcnExactControl = exactMinimumControl;
    const std::function<double(double, double)> fcnDerivativeExact = [](double f, double x)
    {
        return problemOde(f, exactMinimumControl(f));
    };

    EulerMethod euler(startingBound, exitingBound, euler_grid_size);
    euler.solve(fcnDerivativeExact, initGuess);

    // Create file with exact data
    ofstream exactEulerFile;
    exactEulerFile.open("ExactEulerResult.dat", std::ofstream::out);
    if (exactEulerFile.good())
    {
        exactEulerFile << "t f\n";
        double** result = euler.getSolutionPtr();
        for (int i = 0; i < euler.getGridLength(); ++i)
        {
            // Fill control variable for plotting
            exactEulerFile << result[i][GRID] << " " << result[i][FUNC] << endl;
        }
        exactEulerFile.close();
    }
    // Also create file with control data
    ofstream exactControl;
    exactControl.open("ExactControl.dat", std::ofstream::out);
    if (exactControl.good())
    {
        exactControl << "t u\n";
        double** result = euler.getSolutionPtr();
        for (int i = 0; i < euler.getGridLength(); ++i)
        {
            // Fill control variable for plotting
            exactControl << result[i][GRID] << " " << exactMinimumControl(result[i][FUNC]) << endl;
        }
        exactControl.close();
    }

    // ------------------------------------------------------------------------------------------------------- MARKOV //
    // Now test the Markov Approximation method ----------------------------------------- //
    auto markov_grid_size = 3*static_cast<unsigned int>(pow(10, 1)) + 1;
    auto markov_alpha_discretisation_size = 2*static_cast<unsigned int>(pow(10,2)) + 1;
    double alphaStartingBound = -2.5;
    double alphaExitingBound = 0.5;
    double h = pow(10.0, -3.0);
    double markovInitGuess = 2.0;
    const std::function<double(double,double)> fcnCost = costFunction;

    MarkovChainApproximation markovCA(alphaStartingBound,
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

    const std::function<double(double, double)> fcnDerivativeMarkov = problemOde;
    euler.solve(fcnDerivativeMarkov, markovCA, initGuess);

    // Create file with exact data
    ofstream markovEulerFile;
    markovEulerFile.open("MarkovEulerResult.dat", std::ofstream::out);
    if (markovEulerFile.good())
    {
        markovEulerFile << "t f\n";
        double** result = euler.getSolutionPtr();
        for (int i = 0; i < euler.getGridLength(); ++i)
        {
            // Control data is already set in the Euler solve
            markovEulerFile << result[i][GRID] << " " << result[i][FUNC] << endl;
        }
        markovEulerFile.close();
    }
    // Also create file with markov control data
    ofstream markovControlFile;
    markovControlFile.open("MarkovControl.dat", std::ofstream::out);
    if (markovControlFile.good())
    {
        markovControlFile << "t u\n";
        double** result = euler.getSolutionPtr();
        for (int i = 0; i < euler.getGridLength(); ++i)
        {
            // Control data is already set in the Euler solve
            markovControlFile << result[i][GRID] << " " << markovCA.getMarkovControlFunction(result[i][FUNC]) << endl;
        }
        markovControlFile.close();
    }

    // ----------------------------------------------------------------------------------------------------------- C3 //
    // Now test the C3 (Tensor decomposition method) ------------------------------------ //
    /*size_t dim = 2;
    size_t dx = dim;
    size_t dw = dim;
    size_t du = 1;
    auto lb = new double(dx);
    auto ub = new double(dx);
    auto Narr = new size_t(dx);
    size_t N = 20;
    for (int i = 0; i < dx; ++i)
    {
        lb[i] = -2.0;
        ub[i] = 2.0;
        Narr[i] = N;
    }
    double beta = 0.1;
    struct C3Control* c3c = c3control_create(dx, du, dw, lb, ub, Narr, beta);
    c3control_destroy(c3c);*/

    // Plot the data using gnuplot ------------------------------------------------------ //
    system("cd .. && ./viewplots.gla");

    // Finished with no errors ---------------------------------------------------------- //
    cout << " ~~~ Finished executing ~~~ ";
    return 0;
}

