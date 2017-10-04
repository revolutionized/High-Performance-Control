/* Script example and test case for 1D problem.
 * Written by: David Dos Santos (student id - 7521626)
 *
 * This script is given as an example on how to access the Markov Chain Approximation library (the 1d version), and how
 * to access the continuous decomposition library. This script also shows how we can compare the two files.
 */

#include <cmath>
#include <cstring>
#include <iostream>
#include <functional>
#include <fstream>
#include <cassert>

#include "Functions.h"
#include "../../src/mca_1Donly/MarkovChainApproximation1D.h"
#include "../../src/eulersmethod_1Donly/EulersMethod1D.h"
extern "C"
{
#include "c3/c3.h" // Already included in the c3sc file
#include "c3sc/c3sc.h"
#include "c3sc/bellman.h"
#include "c3sc/valuefunc.h"
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

    // Here we create a std::function for the ODE derivative. This is made from the problem ODE given in "functions.h"
    // and matched with the exact minimum control value needed for optimal control (that is the analytical exact value)
    std::function<double(double)> fcnExactControl = exactMinimumControl;
    const std::function<double(double, double)> fcnDerivativeExact = [](double f, double x)
    {
        return problemOde(f, exactMinimumControl(f));
    };

    // We use the Euler Method class (giving it it's bounding grid)
    EulersMethod1D euler(startingBound, exitingBound, euler_grid_size);
    // And then solve the ODE using Euler's method and the function for the derivative we put together before
    euler.solve(fcnDerivativeExact, initGuess);

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
    auto markov_alpha_discretisation_size = 2*static_cast<unsigned int>(pow(10,2)) + 1;
    double alphaStartingBound = -2.5;
    double alphaExitingBound = 0.5;
    double h = pow(10.0, -3.0);
    double markovInitGuess = 2.0;
    const std::function<double(double,double)> fcnCost = costFunction;

    // Markov Chain Approximation (MCA) will solve the optimal cost values and ODE values at each point
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

    // The EulersMethod1D function has a solve that allows you to pass it the MCA object, and thus it utilises the MCA
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

    // ----------------------------------------------------------------------------------------------------------- C3 //
    // Now test the C3 (Tensor decomposition method) ------------------------------------ //

    // SO THE C3 CODE DOES EXECUTE BUT HAVEN'T QUITE FIGURED OUT HOW TO USE THE PRODUCED RESULTS FOR PLOTTING OR
    // CONFIRMATION ETC.

    // First we set all the "scene" parameters
    size_t dim = 2;
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
    // Now with those parameters we "set the scene"
    struct C3Control* c3c = c3control_create(dx, du, dw, lb, ub, Narr, beta);

    // We have created the C3Control scene, but now we need to add a diffusion equation (our ODE) to it. In this 1-d
    // example there is no diffusion so it's pointless to add it
    c3control_add_diff(c3c, diff, nullptr);

    // We also have to add drift now in the same way
    c3control_add_drift(c3c, drift, nullptr);

    // Now add the cost function to the scene. Other costs can be added, such as boundary costs and obstacle costs.
    c3control_add_stagecost(c3c, stagecost);

    // Add the cost function and boundary function
    c3control_add_stagecost(c3c,stagecost);
    c3control_add_boundcost(c3c,boundcost);

    // ??? Copied from the Dubins_car_new example, so don't fully understand it all ???
    char filename[256];
    sprintf(filename, "%s.c3", "cost");
    double ** xgrid = c3control_get_xgrid(c3c);
    struct ValueF * cost = valuef_load(filename, Narr, xgrid);
    struct ApproxArgs * aargs = approx_args_init();
    if (cost == nullptr){
        cost = c3control_init_value(c3c, startcost, nullptr, aargs, 0);
    }

    size_t maxiter_vi = 10;
    double abs_conv_vi = 1e-3;
    size_t maxiter_pi = 10;
    double abs_conv_pi = 1e-2;
    struct Diag * diag = NULL;
    char filename_diag[256];
    struct c3Opt * opt = c3opt_alloc(BRUTEFORCE,du);
    int verbose = 0;
    sprintf(filename_diag,"n%zu_%s", N,"diagnostic.dat");

    cout << "==== Starting C3 Iteration technique ====" << endl;
    for (size_t ii = 0; ii < maxiter_vi; ii++){
        /* for (size_t ii = 0; ii < 0; ii++){ */
        struct ValueF * next = c3control_pi_solve(c3c,maxiter_pi,abs_conv_pi,
                                                  cost,aargs,opt,verbose,&diag);

        valuef_destroy(cost); cost = NULL;
        cost = c3control_vi_solve(c3c,1,abs_conv_vi,next,aargs,opt,verbose,&diag);
        valuef_destroy(next); next = NULL;

        sprintf(filename, "%s.c3", "cost");
        int saved = valuef_save(cost,filename);
        assert (saved == 0);

        if (verbose != 0){
            printf("ii=%zu ranks=",ii);
            size_t * ranks = valuef_get_ranks(cost);
            iprint_sz(4,ranks);
        }

        saved = diag_save(diag,filename_diag);
        assert (saved == 0);
    }

    // Finally, we destroy our C3Control data
    valuef_destroy(cost); cost = NULL;
    c3control_destroy(c3c); c3c = NULL;
    diag_destroy(&diag); diag = NULL;
    c3opt_free(opt); opt = NULL;
    approx_args_free(aargs); aargs = NULL;

    // Plot the data using gnuplot ------------------------------------------------------ //
    // The following command just requests to run a shell script with the commands in the file ../viewplots.gla
//    system("cd .. && ./viewplots.gla");

    // Finished with no errors ---------------------------------------------------------- //
    cout << " ~~~ Finished executing ~~~ ";
    return 0;
}

