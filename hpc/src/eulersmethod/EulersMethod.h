//
// Created by David on 26/09/2017.
//

#pragma once

#include <functional>
#include "multidim/GridParameters.h"
#include "mca/MarkovChainApproximation.h"

#ifdef _WIN32
// New line for Windows
#define NEWL "\r\n"
#else
// New line for Unix
#define NEWL "\n"
#endif

// Pre-declaration - needed for correct compilation
namespace utils
{
    void printProgress(int progressComplete);
}

/// Euler's method for approximating a differential system of equations. To be used with the MarkovChainApproximation
/// library to ensure the minimum control is used.
class EulersMethod
{
public:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

	/// Constructor to Euler's Method which allows you to specify the delta between time steps
	/// \param timeLeftBound The starting point of Euler's method. Can be less than 0.0 if that makes sense to the
	/// problem.
	/// \param timeRightBound The ending point
	/// \param deltaTime
	/// \param dimensions
    EulersMethod(double timeLeftBound,
                 double timeRightBound,
                 double deltaTime,
                 unsigned int dimensions);

	/// Constructor to Euler's Method which acts like MATLAB's linspace command, and you give it the length of the
	/// time vector rather than the delta between time steps.
	/// \param timeLeftBound
	/// \param timeRightBound
	/// \param timeLength
	/// \param dimensions
    EulersMethod(double timeLeftBound,
                 double timeRightBound,
                 unsigned int timeLength,
                 unsigned int dimensions);

    explicit EulersMethod(GridParameters& gpa);

    ~EulersMethod();

    void solve(v_fcn_pddpd& fcnDerivative,
               const double* initGuess,
               MarkovChainApproximation* mca);

    void saveSolution(std::ofstream& stream);

    void getSolutionAt(int index, double* out) const;

private:

    void setupSolution();

    // FIELDS --------------------------------------------------------------------------------------------------- FIELDS
    std::vector<double*>* solution_ = nullptr;
    double deltaTime_;
    unsigned int timeLength_;
    unsigned int dimensions_;
};


