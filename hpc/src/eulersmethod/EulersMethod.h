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
/// library to ensure the minimum control is used. Used in this project to calculate trajectory.
class EulersMethod
{
public:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

	/// Constructor to Euler's Method which allows you to specify the delta between time steps
	/// \param timeLeftBound The starting point of Euler's method. Can be less than 0.0 if that makes sense to the
	/// problem.
	/// \param timeRightBound The ending point. Must be greater than the left bound
	/// \param deltaTime Distance between time steps (for Euler's method, the smaller this is the greater accuracy
	/// you will have)
	/// \param dimensions The output will of Euler's will be each dimension against time
    EulersMethod(double timeLeftBound,
                 double timeRightBound,
                 double deltaTime,
                 unsigned int dimensions);

	/// Constructor to Euler's Method which acts like MATLAB's linspace command, and you give it the length of the
	/// time vector rather than the delta between time steps.
	/// \param timeLeftBound The starting point of Euler's method. Can be less than 0.0 if that makes sense to the
	/// problem.
	/// \param timeRightBound The ending point. Must be greater than the left bound
	/// \param timeLength Total number of elements in the time vector (for Euler's method, the larger this is the
	/// greater accuracy you will have - since large timeLength equates to smaller steps between time)
	/// \param dimensions The output will of Euler's will be each dimension against time
    EulersMethod(double timeLeftBound,
                 double timeRightBound,
                 unsigned int timeLength,
                 unsigned int dimensions);

	/// This is the preferred method for constructing Euler's method.
	/// \param gpa Contains all the parameters needed for a successful Euler's execution
    explicit EulersMethod(GridParameters& gpa);

    ~EulersMethod();

	/// This will solve the given derivative using the "minimum/optimised" control found using the
	/// MarkovChainApproximation.
	/// \param fcnDerivative A std::function that represents the Differential system trying to be solved. It returns
	/// void, and takes as its values:
	/// Parameter (in order)[type] | Description
	/// -------------------------- | -----------
	/// x [double*]				   | Current grid value / location. Eg. x[0]=x, x[1]=y, and so on
	/// alpha [double]			   | Control value
	/// out [double*]			   | This will be filled with the result (and hence should be the same size as x)
	/// \param initGuess The initial guess of each component in the system (hence should be the length of the number
	/// of dimensions)
	/// \param mca A pointer to the MarkovChainApproximation used to minimise the control value of the differential
	/// system.
    void solve(v_fcn_pddpd& fcnDerivative,
               const double* initGuess,
               MarkovChainApproximation* mca);

	/// Will save the solution (along with the time it occurs)
	/// \param stream The stream to output the solution to
    void saveSolution(std::ofstream& stream);

	/// Returns the solution at a specific point, note it uses an index instead of "closest point to this value"
	/// \param index The index should be a valid range, hence it should be > 0 and < the length of the time vector
	/// \param out The result will be stored in here (and hence should be the same dimension length as the system
	/// originally solved in the solve function.
    void getSolutionAt(int index, double* out) const;

private:
	// METHODS ------------------------------------------------------------------------------------------------- METHODS

	/// Tidy-up function that sets up the memory needed for the solution space
    void setupSolution();

    // FIELDS --------------------------------------------------------------------------------------------------- FIELDS

	/// The solution is a vector for each time step, and a double for each dimension
    std::vector<double*>* solution_ = nullptr;
	/// The difference between time steps
    double deltaTime_;
	/// The length of the time vector
    unsigned int timeLength_;
	/// The number of dimensions the system has
    unsigned int dimensions_;
};


