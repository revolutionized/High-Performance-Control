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

// Pre-declaration
namespace utils
{
    void printProgress(int progressComplete);
}

class EulersMethod
{
public:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

    EulersMethod(double timeLeftBound,
                 double timeRightBound,
                 double deltaTime,
                 unsigned int dimensions);

    EulersMethod(double timeLeftBound,
                 double timeRightBound,
                 unsigned int timeLength,
                 unsigned int dimensions);

    explicit EulersMethod(GridParameters& gpa);

    ~EulersMethod();

    void solve(v_fcn2dep& fcnDerivative,
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


