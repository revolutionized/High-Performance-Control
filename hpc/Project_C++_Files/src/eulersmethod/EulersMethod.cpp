//

#include "EulersMethod.h"
#include <cmath>
#include <iostream>

EulersMethod::EulersMethod(double timeLeftBound,
                           double timeRightBound,
                           double deltaTime,
                           unsigned int dimensions)
    : timeLeftBound_(timeLeftBound),
      deltaTime_(deltaTime),
      dimensions_(dimensions)
{
    assert(timeRightBound > timeLeftBound && timeLength > 0 && dimensions > 0);
    timeLength_ = static_cast<uint>(floor((timeRightBound - timeLeftBound)/deltaTime));

    setupSolution();
}

EulersMethod::EulersMethod(double timeLeftBound,
                           double timeRightBound,
                           unsigned int timeLength,
                           unsigned int dimensions)
        : timeLeftBound_(timeLeftBound),
          timeLength_(timeLength),
          dimensions_(dimensions)
{
    assert(timeRightBound > timeLeftBound && timeLength > 0 && dimensions > 0);
    deltaTime_ = (timeRightBound - timeLeftBound) / (timeLength - 1);

    setupSolution();
}

EulersMethod::~EulersMethod()
{
    for (auto& it : *solution_)
    {
        delete it;
    }
    delete solution_;
}

template <typename double>
void EulersMethod::solve(fcn2dep& fcnDerivative,
                         double* initGuess,
                         MarkovChainApproximation* mca)
{
    for (unsigned int ii = 0; ii < dimensions_; ++ii)
    {
        (*solution_)[0][ii] = initGuess[ii];
    }

    for (unsigned int ii = 1; ii < timeLength_; ++ii)
    {
        // Solve derivative
        double* xDash;
        if (mca != nullptr)
        {
            xDash = fcnDerivative((*solution_)[ii-1], mca->getMarkovControlFunction((*solution_)[ii-1]));
        }
        else
        {
            xDash = fcnDerivative((*solution_)[ii-1], 0);
        }

        // Add derivative by deltaTime to the previous solution
        for (int jj = 0; jj < dimensions_; ++jj)
        {
            (*solution_)[ii][jj] = (*solution_)[ii-1][jj] + xDash[jj]*deltaTime_;
        }

        // Since each of the functionDerivatives return a pointer we need to deallocate that pointer here
        delete xDash;

        // Print progress to screen
        printProgress(static_cast<float>(ii/(timeLength_-1)));
    }
}

void EulersMethod::saveSolution(std::ofstream& stream)
{
    for (auto& sol : *solution_)
    {
        for (int jj = 0; jj < dimensions_; jj++)
        {
            stream << sol[jj] << " ";
        }
        stream << NEWL;
    }
}


// PRIVATE METHODS -------------------------------------------------------- PRIVATE METHODS //

void EulersMethod::setupSolution()
{
    // Setup solution memory and insert initial guess
    solution_ = new std::vector<double*>(timeLength_);
    for (auto& ii : *solution_)
    {
        ii = new double[dimensions_];
    }
}

void EulersMethod::printProgress(float progress)
{
    int barWidth = 70;

    std::cout << "[";
    int pos = static_cast<int>(barWidth * progress);
    for (int ii = 0; ii < barWidth; ++ii) {
        if (ii < pos)
        {
            std::cout << "=";
        }
        else if (ii == pos)
        {
            std::cout << ">";
        }
        else
        {
            std::cout << " ";
        }
    }
    std::cout << "] " << int(progress * 100.0) << NEWL;
    std::cout.flush();
}
