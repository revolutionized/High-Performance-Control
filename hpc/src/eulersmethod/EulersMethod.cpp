//

#include "EulersMethod.h"
#include <cmath>
#include <iostream>

EulersMethod::EulersMethod(double timeLeftBound,
                           double timeRightBound,
                           double deltaTime,
                           unsigned int dimensions)
    : deltaTime_(deltaTime),
      dimensions_(dimensions)
{
    assert(timeRightBound > timeLeftBound && deltaTime > 0.0 && dimensions > 0);
    timeLength_ = static_cast<uint>(floor((timeRightBound - timeLeftBound)/deltaTime));

    setupSolution();
}

EulersMethod::EulersMethod(double timeLeftBound,
                           double timeRightBound,
                           unsigned int timeLength,
                           unsigned int dimensions)
        : timeLength_(timeLength),
          dimensions_(dimensions)
{
    assert(timeRightBound > timeLeftBound && timeLength > 0 && dimensions > 0);
    deltaTime_ = (timeRightBound - timeLeftBound) / (timeLength - 1);

    setupSolution();
}

EulersMethod::EulersMethod(GridParameters& gpa)
        : timeLength_(gpa.getGridLength(0)),
          deltaTime_(gpa.getDeltaGrid(0)),
          dimensions_(gpa.getNumOfGrids())
{
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

void EulersMethod::solve(v_fcn2dep& fcnDerivative,
                         const double* initGuess,
                         MarkovChainApproximation* mca)
{
    std::cout << "==== Starting Euler's Method Approximation ====" << std::endl;

    for (unsigned int ii = 0; ii < dimensions_; ++ii)
    {
        (*solution_)[0][ii] = initGuess[ii];
    }

    for (unsigned int ii = 1; ii < timeLength_; ++ii)
    {
        // Solve derivative
        double xDash[dimensions_];
        if (mca != nullptr)
        {
            fcnDerivative((*solution_)[ii-1], mca->getMarkovControlFunction((*solution_)[ii-1]), xDash);
        }
        else
        {
            fcnDerivative((*solution_)[ii-1], 0, xDash);
        }

        // Add derivative by deltaTime to the previous solution
        for (int jj = 0; jj < dimensions_; ++jj)
        {
            (*solution_)[ii][jj] = (*solution_)[ii-1][jj] + xDash[jj]*deltaTime_;
        }

        // Print progress to screen
        printProgress(static_cast<float>(ii) / static_cast<float>(timeLength_-1));
    }

    // Print one new line since our progress bar uses "\r" which returns to the beginning, and then
    // two new lines at the end to have a gap after euler conclusion
    std::cout << NEWL << "==== Finished Euler's Method Approximation ====" << NEWL << NEWL;
}


void EulersMethod::saveSolution(std::ofstream& stream)
{
    int counter = 0;
    for (auto& sol : *solution_)
    {
        // Print the time component
        stream << counter*deltaTime_ << " ";
        counter++;

        // Print the value that the solution has at each dimension
        for (int jj = 0; jj < dimensions_; jj++)
        {
            stream << sol[jj] << " ";
        }

        // Create a new line
        stream << NEWL;
    }
}

void EulersMethod::getSolutionAt(int index, double* out) const
{
    for (int ii = 0; ii < dimensions_; ++ii)
    {
        out[ii] = (*solution_)[index][ii];
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

void EulersMethod::printProgress(float progress) const
{
    int barWidth = 50;

    std::cout << "[";
    auto pos = static_cast<int>(barWidth * progress);
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
    std::cout << "] " << "%" << int(progress * 100.0) << "\r";
    std::cout.flush();
}
