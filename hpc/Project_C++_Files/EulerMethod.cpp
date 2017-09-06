//
// Created by David Dos Santos on 16/08/17.
//

#include <iostream>
#include <cassert>

#include "MarkovChainApproximation.h"
#include "EulerMethod.h"


EulerMethod::EulerMethod(double leftBound, double rightBound, unsigned int gridSize)
: leftBound_(leftBound),
  rightBound_(rightBound),
  gridLength_(gridSize)
{
	assert(leftBound_ < rightBound_);

	// Set up grid spacing
	dx_ = (rightBound_ - leftBound_) / (gridLength_ - 1);
}

EulerMethod::~EulerMethod() 
{
    clearSolution();
}

void EulerMethod::solve(const std::function<double(double, double)> &fcnDerivative, double initGuess)
{
    setUpSolutionDuringSolve(initGuess);

    // Perform Euler's method using control function derivative only
    double df;
    for (unsigned int i = 1; i < gridLength_; ++i)
    {
        df = fcnDerivative(solution_[i-1][FUNC], getGridAtIndex(i));
        solution_[i][FUNC] = solution_[i-1][FUNC] + dx_*df;
    }
}

double EulerMethod::getGridAtIndex(unsigned int index)
{
    assert(index < gridLength_);

    double gridPosition = leftBound_ + index*dx_;
    return gridPosition;
}

unsigned int EulerMethod::getGridLength()
{
    return gridLength_;
}

double** EulerMethod::getSolutionPtr()
{
    return solution_;
}

void EulerMethod::solve(const std::function<double(double, double)>& fcnDerivative,
                        MarkovChainApproximation& mca,
                        double initGuess)
{
    setUpSolutionDuringSolve(initGuess);

    // Perform Euler's method using markov control value
    double df;
    for (unsigned int i = 1; i < gridLength_; ++i)
    {
        df = fcnDerivative(solution_[i-1][FUNC], mca.getMarkovControlFunction(solution_[i-1][FUNC]));
        solution_[i][FUNC] = solution_[i-1][FUNC] + dx_*df;
    }
}

void EulerMethod::setUpSolutionDuringSolve(double initGuess)
{
    // Setup solution memory and insert initial guess
    if (solution_ == nullptr)
    {
        solution_ = new double*[gridLength_];
        for (int i = 0; i < gridLength_; ++i)
        {
            solution_[i] = new double[FUNC];
        }
    }
    solution_[0][FUNC] = initGuess;

    // Fill the grid column of the solution_ with the grid
    for (unsigned int i = 0; i < gridLength_; ++i)
    {
        solution_[i][GRID] = getGridAtIndex(i);
    }
}

void EulerMethod::clearSolution()
{
    if (solution_ != nullptr)
    {
        for (int i = 0; i < gridLength_; ++i)
        {
            delete[] solution_[i];
        }
        delete[] solution_;
    }
    solution_ = nullptr;
}

