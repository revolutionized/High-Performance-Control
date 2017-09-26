//
// Created by David Dos Santos on 16/08/17.
//

#include <iostream>
#include <cassert>

#include "MarkovChainApproximation1D.h"
#include "EulersMethod1D.h"


EulersMethod1D::EulersMethod1D(double leftBound, double rightBound, unsigned int gridSize)
: leftBound_(leftBound),
  rightBound_(rightBound),
  gridLength_(gridSize)
{
	assert(leftBound_ < rightBound_);

	// Set up grid spacing
	dx_ = (rightBound_ - leftBound_) / (gridLength_ - 1);
}

EulersMethod1D::~EulersMethod1D()
{
    clearSolution();
}

void EulersMethod1D::solve(const std::function<double(double, double)> &fcnDerivative, double initGuess)
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

double EulersMethod1D::getGridAtIndex(unsigned int index)
{
    assert(index < gridLength_);

    double gridPosition = leftBound_ + index*dx_;
    return gridPosition;
}

unsigned int EulersMethod1D::getGridLength()
{
    return gridLength_;
}

double** EulersMethod1D::getSolutionPtr()
{
    return solution_;
}

void EulersMethod1D::solve(const std::function<double(double, double)>& fcnDerivative,
                        MarkovChainApproximation1D& mca,
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

void EulersMethod1D::setUpSolutionDuringSolve(double initGuess)
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

void EulersMethod1D::clearSolution()
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

