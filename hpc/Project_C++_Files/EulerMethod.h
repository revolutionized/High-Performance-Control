//
// Created by David Dos Santos on 16/08/17.
//
#pragma once

#include <functional>
#include "MarkovChainApproximation.h"

#define GRID 1
#define FUNC 2

///
class EulerMethod
{
public:
    ///
    /// \param leftBound
    /// \param rightBound
    /// \param gridSize
	explicit EulerMethod(double leftBound, double rightBound, unsigned int gridSize);

    ///
    ~EulerMethod();

    ///
    /// \param fcnDerivative
    /// \param initGuess
	void solve(const std::function<double(double, double)> &fcnDerivative, double initGuess);

    ///
    /// \param fcnDerivative
    /// \param mca
    /// \param initGuess
    void solve(const std::function<double(double, double)>& fcnDerivative,
			   MarkovChainApproximation& mca,
			   double initGuess);

    ///
    /// \return
    double** getSolutionPtr();

    ///
    /// \return
	unsigned int getGridLength();

    ///
	void clearSolution();
private:
    ///
    double getGridAtIndex(unsigned int index);
    ///
    void setUpSolutionDuringSolve(double initGuess);
    ///
    double leftBound_;
    ///
	double rightBound_;
    ///
	double dx_;
    ///
	unsigned int gridLength_;
    ///
    double** solution_ = nullptr;
};