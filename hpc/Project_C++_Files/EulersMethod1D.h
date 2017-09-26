//
// Created by David Dos Santos on 16/08/17.
//
#pragma once

#include <functional>
#include "MarkovChainApproximation1D.h"

#define GRID 1
#define FUNC 2

///
class EulersMethod1D
{
public:
    ///
    /// \param leftBound
    /// \param rightBound
    /// \param gridSize
	explicit EulersMethod1D(double leftBound, double rightBound, unsigned int gridSize);

    ///
    ~EulersMethod1D();

    ///
    /// \param fcnDerivative
    /// \param initGuess
	void solve(const std::function<double(double, double)> &fcnDerivative, double initGuess);

    ///
    /// \param fcnDerivative
    /// \param mca
    /// \param initGuess
    void solve(const std::function<double(double, double)>& fcnDerivative,
			   MarkovChainApproximation1D& mca,
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