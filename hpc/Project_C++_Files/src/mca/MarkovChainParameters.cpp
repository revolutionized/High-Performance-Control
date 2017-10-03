//
// Created by david on 23/09/17.
//

#include "MarkovChainParameters.h"
#include <iostream>
#include <cmath>

// CONSTRUCTORS + DESTRUCTOR ------------------------------------ CONSTRUCTORS + DESTRUCTOR //

MarkovChainParameters::MarkovChainParameters(const double* gridLeftBound,
                                             const double* gridRightBound,
                                             uint* gridLength,
                                             double alphaLeftBound,
                                             double alphaRightBound,
                                             uint alphaLength,
                                             double h,
                                             uint dimensions)
    : GridParameters(gridLeftBound, gridRightBound, gridLength, dimensions),
      h_(h),
      alphaLeftBound_(alphaLeftBound),
      alphaLength_(alphaLength)
{
    // Assert appropriate parameters have been provided
    assert(alphaLeftBound < alphaRightBound);
    assert(alphaLength > 1);

    // Set up alpha discretisation
    deltaAlpha_ = (alphaRightBound - alphaLeftBound) / (alphaLength - 1);
}

MarkovChainParameters::MarkovChainParameters(const double* gridLeftBound,
                                             const double* gridRightBound,
                                             const double* deltaGrid,
                                             double alphaLeftBound,
                                             double alphaRightBound,
                                             double deltaAlpha,
                                             double h,
                                             uint dimensions)
    : GridParameters(gridLeftBound, gridRightBound, deltaGrid, dimensions),
      h_(h),
      alphaLeftBound_(alphaLeftBound),
      deltaAlpha_(deltaAlpha)
{
    // Assert appropriate parameters have been provided
    assert(alphaLeftBound < alphaRightBound);
    assert(deltaAlpha > 0.0);

    // Set up alpha discretisation
    alphaLength_ = static_cast<uint>(floor((alphaRightBound - alphaLeftBound)/deltaAlpha));
}


MarkovChainParameters::MarkovChainParameters(const MarkovChainParameters& mcp)
        : GridParameters(static_cast<GridParameters>(mcp)),
          h_(mcp.h_),
          alphaLeftBound_(mcp.alphaLeftBound_),
          alphaLength_(mcp.alphaLength_),
          deltaAlpha_(mcp.deltaAlpha_)
{}

// SETTERS ------------------------------------------------------------------------------------------------- SETTERS

void MarkovChainParameters::setMaxIterations(uint maxIters)
{
    maxIterations_ = maxIters;
}


void MarkovChainParameters::setMinError(double epsErr)
{
    epsErr_ = epsErr;
}

// GETTERS ------------------------------------------------------------------------ GETTERS //

double MarkovChainParameters::getAlphaAtIndex(uint index)
{
    assert(index < alphaLength_);

    double alphaPosition = alphaLeftBound_ + index * deltaAlpha_;
    return alphaPosition;
}

uint MarkovChainParameters::getAlphaLength()
{
    return alphaLength_;
}

uint MarkovChainParameters::getMaxIterations()
{
    return maxIterations_;
}

double MarkovChainParameters::getRelativeError()
{
    return epsErr_;
}

double MarkovChainParameters::getH()
{
    return h_;
}

unsigned int MarkovChainParameters::getGridLengthAccumulation()
{
    uint total = 1;
    for (int ii = 0; ii < getNumOfGrids(); ++ii)
    {
        total *= getGridLength(ii);
    }
    return total;
}

unsigned int MarkovChainParameters::getThreadCount() const
{
    return threadCount_;
}

void MarkovChainParameters::setThreadCount(unsigned int threadCount)
{
    MarkovChainParameters::threadCount_ = threadCount;
}









