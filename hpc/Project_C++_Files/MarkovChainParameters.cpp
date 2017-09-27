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
    : h_(h),
      alphaLeftBound_(alphaLeftBound),
      alphaLength_(alphaLength),
      numOfGrids_(dimensions)
{
    // Assert appropriate parameters have been provided
    assertParameters(gridLeftBound, gridRightBound, gridLength, alphaLeftBound, alphaRightBound, alphaLength);

    // Set up grid spacing and assign memory for own copies of lower bound and lengths
    uint numOfGridDimensions = 0;
    while (gridLeftBound[numOfGridDimensions] != nullptr)
    {
        numOfGridDimensions++;
    };
    deltaGrid_ = new double[numOfGridDimensions];
    gridLeftBound_ = new double[numOfGridDimensions];
    gridLength_ = new uint[numOfGridDimensions];
    for (int ii = 0; ii < numOfGridDimensions; ++ii)
    {
        deltaGrid_[ii] = (gridRightBound[ii] - gridLeftBound[ii]) / (gridLength[ii] - 1);
        gridLeftBound_[ii] = gridLeftBound[ii];
        gridLength_[ii] = gridLength[ii];
    }

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
    : h_(h),
      alphaLeftBound_(alphaLeftBound),
      deltaAlpha_(deltaAlpha),
      numOfGrids_(dimensions)
{
    // Assert appropriate parameters have been provided
    assertParameters(gridLeftBound, gridRightBound, deltaGrid, alphaLeftBound, alphaRightBound, deltaAlpha);

    deltaGrid_ = new double[numOfGrids_];
    gridLeftBound_ = new double[numOfGrids_];
    gridLength_ = new uint[numOfGrids_];
    for (int ii = 0; ii < numOfGrids_; ++ii)
    {
        deltaGrid_[ii] = deltaGrid[ii];
        gridLeftBound_[ii] = gridLeftBound[ii];
        gridLength_[ii] = static_cast<uint>(floor((gridRightBound[ii] - gridLeftBound[ii])/deltaGrid[ii]));
    }

    // Set up alpha discretisation
    alphaLength_ = static_cast<uint>(floor((alphaRightBound - alphaLeftBound)/deltaAlpha));
}

MarkovChainParameters::~MarkovChainParameters()
{
    delete gridLeftBound_;
    delete gridLength_;
    delete deltaGrid_;
}


// SETTERS ------------------------------------------------------------------------------------------------- SETTERS

void MarkovChainParameters::setMaxIterations(uint maxIters)
{
    maxIterations_ = maxIters;
}

void MarkovChainParameters::setRelativeError(double epsErr)
{
    epsErr_ = epsErr;
}


// GETTERS ------------------------------------------------------------------------ GETTERS //

double MarkovChainParameters::getGridAtIndex(uint index, uint gridNum)
{
    assert(gridNum < numOfGrids_ && index < gridLength_[gridNum]);

    double gridPosition = gridLeftBound_[gridNum] + index * deltaGrid_[gridNum];
    return gridPosition;
}

void MarkovChainParameters::getGridAtIndex(uint* index, double* outGrid)
{
    for (uint ii = 0; ii < numOfGrids_; ++ii)
    {
        outGrid[ii] = getGridAtIndex(index[ii], ii);
    }
}

double MarkovChainParameters::getAlphaAtIndex(uint index)
{
    assert(index < alphaLength_);

    double alphaPosition = alphaLeftBound_ + index * deltaAlpha_;
    return alphaPosition;
}

uint MarkovChainParameters::getNumOfGrids()
{
    return numOfGrids_;
}

uint MarkovChainParameters::getGridLength(uint gridNum)
{
    assert(gridNum < numOfGrids_);
    return gridLength_[gridNum];
}

uint MarkovChainParameters::getAlphaLength()
{
    return alphaLength_;
}

double MarkovChainParameters::getDeltaGrid(uint gridNum)
{
    return deltaGrid_[gridNum];
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

// PRIVATE METHODS -------------------------------------------------------- PRIVATE METHODS //

void MarkovChainParameters::assertParameters(const double* gridLeftBound,
                                             const double* gridRightBound,
                                             const uint* gridLength,
                                             double alphaLeftBound,
                                             double alphaRightBound,
                                             uint alphaLength)
{
    // First we check that both grid bound arrays and grid length are all of the same array size
    if (gridLeftBound[numOfGrids_-1] == nullptr
        || gridRightBound[numOfGrids_-1] == nullptr
        || gridLength[numOfGrids_-1] == nullptr)
    {
        assert(false);
    }

    // Now we ensure appropriate bounds have been given
    assert(alphaLeftBound < alphaRightBound);
    assert(gridLeftBound < gridRightBound);

    // Now we ensure appropriate grid lengths and alpha lengths have been given
    for (int ii = 0; ii < numOfGrids_; ++ii)
    {
        assert(gridLength[ii] > 1);
    }
    assert(alphaLength > 1);
}

void MarkovChainParameters::assertParameters(const double* gridLeftBound,
                                             const double* gridRightBound,
                                             const double* deltaGrid,
                                             double alphaLeftBound,
                                             double alphaRightBound,
                                             double deltaAlpha)
{
    // First we check that both grid bound arrays and grid length are all of the same array size
    if (gridLeftBound[numOfGrids_-1] == nullptr
        || gridRightBound[numOfGrids_-1] == nullptr
        || deltaGrid[numOfGrids_-1] == nullptr)
    {
        assert(false);
    }

    // Now we ensure appropriate bounds have been given
    assert(alphaLeftBound < alphaRightBound);
    assert(gridLeftBound < gridRightBound);

    // Now we ensure appropriate grid lengths and alpha lengths have been given
    for (int ii = 0; ii < numOfGrids_; ++ii)
    {
        assert(deltaGrid[ii] > 1);
    }
    assert(deltaAlpha > 0);
}









