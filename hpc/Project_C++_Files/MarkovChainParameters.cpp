//
// Created by david on 23/09/17.
//

#include "MarkovChainParameters.h"
#include <iostream>
#include <typeinfo>
#include <cmath>

// CONSTRUCTORS + DESTRUCTOR ------------------------------------ CONSTRUCTORS + DESTRUCTOR //

MarkovChainParameters::MarkovChainParameters(const double* gridLeftBound,
                                             const double* gridRightBound,
                                             unsigned int* gridLength,
                                             double alphaLeftBound,
                                             double alphaRightBound,
                                             unsigned int alphaLength,
                                             double h)
    : h_(h),
      alphaLeftBound_(alphaLeftBound),
      alphaLength_(alphaLength)
{
    // Assert appropriate parameters have been provided
    assertParameters(gridLeftBound, gridRightBound, gridLength, alphaLeftBound, alphaRightBound, alphaLength);

    // Set up grid spacing and assign memory for own copies of lower bound and lengths
    unsigned int numOfGridDimensions = 0;
    while (gridLeftBound[numOfGridDimensions] != nullptr)
    {
        numOfGridDimensions++;
    };
    deltaGrid_ = new double[numOfGridDimensions];
    gridLeftBound_ = new double[numOfGridDimensions];
    gridLength_ = new unsigned int[numOfGridDimensions];
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
                                             double h)
    : h_(h),
      alphaLeftBound_(alphaLeftBound),
      deltaAlpha_(deltaAlpha)
{
    // Assert appropriate parameters have been provided
    assertParameters(gridLeftBound, gridRightBound, deltaGrid, alphaLeftBound, alphaRightBound, deltaAlpha);

    // Set up grid spacing and assign memory for own copies of lower bound and lengths
    unsigned int numOfGridDimensions = 0;
    while (gridLeftBound[numOfGridDimensions] != nullptr)
    {
        numOfGridDimensions++;
    };
    deltaGrid_ = new double[numOfGridDimensions];
    gridLeftBound_ = new double[numOfGridDimensions];
    gridLength_ = new unsigned int[numOfGridDimensions];
    for (int ii = 0; ii < numOfGridDimensions; ++ii)
    {
        deltaGrid_[ii] = deltaGrid[ii];
        gridLeftBound_[ii] = gridLeftBound[ii];
        gridLength_[ii] = static_cast<unsigned int>(floor((gridRightBound[ii] - gridLeftBound[ii])/deltaGrid[ii]));
    }

    // Set up alpha discretisation
    alphaLength_ = static_cast<unsigned int>(floor((alphaRightBound - alphaLeftBound)/deltaAlpha));
}

MarkovChainParameters::~MarkovChainParameters()
{
    delete gridLeftBound_;
    delete gridLength_;
    delete deltaGrid_;
    delete alphaLeftBound_;
    delete alphaLength_;
    delete deltaAlpha_;
}


// SETTERS ------------------------------------------------------------------------------------------------- SETTERS

void MarkovChainParameters::setMaxIterations(unsigned int maxIters)
{
    maxIterations_ = maxIters;
}

void MarkovChainParameters::setRelativeError(double epsErr)
{
    epsErr_ = epsErr;
}


// GETTERS ------------------------------------------------------------------------ GETTERS //

double MarkovChainParameters::getGridAtIndex(unsigned int index, unsigned int gridNum)
{
    assert(gridNum < numOfGrids_ && index < gridLength_[gridNum]);

    double gridPosition = gridLeftBound_[gridNum] + index * deltaGrid_[gridNum];
    return gridPosition;
}

double* MarkovChainParameters::getGridAtIndex(unsigned int* index)
{
    double grid[numOfGrids_];
    for (int ii = 0; ii < numOfGrids_; ++ii)
    {
        grid[ii] = getGridAtIndex(index[ii], ii);
    }
    return grid;
}

double MarkovChainParameters::getAlphaAtIndex(unsigned int index)
{
    assert(index < alphaLength_);

    double alphaPosition = alphaLeftBound_ + index * deltaAlpha_;
    return alphaPosition;
}

unsigned int MarkovChainParameters::getNumOfGrids()
{
    return numOfGrids_;
}

unsigned int MarkovChainParameters::getGridLength(unsigned int gridNum)
{
    assert(gridNum < numOfGrids_);
    return gridLength_[gridNum];
}

unsigned int MarkovChainParameters::getAlphaLength()
{
    return alphaLength_;
}

double MarkovChainParameters::getDeltaGrid(unsigned int gridNum)
{
    return deltaGrid_[gridNum];
}

unsigned int MarkovChainParameters::getMaxIterations()
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
                                             unsigned int* gridLength,
                                             double alphaLeftBound,
                                             double alphaRightBound,
                                             unsigned int alphaLength)
{
    // First we check that both grid bound arrays and grid length are all of the same array size
    int counter = 0;
    do
    {
        if (gridLeftBound[counter] != nullptr
            && gridRightBound[counter] != nullptr
            && gridLength[counter] != nullptr)
        {
            counter++;
        }
        else if (gridLeftBound[counter] == nullptr
                && gridRightBound[counter] == nullptr
                && gridLength[counter] == nullptr)
        {
            break;
        }
        else
        {
            // An error here is thrown signifying that the gridLeftBound, gridRightBound and gridLength have differing
            // array sizes.
            std::cout << "ERROR: In " << typeid(this).name() << " on line " << __LINE__ << std::endl;
            assert(false);
        }
    } while(true);

    // Now we ensure appropriate bounds have been given
    assert(alphaLeftBound < alphaRightBound);
    assert(gridLeftBound < gridRightBound);

    // Now we ensure appropriate grid lengths and alpha lengths have been given
    counter = 0;
    while (gridLength[counter] != nullptr)
    {
        assert(gridLength[counter] > 1);
    };
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
    int counter = 0;
    do
    {
        if (gridLeftBound[counter] != nullptr
            && gridRightBound[counter] != nullptr)
        {
            counter++;
        }
        else if (gridLeftBound[counter] == nullptr
                 && gridRightBound[counter] == nullptr)
        {
            break;
        }
        else
        {
            // An error here is thrown signifying that the gridLeftBound and gridRightBound have differing
            // array sizes.
            std::cout << "ERROR: In " << typeid(this).name() << " on line " << __LINE__ << std::endl;
            assert(false);
        }
    } while(true);

    // Now we ensure appropriate bounds have been given
    assert(alphaLeftBound < alphaRightBound);
    assert(gridLeftBound < gridRightBound);

    // Now we ensure appropriate grid deltas and alpha deltas have been given (can't be negative)
    counter = 0;
    while (deltaGrid[counter] != nullptr)
    {
        assert(deltaGrid[counter] > 0);
    };
    assert(deltaAlpha > 0);
}









