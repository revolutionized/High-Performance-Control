//
// Created by David on 26/09/2017.
//

#include <cmath>
#include <cassert>
#include "EulerParameters.h"

typedef unsigned int uint;

EulerParameters::EulerParameters(const double* gridLeftBound,
                                 const double* gridRightBound,
                                 uint* gridLength,
                                 uint dimensions)
    : numOfGridDimensions_(dimensions)
{
    assertParameters(gridLeftBound, gridRightBound, gridLength);
    deltaGrid_ = new double[numOfGridDimensions_];
    gridLeftBound_ = new double[numOfGridDimensions_];
    gridLength_ = new uint[numOfGridDimensions_];
    for (int ii = 0; ii < numOfGridDimensions_; ++ii)
    {
        deltaGrid_[ii] = (gridRightBound[ii] - gridLeftBound[ii]) / (gridLength[ii] - 1);
        gridLeftBound_[ii] = gridLeftBound[ii];
        gridLength_[ii] = gridLength[ii];
    }
}

EulerParameters::EulerParameters(const double* gridLeftBound,
                                 const double* gridRightBound,
                                 const double* deltaGrid,
                                 uint dimensions)
    : numOfGridDimensions_(dimensions)
{
    assertParameters(gridLeftBound, gridRightBound, deltaGrid);
    deltaGrid_ = new double[numOfGridDimensions_];
    gridLeftBound_ = new double[numOfGridDimensions_];
    gridLength_ = new uint[numOfGridDimensions_];
    for (int ii = 0; ii < numOfGridDimensions_; ++ii)
    {
        deltaGrid_[ii] = deltaGrid[ii];
        gridLeftBound_[ii] = gridLeftBound[ii];
        gridLength_[ii] = static_cast<uint>(floor((gridRightBound[ii] - gridLeftBound[ii])/deltaGrid[ii]));
    }
}

EulerParameters::~EulerParameters()
{
    delete gridLeftBound_;
    delete gridLength_;
    delete deltaGrid_;
}


// GETTERS ------------------------------------------------------------------------ GETTERS //

double EulerParameters::getGridAtIndex(unsigned int index, unsigned int gridNum)
{
    assert(gridNum < numOfGridDimensions_ && index < gridLength_[gridNum]);

    double gridPosition = gridLeftBound_[gridNum] + index * deltaGrid_[gridNum];
    return gridPosition;
}

void EulerParameters::getGridAtIndex(unsigned int* index, double* outGrid)
{
    for (unsigned int ii = 0; ii < numOfGridDimensions_; ++ii)
    {
        outGrid[ii] = getGridAtIndex(index[ii], ii);
    }
}

unsigned int EulerParameters::getNumOfGrids()
{
    return numOfGridDimensions_;
}

unsigned int EulerParameters::getGridLength(unsigned int gridIndex)
{
    return gridLength_[gridIndex];
}

double EulerParameters::getDeltaGrid(unsigned int gridIndex)
{
    return deltaGrid_[gridIndex];
}


// PRIVATE METHODS -------------------------------------------------------- PRIVATE METHODS //

void EulerParameters::assertParameters(const double* gridLeftBound,
                                       const double* gridRightBound,
                                       uint* gridLength)
{
    // First we check that both grid bound arrays and grid length are all of the same array size
    if (gridLeftBound[numOfGridDimensions_-1] == nullptr
        || gridRightBound[numOfGridDimensions_-1] == nullptr
        || gridLength[numOfGridDimensions_-1] == nullptr)
    {
        assert(false);
    }

    // Now we ensure appropriate bounds have been given
    assert(gridLeftBound < gridRightBound);

    // Now we ensure appropriate grid lengths have been given
    for (int ii = 0; ii < numOfGridDimensions_; ++ii)
    {
        assert(gridLength[ii] > 1);
    }
}

void EulerParameters::assertParameters(const double* gridLeftBound,
                                       const double* gridRightBound,
                                       const double* deltaGrid)
{
    // First we check that both grid bound arrays and grid length are all of the same array size
    if (gridLeftBound[numOfGridDimensions_-1] == nullptr
        || gridRightBound[numOfGridDimensions_-1] == nullptr
        || deltaGrid[numOfGridDimensions_-1] == nullptr)
    {
        assert(false);
    }

    // Now we ensure appropriate bounds have been given
    assert(gridLeftBound < gridRightBound);

    // Now we ensure appropriate grid lengths have been given
    for (int ii = 0; ii < numOfGridDimensions_; ++ii)
    {
        assert(deltaGrid[ii] > 1);
    }
}






