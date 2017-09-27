//
// Created by david on 28/09/17.
//

#include "GridParameters.h"
#include <cmath>
#include <cassert>


typedef unsigned int uint;

GridParameters::GridParameters(const double* gridLeftBound,
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

GridParameters::GridParameters(const double* gridLeftBound,
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

GridParameters::~GridParameters()
{
    delete gridLeftBound_;
    delete gridLength_;
    delete deltaGrid_;
}


// GETTERS ------------------------------------------------------------------------ GETTERS //

double GridParameters::getGridAtIndex(unsigned int index, unsigned int gridNum)
{
    assert(gridNum < numOfGridDimensions_ && index < gridLength_[gridNum]);

    double gridPosition = gridLeftBound_[gridNum] + index * deltaGrid_[gridNum];
    return gridPosition;
}

void GridParameters::getGridAtIndex(unsigned int* index, double* outGrid)
{
    for (unsigned int ii = 0; ii < numOfGridDimensions_; ++ii)
    {
        outGrid[ii] = getGridAtIndex(index[ii], ii);
    }
}

unsigned int GridParameters::getNumOfGrids()
{
    return numOfGridDimensions_;
}

unsigned int GridParameters::getGridLength(unsigned int gridIndex)
{
    return gridLength_[gridIndex];
}

double GridParameters::getDeltaGrid(unsigned int gridIndex)
{
    return deltaGrid_[gridIndex];
}


// PRIVATE METHODS -------------------------------------------------------- PRIVATE METHODS //

void GridParameters::assertParameters(const double* gridLeftBound,
                                       const double* gridRightBound,
                                       const uint* gridLength)
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

void GridParameters::assertParameters(const double* gridLeftBound,
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






