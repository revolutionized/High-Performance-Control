//
// Created by David on 26/09/2017.
//

#include "EulersMethod.h"


EulersMethod::EulersMethod(EulerParameters& epm)
    : epm_(epm)
{
    setUpSolution();
}

EulersMethod::~EulersMethod()
{
    delete solution_;
}

void EulersMethod::solve(fcn1dep& fcnDerivative, double initGuess)
{
    uint gridIndices[epm_.getNumOfGrids()];
    uint previousGridIndices[epm_.getNumOfGrids()];
    resetIndices(gridIndices, 0);
    resetIndices(previousGridIndices, 0);

    (*solution_)[gridIndices] = initGuess;

    for (uint ii = 0; ii < epm_.getNumOfGrids(); ++ii)
    {
        recursiveSolve(ii, gridIndices, previousGridIndices, fcnDerivative);
    }
}

void EulersMethod::solve(fcn2dep& fcnDerivative, MarkovChainApproximation& mca, double initGuess)
{
    uint gridIndices[epm_.getNumOfGrids()];
    uint previousGridIndices[epm_.getNumOfGrids()];
    resetIndices(gridIndices, 0);
    resetIndices(previousGridIndices, 0);

    (*solution_)[gridIndices] = initGuess;

    for (uint ii = 0; ii < epm_.getNumOfGrids(); ++ii)
    {
        recursiveSolve(ii, gridIndices, previousGridIndices, fcnDerivative, mca);
    }
}


void EulersMethod::saveSolution(std::ofstream& stream)
{
    uint gridIndices[epm_.getNumOfGrids()];
    resetIndices(gridIndices, 0);

    // Fill solution at each point of the grid
    do
    {
        for (uint ii = 0; ii < epm_.getNumOfGrids(); ++ii)
        {
            stream << epm_.getGridAtIndex(gridIndices[ii], ii) << " ";
        }
        stream << (*solution_)[gridIndices] << NEWL;
    } while (nextRecursiveGrid(gridIndices, nullptr, 0));
}


void EulersMethod::saveGrid(std::ofstream& stream)
{
    uint gridIndices[epm_.getNumOfGrids()];
    resetIndices(gridIndices, 0);

    // Fill each grid point
    do
    {
        for (uint ii = 0; ii < epm_.getNumOfGrids(); ++ii)
        {
            stream << epm_.getGridAtIndex(gridIndices[ii], ii) << " ";
        }
        stream << NEWL;
    } while (nextRecursiveGrid(gridIndices, nullptr, 0));
}

// PRIVATE METHODS -------------------------------------------------------- PRIVATE METHODS //

uint EulersMethod::getSolutionIndexRelative(const uint* gridIndices)
{
    uint relIndex = 0;
    for (int ii = epm_.getNumOfGrids() - 1; ii >= 0; --ii)
    {
        relIndex += (epm_.getNumOfGrids() - 1 - ii)*gridIndices[ii];
    }
}

void EulersMethod::recursiveSolve(uint currentGrid,
                                  uint* gridIndices,
                                  uint* previousIndices,
                                  fcn1dep& fcnDerivative)
{
    for (unsigned int jj = currentGrid; jj < epm_.getGridLength(currentGrid); ++jj)
    {
        recursiveSolve(currentGrid, gridIndices, previousIndices, fcnDerivative);
        double gridLocation[epm_.getNumOfGrids()];
        epm_.getGridAtIndex(gridIndices, gridLocation);
        double df = fcnDerivative(gridLocation);
        (*solution_)[gridIndices] = (*solution_)[previousIndices] + epm_.getDeltaGrid(currentGrid)*df;
    }
}

void EulersMethod::recursiveSolve(uint currentGrid,
                                  uint* gridIndices,
                                  uint* previousIndices,
                                  fcn2dep& fcnDerivative,
                                  MarkovChainApproximation& mca)
{
    for (unsigned int jj = currentGrid; jj < epm_.getGridLength(currentGrid); ++jj)
    {
        recursiveSolve(currentGrid, gridIndices, previousIndices, fcnDerivative, mca);
        double gridLocation[epm_.getNumOfGrids()];
        epm_.getGridAtIndex(gridIndices, gridLocation);
        double df = fcnDerivative(gridLocation, mca.getMarkovControlFunction(x));
        (*solution_)[gridIndices] = (*solution_)[previousIndices] + epm_.getDeltaGrid(currentGrid)*df;
    }
}


uint EulersMethod::recursionCount(uint dimIndex, uint* indices, uint padding, bool& stillInGrid)
{
    if (dimIndex >= 0)
    {
        if (indices[dimIndex] < epm_.getGridLength(dimIndex) - padding - 1)
        {
            indices[dimIndex]++;
        }
        else
        {
            indices[dimIndex] = padding;
            dimIndex--;
            recursionCount(dimIndex, indices, padding, stillInGrid);
        }
    }
    else
    {
        stillInGrid = false;
    }
}

void EulersMethod::setUpSolution()
{
    delete solution_;

    // Setup solution memory and insert initial guess
    uint gridIndices[epm_.getNumOfGrids()];
    solution_ = new std::map<uint[epm_.getNumOfGrids()], double>;
    resetIndices(gridIndices, 0);

    do
    {
        solution_->insert(std::make_pair(gridIndices, 0));
    } while (nextRecursiveGrid(gridIndices, nullptr, 0));

}

void EulersMethod::resetIndices(uint* currentIndices, uint padding)
{
    // We ignore boundaries
    for (uint ii = 0; ii < epm_.getNumOfGrids(); ++ii)
    {
        currentIndices[ii] = padding;
    }
}

bool EulersMethod::nextRecursiveGrid(uint* currentIndices, uint* previousIndices, uint padding)
{
    if (previousIndices != nullptr)
    {
        for (int ii = 0; ii < epm_.getNumOfGrids(); ++ii)
        {
            previousIndices[ii] = currentIndices[ii];
        }
    }

    bool stillInGrid = true;
    recursionCount(epm_.getNumOfGrids() - 1, currentIndices, padding, stillInGrid);
    return stillInGrid;
}
