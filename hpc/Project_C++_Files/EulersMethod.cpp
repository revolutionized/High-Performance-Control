//
// Created by David on 26/09/2017.
//

#include "EulersMethod.h"

EulersMethod::EulersMethod(EulerParameters& epm)
    : epm_(epm)
{
    setUpSolution();
}

EulersMethod::~EulerMethod()
{
    delete solution_;
}

void EulersMethod::solveExact(fcn1dep& fcnDerivative, double initGuess)
{
    double df;
    unsigned int gridIndices[epm_.getNumOfGrids()];
    unsigned int previousGridIndices[epm_.getNumOfGrids()];

    resetIndices(gridIndices, 0);
    resetIndices(previousGridIndices, 0);
    (*solution_)[gridIndices] = initGuess;

    for (int ii = 0; ii < epm_.getNumOfGrids(); ++ii)
    {
        recursiveSolve(ii, gridIndices, previousGridIndices, fcnDerivative);
    }

    do
    {

    } while (nextRecursiveGrid(gridIndices, previousGridIndices, 0));
}

void EulersMethod::setUpSolution()
{
    // Setup solution memory and insert initial guess
    if (solution_ != nullptr)
    {
        delete solution_;
    }
    else
    {
        unsigned int gridIndices[epm_.getNumOfGrids()];
        solution_ = new std::map<unsigned int[epm_.getNumOfGrids()], double>;
        resetIndices(gridIndices, 0);

        do
        {
            solution_->insert(std::make_pair(gridIndices, 0));
        } while (nextRecursiveGrid(gridIndices, nullptr, 0));
    }
}

void EulersMethod::resetIndices(uint* currentIndices, uint padding)
{
    // We ignore boundaries
    for (unsigned int ii = 0; ii < epm_.getNumOfGrids(); ++ii)
    {
        currentIndices[ii] = padding;
    }
}

bool EulersMethod::nextRecursiveGrid(uint* currentIndices, uint* previousIndices, uint padding)
{
    for (int ii = 0; ii < epm_.getNumOfGrids(); ++ii)
    {
        previousIndices[ii] = currentIndices[ii];
    }

    bool stillInGrid = true;
    recursionCount(epm_.getNumOfGrids() - 1, currentIndices, padding, stillInGrid);
    return stillInGrid;
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

unsigned int EulersMethod::getSolutionIndexRelative(unsigned int* gridIndices)
{
    unsigned int relIndex = 0;
    for (int ii = epm_.getNumOfGrids() - 1; ii >= 0; --ii)
    {
        relIndex += (epm_.getNumOfGrids() - 1 - ii)*gridIndices[ii];
    }
}

void EulersMethod::recursiveSolve(int currentGrid, unsigned int* gridIndices, unsigned int* previousIndices,
                                  fcn1dep& fcnDerivative)
{
    for (int jj = currentGrid; jj < epm_.getGridLength(currentGrid); ++jj)
    {
        recursiveSolve(currentGrid, gridIndices, previousIndices, fcnDerivative);
        double gridLocation[epm_.getNumOfGrids()];
        epm_.getGridAtIndex(gridIndices, gridLocation);
        double df = fcnDerivative(gridLocation);
        (*solution_)[gridIndices] = (*solution_)[previousIndices] + epm_.getDeltaGrid(currentGrid);
    }
}
