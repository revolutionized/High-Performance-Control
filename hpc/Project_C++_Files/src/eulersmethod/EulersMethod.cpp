//
// Created by David on 26/09/2017.
//

#include "EulersMethod.h"


EulersMethod::EulersMethod(GridParameters& epm)
    : epm_(epm)
{
    // Setup solution memory and insert initial guess
    GridIndex mainGridIndices(epm_.getNumOfGrids());
    solution_ = new GridValue;
    mainGridIndices.resetToOrigin();

    do
    {
        GridIndex newGridIndices(mainGridIndices);
        solution_->insert(std::make_pair(newGridIndices, 0));
    } while (mainGridIndices.nextGridElement(epm_));
}

EulersMethod::~EulersMethod()
{
    delete solution_;
}

void EulersMethod::solve(fcn1dep& fcnDerivative, double initGuess)
{
    GridIndex gridIndices(epm_.getNumOfGrids());
    gridIndices.resetToOrigin();
    GridIndex previousGridIndices(gridIndices);

    (*solution_)[gridIndices] = initGuess;

    for (uint ii = 0; ii < epm_.getNumOfGrids(); ++ii)
    {
        recursiveSolve(ii, gridIndices, previousGridIndices, fcnDerivative);
    }
}

void EulersMethod::solve(fcn2dep& fcnDerivative, MarkovChainApproximation& mca, double initGuess)
{
    GridIndex gridIndices(epm_.getNumOfGrids());
    gridIndices.resetToOrigin();
    GridIndex previousGridIndices(gridIndices);

    (*solution_)[gridIndices] = initGuess;

    for (uint ii = 0; ii < epm_.getNumOfGrids(); ++ii)
    {
        recursiveSolve(ii, gridIndices, previousGridIndices, fcnDerivative, mca);
    }
}


void EulersMethod::saveSolution(std::ofstream& stream)
{
    GridIndex gridIndices(epm_.getNumOfGrids());
    gridIndices.resetToOrigin();

    // Fill solution at each point of the grid
    do
    {
        for (uint ii = 0; ii < epm_.getNumOfGrids(); ++ii)
        {
            stream << epm_.getGridAtIndex(gridIndices.getIndexOfDim(ii), ii) << " ";
        }
        stream << (*solution_)[gridIndices] << NEWL;
    } while (gridIndices.nextGridElement(epm_));
}


// PRIVATE METHODS -------------------------------------------------------- PRIVATE METHODS //

void EulersMethod::recursiveSolve(uint currentGrid,
                                  GridIndex& gridIndices,
                                  GridIndex& previousIndices,
                                  fcn1dep& fcnDerivative)
{
    for (unsigned int jj = currentGrid; jj < epm_.getGridLength(currentGrid); ++jj)
    {
        recursiveSolve(currentGrid+1, gridIndices, previousIndices, fcnDerivative);
        double gridLocation[epm_.getNumOfGrids()];
        for (uint kk = 0; kk < epm_.getNumOfGrids(); ++kk)
        {
            gridLocation[kk] = epm_.getGridAtIndex(gridIndices.getIndexOfDim(kk), kk);
        }
        double df = fcnDerivative(gridLocation);
        (*solution_)[gridIndices] = (*solution_)[previousIndices] + epm_.getDeltaGrid(currentGrid)*df;
    }
}


void EulersMethod::recursiveSolve(uint currentGrid,
                                  GridIndex& gridIndices,
                                  GridIndex& previousIndices,
                                  fcn2dep& fcnDerivative,
                                  MarkovChainApproximation& mca)
{
    for (unsigned int jj = currentGrid; jj < epm_.getGridLength(currentGrid); ++jj)
    {
        recursiveSolve(currentGrid+1, gridIndices, previousIndices, fcnDerivative, mca);
        double gridLocation[epm_.getNumOfGrids()];
        for (uint kk = 0; kk < epm_.getNumOfGrids(); ++kk)
        {
            gridLocation[kk] = epm_.getGridAtIndex(gridIndices.getIndexOfDim(kk), kk);
        }
        double df = fcnDerivative(gridLocation, mca.getMarkovControlFunction(gridLocation));
        (*solution_)[gridIndices] = (*solution_)[previousIndices] + epm_.getDeltaGrid(currentGrid)*df;
    }
}

