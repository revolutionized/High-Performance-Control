//
// Created by David on 26/09/2017.
//

#pragma once

#include <functional>
#include "multidim/GridParameters.h"
#include "mca/MarkovChainApproximation.h"

#ifdef _WIN32
// New line for Windows
#define NEWL "\r\n"
#else
// New line for Unix
#define NEWL "\n"
#endif

typedef const std::function<double(double*, double)> fcn2dep;
typedef const std::function<double(double*)> fcn1dep;

class EulersMethod
{
public:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

    explicit EulersMethod(GridParameters& epm);

    ~EulersMethod();

    void solve(fcn1dep& fcnDerivative, double initGuess);

    void solve(fcn2dep& fcnDerivative,
               MarkovChainApproximation& mca,
               double initGuess);

    void saveSolution(std::ofstream& stream);

private:


    // FIELDS --------------------------------------------------------------------------------------------------- FIELDS

    GridValue* solution_ = nullptr;
    GridParameters epm_;

    void recursiveSolve(uint currentGrid,
                        GridIndex& gridIndices,
                        GridIndex& previousIndices,
                        fcn1dep& fcnDerivative);

    void recursiveSolve(uint currentGrid,
                        GridIndex& gridIndices,
                        GridIndex& previousIndices,
                        fcn2dep& fcnDerivative,
                        MarkovChainApproximation& mca);

    void test();
};


