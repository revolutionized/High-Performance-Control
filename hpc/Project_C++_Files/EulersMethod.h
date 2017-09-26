//
// Created by David on 26/09/2017.
//

#pragma once

#include <functional>
#include "EulerParameters.h"
#include "MarkovChainApproximation.h"


typedef const std::function<double(double*, double)> fcn2dep;
typedef const std::function<double(double*)> fcn1dep;

class EulersMethod
{
public:
    explicit EulersMethod(EulerParameters& epm);

    ~EulerMethod();

    void solveExact(fcn1dep& fcnDerivative, double initGuess);

    void solve(fcn2dep& fcnDerivative,
               MarkovChainApproximation& mca,
               double initGuess);

private:
    void resetIndices(uint* currentIndices, uint padding);
    unsigned int getSolutionIndexRelative(unsigned int* gridIndices);
    void setUpSolution();
    unsigned int recursionCount(uint dimIndex,
                                uint* indices,
                                uint padding,
                                bool& stillInGrid);
    bool nextRecursiveGrid(uint* currentIndices, uint* previousIndices, uint padding);

    void setUpSolutionDuringSolve(double initGuess);

    std::map<unsigned int[], double>* solution_ = nullptr;
    EulerParameters epm_;

    void recursiveSolve(int currentGrid, unsigned int* gridIndices, unsigned int* previousIndices,
                            fcn1dep& fcnDerivative);
};


