//
// Created by David on 26/09/2017.
//

#pragma once

#include <functional>
#include "EulerParameters.h"
#include "MarkovChainApproximation.h"

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

    explicit EulersMethod(EulerParameters& epm);

    ~EulerMethod();

    void solve(fcn1dep& fcnDerivative, double initGuess);

    void solve(fcn2dep& fcnDerivative,
               MarkovChainApproximation& mca,
               double initGuess);

    void saveSolution(std::ofstream& stream);

    void saveGrid(std::ofstream& stream);

    bool nextRecursiveGrid(uint* currentIndices, uint* previousIndices, uint padding);
private:

    // METHODS ------------------------------------------------------------------------------------------------- METHODS
    void resetIndices(uint* currentIndices, uint padding);
    unsigned int getSolutionIndexRelative(const uint* gridIndices);
    void setUpSolution();
    unsigned int recursionCount(uint dimIndex,
                                uint* indices,
                                uint padding,
                                bool& stillInGrid);


    // FIELDS --------------------------------------------------------------------------------------------------- FIELDS

    std::map<unsigned int[], double>* solution_ = nullptr;
    EulerParameters epm_;

    void recursiveSolve(uint currentGrid,
                        uint* gridIndices,
                        uint* previousIndices,
                        fcn1dep& fcnDerivative);
};


