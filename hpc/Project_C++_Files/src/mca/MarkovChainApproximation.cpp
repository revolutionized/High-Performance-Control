//
// Created by david on 23/09/17.
//

#include "MarkovChainApproximation.h"
#include <cassert>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;

#ifdef _WIN32
// New line for Windows
#define NEWL "\r\n"
#else
// New line for Unix
#define NEWL "\n"
#endif

MarkovChainApproximation::MarkovChainApproximation(MarkovChainParameters& mcp,
                                                   double initStateGuess,
                                                   unsigned int precision,
                                                   bool memoryModeRAM)
        : memoryModeRAM_(memoryModeRAM),
          mcp_(mcp)
{
    // Set precision of std::out
    std::setprecision(precision);

    if (!memoryModeRAM_)
    {
        // Instead of using in-RAM arrays, we write to a fstreams (files),
        // since there may not be enough space on the memory
        for (uint ii = 0; ii < mcp.getNumOfGrids(); ++ii)
        {
            file aNewFile;
            aNewFile.stream = new std::fstream;
            // Set precision
            aNewFile.stream->precision(precision);
            // Set filename
            aNewFile.filename = ".tmpMCA_Old~" + std::to_string(ii);
            aNewFile.stream->open(aNewFile.filename);
            // Apply initial guess
            for (uint jj = 0; jj < mcp.getGridLength(ii); ++jj)
            {
                (*aNewFile.stream) << std::to_string(initStateGuess) << NEWL;
            }
            oldVFile_->push_back(aNewFile);


            file anotherNewFile;
            anotherNewFile.stream = new std::fstream;
            // Set precision
            anotherNewFile.stream->precision(precision);
            // Set filename
            anotherNewFile.filename = ".tmpMCA_New~" + std::to_string(ii);
            anotherNewFile.stream->open(anotherNewFile.filename);
            // Apply initial guess
            for (uint jj = 0; jj < mcp.getGridLength(ii); ++jj)
            {
                (*anotherNewFile.stream) << std::to_string(initStateGuess) << NEWL;
            }
            newVFile_->push_back(anotherNewFile);
        }
        minAlphaFile_.stream = new std::fstream;
        // Set precision
        minAlphaFile_.stream->precision(precision);
        // Set filename
        minAlphaFile_.filename = ".tmpMCA_alpha";
        minAlphaFile_.stream->open(minAlphaFile_.filename);

        // Assert that all files were correctly opened
        assert(allStreamsOpen());
    }
    else
    {
        // Set allocation sizes for new vectors
        oldV_ = new GridValue;
        newV_ = new GridValue;
        minAlpha_ = new GridValue;

        GridIndex mainGridIndices(mcp.getNumOfGrids());

        do
        {
            GridIndex oldVGridIndices(mainGridIndices);
            oldV_->insert(std::make_pair(oldVGridIndices, initStateGuess));
            GridIndex newVGridIndices(mainGridIndices);
            minAlpha_->insert(std::make_pair(newVGridIndices,initStateGuess));
        } while (mainGridIndices.nextGridElement(mcp));
    }
}

bool MarkovChainApproximation::allStreamsOpen()
{
    for (auto&& oldVFile : *oldVFile_)
    {
        if (!oldVFile.stream->is_open())
        {
            return false;
        }
    }
    for (auto&& newVFile : *newVFile_)
    {
        if (!newVFile.stream->is_open())
        {
            return false;
        }
    }
    return minAlphaFile_.stream->is_open();
}

MarkovChainApproximation::~MarkovChainApproximation()
{
    if (memoryModeRAM_)
    {
        /*
        for (uint ii = 0; ii < oldV_->size(); ++ii)
        {
            delete oldV_[ii];
        }
        delete oldV_;
        for (uint ii = 0; ii < newV_->size(); ++ii)
        {
            delete newV_[ii];
        }
        delete newV_;
        for (uint ii = 0; ii < minAlpha_->size(); ++ii)
        {
            delete minAlpha_[ii];
        }
         */
        delete oldV_;
        delete newV_;
        delete minAlpha_;
    }
    else
    {
        // Close and remove files used for dynamic programming equations
        for (auto it = oldVFile_->begin() ; it != oldVFile_->end(); ++it)
        {
            it->deallocate();
        }
        oldVFile_->clear();
        for (auto it = newVFile_->begin() ; it != newVFile_->end(); ++it)
        {
            it->deallocate();
        }
        newVFile_->clear();

        minAlphaFile_.deallocate();
    }
}

void MarkovChainApproximation::computeMarkovApproximation(fcn2dep& costFunction,
                                                          fcn2dep& driftFunction,
                                                          fcn1dep& diffFunction)
{
    // Assign memory for v_summed now to avoid having to reassign the memory in loops
    vector<double> v_summed(mcp_.getAlphaLength());
    // Current index of grid (e.g. for a 3x3 system we are at (0,1,3) or something
    GridIndex gridIndices(mcp_.getNumOfGrids());

    // These are the two parameters we compare against to know if we should keep looping or not
    double relErr;
    uint iterations = 0;

    std::setprecision(4); // For showing relative error approximate to 4 decimals
    cout << "==== Starting Markov Chain Approximation ====" << endl;

    if (memoryModeRAM_)
    {
        do
        {
            oldV_->swap(*newV_);
            gridIndices.resetToOrigin(1);
            
            do
            {
                // Solve all summations for different alpha values
                solveTransitionSummations(v_summed, gridIndices, costFunction, driftFunction, diffFunction);

                // Find the minimum alpha term for the DPE
                determineMinimumAlpha(v_summed, gridIndices);

            } while (gridIndices.nextGridElement(mcp_, 1));
            
            relErr = getRelativeError();
            
            cout << "Markov chain iteration complete. Relative error: " << relErr;
            cout << ". Iteration: " << iterations << endl;
            iterations++;
        } while (relErr > mcp_.getRelativeError() && iterations < mcp_.getMaxIterations());
    }
    else
    {

    }

    cout << "==== Finished Markov Chain Approximation ====" << endl;
}

void MarkovChainApproximation::solveTransitionSummations(std::vector<double>& v_summed,
                                                         GridIndex& gridIndices,
                                                         fcn2dep& costFunctionK,
                                                         fcn2dep& driftFunction,
                                                         fcn1dep& diffFunction)
{
    // Nodes can jump either side of a plane, plus the possibility of not moving
    uint numOfPossibilities = mcp_.getNumOfGrids()*2 + 1;
    vector<double> v_probabilities(numOfPossibilities);

    double gridLocation[mcp_.getNumOfGrids()];
    for (uint ii = 0; ii < mcp_.getNumOfGrids(); ++ii)
    {
        gridLocation[ii] = mcp_.getGridAtIndex(gridIndices.getIndexOfDim(ii), ii);
    }

    for (uint ai = 0; ai < mcp_.getAlphaLength(); ++ai)
    {
        double alpha = mcp_.getAlphaAtIndex(ai);
        double den = pow(diffFunction(gridLocation), 2.0) + mcp_.getH() * B_func(driftFunction, gridLocation, alpha);
        double delta_t = pow(mcp_.getH(), 2.0) / den;
        double k = costFunctionK(gridLocation, alpha);

        solveTransitionProbabilities(v_probabilities, gridIndices, alpha, den, driftFunction, diffFunction);

        // Solve dynamic programming equation
        v_summed[ai] = kahanSum(v_probabilities, numOfPossibilities) + delta_t * k;
    }
}

double MarkovChainApproximation::B_func(fcn2dep& driftFunction,
                                         double* gridLocation,
                                         double alpha)
{
    double result = driftFunction(gridLocation, alpha);
    if (result < 0)
    {
        result = -result;
    }
    return result;
}

void MarkovChainApproximation::solveTransitionProbabilities(std::vector<double>& v_probabilities,
                                                            GridIndex& gridIndices,
                                                            double alpha,
                                                            double den,
                                                            fcn2dep& driftFunction,
                                                            fcn1dep& diffFunction)
{
    vector<double> p(v_probabilities.size());
    // Begin with probability of getting anything but itself as 0
    for (uint ii = 0; ii < p.size() - 1; ++ii)
    {
        p[ii] = 0;
    }
    p[p.size() - 1] = 1;
    // p[end][1] is unused

    // Determine 'y' -> the positions that this node can move to
    vector<double*> y(2);
    y[0] = new double[mcp_.getNumOfGrids()];
    y[1] = new double[mcp_.getNumOfGrids()];

    for (uint ii = 0; ii < mcp_.getNumOfGrids(); ++ii)
    {
        // One grid position moved up on the current dimension
        GridIndex upperGridIndex(gridIndices);
        upperGridIndex.setGridIndexOfDim(ii, upperGridIndex.getIndexOfDim(ii) + 1);
        // One grid position moved down on the current dimension
        GridIndex lowerGridIndex(gridIndices);
        lowerGridIndex.setGridIndexOfDim(ii, lowerGridIndex.getIndexOfDim(ii) - 1);

        for (uint jj = 0; jj < mcp_.getNumOfGrids(); ++jj)
        {
            if (ii == jj)
            {
                continue;
            }
            y[0][jj] = mcp_.getGridAtIndex(gridIndices.getIndexOfDim(jj), jj);
            y[1][jj] = mcp_.getGridAtIndex(gridIndices.getIndexOfDim(jj), jj);
        }
        y[0][ii] = mcp_.getGridAtIndex(upperGridIndex.getIndexOfDim(ii), ii) - mcp_.getH();
        y[1][ii] = mcp_.getGridAtIndex(lowerGridIndex.getIndexOfDim(ii), ii) + mcp_.getH();


        // Need to match the y with same value/position in the old dynamic equation
        uint yi = ii*2;

        // Now use the values of y to determine transition probabilities (for both the negative and positve part
        p[yi] = transitionProb(y[0], true, alpha, den, driftFunction, diffFunction);
        p[yi+1] = transitionProb(y[1], false, alpha, den, driftFunction, diffFunction);

        // Calculate probability based on old dynamic equation values at that grid position
        v_probabilities[yi] = p[yi]*(*oldV_)[upperGridIndex];
        v_probabilities[yi+1] = p[yi+1]*(*oldV_)[lowerGridIndex];
    }
    // Get probability of staying and not moving (complement to moving)
    for (unsigned ii = 0; ii < mcp_.getNumOfGrids(); ii++)
    {
        p[p.size()-1] = 1;
        for (uint jj = 0; jj < p.size() - 1; ++jj)
        {
            p[p.size()-1] -= p[jj];
        }

        // TODO: Remove this if no problems occur, this is just for error checking
        if (p[p.size()-1] < 0)
        {
            cout << "ERROR: Remaining probability less than 0" << endl;
        }
    }
    v_probabilities[v_probabilities.size()-1] = p[p.size()-1]*(*oldV_)[gridIndices];

    // Clear up memory allocated
    delete y[0];
    delete y[1];
}


double MarkovChainApproximation::transitionProb(double* x,
                                                bool upperJump,
                                                double alpha,
                                                double den,
                                                fcn2dep& driftFunction,
                                                fcn1dep& diffFunction)
{
    double num = pow(diffFunction(x), 2.0);
    double b_part = driftFunction(x, alpha);

    if (upperJump)
    {
        b_part = b_part > 0 ? b_part : 0;
    } else
    {
        b_part = -b_part > 0 ? -b_part : 0;
    }
    num = num + mcp_.getH() * b_part;

    double result = num / den;
    if (result > 1)
    {
        cout << "Error: Transition probability is > 1" << std::endl;
    } else if (result < 0)
    {
        cout << "Error: Transition probability is < 0 " << std::endl;
    }
    return result;
}

double MarkovChainApproximation::kahanSum(const vector<double>& doubleArray, size_t size_n)
{
    // See https://en.wikipedia.org/wiki/Kahan_summation_algorithm

    double sum = 0.0;
    double c = 0.0;
    for (uint i = 0; i < size_n; ++i)
    {
        double y = doubleArray[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    return sum;
}

void MarkovChainApproximation::determineMinimumAlpha(const vector<double>& v_summed, GridIndex& gridIndices)
{
    // Find minimum v_summed and its index
    uint minIndex = 0;
    for (uint i = 1; i < mcp_.getAlphaLength(); ++i)
    {
        if (v_summed[i] < v_summed[minIndex])
        {
            minIndex = i;
        }
    }

    // Set the minimum value of V for that current index
    (*newV_)[gridIndices] = v_summed[minIndex];
    // Store the minimum alpha value
    (*minAlpha_)[gridIndices] = mcp_.getAlphaAtIndex(minIndex);
}

double MarkovChainApproximation::getRelativeError()
{
    double temp = 0.0;

    GridIndex gridIndices(mcp_.getNumOfGrids());
    gridIndices.resetToOrigin();
    
    do
    {
        // Perform Frobenius tensor norm
        temp += pow(fabs((*oldV_)[gridIndices]) - fabs((*newV_)[gridIndices]), 2.0);
    } while (gridIndices.nextGridElement(mcp_));
        
    return sqrt(temp);
}

double MarkovChainApproximation::getMarkovControlFunction(double* gridLocation)
{
    return (*minAlpha_)[getGridIndicesClosestTo(gridLocation)];
}

GridIndex MarkovChainApproximation::getGridIndicesClosestTo(double* gridLocation)
{
    // Set initial guess to the origin
    GridIndex currentGridIndex(mcp_.getNumOfGrids());
    currentGridIndex.resetToOrigin();
    GridIndex closestGridIndex(currentGridIndex);

    double minDistance = 0;
    for (uint ii = 0; ii < mcp_.getNumOfGrids(); ++ii)
    {
        minDistance += pow(mcp_.getDeltaGrid(ii) * 1.5, 2.0);
    }
    minDistance = sqrt(minDistance);

    // Basically we do a relative distance check between each point (just using Pythagoras theorem)
    do
    {
        double currentGridLoc[mcp_.getNumOfGrids()];
        double currentDistance = 0;
        for (uint ii = 0; ii < mcp_.getNumOfGrids(); ++ii)
        {
            currentGridLoc[ii] = mcp_.getGridAtIndex(currentGridIndex.getIndexOfDim(ii), ii);
            if (currentGridLoc[ii] > 0 && gridLocation[ii] < 0)
            {
                currentDistance += pow(gridLocation[ii] - currentGridLoc[ii], 2.0);
            } else
            {
                currentDistance += pow(currentGridLoc[ii] - gridLocation[ii], 2.0);
            }
        }
        currentDistance = sqrt(currentDistance);

        if (currentDistance < minDistance)
        {
            for (uint ii = 0; ii < mcp_.getNumOfGrids(); ++ii)
            {
                closestGridIndex.setGridIndexOfDim(ii, currentGridIndex.getIndexOfDim(ii));
            }
            minDistance = currentDistance;
        }

    } while (currentGridIndex.nextGridElement(mcp_));

    return closestGridIndex;
}
