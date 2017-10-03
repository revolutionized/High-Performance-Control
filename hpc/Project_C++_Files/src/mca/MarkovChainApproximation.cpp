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

    cout << "Allocating needed memory" << NEWL;

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
        oldV_ = new GridIndices;
        newV_ = new GridIndices;
        minAlpha_ = new GridIndices;

        GridIndex mainGridIndices(mcp.getNumOfGrids());

        uint progressCount = 0;

        do
        {
            GridIndex oldVGridIndices(mainGridIndices);
            oldV_->insert(std::make_pair(oldVGridIndices, initStateGuess));
            GridIndex newVGridIndices(mainGridIndices);
            newV_->insert(std::make_pair(newVGridIndices, initStateGuess));
            GridIndex alphaGridIndices(mainGridIndices);
            minAlpha_->insert(std::make_pair(alphaGridIndices,0.0));

            // We say its the total number of grid elements + 1 since the 1% remaining is for vectors below
            printProgress(static_cast<float>(++progressCount)/ static_cast<float>(mcp.getGridLengthAccumulation() + 1));
        } while (mainGridIndices.nextGridElement(mcp));
    }

    // Set allocation for computation convenience vectors to avoid having to reassign the memory in loops
    pHat_ = new vector<double>(3); // 3 because always probability to move up one, down one, or stay in same place
    vProbs_ = new vector<double>(3); // Same as for pHat.
    vSummed_ = new vector<vector<double>>(mcp.getAlphaLength());
    for (auto& ii : *vSummed_)
    {
        ii.reserve(mcp_.getNumOfGrids());
    }

    // Show that memory allocation has been finished
    cout << "\r" << "Finished allocating memory." << "\r";
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
        delete oldV_;
        delete newV_;
        delete minAlpha_;
    }
    else
    {
        // Close and remove files used for dynamic programming equations
        for (auto& it : *oldVFile_)
        {
            it.deallocate();
        }
        oldVFile_->clear();
        for (auto& it : *newVFile_)
        {
            it.deallocate();
        }
        newVFile_->clear();

        minAlphaFile_.deallocate();
    }

    // All modes utilise these vectors during computation
    delete pHat_;
    delete vProbs_;
    delete vSummed_;
}

void MarkovChainApproximation::computeMarkovApproximation(d_fcn2dep& costFunction,
                                                          v_fcn2dep& driftFunction,
                                                          v_fcn1dep& diffFunction)
{
    // Current index of grid (e.g. for a 3x3 system we are at (0,1,3) or something
    GridIndex gridIndices(mcp_.getNumOfGrids());

    // These are the two parameters we compare against to know if we should keep looping or not
    double relErr;
    uint iterations = 0;

    std::setprecision(4); // For showing relative error approximate to 4 decimals
    cout << "==== Starting Markov Chain Approximation ====" << NEWL;

    if (memoryModeRAM_)
    {
        do
        {
            // The newV always relies on the oldV so we just replace the two pointers
            oldV_->swap(*newV_);

            // Add a padding of 1 since we don't want to be on boundaries
            gridIndices.resetToOrigin(1);

            uint progressCount = 0;

            do
            {
                // Solve all summations for different alpha values
                solveTransitionSummations(gridIndices, costFunction, driftFunction, diffFunction);

                // Find the minimum alpha term for the DPE
                determineMinimumAlpha(gridIndices);

                // Here we print progress of a single iteration
                printProgress(++progressCount/static_cast<float>(mcp_.getGridLengthAccumulation()));

            } while (gridIndices.nextGridElement(mcp_, 1));

            // Calculate relative error between newV and oldV
            relErr = getRelativeError();

            // Display progress
            cout << "Markov chain iteration complete. Relative error: " << relErr;
            cout << ". Iteration: " << iterations + 1 << NEWL;
            iterations++;

        } while (relErr > mcp_.getRelativeError() && iterations < mcp_.getMaxIterations());
    }
    else
    {
        // TODO: File writing mode
    }

    // Print one new line since our iterations use "\r" which returns to the beginning, and then
    // two new lines at the end to have a gap after MCA conclusion
    cout << NEWL << "==== Finished Markov Chain Approximation ====" << NEWL << NEWL;
}

void MarkovChainApproximation::solveTransitionSummations(const GridIndex& gridIndices,
                                                         d_fcn2dep& costFunctionK,
                                                         v_fcn2dep& driftFunction,
                                                         v_fcn1dep& diffFunction)
{
    double gridLocation[mcp_.getNumOfGrids()];
    for (uint ii = 0; ii < mcp_.getNumOfGrids(); ++ii)
    {
        gridLocation[ii] = mcp_.getGridAtIndex(gridIndices.getIndexOfDim(ii), ii);
    }

    for (uint ai = 0; ai < mcp_.getAlphaLength(); ++ai)
    {
        double alpha = mcp_.getAlphaAtIndex(ai);

        // Denominator for each dimension
        double den[mcp_.getNumOfGrids()];
        diffFunction(gridLocation, den);
        double bFunc[mcp_.getNumOfGrids()];
        B_func(driftFunction, gridLocation, alpha, bFunc);
        for (uint ii = 0; ii < mcp_.getNumOfGrids(); ++ii)
        {
            den[ii] = pow(den[ii], 2.0) + mcp_.getH()*bFunc[ii];
        }

        // Delta_t for each dimension
        double delta_t[mcp_.getNumOfGrids()];
        for (uint ii = 0; ii < mcp_.getNumOfGrids(); ++ii)
        {
            delta_t[ii] = pow(mcp_.getH(), 2.0) / den[ii];
        }

        // k for each dimension
        double k = costFunctionK(gridLocation, alpha);

        for (uint ii = 0; ii < mcp_.getNumOfGrids(); ++ii)
        {
            solveTransitionProbabilities(gridIndices, ii, alpha, den[ii], driftFunction, diffFunction);
            // Solve dynamic programming equation
            (*vSummed_)[ai][ii] = kahanSum(vProbs_) + delta_t[ii] * k;
        }
    }
}

void MarkovChainApproximation::B_func(v_fcn2dep& driftFunction,
                                      double* gridLocation,
                                      double alpha,
                                      double* out)
{
    driftFunction(gridLocation, alpha, out);
    for (uint ii = 0; ii < mcp_.getNumOfGrids(); ++ii)
    {
        out[ii] = fabs(out[ii]);
    }
}

void MarkovChainApproximation::solveTransitionProbabilities(const GridIndex& currentIndices,
                                                            uint currentDimension,
                                                            double alpha,
                                                            double den,
                                                            v_fcn2dep& driftFunction,
                                                            v_fcn1dep& diffFunction)
{
    // Begin with probability of getting anything but itself as 0 (meaning itself is 1)
    int stay = 2;
    (*pHat_)[stay] = 1;

    double gridLoc[mcp_.getNumOfGrids()];
    double num[mcp_.getNumOfGrids()];
    double b_part[mcp_.getNumOfGrids()];
    // Get
    mcp_.getGridAtIndices(currentIndices, gridLoc);

    // Solve numerator
    diffFunction(gridLoc, num);
    // Solve b_part
    driftFunction(gridLoc, alpha, b_part);

    // First consider the case it moves up the current dimension plane
    int lower = 0;
    (*pHat_)[lower] = transitionProb(false, num[currentDimension], den, b_part[currentDimension]);

    // Now consider the case it moves down the current dimension plane
    int upper = 1;
    (*pHat_)[upper] = transitionProb(true, num[currentDimension], den, b_part[currentDimension]);

    // Finally we consider the case that it does not move at all (which is just the complement)
    (*pHat_)[stay] -= (*pHat_)[lower] + (*pHat_)[upper];

    // TODO: Remove this if no problems occur, this is just for error checking
    if ((*pHat_)[stay] < 0)
    {
        cout << "ERROR: Remaining probability less than 0" << endl;
    }

    // We need to match the current index to the same that y is pointing to (see that y represents the grid location,
    // but we need the indices that represent that location
    GridIndex lowerIndices = GridIndex(currentIndices);
    GridIndex upperIndices = GridIndex(currentIndices);
    lowerIndices.setGridIndexOfDim(currentDimension, currentIndices.getIndexOfDim(currentDimension) - 1);
    upperIndices.setGridIndexOfDim(currentDimension, currentIndices.getIndexOfDim(currentDimension) + 1);

    // Now we can calculate the dynamic probabilities
    (*vProbs_)[lower] = (*pHat_)[lower]*(*oldV_)[lowerIndices];
    (*vProbs_)[upper] = (*pHat_)[upper]*(*oldV_)[upperIndices];
    (*vProbs_)[stay] = (*pHat_)[stay]*(*oldV_)[currentIndices];
}


double MarkovChainApproximation::transitionProb(bool upperJump, double num, double den, double b_part)
{
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

double MarkovChainApproximation::kahanSum(const vector<double>* array)
{
    // See https://en.wikipedia.org/wiki/Kahan_summation_algorithm

    double sum = 0.0;
    double c = 0.0;
    for (double i : *array)
    {
        double y = i - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    return sum;
}

void MarkovChainApproximation::determineMinimumAlpha(const GridIndex& gridIndices)
{
    // We search for the alpha that minimises the dynamic programming equation, and return it's index
    uint minIndex = 0;
    // For every alpha considered, we need to compare which creates the minimum tensor norm
    double minNorm = pow(10.0, 3.0);
    double currentNorm;


    for (uint ii = 0; ii < mcp_.getAlphaLength(); ++ii)
    {
        currentNorm = 0.0;

        for (int jj = 0; jj < mcp_.getNumOfGrids(); ++jj)
        {
            currentNorm += pow((*vSummed_)[ii][jj], 2.0);
        }
        currentNorm = sqrt(currentNorm);


        if (currentNorm < minNorm)
        {
            minIndex = ii;
            minNorm = currentNorm;
        }
    }

    // Set the minimum value of V for that current index
    (*newV_)[gridIndices] = 0.0;
    for (int ii = 0; ii < mcp_.getNumOfGrids(); ++ii)
    {
        (*newV_)[gridIndices] += (*vSummed_)[minIndex][ii];
    }

    // Store the minimum alpha value
    (*minAlpha_)[gridIndices] = mcp_.getAlphaAtIndex(minIndex);
}

double MarkovChainApproximation::getRelativeError()
{
    double temp = 0.0;

    GridIndex gridIndices(mcp_.getNumOfGrids());
    gridIndices.resetToOrigin();

    for (auto& oldV : *oldV_)
    {
        // Perform Frobenius tensor norm
        temp += pow(fabs(oldV.second) - fabs((*newV_)[oldV.first]), 2.0);
    }

//    do
//    {
//        // Perform Frobenius tensor norm
//        temp += pow(fabs((*oldV_)[gridIndices]) - fabs((*newV_)[gridIndices]), 2.0);
//    } while (gridIndices.nextGridElement(mcp_));

    temp /= 10;
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

    double minDistance = pow(10.0, 3.0);

    // Basically we do a relative distance check between each point (just using Pythagoras theorem /
    // tensor Frobenius norm)
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

void MarkovChainApproximation::printAlpha(std::ofstream& stream)
{
    GridIndex newGrid(1);
    newGrid.resetToOrigin();

    do
    {
        stream << (*minAlpha_)[newGrid] << NEWL;
    } while (newGrid.nextGridElement(mcp_));
}

void MarkovChainApproximation::printProgress(float progress)
{
    int barWidth = 50;

    std::cout << "[";
    auto pos = static_cast<int>(barWidth * progress);
    for (int ii = 0; ii < barWidth; ++ii) {
        if (ii < pos)
        {
            std::cout << "=";
        }
        else if (ii == pos)
        {
            std::cout << ">";
        }
        else
        {
            std::cout << " ";
        }
    }
    std::cout << "] " << "%" << int(progress * 100.0) << "\r";
    std::cout.flush();
}
