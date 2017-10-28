//
// Created by david on 23/09/17.
//

#include "MarkovChainApproximation.h"
#include <cassert>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <utilities/ProgressBar.h>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;


#ifdef _WIN32
// New line for Windows
#ifndef NEWL
#define NEWL "\r\n"
#endif

#else

// New line for Unix
#ifndef NEWL
#define NEWL "\n"
#endif // NEWL

#endif //_WIN32

// Also MAXFLOAT not defined in Unix
#ifndef MAXFLOAT
#define MAXFLOAT 1.79769e+308
#endif

#define PHAT_TEST(phat) { if ((phat) > 1.0 || (phat) < 0.0) cout << "ERROR Phat > 1.0 or < 0.0 " << NEWL; }


MarkovChainApproximation::MarkovChainApproximation(MarkovChainParameters& mcp,
                                                   double initStateGuess,
                                                   unsigned int precision,
                                                   bool memoryModeRAM)
        : memoryModeRAM_(memoryModeRAM),
          mcp_(mcp)
{
    // Set precision of std::out
    std::setprecision(precision);

    cout << "Allocating needed memory";

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

        do
        {
            GridIndex oldVGridIndices(mainGridIndices);
            oldV_->insert(std::make_pair(oldVGridIndices, initStateGuess));
            GridIndex newVGridIndices(mainGridIndices);
            newV_->insert(std::make_pair(newVGridIndices, initStateGuess));
            GridIndex alphaGridIndices(mainGridIndices);
            minAlpha_->insert(std::make_pair(alphaGridIndices,0.0));
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
    cout << "\r" << "Finished allocating memory.";
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
                                                          v_fcn1dep& diffFunction,
                                                          v_fcn1dep& diffMatrix)
{
    // Current index of grid (e.g. for a 3x3 system we are at (0,1,3) or something
    GridIndex gridIndex(mcp_.getNumOfGrids());

    // These are the two parameters we compare against to know if we should keep looping or not
    double relErr;
    uint iterations = 0;

    std::setprecision(4); // For showing relative error approximate to 4 decimals
    cout << "\r==== Performing Markov Chain Approximation ====" << NEWL;

    float totalCountToComplete = mcp_.getGridLengthAccumulation() + 1;

    uint dimensions = mcp_.getNumOfGrids();
    // Define transition probabilities
//    vector<double> phat;
//    phat.reserve(static_cast<unsigned long>(int(pow(3.0, dimensions))));
//    auto test = phat.size();

    if (memoryModeRAM_)
    {
        do
        {
            // The newV always relies on the oldV so we just replace the two pointers
            oldV_->swap(*newV_);

            // Add a padding of 1 since we don't want to be on boundaries
            gridIndex.resetToOrigin(1);

            // Monitor progress
            uint progressCount = 0;
            int percentage_complete = 0;

            bool weakDiagonals = false;

            do
            {
                // Get covariance matrix
                double location[dimensions];
                mcp_.getGridAtIndices(gridIndex, location);
                double sigmaMatrix[dimensions*dimensions];
                diffMatrix(location, sigmaMatrix);
                double covariance[dimensions][dimensions];

                // Perform covariance * covariance' <- where ' signifies transpose
                for (int ii = 0; ii < dimensions; ++ii)
                {
                    for (int jj = 0; jj < dimensions; ++jj)
                    {
                        covariance[ii][jj] = 0.0;

                        for (int kk = 0; kk < dimensions; ++kk)
                        {
                            covariance[ii][jj] += sigmaMatrix[dimensions*kk + ii] * sigmaMatrix[dimensions*jj + kk];
                            auto test = covariance[ii][jj];
                            auto stop = false;
                        }
                    }
                }

                // Test for weak diagonals of covariance matrix
                for (int ii = 0; ii < dimensions; ++ii)
                {
                    double diag = covariance[ii][ii];
                    for (int jj = 0; jj < dimensions; ++jj)
                    {
                        if (jj == ii)
                        {
                            continue;
                        }
                        diag -= covariance[ii][jj];
                    }
                    if (diag < 0)
                    {
                        weakDiagonals = true;
                    }
                }


                auto test = mcp_.getAlphaLength();
                // Loops through each alpha value
                for (int ai = 0; ai < mcp_.getAlphaLength(); ++ai)
                {
                    double alpha = mcp_.getAlphaAtIndex(static_cast<uint>(ai));


                    double k = costFunction(location, alpha);

                    double deltaT[dimensions]; // interpolation interval

                    double b[dimensions]; // drift
                    driftFunction(location, alpha, b);

                    double Qhat[dimensions]; // normalising coefficient
                    double diagSum = 0.0;
                    for (int ii = 0; ii < dimensions; ++ii)
                    {
                        diagSum += covariance[ii][ii];
                    }

                    for (int ii = 0; ii < dimensions; ++ii)
                    {
                        Qhat[ii] = diagSum + mcp_.getH()*fabs(b[ii]);
                        // Get off-diagonal components
                        for (int jj = 0; jj < dimensions && jj != ii; ++jj)
                        {
                            Qhat[ii] -= covariance[ii][jj]/2.0;
                        }

                        deltaT[ii] = pow(mcp_.getH(), 2.0) / Qhat[ii];
                    }

                    // Transition probabilities
                    for (int ii = 0; ii < dimensions; ++ii)
                    {
                        double phat[2];

                        // Probability of transitioning either side of current dimension
//                        int thisDimensionPos = ii*int(pow(3.0, dimensions - 1));
                        int thisDimensionPos = 0;
//                        int thisDimensionNeg = thisDimensionPos + 1;
                        int thisDimensionNeg = 1;

                        double numPos = covariance[ii][ii]/2.0;
                        double numNeg = numPos;
                        for (int kk = 0; kk < dimensions && kk != ii; ++kk)
                        {
                            numPos -= fabs(covariance[ii][kk])/2.0;
                            numNeg -= fabs(covariance[ii][kk])/2.0;
                        }
                        double bPos = b[ii] > 0 ? b[ii] : 0;
                        double bNeg = -b[ii] > 0 ? -b[ii] : 0;
                        numPos += mcp_.getH()*bPos;
                        numNeg += mcp_.getH()*bNeg;


                        phat[thisDimensionPos] = numPos/Qhat[ii]; PHAT_TEST(phat[thisDimensionPos]);
                        phat[thisDimensionNeg] = numNeg/Qhat[ii]; PHAT_TEST(phat[thisDimensionNeg]);

                        // Need to consider indices related to the current transition
                        GridIndex posGridIndex(gridIndex);
                        posGridIndex.setGridIndexOfDim(static_cast<uint>(ii),
                                                       posGridIndex.getIndexOfDim(static_cast<uint>(ii)) + 1);
                        GridIndex negGridIndex(gridIndex);
                        negGridIndex.setGridIndexOfDim(static_cast<uint>(ii),
                                                       negGridIndex.getIndexOfDim(static_cast<uint>(ii)) - 1);

                        // Need to find the related gridIndex of moving up or down
                        (*vSummed_)[ai][ii] = k*deltaT[ii];
                        (*vSummed_)[ai][ii] += phat[thisDimensionPos]*(*oldV_)[posGridIndex];
                        (*vSummed_)[ai][ii] += phat[thisDimensionNeg]*(*oldV_)[negGridIndex];

                        // Now the remaining probability is that it stays
                        double phat_stay = 1.0 - phat[thisDimensionNeg] - phat[thisDimensionNeg];

                        // Probability of transitioning either both +ve or both -ve along two axes
                        for (int jj = 0; jj < dimensions; ++jj)
                        {
                            if (jj == ii)
                            {
                                continue;
                            }
                            // Positive along both grid axes
                            double num = covariance[ii][jj] > 0.0 ? covariance[ii][jj] : 0.0;
                            double den = 2*Qhat[ii];
                            double posposPhat = num/den; PHAT_TEST(posposPhat);
                            GridIndex posposIndex(gridIndex);
                            posposIndex.setGridIndexOfDim(static_cast<uint>(ii),
                                                               posposIndex.getIndexOfDim(static_cast<uint>(ii)) + 1);
                            posposIndex.setGridIndexOfDim(static_cast<uint>(jj),
                                                          posposIndex.getIndexOfDim(static_cast<uint>(jj)) + 1);
                            (*vSummed_)[ai][ii] += posposPhat*(*oldV_)[posposIndex];

                            // Negative along both grid axes
                            double negnegPhat = num/den; PHAT_TEST(negnegPhat);
                            GridIndex negnegIndex(gridIndex);
                            negnegIndex.setGridIndexOfDim(static_cast<uint>(ii),
                                                          negnegIndex.getIndexOfDim(static_cast<uint>(ii)) - 1);
                            negnegIndex.setGridIndexOfDim(static_cast<uint>(jj),
                                                          negnegIndex.getIndexOfDim(static_cast<uint>(jj)) - 1);
                            (*vSummed_)[ai][ii] += negnegPhat*(*oldV_)[negnegIndex];

                            // Positive along one axis and negative along the other (first way)
                            num = -covariance[ii][jj] > 0.0 ? -covariance[ii][jj] : 0.0;
                            double posnegPhat = num/den; PHAT_TEST(posnegPhat);
                            GridIndex posnegIndex(gridIndex);
                            posnegIndex.setGridIndexOfDim(static_cast<uint>(ii),
                                                          posnegIndex.getIndexOfDim(static_cast<uint>(ii)) + 1);
                            posnegIndex.setGridIndexOfDim(static_cast<uint>(jj),
                                                          posnegIndex.getIndexOfDim(static_cast<uint>(jj)) - 1);
                            (*vSummed_)[ai][ii] += posnegPhat*(*oldV_)[posnegIndex];

                            // Positive along one axis and negative along the other (second way)
                            double negposPhat = num/den; PHAT_TEST(negposPhat);
                            GridIndex negposIndex(gridIndex);
                            negposIndex.setGridIndexOfDim(static_cast<uint>(ii),
                                                          negposIndex.getIndexOfDim(static_cast<uint>(ii)) - 1);
                            negposIndex.setGridIndexOfDim(static_cast<uint>(jj),
                                                          negposIndex.getIndexOfDim(static_cast<uint>(jj)) + 1);
                            (*vSummed_)[ai][ii] += negposPhat*(*oldV_)[negposIndex];

                            // Calculate remainder
                            phat_stay -= posposPhat;
                            phat_stay -= negnegPhat;
                            phat_stay -= posnegPhat;
                            phat_stay -= negposPhat;
                        }

                        (*vSummed_)[ai][ii] += phat_stay*(*oldV_)[gridIndex];
                    }
                }




                // Solve all summations for different alpha values
//                solveTransitionSummations(gridIndex, costFunction, driftFunction, diffFunction);

                // Find the minimum alpha term for the DPE
                determineMinimumAlpha(gridIndex);


                // If progress has jumped up more than 1% then redraw the progress bar
                // To reduce multiple prints of same percentage, just skip until a different percentage is reached
                auto prog = int(static_cast<float>(++progressCount*100.0/totalCountToComplete));
                if (prog != percentage_complete)
                {
                    percentage_complete = prog;
                    // We say its the total number of grid elements + 1 since the 1% remaining is for vectors below
                    utils::printProgress(percentage_complete);
                }

            } while (gridIndex.nextGridElement(mcp_, 1));

            // Calculate relative error between newV and oldV
            relErr = getRelativeError();

            // Display progress
            cout << "Relative error: " << relErr;
            cout << "  Iteration: " << iterations + 1 << NEWL;
            iterations++;
            // Display errors if any
            if (weakDiagonals)
            {
                cout << "  warning: Weak covariance diagonals were found in last iteration. "
                        "Approximation may be be poor." << NEWL;
            }

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
    // gridLocation basically represents the actual x, y, z (etc.) value that we are at
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

        double gridLoc[mcp_.getNumOfGrids()];
        mcp_.getGridAtIndices(gridIndices, gridLoc);

        // Solve initial component of numerator
        double num[mcp_.getNumOfGrids()];
        diffFunction(gridLoc, num);
        // Solve 'b_part'
        double b_part[mcp_.getNumOfGrids()];
        driftFunction(gridLoc, alpha, b_part);


        for (uint ii = 0; ii < mcp_.getNumOfGrids(); ++ii)
        {
            double currentNum = num[ii];
            double currentbPart = b_part[ii];
            solveTransitionProbabilities(gridIndices,
                                         ii,
                                         alpha,
                                         den[ii],
                                         driftFunction,
                                         diffFunction,
                                         currentNum,
                                         currentbPart);
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

void
MarkovChainApproximation::solveTransitionProbabilities(const GridIndex& currentIndices,
                                                       uint currentDimension,
                                                       double alpha,
                                                       double den,
                                                       v_fcn2dep& driftFunction,
                                                       v_fcn1dep& diffFunction,
                                                       double currentNumerator,
                                                       double currentBpart)
{
    // Begin with probability of getting anything but itself as 0 (meaning itself is 1)
    int stay = 2;
    (*pHat_)[stay] = 1;

    // First consider the case it moves up the current dimension plane
    int lower = 0;
    (*pHat_)[lower] = transitionProb(false, currentNumerator, den, currentBpart);

    // Now consider the case it moves down the current dimension plane
    int upper = 1;
    (*pHat_)[upper] = transitionProb(true, currentNumerator, den, currentBpart);

    // Finally we consider the case that it does not move at all (which is just the complement)
    (*pHat_)[stay] -= (*pHat_)[lower] + (*pHat_)[upper];

    // TODO: Remove this if no problems occur, this is just for error checking
    if ((*pHat_)[stay] < 0)
    {
        double val = currentIndices.getIndexOfDim(currentDimension);

        auto test_lower = (*pHat_)[lower];
        auto test_upper = (*pHat_)[upper];
        auto test_stay = (*pHat_)[stay];
        cout << "ERROR: Remaining probability less than 0" << endl;

        auto test = transitionProb(false, currentNumerator, den, currentBpart);
        auto stop = false;
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
    double temp;
    double maxInf = 0.0;

    GridIndex gridIndices(mcp_.getNumOfGrids());
    gridIndices.resetToOrigin();

    for (auto& oldV : *oldV_)
    {
        // Perform Frobenius tensor norm

        temp = fabs(oldV.second - (*newV_)[oldV.first]);
        maxInf = temp > maxInf ? temp : maxInf;
    }

    return maxInf;
}

double MarkovChainApproximation::getMarkovControlFunction(double* gridLocation)
{
    return (*minAlpha_)[getGridIndicesClosestTo(gridLocation)];
}

GridIndex MarkovChainApproximation::getGridIndicesClosestTo(double* gridLocation)
{
    // Set initial guess to the origin
//    GridIndex currentGridIndex(mcp_.getNumOfGrids());
//    currentGridIndex.resetToOrigin();
    GridIndex closestGridIndex(mcp_.getNumOfGrids());
    closestGridIndex.resetToOrigin();

    double minDistance = MAXFLOAT;

    // Basically we do a relative distance check between each point (just using Pythagoras theorem /
    // tensor Frobenius norm)
    for (auto& minAlpha : *minAlpha_)
    {
        double euclidianDistance = 0.0;

        // Find current grid location
        double currentGridLocation[mcp_.getNumOfGrids()];
        for (uint ii = 0; ii < mcp_.getNumOfGrids(); ++ii)
        {
            currentGridLocation[ii] = mcp_.getGridAtIndex(minAlpha.first.getIndexOfDim(ii), ii);
            // Compare it to the location we are given
            euclidianDistance += pow(currentGridLocation[ii] - gridLocation[ii], 2.0);
        }
        // Find total euclidean distance
        euclidianDistance = sqrt(euclidianDistance);

        if (euclidianDistance < minDistance)
        {
            minDistance = euclidianDistance;
            closestGridIndex = minAlpha.first;
        }

    }
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


