//
// Created by david on 23/09/17.
//

#include "MarkovChainApproximation2.h"
#include "Functions2.h"
#include <cassert>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>

#ifdef _WIN32
// New line for Windows
#define NEWL "\r\n"
#else
// New line for Unix
#define NEWL "\n"
#endif

using std::cout;
using std::string;
using std::fstream;

MarkovChainApproximation2::MarkovChainApproximation2(MarkovChainParameters& mcp,
                                                     double* initGridGuess,
                                                     double* initAlphaGuess,
                                                     unsigned int precision)
        :
{
    // Set precision of std::out
    std::setprecision(precision);

    // Instead of using in-RAM arrays, we write to a fstreams (files), since there may not be enough space on the memory
    for (int ii = 0; ii < mcp.getNumOfGrids(); ++ii)
    {
        file aNewFile;
        aNewFile.stream = new fstream;
        // Set precision
        aNewFile.stream->precision(precision);
        // Set filename
        aNewFile.filename = ".tmpMCA_Old~" + std::to_string(ii);
        aNewFile.stream->open(aNewFile.filename);
        // Apply initial guess
        for (int jj = 0; jj < mcp.getGridLength(ii); ++jj)
        {
            (*aNewFile.stream) << std::to_string(initGridGuess[ii]) << NEWL;
        }
        oldVFile_.push_back(aNewFile);


        file anotherNewFile;
        anotherNewFile.stream = new fstream;
        // Set precision
        anotherNewFile.stream->precision(precision);
        // Set filename
        anotherNewFile.filename = ".tmpMCA_New~" + std::to_string(ii);
        anotherNewFile.stream->open(anotherNewFile.filename);
        // Apply initial guess
        for (int jj = 0; jj < mcp.getGridLength(ii); ++jj)
        {
            (*anotherNewFile.stream) << std::to_string(initGridGuess[ii]) << NEWL;
        }
        newVFile_.push_back(anotherNewFile);
    }
    for (int ii = 0; ii < mcp.getNumOfAlphas(); ++ii)
    {
        file aNewFile;
        aNewFile.stream = new fstream;
        // Set precision
        aNewFile.stream->precision(precision);
        // Set filename
        aNewFile.filename = ".tmpMCA_alpha~" + std::to_string(ii);
        aNewFile.stream->open(aNewFile.filename);
        // Apply initial guess
        for (int jj = 0; jj < mcp.getAlphaLength(ii); ++jj)
        {
            (*aNewFile.stream) << std::to_string(initAlphaGuess[ii]) << NEWL;
        }
        minAlphaFile_.push_back(aNewFile);
    }

    // Assert that all files were correctly opened
    assert(allStreamsOpen());

    /*
    // Apply initial guess' for dynamic programming equations, and use same values for minimum alpha values
    for (int ii = 0; ii < mcp.getNumOfGrids(); ++ii)
    {
        // Start back at top and create new column
        newVFile_.seekg(0);

        for (int jj = 0; jj < mcp.getGridLength(ii); ++jj)
        {
            std::string line;
            // Remove the end newline component
            getline(newVFile_, line);
            line.erase(std::remove(line.begin(), line.end(), NEWL), line.end());
            // Add on to the end of the line the next guess for that grid dynamic equation
            line += " " + std::to_string(initGridGuess[ii]);
            newVFile_ << line << NEWL;
        }
    }
    for (int ii = 0; ii < mcp.getNumOfAlphas(); ++ii)
    {
        // Start back at top and create new column
        minAlphaFile_.seekg(0);

        // Cause BC's aren't really considered, just say the minimum alpha boundaries are zero
        std::string line;
        // Remove the end newline component
        getline(minAlphaFile_, line);
        line.erase(std::remove(line.begin(), line.end(), NEWL), line.end());
        line += " " + std::to_string(0.0);
        minAlphaFile_ << line;

        for (int jj = 1; jj < mcp.getAlphaLength(ii) - 1; ++jj)
        {
            std::string line;
            // Remove the end newline component
            getline(minAlphaFile_, line);
            line.erase(std::remove(line.begin(), line.end(), NEWL), line.end());
            // Add on to the end of the line the next init guess
            line += " " + std::to_string(initAlphaGuess[ii]);
            newVFile_ << line << NEWL;
        }

        // end boundary same as before
        getline(minAlphaFile_, line);
        line.erase(std::remove(line.begin(), line.end(), NEWL), line.end());
        line += " " + std::to_string(0.0);
        minAlphaFile_ << line;
    }
    */
}

MarkovChainApproximation2::~MarkovChainApproximation2()
{
    delete oldV_;
    delete newV_;
    delete minAlpha_;

    // Close and remove files used for dynamic programming equations
    for (int ii = 0; ii < oldVFile_.size(); ++ii)
    {
        oldVFile_[ii].deallocate();
        delete oldVFile_[ii];
    }
    for (int ii = 0; ii < newVFile_.size(); ++ii)
    {
        newVFile_[ii].deallocate();
        delete newVFile_[ii];
    }
    for (int ii = 0; ii < minAlphaFile_.size(); ++ii)
    {
        minAlphaFile_[ii].deallocate();
        delete minAlphaFile_[ii];
    }
}

void MarkovChainApproximation2::computeMarkovApproximation(const std::function<double(double, double)>& costFunction,
                                                          const std::function<double(double)>& sigmaFunction)
{
    vector<double> v_probabilities(3);
    vector<double> v_summed(alphaLength_);
    double relErr;
    unsigned int iterations = 0;

    std::setprecision(4); // For showing relative error approximate to 4 decimals
    cout << "==== Starting Markov Chain Approximation technique ====" << endl;
    do
    {
        oldV_->swap(*newV_);

        for (unsigned int xi = 1; xi < newV_->size() - 1; ++xi)
        {
            // Solve all summations for different alpha values
            determineTransitionSummations(v_probabilities, xi, v_summed, costFunction, sigmaFunction);

            // Find the minimum alpha term for the DPE
            determineMinimumAlpha(v_summed, xi);
        }

        relErr = getRelativeError(newV_, oldV_, static_cast<unsigned int>(newV_->size()));
        cout << "Markov chain iteration complete. Relative error: " << relErr;
        cout << ". Iteration: " << iterations << std::endl;
        iterations++;

    } while (relErr > epsErr_ && iterations < maxIterations_);
}

double MarkovChainApproximation2::B_func(double x, double alpha)
{
    // B_func results already produced, just need to find max
//    return max(abs(b_result));
    double result = A * x + B * alpha;
    if (result < 0)
    {
        result = -result;
    }
    return result;
}



double MarkovChainApproximation2::transitionProb(double x,
                                                double y,
                                                double alpha,
                                                double den,
                                                const std::function<double(double)>& sigmaFunction)
{
    double num = pow(sigmaFunction(x), 2.0);

    double b_part = A * x + B * alpha;
    if (y > x)
    {
        b_part = b_part > 0 ? b_part : 0;
    } else
    {
        b_part = -b_part > 0 ? -b_part : 0;
    }
    num = num + h_ * b_part;

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

double MarkovChainApproximation2::kahanSum(const vector<double>& doubleArray, size_t size_n)
{
    // See https://en.wikipedia.org/wiki/Kahan_summation_algorithm

    double sum = 0.0;
    double c = 0.0;
    for (int i = 0; i < size_n; ++i)
    {
        double y = doubleArray[i] - c;
        double t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    return sum;
}

void MarkovChainApproximation2::determineTransitionProbabilities(vector<double>& v_probabilities,
                                                                unsigned int xi,
                                                                double alpha,
                                                                double den,
                                                                const std::function<double(double)>& sigmaFunction)
{
    std::array<double, 3> p{0, 0, 1};
    double x = getGridAtIndex(xi);
    std::array<double, 3> y{x - h_, x + h_, x};

    unsigned int yi;
    for (yi = 0; yi < 2; ++yi)
    {
        p[yi] = transitionProb(x, y[yi], alpha, den, sigmaFunction);
        // Need to match y iteration to the same x iteration value
        unsigned int vyi;
        if (yi == 0)
        {
            vyi = xi - 1;
        } else
        {
            vyi = xi + 1;
        }

        v_probabilities[yi] = p[yi] * (*oldV_)[vyi];
    }
    yi = 2;
    p[yi] = 1 - p[0] - p[1];
    v_probabilities[yi] = p[yi] * (*oldV_)[xi];
}

void MarkovChainApproximation2::determineTransitionSummations(vector<double>& v_probabilities,
                                                             unsigned int xi,
                                                             vector<double>& v_summed,
                                                             const std::function<double(double, double)>& costFunctionK,
                                                             const std::function<double(double)>& sigmaFunction)
{
    double x = getGridAtIndex(xi);

    for (unsigned int ai = 0; ai < alphaLength_; ++ai)
    {
        double alpha = getAlphaAtIndex(ai);
        double den = pow(sigmaFunction(x), 2.0) + h_ * B_func(x, alpha);
        double delta_t = pow(h_, 2.0) / den;
        double k = costFunctionK(x, alpha);

        determineTransitionProbabilities(v_probabilities, xi, alpha, den, sigmaFunction);

        // Solve dynamic programming equation
        v_summed[ai] = kahanSum(v_probabilities, 3) + delta_t * k;
    }
}

void MarkovChainApproximation2::determineMinimumAlpha(const vector<double>& v_summed, unsigned int x_index)
{
    // Find minimum v_summed and its index
    int minIndex = 0;
    for (int i = 1; i < alphaLength_; ++i)
    {
        if (v_summed[i] < v_summed[minIndex])
        {
            minIndex = i;
        }
    }
    // Set the minimum value of V for that x-index
    (*newV_)[x_index] = v_summed[minIndex];
    // Store the minimum alpha value
    (*minAlpha_)[x_index] = getAlphaAtIndex(static_cast<unsigned int>(minIndex));
}

double MarkovChainApproximation2::getRelativeError(const vector<double>* firstArr,
                                                  const vector<double>* secondArr,
                                                  unsigned int arrLength)
{
    double temp;
    double maxValue = 0.0;

    for (int i = 0; i < arrLength; ++i)
    {
        // Perform infinity norm
        temp = (*firstArr)[i] - (*secondArr)[i];
        if (abs(temp) > maxValue)
        {
            maxValue = abs(temp);
        }
    }

    return maxValue;
}

double MarkovChainApproximation2::getMarkovControlFunction(double x)
{
    return (*minAlpha_)[getGridIndexClosestTo(x)];
}

unsigned int MarkovChainApproximation2::getGridIndexClosestTo(double x)
{
    unsigned int closestGridIndex = 0;
    double minDistance = abs(getGridAtIndex(0) - x);

    for (unsigned int i = 1; i < newV_->size(); ++i)
    {
        double gridLoc = getGridAtIndex(i);
        if (abs(gridLoc - x) < minDistance)
        {
            closestGridIndex = i;
            minDistance = abs(gridLoc - x);
        }
    }

    return closestGridIndex;
}

bool MarkovChainApproximation2::allStreamsOpen()
{
    for (auto&& oldVFile : oldVFile_)
    {
        if (!oldVFile.stream->is_open())
        {
            return false;
        }
    }
    for (auto&& newVFile : newVFile_)
    {
        if (!newVFile.stream->is_open())
        {
            return false;
        }
    }
    for (auto&& minAlphaFile : minAlphaFile_)
    {
        if (!minAlphaFile.stream->is_open())
        {
            return false;
        }
    }
    return true;
}

