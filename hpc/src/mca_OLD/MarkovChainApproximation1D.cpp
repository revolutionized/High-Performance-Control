//
// Created by david on 22/08/17.
//

#include "MarkovChainApproximation1D.h"
#include "mca_1Donly/Functions_OLD.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>

//extern "C" {
//// Get c3 library which is in C-code not C++
//#include "c3.h"
//}

using namespace std;

MarkovChainApproximation1D::MarkovChainApproximation1D(double alphaLeftBound,
                                                   double alphaRightBound,
                                                   unsigned int alphaLength,
                                                   double gridLeftBound,
                                                   double gridRightBound,
                                                   unsigned int gridLength,
                                                   double h,
                                                   double initGuess)
        : alphaLength_(alphaLength),
          alphaLeftBound_(alphaLeftBound),
          gridLeftBound_(gridLeftBound),
          h_(h)
{
    assert(alphaLeftBound_ < alphaRightBound);
    assert(gridLeftBound_ < gridRightBound);
    assert(gridLength > 0);
    assert(alphaLength > 0);

    // Set up grid spacing
    deltaGrid_ = (gridRightBound - gridLeftBound_) / (gridLength - 1);
    // Set up alpha discretisation
    deltaAlpha_ = (alphaRightBound - alphaLeftBound_) / (alphaLength - 1);

    // Set allocation sizes for new vectors
    oldV_ = new std::vector<double>(gridLength);
    newV_ = new std::vector<double>(gridLength);
    minAlpha_ = new std::vector<double>(gridLength);

    // Apply initial guess
    for (int i = 0; i < gridLength; ++i)
    {
        (*newV_)[i] = initGuess;
        (*minAlpha_)[i] = initGuess;
    }
    *oldV_ = *newV_; // Copy the vectors

    // Cause BC's aren't really considered, just say the minimum alpha boundaries are zero
    (*minAlpha_)[0] = 0;
    (*minAlpha_)[minAlpha_->size() - 1] = 0;

    // Set relative error to test against and max iterations
    epsErr_ = pow(10.0, -3.0);
    maxIterations_ = 100;
}

MarkovChainApproximation1D::~MarkovChainApproximation1D()
{
    delete oldV_;
    delete newV_;
    delete minAlpha_;
}

void MarkovChainApproximation1D::computeMarkovApproximation(const std::function<double(double, double)>& costFunction,
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

double MarkovChainApproximation1D::B_func(double x, double alpha)
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


double MarkovChainApproximation1D::getGridAtIndex(unsigned int index)
{
//    assert(index < newV_.size());

    double gridPosition = gridLeftBound_ + index * deltaGrid_;
    return gridPosition;
}

double MarkovChainApproximation1D::getAlphaAtIndex(unsigned int index)
{
    assert(index < alphaLength_);

    double alphaPosition = alphaLeftBound_ + index * deltaAlpha_;
    return alphaPosition;
}

double MarkovChainApproximation1D::transitionProb(double x,
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

double MarkovChainApproximation1D::kahanSum(const vector<double>& doubleArray, size_t size_n)
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

void MarkovChainApproximation1D::determineTransitionProbabilities(vector<double>& v_probabilities,
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

void MarkovChainApproximation1D::determineTransitionSummations(vector<double>& v_probabilities,
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

void MarkovChainApproximation1D::determineMinimumAlpha(const vector<double>& v_summed, unsigned int x_index)
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

double MarkovChainApproximation1D::getRelativeError(const vector<double>* firstArr,
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

double MarkovChainApproximation1D::getMarkovControlFunction(double x)
{
    return (*minAlpha_)[getGridIndexClosestTo(x)];
}

unsigned int MarkovChainApproximation1D::getGridIndexClosestTo(double x)
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
