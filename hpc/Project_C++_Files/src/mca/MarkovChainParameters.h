//
// Created by david on 23/09/17.
//

#pragma once

#include "multidim/GridParameters.h"
#include <cassert>

struct MarkovChainParameters : public GridParameters
{
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

    /// \brief The class must be started with all the appropriate bounds.
    ///
    /// Incorrect set bounds will cause assert errors. The arrays (pointers) passed to this class are not managed by
    /// this class, but are rather copied and the copies are managed internally by this class.
    /// \param alphaLeftBound The lower bound for the range of alpha values we want to test
    /// \param alphaRightBound The upper bound for the range of alpha values we want to test
    /// \param alphaLength The total length of alphas to test (ie. if lower bound = 0.0, upper bound = 1.0, and the
    /// length is 5, then we would expect to see the alpha vector look something like { 0.0, 0.25, 0.5, 0.75, 1.0 }).
    /// This constructor is similar to calling in MATLAB: linspace(alphaLeftBound, alphaRightBound, alphaLength).
    /// \param gridLeftBound The lower bound for the arbitrary range that we are covering in this 1-dimensional problem
    /// \param gridRightBound The upper bound for the arbitrary range that we are covering in this 1-dimensional problem
    /// \param gridLength Same as the alphaLength variable, but applies to the 1-dimensional grid
    /// \param h This is the level of detail that the Markov Approximation considers (for convergence we want h <
    /// grid spacing (delta x)
    /// \param initGuess The initial guess for the value of the function we are approximating (at x = 0)
        explicit MarkovChainParameters(const double* gridLeftBound,
                                       const double* gridRightBound,
                                       unsigned int* gridLength,
                                       double alphaLeftBound,
                                       double alphaRightBound,
                                       unsigned int alphaLength,
                                       double h,
                                       unsigned int dimensions);

    /// \brief The class must be started with all the appropriate bounds.
    ///
    /// Incorrect set bounds will cause assert errors. The arrays (pointers) passed to this class are not managed by
    /// this class, but are rather copied and the copies are managed internally by this class.
    /// \param alphaLeftBound The lower bound for the range of alpha values we want to test
    /// \param alphaRightBound The upper bound for the range of alpha values we want to test
    /// \param deltaAlpha The size of spacing between alpha guesses (ie. if lower bound = 0.0, upper bound = 1.0, and
    /// the deltaAlpha is 0.25, then we would expect to see the alpha vector look something like { 0.0, 0.25, 0.5, 0.75,
    /// 1.0 }). This constructor is similar to calling in MATLAB: alphaLeftBound:deltaAlpha:alphaRightBound.
    /// \param gridLeftBound The lower bound for the arbitrary range that we are covering in this 1-dimensional problem
    /// \param gridRightBound The upper bound for the arbitrary range that we are covering in this 1-dimensional problem
    /// \param deltaGrid Same as the deltaAlpha variable, but applies to the n-dimensional grid
    /// \param h This is the level of detail that the Markov Approximation considers (for convergence we want h <
    /// grid spacing (delta x)
    /// \param initGuess The initial guess for the value of the function we are approximating (at x = 0)
    MarkovChainParameters(const double* gridLeftBound, const double* gridRightBound,
                          const double* deltaGrid, double alphaLeftBound, double alphaRightBound,
                          double deltaAlpha, double h, unsigned int dimensions);


    // SETTERS ------------------------------------------------------------------------------------------------- SETTERS

    void setMaxIterations(unsigned int maxIters);

    void setMinError(double epsErr);

    // GETTERS ------------------------------------------------------------------------------------------------- GETTERS

    /// \brief Returns value of alpha array at particular index.
    ///
    /// Technically there is no array for alpha. Instead of storing the alpha as another large array and obvious (in
    /// the sense that it's one increment after the other), this function just calculates and returns the alpha value at
    /// the specified index.
    /// \param index The index along the alpha array. Must be >= 0 and less than the total length of the alpha array.
    double getAlphaAtIndex(unsigned int index);

    unsigned int getAlphaLength();

    unsigned int getMaxIterations();

    double getRelativeError();

    double getH();

private:
    // FIELDS --------------------------------------------------------------------------------------------------- FIELDS

    /// Length of the alpha array (could also use sizeof, but this is for convenience)
    unsigned int alphaLength_;
    /// \brief Lower bound of the alpha array.
    ///
    /// Calculating alpha points starts from the lower bound and works it way up (so no
    /// need to store the upper bound)
    double alphaLeftBound_;
    /// Size of steps between alpha array values
    double deltaAlpha_;
    /// Scalar approximation parameter that determines
    double h_;
    /// \brief Smallest relative error before markov approximation exits.
    ///
    /// The markov chain continues to repeat until its relative error (norm_inf |oldV - newV|) is less than this value.
    /// The smaller this value, the higher the precision, but the greater number of iterations it will take.
    double epsErr_ = 0.001;
    /// \brief Total number of iterations before markov approximation exits.
    ///
    /// If the max number of iterations has been reached and the relative error still isn't less than epsErr_, then the
    /// markov approximation method will finish
    unsigned int maxIterations_ = 100;
};



