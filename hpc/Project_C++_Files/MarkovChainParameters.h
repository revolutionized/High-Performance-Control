//
// Created by david on 23/09/17.
//

#pragma once

#include <cassert>

struct MarkovChainParameters
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
    MarkovChainParameters(double* gridLeftBound,
                          double* gridRightBound,
                          unsigned int* gridLength,
                          double* alphaLeftBound,
                          double* alphaRightBound,
                          unsigned int* alphaLength,
                          double h);

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
    MarkovChainParameters(double* gridLeftBound,
                          double* gridRightBound,
                          double* deltaGrid,
                          double* alphaLeftBound,
                          double* alphaRightBound,
                          double* deltaAlpha,
                          double h);

    ~MarkovChainParameters();

    // GETTERS ------------------------------------------------------------------------------------------------- GETTERS


private:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

    /// Asserts pro
    void assertParameters(double* gridLeftBound,
                          double* gridRightBound,
                          double* deltaGrid,
                          double* alphaLeftBound,
                          double* alphaRightBound,
                          double* deltaAlpha);

    void assertParameters(double* gridLeftBound,
                          double* gridRightBound,
                          unsigned int* gridLength,
                          double* alphaLeftBound,
                          double* alphaRightBound,
                          unsigned int* alphaLength);


    // FIELDS --------------------------------------------------------------------------------------------------- FIELDS

    /// Length of the alpha array (could also use sizeof, but this is for convenience)
    unsigned int* alphaLength_;

    unsigned int* gridLength_;

    /// \brief Lower bound of the grid array.
    ///
    /// Calculating grid points starts from the lower bound and works it way up (so no
    /// need to store the upper bound)
    double* gridLeftBound_;
    /// \brief Lower bound of the alpha array.
    ///
    /// Calculating alpha points starts from the lower bound and works it way up (so no
    /// need to store the upper bound)
    double* alphaLeftBound_;
    /// Size of steps between grid array values
    double* deltaGrid_;
    /// Size of steps between alpha array values
    double* deltaAlpha_;
    /// Scalar approximation parameter that determines
    double h_;
};



