//
// Created by David Dos Santos on 22/08/17.
//
#pragma once


#include <functional>
#include <vector>

/// \brief Markov Approximation Method (MCA) for 1-dimensional problems.
/// 
/// This is the so called 'library' for solving a nonlinear stochastic problem using the Markov Approximation Method. It
/// does not solve anything more than a one-dimensional problem.
class MarkovChainApproximation1D
{
public:
    /// \brief The class must be started with all the appropriate bounds. 
    ///
    /// Incorrect set bounds will cause assert errors.
    /// \param alphaLeftBound The lower bound for the range of alpha values we want to test
    /// \param alphaRightBound The upper bound for the range of alpha values we want to test
    /// \param alphaLength The total length of alphas to test (ie. if lower bound = 0.0, upper bound = 1.0, and the
    /// length is 10, then we would expect to see the alpha vector look something like { 0.0, 0.25, 0.5, 0.75, 1.0 })
    /// \param gridLeftBound The lower bound for the arbitrary range that we are covering in this 1-dimensional problem
    /// \param gridRightBound The upper bound for the arbitrary range that we are covering in this 1-dimensional problem
    /// \param gridLength Same as the alphaLength variable, but applies to the 1-dimensional grid
    /// \param h This is the level of detail that the Markov Approximation considers (for convergence we want h <
    /// grid spacing (delta x)
    /// \param initGuess The initial guess for the value of the function we are approximating (at x = 0)
    explicit MarkovChainApproximation1D(double alphaLeftBound,
                                      double alphaRightBound,
                                      unsigned int alphaLength,
                                      double gridLeftBound,
                                      double gridRightBound,
                                      unsigned int gridLength,
                                      double h,
                                      double initGuess);

    /// \brief Deallocates all memory that this class needed
    ~MarkovChainApproximation1D();

    /// \brief Main function to call when ready to compute the function approximation using the MCA method.
    /// \param costFunction The cost function is given as a std::function and takes two doubles (x, alpha) and returns 
    /// the value as a double.
    /// \param sigmaFunction The sigma function is often called the ?advection term?. For this one-dimensional case it
    /// represents the noise / stochastic disturbance being introduced to the scenario.
    void computeMarkovApproximation(const std::function<double(double, double)>& costFunction,
                                    const std::function<double(double)>& sigmaFunction);

    /// \brief Finds the most optimal control value at a specific location along arbitrary one-dimensional grid.
    ///
    /// \pre Must use computeMarkovApproximation prior to using this function or your return values will be incorrect
    /// \param x The location along the arbitrary grid that you wish to find the most optimal control at this point. 
    /// This given x location can be anywhere, but for proper results it should exist within the bounds given at the
    /// construction of this class.
    /// \return The value of the optimal control parameter at the given x location. Or rather, at the location closest 
    /// to x.
    double getMarkovControlFunction(double x);

private:
    /// \brief Returns value of grid array at particular index
    ///
    /// Technically there is no array for the grid. Instead of storing the grid as another large array and obvious (in 
    /// the sense that it's one increment after the other), this function just calculates and returns the grid at the 
    /// specified index.
    /// \param index The index along the grid array. Must be >= 0 and less than the total length of the grid array.
    double getGridAtIndex(unsigned int index);

    /// \brief Returns value of alpha array at particular index.
    ///
    /// Technically there is no array for alpha. Instead of storing the alpha as another large array and obvious (in 
    /// the sense that it's one increment after the other), this function just calculates and returns the alpha value at 
    /// the specified index.
    /// \param index The index along the alpha array. Must be >= 0 and less than the total length of the alpha array.
    double getAlphaAtIndex(unsigned int index);

    /// \brief B function with no control dependency
    /// See "Numerical Methods for Stochastic Control Problems in Continuous Time" (H. J. Kushner) pg. 99 where it gives
    /// B(x) = max_{\alpha\el U} \abs{b(x,\alpha)}.
    double B_func(double x, double alpha);

    /// \brief Returns the value of $ p^{h}(x, y | \alpha) $
    /// \param x Current value of grid location
    /// \param y Transition state (x + h or x - h)
    /// \param alpha Current value of alpha
    /// \param den Denominator $ \sigma(x)^2 + hB(x) $
    /// \param sigmaFunction Same sigma function given at the computeMarkovApproximation call
    double transitionProb(double x,
                          double y,
                          double alpha,
                          double den,
                          const std::function<double(double)>& sigmaFunction);

    /// \brief Summation of double array elements using Kahan formula (which ensures rounding errors aren't lost)
    /// \param doubleArray The array whose elements are to be summed
    /// \param size_n Total size of the given array
    /// \return The total value of all the elements summed together
    double kahanSum(const std::vector<double>& doubleArray, size_t size_n);

    /// \brief Computes the transition values. This is a convenience function for splitting the code into smaller chunks
    /// \param v_probabilities This is the vector that gets filled with the transition values (being the transition
    /// probability multiplied by the previous value of the dynamic programming equation at the position).
    /// \param xi Current index of the grid
    /// \param alpha Current value of alpha 
    /// \param den The denominator for $ \Delta t $ is the same here as it is for $ p^{h} $, so we pass the value of the
    /// denominator in to save re-computation over each iteration.
    /// \param sigmaFunction Same sigma function given at the computeMarkovApproximation call
    void determineTransitionProbabilities(std::vector<double>& v_probabilities,
                                          unsigned int xi,
                                          double alpha,
                                          double den,
                                          const std::function<double(double)>& sigmaFunction);

    /// \brief Loops through each alpha and caluclates dyn prog eqn summation for that particular alpha. 
    /// \param v_probabilities This is the possible transition states that can be taken on (so its of length 3 for the 
    /// one-dimensional case). We pass this in as a parameter instead of constructing it in the function to avoid the 
    /// time needed to allocate and deallocate memory for this vector (instead it is already allocated prior and we just
    /// reuse the vector as many times as needed).
    /// \param xi Current index along the grid
    /// \param v_summed This is where all the summations are stored
    /// \param costFunctionK Same cost function given at the computeMarkovApproximation call
    /// \param sigmaFunction Same sigma function given at the computeMarkovApproximation call
    void determineTransitionSummations(std::vector<double>& v_probabilities,
                                       unsigned int xi,
                                       std::vector<double>& v_summed,
                                       const std::function<double(double, double)>& costFunctionK,
                                       const std::function<double(double)>& sigmaFunction);

    /// \brief Picks the alpha that minimises the dynamic programming equation at the specific grid location. 
    ///
    /// This is a convenience function for splitting the code into smaller chunks.
    /// \param v_summed The values of the dynamic programming equation at each alpha value
    /// \param x_index The current index of the grid
    void determineMinimumAlpha(const std::vector<double>& v_summed, unsigned int x_index);

    /// \brief Returns relative error between two arrays.
    ///
    /// The relative error is calculated as the infinity norm of the absolute difference between the two arrays - thus 
    /// the two arrays must be the same length.
    /// \param arrLength The length of the two arrays (both arrays should be length arrLength)
    double getRelativeError(const std::vector<double>* firstArr,
                            const std::vector<double>* secondArr,
                            unsigned int arrLength);

    /// \brief Searches and finds the index of the grid array that most closely matches with the value of x given
    /// \param x Can be anything, but for proper use should be between the lower and upper bounds of the grid
    unsigned int getGridIndexClosestTo(double x);


    // FIELDS --------------------------------------------------------------------------------------------------- FIELDS

    /// Length of the alpha array (could also use sizeof, but this is for convenience)
    unsigned int alphaLength_;
    /// \brief Lower bound of the grid array. 
    ///
    /// Calculating grid points starts from the lower bound and works it way up (so no 
    /// need to store the upper bound)
    double gridLeftBound_;
    /// \brief Lower bound of the alpha array.
    ///
    /// Calculating alpha points starts from the lower bound and works it way up (so no 
    /// need to store the upper bound)
    double alphaLeftBound_;
    /// Size of steps between grid array values
    double deltaGrid_;
    /// Size of steps between alpha array values
    double deltaAlpha_;
    /// Scalar approximation parameter that determines 
    double h_;
    /// The last iterations values of the dynamic programming equation at each grid index
    std::vector<double>* oldV_;
    /// The current iterations values of the dynamic programming equation at each grid index
    std::vector<double>* newV_;
    /// The values for the optimal alpha (or control paramater) for each index along the grid is stored in this vector
    std::vector<double>* minAlpha_;
    /// \brief Smallest relative error before markov approximation exits.
    ///
    /// The markov chain continues to repeat until its relative error (norm_inf |oldV - newV|) is less than this value.
    /// The smaller this value, the higher the precision, but the greater number of iterations it will take.
    double epsErr_;
    /// \brief Total number of iterations before markov approximation exits.
    ///
    /// If the max number of iterations has been reached and the relative error still isn't less than epsErr_, then the
    /// markov approximation method will finish
    unsigned int maxIterations_;
};

