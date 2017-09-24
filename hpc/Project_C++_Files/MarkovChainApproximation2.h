//
// Created by david on 23/09/17.
//

#pragma once

#include <functional>
#include <vector>
#include <fstream>
#include "MarkovChainParameters.h"

class MarkovChainApproximation2
{
public:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS
    explicit MarkovChainApproximation2(MarkovChainParameters& mcp,
                                       double* initGuess,
                                       double* initAlphaGuess,
                                       unsigned int precision);

    /// \brief Deallocates all memory that this class needed
    ~MarkovChainApproximation2();

    /// \brief Main function to call when ready to compute the function approximation using the MCA method.
    /// \param costFunction The cost function is given as a std::function and takes two double arrays (x, alpha) and
    /// returns the value as a double.
    /// \param sigmaFunction The sigma function is often called the diffusion term and often represents represents the
    /// noise / stochastic disturbance being introduced to the system. It is given as a std::function and takes a double
    /// array representing the diffusion applied to each subsystem
    void computeMarkovApproximation(const std::function<double(double*, double*)>& costFunction,
                                    const std::function<double(double*)>& sigmaFunction);

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
    // METHODS ------------------------------------------------------------------------------------------------- METHODS



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

    /// \brief Loops through each alpha and calculates dynamic programming equation summation for that particular alpha.
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

    bool allStreamsOpen();

    // FIELDS --------------------------------------------------------------------------------------------------- FIELDS

    struct file
    {
        std::fstream* stream = nullptr;
        std::string filename;
        void deallocate()
        {
            if (stream != nullptr)
            {
                stream->close();
                delete stream;
                std::remove(filename.c_str());
            }
        }
    };
    std::vector<file> oldVFile_;
    std::vector<file> newVFile_;
    std::vector<file> minAlphaFile_;

    /// The last iterations values of the dynamic programming equation at each grid index
    std::vector<double>* oldV_ = nullptr;
    /// The current iterations values of the dynamic programming equation at each grid index
    std::vector<double>* newV_ = nullptr;
    /// The values for the optimal alpha (or control parameter) for each index along the grid is stored in this vector
    std::vector<double>* minAlpha_ = nullptr;
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

