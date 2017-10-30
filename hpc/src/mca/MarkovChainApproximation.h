//
// Created by david on 23/09/17.
//

#pragma once

#include <functional>
#include <vector>
#include <map>
#include <fstream>
#include "MarkovChainParameters.h"
#include "multidim/GridIndex.h"

typedef const std::function<void(double*, double, double*)> v_fcn_pddpd;
typedef const std::function<double(double*, double)> d_fcn_pdd;
typedef const std::function<void(double*, double**)> v_fcn_pdppd;

// region Old Code
/*
typedef const std::function<void(double*, double*)> v_fcn_pdpd;
*/
// endregion

/// This is the computational part of the Markov Chain Approximation method. It relies on having a pre-set
/// MarkovChainParameters to know what range of values to test against.
class MarkovChainApproximation
{
public:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

    /// Must initialise the Markov chain approximation with parameters, and an initial guess.
    /// \param mcp Parameters for this specific Markov Chain Approximation. Any changes to these parameters outside
    /// of this class will not be reflected in this class (because the MarkovChainParameters is basically copied).
    /// \param initStateGuess Initial guess is the initial guess of the dynamic programming equation. If you find
    /// that the process is converging slowly or diverging, consider adjusting this guess value)
    /// \param precision The precision is how many decimals are shown from the output window. The only output is the
    /// relative error between iterations.
    /// \param memoryModeRAM If RAM mode is used, then all variables are saved in memory. If working on a complex
    /// problem this may not be viable and instead the large memory components are written to a file. This will
    /// obviously be slower but will take advantage of the fact that there is more space on a drive than there is in
    /// memory.
    explicit MarkovChainApproximation(MarkovChainParameters& mcp,
                                      double initStateGuess,
                                      unsigned int precision = 4,
                                      bool memoryModeRAM = false);

    /// \brief Deallocates all memory that this class needed
    ~MarkovChainApproximation();

    /// \brief Main function to call when ready to compute the function approximation using the MCA method.
    /// \param costFunction The cost function is given as a std::function. It returns a double, and takes:
	/// Parameter (in order)[type] | Description
	/// -------------------------- | -----------
	/// x [double*]				   | Current grid value / location. Eg. x[0]=x, x[1]=y, and so on
	/// alpha [double]			   | Control value
	/// \param driftFunction The drift function is given as a std::function and takes:
	/// Parameter (in order)[type] | Description
	/// -------------------------- | -----------
	/// x [double*]				   | Current grid value / location. Eg. x[0]=x, x[1]=y, and so on
	/// alpha [double]			   | Control value
	/// out [double*]			   | This will be filled with the result (and hence should be the same size as x)
	/// \param diffMatrix The diffusion matrix is the given as a std::function and takes:
	/// Parameter (in order)[type] | Description
	/// -------------------------- | -----------
	/// x [double*]				   | Current grid value / location. Eg. x[0]=x, x[1]=y, and so on
	/// out [double**]			   | This will be filled with the diffusion matrix (and hence should be the same size
	/// as x in both directions - i.e. if x is of length n, then out should of size n x n)
	/// Note that this does not take a control (as the current Markov Chain Approximation is not set up to take a
	/// control parameter in the diffusion component).
    void computeMarkovApproximation(d_fcn_pdd& costFunction,
                                    v_fcn_pddpd& driftFunction,
									v_fcn_pdppd& diffMatrix);

    /// \brief Finds the most optimal control value at a specific location along arbitrary one-dimensional grid.
    ///
    /// \pre Must use computeMarkovApproximation prior to using this function or your return values will be incorrect
    /// \param gridLocation The location along the arbitrary grid that you wish to find the most optimal control at this
    /// point. This given x location can be anywhere, but for proper results it should exist within the bounds given
    /// at the construction of this class.
    /// \return The value of the optimal control parameter at the given x location. Or rather, at the location closest
    /// to x.
    double getMarkovControlFunction(double* gridLocation);

	/// Convenience function if you just want to view the control values for each grid
	/// \param stream The stream to write the control values to
    void printAlpha(std::ofstream& stream);

private:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

	/// Calculates the covariance matrix which is just the diffMatrix * diffMatrix', where ' means transpose
	/// \param sigmaMatrix The diffusion matrix
	/// \param dimensions Number of dimensions, n. Both sigmaMatrix and out should be n x n in size
	/// \param out This will become the covariance matrix, so make sure its the same size as sigmaMatrix
	/// \return
	bool getCovariance(double** sigmaMatrix, double** out);

	/// Tidy-up function to make code more manageable. Solves the summation of a specific alpha (filling up
	/// MarkovChainApproximation::vSummed_)
	/// \param covariance Covariance matrix
	/// \param k Cost function evaluation (at this location)
	/// \param b Drift function evaluation (at this location)
	/// \param deltaT Interpolation interval
	/// \param Qhat Normalising coefficient
	/// \param gridIndex Current index
	/// \param ai Current alpha iteration
	void solveTransitionSummations(double** covariance,
								   double k,
								   const double* b,
								   const double* deltaT,
								   const double* Qhat,
								   const GridIndex& gridIndex,
								   int ai);

	// region Old Code
/*
    /// \brief Loops through each alpha and calculates dynamic programming equation summation for that particular alpha.
    /// \param v_probabilities This is the possible transition states that can be taken on (so its of length 3 for the
    /// one-dimensional case). We pass this in as a parameter instead of constructing it in the function to avoid the
    /// time needed to allocate and deallocate memory for this vector (instead it is already allocated prior and we just
    /// reuse the vector as many times as needed).
    /// \param x_index Current index along the grid
    /// \param v_summed This is where all the summations are stored
    /// \param costFunctionK Same cost function given at the computeMarkovApproximation call
    /// \param diffFunction Same sigma function given at the computeMarkovApproximation call
    void solveTransitionSummations(const GridIndex& gridIndices,
                                   d_fcn_pdd& costFunctionK,
                                   v_fcn_pddpd& driftFunction,
                                   v_fcn_pdpd& diffFunction);

    /// \brief B function with no control dependency
    /// See "Numerical Methods for Stochastic Control Problems in Continuous Time" (H. J. Kushner) pg. 99 where it gives
    /// B(x) = max_{\alpha\el U} \abs{b(x,\alpha)}.
    void B_func(v_fcn_pddpd& driftFunction,
                double* gridLocation,
                double alpha,
                double* out);

    /// \brief Computes the transition values. This is a convenience function for splitting the code into smaller chunks
    /// \param v_probabilities This is the vector that gets filled with the transition values (being the transition
    /// probability multiplied by the previous value of the dynamic programming equation at the position).
    /// \param currentGridIndex Current index of the grid
    /// \param alpha Current value of alpha
    /// \param den The denominator for $ \Delta t $ is the same here as it is for $ p^{h} $, so we pass the value of the
    /// denominator in to save re-computation over each iteration.
    /// \param diffFunction Same sigma function given at the computeMarkovApproximation call
    void
    solveTransitionProbabilities(const GridIndex& currentIndices,
                                 uint currentDimension,
                                 double alpha,
                                 double den,
                                 v_fcn_pddpd& driftFunction,
                                 v_fcn_pdpd& diffFunction,
                                 double currentNumerator,
                                 double currentBpart);

    /// \brief Returns the value of $ p^{h}(x, y | \alpha) $
    /// \param currentLocation Current value of grid location
    /// \param y Transition state (x + h or x - h)
    /// \param alpha Current value of alpha
    /// \param den Denominator $ \sigma(x)^2 + hB(x) $
    /// \param b_part Same sigma function given at the computeMarkovApproximation call
    double transitionProb(bool upperJump, double num, double den, double b_part);

    /// \brief Summation of double array elements using Kahan formula (which ensures rounding errors aren't lost)
    /// \param array The array whose elements are to be summed
    /// \param size_n Total size of the given array
    /// \return The total value of all the elements summed together
    double kahanSum(const std::vector<double>* array);
*/
	// endregion

    /// \brief Picks the alpha that minimises the dynamic programming equation at the specific grid location.
    ///
    /// This is a convenience function for splitting the code into smaller chunks.
    /// \param gridIndices The current index of the grid
    void determineMinimumAlpha(const GridIndex& gridIndices);

    /// \brief Returns relative error between two arrays.
    ///
    /// The relative error is calculated as the infinity norm of the absolute difference between the two arrays - thus
    /// the two arrays must be the same length.
    /// \param arrLength The length of the two arrays (both arrays should be length arrLength)
    double getRelativeError();

    /// \brief Searches and finds the index of the grid array that most closely matches with the value of x given
    /// \param gridLocation Can be anything, but for proper use should be between the lower and upper bounds of the grid
    GridIndex getGridIndicesClosestTo(double* gridLocation);

	/// Tidy-up function to make code more manageable. Checks if all file streams are open. Only applicable for
	/// non-memory usage mode, and should not be called if steams have not been initialised.
    bool allStreamsOpen();

	/// Returns the trace of an n x n matrix (the sum of the diagonals)
	/// \param covariance The n x n matrix
	/// \param dimensions The 'n' part
	/// \return The value of the summation of the diagonals
	double getTrace(double** covariance, int dimensions);

    // FIELDS --------------------------------------------------------------------------------------------------- FIELDS

	/// A simple struct which represents a file being "read" and "written" to
    struct file
    {
		/// The file stream
        std::fstream* stream = nullptr;
		/// The name of the file
        std::string filename;
		/// Destructor calls file::deallocate
        ~file()
        {
            deallocate();
        }
		/// Checks if the stream has been initialised, and if so, closes it and deallocates its memory usage (but
		/// does not delete the file)
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
	/// The file where the minimum (optimal) alpha / control value is stored - Only if not using memory mode
    file minAlphaFile_;
	/// The file where the old dynamic programming evaluation is stored - Only if not using memory mode
    file oldVFile_;
	/// The file where the new dynamic programming evaluations are stored - Only if not using memory mod
    file newVFile_;
	/// If true, will use RAM for any storage of large arrays / maps. If false, will write large arrays / maps to
	/// files (NOT CURRENTLY IMPLEMENTED)
    bool memoryModeRAM_ = false;
	/// The parameters that define where the Markov Chain Approximation should look (e.g. the range of values)
    MarkovChainParameters mcp_;
    /// The last iterations values of the dynamic programming equation at each grid index
    GridIndices* oldV_ = nullptr;
    /// The current iterations values of the dynamic programming equation at each grid index
    GridIndices* newV_ = nullptr;
    /// The values for the optimal alpha (or control parameter) for each index along the grid is stored in this vector
    GridIndices* minAlpha_ = nullptr;
	/// The summation of all transition probabilities for each different alpha. Basically the part before finding the
	/// minimum control.
    std::vector<std::vector<double>>* vSummed_ = nullptr;

	// region Old Code
	/*
	///
    std::vector<double>* pHat_ = nullptr;
	///
    std::vector<double>* vProbs_ = nullptr;
    */
	// endregion
};

