//
// Created by david on 28/09/17.
//

#pragma once

// Foreword declaration
class GridIndex;

/// Contains the needed parameters that define a grid. Rather than storing all the grid values in memory, this class
/// takes into account that the grids being considered are uniform and hence we can calculate a current grid value
/// based on the distance between positions, and the current index value.
class GridParameters
{
public:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

	/// Allows you to define a multi dimensional grid based on the size of the vectors of each grid
	/// \param gridLeftBound An array of values that signify where that current dimension or 'grid' starts
	/// \param gridRightBound An array of values that signify where that current dimension or 'grid' ends, and hence
	/// should be greater than its counterpart left bound
	/// \param gridLength The length of each vector representing the discretised grid
	/// \param dimensions The total number of dimensions (hence all arrays should have this many dimensions)
    explicit GridParameters(const double* gridLeftBound,
							const double* gridRightBound,
							unsigned int* gridLength,
							unsigned int dimensions);

	/// Allows you to define a multi dimensional grid based on the distances between steps of each grid
	/// \param gridLeftBound An array of values that signify where that current dimension or 'grid' starts
	/// \param gridRightBound An array of values that signify where that current dimension or 'grid' ends, and hence
	/// should be greater than its counterpart left bound
	/// \param deltaGrid The distance between steps of each vector representing the discretised grid
	/// \param dimensions The total number of dimensions (hence all arrays should have this many dimensions)
    GridParameters(const double* gridLeftBound,
                   const double* gridRightBound,
                   const double* deltaGrid,
                   unsigned int dimensions);

	/// Copy constructor - all values are copied from the supplied GridParameters
	/// \param gp
    GridParameters(const GridParameters& gp);

	///
    virtual ~GridParameters();

    // GETTERS ------------------------------------------------------------------------------------------------- GETTERS

    /// \brief Returns continuous value of a grid array at a particular index and along a certain dimension
    ///
    /// Technically there is no array for the grid. Instead of storing the grid as another large array and obvious (in
    /// the sense that it's one increment after the other), this function just calculates and returns the grid at the
    /// specified index.
    /// \param index The index along the grid array. Must be >= 0 and less than the total length of the grid array.
	/// \param gridNum The index dimension being considered
    double getGridAtIndex(unsigned int index, unsigned int gridNum) const;

	/// Returns number of dimensions that this class accounts for
    unsigned int getNumOfGrids() const;

	/// Gets the total length of the vector for a certain dimension
	/// \param dimension The index of the dimension being considered
    unsigned int getGridLength(unsigned int dimension) const;

	/// Gets the distance between discretised points of a certain dimension
	/// \param dimension The index of the dimension being considered
    double getDeltaGrid(unsigned int dimension) const;

	/// Fills an array with the continuous values of each dimension based on the given indices
	/// \param gridIndices Contains the current index of each dimension
	/// \param gridLocOut This is the array that will be filled with the continuous values of each dimension, and
	/// thus should be the same size as the number of dimensions
    void getGridAtIndices(const GridIndex& gridIndices, double* gridLocOut) const;
private:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

	/// Tidy-up code used to confirm values given were appropriate
    void assertParameters(const double* gridLeftBound,
                          const double* gridRightBound,
                          const unsigned int* gridLength);

	/// Tidy-up code used to confirm values given were appropriate
    void assertParameters(const double* gridLeftBound,
                          const double* gridRightBound,
                          const double* deltaGrid);

    // FIELDS --------------------------------------------------------------------------------------------------- FIELDS

    /// \brief Lower bound of the grid array.
    ///
    /// Calculating grid points starts from the lower bound and works it way up (so no
    /// need to store the upper bound)
    double* gridLeftBound_ = nullptr;
    /// The length of each vector for each dimensions
    unsigned int* gridLength_ = nullptr;
    /// Size of steps between grid array values
    double* deltaGrid_ = nullptr;
    /// Number of dimensions this class considers
    unsigned int numOfGridDimensions_;
};


