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
    GridParameters(const double* gridLeftBound,
                    const double* gridRightBound,
                    unsigned int* gridLength,
                    unsigned int dimensions);

	///
	/// \param gridLeftBound
	/// \param gridRightBound
	/// \param deltaGrid
	/// \param dimensions
    GridParameters(const double* gridLeftBound,
                    const double* gridRightBound,
                    const double* deltaGrid,
                    unsigned int dimensions);

	///
	/// \param gp
    GridParameters(const GridParameters& gp);

	///
    virtual ~GridParameters();

    // GETTERS ------------------------------------------------------------------------------------------------- GETTERS

    /// \brief Returns value of grid array at particular index
    ///
    /// Technically there is no array for the grid. Instead of storing the grid as another large array and obvious (in
    /// the sense that it's one increment after the other), this function just calculates and returns the grid at the
    /// specified index.
    /// \param index The index along the grid array. Must be >= 0 and less than the total length of the grid array.
    double getGridAtIndex(unsigned int index, unsigned int gridNum) const;

	///
	/// \return
    unsigned int getNumOfGrids() const;

	///
	/// \param dimension
	/// \return
    unsigned int getGridLength(unsigned int dimension) const;

	///
	/// \param gridIndex
	/// \return
    double getDeltaGrid(unsigned int gridIndex) const;

	///
	/// \param gridIndices
	/// \param gridLocOut
    void getGridAtIndices(const GridIndex& gridIndices, double* gridLocOut) const;
private:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

	///
	/// \param gridLeftBound
	/// \param gridRightBound
	/// \param gridLength
    void assertParameters(const double* gridLeftBound,
                          const double* gridRightBound,
                          const unsigned int* gridLength);

	///
	/// \param gridLeftBound
	/// \param gridRightBound
	/// \param deltaGrid
    void assertParameters(const double* gridLeftBound,
                          const double* gridRightBound,
                          const double* deltaGrid);

    // FIELDS --------------------------------------------------------------------------------------------------- FIELDS

    /// \brief Lower bound of the grid array.
    ///
    /// Calculating grid points starts from the lower bound and works it way up (so no
    /// need to store the upper bound)
    double* gridLeftBound_ = nullptr;
    ///
    unsigned int* gridLength_ = nullptr;
    /// Size of steps between grid array values
    double* deltaGrid_ = nullptr;
    ///
    unsigned int numOfGridDimensions_;
};


