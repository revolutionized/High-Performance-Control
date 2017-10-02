//
// Created by david on 28/09/17.
//

#pragma once

class GridIndex;

class GridParameters
{
public:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

    GridParameters(const double* gridLeftBound,
                    const double* gridRightBound,
                    unsigned int* gridLength,
                    unsigned int dimensions);

    GridParameters(const double* gridLeftBound,
                    const double* gridRightBound,
                    const double* deltaGrid,
                    unsigned int dimensions);

    GridParameters(const GridParameters& gp);

    virtual ~GridParameters();

    // GETTERS ------------------------------------------------------------------------------------------------- GETTERS

    /// \brief Returns value of grid array at particular index
    ///
    /// Technically there is no array for the grid. Instead of storing the grid as another large array and obvious (in
    /// the sense that it's one increment after the other), this function just calculates and returns the grid at the
    /// specified index.
    /// \param index The index along the grid array. Must be >= 0 and less than the total length of the grid array.
    double getGridAtIndex(unsigned int index, unsigned int gridNum) const;

    unsigned int getNumOfGrids() const;

    unsigned int getGridLength(unsigned int gridIndex) const;

    double getDeltaGrid(unsigned int gridIndex) const;

    void getGridAtIndices(const GridIndex& gridIndices, double* gridLocOut) const;
private:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

    void assertParameters(const double* gridLeftBound,
                          const double* gridRightBound,
                          const unsigned int* gridLength);

    void assertParameters(const double* gridLeftBound,
                          const double* gridRightBound,
                          const double* deltaGrid);

    // FIELDS --------------------------------------------------------------------------------------------------- FIELDS

    /// \brief Lower bound of the grid array.
    ///
    /// Calculating grid points starts from the lower bound and works it way up (so no
    /// need to store the upper bound)
    double* gridLeftBound_ = nullptr;


    unsigned int* gridLength_ = nullptr;

    /// Size of steps between grid array values
    double* deltaGrid_ = nullptr;

    unsigned int numOfGridDimensions_;
};


