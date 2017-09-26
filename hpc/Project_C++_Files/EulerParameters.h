//
// Created by David on 26/09/2017.
//

#pragma once


class EulerParameters
{
public:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

    EulerParameters(const double* gridLeftBound,
                    const double* gridRightBound,
                    unsigned int* gridLength,
                    unsigned int dimensions);

    EulerParameters(const double* gridLeftBound,
                    const double* gridRightBound,
                    const double* deltaGrid,
                    unsigned int dimensions);

    ~EulerParameters();

    // GETTERS ------------------------------------------------------------------------------------------------- GETTERS

    /// \brief Returns value of grid array at particular index
    ///
    /// Technically there is no array for the grid. Instead of storing the grid as another large array and obvious (in
    /// the sense that it's one increment after the other), this function just calculates and returns the grid at the
    /// specified index.
    /// \param index The index along the grid array. Must be >= 0 and less than the total length of the grid array.
    double getGridAtIndex(unsigned int index, unsigned int gridNum);

    void getGridAtIndex(unsigned int* index, double* outGrid);

    unsigned int getNumOfGrids();

    unsigned int getGridLength(unsigned int gridIndex);

    double getDeltaGrid(unsigned int gridIndex);

private:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

    void assertParameters(const double* gridLeftBound,
                          const double* gridRightBound,
                          unsigned int* gridLength);

    void assertParameters(const double* gridLeftBound,
                          const double* gridRightBound,
                          const double* deltaGrid);

    // FIELDS --------------------------------------------------------------------------------------------------- FIELDS

    double* gridLeftBound_ = nullptr;
    unsigned int* gridLength_ = nullptr;
    double* deltaGrid_ = nullptr;
    unsigned int numOfGridDimensions_;


};


