//
// Created by david on 28/09/17.
//

#pragma once

#include <map>

// Foreword declaration
class GridParameters;

///
class GridIndex
{
public:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

	///
	/// \param numOfDimensions_
    explicit GridIndex(unsigned int numOfDimensions_);

    /// Copy constructor  (there is no = operator overload, to copy you must use the copy constructor)
    /// \param gi
    GridIndex(const GridIndex& gi);

	///
    ~GridIndex();

	///
	/// \return
    GridIndex& operator= (const GridIndex&);

	///
	/// \param padding
    void resetToOrigin(unsigned int padding = 0);

	///
	/// \param gp
	/// \param padding
	/// \return
    bool nextGridElement(const GridParameters& gp, unsigned int padding = 0);

	///
	/// \param gridNum
	/// \return
    unsigned int getIndexOfDim(unsigned int gridNum) const;

	///
	/// \param gridIndices
    void setGridIndices(const unsigned int* gridIndices);

	///
	/// \param gridNum
	/// \param index
    void setGridIndexOfDim(unsigned int gridNum, unsigned int index);

	///
	/// \return
    unsigned int id() const;
private:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

    ///
    /// \param padding
    /// \param curDimension
    /// \param gp
    /// \return
    bool recursionCount(unsigned int padding, unsigned int curDimension, const GridParameters& gp);

    // FIELDS --------------------------------------------------------------------------------------------------- FIELDS
    unsigned int numOfDimensions_;
    unsigned int* gridIndices_ = nullptr;

};

///
struct GridIndexCompare
{
    ///
    /// \param lhs
    /// \param rhs
    /// \return
    bool operator() (const GridIndex& lhs, const GridIndex& rhs) const
    {
        return lhs.id() < rhs.id();
    }
};

///
typedef std::map<GridIndex, double, GridIndexCompare> GridIndices;

