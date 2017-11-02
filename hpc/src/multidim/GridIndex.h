//
// Created by david on 28/09/17.
//

#pragma once

#include <map>

// Foreword declaration
class GridParameters;

/// The GridIndex is used to keep record of "current" grid location based on index values rather than continuous
/// values. It works hand-in-hand with the GridParameters when ever the continuous values are needed.
class GridIndex
{
public:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

	/// The GridIndex must be constructed with a specified number of dimensions and this number of dimensions will
	/// not change.
	/// \param numOfDimensions_
    explicit GridIndex(unsigned int numOfDimensions_);
	// todo: change so it takes GridParameters as input, and then nextGridElement should not have a need for it

    /// Copy constructor  (alteernatively you can use the = operator overload)
    /// \param gi The other GridIndex you wish to copy from
    GridIndex(const GridIndex& gi);

	/// Destructor - deallocates used memory
    ~GridIndex();

	/// Equals overload: Like a copy constructor, but requires the supplied GridIndex to be exist already
    GridIndex& operator= (const GridIndex&);

	/// Resets all indices to the value of padding - origin referring to the value of the "left bound" for all
	/// dimensions
	/// \param padding Value to reset all dimension indicies to (and hence is an unsigned int so that it is always
	/// greater than 0)
    void resetToOrigin(unsigned int padding = 0);

	/// Moves the Grid indices up along one dimension, if it goes beyond end - padding then it wraps back to the
	/// beginning and the next dimension goes up one. Much like a counter!
	/// \param gp The parameters that define length of each dimension
	/// \param padding This padding is applied to either side; so padding of 1 means that you "hit a wall" at end - 1
	/// index, and return to the 1 position and not the 0 position. This applies to all dimensions.
	/// \return False if the grid has reached it's maximum (i.e. it cannot step any further), else it returns true
	/// signifying that it successfully moved one step on some dimension.
    bool nextGridElement(const GridParameters& gp, unsigned int padding = 0);

	/// Returns the current index value of a certain dimension
	/// \param gridNum The index of the dimension being considered
    unsigned int getIndexOfDim(unsigned int gridNum) const;

	/// Changes the index value of all indices
	/// \param gridIndices Should be the same length as the number of dimensions, and all indices should not be
	/// greater than the number of elements in each vector (as set by the GridParameters being used in conjuction
	/// with this object)
    void setGridIndices(const unsigned int* gridIndices);

	/// Changes the index value of a single dimension
	/// \param gridNum The index of the dimension being considered
	/// \param index The value to change the index to
    void setGridIndexOfDim(unsigned int gridNum, unsigned int index);

	/// This id class is mainly intended to improve the ordering efficiency of GridIndices, thus when you loop
	/// through a GridIndices by moving along each element, it has the same (even less really) efficiency as calling
	/// nextGridElement.
	/// \return Basically adds up all dimensions by the weighting of their dimension value, e.g. if there were two
	/// dimensions and the current index values were [2,3] then id would return 2*2 + 1*3 = 7.
    unsigned int id() const;
private:
    // METHODS ------------------------------------------------------------------------------------------------- METHODS

    /// Really this is the nextGridElement function, so see the comments for GridIndex::nextGridElement
    bool recursionCount(unsigned int padding, unsigned int curDimension, const GridParameters& gp);

    // FIELDS --------------------------------------------------------------------------------------------------- FIELDS

	/// number of dimensions that this GridIndex accounts for (will become redundant in future releases)
    unsigned int numOfDimensions_;
	/// Array of the index values along each dimension (often called grids)
    unsigned int* gridIndices_ = nullptr;
};

/// This structure is used to properly order the GridIndices std::map, making it efficient to search through, and
/// making each element in the map act just like GridIndex::nextGridElement
struct GridIndexCompare
{
    /// Overloaded operator - see std::less if this doesn't make sense
    bool operator() (const GridIndex& lhs, const GridIndex& rhs) const
    {
        return lhs.id() < rhs.id();
    }
};

/// This map effectively defines the entire state space, but on discrete index values rather than continuous values
typedef std::map<GridIndex, double, GridIndexCompare> GridIndices;

