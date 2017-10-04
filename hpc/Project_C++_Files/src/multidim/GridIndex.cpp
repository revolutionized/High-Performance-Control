#include "GridIndex.h"
#include "GridParameters.h"

GridIndex::GridIndex(unsigned int numOfDimensions_)
        : numOfDimensions_(numOfDimensions_)
{
    gridIndices_ = new unsigned int[numOfDimensions_];
    resetToOrigin();
}

GridIndex::GridIndex(const GridIndex& gi)
    : numOfDimensions_(gi.numOfDimensions_)
{
    gridIndices_ = new unsigned int[numOfDimensions_];
    for (int ii = 0; ii < numOfDimensions_; ++ii)
    {
        gridIndices_[ii] = gi.gridIndices_[ii];
    }
}

GridIndex::~GridIndex()
{
    delete gridIndices_;
}

GridIndex& GridIndex::operator=(const GridIndex& rhs)
{
    if (this != &rhs)
    {
        GridIndex::numOfDimensions_ = rhs.numOfDimensions_;
        for (int ii = 0; ii < numOfDimensions_; ++ii)
        {
            GridIndex::gridIndices_[ii] = rhs.gridIndices_[ii];
        }
    }

    return *this;
}

void GridIndex::resetToOrigin(unsigned int padding)
{
    for (unsigned int ii = 0; ii < numOfDimensions_; ++ii)
    {
        gridIndices_[ii] = padding;
    }
}

bool GridIndex::nextGridElement(const GridParameters& gp, unsigned int padding)
{
    return recursionCount(padding, numOfDimensions_, gp);
}

bool GridIndex::recursionCount(unsigned int padding, unsigned int curDimension, const GridParameters& gp)
{
    if (curDimension > 0)
    {
        if (gridIndices_[curDimension-1] < gp.getGridLength(curDimension-1) - padding - 1)
        {
            gridIndices_[curDimension-1]++;
        }
        else
        {
            gridIndices_[curDimension-1] = padding;
            curDimension--;
            if (curDimension > 0)
            {
                return recursionCount(padding, curDimension, gp);
            }
            else
            {
                return false;
            }
        }
        return true;
    }
    return false;
}

unsigned int GridIndex::getIndexOfDim(unsigned int gridNum) const
{
    return gridIndices_[gridNum];
}

void GridIndex::setGridIndices(const unsigned int* gridIndices)
{
    for (int ii = 0; ii < numOfDimensions_; ++ii)
    {
        gridIndices_[ii] = gridIndices[ii];
    }
}

void GridIndex::setGridIndexOfDim(unsigned int gridNum, unsigned int index)
{
    gridIndices_[gridNum] = index;
}

unsigned int GridIndex::id() const
{
    unsigned int full_index = 0;
    for (int ii = static_cast<int>(numOfDimensions_) - 1; ii >= 0; --ii)
    {
        full_index += (numOfDimensions_ - ii)*getIndexOfDim(static_cast<unsigned int>(ii));
    }
    return full_index;
}

