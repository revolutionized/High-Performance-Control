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
        auto test = gp.getGridLength(0);
        auto test1 = gp.getGridLength(1);
        auto test2 = gp.getGridLength(2);
        auto test3 = curDimension;
        auto test4 = gridIndices_[2];
        auto test5 = gridIndices_[1];
        auto test6 = gridIndices_[0];
        if (test6 > 17 && test5 > 23 && test4 > 23)
        {
            auto test7 = true;
        }
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

