//

#include "GridIndex.h"

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

bool GridIndex::nextGridElement(unsigned int padding, GridParameters& gp)
{
    return recursionCount(padding, numOfDimensions_, gp);
}

bool GridIndex::recursionCount(unsigned int padding, unsigned int curDimension, GridParameters& gp)
{
    if (curDimension > 0)
    {
        if (gridIndices_[curDimension-1] < gp.getGridLength(curDimension) - padding - 1)
        {
            gridIndices_[curDimension-1]++;
        }
        else
        {
            gridIndices_[curDimension] = padding;
            curDimension--;
            recursionCount(padding, curDimension, gp);
        }
        return true;
    }
    return false;
}

unsigned int GridIndex::getIndexOfDim(unsigned int gridNum)
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

void GridIndex::setGridIndexOfDum(unsigned int gridNum, unsigned int index)
{
    gridIndices_[gridNum] = index;
}

