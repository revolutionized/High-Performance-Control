//
// Created by david on 28/09/17.
//

#pragma once

#include "GridParameters.h"

class GridIndex
{
public:
    GridIndex(unsigned int numOfDimensions_);

    // Copy constrcutor
    GridIndex(const GridIndex& gi);

    ~GridIndex();

    void resetToOrigin(unsigned int padding = 0);

    bool nextGridElement(GridParameters& gp, unsigned int padding = 0);

    unsigned int getIndexOfDim(unsigned int gridNum);

    void setGridIndices(const unsigned int* gridIndices);

    void setGridIndexOfDim(unsigned int gridNum, unsigned int index);

private:
    bool recursionCount(unsigned int padding, unsigned int curDimension, GridParameters& gp);
    unsigned int numOfDimensions_;
    unsigned int* gridIndices_ = nullptr;

};


