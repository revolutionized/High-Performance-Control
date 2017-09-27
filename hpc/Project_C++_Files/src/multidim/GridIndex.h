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

    bool nextGridElement(unsigned int padding, GridParameters& gp);

    unsigned int getIndexOfDim(unsigned int gridNum);

    void setGridIndices(const unsigned int* gridIndices);

    void setGridIndexOfDum(unsigned int gridNum, unsigned int index);

private:
    bool recursionCount(unsigned int padding, unsigned int curDimension, GridParameters& gp);
    unsigned int numOfDimensions_;
    unsigned int* gridIndices_ = nullptr;

};


