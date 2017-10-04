//
// Created by david on 28/09/17.
//

#pragma once

#include <map>
//#include "GridParameters.h"

class GridParameters;

class GridIndex
{
public:
    GridIndex(unsigned int numOfDimensions_);

    // Copy constrcutor
    GridIndex(const GridIndex& gi);

    ~GridIndex();

    GridIndex& operator= (const GridIndex&);

    void resetToOrigin(unsigned int padding = 0);

    bool nextGridElement(const GridParameters& gp, unsigned int padding = 0);

    unsigned int getIndexOfDim(unsigned int gridNum) const;

    void setGridIndices(const unsigned int* gridIndices);

    void setGridIndexOfDim(unsigned int gridNum, unsigned int index);

    unsigned int id() const;
private:
    bool recursionCount(unsigned int padding, unsigned int curDimension, const GridParameters& gp);
    unsigned int numOfDimensions_;
    unsigned int* gridIndices_ = nullptr;

};

struct GridIndexCompare
{
    bool operator() (const GridIndex& lhs, const GridIndex& rhs) const
    {
        return lhs.id() < rhs.id();
    }
};

typedef std::map<GridIndex, double, GridIndexCompare> GridIndices;

