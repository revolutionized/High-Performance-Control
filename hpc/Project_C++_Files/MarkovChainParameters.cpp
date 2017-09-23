//
// Created by david on 23/09/17.
//

#include "MarkovChainParameters.h"
#include <iostream>
#include <typeinfo>
#include <cmath>

// CONSTRUCTORS + DESTRUCTOR ------------------------------------ CONSTRUCTORS + DESTRUCTOR //

MarkovChainParameters::MarkovChainParameters(double* gridLeftBound,
                                             double* gridRightBound,
                                             unsigned int* gridLength,
                                             double* alphaLeftBound,
                                             double* alphaRightBound,
                                             unsigned int* alphaLength,
                                             double h,)
    : h_(h)
{
    // Assert appropriate parameters have been provided
    assertParameters(gridLeftBound, gridRightBound, gridLength, alphaLeftBound, alphaRightBound, alphaLength);

    // Set up grid spacing and assign memory for own copies of lower bound and lengths
    unsigned int numOfGridDimensions = 0;
    while (gridLeftBound[numOfGridDimensions] != nullptr)
    {
        numOfGridDimensions++;
    };
    deltaGrid_ = new double[numOfGridDimensions];
    gridLeftBound_ = new double[numOfGridDimensions];
    gridLength_ = new unsigned int[numOfGridDimensions];
    for (int ii = 0; ii < numOfGridDimensions; ++ii)
    {
        deltaGrid_[ii] = (gridRightBound[ii] - gridLeftBound[ii]) / (gridLength[ii] - 1);
        gridLeftBound_[ii] = gridLeftBound[ii];
        gridLength_[ii] = gridLength[ii];
    }
    numOfAlphas_ = numOfGridDimensions;

    // Set up alpha discretisation  and assign memory for own copies of lower bound and lengths
    unsigned int numOfAlphaDimensions = 0;
    while (alphaLeftBound[numOfAlphaDimensions] != nullptr)
    {
        numOfAlphaDimensions++;
    };
    deltaAlpha_ = new double[numOfAlphaDimensions];
    alphaLeftBound_ = new double[numOfAlphaDimensions];
    alphaLength_ = new unsigned int[numOfAlphaDimensions];
    for (int ii = 0; ii < numOfAlphaDimensions; ++ii)
    {
        deltaAlpha_[ii] = (alphaRightBound[ii] - alphaLeftBound[ii]) / (alphaLength[ii] - 1);
        alphaLeftBound_[ii] = alphaLeftBound[ii];
        alphaLength_[ii] = alphaLength[ii];
    }
    numOfAlphas_ = numOfAlphaDimensions;
}

MarkovChainParameters::MarkovChainParameters(double* gridLeftBound,
                                             double* gridRightBound,
                                             double* deltaGrid,
                                             double* alphaLeftBound,
                                             double* alphaRightBound,
                                             double* deltaAlpha,
                                             double h)
    : h_(h)
{
    // Assert appropriate parameters have been provided
    assertParameters(gridLeftBound, gridRightBound, deltaGrid, alphaLeftBound, alphaRightBound, deltaAlpha);

    // Set up grid spacing and assign memory for own copies of lower bound and lengths
    unsigned int numOfGridDimensions = 0;
    while (gridLeftBound[numOfGridDimensions] != nullptr)
    {
        numOfGridDimensions++;
    };
    deltaGrid_ = new double[numOfGridDimensions];
    gridLeftBound_ = new double[numOfGridDimensions];
    gridLength_ = new unsigned int[numOfGridDimensions];
    for (int ii = 0; ii < numOfGridDimensions; ++ii)
    {
        deltaGrid_[ii] = deltaGrid[ii];
        gridLeftBound_[ii] = gridLeftBound[ii];
        gridLength_[ii] = static_cast<unsigned int>(floor((gridRightBound[ii] - gridLeftBound[ii])/deltaGrid[ii]));
    }

    // Set up alpha discretisation  and assign memory for own copies of lower bound and lengths
    unsigned int numOfAlphaDimensions = 0;
    while (alphaLeftBound[numOfAlphaDimensions] != nullptr)
    {
        numOfAlphaDimensions++;
    };
    deltaAlpha_ = new double[numOfAlphaDimensions];
    alphaLeftBound_ = new double[numOfAlphaDimensions];
    alphaLength_ = new unsigned int[numOfAlphaDimensions];
    for (int ii = 0; ii < numOfAlphaDimensions; ++ii)
    {
        deltaAlpha_[ii] = deltaAlpha[ii];
        alphaLeftBound_[ii] = alphaLeftBound[ii];
        alphaLength_[ii] = static_cast<unsigned int>(floor((alphaRightBound[ii] - alphaLeftBound[ii])/deltaAlpha[ii]));
    }

    // Set relative error to test against and max iterations
    epsErr_ = pow(10.0, -3.0);
    maxIterations_ = 100;
}

MarkovChainParameters::~MarkovChainParameters()
{
    delete gridLeftBound_;
    delete gridLength_;
    delete deltaGrid_;
    delete alphaLeftBound_;
    delete alphaLength_;
    delete deltaAlpha_;
}


// GETTERS ------------------------------------------------------------------------ GETTERS //

double MarkovChainParameters::getGridAtIndex(unsigned int index, unsigned int gridNum)
{
    assert(gridNum < numOfGrids_ && index < gridLength_[gridNum]);

    double gridPosition = gridLeftBound_[gridNum] + index * deltaGrid_[gridNum];
    return gridPosition;
}

double MarkovChainParameters::getAlphaAtIndex(unsigned int index, unsigned int alphaNum)
{
    assert(alphaNum < numOfAlphas_ && index < alphaLength_[alphaNum]);

    double alphaPosition = alphaLeftBound_[alphaNum] + index * deltaAlpha_[alphaNum];
    return alphaPosition;
}

unsigned int MarkovChainParameters::getNumOfGrids()
{
    return numOfGrids_;
}

unsigned int MarkovChainParameters::getNumOfAlphas()
{
    return numOfAlphas_;
}

unsigned int MarkovChainParameters::getGridLength(unsigned int gridNum)
{
    assert(gridNum < numOfGrids_);
    return gridLength_[gridNum];
}

unsigned int MarkovChainParameters::getAlphaLength(unsigned int alphaNum)
{
    assert(alphaNum < numOfAlphas_);
    return alphaLength_[alphaNum];
}


// PRIVATE METHODS -------------------------------------------------------- PRIVATE METHODS //

void MarkovChainParameters::assertParameters(double* gridLeftBound,
                                             double* gridRightBound,
                                             unsigned int* gridLength,
                                             double* alphaLeftBound,
                                             double* alphaRightBound,
                                             unsigned int* alphaLength)
{
    // First we check that both grid bound arrays and grid length are all of the same array size
    int counter = 0;
    do
    {
        if (gridLeftBound[counter] != nullptr
            && gridRightBound[counter] != nullptr
            && gridLength[counter] != nullptr)
        {
            counter++;
        }
        else if (gridLeftBound[counter] == nullptr
                && gridRightBound[counter] == nullptr
                && gridLength[counter] == nullptr)
        {
            break;
        }
        else
        {
            // An error here is thrown signifying that the gridLeftBound, gridRightBound and gridLength have differing
            // array sizes.
            std::cout << "ERROR: In " << typeid(this).name() << " on line " << __LINE__ << std::endl;
            assert(false);
        }
    } while(true);
    counter = 0;
    do
    {
        if (alphaLeftBound[counter] != nullptr
            && alphaRightBound[counter] != nullptr
            && alphaLength[counter] != nullptr)
        {
            counter++;
        }
        else if (alphaLeftBound[counter] == nullptr
                 && alphaRightBound[counter] == nullptr
                 && alphaLength[counter] == nullptr)
        {
            break;
        }
        else
        {
            // An error here is thrown signifying that the alphaLeftBound, alphaRightBound and alphaLength have
            // differing array sizes.
            std::cout << "ERROR: In " << typeid(this).name() << " on line " << __LINE__ << std::endl;
            assert(false);
        }
    } while(true);


    // Now we ensure appropriate bounds have been given
    while (alphaLeftBound[counter] != nullptr && alphaRightBound[counter] != nullptr)
    {
        assert(alphaLeftBound[counter] < alphaRightBound[counter]);
        counter++;
    };
    counter = 0;
    while (gridLeftBound[counter] != nullptr && gridRightBound[counter] != nullptr)
    {
        assert(gridLeftBound[counter] < gridRightBound[counter]);
        counter++;
    };

    // Now we ensure appropriate grid lengths and alpha lengths have been given
    counter = 0;
    while (gridLength[counter] != nullptr)
    {
        assert(gridLength[counter] > 1);
    };
    counter = 0;
    while (alphaLength[counter] != nullptr)
    {
        assert(alphaLength[counter] > 1);
    };
}

void MarkovChainParameters::assertParameters(double* gridLeftBound,
                                             double* gridRightBound,
                                             double* deltaGrid,
                                             double* alphaLeftBound,
                                             double* alphaRightBound,
                                             double* deltaAlpha)
{
    // First we check that both grid bound arrays and grid length are all of the same array size
    int counter = 0;
    do
    {
        if (gridLeftBound[counter] != nullptr
            && gridRightBound[counter] != nullptr)
        {
            counter++;
        }
        else if (gridLeftBound[counter] == nullptr
                 && gridRightBound[counter] == nullptr)
        {
            break;
        }
        else
        {
            // An error here is thrown signifying that the gridLeftBound and gridRightBound have differing
            // array sizes.
            std::cout << "ERROR: In " << typeid(this).name() << " on line " << __LINE__ << std::endl;
            assert(false);
        }
    } while(true);
    counter = 0;
    do
    {
        if (alphaLeftBound[counter] != nullptr
            && alphaRightBound[counter] != nullptr)
        {
            counter++;
        }
        else if (alphaLeftBound[counter] == nullptr
                 && alphaRightBound[counter] == nullptr)
        {
            break;
        }
        else
        {
            // An error here is thrown signifying that the alphaLeftBound and alphaRightBound have
            // differing array sizes.
            std::cout << "ERROR: In " << typeid(this).name() << " on line " << __LINE__ << std::endl;
            assert(false);
        }
    } while(true);


    // Now we ensure appropriate bounds have been given
    while (alphaLeftBound[counter] != nullptr && alphaRightBound[counter] != nullptr)
    {
        assert(alphaLeftBound[counter] < alphaRightBound[counter]);
        counter++;
    };
    counter = 0;
    while (gridLeftBound[counter] != nullptr && gridRightBound[counter] != nullptr)
    {
        assert(gridLeftBound[counter] < gridRightBound[counter]);
        counter++;
    };

    // Now we ensure appropriate grid deltas and alpha deltas have been given (can't be negative)
    counter = 0;
    while (deltaGrid[counter] != nullptr)
    {
        assert(deltaGrid[counter] > 0);
    };
    counter = 0;
    while (deltaAlpha[counter] != nullptr)
    {
        assert(deltaAlpha[counter] > 0);
    };
}



