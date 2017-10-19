//
// Created by david on 24/08/17.
//

#include <iostream>
#include "DotMatCreator.h"
//#include <cstring>

using namespace std;

DotMatCreator::DotMatCreator(const char *filename, const char *MATvarname)
: filename_(filename),
  MATVarname_(MATvarname)
{
    cout << "Creating file " << filename_ << endl;
    pmat_ = matOpen(filename_, "w");

    if (pmat_ == nullptr)
    {
        cout << "error creating file " << filename_ << endl;
        cout << "(Do you have write permission in this directory?" << endl;
        fileStatus_ = FileStatus::NO_WRITE_PERMISSION;
    }
}

DotMatCreator::~DotMatCreator()
{
    // Clean up before exit
    if (matClose(pmat_) != 0)
    {
        cout << "Error closing file " << filename_ << endl;
    }
}

void DotMatCreator::fillMatFile(double **arr, size_t arrLongLength)
{
    // Create a 'MATLAB' styled array (reversed order since C stores opposite to matlab)
    mxArray* pa = mxCreateNumericMatrix(2, arrLongLength, mxDOUBLE_CLASS, mxREAL);
    if (pa == nullptr) {
        cout << __FILE__ << " : Out of memory on line " << __LINE__ << endl;
        cout << "Unable to create mxArray" << endl;
        fileStatus_ = FileStatus::OUT_OF_MEMORY;
    }

    // Copy solution to MATLAB matrix
    // (memcpy is faster but not working)
    double* mxPtr = mxGetPr(pa);
//    memcpy(mxGetData(pa), *(arr), sizeof(*(arr)));
    for (int i = 0; i < arrLongLength; ++i)
    {
        mxPtr[i] = arr[i][1];
        mxPtr[i+1] = arr[i][2];
    }

    // Now place this MATLAB array into the .mat file
    int status = matPutVariable(pmat_, MATVarname_, pa);
    if (status != 0) {
        cout << __FILE__ << " : Error using matPutVariable on line " << __LINE__ << endl;
        fileStatus_ = FileStatus::UNDEFINED_ERROR;
    }

    // Clear up memory
    mxDestroyArray(pa);
}

DotMatCreator::FileStatus DotMatCreator::getFilestatus()
{
    return fileStatus_;
}
