//
// Created by david on 24/08/17.
//
#pragma once
#include "mat.h"
#include "matrix.h"
#include <string>

///
class DotMatCreator
{
public:
    ///
    enum class FileStatus
    {
        NO_WRITE_PERMISSION,
        OUT_OF_MEMORY,
        UNDEFINED_ERROR,
        NO_ERROR
    };

    ///
    /// \param filename
    /// \param MATvarname
    explicit DotMatCreator(const char *filename, const char *MATvarname);

    ///
    ~DotMatCreator();

    ///
    /// \param arr
    /// \param arrLongLength
    void fillMatFile(double **arr, size_t arrLongLength);

    ///
    /// \return
    FileStatus getFilestatus();
private:
    ///
    MATFile* pmat_ = nullptr;
    ///
    const char* filename_;
    ///
    const char* MATVarname_;
    ///
    FileStatus fileStatus_ = FileStatus::NO_ERROR;
};
