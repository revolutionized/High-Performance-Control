//
// Created by David on 20/10/2017.
//

#pragma once

#include <iostream>

namespace utils
{
    void printProgress(int progressComplete)
    {
        int barWidth = 40;
        std::cout << "[";
        int pos = int(barWidth * progressComplete/100.0);
        for (int ii = 0; ii < barWidth; ++ii) {
            if (ii < pos)
            {
                std::cout << "=";
            }
            else if (ii == pos)
            {
                std::cout << ">";
            }
            else
            {
                std::cout << " ";
            }
        }
        std::cout << "] " << "%" << progressComplete << "\r";
        std::cout.flush();
    }
}

