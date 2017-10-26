//
// Created by David on 20/10/2017.
//

#include <iostream>
// Include tclap for command line parsing
#include "tclap/CmdLine.h"
// Include examples files to check which user has selected to use
#include "dubinscar/dubinscar.h"
#include "onedimensional/OneDimensionPlus_Script.h"

// Want to hide cursor for console
#ifdef _WIN32
#include <windows.h>
void hidecursor(bool hide = true)
{
    HANDLE consoleHandle = GetStdHandle(STD_OUTPUT_HANDLE);
    CONSOLE_CURSOR_INFO info{};
    info.dwSize = 100;
    if (!hide)
    {
        info.bVisible = FALSE;
    }
    else
    {
        info.bVisible = TRUE;
    }
    SetConsoleCursorInfo(consoleHandle, &info);
}

#else

void hidecursor(bool hide = true)
{
    fputs(hide ? "\e[?25l" : "\e[?25h", stdout);
}

#endif

using std::cout;

int main(int argc, const char** argv)
{
    // Hide cursor
    hidecursor();

    // Wrap everything in a try block.  Do this every time,
    // because exceptions will be thrown for problems.
    try {

        // Definitioon of program
        TCLAP::CmdLine cmd("High Performance Control - Written by David Dos Santos, QUT", ' ', "1.1");

        // Define an argument for which example to run and add it to the command line.
        TCLAP::ValueArg<std::string> caseArg("c","case","Case / scenario / example to execute",
                                                true,"list","case_name");
        cmd.add(caseArg);

        // Define an argument for max iterations
        std::string dsc = "Max iterations to run Markov Chain Approximation for. Overrides any config file values.";
        TCLAP::ValueArg<int> maxItersArg("i", "iterations", dsc, false, 100, "unsigned int");
        cmd.add(maxItersArg);

        // Define a switch for utilising multi-threading TODO: NOT IMPLEMENTED YET
        TCLAP::SwitchArg parallelSwitchArg("p","parallel","Parallel computation / multithreading", cmd, false);

        // Parse the argv array.
        cmd.parse( argc, argv );

        // Parallel processing
        bool parallel = parallelSwitchArg.getValue();
        if ( parallel )
        {
            cout << "Parallel processing not implemented yet sorry... " << std::endl;
        }

        // Determine max iterations
        int testIterations = maxItersArg.getValue();
        if (testIterations <= 0)
        {
            cout << "Invalid iteration entered. Must be > 0. Using default of 100..." << std::endl;
            testIterations = 100;
        }
        auto maxIterations = static_cast<unsigned int>(testIterations);

        // Example to run
        std::string& example = caseArg.getValue();
        if (example == "list")
        {
            cout << "Please enter a valid case. Current cases/examples are: " << DUBINSCAR << "," << ONEDIMENSION;
            // Show cursor before exiting
            hidecursor(false);
            return 0;
        }
        else if (example == DUBINSCAR)
        {
            ExecuteDubinsCar(maxIterations);
        }
        else if (example == ONEDIMENSION)
        {
            ExecuteOneDimension(maxIterations);
        }


        // Finished with no errors ---------------------------------------------------------- //
        cout << " ~~~ Finished executing ~~~ ";
    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
        // Show cursor before exiting
        hidecursor(false);
        return 1;
    }

    // Show cursor before exiting
    hidecursor(false);
    return 0;
}