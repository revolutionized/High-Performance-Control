#!/bin/bash -- Bash script
echo "~~~ INSTALLING COMPONENT DEPENDENCIES FOR HIGH-PERFORMANCE-CONTROL ~~~"
echo " -------------------- Installing C3 Library"
git clone https://github.com/goroda/Compressed-Continuous-Computation.git c3
cd c3
mkdir build
cd build
cmake ..
make
cd ../..
echo " -------------------- Installing CDYN Library"
git clone https://github.com/goroda/cdyn.git cdyn
cd cdyn
mkdir build
cd build
cmake ..
make
cd ../..
echo " -------------------- Installing C3SC Library"
git clone https://github.com/goroda/c3sc.git c3sc
cd c3sc
mkdir build
cd build
cmake ..
make
cd ../..
echo "~~~ INSTALLING HIGH-PERFORMANCE-CONTROL ~~~"
cd hpc
mkdir build
cd build
cmake ../Project_C++_Files
make
echo "~~~ FINISHED ~~~"
