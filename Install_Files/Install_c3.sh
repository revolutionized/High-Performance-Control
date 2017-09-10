RED='\033[1;31m'
NC='\033[0m' # No color
echo -e " ${RED}-------------------- Installing C3 Library ${NC}"
git clone https://github.com/goroda/Compressed-Continuous-Computation.git c3
cd c3
mkdir build
cd build
cmake ..
make
make install
cd ../..