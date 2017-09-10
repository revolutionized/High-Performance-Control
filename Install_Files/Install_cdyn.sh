RED='\033[1;31m'
NC='\033[0m' # No color
echo -e " ${RED}-------------------- Installing CDYN Library ${NC}"
git clone https://github.com/goroda/cdyn.git cdyn
cd cdyn
mkdir build
cd build
cmake ..
make
make install
cd ../..