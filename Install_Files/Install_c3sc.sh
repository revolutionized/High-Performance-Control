RED='\033[1;31m'
NC='\033[0m' # No color
echo -e " ${RED}-------------------- Installing C3SC Library ${NC}"
git clone https://github.com/goroda/c3sc.git c3sc
cd c3sc
mkdir build
cd build
cmake ..
make
make install
cd ../..