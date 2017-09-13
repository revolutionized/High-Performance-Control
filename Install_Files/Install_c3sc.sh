# Assumes it is in the 
RED='\033[1;31m'
NC='\033[0m' # No color
echo -e " ${RED}-------------------- Installing C3SC Library ${NC}"
#git clone https://github.com/goroda/c3sc.git c3sc
cd c3sc
#mkdir build
cd build
#cmake ..

# There is a case in linux where it c3sc can't find the cdyn library, so I've created another CMakeLists to replace it
# with.
case "$OSTYPE" in
  linux*) . ../../Install_Files/Linux_replacement.sh;;
esac

make
make install
cd ../..