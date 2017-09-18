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
# Now need to add the installed library to a local directory since the cdyn install only goes to the bin/user path
case "$OSTYPE" in
  linux*)   mkdir ../lib && cp /usr/local/lib/cdyn/libcdyn.so ../lib/libcdyn.so ;;
  darwin*)  mkdir ../lib && echo "${RED}Un-implemented code - please copy the libcdyn.so from the OSX bin directory into cdyn/lib${NC}" ;; 
  cygwin*)  mkdir ../lib && echo "${RED}Un-implemented code - please copy the libcdyn.so from the Cygwin bin directory into cdyn/lib${NC}" ;;
esac
cd ../..
