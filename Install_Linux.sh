#!/bin/bash 
# -- Bash script
YELLOW='\033[1;33m'
NC='\033[0m' # No color
echo -e "${YELLOW}~~~ INSTALLING COMPONENT DEPENDENCIES FOR HIGH-PERFORMANCE-CONTROL ~~~ ${NC}"
#. Install_Files/Install_c3.sh
#. Install_Files/Install_cdyn.sh
. Install_Files/Install_c3sc.sh

echo -e "${YELLOW}~~~ INSTALLING HIGH-PERFORMANCE-CONTROL ~~~ ${NC}"
#. Install_Files/Install_hpc.sh

echo -e "${YELLOW}~~~ FINISHED ~~~ ${NC}"