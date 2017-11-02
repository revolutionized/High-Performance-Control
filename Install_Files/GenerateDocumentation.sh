# Uses doxygen to compile documentation.
# IMPORTANT: Needs doxygen to be installed - else it will fail, it will all fail
# Also needs to be run from the main directory (the one with README.md in it)
cd hpc/documentation
doxygen DoxyConfig
echo 
echo ~~~ Succesfully generated code documentation. View it by clicking on hpc\documentation\html\index.html ~~~
