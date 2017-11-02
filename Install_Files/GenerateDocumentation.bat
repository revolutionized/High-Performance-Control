@echo off
REM Uses doxygen to compile documentation.
REM IMPORTANT: Needs doxygen to be installed - else it will fail, it will all fail
REM Also needs to be run from the main directory (the one with README.md in it)
cd hpc\documentation
doxygen DoxyConfig
@echo on
echo 
echo ~~~ Succesfully generated code documentation. View it by clicking on hpc\documentation\html\index.html ~~~