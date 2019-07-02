#!/bin/bash

# Root folder
PROGRAM_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo "  Program will be installed in $PROGRAM_ROOT/Install"

# Create install directory if it does not exist
if [ ! -d $PROGRAM_ROOT/Install ]
then
  mkdir -p $PROGRAM_ROOT/Install
  echo "  Install directory did not exist. Install directory has now been created."
else
  echo "  Install directory already exists. "
fi

# Run the CMake configuration script
echo -e "\n  Running CMake configuration script: "
cd $PROGRAM_ROOT/Build
./configure*
echo -e " \n !! If the CMake script compiled successfully go to the build folder and run make. \n"
