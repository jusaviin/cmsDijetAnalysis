#!/bin/bash
echo "This is job number $1"

# Read the script arguments in format name=value
# The line below only reads the value from the above format
CARD=${2#*=}
OUTPUT=${3#*=}

# Untar the input file list
tar xf input_files.tar.gz

# Unzip tar ball
tar -xvzf dijet5TeV.tar.gz

# Compile the code
make

# Run the code
./dijetAnalysis $1 $CARD $OUTPUT
