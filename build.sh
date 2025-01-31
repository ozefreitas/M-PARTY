#!/bin/bash

set -eux 

INSTALL_FOLDER=$PREFIX/share
BIN_FOLDER=$PREFIX/bin
mkdir -p $INSTALL_FOLDER
mkdir -p $BIN_FOLDER/workflow

# cp -r workflow/scripts $BIN_FOLDER/
cp -r workflow $BIN_FOLDER/
# cp workflow/Snakefile $BIN_FOLDER/workflow
cp -r config $BIN_FOLDER/
# cp -r results $BIN_FOLDER/
cp m_party.py $BIN_FOLDER/
# ls -lt $INSTALL_FOLDER

echo "print bin folder"
ls -la $BIN_FOLDER

echo "print workflow folder"
ls -la $BIN_FOLDER/workflow

echo "print config folder"
ls -la $BIN_FOLDER/config

chmod u+x $BIN_FOLDER/m_party.py
ln -s $BIN_FOLDER/m_party.py $BIN_FOLDER/m-party
