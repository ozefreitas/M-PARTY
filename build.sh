#!/bin/bash

INSTALL_FOLDER=$PREFIX/share
BIN_FOLDER=$PREFIX/bin
mkdir -p $INSTALL_FOLDER
# mkdir -p $INSTALL_FOLDER/resources/Data/HMMs/PE/After_tcoffee_UPI
# mkdir -p $INSTALL_FOLDER/resources/Data/FASTA/PE/CDHIT
mkdir -p $BIN_FOLDER

cp -r workflow/scripts $BIN_FOLDER/
cp -r config $BIN_FOLDER/
cp -r results $BIN_FOLDER/
cp m-party.py $BIN_FOLDER/
# ls -lt $INSTALL_FOLDER

# cp $BIN_FOLDER/scripts/* $BIN_FOLDER
# for testing
# cp resources/Data/FASTA/human_gut_metagenome.fasta $INSTALL_FOLDER/resources/Data/FASTA/

ls -l $INSTALL_FOLDER
ls -l $BIN_FOLDER
chmod u+x $BIN_FOLDER/m-party.py
ln -s $BIN_FOLDER/m-party.py $BIN_FOLDER/m-party
