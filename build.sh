#!/bin/bash

INSTALL_FOLDER=$PREFIX/share
BIN_FOLDER=$PREFIX/bin
mkdir -p $INSTALL_FOLDER
mkdir -p $INSTALL_FOLDER/resources/Data/HMMs/After_tcoffee_UPI
mkdir -p $INSTALL_FOLDER/resources/Data/FASTA/CDHIT
mkdir -p $BIN_FOLDER

cp -r workflow/scripts $INSTALL_FOLDER/
cp -r config $INSTALL_FOLDER/
cp -r resources/Data/HMMs/After_tcoffee_UPI/* $INSTALL_FOLDER/
cp -r resources/Data/FASTA/CDHIT/* $INSTALL_FOLDER/
cp -r results $INSTALL_FOLDER/
cp plastedma.py $INSTALL_FOLDER/
# ls -lt $INSTALL_FOLDER

cp $INSTALL_FOLDER/scripts/* $BIN_FOLDER
# for testing
cp resources/Data/FASTA/literature_seq/lit_sequences.fasta $INSTALL_FOLDER
# pre-built models for PE
# cp -r $INSTALL_FOLDER/resources/Data/HMMs/After_tcoffee_UPI/* $BIN_FOLDER/resources/Data/HMMs/After_tcoffee_UPI/
# sequences that give the models
# cp -r $INSTALL_FOLDER/resources/Data/FASTA/CDHIT/* $BIN_FOLDER/resources/Data/FASTA/CDHIT/
# default for negative control
cp resources/Data/FASTA/human_gut_metagenome.fasta $INSTALL_FOLDER/resources/Data/FASTA/

ls -l $INSTALL_FOLDER
#ls -l $INSTALL_FOLDER/resources/Data/FASTA/CDHIT
#ls -l $INSTALL_FOLDER/resources/Data/HMMs&/After_tcoffee_UPI
ls -l $BIN_FOLDER
chmod u+x $INSTALL_FOLDER/plastedma.py
ln -s $INSTALL_FOLDER/plastedma.py $BIN_FOLDER/plastedma
# chmod u+x $BIN_FOLDER/plastedma.py
