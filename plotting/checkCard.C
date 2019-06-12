#include "DijetCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)

/*
 * Macro for printing out DijetCard information from a given file
 */ 
 void checkCard(const char *fileName){
  TFile *file = TFile::Open(fileName);
  DijetCard *card = new DijetCard(file);
  card->Print();
  file->Close();
  delete card;
 }
