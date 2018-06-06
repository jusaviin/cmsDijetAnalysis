#ifndef DIJETCARD_H
#define DIJETCARD_H

// C++ includes
#include <iostream>

// Root includes
#include <TFile.h>
#include <TString.h>
#include <TVectorT.h>

/*
 * DijetCard class
 *
 * This class reads the ConfigurationCard from the input root file and decodes the
 * necessary information for the analysis.
 */
class DijetCard {
  
private:
  
  TFile *fInputFile;         // Input file from which all the data is read
  TString fCardDirectory;    // Path to the ConfigurationCard directory
  int fDataType;             // Total number of centrality bins in the analysis
  int fMonteCarloType;       // Type of Monte Carlo used for jet-track correlations
  TString fDataTypeString;   // Total number of eta gaps in the analysis
  
  void FindDataTypeString(); // Construct a data type string based on information on the card
  
public:
  
  DijetCard(TFile *inFile); // Contructor with input file
  ~DijetCard();             // Destructor
  
  TString GetDataType() const; // Getter for data type string
  
};

#endif
