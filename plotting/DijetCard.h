#ifndef DIJETCARD_H
#define DIJETCARD_H

// Root includes
#include <TFile.h>
#include <TString.h>
#include <TVectorT.h>

/*
 * Implementation of the DijetCard class
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
  
  void FindDataTypeString(){  // Construct a data type string based on information on the card
    
    // Define the different data types corresponding to certain indices
    TString dataTypes[5] = {"pp","PbPb","pp MC","PbPb MC","localTest"};
    if(fDataType < 0 || fDataType > 4){
      fDataTypeString = "Unknown";
      return;
    }
    
    // Define Monte Carlo types and add them to MC productions, which are data types 2 and 3
    TString monteCarloString[4] = {" RecoReco"," RecoGen"," GenReco"," GenGen"};
    if(fMonteCarloType < 0 || fMonteCarloType > 3){
      fDataTypeString = "Unknown";
      return;
    }
    
    for(int iDataType = 2; iDataType <=3; iDataType++){
      dataTypes[iDataType].Append(monteCarloString[fMonteCarloType]);
    }
    
    // Remember the constructed data type string
    fDataTypeString = dataTypes[fDataType];
  }
  
public:
  
  /*
   * Contructor with input file
   *
   *  TFile *inFile = Input file
   */
  DijetCard(TFile *inFile):
    fInputFile(inFile),
    fCardDirectory("JCard"),
    fDataType(-1),
    fDataTypeString("")
  {
    fInputFile->cd(fCardDirectory.Data());
    TVectorT<float> *reader = (TVectorT<float>*) gDirectory->Get("DataType");
    fDataType = (*reader)[1];
    TVectorT<float> *monteCarloType = (TVectorT<float>*) gDirectory->Get("McCorrelationType");
    fMonteCarloType = (*monteCarloType)[1];
    FindDataTypeString();
  }
  
  /*
   *  Getter for data type string
   */
  TString GetDataType() const{
    return fDataTypeString;
  }
  
};

#endif
