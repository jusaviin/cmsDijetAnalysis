#ifndef DIJETCARD_H
#define DIJETCARD_H

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
  TString fDataTypeString;   // Total number of eta gaps in the analysis
  
  void FindDataTypeString(){
    TString dataTypes[4] = {"pp","PbPb","pp MC","PbPb MC"};
    if(fDataType < 0 || fDataType > 3){
      fDataTypeString = "Unknown";
      return;
    }
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
    TVector *reader = (TVector*) gDirectory->Get("DataType");
    fDataType = (*reader)[1];
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
