/*
 * Implementation of the DijetCard class
 */

// Own includes
#include "DijetCard.h"

/*
 * Contructor with input file
 *
 *  TFile *inFile = Input file
 */
DijetCard::DijetCard(TFile *inFile):
  fInputFile(inFile),
  fCardDirectory("JCard"),
  fDataType(-1),
  fDataTypeString("")
{
  fInputFile->cd(fCardDirectory.Data());
  TVectorT<float> *reader = (TVectorT<float>*) gDirectory->Get("DataType");
  fDataType = (*reader)[1];
  TVectorT<float> *monteCarloType = (TVectorT<float>*) gDirectory->Get("McCorrelationType");
  if(monteCarloType) {
    fMonteCarloType = (*monteCarloType)[1];
  } else {
    fMonteCarloType = 0;
  }
  FindDataTypeString();
}

/*
 * Destructor
 */
DijetCard::~DijetCard(){

}

/*
 * Construct a data type string based on information on the card
 */
void DijetCard::FindDataTypeString(){ 
  
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

/*
 *  Getter for data type string
 */
TString DijetCard::GetDataType() const{
  return fDataTypeString;
}

