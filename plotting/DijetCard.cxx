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
  ReadVectors();
  fDataType = (*fCardEntries[kDataType])[1];
  
  // Read the Monte Carlo type from the vector: 0 = RecoReco, 1 = RecoGen, 2 = GenReco, 3 = GenGen
  if(fCardEntries[kMcCorrelationType]) {
    fMonteCarloType = (*fCardEntries[kMcCorrelationType])[1];
  } else {
    fMonteCarloType = 0;
  }
  
  FindDataTypeString();
  
  // Read the asymmetry bin type from the vector: 0 = AJ, 1 = xJ
  if(fCardEntries[kAsymmetryBinType]){
    fAsymmetryBinType = (*fCardEntries[kAsymmetryBinType])[1];
  } else {
    fAsymmetryBinType = 0;
  }
  
}

/*
 * Destructor
 */
DijetCard::~DijetCard(){

}

/*
 * Reader for all the vectors from the input file
 */
void DijetCard::ReadVectors(){
  
  // Read the TVectorT<float>:s
  for(int iEntry = 0; iEntry < knEntries; iEntry++){
    fCardEntries[iEntry] = (TVectorT<float>*) gDirectory->Get(fCardEntryNames[iEntry]);
  }
  
  // Read the file names
  for(int iFileName = 0; iFileName < knFileNames; iFileName++){
    fFileNames[iFileName] = (TObjString*) gDirectory->Get(fFileNameSaveName[iFileName]);
  }

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

/*
 *  Getter for maximum deltaEta
 */
double DijetCard::GetMaxDeltaEta() const{
  return (*fCardEntries[kJetEtaCut])[1] + (*fCardEntries[kTrackEtaCut])[1];
}

/*
 * Get the number of bins for internal index
 * If no vector is found in the index, return 0.
 */
int DijetCard::GetNBins(const int index) const{
  if(fCardEntries[index]) return fCardEntries[index]->GetNoElements()-1;
  return 0;
}

// Get the number of centrality bins
int DijetCard::GetNCentralityBins() const{
  return GetNBins(kCentralityBinEdges);
}

// Get the number of track pT bins
int DijetCard::GetNTrackPtBins() const{
  return GetNBins(kTrackPtBinEdges);
}

// Get the number of asymmetry bins
int DijetCard::GetNAsymmetryBins() const{
  return GetNBins(kAsymmetryBinEdges);
}

// Get the low border of i:th bin from internal index
double DijetCard::GetLowBinBorder(const int index, const int iBin) const{
  
  // Sanity check for the input
  if(iBin < 0) return -1;
  if(iBin > GetNBins(index)) return -1;
  
  // Return the asked bin index
  if(fCardEntries[index]) return (*fCardEntries[index])[iBin+1];
  return -1;
}

// Get the low border of i:th centrality bin
double DijetCard::GetLowBinBorderCentrality(const int iBin) const{
  return GetLowBinBorder(kCentralityBinEdges,iBin);
}

// Get the low border of i:th track pT bin
double DijetCard::GetLowBinBorderTrackPt(const int iBin) const{
  return GetLowBinBorder(kTrackPtBinEdges,iBin);
}

// Get the low border of i:th asymmetry bin
double DijetCard::GetLowBinBorderAsymmetry(const int iBin) const{
  return GetLowBinBorder(kAsymmetryBinEdges,iBin);
}

// Get the high border of i:th bin from internal index
double DijetCard::GetHighBinBorder(const int index, const int iBin) const{
  
  // Sanity check for the input
  if(iBin < 0) return -1;
  if(iBin >= GetNBins(index)) return -1;
  
  // Return the asked bin index
  if(fCardEntries[index]) return (*fCardEntries[index])[iBin+2];
  return -1;
}

// Get the high border of i:th centrality bin
double DijetCard::GetHighBinBorderCentrality(const int iBin) const{
  return GetHighBinBorder(kCentralityBinEdges,iBin);
}

// Get the high border of i:th track pT bin
double DijetCard::GetHighBinBorderTrackPt(const int iBin) const{
  return GetHighBinBorder(kTrackPtBinEdges,iBin);
}

// Get the high border of i:th asymmetry bin
double DijetCard::GetHighBinBorderAsymmetry(const int iBin) const{
  return GetHighBinBorder(kAsymmetryBinEdges,iBin);
}

 // Get a description of the used asymmetry bin type
const char* DijetCard::GetAsymmetryBinType(TString latexIt) const{
  const char *asymmetryNames[] = {"Aj","Xj"};
  const char *asymmetryLatex[] = {"A_{J}","x_{J}"};
  if(latexIt.Contains("tex"),TString::kIgnoreCase) return asymmetryLatex[fAsymmetryBinType];
  return asymmetryNames[fAsymmetryBinType];
}

/*
 * Get the number of deltaPhi bins bins
 * Note that this needs to be separate from the general GetNBins as low bin borders and high bin borders
 * are in different arrays, to there is a displacement of one bin with respect to the regular case.
 */
int DijetCard::GetNDeltaPhiBins() const{
  if(fCardEntries[kLowDeltaPhiBinBorders]) return fCardEntries[kLowDeltaPhiBinBorders]->GetNoElements();
  return 0;
}

// Get the low border of i:th deltaPhi bin
double DijetCard::GetLowBinBorderDeltaPhi(const int iBin) const{
  
  // Sanity check for the input
  if(iBin < 0) return -1;
  if(iBin >= GetNDeltaPhiBins()) return -1;
  
  // Return the asked bin index
  return (*fCardEntries[kLowDeltaPhiBinBorders])[iBin+1];
}

// Get the high border of i:th deltaPhi bin
double DijetCard::GetHighBinBorderDeltaPhi(const int iBin) const{
  
  // Sanity check for the input
  if(iBin < 0) return -1;
  if(iBin >= GetNDeltaPhiBins()) return -1;
  
  // Return the asked bin index
  return (*fCardEntries[kHighDeltaPhiBinBorders])[iBin+1];
}

/*
 * Add one-dimensional vector to the card
 *
 * Arguments:
 *  int entryIndex = Internal index for the added vector in entry array
 *  float entryContent = Content to be put into the vector
 */
void DijetCard::AddOneDimensionalVector(int entryIndex, float entryContent){
  
  // Only allow addition to postprocessing vectors
  if(entryIndex < kJffCorrection) return;
  
  // Make a new one dimensional vector to the desired index with given content
  float contents[1] = {entryContent};
  fCardEntries[entryIndex] = new TVectorT<float>(1,1,contents);

}

/*
 * Add one-dimensional vector to the card
 *
 * Arguments:
 *  int entryIndex = Internal index for the added vector in entry array
 *  int dimension = Number of entries in the given array
 *  float *contents = Content to be put into the vector
 */
void DijetCard::AddVector(int entryIndex, int dimension, double *contents){
  
  // Convert double pointer to float pointer
  float* convertedContents = new float[dimension];
  for(int i = 0; i < dimension; i++){
    convertedContents[i] = contents[i];
  }
  
  // Make a new vector to the desired index with given content
  fCardEntries[entryIndex] = new TVectorT<float>(1,dimension,convertedContents);
  
  // Delete the converted contents array
  delete[] convertedContents;
  
}

/*
 * Add file name to the card
 *
 * Arguments:
 *  int entryIndex = Internal index for the added file name in file name array
 *  TString fileName = Added file name
 */
void DijetCard::AddFileName(int entryIndex, TString fileName){
  
  // Make convert the string to TObjString and add it to the file name array
  fFileNames[entryIndex] = new TObjString(fileName.Data());
  
}

/*
 * Print the contents of the card to the console
 */
void DijetCard::Print() const{

  std::cout<<std::endl<<"========================= DijetCard =========================="<<std::endl;
  for(int iEntry = 0; iEntry < knEntries; iEntry++){
    if(fCardEntries[iEntry]){
      std::cout << Form("%25s",fCardEntryNames[iEntry]); //print keyword
      std::cout << " (dim = "<<fCardEntries[iEntry]->GetNoElements() << ") ";//print size of TVector
      for(int iElement = 1; iElement <= fCardEntries[iEntry]->GetNoElements(); iElement++){
        std::cout << (*fCardEntries[iEntry])[iElement] << " ";//TVector components
      }
      std::cout << std::endl;
    }
  }
  std::cout << std::endl;
  for(int iFileName = 0; iFileName < knFileNames; iFileName++){
    if(fFileNames[iFileName]){
      std::cout << "Used " << fFileNameType[iFileName] << " file: " << fFileNames[iFileName]->String().Data() << std::endl;
    }
  }
}

/*
 * Reader for all the vectors from the input file
 */
void DijetCard::Write(TDirectory *file){
  
  // Create a directory to store the card parameters
  if(!file->GetDirectory("JCard")) file->mkdir("JCard");
  file->cd("JCard");
  
  // Write all the vectors to the file. Not all of these exist in older versions of cards, thus check if exists before writing.
  for(int iEntry = 0; iEntry < knEntries; iEntry++){
    if(fCardEntries[iEntry])  fCardEntries[iEntry]->Write(fCardEntryNames[iEntry]);
  }
  
  // Write all the data names to the file.
  for(int iFileName = 0; iFileName < knFileNames; iFileName++){
    if(fFileNames[iFileName]) fFileNames[iFileName]->Write(fFileNameSaveName[iFileName]);
  }
   
  // Return back to the main directory
  file->cd("../");
}
