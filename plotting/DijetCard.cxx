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
  fDataType = (*fDataTypeVector)[1];
  
  // Read the Monte Carlo type from the vector: 0 = RecoReco, 1 = RecoGen, 2 = GenReco, 3 = GenGen
  if(fMcCorrelationTypeVector) {
    fMonteCarloType = (*fMcCorrelationTypeVector)[1];
  } else {
    fMonteCarloType = 0;
  }
  
  FindDataTypeString();
  
  // Read the asymmetry bin type from the vector: 0 = AJ, 1 = xJ
  if(fAsymmetryBinTypeVector){
    fAsymmetryBinType = (*fAsymmetryBinTypeVector)[1];
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
  fDataTypeVector = (TVectorT<float>*) gDirectory->Get("DataType");
  fMcCorrelationTypeVector = (TVectorT<float>*) gDirectory->Get("McCorrelationType");
  fMatchJetsVector = (TVectorT<float>*) gDirectory->Get("MatchJets");
  fForestTypeVector = (TVectorT<float>*) gDirectory->Get("ForestType");
  fReadModeVector = (TVectorT<float>*) gDirectory->Get("ReadMode");
  fJetTypeVector = (TVectorT<float>*) gDirectory->Get("JetType");
  fJetAxisVector = (TVectorT<float>*) gDirectory->Get("JetAxis");
  fJetEtaCutVector = (TVectorT<float>*) gDirectory->Get("JetEtaCut");
  fSearchEtaCutVector = (TVectorT<float>*) gDirectory->Get("SearchEtaCut");
  fMaxPtCutVector = (TVectorT<float>*) gDirectory->Get("MaxPtCut");
  fMinPtCutVector = (TVectorT<float>*) gDirectory->Get("MinPtCut");
  fSubleadingPtCutVector = (TVectorT<float>*) gDirectory->Get("SubleadingPtCut");
  fDeltaPhiCutVector = (TVectorT<float>*) gDirectory->Get("DeltaPhiCut");
  fMinMaxTrackPtFractionVector = (TVectorT<float>*) gDirectory->Get("MinMaxTrackPtFraction");
  fMaxMaxTrackPtFractionVector = (TVectorT<float>*) gDirectory->Get("MaxMaxTrackPtFraction");
  fTrackEtaCutVector = (TVectorT<float>*) gDirectory->Get("TrackEtaCut");
  fMinTrackPtCutVector = (TVectorT<float>*) gDirectory->Get("MinTrackPtCut");
  fMaxTrackPtRelativeErrorVector = (TVectorT<float>*) gDirectory->Get("MaxTrackPtRelativeError");
  fVertexMaxDistanceVector = (TVectorT<float>*) gDirectory->Get("VertexMaxDistance");
  fCalorimeterSignalLimitPtVector = (TVectorT<float>*) gDirectory->Get("CalorimeterSignalLimitPt");
  fHighPtEtFractionVector = (TVectorT<float>*) gDirectory->Get("HighPtEtFraction");
  fChi2QualityCutVector = (TVectorT<float>*) gDirectory->Get("Chi2QualityCut");
  fMinimumTrackHitsVector = (TVectorT<float>*) gDirectory->Get("MinimumTrackHits");
  fSubeventCutVector = (TVectorT<float>*) gDirectory->Get("SubeventCut");
  fZVertexCutVector = (TVectorT<float>*) gDirectory->Get("ZVertexCut");
  fLowPtHatCutVector = (TVectorT<float>*) gDirectory->Get("LowPtHatCut");
  fHighPtHatCutVector = (TVectorT<float>*) gDirectory->Get("HighPtHatCut");
  fAsymmetryBinTypeVector = (TVectorT<float>*) gDirectory->Get("AsymmetryBinType");
  fCentralityBinEdgesVector = (TVectorT<float>*) gDirectory->Get("CentralityBinEdges");
  fTrackPtBinEdgesVector = (TVectorT<float>*) gDirectory->Get("TrackPtBinEdges");
  fAsymmetryBinEdgesVector = (TVectorT<float>*) gDirectory->Get("AsymmetryBinEdges");
  fPtHatBinEdgesVector = (TVectorT<float>*) gDirectory->Get("PtHatBinEdges");
  fDoEventMixingVector = (TVectorT<float>*) gDirectory->Get("DoEventMixing");
  fMixWithPoolVector = (TVectorT<float>*) gDirectory->Get("MixWithPool");
  fNMixedEventsPerDijetVector = (TVectorT<float>*) gDirectory->Get("NMixedEventsPerDijet");
  fVzToleranceVector = (TVectorT<float>*) gDirectory->Get("VzTolerance");
  fMixingVzBinWidthVector = (TVectorT<float>*) gDirectory->Get("MixingVzBinWidth");
  fMixingHiBinWidthVector = (TVectorT<float>*) gDirectory->Get("MixingHiBinWidth");
  fMixingPoolDepthVector = (TVectorT<float>*) gDirectory->Get("MixingPoolDepth");
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
  return (*fJetEtaCutVector)[1] + (*fTrackEtaCutVector)[1];
}

// Get the number of asymmetry bins
int DijetCard::GetNAsymmetryBins() const{
  return fAsymmetryBinEdgesVector->GetNoElements()-1;
}

// Get the low border of i:th asymmetry bin
double DijetCard::GetLowBinBorderAsymmetry(const int iBin) const{
  
  // Sanity check for the input
  if(iBin < 0) return -1;
  if(iBin >= GetNAsymmetryBins()) return -1;
  
  // Return the asked bin index
  return (*fAsymmetryBinEdgesVector)[iBin+1];
}

// Get the high border of i:th asymmetry bin
double DijetCard::GetHighBinBorderAsymmetry(const int iBin) const{
  
  // Sanity check for the input
  if(iBin < 0) return -1;
  if(iBin >= GetNAsymmetryBins()) return -1;
  
  // Return the asked bin index
  return (*fAsymmetryBinEdgesVector)[iBin+2];
}

 // Get a description of the used asymmetry bin type
const char* DijetCard::GetAsymmetryBinType(TString latexIt) const{
  const char *asymmetryNames[] = {"Aj","Xj"};
  const char *asymmetryLatex[] = {"A_{J}","x_{J}"};
  if(latexIt.Contains("tex"),TString::kIgnoreCase) return asymmetryLatex[fAsymmetryBinType];
  return asymmetryNames[fAsymmetryBinType];
}

/*
 * Reader for all the vectors from the input file
 */
void DijetCard::Write(TDirectory *file){
  
  // Create a directory to store the card parameters
  if(!file->GetDirectory("JCard")) file->mkdir("JCard");
  file->cd("JCard");
  
  // Write all the vectors to the file. Not all of these exist in older versions of cards, thus check if exists before writing.
  if(fDataTypeVector) fDataTypeVector->Write("DataType");
  if(fMcCorrelationTypeVector) fMcCorrelationTypeVector->Write("McCorrelationType");
  if(fMatchJetsVector) fMatchJetsVector->Write("MatchJets");
  if(fForestTypeVector) fForestTypeVector->Write("ForestType");
  if(fReadModeVector) fReadModeVector->Write("ReadMode");
  if(fJetTypeVector) fJetTypeVector->Write("JetType");
  if(fJetAxisVector) fJetAxisVector->Write("JetAxis");
  if(fJetEtaCutVector) fJetEtaCutVector->Write("JetEtaCut");
  if(fSearchEtaCutVector) fSearchEtaCutVector->Write("SearchEtaCut");
  if(fMaxPtCutVector) fMaxPtCutVector->Write("MaxPtCut");
  if(fMinPtCutVector) fMinPtCutVector->Write("MinPtCut");
  if(fSubleadingPtCutVector) fSubleadingPtCutVector->Write("SubleadingPtCut");
  if(fDeltaPhiCutVector) fDeltaPhiCutVector->Write("DeltaPhiCut");
  if(fMinMaxTrackPtFractionVector) fMinMaxTrackPtFractionVector->Write("MinMaxTrackPtFraction");
  if(fMaxMaxTrackPtFractionVector) fMaxMaxTrackPtFractionVector->Write("MaxMaxTrackPtFraction");
  if(fTrackEtaCutVector) fTrackEtaCutVector->Write("TrackEtaCut");
  if(fMinTrackPtCutVector) fMinTrackPtCutVector->Write("MinTrackPtCut");
  if(fMaxTrackPtRelativeErrorVector) fMaxTrackPtRelativeErrorVector->Write("MaxTrackPtRelativeError");
  if(fVertexMaxDistanceVector) fVertexMaxDistanceVector->Write("VertexMaxDistance");
  if(fCalorimeterSignalLimitPtVector) fCalorimeterSignalLimitPtVector->Write("CalorimeterSignalLimitPt");
  if(fHighPtEtFractionVector) fHighPtEtFractionVector->Write("HighPtEtFraction");
  if(fChi2QualityCutVector) fChi2QualityCutVector->Write("Chi2QualityCut");
  if(fMinimumTrackHitsVector) fMinimumTrackHitsVector->Write("MinimumTrackHits");
  if(fSubeventCutVector) fSubeventCutVector->Write("SubeventCut");
  if(fZVertexCutVector) fZVertexCutVector->Write("ZVertexCut");
  if(fLowPtHatCutVector) fLowPtHatCutVector->Write("LowPtHatCut");
  if(fHighPtHatCutVector) fHighPtHatCutVector->Write("HighPtHatCut");
  if(fAsymmetryBinTypeVector) fAsymmetryBinTypeVector->Write("AsymmetryBinType");
  if(fCentralityBinEdgesVector) fCentralityBinEdgesVector->Write("CentralityBinEdges");
  if(fTrackPtBinEdgesVector) fTrackPtBinEdgesVector->Write("TrackPtBinEdges");
  if(fAsymmetryBinEdgesVector) fAsymmetryBinEdgesVector->Write("AsymmetryBinEdges");
  if(fPtHatBinEdgesVector) fPtHatBinEdgesVector->Write("PtHatBinEdges");
  if(fDoEventMixingVector) fDoEventMixingVector->Write("DoEventMixing");
  if(fMixWithPoolVector) fMixWithPoolVector->Write("MixWithPool");
  if(fNMixedEventsPerDijetVector) fNMixedEventsPerDijetVector->Write("NMixedEventsPerDijet");
  if(fVzToleranceVector) fVzToleranceVector->Write("VzTolerance");
  if(fMixingVzBinWidthVector) fMixingVzBinWidthVector->Write("MixingVzBinWidth");
  if(fMixingHiBinWidthVector) fMixingHiBinWidthVector->Write("MixingHiBinWidth");
  if(fMixingPoolDepthVector) fMixingPoolDepthVector->Write("MixingPoolDepth");
  
  // Return back to the main directory
  file->cd("../");
}
