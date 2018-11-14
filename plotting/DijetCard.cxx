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
  if(fMcCorrelationTypeVector) {
    fMonteCarloType = (*fMcCorrelationTypeVector)[1];
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
 * Reader for all the vectors from the input file
 */
void DijetCard::ReadVectors(){
  fDataTypeVector = (TVectorT<float>*) gDirectory->Get("DataType");
  fMcCorrelationTypeVector = (TVectorT<float>*) gDirectory->Get("McCorrelationType");
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
 * Reader for all the vectors from the input file
 */
void DijetCard::Write(TDirectory *file){
  
  // Create a directory to store the card parameters
  if(!file->GetDirectory("JCard")) file->mkdir("JCard");
  file->cd("JCard");
  
  // Write all the vectors to the file
  fDataTypeVector->Write("DataType");
  fMcCorrelationTypeVector->Write("McCorrelationType");
  fForestTypeVector->Write("ForestType");
  fReadModeVector->Write("ReadMode");
  fJetTypeVector->Write("JetType");
  fJetAxisVector->Write("JetAxis");
  fJetEtaCutVector->Write("JetEtaCut");
  fSearchEtaCutVector->Write("SearchEtaCut");
  fMaxPtCutVector->Write("MaxPtCut");
  fMinPtCutVector->Write("MinPtCut");
  fSubleadingPtCutVector->Write("SubleadingPtCut");
  fDeltaPhiCutVector->Write("DeltaPhiCut");
  fMinMaxTrackPtFractionVector->Write("MinMaxTrackPtFraction");
  fMaxMaxTrackPtFractionVector->Write("MaxMaxTrackPtFraction");
  fTrackEtaCutVector->Write("TrackEtaCut");
  fMinTrackPtCutVector->Write("MinTrackPtCut");
  fMaxTrackPtRelativeErrorVector->Write("MaxTrackPtRelativeError");
  fVertexMaxDistanceVector->Write("VertexMaxDistance");
  fCalorimeterSignalLimitPtVector->Write("CalorimeterSignalLimitPt");
  fHighPtEtFractionVector->Write("HighPtEtFraction");
  fChi2QualityCutVector->Write("Chi2QualityCut");
  fMinimumTrackHitsVector->Write("MinimumTrackHits");
  fSubeventCutVector->Write("SubeventCut");
  fZVertexCutVector->Write("ZVertexCut");
  fLowPtHatCutVector->Write("LowPtHatCut");
  fHighPtHatCutVector->Write("HighPtHatCut");
  fCentralityBinEdgesVector->Write("CentralityBinEdges");
  fTrackPtBinEdgesVector->Write("TrackPtBinEdges");
  fAsymmetryBinEdgesVector->Write("AsymmetryBinEdges");
  fPtHatBinEdgesVector->Write("PtHatBinEdges");
  fDoEventMixingVector->Write("DoEventMixing");
  fMixWithPoolVector->Write("MixWithPool");
  fNMixedEventsPerDijetVector->Write("NMixedEventsPerDijet");
  fVzToleranceVector->Write("VzTolerance");
  fMixingVzBinWidthVector->Write("MixingVzBinWidth");
  fMixingHiBinWidthVector->Write("MixingHiBinWidth");
  fMixingPoolDepthVector->Write("MixingPoolDepth");
  
  // Return back to the main directory
  file->cd("../");
}
