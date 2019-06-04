/*
 * Implementation of JffCorrector
 */

#include "JffCorrector.h"

/*
 * Default constructor
 */
JffCorrector::JffCorrector() :
  fFileLoaded(false),
  fSpilloverLoaded(false)
{
  
  // JFF correction histograms for jet shape
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iCentrality = 0; iCentrality < DijetHistogramManager::kMaxCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < DijetHistogramManager::kMaxTrackPtBins; iTrackPt++){
        fhJetShapeCorrection[iJetTrack][iCentrality][iTrackPt] = NULL;
        fhDeltaEtaDeltaPhiCorrection[iJetTrack][iCentrality][iTrackPt] = NULL;
        fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrack][iCentrality][iTrackPt] = NULL;
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation type loop
}

/*
 * Constructor
 */
JffCorrector::JffCorrector(TFile *inputFile)
{
  ReadInputFile(inputFile);
}

/*
 * Constructor
 */
JffCorrector::JffCorrector(TFile *inputFile, TFile *spilloverFile)
{
  ReadInputFile(inputFile);
  ReadSpilloverFile(spilloverFile);
}

/*
 * Copy constructor
 */
JffCorrector::JffCorrector(const JffCorrector& in) :
  fFileLoaded(in.fFileLoaded),
  fSpilloverLoaded(in.fSpilloverLoaded)
{
  // Copy constructor
  
  // JFF correction histograms for jet shape
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iCentrality = 0; iCentrality < DijetHistogramManager::kMaxCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < DijetHistogramManager::kMaxTrackPtBins; iTrackPt++){
        fhJetShapeCorrection[iJetTrack][iCentrality][iTrackPt] = in.fhJetShapeCorrection[iJetTrack][iCentrality][iTrackPt];
        fhDeltaEtaDeltaPhiCorrection[iJetTrack][iCentrality][iTrackPt] = in.fhDeltaEtaDeltaPhiCorrection[iJetTrack][iCentrality][iTrackPt];
        fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrack][iCentrality][iTrackPt] = in.fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrack][iCentrality][iTrackPt];
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation type loop
}

/*
 * Destructor
 */
JffCorrector::~JffCorrector(){
  
}

// Setter for input file
void JffCorrector::ReadInputFile(TFile *inputFile){
  
  // Create histogram manager to find correct histogram naming in the input file
  DijetHistogramManager *namerHelper = new DijetHistogramManager();
  DijetCard *card = new DijetCard(inputFile);
  
  // Load the correction histograms from the file
  char histogramNamer[200];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    //for(int iCentrality = 0; iCentrality < card->GetNCentralityBins(); iCentrality++){ // TODO: Uncomment after correction files have been reproduced
    for(int iCentrality = 0; iCentrality < 4; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < card->GetNTrackPtBins(); iTrackPt++){
        sprintf(histogramNamer,"%s_%s/jffCorrection_%s_%s_C%dT%d",namerHelper->GetJetShapeHistogramName(DijetHistogramManager::kJetShape), namerHelper->GetJetTrackHistogramName(iJetTrack),namerHelper->GetJetShapeHistogramName(DijetHistogramManager::kJetShape), namerHelper->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        fhJetShapeCorrection[iJetTrack][iCentrality][iTrackPt] = (TH1D*) inputFile->Get(histogramNamer);
        
        sprintf(histogramNamer,"%sDeltaEtaDeltaPhi/jffCorrection_%sDeltaEtaDeltaPhi_C%dT%d",namerHelper->GetJetTrackHistogramName(iJetTrack), namerHelper->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        fhDeltaEtaDeltaPhiCorrection[iJetTrack][iCentrality][iTrackPt] = (TH2D*) inputFile->Get(histogramNamer);
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation type loop
  
  // Raise the flag that input file has been loaded
  fFileLoaded = true;
  
  // Delete the helper objects
  delete card;
  delete namerHelper;
}

// Setter for spillover file
void JffCorrector::ReadSpilloverFile(TFile *spilloverFile){
  
  // Create histogram manager to find correct histogram naming in the input file
  DijetHistogramManager *namerHelper = new DijetHistogramManager();
  //DijetCard *card = new DijetCard(spilloverFile); TODO: Uncomment after spillover correction files have been reproduced
  
  // Load the correction histograms from the file
  char histogramNamer[200];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    //for(int iCentrality = 0; iCentrality < card->GetNCentralityBins(); iCentrality++){ // TODO: Uncomment after correction files have been reproduced
    for(int iCentrality = 0; iCentrality < 4; iCentrality++){
      //for(int iTrackPt = 0; iTrackPt < card->GetNTrackPtBins(); iTrackPt++){ // TODO: Uncomment after correction files have been reproduced
      for(int iTrackPt = 0; iTrackPt < 6; iTrackPt++){
        
        sprintf(histogramNamer,"%sDeltaEtaDeltaPhi/fittedSpilloverCorrection_%sDeltaEtaDeltaPhi_C%dT%d", namerHelper->GetJetTrackHistogramName(iJetTrack),namerHelper->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrack][iCentrality][iTrackPt] = (TH2D*) spilloverFile->Get(histogramNamer);
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation type loop
  
  // Raise the flag that input file has been loaded
  fSpilloverLoaded = true;
  
  // Delete the helper objects
  //delete card; // TODO: Uncomment after correction files have been reproduced
  delete namerHelper;
}

// Getter for JFF correction histograms for jet shape
TH1D* JffCorrector::GetJetShapeJffCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt) const{
  return fhJetShapeCorrection[iJetTrackCorrelation][iCentrality][iTrackPt];
}

// Getter for deltaEta-deltaPhi JFF correction histograms
TH2D* JffCorrector::GetDeltaEtaDeltaPhiJffCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt) const{
  return fhDeltaEtaDeltaPhiCorrection[iJetTrackCorrelation][iCentrality][iTrackPt];
}

// Getter for deltaEta-deltaPhi spillover correction histograms
TH2D* JffCorrector::GetDeltaEtaDeltaPhiSpilloverCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt) const{
  return fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrackCorrelation][iCentrality][iTrackPt];
}

// Return information, if correction is ready to be obtained
bool JffCorrector::CorrectionReady(){
  return fFileLoaded;
}

// Return information, if correction is ready to be obtained
bool JffCorrector::SpilloverReady(){
  return fSpilloverLoaded;
}
