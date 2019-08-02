/*
 * Implementation of JffCorrector
 */

#include "JffCorrector.h"

/*
 * Default constructor
 */
JffCorrector::JffCorrector() :
  fFileLoaded(false),
  fJffAsymmetryBins(0),
  fSpilloverLoaded(false),
  fSpilloverAsymmetryBins(0),
  fSystematicErrorLoaded(false),
  fSystematicAsymmetryBins(0)
{
  
  // JFF correction histograms for jet shape
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iCentrality = 0; iCentrality < DijetHistogramManager::kMaxCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < DijetHistogramManager::kMaxTrackPtBins; iTrackPt++){
        for(int iAsymmetry = 0; iAsymmetry <= DijetHistogramManager::kMaxAsymmetryBins; iAsymmetry++){
          fhJetShapeCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          fhDeltaEtaDeltaPhiCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
            fhJetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty] = NULL;
          } // Uncertainty source loop
        } // Asymmetry loop
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
  fJffAsymmetryBins(in.fJffAsymmetryBins),
  fSpilloverLoaded(in.fSpilloverLoaded),
  fSpilloverAsymmetryBins(in.fSpilloverAsymmetryBins),
  fSystematicErrorLoaded(in.fSystematicErrorLoaded),
  fSystematicAsymmetryBins(in.fSystematicAsymmetryBins)
{
  // Copy constructor
  
  // JFF correction histograms for jet shape
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iCentrality = 0; iCentrality < DijetHistogramManager::kMaxCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < DijetHistogramManager::kMaxTrackPtBins; iTrackPt++){
        for(int iAsymmetry = 0; iAsymmetry <= DijetHistogramManager::kMaxAsymmetryBins; iAsymmetry++){
          fhJetShapeCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = in.fhJetShapeCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt];
          fhDeltaEtaDeltaPhiCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = in.fhDeltaEtaDeltaPhiCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt];
          fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = in.fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt];
          for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
            fhJetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty] = in.fhJetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty];
          } // Uncertainty source loop
        } // Asymmetry loop
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
  fJffAsymmetryBins = card->GetNAsymmetryBins();
  
  // Load the correction histograms from the file
  char histogramNamer[200];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iCentrality = 0; iCentrality < card->GetNCentralityBins(); iCentrality++){
    //for(int iCentrality = 0; iCentrality < 4; iCentrality++){ // Needed for old files
      for(int iTrackPt = 0; iTrackPt < card->GetNTrackPtBins(); iTrackPt++){
        
        sprintf(histogramNamer, "%s_%s/jffCorrection_%s_%s_C%dT%d", namerHelper->GetJetShapeHistogramName(DijetHistogramManager::kJetShape), namerHelper->GetJetTrackHistogramName(iJetTrack), namerHelper->GetJetShapeHistogramName(DijetHistogramManager::kJetShape), namerHelper->GetJetTrackHistogramName(iJetTrack), iCentrality, iTrackPt);
        fhJetShapeCorrection[iJetTrack][fJffAsymmetryBins][iCentrality][iTrackPt] = (TH1D*) inputFile->Get(histogramNamer);
        
        sprintf(histogramNamer, "%sDeltaEtaDeltaPhi/jffCorrection_%sDeltaEtaDeltaPhi_C%dT%d", namerHelper->GetJetTrackHistogramName(iJetTrack), namerHelper->GetJetTrackHistogramName(iJetTrack), iCentrality, iTrackPt);
        fhDeltaEtaDeltaPhiCorrection[iJetTrack][fJffAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) inputFile->Get(histogramNamer);
        
        for(int iAsymmetry = 0; iAsymmetry < fJffAsymmetryBins; iAsymmetry++){
          
          sprintf(histogramNamer, "%s_%s/jffCorrection_%s_%s_A%dC%dT%d", namerHelper->GetJetShapeHistogramName(DijetHistogramManager::kJetShape), namerHelper->GetJetTrackHistogramName(iJetTrack), namerHelper->GetJetShapeHistogramName(DijetHistogramManager::kJetShape), namerHelper->GetJetTrackHistogramName(iJetTrack), iAsymmetry, iCentrality, iTrackPt);
          fhJetShapeCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH1D*) inputFile->Get(histogramNamer);
          
          sprintf(histogramNamer, "%sDeltaEtaDeltaPhi/jffCorrection_%sDeltaEtaDeltaPhi_A%dC%dT%d", namerHelper->GetJetTrackHistogramName(iJetTrack), namerHelper->GetJetTrackHistogramName(iJetTrack), iAsymmetry, iCentrality, iTrackPt);
          fhDeltaEtaDeltaPhiCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH2D*) inputFile->Get(histogramNamer);
          
        } // Asymmetry loop
        
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
  DijetCard *card = new DijetCard(spilloverFile);
  fSpilloverAsymmetryBins = card->GetNAsymmetryBins();
  
  // Load the correction histograms from the file
  char histogramNamer[200];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iCentrality = 0; iCentrality < card->GetNCentralityBins(); iCentrality++){
      for(int iTrackPt = 0; iTrackPt < card->GetNTrackPtBins(); iTrackPt++){
        
        sprintf(histogramNamer,"%sDeltaEtaDeltaPhi/nofitSpilloverCorrection_%sDeltaEtaDeltaPhi_C%dT%d", namerHelper->GetJetTrackHistogramName(iJetTrack),namerHelper->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        
        fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrack][fSpilloverAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) spilloverFile->Get(histogramNamer);
        
        for(int iAsymmetry = 0; iAsymmetry < fSpilloverAsymmetryBins; iAsymmetry++){
          
          sprintf(histogramNamer,"%sDeltaEtaDeltaPhi/nofitSpilloverCorrection_%sDeltaEtaDeltaPhi_A%dC%dT%d", namerHelper->GetJetTrackHistogramName(iJetTrack), namerHelper->GetJetTrackHistogramName(iJetTrack), iAsymmetry, iCentrality, iTrackPt);
          fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH2D*) spilloverFile->Get(histogramNamer);
          
        } // Asymmetry loop
        
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation type loop
  
  // Raise the flag that input file has been loaded
  fSpilloverLoaded = true;
  
  // Delete the helper objects
  delete card;
  delete namerHelper;
}

// Setter for spillover file
void JffCorrector::ReadSystematicFile(TFile *systematicFile){
  
  // Create histogram manager to find correct histogram naming in the input file
  DijetHistogramManager *namerHelper = new DijetHistogramManager();
  DijetCard *card = new DijetCard(systematicFile); //TODO: Uncomment after spillover correction files have been reproduced
  fSystematicAsymmetryBins = card->GetNAsymmetryBins();
  
  // Read the histograms from the file
  TString asymmetryString;
  TString histogramName;
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    
    for(int iAsymmetry = 0; iAsymmetry <= fSystematicAsymmetryBins; iAsymmetry++){
      
      // No asymmetry bins for inclusive jets
      if(iJetTrack >= DijetHistogramManager::kTrackInclusiveJet && iAsymmetry != fSystematicAsymmetryBins) continue;
      
      // Define a string for asymmetry
      if(iAsymmetry == fSystematicAsymmetryBins) {
        asymmetryString = "";
      } else {
        asymmetryString = Form("A%d",iAsymmetry);
      }
      
      for(int iCentrality = 0; iCentrality < card->GetNCentralityBins(); iCentrality++){
        for(int iTrackPt = 0; iTrackPt < card->GetNTrackPtBins(); iTrackPt++){
          for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
            histogramName = Form("%sUncertainty/jetShapeUncertainty_%s_%sC%dT%d_%s", namerHelper->GetJetTrackHistogramName(iJetTrack), namerHelper->GetJetTrackHistogramName(iJetTrack), asymmetryString.Data(), iCentrality, iTrackPt, uncertaintyName[iUncertainty].Data());
            fhJetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty] = (TH1D*) systematicFile->Get(histogramName.Data());
          } // Uncertainty source loop
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
  } // Jet-track correlation type loop
  
  // Raise the flag that input file has been loaded
  fSystematicErrorLoaded = true;
  
  // Delete the helper objects
  delete card;
  delete namerHelper;
  
}

// Read the histograms related to systematic uncertainties of long range correlations
void JffCorrector::ReadLongRangeSystematicFile(const char *systematicFile){
  
  // Create a stream to read the input file
  std::string lineInFile;
  std::ifstream systematicUncertainties(systematicFile);
  
  // The first line contains binning information
  std::getline(systematicUncertainties, lineInFile);
  int nCentralityBins, nFlowComponents, nTrackPtBins;
  
  std::istringstream lineStream(lineInFile);
  lineStream >> nCentralityBins;
  lineStream >> nFlowComponents;
  lineStream >> fSystematicAsymmetryBins;
  lineStream >> nTrackPtBins;
  
  // Loop over the file and read all the uncertainties to the master table
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    for(int iFlow = 0; iFlow < nFlowComponents; iFlow++){
      for(int iAsymmetry = 0; iAsymmetry <= fSystematicAsymmetryBins; iAsymmetry++){
        std::getline(systematicUncertainties, lineInFile);
        std::istringstream lineStream(lineInFile);
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          lineStream >> fLongRangeUncertaintyTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
        } // Track pT loop
      } // Asymmetry loop
    } // Flow component loop
  } // Centrality loop
  
}

// Getter for JFF correction histograms for jet shape
TH1D* JffCorrector::GetJetShapeJffCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated correction
  if(iAsymmetry < 0 || iAsymmetry >= fJffAsymmetryBins) iAsymmetry = fJffAsymmetryBins;
  
  return fhJetShapeCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt];
}

// Getter for deltaEta-deltaPhi JFF correction histograms
TH2D* JffCorrector::GetDeltaEtaDeltaPhiJffCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated correction
  if(iAsymmetry < 0 || iAsymmetry >= fJffAsymmetryBins) iAsymmetry = fJffAsymmetryBins;
  
  return fhDeltaEtaDeltaPhiCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt];
}

// Getter for deltaEta-deltaPhi spillover correction histograms
TH2D* JffCorrector::GetDeltaEtaDeltaPhiSpilloverCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated correction
  if(iAsymmetry < 0 || iAsymmetry >= fSpilloverAsymmetryBins) return fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrackCorrelation][fSpilloverAsymmetryBins][iCentrality][iTrackPt];
  
  // Return the correction in the selected bin
  return fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt];
}

// Getter systematic uncertainty histogram for jet shapes
TH1D* JffCorrector::GetJetShapeSystematicUncertainty(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry, int iUncertainty) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated uncertainty
  if(iAsymmetry < 0 || iAsymmetry >= fSystematicAsymmetryBins) iAsymmetry = fSystematicAsymmetryBins;
  
  // If uncertainty bin is outside of the uncertainty bin range, return total systematic uncertainty
  if(iUncertainty < 0 || iUncertainty > kTotal) iUncertainty = kTotal;
  
  // Return the uncertainty in the selected bin
  return fhJetShapeUncertainty[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt][iUncertainty];
}

// Getter for absolute systematic uncertainty for long range correlations
double JffCorrector::GetLongRangeSystematicUncertainty(const int iFlow, const int iCentrality, const int iTrackPt, int iAsymmetry) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated uncertainty
  if(iAsymmetry < 0 || iAsymmetry >= fSystematicAsymmetryBins) iAsymmetry = fSystematicAsymmetryBins;
  
  // Return the uncertainty in the selected bin
  return fLongRangeUncertaintyTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
}

// Getter for a name for the source of uncertainty
TString JffCorrector::GetUncertaintyName(const int iUncertainty) const{
  return uncertaintyName[iUncertainty];
}

// Getter for a name for the source of long range uncertainty
TString JffCorrector::GetLongRangeUncertaintyName(const int iUncertainty) const{
  return longRangeUncertaintyName[iUncertainty];
}

// Return information, if correction is ready to be obtained
bool JffCorrector::CorrectionReady(){
  return fFileLoaded;
}

// Return information, if correction is ready to be obtained
bool JffCorrector::SpilloverReady(){
  return fSpilloverLoaded;
}

// Return information, if correction is ready to be obtained
bool JffCorrector::SystematicsReady(){
  return fSystematicErrorLoaded;
}
