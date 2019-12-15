#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "JDrawer.h"

/*
 * Macro for printing the relative uncertainties in latex slides
 */
void printDeltaEtaUncertainties(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString uncertaintyFileName = "uncertainties/systematicUncertaintyForPbPb_15percentSpill20Jff_2019-10-01.root";
  TString ppUncertainryFileName = "uncertainties/systematicUncertaintyForPp_20percentSpillJff_2019-09-30.root";
  
  // Open the input files
  TFile *uncertaintyFile = TFile::Open(uncertaintyFileName);
  TFile *ppUncertainryFile = TFile::Open(ppUncertainryFileName);
  
  // Create histogram manager to manage the names for different jet collections
  DijetHistogramManager *namerHelper = new DijetHistogramManager();
  
  // Make readers to read the files
  JffCorrector *uncertaintyManager = new JffCorrector();
  uncertaintyManager->ReadSystematicFile(uncertaintyFile);
  
  JffCorrector *ppUncertaintyManager = new JffCorrector();
  ppUncertaintyManager->ReadSystematicFile(ppUncertainryFile);
  
  const int nCentralityBins = 4;
  const int nTrackPtBins = 7;
  const int nAsymmetryBins = 3;
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  TString asymmetryString[nAsymmetryBins+1];
  for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
    asymmetryString[iAsymmetry] = Form("%.1f < xj < %.1f",xjBinBorders[iAsymmetry],xjBinBorders[iAsymmetry+1]);
  }
  asymmetryString[nAsymmetryBins] = "xj integrated";
  
  // Define arrays for the jet shapes
  TH1D *uncertaintyAcceptanceCorrection[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  TH1D *uncertaintyBackgroundSubtraction[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = true;    // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = true;     // Produce the correction for inclusive jet-track correlatios
    
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,inclusiveJetTrack,inclusiveJetTrack};
  
  // Read PbPb and pp histograms from files
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          // Different data reader for pp and PbPb files
          if(iCentrality == nCentralityBins){
            uncertaintyAcceptanceCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = ppUncertaintyManager->GetDeltaEtaSystematicUncertainty(iJetTrack, 0, iTrackPt, iAsymmetry, JffCorrector::kPairAcceptance);
            uncertaintyBackgroundSubtraction[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = ppUncertaintyManager->GetDeltaEtaSystematicUncertainty(iJetTrack, 0, iTrackPt, iAsymmetry, JffCorrector::kBackgroundSubtraction);
          } else {
            uncertaintyAcceptanceCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = uncertaintyManager->GetDeltaEtaSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kPairAcceptance);
            uncertaintyBackgroundSubtraction[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = uncertaintyManager->GetDeltaEtaSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kBackgroundSubtraction);
          }
          
          //if(uncertaintyAcceptanceCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] == NULL){
          //  cout << "NULL!!!  iJetTrack: " << iJetTrack << " iAsymmetry: " << iAsymmetry << " iCentrality: " << iCentrality << " iTrackpT" << iTrackPt << endl;
          //}
          
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
  } // Jet-track correlation type loop
  
  // Print the uncertainties as a function of centrality and track pT
  char namer[100];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      if(iJetTrack >= DijetHistogramManager::kTrackInclusiveJet && iAsymmetry < nAsymmetryBins) continue;
      
      for(int iUncertainty = JffCorrector::kPairAcceptance; iUncertainty <= JffCorrector::kBackgroundSubtraction; iUncertainty++){
        
        
        cout << endl;
        cout << Form("Uncertainty for %s from %s in the bin %s", namerHelper->GetJetTrackHistogramName(iJetTrack), uncertaintyManager->GetUncertaintyName(iUncertainty).Data(), asymmetryString[iAsymmetry].Data()) << endl;

        cout << "pT C:0-10 C:10-30 C:30-50 C:50-90 pp " << endl;
        
        // Set the correct precision for printing floating point numbers
        cout << fixed << setprecision(5);
        
        // Print one line for each track pT bin
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          cout << Form("%.1f<pT<%.1f ", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
          for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
            if(iUncertainty == JffCorrector::kPairAcceptance){
              cout << " " << uncertaintyAcceptanceCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(1);
            }
            if(iUncertainty == JffCorrector::kBackgroundSubtraction){
              cout << " " << uncertaintyBackgroundSubtraction[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(1);
            }
          }
          cout << endl;
        }
        
        cout << endl;
      } // Uncertainty loop
    } // Asymmetry loop
  } // Jet-track correlation loop
}

