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
  fJffTrackPtBins(0),
  fSpilloverLoaded(false),
  fSpilloverDeltaRLoaded(false),
  fSpilloverAsymmetryBins(0),
  fSpilloverTrackPtBins(0),
  fSystematicErrorLoaded(false),
  fSystematicAsymmetryBins(0),
  fSystematicTrackPtBins(0),
  fLongRangeAsymmetryBins(0),
  fTrackingCorrectionLoaded(false),
  fTrackingAsymmetryBins(0),
  fTrackingPtBins(0),
  fSmoothUncertainty(false),
  fSymmetrizeDeltaEta(true)
{
  
  // JFF correction histograms for jet shape
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iCentrality = 0; iCentrality < DijetHistogramManager::kMaxCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < DijetHistogramManager::kMaxTrackPtBins; iTrackPt++){
        for(int iAsymmetry = 0; iAsymmetry <= DijetHistogramManager::kMaxAsymmetryBins; iAsymmetry++){
          fhJetShapeCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          fhDeltaEtaDeltaPhiCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          fhJetShapeSpilloverCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          fhJetShapeSpilloverCorrectionManualTune[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          fhDeltaEtaDeltaPhiTrackingDeltaRCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          fhTrackDeltaRResidualScale[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
            fhJetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty] = NULL;
            fhDeltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty] = NULL;
          } // Uncertainty source loop
        } // Asymmetry loop
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation type loop
}

/*
 * Constructor
 */
JffCorrector::JffCorrector(TFile *inputFile) :
  JffCorrector()
{
  ReadInputFile(inputFile);
}

/*
 * Constructor
 */
JffCorrector::JffCorrector(TFile *inputFile, TFile *spilloverFile) :
  JffCorrector()
{
  ReadInputFile(inputFile);
  ReadSpilloverFile(spilloverFile);
}

/*
 * Constructor
 */
JffCorrector::JffCorrector(TFile *inputFile, TFile *spilloverFile, TFile *trackingFile) :
  JffCorrector()
{
  ReadInputFile(inputFile);
  ReadSpilloverFile(spilloverFile);
  ReadTrackDeltaRFile(trackingFile);
}

/*
 * Copy constructor
 */
JffCorrector::JffCorrector(const JffCorrector& in) :
  fFileLoaded(in.fFileLoaded),
  fJffAsymmetryBins(in.fJffAsymmetryBins),
  fJffTrackPtBins(in.fJffTrackPtBins),
  fSpilloverLoaded(in.fSpilloverLoaded),
  fSpilloverDeltaRLoaded(in.fSpilloverDeltaRLoaded),
  fSpilloverAsymmetryBins(in.fSpilloverAsymmetryBins),
  fSpilloverTrackPtBins(in.fSpilloverTrackPtBins),
  fSystematicErrorLoaded(in.fSystematicErrorLoaded),
  fSystematicAsymmetryBins(in.fSystematicAsymmetryBins),
  fSystematicTrackPtBins(in.fSystematicTrackPtBins),
  fLongRangeAsymmetryBins(in.fLongRangeAsymmetryBins),
  fTrackingCorrectionLoaded(in.fTrackingCorrectionLoaded),
  fTrackingAsymmetryBins(in.fTrackingAsymmetryBins),
  fTrackingPtBins(in.fTrackingPtBins),
  fSmoothUncertainty(in.fSmoothUncertainty),
  fSymmetrizeDeltaEta(in.fSymmetrizeDeltaEta)
{
  // Copy constructor
  
  // JFF correction histograms for jet shape
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iCentrality = 0; iCentrality < DijetHistogramManager::kMaxCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < DijetHistogramManager::kMaxTrackPtBins; iTrackPt++){
        for(int iAsymmetry = 0; iAsymmetry <= DijetHistogramManager::kMaxAsymmetryBins; iAsymmetry++){
          fhJetShapeCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = in.fhJetShapeCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt];
          fhDeltaEtaDeltaPhiCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = in.fhDeltaEtaDeltaPhiCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt];
          fhJetShapeSpilloverCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = in.fhJetShapeSpilloverCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt];
          fhJetShapeSpilloverCorrectionManualTune[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = in.fhJetShapeSpilloverCorrectionManualTune[iJetTrack][iAsymmetry][iCentrality][iTrackPt];
          fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = in.fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt];
          fhDeltaEtaDeltaPhiTrackingDeltaRCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = in.fhDeltaEtaDeltaPhiTrackingDeltaRCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt];
          fhTrackDeltaRResidualScale[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = in.fhTrackDeltaRResidualScale[iJetTrack][iAsymmetry][iCentrality][iTrackPt];
          for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
            fhJetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty] = in.fhJetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty];
            fhDeltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty] = in.fhDeltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty];
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

// Read the JFF correction file
void JffCorrector::ReadInputFile(TFile *inputFile){
  
  // Create histogram manager to find correct histogram naming in the input file
  DijetHistogramManager *namerHelper = new DijetHistogramManager();
  DijetCard *card = new DijetCard(inputFile);
  fJffAsymmetryBins = card->GetNAsymmetryBins();
  fJffTrackPtBins = card->GetNTrackPtBins();
  
  // Load the correction histograms from the file
  char histogramNamer[200];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iCentrality = 0; iCentrality < card->GetNCentralityBins(); iCentrality++){
      for(int iTrackPt = 0; iTrackPt < fJffTrackPtBins; iTrackPt++){
        
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

// Read the spillover file
void JffCorrector::ReadSpilloverFile(TFile *spilloverFile){
  
  // Create histogram manager to find correct histogram naming in the input file
  DijetHistogramManager *namerHelper = new DijetHistogramManager();
  DijetCard *card = new DijetCard(spilloverFile);
  fSpilloverAsymmetryBins = card->GetNAsymmetryBins();
  fSpilloverTrackPtBins = card->GetNTrackPtBins();
  
  // Load the correction histograms from the file
  char histogramNamer[200];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iCentrality = 0; iCentrality < card->GetNCentralityBins(); iCentrality++){
      for(int iTrackPt = 0; iTrackPt < fSpilloverTrackPtBins; iTrackPt++){
        
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

// Read the spillover file with histograms as a function of deltaR
void JffCorrector::ReadSpilloverDeltaRFile(TFile *spilloverFile){
  
  // Create histogram manager to find correct histogram naming in the input file
  DijetHistogramManager *namerHelper = new DijetHistogramManager();
  DijetCard *card = new DijetCard(spilloverFile);
  fSpilloverAsymmetryBins = card->GetNAsymmetryBins();
  fSpilloverTrackPtBins = card->GetNTrackPtBins();
  
  // Load the correction histograms from the file
  char histogramNamer[200];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iCentrality = 0; iCentrality < card->GetNCentralityBins(); iCentrality++){
      for(int iTrackPt = 0; iTrackPt < fSpilloverTrackPtBins; iTrackPt++){
        
        sprintf(histogramNamer,"spilloverDeltaR_%s/spilloverDeltaR_%s_C%dT%d", namerHelper->GetJetTrackHistogramName(iJetTrack),namerHelper->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        
        fhJetShapeSpilloverCorrection[iJetTrack][fSpilloverAsymmetryBins][iCentrality][iTrackPt] = (TH1D*) spilloverFile->Get(histogramNamer);
        
        sprintf(histogramNamer,"spilloverDeltaRManualTune_%s/spilloverDeltaRManualTune_%s_C%dT%d", namerHelper->GetJetTrackHistogramName(iJetTrack),namerHelper->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        
        fhJetShapeSpilloverCorrectionManualTune[iJetTrack][fSpilloverAsymmetryBins][iCentrality][iTrackPt] = (TH1D*) spilloverFile->Get(histogramNamer);
        
        for(int iAsymmetry = 0; iAsymmetry < fSpilloverAsymmetryBins; iAsymmetry++){
          
          sprintf(histogramNamer,"spilloverDeltaR_%s/spilloverDeltaR_%s_A%dC%dT%d", namerHelper->GetJetTrackHistogramName(iJetTrack), namerHelper->GetJetTrackHistogramName(iJetTrack), iAsymmetry, iCentrality, iTrackPt);
          fhJetShapeSpilloverCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH1D*) spilloverFile->Get(histogramNamer);
          
          sprintf(histogramNamer,"spilloverDeltaRManualTune_%s/spilloverDeltaRManualTune_%s_A%dC%dT%d", namerHelper->GetJetTrackHistogramName(iJetTrack), namerHelper->GetJetTrackHistogramName(iJetTrack), iAsymmetry, iCentrality, iTrackPt);
          fhJetShapeSpilloverCorrectionManualTune[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH1D*) spilloverFile->Get(histogramNamer);
          
        } // Asymmetry loop
        
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation type loop
  
  // Raise the flag that input file has been loaded
  fSpilloverDeltaRLoaded = true;
  
  // Delete the helper objects
  delete card;
  delete namerHelper;
}

// Read the residual tracking correction file
void JffCorrector::ReadTrackDeltaRFile(TFile *trackFile){
  
  // Create histogram manager to find correct histogram naming in the input file
  DijetHistogramManager *namerHelper = new DijetHistogramManager();
  DijetCard *card = new DijetCard(trackFile);
  fTrackingAsymmetryBins = card->GetNAsymmetryBins();
  fTrackingPtBins = card->GetNTrackPtBins();
  
  // Load the correction histograms from the file
  char histogramNamer[200];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    for(int iCentrality = 0; iCentrality < card->GetNCentralityBins(); iCentrality++){
      for(int iTrackPt = 0; iTrackPt < fTrackingPtBins; iTrackPt++){
        
        sprintf(histogramNamer,"%sDeltaEtaDeltaPhi/trackDeltaRCorrection_%sDeltaEtaDeltaPhi_C%dT%d", namerHelper->GetJetTrackHistogramName(iJetTrack), namerHelper->GetJetTrackHistogramName(iJetTrack), iCentrality, iTrackPt);
        
        fhDeltaEtaDeltaPhiTrackingDeltaRCorrection[iJetTrack][fTrackingAsymmetryBins][iCentrality][iTrackPt] = (TH2D*) trackFile->Get(histogramNamer);
        
        sprintf(histogramNamer,"%sResidualScale/residualScale_%s_C%dT%d", namerHelper->GetJetTrackHistogramName(iJetTrack), namerHelper->GetJetTrackHistogramName(iJetTrack), iCentrality, iTrackPt);
        
        fhTrackDeltaRResidualScale[iJetTrack][fTrackingAsymmetryBins][iCentrality][iTrackPt] = (TH1D*) trackFile->Get(histogramNamer);
        
        
        for(int iAsymmetry = 0; iAsymmetry < fTrackingAsymmetryBins; iAsymmetry++){
          
          sprintf(histogramNamer,"%sDeltaEtaDeltaPhi/trackDeltaRCorrection_%sDeltaEtaDeltaPhi_A%dC%dT%d", namerHelper->GetJetTrackHistogramName(iJetTrack), namerHelper->GetJetTrackHistogramName(iJetTrack), iAsymmetry, iCentrality, iTrackPt);
          fhDeltaEtaDeltaPhiTrackingDeltaRCorrection[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH2D*) trackFile->Get(histogramNamer);
          
          sprintf(histogramNamer,"%sResidualScale/residualScale_%s_A%dC%dT%d", namerHelper->GetJetTrackHistogramName(iJetTrack), namerHelper->GetJetTrackHistogramName(iJetTrack), iAsymmetry, iCentrality, iTrackPt);
          fhTrackDeltaRResidualScale[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH1D*) trackFile->Get(histogramNamer);
          
        } // Asymmetry loop
        
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation type loop
  
  // Raise the flag that input file has been loaded
  fTrackingCorrectionLoaded = true;
  
  // Delete the helper objects
  delete card;
  delete namerHelper;
}

// Read the systematic uncertainty file
void JffCorrector::ReadSystematicFile(TFile *systematicFile){
  
  // Create histogram manager to find correct histogram naming in the input file
  DijetHistogramManager *namerHelper = new DijetHistogramManager();
  DijetCard *card = new DijetCard(systematicFile);
  fSystematicAsymmetryBins = card->GetNAsymmetryBins();
  fSystematicTrackPtBins = card->GetNTrackPtBins();
  
  // Read the histograms from the file
  TString asymmetryString;
  TString histogramName;
  TH1D *jetShapeSum;
  TH1D *deltaEtaSum;
  double oldContent, addedContent, newContent;
  
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
        
        for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
          
          // Set the sums to NULL. If still NULL after loop, we know that we do not assign this to pT sum.
          jetShapeSum = NULL;
          deltaEtaSum = NULL;
          
          for(int iTrackPt = 0; iTrackPt < fSystematicTrackPtBins; iTrackPt++){
            
            histogramName = Form("%sUncertainty/jetShapeUncertainty_%s_%sC%dT%d_%s", namerHelper->GetJetTrackHistogramName(iJetTrack), namerHelper->GetJetTrackHistogramName(iJetTrack), asymmetryString.Data(), iCentrality, iTrackPt, uncertaintyName[iUncertainty].Data());
            fhJetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty] = (TH1D*) systematicFile->Get(histogramName.Data());
            
            histogramName = Form("%sUncertainty/deltaEtaUncertainty_%s_%sC%dT%d_%s", namerHelper->GetJetTrackHistogramName(iJetTrack), namerHelper->GetJetTrackHistogramName(iJetTrack), asymmetryString.Data(), iCentrality, iTrackPt, uncertaintyName[iUncertainty].Data());
              
            fhDeltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty] = (TH1D*) systematicFile->Get(histogramName.Data());
            
            // If the histogram does not exist, do not do the summing
            if(fhJetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty] == NULL) continue;
            
            // For pT summed histogram, in the first pT bin clone the histogram and add the histograms in quadrature in the other bins
            if(iTrackPt == 0){
              histogramName = Form("jetShapeUncertainty_%s_%sC%dT%d_%s", namerHelper->GetJetTrackHistogramName(iJetTrack), asymmetryString.Data(), iCentrality, fSystematicTrackPtBins, uncertaintyName[iUncertainty].Data());
              jetShapeSum = (TH1D*) fhJetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty]->Clone(histogramName);
              
              histogramName = Form("deltaEtaUncertainty_%s_%sC%dT%d_%s", namerHelper->GetJetTrackHistogramName(iJetTrack), asymmetryString.Data(), iCentrality, fSystematicTrackPtBins, uncertaintyName[iUncertainty].Data());
              deltaEtaSum = (TH1D*) fhDeltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty]->Clone(histogramName);
              
            } else {
              
              for(int iBin = 1; iBin <= jetShapeSum->GetNbinsX(); iBin++){
                oldContent = jetShapeSum->GetBinContent(iBin);
                addedContent = fhJetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty]->GetBinContent(iBin);
                
                // For tracking efficiency and residual tracking, add the uncertainty without quadratures
                if(iUncertainty == kTrackingEfficiency || iUncertainty == kResidualTracking){
                  newContent = oldContent+addedContent;
                } else {
                  newContent = TMath::Sqrt(oldContent*oldContent+addedContent*addedContent);
                }
                
                jetShapeSum->SetBinContent(iBin, newContent);
              }
              
              for(int iBin = 1; iBin <= deltaEtaSum->GetNbinsX(); iBin++){
                oldContent = deltaEtaSum->GetBinContent(iBin);
                addedContent = fhDeltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iUncertainty]->GetBinContent(iBin);
                
                // For tracking efficiency and residual tracking, add the uncertainty without quadratures
                if(iUncertainty == kTrackingEfficiency || iUncertainty == kResidualTracking){
                  newContent = oldContent+addedContent;
                } else {
                  newContent = TMath::Sqrt(oldContent*oldContent+addedContent*addedContent);
                }
                
                deltaEtaSum->SetBinContent(iBin, newContent);
              }
            }
            
          } // Track pT loop

          if(jetShapeSum) fhJetShapeUncertainty[iJetTrack][iAsymmetry][iCentrality][fSystematicTrackPtBins][iUncertainty] = (TH1D*) jetShapeSum->Clone();
          if(deltaEtaSum) fhDeltaEtaUncertainty[iJetTrack][iAsymmetry][iCentrality][fSystematicTrackPtBins][iUncertainty] = (TH1D*) deltaEtaSum->Clone();
          
        } // Uncertainty source loop
      } // Centrality loop
    } // Asymmetry loop
  } // Jet-track correlation type loop
  
  // Raise the flag that input file has been loaded
  fSystematicErrorLoaded = true;
  
  // Delete the helper objects
  delete card;
  delete namerHelper;
  
}

// Read the correction to v2 due to jet reconstruction bias
void JffCorrector::ReadJetReconstructionBiasFile(const char *fileName){
  
  // Create a stream to read the input file
  std::string lineInFile;
  std::ifstream fourierCorrections(fileName);
  
  // The first line contains binning information
  std::getline(fourierCorrections, lineInFile);
  int nCentralityBins, nFlowComponents, nTrackPtBins;
  
  std::istringstream lineStream(lineInFile);
  lineStream >> nCentralityBins;
  lineStream >> nFlowComponents;
  lineStream >> fLongRangeAsymmetryBins;
  lineStream >> nTrackPtBins;
  
  // Loop over the file and read all the uncertainties to the master table
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iFlow = 0; iFlow < nFlowComponents; iFlow++){
      for(int iAsymmetry = 0; iAsymmetry <= fLongRangeAsymmetryBins; iAsymmetry++){
        std::getline(fourierCorrections, lineInFile);
        std::istringstream lineStream(lineInFile);
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          lineStream >> fJetReconstructionBiasCorrection[iAsymmetry][iCentrality][iTrackPt][iFlow];
        } // Track pT loop
      } // Asymmetry loop
    } // Flow component loop
  } // Centrality loop
  
}

// Read the correction to v2 due to jet reconstruction bias
void JffCorrector::ReadJetReconstructionBiasFileForJetVn(const char *fileName){
  
  // Create a stream to read the input file
  std::string lineInFile;
  std::ifstream fourierCorrections(fileName);
  
  // The first line contains binning information
  std::getline(fourierCorrections, lineInFile);
  int nCentralityBins, nFlowComponents, nTrackPtBins;
  
  std::istringstream lineStream(lineInFile);
  lineStream >> nCentralityBins;
  lineStream >> nFlowComponents;
  lineStream >> fLongRangeAsymmetryBins;
  lineStream >> nTrackPtBins;
  
  // Loop over the file and read all the uncertainties to the master table
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iFlow = 0; iFlow < nFlowComponents; iFlow++){
      for(int iAsymmetry = 0; iAsymmetry <= fLongRangeAsymmetryBins; iAsymmetry++){
        std::getline(fourierCorrections, lineInFile);
        std::istringstream lineStream(lineInFile);
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          lineStream >> fJetReconstructionBiasCorrectionForJetVn[iAsymmetry][iCentrality][iTrackPt][iFlow];
        } // Track pT loop
      } // Asymmetry loop
    } // Flow component loop
  } // Centrality loop
  
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
  if(iAsymmetry < 0 || iAsymmetry > fJffAsymmetryBins) iAsymmetry = fJffAsymmetryBins;
  
  // If number of pT bins is given as an argument, return a sum of all pT bins
  if(iTrackPt == fJffTrackPtBins){
    TH2D *ptSumHistogram = (TH2D*) fhDeltaEtaDeltaPhiCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][0]->Clone(Form("%s_ptSum", fhDeltaEtaDeltaPhiCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][0]->GetName()));
    for(int iTrackPt = 1; iTrackPt < fJffTrackPtBins; iTrackPt++){
      ptSumHistogram->Add(fhDeltaEtaDeltaPhiCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt]);
    }
    return ptSumHistogram;
  }
  
  return fhDeltaEtaDeltaPhiCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt];
}

// Getter for deltaEta-deltaPhi spillover correction histograms
TH2D* JffCorrector::GetDeltaEtaDeltaPhiSpilloverCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated correction
  if(iAsymmetry < 0 || iAsymmetry > fSpilloverAsymmetryBins) iAsymmetry = fSpilloverAsymmetryBins;
  
  // If number of pT bins is given as an argument, return a sum of all pT bins
  if(iTrackPt == fSpilloverTrackPtBins){
    TH2D *ptSumHistogram = (TH2D*) fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][0]->Clone(Form("%s_ptSum", fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][0]->GetName()));
    for(int iTrackPt = 1; iTrackPt < fSpilloverTrackPtBins; iTrackPt++){
      ptSumHistogram->Add(fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt]);
    }
    return ptSumHistogram;
  }
  
  // Return the correction in the selected bin
  return fhDeltaEtaDeltaPhiSpilloverCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt];
}

// DeltaEta-DeltaPhi spillover correction histograms. Use scaled xj integrated distribution to calculate correction in different xj bins
TH2D* JffCorrector::GetDeltaEtaDeltaPhiSpilloverCorrectionAsymmetryScale(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, const int iAsymmetry, int usedBin) const{
  
  if(usedBin < 0) usedBin = fSpilloverAsymmetryBins;
  
  TH2D *correctionHistogram = (TH2D*) GetDeltaEtaDeltaPhiSpilloverCorrection(iJetTrackCorrelation, iCentrality, iTrackPt, usedBin)->Clone();
  
  // The spillover correction in different asymmetry bins is for a very good approximation a constant times xj integrated
  // This method can be used to suppress fluctuations for xj bins.
  double scale[] = {1.3, 1, 0.7, 1};
  correctionHistogram->Scale(scale[iAsymmetry]/scale[usedBin]);
  
  return correctionHistogram;
  
}

// Getter for spillover correction histograms as a function of deltaR
TH1D* JffCorrector::GetJetShapeSpilloverCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated correction
  if(iAsymmetry < 0 || iAsymmetry > fSpilloverAsymmetryBins) iAsymmetry = fSpilloverAsymmetryBins;
  
  // If number of pT bins is given as an argument, return a sum of all pT bins
  if(iTrackPt == fSpilloverTrackPtBins){
    TH1D *ptSumHistogram = (TH1D*) fhJetShapeSpilloverCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][0]->Clone(Form("%s_ptSum", fhJetShapeSpilloverCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][0]->GetName()));
    for(int iTrackPt = 1; iTrackPt < fSpilloverTrackPtBins; iTrackPt++){
      ptSumHistogram->Add(fhJetShapeSpilloverCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt]);
    }
    return ptSumHistogram;
  }
  
  // Return the correction in the selected bin
  return fhJetShapeSpilloverCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt];
}

// Getter for manually tuned spillover correction histograms as a function of deltaR
TH1D* JffCorrector::GetJetShapeSpilloverCorrectionManualTune(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated correction
  if(iAsymmetry < 0 || iAsymmetry > fSpilloverAsymmetryBins) iAsymmetry = fSpilloverAsymmetryBins;
  
  // If number of pT bins is given as an argument, return a sum of all pT bins
  if(iTrackPt == fSpilloverTrackPtBins){
    TH1D *ptSumHistogram = (TH1D*) fhJetShapeSpilloverCorrectionManualTune[iJetTrackCorrelation][iAsymmetry][iCentrality][0]->Clone(Form("%s_ptSum", fhJetShapeSpilloverCorrectionManualTune[iJetTrackCorrelation][iAsymmetry][iCentrality][0]->GetName()));
    for(int iTrackPt = 1; iTrackPt < fSpilloverTrackPtBins; iTrackPt++){
      ptSumHistogram->Add(fhJetShapeSpilloverCorrectionManualTune[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt]);
    }
    return ptSumHistogram;
  }
  
  // Return the correction in the selected bin
  return fhJetShapeSpilloverCorrectionManualTune[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt];
}

// Getter for deltaEta-deltaPhi spillover correction histograms
TH2D* JffCorrector::GetDeltaEtaDeltaPhiTrackDeltaRCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated correction
  if(iAsymmetry < 0 || iAsymmetry > fTrackingAsymmetryBins) iAsymmetry = fTrackingAsymmetryBins;
  
  // If number of pT bins is given as an argument, return a sum of all pT bins
  if(iTrackPt == fTrackingPtBins){
    TH2D *ptSumHistogram = (TH2D*) fhDeltaEtaDeltaPhiTrackingDeltaRCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][0]->Clone(Form("%s_ptSum", fhDeltaEtaDeltaPhiTrackingDeltaRCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][0]->GetName()));
    for(int iTrackPt = 1; iTrackPt < fTrackingPtBins; iTrackPt++){
      ptSumHistogram->Add(fhDeltaEtaDeltaPhiTrackingDeltaRCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt]);
    }
    return ptSumHistogram;
  }
  
  // Return the correction in the selected bin
  return fhDeltaEtaDeltaPhiTrackingDeltaRCorrection[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt];
}

// Getter for residual scale between reconstructed and generated tracks
double JffCorrector::GetTrackDeltaRResidualScale(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated correction
  if(iAsymmetry < 0 || iAsymmetry > fTrackingAsymmetryBins) iAsymmetry = fTrackingAsymmetryBins;
  
  // Old files do not have this histogram. Do not provide any scaling in that case
  if(fhTrackDeltaRResidualScale[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt] == NULL) return 1;
  
  return fhTrackDeltaRResidualScale[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(1);
  
}

// Getter systematic uncertainty histogram for jet shapes
TH1D* JffCorrector::GetJetShapeSystematicUncertainty(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry, int iUncertainty) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated uncertainty
  if(iAsymmetry < 0 || iAsymmetry >= fSystematicAsymmetryBins) iAsymmetry = fSystematicAsymmetryBins;
  
  // If uncertainty bin is outside of the uncertainty bin range, return total systematic uncertainty
  if(iUncertainty < 0 || iUncertainty > kTotal) iUncertainty = kTotal;
  
  if(iUncertainty != kPairAcceptance && iUncertainty != kBackgroundSubtraction && iUncertainty != kResidualTracking && iUncertainty != kTrackingEfficiency && fSmoothUncertainty){
    fhJetShapeUncertainty[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt][iUncertainty]->GetXaxis()->SetRangeUser(0.06,1);
    fhJetShapeUncertainty[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt][iUncertainty]->Smooth(1,"R");
    fhJetShapeUncertainty[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt][iUncertainty]->GetXaxis()->SetRangeUser(0,1);
  }
  
  // Return the uncertainty in the selected bin
  return fhJetShapeUncertainty[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt][iUncertainty];
}

// Getter systematic uncertainty histogram for deltaEta
TH1D* JffCorrector::GetDeltaEtaSystematicUncertainty(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry, int iUncertainty) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated uncertainty
  if(iAsymmetry < 0 || iAsymmetry >= fSystematicAsymmetryBins) iAsymmetry = fSystematicAsymmetryBins;
  
  // If uncertainty bin is outside of the uncertainty bin range, return total systematic uncertainty
  if(iUncertainty < 0 || iUncertainty > kTotal) iUncertainty = kTotal;
  
  // If symmetrization is set, symmetrize the uncertainty histogram before returning it
  if(fSymmetrizeDeltaEta){
    DijetMethods *symmetrizer = new DijetMethods();
    symmetrizer->SymmetrizeDeltaEta(fhDeltaEtaUncertainty[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt][iUncertainty]);
    delete symmetrizer;
  }
  
  // Return the uncertainty in the selected bin
  return fhDeltaEtaUncertainty[iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt][iUncertainty];
}

// Getter for jet reconstruction bias correction for jet-hadron Vn
double JffCorrector::GetJetReconstructionBiasCorrection(const int iFlow, const int iCentrality, const int iTrackPt, int iAsymmetry) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated uncertainty
  if(iAsymmetry < 0 || iAsymmetry >= fLongRangeAsymmetryBins) iAsymmetry = fLongRangeAsymmetryBins;
  
  // Return the uncertainty in the selected bin
  return fJetReconstructionBiasCorrection[iAsymmetry][iCentrality][iTrackPt][iFlow];
}

// Getter for jet reconstruction bias correction for jet vn
double JffCorrector::GetJetReconstructionBiasCorrectionForJetVn(const int iFlow, const int iCentrality, const int iTrackPt, int iAsymmetry) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated uncertainty
  if(iAsymmetry < 0 || iAsymmetry >= fLongRangeAsymmetryBins) iAsymmetry = fLongRangeAsymmetryBins;
  
  // Return the uncertainty in the selected bin
  return fJetReconstructionBiasCorrectionForJetVn[iAsymmetry][iCentrality][iTrackPt][iFlow];
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

// Return information, if JFF correction is ready to be obtained
bool JffCorrector::CorrectionReady(){
  return fFileLoaded;
}

// Return information, if spillover correction is ready to be obtained
bool JffCorrector::SpilloverReady(){
  return fSpilloverLoaded;
}

// Return information, if spillover correction as a function of deltaR is ready to be obtained
bool JffCorrector::SpilloverDeltaRReady(){
  return fSpilloverDeltaRLoaded;
}

// Return information, if residual tracking correction is ready to be obtained
bool JffCorrector::TrackingCorrectionReady(){
  return fTrackingCorrectionLoaded;
}

// Return information, if systematic uncertainties are ready to be obtained
bool JffCorrector::SystematicsReady(){
  return fSystematicErrorLoaded;
}

// Setter for smoothing the uncertainties
void JffCorrector::SetUncertaintySmooth(const bool smooth){
  fSmoothUncertainty = smooth;
}

// Setter for smoothing the uncertainties
void JffCorrector::SetDeltaEtaSymmetrization(const bool symmetrize){
  fSymmetrizeDeltaEta = symmetrize;
}
