#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "DijetMethods.h"
#include "JDrawer.h"

/*
 * Macro for estimating systematic uncertainties for long range correlation results
 *
 * The following sources of systematic uncertainty are estimated:
 *  - Adjustment on the gluing point of leading/subleading distributions
 *  - Fit projections over regions -2.5 < deltaEta < -1.5 and 1.5 < deltaEta < 2.5 and compare to nominal
 *  - Fit projections over regions 1.5 < |deltaEta| < 2 and 2 < |deltaEta| < 2.5 and compare to nominal
 *  - Fit projections before and after mixed event correction
 *  - Fir different v_{z} selection (lika positive and negative)
 *
 *  The results will be saved to a file and also a slide with different contribution can be printed
 */
void estimateLongRangeSystematics(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Main files from which the long range asymmetries are obtained
  TString pbpbFileName = "data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_adjustedBackground_processed_2019-07-05.root";
  TString ppFileName = "data/dijet_pp_highForest_pfJets_noUncOrInc_allCorrections_adjustedBackground_wtaAxis_processed_2019-07-13.root";
  
  // For systematic uncertainty estimation, need files in which the background is not adjusted between leading and subleading side
  TString pbpbUnadjustedFileName = "data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root";
  TString ppUnadjustedFileName = "data/dijet_pp_highForest_pfJets_noUncOrInc_allCorrections_wtaAxis_processed_2019-07-13.root";
  
  const bool printUncertainties = true;
  
  const int nRefit = 4;    // Number of vn:s included in the refit
  int backgroundRebin = 4; // Rebin applied to the long range distributions
  
  // ==================================================================
  // ====================== Configuration done ========================
  // ==================================================================
  
  // Open the input files
  TFile *pbpbDataFile = TFile::Open(pbpbFileName);
  TFile *ppDataFile = TFile::Open(ppFileName);
  
  // Open the files needed for systematic uncertainty estimation
  TFile *pbpbUnadjustedFile = TFile::Open(pbpbUnadjustedFileName);
  TFile *ppUnadjustedFile = TFile::Open(ppUnadjustedFileName);
  
  // Create readers for the histograms and define binning
  DijetHistogramManager *pbpbReader = new DijetHistogramManager(pbpbDataFile);
  DijetHistogramManager *ppReader = new DijetHistogramManager(ppDataFile);
  DijetHistogramManager *pbpbUnadjustedReader = new DijetHistogramManager(pbpbUnadjustedFile);
  DijetHistogramManager *ppUnadjustedReader = new DijetHistogramManager(ppUnadjustedFile);
  
  const int nCentralityBins = pbpbReader->GetNCentralityBins();
  const int nTrackPtBins = pbpbReader->GetNTrackPtBins();
  const int nAsymmetryBins = pbpbReader->GetNAsymmetryBins();
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  // Define a fitter to test fitting the background distribution with different amount of vn:s
  DijetMethods *refitter = new DijetMethods();
  
  // Define arrays for the histograms
  TH1D *background[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  TF1 *backgroundFit[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  
  TH1D *sameEventLongRange[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  TF1 *sameEventFit[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  
  // Two dimensional histogram needed for projections over different regions
  TString typeString[2] = {"SameEvent","Corrected"};
  TH2D *deltaPhiDeltaEtaDistribution[2][2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins]; // First bin: leading/subleading. Second bin: same/corrected
  
  // Custom projections from different regions
  TH1D *positiveEtaProjection[2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];  // First bin: same/corrected
  TH1D *negativeEtaProjection[2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];  // First bin: same/corrected
  TH1D *closeEtaProjection[2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];     // First bin: same/corrected
  TH1D *farEtaProjection[2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];       // First bin: same/corrected
  
  // Fits to projections from different regions
  TF1 *positiveEtaFit[2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];  // First bin: same/corrected
  TF1 *negativeEtaFit[2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];  // First bin: same/corrected
  TF1 *closeEtaFit[2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];     // First bin: same/corrected
  TF1 *farEtaFit[2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];       // First bin: same/corrected
  
  // Define arrays needed for systematic uncertainty estimation
  TH1D *leadingUnadjustedBackground[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  TH1D *leadingUnadjustedBackgroundOverlap[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  
  // We will fit also the unadjusted histogram to estimate systematic uncertainties
  TF1 *unadjustedFit[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  
  // Define arrays for extracted vn numbers
  double masterFlowTable[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double masterFlowError[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // Read the numbers also for unadjusted fits to estimate systematic uncertainties
  double unadjustedFlowTable[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double unadjustedFlowError[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // Read the numbers from same event long range correlation
  double sameEventFlowTable[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double sameEventFlowError[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  double positiveEtaTable[2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double negativeEtaTable[2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double closeEtaTable[2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double farEtaTable[2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // Make a big table with all relative uncertainties
  int nUncertainties = 4;
  if(nUncertainties != 4){
    cout << "Check the number of uncertainty sources!. Expected = " << nUncertainties << endl;
    return;
  }
  const char *uncertaintyNames[4] = {"Same event","Relative shift","$\\eta$ side","$\\eta$ region"};
  double relativeUncertaintyTable[nUncertainties][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double absoluteUncertaintyTable[nUncertainties][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // To pass the array as an argument, the numbers will need to be given by hand.
  if(nAsymmetryBins != 3 || nCentralityBins != 4 || nTrackPtBins != 7){
    cout << "Wrong binning for chi2/ndf and background uncertainty histograms!" << endl;
    cout << "nAsymmetry = " << nAsymmetryBins << " expected = 3" << endl;
    cout << "nCentralityBins = " << nCentralityBins << " expected = 4" << endl;
    cout << "nTrackPtBins = " << nTrackPtBins << " expected = 7" << endl;
    cout << "Fix the binning and the code will run again!" << endl;
    return;
  }
  double chi2PerNdf[4][5][7];
  double backgroundAdjustmentUncertainty[4][5][7];
  
  // Read the histograms from the input files
  pbpbReader->SetLoadTrackLeadingJetCorrelations(true);
  pbpbReader->SetAsymmetryBinRange(0,nAsymmetryBins);
  pbpbReader->SetLoad2DHistograms(true);
  pbpbReader->LoadProcessedHistograms();
  
  ppReader->SetLoadTrackLeadingJetCorrelations(true);
  ppReader->SetAsymmetryBinRange(0,nAsymmetryBins);
  ppReader->SetLoad2DHistograms(true);
  ppReader->LoadProcessedHistograms();
  
  pbpbUnadjustedReader->SetLoadTrackLeadingJetCorrelations(true);
  pbpbUnadjustedReader->SetAsymmetryBinRange(0,nAsymmetryBins);
  pbpbUnadjustedReader->LoadProcessedHistograms();
  
  ppUnadjustedReader->SetLoadTrackLeadingJetCorrelations(true);
  ppUnadjustedReader->SetAsymmetryBinRange(0,nAsymmetryBins);
  ppUnadjustedReader->LoadProcessedHistograms();
  
  double backgroundLevel;
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        
        // Read the histograms for pp
        if(iCentrality == nCentralityBins){
          
          // Regular background histogram
          background[iAsymmetry][nCentralityBins][iTrackPt] = ppReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, iAsymmetry, 0, iTrackPt, DijetHistogramManager::kWholeEta);
          
          // Track-leading jet unadjusted background histogram for systematic uncertainty estimation
          leadingUnadjustedBackground[iAsymmetry][nCentralityBins][iTrackPt] = ppUnadjustedReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, iAsymmetry, 0, iTrackPt, DijetHistogramManager::kWholeEta);
          
          // Track-subleading jet unadjusted background histogram for systematic uncertainty estimation
          leadingUnadjustedBackgroundOverlap[iAsymmetry][nCentralityBins][iTrackPt] = ppUnadjustedReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackgroundOverlap, iAsymmetry, 0, iTrackPt, DijetHistogramManager::kWholeEta);
          
          // Same event histograms
          deltaPhiDeltaEtaDistribution[0][0][iAsymmetry][nCentralityBins][iTrackPt] = ppReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kSameEvent, iAsymmetry, 0, iTrackPt);
          
          deltaPhiDeltaEtaDistribution[1][0][iAsymmetry][nCentralityBins][iTrackPt] = ppReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackSubleadingJet, DijetHistogramManager::kSameEvent, iAsymmetry, 0, iTrackPt);
          
          // Mixed event corrected histograms
          deltaPhiDeltaEtaDistribution[0][1][iAsymmetry][nCentralityBins][iTrackPt] = ppReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kCorrected, iAsymmetry, 0, iTrackPt);
          
          deltaPhiDeltaEtaDistribution[1][1][iAsymmetry][nCentralityBins][iTrackPt] = ppReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackSubleadingJet, DijetHistogramManager::kCorrected, iAsymmetry, 0, iTrackPt);
          
        // Read the histograms for PbPb
        } else {
          
          // Regular background histogram
          background[iAsymmetry][iCentrality][iTrackPt] = pbpbReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta);
          
          // Track-leading jet unadjusted background histogram for systematic uncertainty estimation
          leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt] = pbpbUnadjustedReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta);
          
          // Track-subleading jet unadjusted background histogram for systematic uncertainty estimation
          leadingUnadjustedBackgroundOverlap[iAsymmetry][iCentrality][iTrackPt] = pbpbUnadjustedReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackgroundOverlap, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta);
          
          cout << "iAsymmetry: " << iAsymmetry << " iCerntrality: " << iCentrality << " iTrackPt: " << iTrackPt << endl;
          
          // Same event histograms
          deltaPhiDeltaEtaDistribution[0][0][iAsymmetry][iCentrality][iTrackPt] = pbpbReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt);
          
          deltaPhiDeltaEtaDistribution[1][0][iAsymmetry][iCentrality][iTrackPt] = pbpbReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackSubleadingJet, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt);
          
          // Mixed event corrected histograms
          deltaPhiDeltaEtaDistribution[0][1][iAsymmetry][iCentrality][iTrackPt] = pbpbReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kCorrected, iAsymmetry, iCentrality, iTrackPt);
          
          deltaPhiDeltaEtaDistribution[1][1][iAsymmetry][iCentrality][iTrackPt] = pbpbReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackSubleadingJet, DijetHistogramManager::kCorrected, iAsymmetry, iCentrality, iTrackPt);
          
        }
        
        // Do the projection for different eta regions from same event and corrected distributions
        for(int iType = 0; iType < 2; iType++){
        
          positiveEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt] = refitter->CombineDeltaPhi(deltaPhiDeltaEtaDistribution[0][iType][iAsymmetry][iCentrality][iTrackPt], deltaPhiDeltaEtaDistribution[1][iType][iAsymmetry][iCentrality][iTrackPt], 1.5, 2.5, Form("positiveProjection%s%d%d%d",typeString[iType].Data(),iAsymmetry,iCentrality,iTrackPt), true);
          
          negativeEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt] = refitter->CombineDeltaPhi(deltaPhiDeltaEtaDistribution[0][iType][iAsymmetry][iCentrality][iTrackPt], deltaPhiDeltaEtaDistribution[1][iType][iAsymmetry][iCentrality][iTrackPt], -2.5, -1.5, Form("negativeProjection%s%d%d%d",typeString[iType].Data(),iAsymmetry,iCentrality,iTrackPt), true);
          
          closeEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt] = refitter->CombineDeltaPhi(deltaPhiDeltaEtaDistribution[0][iType][iAsymmetry][iCentrality][iTrackPt], deltaPhiDeltaEtaDistribution[1][iType][iAsymmetry][iCentrality][iTrackPt], 1.5, 2, Form("closeProjection%s%d%d%d",typeString[iType].Data(),iAsymmetry,iCentrality,iTrackPt));
          
          farEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt] = refitter->CombineDeltaPhi(deltaPhiDeltaEtaDistribution[0][iType][iAsymmetry][iCentrality][iTrackPt], deltaPhiDeltaEtaDistribution[1][iType][iAsymmetry][iCentrality][iTrackPt], 2, 2.5, Form("farProjection%s%d%d%d",typeString[iType].Data(),iAsymmetry,iCentrality,iTrackPt));
          
          // Fit the projections with a Fourier fit
          positiveEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
          positiveEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
          refitter->FourierFit(positiveEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt], nRefit);
          positiveEtaFit[iType][iAsymmetry][iCentrality][iTrackPt] = positiveEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
          
          negativeEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
          negativeEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
          refitter->FourierFit(negativeEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt], nRefit);
          negativeEtaFit[iType][iAsymmetry][iCentrality][iTrackPt] = negativeEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
          
          closeEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
          closeEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
          refitter->FourierFit(closeEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt], nRefit);
          closeEtaFit[iType][iAsymmetry][iCentrality][iTrackPt] = closeEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
          
          farEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
          farEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
          refitter->FourierFit(farEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt], nRefit);
          farEtaFit[iType][iAsymmetry][iCentrality][iTrackPt] = farEtaProjection[iType][iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        }
        
        // Get the long range correlations directly from the raw distribution
        sameEventLongRange[iAsymmetry][iCentrality][iTrackPt] = refitter->CombineDeltaPhi(deltaPhiDeltaEtaDistribution[0][0][iAsymmetry][iCentrality][iTrackPt], deltaPhiDeltaEtaDistribution[1][0][iAsymmetry][iCentrality][iTrackPt], 1.5, 2.5, Form("sameProjection%d%d%d",iAsymmetry,iCentrality,iTrackPt));
        
        sameEventLongRange[iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
        sameEventLongRange[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
        
        refitter->FourierFit(sameEventLongRange[iAsymmetry][iCentrality][iTrackPt], nRefit);
        sameEventFit[iAsymmetry][iCentrality][iTrackPt] = sameEventLongRange[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        
        // Find the fourier fit function from the histogram
        backgroundFit[iAsymmetry][iCentrality][iTrackPt] = background[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        unadjustedFit[iAsymmetry][iCentrality][iTrackPt] = leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        
        // Fit the standard background distribution
        background[iAsymmetry][iCentrality][iTrackPt]->RecursiveRemove(backgroundFit[iAsymmetry][iCentrality][iTrackPt]);
        background[iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
        background[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
        
        refitter->FourierFit(background[iAsymmetry][iCentrality][iTrackPt], nRefit);
        backgroundFit[iAsymmetry][iCentrality][iTrackPt] = background[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        
        // Fit also the unadjusted background to have same exactly the settings as for adjusted background
        leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt]->RecursiveRemove(unadjustedFit[iAsymmetry][iCentrality][iTrackPt]);
        leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
        leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
        refitter->FourierFit(leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt], nRefit);
        unadjustedFit[iAsymmetry][iCentrality][iTrackPt] = leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        
        
        
      } // asymmetry loop
    } // track pT loop
  } // centrality loop
  
  // Extract the vn parameters from the fits
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = 0; iAsymmetry <=nAsymmetryBins; iAsymmetry++){
        for(int iFlow = 0; iFlow < nRefit; iFlow++){
          masterFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = backgroundFit[iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1);
          masterFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = backgroundFit[iAsymmetry][iCentrality][iTrackPt]->GetParError(iFlow+1);
          unadjustedFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = unadjustedFit[iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1);
          unadjustedFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = unadjustedFit[iAsymmetry][iCentrality][iTrackPt]->GetParError(iFlow+1);
          
          sameEventFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = sameEventFit[iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1);
          sameEventFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = sameEventFit[iAsymmetry][iCentrality][iTrackPt]->GetParError(iFlow+1);
          
          for(int iType = 0; iType < 2; iType++){
            
            if(positiveEtaFit[iType][iAsymmetry][iCentrality][iTrackPt]){
              positiveEtaTable[iType][iAsymmetry][iCentrality][iTrackPt][iFlow] = positiveEtaFit[iType][iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1);
            } else {
              positiveEtaTable[iType][iAsymmetry][iCentrality][iTrackPt][iFlow] = -999;
            }
            
            if(negativeEtaFit[iType][iAsymmetry][iCentrality][iTrackPt]){
              negativeEtaTable[iType][iAsymmetry][iCentrality][iTrackPt][iFlow] = negativeEtaFit[iType][iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1);
            } else {
              negativeEtaTable[iType][iAsymmetry][iCentrality][iTrackPt][iFlow] = -999;
            }
            
            if(closeEtaFit[iType][iAsymmetry][iCentrality][iTrackPt]){
              closeEtaTable[iType][iAsymmetry][iCentrality][iTrackPt][iFlow] = closeEtaFit[iType][iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1);
            } else {
              closeEtaTable[iType][iAsymmetry][iCentrality][iTrackPt][iFlow] = -999;
            }
            
            if(farEtaFit[iType][iAsymmetry][iCentrality][iTrackPt]){
              farEtaTable[iType][iAsymmetry][iCentrality][iTrackPt][iFlow] = farEtaFit[iType][iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1);
            } else {
              farEtaTable[iType][iAsymmetry][iCentrality][iTrackPt][iFlow] = -999;
            }
          } // Loop over same event and corrected
          
        } // flow components
      } // Asymmetry loop
    } // track pT
  } // centrality
  
  double firstValue, secondValue;
  
  // Collect all the numbers in one big table
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        for(int iFlow = 0; iFlow < nRefit; iFlow++){
          relativeUncertaintyTable[0][iAsymmetry][iCentrality][iTrackPt][iFlow] = TMath::Abs(1-sameEventFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow]/unadjustedFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow]);
          absoluteUncertaintyTable[0][iAsymmetry][iCentrality][iTrackPt][iFlow] = TMath::Abs(sameEventFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] - unadjustedFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow]);
          
          relativeUncertaintyTable[1][iAsymmetry][iCentrality][iTrackPt][iFlow] = TMath::Abs(1-masterFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow]/unadjustedFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow]);
          absoluteUncertaintyTable[1][iAsymmetry][iCentrality][iTrackPt][iFlow] = TMath::Abs(masterFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] - unadjustedFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow]);
          
          firstValue = TMath::Abs(1-positiveEtaTable[1][iAsymmetry][iCentrality][iTrackPt][iFlow]/unadjustedFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow]);
          secondValue = TMath::Abs(1-negativeEtaTable[1][iAsymmetry][iCentrality][iTrackPt][iFlow]/unadjustedFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow]);
          relativeUncertaintyTable[2][iAsymmetry][iCentrality][iTrackPt][iFlow] = TMath::Max(firstValue,secondValue);
          
          firstValue = TMath::Abs(positiveEtaTable[1][iAsymmetry][iCentrality][iTrackPt][iFlow] - unadjustedFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow]);
          secondValue = TMath::Abs(negativeEtaTable[1][iAsymmetry][iCentrality][iTrackPt][iFlow] - unadjustedFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow]);
          absoluteUncertaintyTable[2][iAsymmetry][iCentrality][iTrackPt][iFlow] = TMath::Max(firstValue,secondValue);
          
          firstValue = TMath::Abs(1-closeEtaTable[1][iAsymmetry][iCentrality][iTrackPt][iFlow]/unadjustedFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow]);
          secondValue = TMath::Abs(1-farEtaTable[1][iAsymmetry][iCentrality][iTrackPt][iFlow]/unadjustedFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow]);
          relativeUncertaintyTable[3][iAsymmetry][iCentrality][iTrackPt][iFlow] = TMath::Max(firstValue,secondValue);
          
          firstValue = TMath::Abs(closeEtaTable[1][iAsymmetry][iCentrality][iTrackPt][iFlow] - unadjustedFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow]);
          secondValue = TMath::Abs(farEtaTable[1][iAsymmetry][iCentrality][iTrackPt][iFlow] - unadjustedFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow]);
          absoluteUncertaintyTable[3][iAsymmetry][iCentrality][iTrackPt][iFlow] = TMath::Max(firstValue,secondValue);
        }
      }
    }
  } // Asymmetry loop
  
  if(printUncertainties){
    
    // Print a slide with uncertainties for each source and each centrality for each V
    char namer[100];
    for(int iFlow = 0; iFlow < 3; iFlow++){
      for(int iTrackPt = 0; iTrackPt < 5; iTrackPt++){
        
        sprintf(namer,"\\frametitle{Uncertainty for $V_{%d}$, $%.1f < p_{\\mathrm{T}} < %.1f$}", iFlow+1, trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
        
        cout << endl;
        cout << "\\begin{frame}" << endl;
        cout << namer << endl;
        cout << "\\begin{center}" << endl;
        cout << "  \\begin{tabular}{cccccc}" << endl;
        cout << "    \\toprule" << endl;
        cout << "    Source & C: 0-10 & C: 10-30 & C: 30-50 & pp \\\\" << endl;
        cout << "    \\midrule" << endl;
        
        // Set the correct precision for printing floating point numbers
        cout << fixed << setprecision(3);
        
        // First, print the actual values of vn:s
        cout << "$V_{" << iFlow+1 << "}$ ";
        
        for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
          if(iCentrality == nCentralityBins - 1) continue;
          cout << " & $" << unadjustedFlowTable[nAsymmetryBins][iCentrality][iTrackPt][iFlow] << "$";
        }
        cout << " \\\\" << endl;
        cout << "    \\midrule" << endl;
        
        // Print the relative uncertainties
        
        for(int iUncertainty = 0; iUncertainty < nUncertainties; iUncertainty++){
          cout << uncertaintyNames[iUncertainty];
          for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
            if(iCentrality == nCentralityBins - 1) continue;
            cout << " & $" << relativeUncertaintyTable[iUncertainty][nAsymmetryBins][iCentrality][iTrackPt][iFlow] << "$";
          }
          cout << " \\\\" << endl;
        }
        
        cout << "    \\midrule" << endl;
        
        // Print the absolute uncertainties
        
        for(int iUncertainty = 0; iUncertainty < nUncertainties; iUncertainty++){
          cout << uncertaintyNames[iUncertainty];
          for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
            if(iCentrality == nCentralityBins - 1) continue;
            cout << " & $" << absoluteUncertaintyTable[iUncertainty][nAsymmetryBins][iCentrality][iTrackPt][iFlow] << "$";
          }
          cout << " \\\\" << endl;
        }
        
        
        cout << "    \\bottomrule" << endl;
        cout << "  \\end{tabular}" << endl;
        cout << "\\end{center}" << endl;
        cout << "\\begin{itemize}" << endl;
        cout << "  \\item Top: value. Middle: relative uncertainty. Bottom: absolute uncertainty." << endl;
        cout << "\\end{itemize}" << endl;
        cout << "\\end{frame}" << endl;
        cout << endl;
        
      } // Track pT loop
    } // vn loop*/
    
  } // Printing uncertainties
}

