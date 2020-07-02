#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"
#include "JDrawer.h"

/*
 * Produce a correction for the fact the jet reconstruction biases the long range correlation distribution.
 * We compare the Fourier fits to RecoGen and GenGen to extract the correction.
 */ 
void produceJetReconstructionBiasCorrection(){

  // Name for the output file
  const char *outputFileName = "";
  // corrections/jetReconstructionBiasCorrection_noShiftFitUpToV4_forTestingPurposes.txt
  // corrections/jetReconstructionBiasCorrectionWithoutShift_forTestingPurposes.txt
  // corrections/jetReconstructionBiasCorrection_forTestingPurposes.txt
  
  // Draw the graphs showing the size of correction
  bool drawCorrection = true;
  bool drawRecoGenFits = false;
  bool drawGenGenFits = true;
  bool drawFits = drawRecoGenFits || drawGenGenFits;
  bool onlyNearSideFit = false;
  
  // Save the illustration plots to file
  bool saveFigures = false;
  TString saveComment = "_recoGenIllustration";
  
  // Define the file names for the MC files used in deriving the correction
  TString recoGenFileName = "data/PbPbMC2018_RecoGen_akFlowJet_noUncorr_noCentShift_improvisedMixing_noCorrections_jet100trigger_processed_2020-06-22.root";
  // data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_xjBins_5pShiftedCent_5eveMix_jet100Trigger_allCorrections_tuning_processed_2019-10-21.root
  // data/PbPbMC2018_RecoGen_akFlowJet_noUncorr_noCentShift_improvisedMixing_noCorrections_jet100trigger_processed_2020-06-22.root
  // data/PbPbMC2018_RecoGen_akPfCsJet_noUncorr_5pCentShift_improvisedMixing_jet100trigger_noCorrections_processed_2020-06-22.root
  TString genGenFileName = "data/PbPbMC2018_GenGen_akFlowJet_noUncorr_noCentShift_xjBins_improvisedMixing_noCorrections_sube0_processed_2020-06-30.root";
  // data/PbPbMC2018_GenGen_akFlowJet_noUncorr_noCentShift_improvisedMixing_noTrigger_noCorrections_processed_2020-06-22.root
  // data/PbPbMC2018_GenGen_akPfCsJet_noUncorr_5pCentShift_improvisedMixing_noTrigger_noCorrections_processed_2020-06-22.root
  // data/PbPbMC2018_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_centShift5_noCorrections_reProcess_processed_2019-10-12.root
  
  // Open the files
  TFile *recoGenFile = TFile::Open(recoGenFileName);
  TFile *genGenFile = TFile::Open(genGenFileName);
  
  // Create readers for the histograms and define binning
  DijetHistogramManager *recoGenReader = new DijetHistogramManager(recoGenFile);
  DijetHistogramManager *genGenReader = new DijetHistogramManager(genGenFile);
  
  const int nCentralityBins = recoGenReader->GetNCentralityBins();
  const int nTrackPtBins = recoGenReader->GetNTrackPtBins();
  const int nAsymmetryBins = recoGenReader->GetNAsymmetryBins();
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  const char *xjString[] = {"0 < x_{j} < 0.6", "0.6 < x_{j} < 0.8", "0.8 < x_{j} < 1.0", "x_{j} inclusive"};
  
  // Define how many v:s is fitted for the correction
  const int nRefit = 4; // Number of vn:s included in the refit
  const int backgroundRebin = 4;  // Rebinning for the background histogram before Fourier fit
  
  // Select if you want to use the distributions directly from the same event before mixing correction
  const bool useSameEvent = false;
  
  // Read the histograms from the input files
  recoGenReader->SetLoadTrackLeadingJetCorrelations(true);
  if(useSameEvent) recoGenReader->SetLoadTrackSubleadingJetCorrelations(true);
  recoGenReader->SetAsymmetryBinRange(0,nAsymmetryBins);
  recoGenReader->SetLoad2DHistograms(useSameEvent);
  recoGenReader->LoadProcessedHistograms();
  
  genGenReader->SetLoadTrackLeadingJetCorrelations(true);
  if(useSameEvent) genGenReader->SetLoadTrackSubleadingJetCorrelations(true);
  genGenReader->SetAsymmetryBinRange(0,nAsymmetryBins);
  genGenReader->SetLoad2DHistograms(useSameEvent);
  genGenReader->LoadProcessedHistograms();
  
  
  // Define arrays for extracted vn numbers
  double recoGenFlowTable[nAsymmetryBins+1][nCentralityBins][nTrackPtBins][nRefit];
  double recoGenFlowError[nAsymmetryBins+1][nCentralityBins][nTrackPtBins][nRefit];
  
  double genGenFlowTable[nAsymmetryBins+1][nCentralityBins][nTrackPtBins][nRefit];
  double genGenFlowError[nAsymmetryBins+1][nCentralityBins][nTrackPtBins][nRefit];
  
  // Extra histogram needed in case we use same events
  TH2D *sameEventDeltaEtaDeltaPhiLeadingRecoGen[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *sameEventDeltaEtaDeltaPhiLeadingGenGen[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
  TH2D *sameEventDeltaEtaDeltaPhiSubleadingRecoGen[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *sameEventDeltaEtaDeltaPhiSubleadingGenGen[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
  TH2D *sameEventDeltaEtaDeltaPhiLongRangeRecoGen[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *sameEventDeltaEtaDeltaPhiLongRangeGenGen[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
  // Variables needed to calculate long range from same event
  int nBins;
  double binError, errorScale;
  
  // Define arrays for the histograms
  TH1D *recoGenLongRange[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TF1 *recoGenLongRangeFit[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
  TH1D *genGenLongRange[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TF1 *genGenLongRangeFit[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
  // Get the background histograms from the files and do the Fourier fit
  
  DijetMethods *refitter = new DijetMethods();
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        
        // If we are doing the same event correction, we need to read the two-dimensional histograms
        if(useSameEvent){
          sameEventDeltaEtaDeltaPhiLeadingRecoGen[iAsymmetry][iCentrality][iTrackPt] = recoGenReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt);
          sameEventDeltaEtaDeltaPhiLeadingRecoGen[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0 / recoGenReader->GetPtIntegral(iCentrality, iAsymmetry));
          
          sameEventDeltaEtaDeltaPhiSubleadingRecoGen[iAsymmetry][iCentrality][iTrackPt] = recoGenReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackSubleadingJet, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt);
          sameEventDeltaEtaDeltaPhiSubleadingRecoGen[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0 / recoGenReader->GetPtIntegral(iCentrality, iAsymmetry));
          
          sameEventDeltaEtaDeltaPhiLeadingGenGen[iAsymmetry][iCentrality][iTrackPt] = genGenReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt);
          sameEventDeltaEtaDeltaPhiLeadingGenGen[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0 / genGenReader->GetPtIntegral(iCentrality, iAsymmetry));
          
          sameEventDeltaEtaDeltaPhiSubleadingGenGen[iAsymmetry][iCentrality][iTrackPt] = genGenReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackSubleadingJet, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt);
          sameEventDeltaEtaDeltaPhiSubleadingGenGen[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0 / genGenReader->GetPtIntegral(iCentrality, iAsymmetry));
          
          
          // Get the background histograms from the same event RecoGen distributions
          refitter->SubtractBackground(sameEventDeltaEtaDeltaPhiLeadingRecoGen[iAsymmetry][iCentrality][iTrackPt], sameEventDeltaEtaDeltaPhiSubleadingRecoGen[iAsymmetry][iCentrality][iTrackPt], 4, false);
          
          sameEventDeltaEtaDeltaPhiLongRangeRecoGen[iAsymmetry][iCentrality][iTrackPt] = refitter->GetBackground();
          
          nBins = sameEventDeltaEtaDeltaPhiLongRangeRecoGen[iAsymmetry][iCentrality][iTrackPt]->GetNbinsY();
          recoGenLongRange[iAsymmetry][iCentrality][iTrackPt] = sameEventDeltaEtaDeltaPhiLongRangeRecoGen[iAsymmetry][iCentrality][iTrackPt]->ProjectionX(Form("sameLongRecoGen%d%d%d", iAsymmetry, iCentrality, iTrackPt), 1, nBins);  // Exclude underflow and overflow bins by specifying range
          
          recoGenLongRange[iAsymmetry][iCentrality][iTrackPt]->Scale( sameEventDeltaEtaDeltaPhiLongRangeRecoGen[iAsymmetry][iCentrality][iTrackPt]->GetYaxis()->GetBinWidth(1) );  // For correct normalization, need to divide out deltaEta bin width of the two-dimensional histogram
          
          /*
           * Do error scaling and for the background deltaPhi distribution
           *
           * Error scaling is needed because information only from 1.5 < |deltaEta| < 2.5 is used to determine the background.
           * When doing projection, root by default scaled the histogram with square root of bins projected over
           * The whole range of the analysis is |deltaEta| < 4, which means that four times the bins used to determine
           * the background are used in normalixing the error. We can fix this be scaling the error up by sqrt(4) = 2.
           * This number is provided by the DijetMethods, since if we change the limits described above, the number is
           * automatically adjusted for the new limits inside DijetMethods.
           */
          
          for(int iBin = 1; iBin < recoGenLongRange[iAsymmetry][iCentrality][iTrackPt]->GetNbinsX(); iBin++){
            binError = recoGenLongRange[iAsymmetry][iCentrality][iTrackPt]->GetBinError(iBin);
            errorScale = refitter->GetBackgroundErrorScalingFactor();
            recoGenLongRange[iAsymmetry][iCentrality][iTrackPt]->SetBinError(iBin, binError*errorScale);
          }
          
          // Get the background histograms from the same event GenGen distributions
          refitter->SubtractBackground(sameEventDeltaEtaDeltaPhiLeadingGenGen[iAsymmetry][iCentrality][iTrackPt], sameEventDeltaEtaDeltaPhiSubleadingGenGen[iAsymmetry][iCentrality][iTrackPt], 4, false);
          
          sameEventDeltaEtaDeltaPhiLongRangeGenGen[iAsymmetry][iCentrality][iTrackPt] = refitter->GetBackground();
          
          nBins = sameEventDeltaEtaDeltaPhiLongRangeGenGen[iAsymmetry][iCentrality][iTrackPt]->GetNbinsY();
          genGenLongRange[iAsymmetry][iCentrality][iTrackPt] = sameEventDeltaEtaDeltaPhiLongRangeGenGen[iAsymmetry][iCentrality][iTrackPt]->ProjectionX(Form("sameLongGenGen%d%d%d", iAsymmetry, iCentrality, iTrackPt), 1, nBins);  // Exclude underflow and overflow bins by specifying range
          
          genGenLongRange[iAsymmetry][iCentrality][iTrackPt]->Scale( sameEventDeltaEtaDeltaPhiLongRangeGenGen[iAsymmetry][iCentrality][iTrackPt]->GetYaxis()->GetBinWidth(1) );  // For correct normalization, need to divide out deltaEta bin width of the two-dimensional histogram
          
          /*
           * Do error scaling and for the background deltaPhi distribution
           *
           * Error scaling is needed because information only from 1.5 < |deltaEta| < 2.5 is used to determine the background.
           * When doing projection, root by default scaled the histogram with square root of bins projected over
           * The whole range of the analysis is |deltaEta| < 4, which means that four times the bins used to determine
           * the background are used in normalixing the error. We can fix this be scaling the error up by sqrt(4) = 2.
           * This number is provided by the DijetMethods, since if we change the limits described above, the number is
           * automatically adjusted for the new limits inside DijetMethods.
           */
          
          for(int iBin = 1; iBin < genGenLongRange[iAsymmetry][iCentrality][iTrackPt]->GetNbinsX(); iBin++){
            binError = genGenLongRange[iAsymmetry][iCentrality][iTrackPt]->GetBinError(iBin);
            errorScale = refitter->GetBackgroundErrorScalingFactor();
            genGenLongRange[iAsymmetry][iCentrality][iTrackPt]->SetBinError(iBin, binError*errorScale);
          }
          
        } else {
          // Read the background distributions directly from the file
          
          recoGenLongRange[iAsymmetry][iCentrality][iTrackPt] = recoGenReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta);
          
          genGenLongRange[iAsymmetry][iCentrality][iTrackPt] = genGenReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta);
          
          // Find the fourier fit function from the histogram
          recoGenLongRangeFit[iAsymmetry][iCentrality][iTrackPt] = recoGenLongRange[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
          genGenLongRangeFit[iAsymmetry][iCentrality][iTrackPt] = genGenLongRange[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
          
          // Remove possible previous fits and optionally do rebinning
          recoGenLongRange[iAsymmetry][iCentrality][iTrackPt]->RecursiveRemove(recoGenLongRangeFit[iAsymmetry][iCentrality][iTrackPt]);
          recoGenLongRange[iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
          recoGenLongRange[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
          
          genGenLongRange[iAsymmetry][iCentrality][iTrackPt]->RecursiveRemove(genGenLongRangeFit[iAsymmetry][iCentrality][iTrackPt]);
          genGenLongRange[iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
          genGenLongRange[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
          
        }
        
        // Fit the distributions
        refitter->FourierFit(recoGenLongRange[iAsymmetry][iCentrality][iTrackPt], nRefit, onlyNearSideFit);
        recoGenLongRangeFit[iAsymmetry][iCentrality][iTrackPt] = recoGenLongRange[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        
        refitter->FourierFit(genGenLongRange[iAsymmetry][iCentrality][iTrackPt], nRefit, onlyNearSideFit);
        genGenLongRangeFit[iAsymmetry][iCentrality][iTrackPt] = genGenLongRange[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        
      } // asymmetry loop
    } // track pT loop
  } // centrality loop
  
  // Extract the vn parameters from the fits
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        for(int iFlow = 0; iFlow < nRefit; iFlow++){
          recoGenFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = recoGenLongRangeFit[iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1);
          recoGenFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = recoGenLongRangeFit[iAsymmetry][iCentrality][iTrackPt]->GetParError(iFlow+1);
          genGenFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = genGenLongRangeFit[iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1);
          genGenFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = genGenLongRangeFit[iAsymmetry][iCentrality][iTrackPt]->GetParError(iFlow+1);
        } // flow components
      } // Asymmetry loop
    } // track pT
  } // centrality
  
  
  // If an output file name is given, put the total uncertainties to an output file
  if(strcmp(outputFileName, "") != 0){
    std::ofstream outputFile;
    outputFile.open(outputFileName);
    
    // Print first the information about number of bins
    outputFile << nCentralityBins << " " << nRefit << " " << nAsymmetryBins << " " << nTrackPtBins << endl;
    
    // Next, print all the uncertainties in a specific order
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iFlow = 0; iFlow < nRefit; iFlow++){
        for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
          for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
            outputFile << recoGenFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] - genGenFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] << " ";
          } // Track pT loop
          outputFile << endl;
        } // Asymmetry loop
      } // Flow component loop
    } // Centrality loop
    
    outputFile.close();
  }
  
  // Draw the correction to check it looks nice
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceGraph();
  
  // Legend
  TLegend *legend;
  
  // Naming helper
  TString compactAsymmetryString;
  
  if(drawCorrection){
    
    // Define the needed graphs
    TGraphErrors *correctionGraph[nAsymmetryBins+1][nCentralityBins];
    TGraphErrors *recoGenGraph[nAsymmetryBins+1][nCentralityBins];
    TGraphErrors *genGenGraph[nAsymmetryBins+1][nCentralityBins];
    
    // Default values for x-axis
    double defaultXpoints[] = {0.85, 1.5, 2.5, 3.5, 6, 10, 14};
    double defaultXerrors[] = {0, 0, 0, 0, 0, 0, 0};
    
    double yPoints[3][nTrackPtBins];  // First bin: correction/RecoGen/GenGen
    double yErrors[3][nTrackPtBins];  // First bin: correction/RecoGen/GenGen
    
    // Construct the graphs and draw them
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      
      // Set the asymmetry string based on the selected asymmetry bin
      if(iAsymmetry < nAsymmetryBins){
        compactAsymmetryString = Form("_A=%.1f-%.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
        compactAsymmetryString.ReplaceAll(".","v");
      } else {
        compactAsymmetryString = "";
      }
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        
        // Collect the y-axis information to arrays
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          yPoints[0][iTrackPt] = recoGenFlowTable[iAsymmetry][iCentrality][iTrackPt][1] - genGenFlowTable[iAsymmetry][iCentrality][iTrackPt][1];
          yPoints[1][iTrackPt] = recoGenFlowTable[iAsymmetry][iCentrality][iTrackPt][1];
          yPoints[2][iTrackPt] = genGenFlowTable[iAsymmetry][iCentrality][iTrackPt][1];
          
          yErrors[0][iTrackPt] = TMath::Sqrt(TMath::Power(recoGenFlowError[iAsymmetry][iCentrality][iTrackPt][1],2) + TMath::Power(genGenFlowError[iAsymmetry][iCentrality][iTrackPt][1],2));
          yErrors[1][iTrackPt] = recoGenFlowError[iAsymmetry][iCentrality][iTrackPt][1];
          yErrors[2][iTrackPt] = genGenFlowError[iAsymmetry][iCentrality][iTrackPt][1];
        
        } // Track pT loop
        
        // Create the graphs
        correctionGraph[iAsymmetry][iCentrality] = new TGraphErrors(nTrackPtBins, defaultXpoints, yPoints[0], defaultXerrors, yErrors[0]);
        recoGenGraph[iAsymmetry][iCentrality] = new TGraphErrors(nTrackPtBins, defaultXpoints, yPoints[1], defaultXerrors, yErrors[1]);
        genGenGraph[iAsymmetry][iCentrality] = new TGraphErrors(nTrackPtBins, defaultXpoints, yPoints[2], defaultXerrors, yErrors[2]);
        
        // Set the style for graphs
        correctionGraph[iAsymmetry][iCentrality]->SetMarkerStyle(kFullStar);
        correctionGraph[iAsymmetry][iCentrality]->SetMarkerColor(kBlack);
        
        recoGenGraph[iAsymmetry][iCentrality]->SetMarkerStyle(kFullCircle);
        recoGenGraph[iAsymmetry][iCentrality]->SetMarkerColor(kBlue);
        
        genGenGraph[iAsymmetry][iCentrality]->SetMarkerStyle(kFullSquare);
        genGenGraph[iAsymmetry][iCentrality]->SetMarkerColor(kRed);
        
        drawer->DrawGraph(correctionGraph[iAsymmetry][iCentrality], 0, 12, -0.02, 0.2, "p_{T} (GeV)", "v_{2}", " ", "p");
        recoGenGraph[iAsymmetry][iCentrality]->Draw("p,same");
        genGenGraph[iAsymmetry][iCentrality]->Draw("p,same");
        
        // Draw a legend to the graph
        legend = new TLegend(0.2,0.7,0.5,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
        legend->SetHeader(Form("C = %.0f-%.0f, %s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], xjString[iAsymmetry]));
        legend->AddEntry(correctionGraph[iAsymmetry][iCentrality], "Correction", "p");
        legend->AddEntry(recoGenGraph[iAsymmetry][iCentrality], "RecoGen", "p");
        legend->AddEntry(genGenGraph[iAsymmetry][iCentrality], "GenGen", "p");
        legend->Draw();
        
        // Save the figures to file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/jetFlowcorrectionIllustration%s_v2_C=%.0f-%.0f%s.pdf", saveComment.Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], compactAsymmetryString.Data()));
          
        }
        
      } // Centrality loop
    } // Asymmetry loop
    
    
  } // Draw corrections
  
  // Draw the fits in each bin to see that nothing crazy is happening
  if(drawFits){
    
    drawer->Reset();
    
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      
      // Set the asymmetry string based on the selected asymmetry bin
      if(iAsymmetry < nAsymmetryBins){
        compactAsymmetryString = Form("_A=%.1f-%.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
        compactAsymmetryString.ReplaceAll(".","v");
      } else {
        compactAsymmetryString = "";
      }
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        
        // Collect the y-axis information to arrays
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          if(drawRecoGenFits){
            drawer->DrawHistogram(recoGenLongRange[iAsymmetry][iCentrality][iTrackPt], "#Delta#varphi", "#frac{dN}{d#Delta#varphi}", " ");
            
            legend = new TLegend(0.2,0.7,0.5,0.9);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
            legend->SetHeader("RecoGen");
            legend->AddEntry((TObject*) 0, Form("C = %.0f-%.0f, %s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], xjString[iAsymmetry]), "");
            legend->AddEntry((TObject*) 0, Form("%.1f < pT < %.1f GeV", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]), "");
            
            legend->Draw();
          }
          
          if(drawGenGenFits){
            drawer->DrawHistogram(genGenLongRange[iAsymmetry][iCentrality][iTrackPt], "#Delta#varphi", "#frac{dN}{d#Delta#varphi}", " ");
            
            legend = new TLegend(0.2,0.7,0.5,0.9);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
            legend->SetHeader("GenGen");
            legend->AddEntry((TObject*) 0, Form("C = %.0f-%.0f, %s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], xjString[iAsymmetry]), "");
            legend->AddEntry((TObject*) 0, Form("%.1f < pT < %.1f GeV", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]), "");
            
            legend->Draw();
          }
          
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
    
  }
  
  
}
