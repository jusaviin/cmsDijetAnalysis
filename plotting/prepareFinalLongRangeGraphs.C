#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "DijetMethods.h"
#include "JDrawer.h"
#include "JffCorrector.h"

/*
 * Macro for preparing long range correlation graphs.
 * The results are fully corrected down to jet vn level.
 * Systematic uncertainties can be saved together with data points.
 */
void prepareFinalLongRangeGraphs(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // File for Vn from jet-hadron correlations
  TString dijetFileName = "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_averagePeakMixing_allCorrections_processed_2020-03-13.root";
  
  // File for Vn from dihadron correlations
  TString dihadronFileName = "data/dihadronPbPb2018_sameTriggerAssoc_5eventMixed_onlySeagull_deltaEta2-3v5_processed_2020-06-18_smallStats.root";
  // data/dihadronPbPb2018_sameTriggerAssoc_5eventMixed_onlySeagull_deltaEta2-3v5_processed_2020-06-18_smallStats.root
  // data/dihadronPbPb2018_sameTriggerAssoc_5eventMixed_noCorrections_processed_2020-06-18_smallStats.root

  // Reconstruction bias correction
  const char *jetReconstructionBiasFile = "corrections/jetReconstructionBiasCorrection_noShiftFitUpToV4_forTestingPurposes.txt";
  
  // Systematic uncertainty configuration
  const char *uncertaintyFile = "data/vnUncertaintyPreliminary2018.txt";
  const bool disableSystematicUncertainty = true;
  
  // Open the input file and read bin numbers from it
  TFile *dijetFile = TFile::Open(dijetFileName);
  DijetHistogramManager *dijetReader = new DijetHistogramManager(dijetFile);
  const int nCentralityBins = dijetReader->GetNCentralityBins();
  const int nTrackPtBins = dijetReader->GetNTrackPtBins();
  const int nAsymmetryBins = dijetReader->GetNAsymmetryBins();
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  const bool useSameEventDijet = false;     // True: Prepare final graphs from raw dijet distribution. False: Use mixed event corrected distributions
  const bool useSameEventDihadron = false;  // True: Prepare final graphs from raw dihadron distribution. False: Use mixed event corrected distributions
  
  const bool drawFourierFitDijet = false;   // Draw the fits done to the dijet distributions
  const bool drawFourierFitDihadron = false;   // Draw the fits done to the dihadron distributions
    
  const bool saveFigures = false;
  TString saveComment = "_jetV2";
  
  const int nRefit = 4; // Number of vn:s included in the refit
  bool refitBackground = true; // Refit the background distribution
  int backgroundRebin = 4; // Rebin applied to the background distributions
  
  // To get the single hadron vn from dihadron vn, we need to divide with the trigger bin vn
  const int dihadronNormalizationBin = -1; // Bin used for normalizing dihadron V2 to hadron v2. For -1, each bin is normalized by the square root of that bin
  
  TString outputFileName = "finalGraphTestNew.root";
  
  // ==================================================================
  // ====================== Configuration done ========================
  // ==================================================================
  
  // Read the dihadron file
  TFile *dihadronFile = TFile::Open(dihadronFileName);
  DijetHistogramManager *dihadronReader = new DijetHistogramManager(dihadronFile);
  
  // Read the uncertainties from the uncertainty file
  JffCorrector *uncertaintyProvider = new JffCorrector();
  uncertaintyProvider->ReadLongRangeSystematicFile(uncertaintyFile);
  uncertaintyProvider->ReadJetReconstructionBiasFile(jetReconstructionBiasFile);
  
  // Define a fitter to test fitting the background distribution with different amount of vn:s
  DijetMethods *refitter = new DijetMethods();
  
  // Define arrays for the histograms
  TH1D *longRangeDijet[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  TF1 *longRangeFitDijet[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  
  TH1D *longRangeDihadron[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  TF1 *longRangeFitDihadron[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  
  // Helper histograms to get the background projection from same event histgorams
  TH2D *helperHistogram;
  TH2D *helperHistogramLeading;
  TH2D *helperHistogramSubleading;
  
  // Define track histograms that are needed to find correct place to put point for the graphs
  TH1D *tracksForGraph[nCentralityBins+1];
  
  // Arrays for extracted vn numbers for dijets
  double dijetFlowTable[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double dijetFlowError[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // Arrays for extracted vn numbers for dijets corrected for jet reconstruction bias
  double dijetFlowTableCorrected[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double dijetFlowErrorCorrected[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // Arrays for extracted vn numbers for dihadrons
  double dihadronFlowTable[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double dihadronFlowError[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // Arrays for extracted vn numbers for single hadrons
  double hadronFlowTable[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double hadronFlowError[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // Arrays for extracted vn numbers for jets
  double jetFlowTable[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double jetFlowError[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // Read the histograms from the input files
  dijetReader->SetLoadTracks(true);
  dijetReader->SetLoadTrackLeadingJetCorrelations(true);
  if(useSameEventDijet) dijetReader->SetLoadTrackSubleadingJetCorrelations(true);
  dijetReader->SetAsymmetryBinRange(0,nAsymmetryBins);
  dijetReader->SetLoad2DHistograms(useSameEventDijet);
  dijetReader->LoadProcessedHistograms();
  
  dihadronReader->SetLoadTracks(true);
  dihadronReader->SetLoadTrackLeadingJetCorrelations(true);
  if(useSameEventDihadron) dihadronReader->SetLoadTrackSubleadingJetCorrelations(true);
  dihadronReader->SetAsymmetryBinRange(0,nAsymmetryBins);
  dihadronReader->SetLoad2DHistograms(useSameEventDihadron);
  dihadronReader->LoadProcessedHistograms();

  
  double binCenter, binError, errorScale;
  int nBins;
  
  // Read the long range dijet histograms and do the fitting
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    // Track histogram for PbPb (needed for graph binning)
    tracksForGraph[iCentrality] = dijetReader->GetHistogramTrackPt(DijetHistogramManager::kTrack, DijetHistogramManager::kSameEvent, iCentrality);
    
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = nAsymmetryBins; iAsymmetry <= nAsymmetryBins; iAsymmetry++){  // TODO: Add asymmetry binning
        
        // Regular long range histogram
        longRangeDijet[iAsymmetry][iCentrality][iTrackPt] = dijetReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta);
        
        // Remove earlier fit from the histogram
        longRangeFitDijet[iAsymmetry][iCentrality][iTrackPt] = longRangeDijet[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        longRangeDijet[iAsymmetry][iCentrality][iTrackPt]->RecursiveRemove(longRangeFitDijet[iAsymmetry][iCentrality][iTrackPt]);
        
        // Read the two dimensional distribution from the same event
        if(useSameEventDijet){
          helperHistogramLeading = dijetReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt);
          helperHistogramLeading->Scale(1.0 / dijetReader->GetPtIntegral(iCentrality, iAsymmetry));
          
          helperHistogramSubleading = dijetReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackSubleadingJet, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt);
          helperHistogramSubleading->Scale(1.0 / dijetReader->GetPtIntegral(iCentrality, iAsymmetry));
          
          refitter->SubtractBackground(helperHistogramLeading, helperHistogramSubleading, 4, false);
                    
          helperHistogram = refitter->GetBackground();
          nBins = helperHistogram->GetNbinsY();
          
          longRangeDijet[iAsymmetry][iCentrality][iTrackPt] = helperHistogram->ProjectionX(Form("sameLongDijet%d%d%d", iAsymmetry, iCentrality, iTrackPt), 1, nBins);  // Exclude underflow and overflow bins by specifying range
          
          longRangeDijet[iAsymmetry][iCentrality][iTrackPt]->Scale(helperHistogram->GetYaxis()->GetBinWidth(1));  // For correct normalization, need to divide out deltaEta bin width of the two-dimensional histogram
          
          /*
           * Do error scaling and for the long range deltaPhi distribution
           *
           * Error scaling is needed because information only from 1.5 < |deltaEta| < 2.5 is used to determine the background.
           * When doing projection, root by default scaled the histogram with square root of bins projected over
           * The whole range of the analysis is |deltaEta| < 4, which means that four times the bins used to determine
           * the background are used in normalixing the error. We can fix this be scaling the error up by sqrt(4) = 2.
           * This number is provided by the DijetMethods, since if we change the limits described above, the number is
           * automatically adjusted for the new limits inside DijetMethods.
           */
          
          for(int iBin = 1; iBin < longRangeDijet[iAsymmetry][iCentrality][iTrackPt]->GetNbinsX(); iBin++){
            binError = longRangeDijet[iAsymmetry][iCentrality][iTrackPt]->GetBinError(iBin);
            errorScale = refitter->GetBackgroundErrorScalingFactor();
            longRangeDijet[iAsymmetry][iCentrality][iTrackPt]->SetBinError(iBin, binError*errorScale);
          }
        }
        
        // Fit the background with Fourier fit
        longRangeDijet[iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
        longRangeDijet[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
        
        refitter->FourierFit(longRangeDijet[iAsymmetry][iCentrality][iTrackPt], nRefit);
        longRangeFitDijet[iAsymmetry][iCentrality][iTrackPt] = longRangeDijet[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        
        
        
      } // asymmetry loop
    } // track pT loop
  } // centrality loop
  
  // Read the long range dihadron histograms and do the fitting
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = nAsymmetryBins; iAsymmetry <= nAsymmetryBins; iAsymmetry++){ // TODO: Add asymmetry binning
        
        // Regular long range histogram
        longRangeDihadron[iAsymmetry][iCentrality][iTrackPt] = dihadronReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta);
        
        // Remove earlier fit from the histogram
        longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt] = longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        if(longRangeDihadron[iAsymmetry][iCentrality][iTrackPt] != NULL) longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->RecursiveRemove(longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt]);
        
        // Read the two dimensional distribution from the same event
        if(useSameEventDihadron){
          helperHistogramLeading = dihadronReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt);
          helperHistogramLeading->Scale(1.0 / dijetReader->GetPtIntegral(iCentrality, iAsymmetry));
          
          helperHistogramSubleading = dihadronReader->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackSubleadingJet, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt);
          helperHistogramSubleading->Scale(1.0 / dijetReader->GetPtIntegral(iCentrality, iAsymmetry));
          
          refitter->SubtractBackground(helperHistogramLeading, helperHistogramSubleading, 4, false);
                    
          helperHistogram = refitter->GetBackground();
          nBins = helperHistogram->GetNbinsY();
          
          longRangeDihadron[iAsymmetry][iCentrality][iTrackPt] = helperHistogram->ProjectionX(Form("sameLongDihadron%d%d%d", iAsymmetry, iCentrality, iTrackPt), 1, nBins);  // Exclude underflow and overflow bins by specifying range
          
          longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->Scale(helperHistogram->GetYaxis()->GetBinWidth(1));  // For correct normalization, need to divide out deltaEta bin width of the two-dimensional histogram
          
          /*
           * Do error scaling and for the long range deltaPhi distribution
           *
           * Error scaling is needed because information only from 1.5 < |deltaEta| < 2.5 is used to determine the background.
           * When doing projection, root by default scaled the histogram with square root of bins projected over
           * The whole range of the analysis is |deltaEta| < 4, which means that four times the bins used to determine
           * the background are used in normalixing the error. We can fix this be scaling the error up by sqrt(4) = 2.
           * This number is provided by the DijetMethods, since if we change the limits described above, the number is
           * automatically adjusted for the new limits inside DijetMethods.
           */
          
          for(int iBin = 1; iBin < longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->GetNbinsX(); iBin++){
            binError = longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->GetBinError(iBin);
            errorScale = refitter->GetBackgroundErrorScalingFactor();
            longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->SetBinError(iBin, binError*errorScale);
          }
        }
        
        // Fit the background with Fourier fit
        longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
        longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
        
        refitter->FourierFit(longRangeDihadron[iAsymmetry][iCentrality][iTrackPt], nRefit);
        longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt] = longRangeDihadron[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        
        
        
      } // asymmetry loop
    } // track pT loop
  } // centrality loop
  
  double hadronFlowNormalizer = 1;
  
  // Extract the vn parameters from the fits
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iAsymmetry = nAsymmetryBins; iAsymmetry <= nAsymmetryBins; iAsymmetry++){  // TODO: Add asymmetry binning
      for(int iFlow = 0; iFlow < nRefit; iFlow++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          // For dijet vn, these numbers can be directly read from the Fourier fits
          dijetFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = longRangeFitDijet[iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1);
          dijetFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = longRangeFitDijet[iAsymmetry][iCentrality][iTrackPt]->GetParError(iFlow+1);
          
          // For dihadron vn, these numbers can be directly read from the Fourier fits
          dihadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1);
          dihadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt]->GetParError(iFlow+1);
          
          if(dihadronNormalizationBin < 0) {
            hadronFlowNormalizer = TMath::Sqrt(dihadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow]);
          } else {
            hadronFlowNormalizer = TMath::Sqrt(dihadronFlowTable[iAsymmetry][iCentrality][dihadronNormalizationBin][iFlow]);
          }
          
          // Single hadron vn are obtained from dihadron vn by dividing out the trigger particle component
          hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1) / hadronFlowNormalizer;
          hadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = longRangeFitDihadron[iAsymmetry][iCentrality][iTrackPt]->GetParError(iFlow+1) / hadronFlowNormalizer;
          
          // The dijet vn needs first be corrected for jet reconstruction bias effects
          dijetFlowTableCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow] = longRangeFitDijet[iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1) - uncertaintyProvider->GetJetReconstructionBiasCorrection(iFlow, iCentrality, iTrackPt, iAsymmetry);;
          dijetFlowErrorCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow] = longRangeFitDijet[iAsymmetry][iCentrality][iTrackPt]->GetParError(iFlow+1);
          
          // To get the final jet vn, we need to divide out the hadron vn from corrected dijet vn
          jetFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = dijetFlowTableCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow] / hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
          jetFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = TMath::Sqrt(TMath::Power(dijetFlowErrorCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow] / hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow],2) + TMath::Power(dijetFlowTableCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow] *  hadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] / TMath::Power(hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow], 2),2));
          
        } // track pT
      } // flow components
    } // Asymmetry loop
    
  } // centrality
  
  // Draw all asymmery bins deltaR curves to a single plot to see if the spillover is similar regardless of asymmetry bin
  JDrawer *drawer = new JDrawer();

  // Helper variables for drawing figures
  TLegend *legend;

  char namerY[100];
  TString asymmetryString[] = {"_A=0v0-0v6", "_A=0v6-0v8", "_A=0v8-1v0", ""};
  TH1D *drawnHistogram;
  TF1 *drawnFit;
  
  if(drawFourierFitDijet){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        for(int iAsymmetry = nAsymmetryBins; iAsymmetry <= nAsymmetryBins; iAsymmetry++){   // TODO: Add astmmetry binning
          
          drawer->DrawHistogram(longRangeDijet[iAsymmetry][iCentrality][iTrackPt], "#Delta#phi", "#frac{dN}{D#Delta#phi}", " ");
          
          // Ass a legend to the figure
          legend = new TLegend(0.3,0.6,0.7,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          
          legend->AddEntry((TObject*) 0,"Dijet Fourier fit","");
          legend->AddEntry((TObject*) 0, Form("C: %.0f-%.0f %%",centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "");
          legend->AddEntry((TObject*) 0, Form("%.1f < p_{T} < %.1f GeV",trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]), "");
          if(iAsymmetry < nAsymmetryBins){
            legend->AddEntry((TObject*) 0, Form("%.1f < x_{j} < %.1f",xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]), "");
          }
          
          
          // Draw the legend
          legend->Draw();
          
          // Save the figures to file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/longRangeFitCheckDijet%s%s_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", saveComment.Data(), asymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
            
          }
        } // asymmetry loop
      } // Track pt Loop
    } // Centrality loop
  } // Drawing Fourier fits

  if(drawFourierFitDihadron){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        for(int iAsymmetry = nAsymmetryBins; iAsymmetry <= nAsymmetryBins; iAsymmetry++){  // TODO: Add asymmetry binning
          
          drawer->DrawHistogram(longRangeDihadron[iAsymmetry][iCentrality][iTrackPt], "#Delta#phi", "#frac{dN}{D#Delta#phi}", " ");
          
          // Ass a legend to the figure
          legend = new TLegend(0.3,0.6,0.7,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          
          legend->AddEntry((TObject*) 0,"Dihadron Fourier fit","");
          legend->AddEntry((TObject*) 0, Form("C: %.0f-%.0f %%",centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "");
          legend->AddEntry((TObject*) 0, Form("%.1f < p_{T} < %.1f GeV",trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]), "");
          if(iAsymmetry < nAsymmetryBins){
            legend->AddEntry((TObject*) 0, Form("%.1f < x_{j} < %.1f",xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]), "");
          }
          
          
          // Draw the legend
          legend->Draw();
          
          // Save the figures to file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/longRangeFitCheckDihadron%s%s_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", saveComment.Data(), asymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
            
          }
        } // asymmetry loop
      } // Track pt Loop
    } // Centrality loop
  } // Drawing Fourier fits
  
  // Contruct graphs from the arrays and save them to a file
  double graphPointsX[nTrackPtBins-2];               // x-axis points in flow graphs
  double graphErrorsX[nTrackPtBins-2];               // No errors for x-axis
  double graphPointsYDijet[nTrackPtBins-2];          // Vn values for dijets
  double graphErrorsYDijet[nTrackPtBins-2];          // Statistical errors for dijet Vn
  double graphPointsYDihadron[nTrackPtBins-2];       // Vn values for dihadrons
  double graphErrorsYDihadron[nTrackPtBins-2];       // Statistical errors for dihadron Vn
  double graphPointsYDijetCorrected[nTrackPtBins-2]; // Vn values for corrected dijets
  double graphErrorsYDijetCorrected[nTrackPtBins-2]; // Statistical errors for corrected dijet Vn
  double graphPointsYHadron[nTrackPtBins-2];         // vn values for hadrons
  double graphErrorsYHadron[nTrackPtBins-2];         // Statistical errors for hadron vn
  double graphPointsYJet[nTrackPtBins-2];            // vn values for jets
  double graphErrorsYJet[nTrackPtBins-2];            // Statistical errors for jet vn
  
  TGraphErrors *flowGraphDijet[nAsymmetryBins+1][nCentralityBins+1][nRefit];
  TGraphErrors *flowGraphDijetCorrected[nAsymmetryBins+1][nCentralityBins+1][nRefit];
  TGraphErrors *flowGraphDihadron[nAsymmetryBins+1][nCentralityBins+1][nRefit];
  TGraphErrors *flowGraphHadron[nAsymmetryBins+1][nCentralityBins+1][nRefit];
  TGraphErrors *flowGraphJet[nAsymmetryBins+1][nCentralityBins+1][nRefit];
  
  int lowPtBin, highPtBin;
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){

    for(int iTrackPt = 0; iTrackPt < nTrackPtBins - 2; iTrackPt++){

      // Find a good place to put the track pT points for the graphs
      lowPtBin = tracksForGraph[iCentrality]->FindBin(trackPtBinBorders[iTrackPt]);
      highPtBin = tracksForGraph[iCentrality]->FindBin(trackPtBinBorders[iTrackPt+1]);
      tracksForGraph[iCentrality]->GetXaxis()->SetRange(lowPtBin,highPtBin);
      graphPointsX[iTrackPt] = tracksForGraph[iCentrality]->GetMean();

      // Initialize other arrays to zero
      graphErrorsX[iTrackPt] = 0;
      graphPointsYDijet[iTrackPt] = 0;
      graphErrorsYDijet[iTrackPt] = 0;
      graphPointsYDijetCorrected[iTrackPt] = 0;
      graphErrorsYDijetCorrected[iTrackPt] = 0;
      graphPointsYDihadron[iTrackPt] = 0;
      graphErrorsYDihadron[iTrackPt] = 0;
      graphPointsYHadron[iTrackPt] = 0;
      graphErrorsYHadron[iTrackPt] = 0;
      graphPointsYJet[iTrackPt] = 0;
      graphErrorsYJet[iTrackPt] = 0;

    } // Track pT loop for x-axis array

    // Create an array for the y-axis and make a graph out of vn values
    for(int iAsymmetry = nAsymmetryBins; iAsymmetry <= nAsymmetryBins; iAsymmetry++){  // TODO: Add asymmetry binning
      for(int iFlow = 0; iFlow < nRefit; iFlow++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins-2; iTrackPt++){
          
          graphPointsYDijet[iTrackPt] = dijetFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
          graphErrorsYDijet[iTrackPt] = dijetFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow];
          flowGraphDijet[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsYDijet, graphErrorsX, graphErrorsYDijet);

          graphPointsYDijetCorrected[iTrackPt] = dijetFlowTableCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow];
          graphErrorsYDijetCorrected[iTrackPt] = dijetFlowErrorCorrected[iAsymmetry][iCentrality][iTrackPt][iFlow];
          flowGraphDijetCorrected[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsYDijetCorrected, graphErrorsX, graphErrorsYDijetCorrected);
          
          graphPointsYDihadron[iTrackPt] = dihadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
          graphErrorsYDihadron[iTrackPt] = dihadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow];
          flowGraphDihadron[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsYDihadron, graphErrorsX, graphErrorsYDihadron);
          
          graphPointsYHadron[iTrackPt] = hadronFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
          graphErrorsYHadron[iTrackPt] = hadronFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow];
          flowGraphHadron[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsYHadron, graphErrorsX, graphErrorsYHadron);
          
          graphPointsYJet[iTrackPt] = jetFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
          graphErrorsYJet[iTrackPt] = jetFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow];
          flowGraphJet[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsYJet, graphErrorsX, graphErrorsYJet);

        } // Track pT loop
      } // Flow component loop
    } // Asymmetry loop
  } // Centrality loop
  
  // Save the graphs to a file
  if(outputFileName.EndsWith(".root")){
    
    // Create the output file
    TFile *outputFile = new TFile(outputFileName,"UPDATE");
    char histogramNamer[100];
    
    // Create a directory for dijet Vn graphs
    sprintf(histogramNamer,"dijetVn");
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iFlow = 0; iFlow < nRefit; iFlow++){
      for(int iAsymmetry = nAsymmetryBins; iAsymmetry <= nAsymmetryBins; iAsymmetry++){  // TODO: Add asymmetry binning
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          flowGraphDijet[iAsymmetry][iCentrality][iFlow]->Write(Form("dijetV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality), TObject::kOverwrite);
        } // Centrality loop
      } // Asymmetry loop
    } // Flow component loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for corrected dijet Vn graphs
    sprintf(histogramNamer,"dijetVnCorrected");
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iFlow = 0; iFlow < nRefit; iFlow++){
      for(int iAsymmetry = nAsymmetryBins; iAsymmetry <= nAsymmetryBins; iAsymmetry++){  // TODO: Add asymmetry binning
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          flowGraphDijetCorrected[iAsymmetry][iCentrality][iFlow]->Write(Form("dijetV%dCorrected_A%dC%d", iFlow+1, iAsymmetry, iCentrality), TObject::kOverwrite);
        } // Centrality loop
      } // Asymmetry loop
    } // Flow component loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for dihadron Vn graphs
    sprintf(histogramNamer,"dihadronVn");
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iFlow = 0; iFlow < nRefit; iFlow++){
      for(int iAsymmetry = nAsymmetryBins; iAsymmetry <= nAsymmetryBins; iAsymmetry++){  // TODO: Add asymmetry binning
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          flowGraphDihadron[iAsymmetry][iCentrality][iFlow]->Write(Form("dihadronV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality), TObject::kOverwrite);
        } // Centrality loop
      } // Asymmetry loop
    } // Flow component loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for hadron vn graphs
    sprintf(histogramNamer,"hadronVn");
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iFlow = 0; iFlow < nRefit; iFlow++){
      for(int iAsymmetry = nAsymmetryBins; iAsymmetry <= nAsymmetryBins; iAsymmetry++){  // TODO: Add asymmetry binning
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          flowGraphHadron[iAsymmetry][iCentrality][iFlow]->Write(Form("hadronV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality), TObject::kOverwrite);
        } // Centrality loop
      } // Asymmetry loop
    } // Flow component loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for dijet Vn graphs
    sprintf(histogramNamer,"jetVn");
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iFlow = 0; iFlow < nRefit; iFlow++){
      for(int iAsymmetry = nAsymmetryBins; iAsymmetry <= nAsymmetryBins; iAsymmetry++){  // TODO: Add asymmetry binning
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          flowGraphJet[iAsymmetry][iCentrality][iFlow]->Write(Form("jetV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality), TObject::kOverwrite);
        } // Centrality loop
      } // Asymmetry loop
    } // Flow component loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Close the output file
    outputFile->Close();
    
  } // Save graphs to output file
}

