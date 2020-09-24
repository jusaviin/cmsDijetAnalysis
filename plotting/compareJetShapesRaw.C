#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"
#include "DijetMethods.h"
#include "JffCorrector.h"

/*
 * Macro for comparing jets shapes in different pT bins between pp and PbPb
 */
void compareJetShapesRaw(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString fileNames[] = { "data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_xjBins_5eveMix_allCorrectionsForClosure_tuning_processed_2020-09-21.root", "data/PbPbMC2018_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_wtaAxis_sube0_centShift5_noCorrections_onlyFinalResults_processed_2019-10-12.root" };
  // data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_xjBins_5eveMix_allCorrectionsForClosure_tuning_processed_2020-09-15.root
  // data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_xjBins_5eveMix_allCorrectionsForClosure_tuning_processed_2020-09-21.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_allCorrections_processed_2020-02-17.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_xjBins_20eveMix_trigWeight_allCorrectionsNonWeighted_wtaAxis_processed_2020-04-30.root
  // data/dijetPbPb2018_akFlowJet_improvisedMixing_jet80trigger_preprocessed_2020-05-26.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_20eveMix_20pSmear_wtaAxis_preprocessed_2020-05-13.root

  
  TString uncertaintyFileName = "uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_manualFluctuationReduction_noTriggerBias_jffUpdate_2020-06-03.root";
  TString uncertaintyFileNamePlus2p = "uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_manualFluctuationReduction_addTriggerBias_jffUpdate_2020-06-03.root";
  
  const int nFilesToCompare = 2;
  const bool mcLabel = false;
  
  // Open all the files for the comparison
  TFile *files[nFilesToCompare];
  for(int iFile = 0; iFile < nFilesToCompare; iFile++){
    files[iFile] = TFile::Open(fileNames[iFile]);
  }
  
  // Create histogram managers to read the histograms from the files
  DijetHistogramManager *histograms[nFilesToCompare];
  for(int iFile = 0; iFile < nFilesToCompare; iFile++){
    histograms[iFile] = new DijetHistogramManager(files[iFile]);
  }
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = false;
  
  // Get the number of asymmetry bins
  const int nAsymmetryBins = histograms[0]->GetNAsymmetryBins();
  const int nTrackPtBins = histograms[0]->GetNTrackPtBins();
  const int nCentralityBins = histograms[0]->GetNCentralityBins();
  
  double centralityBinBorders[] = {0,10,30,50,90};
  double asymmetryBinBorders[] = {0,0.6,0.8,1};
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};
  
  int firstDrawnCentralityBin = 3;
  int lastDrawnCentralityBin = 3;
  
  int firstDrawnAsymmetryBin = 0;
  int lastDrawnAsymmetryBin = nAsymmetryBins;
  
  int iJetTrack = DijetHistogramManager::kTrackSubleadingJet; // kTrackLeadingJet kTrackSubleadingJet kPtWeightedTrackLeadingJet kPtWeightedTrackSubleadingJet
  bool doSameEvent = false;   // true = Same event. False = Final corrected distribution
  
  bool normalizeJetShape = false;
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Load leading and subleading correlation histograms from the files
  for(int iFile = 0; iFile < nFilesToCompare; iFile++){
    histograms[iFile]->SetLoadTrackLeadingJetCorrelations(iJetTrack == DijetHistogramManager::kTrackLeadingJet);
    histograms[iFile]->SetLoadTrackSubleadingJetCorrelations(iJetTrack == DijetHistogramManager::kTrackSubleadingJet);
    histograms[iFile]->SetLoadTrackLeadingJetCorrelationsPtWeighted(iJetTrack == DijetHistogramManager::kPtWeightedTrackLeadingJet);
    histograms[iFile]->SetLoadTrackSubleadingJetCorrelationsPtWeighted(iJetTrack == DijetHistogramManager::kPtWeightedTrackSubleadingJet);
    histograms[iFile]->SetLoad2DHistograms(doSameEvent);
    histograms[iFile]->SetAsymmetryBinRange(0,nAsymmetryBins);
    if(iFile == 0) histograms[iFile]->SetAvoidMixingPeak(true);
    histograms[iFile]->LoadProcessedHistograms();
  }

  // Load the systematic uncertainties
  TFile *uncertaintyFile = TFile::Open(uncertaintyFileName);
  JffCorrector *uncertaintyProvider = new JffCorrector();
  uncertaintyProvider->ReadSystematicFile(uncertaintyFile);
  
  TFile *uncertaintyFilePlus2p = TFile::Open(uncertaintyFileNamePlus2p);
  JffCorrector *uncertaintyProviderPlus2p = new JffCorrector();
  uncertaintyProviderPlus2p->ReadSystematicFile(uncertaintyFilePlus2p);
  
  // Methods for doing stuff
  DijetMethods *projector = new DijetMethods();
 
  // Histograms needed to calculate the deltaR distributions
  TH1D *jetShapePbPb[2][nAsymmetryBins+1][nCentralityBins][nTrackPtBins+1];
  TH1D *jetShapeRatio[nAsymmetryBins+1][nCentralityBins][nTrackPtBins+1];
  double jetShapeNormalizerPbPb[2][nAsymmetryBins+1][nCentralityBins];
  TH2D *helperHistogram;
  
  // Histograms for systematic uncertainties
  const int nReleventUncertainties = 9;
  TH1D *jetShapeUncertainty[nReleventUncertainties][nAsymmetryBins+1][nCentralityBins][nTrackPtBins+1];
  TH1D *jetShapeUncertaintyPlus2p[nReleventUncertainties][nAsymmetryBins+1][nCentralityBins][nTrackPtBins+1];
  TH1D *jetShapeRatioUncertaintyAll[nAsymmetryBins+1][nCentralityBins][nTrackPtBins+1];
  TH1D *jetShapeRatioUncertaintyRelevant[nAsymmetryBins+1][nCentralityBins][nTrackPtBins+1];
  TH1D *jetShapeRatioUncertaintyRelevantPlus2[nAsymmetryBins+1][nCentralityBins][nTrackPtBins+1];
  
  // For normalization of the jet shapes, read the pT summed histograms
  for(int iFile = 0; iFile < 2; iFile++){
    for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
      
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        if(doSameEvent){
          helperHistogram = histograms[iFile]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, nTrackPtBins);
          helperHistogram->Scale(1.0/histograms[iFile]->GetPtIntegral(iCentrality,iAsymmetry));
          helperHistogram->SetName(Form("lol%d%d%d",iFile,iAsymmetry,iCentrality));
          jetShapePbPb[iFile][iAsymmetry][iCentrality][nTrackPtBins] = (TH1D*) projector->GetJetShape(helperHistogram)->Clone(Form("PbPbNormalizer%d%d",iFile,iCentrality));
        } else {
          jetShapePbPb[iFile][iAsymmetry][iCentrality][nTrackPtBins] = (TH1D*) histograms[iFile]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, nTrackPtBins )->Clone(Form("PbPbNormalizer%d%d",iFile,iCentrality));
        }
        jetShapeNormalizerPbPb[iFile][iAsymmetry][iCentrality] = 1;
        if(normalizeJetShape) jetShapeNormalizerPbPb[iFile][iAsymmetry][iCentrality] = jetShapePbPb[iFile][iAsymmetry][iCentrality][nTrackPtBins]->Integral(1, jetShapePbPb[iFile][iAsymmetry][iCentrality][nTrackPtBins]->FindBin(0.99), "width");
      }
      
      // Read the jet shape histograms from the file and normalize them to one
      for(int iTrackPt = 0; iTrackPt <= nTrackPtBins; iTrackPt++){
        for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
          if(doSameEvent){
            helperHistogram = histograms[iFile]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kSameEvent, iAsymmetry, iCentrality, iTrackPt);
            if(iTrackPt != nTrackPtBins) helperHistogram->Scale(1.0/histograms[iFile]->GetPtIntegral(iCentrality,iAsymmetry));
            helperHistogram->SetName(Form("lul%d%d%d%d",iFile,iAsymmetry,iCentrality,iTrackPt));
            jetShapePbPb[iFile][iAsymmetry][iCentrality][iTrackPt] = projector->GetJetShape(helperHistogram);
            jetShapePbPb[iFile][iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/jetShapeNormalizerPbPb[iFile][iAsymmetry][iCentrality]);
          } else {
            jetShapePbPb[iFile][iAsymmetry][iCentrality][iTrackPt] = histograms[iFile]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iTrackPt);
            if(iTrackPt != nTrackPtBins) jetShapePbPb[iFile][iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/jetShapeNormalizerPbPb[iFile][iAsymmetry][iCentrality]);
          }
          
          
          if(iFile == 0) jetShapeRatio[iAsymmetry][iCentrality][iTrackPt] = (TH1D*) jetShapePbPb[iFile][iAsymmetry][iCentrality][iTrackPt]->Clone(Form("jetShapeRatio%d%d%d", iAsymmetry, iCentrality, iTrackPt));
          if(iFile == 1) jetShapeRatio[iAsymmetry][iCentrality][iTrackPt]->Divide(jetShapePbPb[iFile][iAsymmetry][iCentrality][iTrackPt]);
        } // Centrality loop
      } // Track pT loop
    } // Asymmetry loop
  } // File loop
    
  // Read the systematic uncertainty histograms and calculate uncertainty for the ratio
  double uncertaintyValue;
  double numerator;
  double denominator;
  double ratioUncertaintyValue;
  
  /*for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      for(int iTrackPt = 0; iTrackPt <= nTrackPtBins; iTrackPt++){
        jetShapeUncertainty[0][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kFragmentationBias);
        jetShapeUncertainty[1][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kBackgroundFluctuation);
        jetShapeUncertainty[2][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kPairAcceptance);
        jetShapeUncertainty[3][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kBackgroundSubtraction);
        jetShapeUncertainty[4][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kJetEnergyScale);
        jetShapeUncertainty[5][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kJetResolution);
        jetShapeUncertainty[6][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kTrackingEfficiency);
        jetShapeUncertainty[7][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kResidualTracking);
        jetShapeUncertainty[8][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kTrackingDeltaR);
        
        jetShapeUncertaintyPlus2p[0][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProviderPlus2p->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kFragmentationBias);
        jetShapeUncertaintyPlus2p[1][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProviderPlus2p->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kBackgroundFluctuation);
        jetShapeUncertaintyPlus2p[2][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProviderPlus2p->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kPairAcceptance);
        jetShapeUncertaintyPlus2p[3][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProviderPlus2p->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kBackgroundSubtraction);
        jetShapeUncertaintyPlus2p[4][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProviderPlus2p->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kJetEnergyScale);
        jetShapeUncertaintyPlus2p[5][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProviderPlus2p->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kJetResolution);
        jetShapeUncertaintyPlus2p[6][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProviderPlus2p->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kTrackingEfficiency);
        jetShapeUncertaintyPlus2p[7][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProviderPlus2p->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kResidualTracking);
        jetShapeUncertaintyPlus2p[8][iAsymmetry][iCentrality][iTrackPt] = uncertaintyProviderPlus2p->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry, JffCorrector::kTrackingDeltaR);
        
        jetShapeRatioUncertaintyAll[iAsymmetry][iCentrality][iTrackPt] = (TH1D*) jetShapeRatio[iAsymmetry][iCentrality][iTrackPt]->Clone(Form("uncertaintyForRatioAll%d%d%d", iAsymmetry, iCentrality, iTrackPt));
        jetShapeRatioUncertaintyRelevant[iAsymmetry][iCentrality][iTrackPt] = (TH1D*) jetShapeRatio[iAsymmetry][iCentrality][iTrackPt]->Clone(Form("uncertaintyForRatioRelevant%d%d%d", iAsymmetry, iCentrality, iTrackPt));
        jetShapeRatioUncertaintyRelevantPlus2[iAsymmetry][iCentrality][iTrackPt] = (TH1D*) jetShapeRatio[iAsymmetry][iCentrality][iTrackPt]->Clone(Form("uncertaintyForRatioRelevantPlus2%d%d%d", iAsymmetry, iCentrality, iTrackPt));
        for(int iBin = 1; iBin <= jetShapeRatio[iAsymmetry][iCentrality][iTrackPt]->GetNbinsX(); iBin++){
          
          numerator = jetShapePbPb[0][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin);
          denominator = jetShapePbPb[1][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin);
          
          // Relevant errors + 2 %
          uncertaintyValue = TMath::Sqrt(TMath::Power(jetShapeUncertaintyPlus2p[0][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin),2) + TMath::Power(jetShapeUncertaintyPlus2p[1][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin),2) + TMath::Power(jetShapeUncertaintyPlus2p[2][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin),2) + TMath::Power(jetShapeUncertaintyPlus2p[3][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin),2));
          
          ratioUncertaintyValue = TMath::Sqrt(TMath::Power(uncertaintyValue/denominator,2) + TMath::Power(numerator*uncertaintyValue/TMath::Power(denominator,2),2));
          
          jetShapeRatioUncertaintyRelevantPlus2[iAsymmetry][iCentrality][iTrackPt]->SetBinError(iBin, ratioUncertaintyValue);
          
          // Relevant errors
          uncertaintyValue = TMath::Sqrt(TMath::Power(jetShapeUncertainty[0][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin)/2.5 ,2)+ TMath::Power(jetShapeUncertainty[1][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin),2) + TMath::Power(jetShapeUncertainty[2][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin),2) + TMath::Power(jetShapeUncertainty[3][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin),2));
          
          ratioUncertaintyValue = TMath::Sqrt(TMath::Power(uncertaintyValue/denominator,2) + TMath::Power(numerator*uncertaintyValue/TMath::Power(denominator,2),2));
          
          jetShapeRatioUncertaintyRelevant[iAsymmetry][iCentrality][iTrackPt]->SetBinError(iBin, ratioUncertaintyValue);
          
          // All errors
          uncertaintyValue = TMath::Sqrt(TMath::Power(jetShapeUncertaintyPlus2p[0][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin),2) + TMath::Power(jetShapeUncertaintyPlus2p[1][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin),2) + TMath::Power(jetShapeUncertaintyPlus2p[2][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin),2) + TMath::Power(jetShapeUncertaintyPlus2p[3][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin),2) + TMath::Power(jetShapeUncertaintyPlus2p[4][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin),2) + TMath::Power(jetShapeUncertaintyPlus2p[5][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin),2) + TMath::Power(jetShapeUncertaintyPlus2p[6][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin),2) + TMath::Power(jetShapeUncertaintyPlus2p[7][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin),2) + TMath::Power(jetShapeUncertaintyPlus2p[8][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin),2));
          
          ratioUncertaintyValue = TMath::Sqrt(TMath::Power(uncertaintyValue/denominator,2) + TMath::Power(numerator*uncertaintyValue/TMath::Power(denominator,2),2));
          
          jetShapeRatioUncertaintyAll[iAsymmetry][iCentrality][iTrackPt]->SetBinError(iBin, ratioUncertaintyValue);
        }
      }
    }
  }*/
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  
  TLegend *legend;
  TLegend *ratioLegend;
  TLegend *ratioLegend2;
  TLegend *ratioLegend3;
  
  TString centralityString;
  TString compactCentralityString;
  TString trackString;
  TString compactTrackString;
  TString asymmetryString;
  TString compactAsymmetryString;
  
  TString figureName;
  int veryNiceColors[] = {kBlue,kRed,kRed,kGreen+3};
  const char *labels[] = {"Leading","Leading","Leading","Subleading","Subleading","Subleading","Inclusive","Inclusive"};
  TLine *oneLine = new TLine(0,1,1,1);
  oneLine->SetLineStyle(2);
  
  TString yAxisName = "P(#Deltar)";
  if(normalizeJetShape) yAxisName = "#rho(#Deltar)";
  
  for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
    
    // Setup asymmetry strings
    if(iAsymmetry < nAsymmetryBins){
      asymmetryString = Form("%.1f < x_{j} < %.1f", asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1]);
      compactAsymmetryString = Form("_A=%.1f-%.1f", asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1]);
      compactAsymmetryString.ReplaceAll(".","v");
    } else {
      asymmetryString = "";
      compactAsymmetryString = "";
    }
    
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("%.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
      compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
      
      for(int iTrackPt = 0; iTrackPt <= nTrackPtBins; iTrackPt++){
        
        drawer->CreateSplitCanvas();
        drawer->SetLogY(true);
        
        if(iTrackPt < nTrackPtBins){
          trackString = Form("%.1f < p_{T} < %.1f GeV",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
          compactTrackString = Form("_T=%.0f-%.0f",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
        } else {
          trackString = "0.7 < p_{T} < 300 GeV";
          compactTrackString = "_T=1-300";
        }
        
        jetShapePbPb[0][iAsymmetry][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
        //jetShapePbPb[0][iAsymmetry][iCentrality][iTrackPt]>GetYaxis()->SetRangeUser(1.2e-2,20); // TODO: Remove this line
        jetShapePbPb[0][iAsymmetry][iCentrality][iTrackPt]->SetMarkerStyle(21);
        jetShapePbPb[0][iAsymmetry][iCentrality][iTrackPt]->SetMarkerColor(kBlack);
        if(iTrackPt == nTrackPtBins) cout << "pp integral: " << jetShapePbPb[0][iAsymmetry][iCentrality][iTrackPt]->Integral(1, jetShapePbPb[0][iAsymmetry][iCentrality][iTrackPt]->FindBin(0.99), "width") << endl;
        drawer->DrawHistogramToUpperPad(jetShapePbPb[0][iAsymmetry][iCentrality][iTrackPt],"#Deltar",yAxisName," ");
        
        jetShapePbPb[1][iAsymmetry][iCentrality][iTrackPt]->SetLineColor(veryNiceColors[1]);
        jetShapePbPb[1][iAsymmetry][iCentrality][iTrackPt]->SetMarkerStyle(20);
        jetShapePbPb[1][iAsymmetry][iCentrality][iTrackPt]->SetMarkerColor(veryNiceColors[1]);
        jetShapePbPb[1][iAsymmetry][iCentrality][iTrackPt]->Draw("same");
        if(iTrackPt == nTrackPtBins) cout << "PbPb integral: " << jetShapePbPb[1][iAsymmetry][iCentrality][iTrackPt]->Integral(1, jetShapePbPb[1][iAsymmetry][iCentrality][iTrackPt]->FindBin(0.99), "width") << endl;
        
        legend = new TLegend(0.45,0.5,0.9,0.8);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0,Form("%s jet %s",labels[iJetTrack],centralityString.Data()),"");
        legend->AddEntry((TObject*) 0,trackString.Data(),"");
        if(iAsymmetry < nAsymmetryBins) legend->AddEntry((TObject*) 0,asymmetryString.Data(),"");
        if(mcLabel){
          legend->AddEntry(jetShapePbPb[0][iAsymmetry][iCentrality][iTrackPt],"Jet100","l");
          legend->AddEntry(jetShapePbPb[1][iAsymmetry][iCentrality][iTrackPt],"Jet80","l");
        } else {
          legend->AddEntry(jetShapePbPb[0][iAsymmetry][iCentrality][iTrackPt],"Nominal","l");
          legend->AddEntry(jetShapePbPb[1][iAsymmetry][iCentrality][iTrackPt],"Smeared","l");
        }
        legend->Draw();
        
        // Draw the ratio to the lower pad
        drawer->SetLogY(false);
        
        jetShapeRatio[iAsymmetry][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
        jetShapeRatio[iAsymmetry][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0.9,1.1);
        jetShapeRatio[iAsymmetry][iCentrality][iTrackPt]->SetMarkerStyle(20);
        jetShapeRatio[iAsymmetry][iCentrality][iTrackPt]->SetMarkerColor(veryNiceColors[1]);
        jetShapeRatio[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(veryNiceColors[1]);
        if(mcLabel){
          drawer->DrawHistogramToLowerPad(jetShapeRatio[iAsymmetry][iCentrality][iTrackPt],"#Deltar","Hydjet / Pythia", " ");
        } else {
          drawer->DrawHistogramToLowerPad(jetShapeRatio[iAsymmetry][iCentrality][iTrackPt],"#Deltar","Nominal / Smear", " ");
        }
        oneLine->Draw();
        
        //jetShapeRatioUncertaintyAll[iAsymmetry][iCentrality][iTrackPt]->SetBarWidth(0.3);
        //jetShapeRatioUncertaintyAll[iAsymmetry][iCentrality][iTrackPt]->SetBarOffset(0.05);
        //jetShapeRatioUncertaintyAll[iAsymmetry][iCentrality][iTrackPt]->SetFillColorAlpha(kRed,0.3);
        //jetShapeRatioUncertaintyAll[iAsymmetry][iCentrality][iTrackPt]->Draw("same,E2");
        
        //jetShapeRatioUncertaintyRelevant[iAsymmetry][iCentrality][iTrackPt]->SetBarWidth(0.3);
        //jetShapeRatioUncertaintyRelevant[iAsymmetry][iCentrality][iTrackPt]->SetBarOffset(0.35);
        //jetShapeRatioUncertaintyRelevant[iAsymmetry][iCentrality][iTrackPt]->SetFillColorAlpha(kBlue,0.3);
        //jetShapeRatioUncertaintyRelevant[iAsymmetry][iCentrality][iTrackPt]->Draw("same,E2");
        
        //jetShapeRatioUncertaintyRelevantPlus2[iAsymmetry][iCentrality][iTrackPt]->SetBarWidth(0.3);
        //jetShapeRatioUncertaintyRelevantPlus2[iAsymmetry][iCentrality][iTrackPt]->SetBarOffset(0.65);
        //jetShapeRatioUncertaintyRelevantPlus2[iAsymmetry][iCentrality][iTrackPt]->SetFillColorAlpha(kGreen+4,0.3);
        //jetShapeRatioUncertaintyRelevantPlus2[iAsymmetry][iCentrality][iTrackPt]->Draw("same,E2");
        
        //ratioLegend = new TLegend(0.25,0.85,0.45,0.95);  // Leading
        //ratioLegend = new TLegend(0.35,0.85,0.55,0.95);  // Subleading
        //ratioLegend->SetFillStyle(0);ratioLegend->SetBorderSize(0);ratioLegend->SetTextSize(0.05);ratioLegend->SetTextFont(62);
        //ratioLegend->AddEntry(jetShapeRatioUncertaintyAll[iAsymmetry][iCentrality][iTrackPt], "All syst", "F");
        //ratioLegend->Draw();
        
        //ratioLegend2 = new TLegend(0.55,0.85,0.85,0.95);
        //ratioLegend2 = new TLegend(0.25,0.85,0.55,0.95); // Centered
        //ratioLegend2->SetFillStyle(0);ratioLegend2->SetBorderSize(0);ratioLegend2->SetTextSize(0.05);ratioLegend2->SetTextFont(62);
        //ratioLegend2->AddEntry(jetShapeRatioUncertaintyRelevant[iAsymmetry][iCentrality][iTrackPt], "Relevant uncertainties", "F");
        //ratioLegend2->Draw();
        
        //ratioLegend3 = new TLegend(0.65,0.85,0.95,0.95);
        //ratioLegend3->SetFillStyle(0);ratioLegend3->SetBorderSize(0);ratioLegend3->SetTextSize(0.05);ratioLegend3->SetTextFont(62);
        //ratioLegend3->AddEntry(jetShapeRatioUncertaintyRelevantPlus2[iAsymmetry][iCentrality][iTrackPt], "Add trigger bias", "F");
        //ratioLegend3->Draw();
        
        // Save the figures into a file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/jetShapeComparisonSmearing%s%s%s%s.png", labels[iJetTrack], compactAsymmetryString.Data(), compactCentralityString.Data(), compactTrackString.Data()));
          //gPad->GetCanvas()->SaveAs(Form("figures/jetShapeComparisonErrorSets%s%s%s%s.pdf", labels[iJetTrack], compactAsymmetryString.Data(), compactCentralityString.Data(), compactTrackString.Data()));
        }
        
      } // Track pT loop
    } // Centrality loop
  } // Asymmetry loop
  
  int ratioColors[] = {kBlack, kBlue, kRed, kGreen+2};
  TH1D *hZeroFiller[4];
  
  // Extra figures for putting different ratios in the same plot
  for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
    
    // Setup asymmetry strings
    if(iAsymmetry < nAsymmetryBins){
      asymmetryString = Form("%.1f < x_{j} < %.1f", asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1]);
      compactAsymmetryString = Form("_A=%.1f-%.1f", asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1]);
      compactAsymmetryString.ReplaceAll(".","v");
    } else {
      asymmetryString = "";
      compactAsymmetryString = "";
    }
    
    drawer->Reset();
    
    jetShapeRatio[iAsymmetry][firstDrawnCentralityBin][nTrackPtBins]->SetMarkerColor(ratioColors[0]);
    jetShapeRatio[iAsymmetry][firstDrawnCentralityBin][nTrackPtBins]->SetLineColor(ratioColors[0]);
    jetShapeRatio[iAsymmetry][firstDrawnCentralityBin][nTrackPtBins]->GetYaxis()->SetRangeUser(0.92,1.1);
    
    drawer->DrawHistogram(jetShapeRatio[iAsymmetry][firstDrawnCentralityBin][nTrackPtBins], "#Deltar", "Nominal / Smeared", " ");
    
    oneLine->Draw();
    
    legend = new TLegend(0.6,0.6,0.9,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry(jetShapeRatio[iAsymmetry][firstDrawnCentralityBin][nTrackPtBins],"0-10 %","l");
    
    ratioLegend = new TLegend(0.18,0.81,0.48,0.91);
    ratioLegend->SetFillStyle(0);ratioLegend->SetBorderSize(0);ratioLegend->SetTextSize(0.05);ratioLegend->SetTextFont(62);
    ratioLegend->AddEntry((TObject*) 0,"PbPb data, all x_{j}","");
    
    for(int iCentrality = firstDrawnCentralityBin+1; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("%.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
      compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
      
      jetShapeRatio[iAsymmetry][iCentrality][nTrackPtBins]->SetMarkerColor(ratioColors[iCentrality]);
      jetShapeRatio[iAsymmetry][iCentrality][nTrackPtBins]->SetLineColor(ratioColors[iCentrality]);
      
      jetShapeRatio[iAsymmetry][iCentrality][nTrackPtBins]->Draw("same");
      legend->AddEntry(jetShapeRatio[iAsymmetry][iCentrality][nTrackPtBins],centralityString,"l");
      
    } // Centrality loop
    
    legend->Draw();
    ratioLegend->Draw();
    
  } // Asymmetry loop
  
  // Extra figures for putting different ratios in the same plot
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    
    centralityString = Form("%.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    
    drawer->Reset();
    
    jetShapeRatio[nAsymmetryBins][iCentrality][nTrackPtBins]->SetMarkerColor(ratioColors[0]);
    jetShapeRatio[nAsymmetryBins][iCentrality][nTrackPtBins]->SetLineColor(ratioColors[0]);
    jetShapeRatio[nAsymmetryBins][iCentrality][nTrackPtBins]->GetYaxis()->SetRangeUser(0.5,1.5);
    
    // Prepare hiding large error points
    hZeroFiller[3] = (TH1D*) jetShapeRatio[nAsymmetryBins][iCentrality][nTrackPtBins]->Clone(Form("zeroFiller3%d",iCentrality));
    hZeroFiller[3]->SetMarkerStyle(kOpenCircle);
    hZeroFiller[3]->GetYaxis()->SetRangeUser(-1,3);
    
    for(int iBin = 1; iBin < hZeroFiller[3]->GetNbinsX(); iBin++){
      hZeroFiller[3]->SetBinError(iBin,0);
      hZeroFiller[3]->SetBinContent(iBin,0);
      //if(iBin > 12) {
      //  hZeroFiller[3]->SetBinContent(iBin,1);
      //  jetShapeRatio[nAsymmetryBins][iCentrality][nTrackPtBins]->SetBinContent(iBin,0);
      //}
    }
    
    drawer->DrawHistogram(jetShapeRatio[nAsymmetryBins][iCentrality][nTrackPtBins], "#Deltar", "Reco / Gen", " ");
    hZeroFiller[3]->Draw("p,same");
    
    oneLine->Draw();
    
    //legend = new TLegend(0.6,0.6,0.9,0.9);
    legend = new TLegend(0.19,0.19,0.49,0.49);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry(jetShapeRatio[nAsymmetryBins][iCentrality][nTrackPtBins],"All x_{j}","lp");
    
    ratioLegend = new TLegend(0.17,0.81,0.48,0.91);
    ratioLegend->SetFillStyle(0);ratioLegend->SetBorderSize(0);ratioLegend->SetTextSize(0.05);ratioLegend->SetTextFont(62);
    ratioLegend->AddEntry((TObject*) 0,Form("Analysis closure, subleading jet, %s", centralityString.Data()),"");
    
    for(int iAsymmetry = 0; iAsymmetry <= 2; iAsymmetry++){
      
      // Setup asymmetry strings
      if(iAsymmetry < nAsymmetryBins){
        asymmetryString = Form("%.1f < x_{j} < %.1f", asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1]);
        compactAsymmetryString = Form("_A=%.1f-%.1f", asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1]);
        compactAsymmetryString.ReplaceAll(".","v");
      } else {
        asymmetryString = "";
        compactAsymmetryString = "";
      }
      
      jetShapeRatio[iAsymmetry][iCentrality][nTrackPtBins]->SetMarkerColor(ratioColors[iAsymmetry+1]);
      jetShapeRatio[iAsymmetry][iCentrality][nTrackPtBins]->SetLineColor(ratioColors[iAsymmetry+1]);
      
      // Prepare hiding large error points
      hZeroFiller[iAsymmetry] = (TH1D*) jetShapeRatio[iAsymmetry][iCentrality][nTrackPtBins]->Clone(Form("zeroFiller%d%d", iAsymmetry, iCentrality));
      hZeroFiller[iAsymmetry]->SetMarkerStyle(kOpenCircle);
      hZeroFiller[iAsymmetry]->GetYaxis()->SetRangeUser(-1,3);
      
      for(int iBin = 1; iBin < hZeroFiller[iAsymmetry]->GetNbinsX(); iBin++){
        hZeroFiller[iAsymmetry]->SetBinError(iBin,0);
        hZeroFiller[iAsymmetry]->SetBinContent(iBin,0);
        //if(iBin > 12) {
        //  hZeroFiller[iAsymmetry]->SetBinContent(iBin,1);
        //  jetShapeRatio[iAsymmetry][iCentrality][nTrackPtBins]->SetBinContent(iBin,0);
        //}
      }
      
      jetShapeRatio[iAsymmetry][iCentrality][nTrackPtBins]->Draw("same");
      hZeroFiller[iAsymmetry]->Draw("p,same");
      legend->AddEntry(jetShapeRatio[iAsymmetry][iCentrality][nTrackPtBins],asymmetryString,"lp");
      
    } // Asymmetry loop
    
    legend->Draw();
    ratioLegend->Draw();
    
  } //  Centrality loop

}
