#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "JDrawer.h"

/*
 * Macro for printing the relative uncertainties in latex slides
 */
void printRelativeUncertainty(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString dataFileName = "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_allCorrections_processed_2020-02-17.root"; // Compare also with uncorrected data
  // data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root
  TString ppFileName = "data/ppData2017_highForest_pfJets_20EveMixed_xjBins_wtaAxis_allCorrections_processed_2020-02-04.root";
  TString uncertaintyFileName[2];
  TString ppUncertaintyFileName[2];
  uncertaintyFileName[0] = "uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_manualFluctuationReduction_smoothedPairBackground_2020-05-22.root";
  ppUncertaintyFileName[0] = "uncertainties/systematicUncertaintyForPp_20eveMix_xjBins_fixJES_2020-02-03.root";
  uncertaintyFileName[1] = "uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_addNewSources_smoothedPairBackground_2020-05-18.root";
  ppUncertaintyFileName[1] = "uncertainties/systematicUncertaintyForPp_20eveMix_newJESestimate_2020-01-13.root";
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_manualFluctuationReduction_smoothedPairBackground_2020-05-22.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_addNewSources_smoothedPairBackground_2020-05-18.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_finalTuning_smoothedPairBackground_2020-03-09.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_includeTrackDeltaR_2020-01-27.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_newSpilloverWithSmoothedBackground_newJES_tunedSeagull_smoothedPairBackground_2020-01-23.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_newSpilloverJES_seagullTuningProcess_2020-01-15.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_newSpilloverJES_tunedSeagull_smoothedPairBackground_2020-01-21.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_newSpilloverJEStunedSeagull_noRebinPairBg_2020-01-14.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_newSpilloverJEStunedSeagull_forComparison_2020-01-14.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_oldJES_15percentSpill10Jff_2019-10-17.root
  
  bool drawComparisonUncertaintyFile = false;
  bool smoothHistograms = true;
  int nFiles = drawComparisonUncertaintyFile ? 2 : 1;
  
  const char* legendLabels[] = {"Dijet analysis", "New estimate", "lul", "Toinen lul"};  // Labels referring to different uncertainty files
  
  TString inclusiveFileName = "uncertainties/inclusiveAnalysis/js_AllSource_syst_err.root";
  
  bool printSlides = false;  // Print slides showing the R-integrated uncertainty in each pT bin
  bool combineTracking = false; // Combine all tracking related uncertainty to one when printing the table
  bool drawUncertaintySourceComparison = true; // Draw all uncertainty sources as a function of R in each pT bin
  bool drawUncertaintySystemComparison = false; // Draw single uncertainty source for all systems as a function of R in each pT bin
  bool drawComparisonToInclusive = false;        // Draw comparison to systematic uncertainty histograms from inclusive analysis
  bool saveFigures = false;   // Save the drawn figures to file
  
  // ==================================================================
  // ====================== Configuration done ========================
  // ==================================================================
  
  // Open the input files
  TFile *dataFile = TFile::Open(dataFileName);
  TFile *uncertaintyFile[2];
  for(int iFile = 0; iFile < nFiles; iFile++) uncertaintyFile[iFile] = TFile::Open(uncertaintyFileName[iFile]);
  TFile *ppFile = TFile::Open(ppFileName);
  TFile *ppUncertaintyFile[2];
  for(int iFile = 0; iFile < nFiles; iFile++) ppUncertaintyFile[iFile] = TFile::Open(ppUncertaintyFileName[iFile]);
  
  TFile *inclusiveFile;
  if(drawComparisonToInclusive) inclusiveFile = TFile::Open(inclusiveFileName);
  
  // Make readers to read the files
  DijetHistogramManager *dataManager = new DijetHistogramManager(dataFile);
  dataManager->SetLoadTrackLeadingJetCorrelationsPtWeighted(true);
  dataManager->SetLoadTrackSubleadingJetCorrelationsPtWeighted(true);
  dataManager->SetAsymmetryBinRange(0,3);
  dataManager->LoadProcessedHistograms();
  
  DijetHistogramManager *ppManager = new DijetHistogramManager(ppFile);
  ppManager->SetLoadTrackLeadingJetCorrelationsPtWeighted(true);
  ppManager->SetLoadTrackSubleadingJetCorrelationsPtWeighted(true);
  ppManager->SetAsymmetryBinRange(0,3);
  ppManager->LoadProcessedHistograms();
  
  JffCorrector *uncertaintyManager[2];
  for(int iFile = 0; iFile < nFiles; iFile++){
    uncertaintyManager[iFile] = new JffCorrector();
    uncertaintyManager[iFile]->ReadSystematicFile(uncertaintyFile[iFile]);
    uncertaintyManager[iFile]->SetUncertaintySmooth(smoothHistograms);
  }
  
  JffCorrector *ppUncertaintyManager[2];
  for(int iFile = 0; iFile < nFiles; iFile++){
    ppUncertaintyManager[iFile] = new JffCorrector();
    ppUncertaintyManager[iFile]->ReadSystematicFile(ppUncertaintyFile[iFile]);
    ppUncertaintyManager[iFile]->SetUncertaintySmooth(smoothHistograms);
  }
  
  const int nCentralityBins = dataManager->GetNCentralityBins();
  const int nTrackPtBins = dataManager->GetNTrackPtBins();
  const int nAsymmetryBins = dataManager->GetNAsymmetryBins();
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  // Select drawn bins
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins-1;
  
  int firstDrawnAsymmetryBin = 3;
  int lastDrawnAsymmetryBin = 3;
    
  // Select which uncertainty sources to draw
  bool drawAllUncertainites = true;
  bool drawSpilloverUncertainty = false;
  bool drawJffUncertainty = false;
  bool drawJetEnergyScaleUncertainty = false;
  bool drawTrackingEfficiencyUncertainty = false;
  bool drawResidualTrackinguncertainty = false;
  bool drawTrackingDeltaRUncertainty = false;
  bool drawPairAcceptanceUncertainty = false;
  bool drawBackgroundSubtractionUncertainty = true;
  bool drawTotalUncertainty = false;
  
  bool drawUncertainty[JffCorrector::knUncertaintySources] = {drawSpilloverUncertainty, drawJffUncertainty, drawJetEnergyScaleUncertainty, drawTrackingEfficiencyUncertainty, drawResidualTrackinguncertainty, drawTrackingDeltaRUncertainty, drawPairAcceptanceUncertainty, drawBackgroundSubtractionUncertainty, drawTotalUncertainty};
  
  if(drawAllUncertainites){
    for(int iUncertainty = 0; iUncertainty < JffCorrector::knUncertaintySources; iUncertainty++){
      drawUncertainty[iUncertainty] = true;
    }
  }
  
  // If we are drawing system comparison or printing slides, we must use all centrality bins
  if(drawUncertaintySystemComparison || printSlides){
    firstDrawnCentralityBin = 0;
    lastDrawnCentralityBin = nCentralityBins-1;
  }
  
  // Define arrays for the jet shapes
  TH1D *dataShape[2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins+1];
  TH1D *uncertaintyShape[2][2][JffCorrector::knUncertaintySources][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins+1];
  
  // Define arrays for yield
  double dataYield[2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins+1];
  double uncertaintyYield[2][2][JffCorrector::knUncertaintySources][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins+1];
  
  // Define histograms for uncertainties from inclusive analysis
  // There are four sources of uncertainty in the inclusive files:
  //   (1) Spillover
  //   (2) JFF
  //   (3) Background subtraction + pair acceptance
  //   (4) All the rest
  TH1D *inclusiveUncertainty[5][nCentralityBins][nTrackPtBins]; // First bin: uncertainty source
  TH1D *helperHistogram;
  TString inclusiveSources[] = {"rel_spill_err","jff_err","rel_bkg_err","rel_rela_err"};
  TString inclusiveCentrality[] = {"Cent0_Cent10","Cent10_Cent30","Cent30_Cent50","Cent50_Cent100"};
  TString inclusiveTrackPt[] = {"TrkPt07_TrkPt1","TrkPt1_TrkPt2","TrkPt2_TrkPt3","TrkPt3_TrkPt4","TrkPt4_TrkPt8","TrkPt8_TrkPt12","TrkPt12_TrkPt16","TrkPt16_TrkPt20","TrkPt20_TrkPt300"};
  const int nInclusiveTrackPtBins = 9;
  double oldContent, addedContent, newContent;
  
  // Read PbPb and pp histograms from files
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    for(int iTrackPt = firstDrawnTrackPtBin; iTrackPt <= lastDrawnTrackPtBin; iTrackPt++){
      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
        
        // Different data reader for pp and PbPb files
        if(iCentrality == nCentralityBins){
          dataShape[0][iAsymmetry][iCentrality][iTrackPt] = ppManager->GetHistogramJetShape(DijetHistogramManager::kJetShape, DijetHistogramManager::kPtWeightedTrackLeadingJet, iAsymmetry, 0, iTrackPt);
          dataShape[1][iAsymmetry][iCentrality][iTrackPt] = ppManager->GetHistogramJetShape(DijetHistogramManager::kJetShape, DijetHistogramManager::kPtWeightedTrackSubleadingJet, iAsymmetry, 0, iTrackPt);
        } else {
          dataShape[0][iAsymmetry][iCentrality][iTrackPt] = dataManager->GetHistogramJetShape(DijetHistogramManager::kJetShape, DijetHistogramManager::kPtWeightedTrackLeadingJet, iAsymmetry, iCentrality, iTrackPt);
          dataShape[1][iAsymmetry][iCentrality][iTrackPt] = dataManager->GetHistogramJetShape(DijetHistogramManager::kJetShape, DijetHistogramManager::kPtWeightedTrackSubleadingJet, iAsymmetry, iCentrality, iTrackPt);
        }
        
        dataYield[0][iAsymmetry][iCentrality][iTrackPt] = dataShape[0][iAsymmetry][iCentrality][iTrackPt]->Integral(1, dataShape[0][iAsymmetry][iCentrality][iTrackPt]->FindBin(0.99), "width");
        dataYield[1][iAsymmetry][iCentrality][iTrackPt] = dataShape[1][iAsymmetry][iCentrality][iTrackPt]->Integral(1, dataShape[1][iAsymmetry][iCentrality][iTrackPt]->FindBin(0.99), "width");
        
        for(int iFile = 0; iFile < nFiles; iFile++){
          for(int iUncertainty = 0; iUncertainty < JffCorrector::knUncertaintySources; iUncertainty++){
            
            // Initialize the uncertainty histograms to NULL
            uncertaintyShape[iFile][0][iUncertainty][iAsymmetry][iCentrality][iTrackPt] = NULL;
            
            // Different reader for pp and PbPb uncertainty files
            if(iCentrality == nCentralityBins){
              uncertaintyShape[iFile][0][iUncertainty][iAsymmetry][iCentrality][iTrackPt] = ppUncertaintyManager[iFile]->GetJetShapeSystematicUncertainty(DijetHistogramManager::kPtWeightedTrackLeadingJet, 0, iTrackPt, iAsymmetry, iUncertainty);
              uncertaintyShape[iFile][1][iUncertainty][iAsymmetry][iCentrality][iTrackPt] = ppUncertaintyManager[iFile]->GetJetShapeSystematicUncertainty(DijetHistogramManager::kPtWeightedTrackSubleadingJet, 0, iTrackPt, iAsymmetry, iUncertainty);
            } else {
              uncertaintyShape[iFile][0][iUncertainty][iAsymmetry][iCentrality][iTrackPt] = uncertaintyManager[iFile]->GetJetShapeSystematicUncertainty(DijetHistogramManager::kPtWeightedTrackLeadingJet, iCentrality, iTrackPt, iAsymmetry, iUncertainty);
              uncertaintyShape[iFile][1][iUncertainty][iAsymmetry][iCentrality][iTrackPt] = uncertaintyManager[iFile]->GetJetShapeSystematicUncertainty(DijetHistogramManager::kPtWeightedTrackSubleadingJet, iCentrality, iTrackPt, iAsymmetry, iUncertainty);
            }
            
            if(uncertaintyShape[iFile][0][iUncertainty][iAsymmetry][iCentrality][iTrackPt]){
              uncertaintyYield[iFile][0][iUncertainty][iAsymmetry][iCentrality][iTrackPt] = uncertaintyShape[iFile][0][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->Integral(1, uncertaintyShape[iFile][0][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->FindBin(0.99), "width");
              uncertaintyYield[iFile][1][iUncertainty][iAsymmetry][iCentrality][iTrackPt] = uncertaintyShape[iFile][1][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->Integral(1, uncertaintyShape[iFile][1][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->FindBin(0.99), "width");
            } else {
              uncertaintyYield[iFile][0][iUncertainty][iAsymmetry][iCentrality][iTrackPt] = 0;
              uncertaintyYield[iFile][1][iUncertainty][iAsymmetry][iCentrality][iTrackPt] = 0;
            }
          } // Uncertainty type loop
        } // Uncertainty file loop
      } // Asymmetry loop
    } // Track pT loop
    
    // Load the histograms for inclusive analysis
    if(drawComparisonToInclusive){
      
      // No pp for inclusive analysis
      if(iCentrality == nCentralityBins) continue;
      
      for(int iTrackPt = 0; iTrackPt < nInclusiveTrackPtBins; iTrackPt++){
        for(int iUncertainty = 0; iUncertainty < 5; iUncertainty++){
          
          // The total uncertainty has different naming convention from individual sources in the uncertainty file
          if(iUncertainty == 4){
            helperHistogram = (TH1D*) inclusiveFile->Get(Form("js_dr_Pb_Syst_Error_%d_%d", iTrackPt, iCentrality));
          } else {
            helperHistogram = (TH1D*) inclusiveFile->Get(Form("dR_Syst_PbPb_%s_%s_%s", inclusiveCentrality[iCentrality].Data(), inclusiveTrackPt[iTrackPt].Data(), inclusiveSources[iUncertainty].Data()));
          }
          
          if(iTrackPt < nTrackPtBins){
            inclusiveUncertainty[iUncertainty][iCentrality][iTrackPt] = helperHistogram;
          } else {

            // Add histograms in quadrature for high pT bins to get same binning as in this analysis
            for(int iBin = 1; iBin <= helperHistogram->GetNbinsX(); iBin++){
              oldContent = inclusiveUncertainty[iUncertainty][iCentrality][nTrackPtBins-1]->GetBinContent(iBin);
              addedContent = helperHistogram->GetBinContent(iBin);
              newContent = TMath::Sqrt(oldContent*oldContent+addedContent*addedContent);
              inclusiveUncertainty[iUncertainty][iCentrality][nTrackPtBins-1]->SetBinContent(iBin,newContent);
            }
            
          }
        } // Inclusive uncertainty loop
      } // Inclusive track pT loop
      
    } // Reading inclusive histograms from the file
  } // Centrality loop
  
  // Helper variables for printing and drawing
  char namer[100];
  const char* jetType[2] = {"Leading","Subleading"};
  
  // Print a slide with uncertainties for each source and each centrality
  if(printSlides){
    for(int iJetType = 0; iJetType < 2; iJetType++){
      for(int iTrackPt = firstDrawnTrackPtBin; iTrackPt <= lastDrawnTrackPtBin; iTrackPt++){
        
        if(iTrackPt == nTrackPtBins){
          sprintf(namer,"\\frametitle{%s jet shape error, $p_{\\mathrm{T}}$} integrated}",jetType[iJetType]);
        } else {
          sprintf(namer,"\\frametitle{%s jet shape error, $%.1f < p_{\\mathrm{T}} < %.1f$}",jetType[iJetType], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
        }
        
        cout << endl;
        cout << "\\begin{frame}" << endl;
        cout << namer << endl;
        cout << "\\begin{center}" << endl;
        cout << "  \\begin{tabular}{ccccc}" << endl;
        cout << "    \\toprule" << endl;
        //cout << "    Source & C: 0-10 & C: 10-30 & C: 30-50 & C: 50-90 & pp \\\\" << endl;
        cout << "    Source & C: 0-10 & C: 10-30 & C: 30-50 & C: 50-90 \\\\" << endl;
        cout << "    \\midrule" << endl;
        
        // Set the correct precision for printing floating point numbers
        cout << fixed << setprecision(3);
        
        // Combine the tracking uncertainties kTrackingEfficiency, kResidualTracking, kTrackingDeltaR
        if(combineTracking){
          for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
            uncertaintyYield[0][iJetType][JffCorrector::kTrackingEfficiency][nAsymmetryBins][iCentrality][iTrackPt] = TMath::Sqrt(TMath::Power(uncertaintyYield[0][iJetType][JffCorrector::kTrackingEfficiency][nAsymmetryBins][iCentrality][iTrackPt],2) + TMath::Power(uncertaintyYield[0][iJetType][JffCorrector::kResidualTracking][nAsymmetryBins][iCentrality][iTrackPt],2) + TMath::Power(uncertaintyYield[0][iJetType][JffCorrector::kTrackingDeltaR][nAsymmetryBins][iCentrality][iTrackPt],2));
          }
        }
        
        // Print one line for each track pT bin
        for(int iUncertainty = 0; iUncertainty < JffCorrector::knUncertaintySources; iUncertainty++){
          if(combineTracking && (iUncertainty == JffCorrector::kResidualTracking || iUncertainty == JffCorrector::kTrackingDeltaR)) continue;
          cout << uncertaintyManager[0]->GetUncertaintyName(iUncertainty);
          for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
            cout << " & $" << uncertaintyYield[0][iJetType][iUncertainty][nAsymmetryBins][iCentrality][iTrackPt] / dataYield[iJetType][nAsymmetryBins][iCentrality][iTrackPt] << "$";
          }
          cout << " \\\\" << endl;
        }
        
        cout << "    \\bottomrule" << endl;
        cout << "  \\end{tabular}" << endl;
        cout << "\\end{center}" << endl;
        cout << "\\end{frame}" << endl;
        cout << endl;
        
      } // Loop over track pT bins
    } // Loop over leading and subleading jets
  } // Printing slides
  
  if(drawUncertaintySourceComparison || drawUncertaintySystemComparison || drawComparisonToInclusive){
    
    JDrawer *drawer = new JDrawer();
    TLegend *legend;
    int colors[] = {kBlue, kRed, kMagenta, kCyan, kGreen+3, kViolet+7, kOrange+7, kSpring, kGray};
    int trueColors[] = {kBlack, kBlue, kMagenta, kCyan, kGreen+3, kViolet+7, kOrange+7, kSpring, kGray};
    double maxValue = 0;
    double localMaximum = 0;
    double histogramMaximum = 0;
    TString centralityString;
    TString compactCentralityString;
    TString trackPtString;
    TString compactTrackPtString;
    TString jetShapeString;
    TString asymmetryString;
    TString compactAsymmetryString;
        
    // Draw all uncertainty sources as a function of R in each pT bin
    if(drawUncertaintySourceComparison){
      
      for(int iJetType = 0; iJetType < 2; iJetType++){
        
        jetShapeString = Form("%s jet shape", jetType[iJetType]);
        
        for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
          
          if(iCentrality == nCentralityBins){
            centralityString = "pp";
            compactCentralityString = "_pp";
          } else {
            centralityString = Form("C = %.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
            compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
          }
          
          for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
            
            if(iAsymmetry < nAsymmetryBins){
              asymmetryString = Form(", %.1f < x_{j} < %.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
              compactAsymmetryString = Form("_A=%.1f-%.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
              compactAsymmetryString.ReplaceAll(".","v");
            } else {
              asymmetryString = "";
              compactAsymmetryString = "";
            }
            
            for(int iTrackPt = firstDrawnTrackPtBin; iTrackPt <= lastDrawnTrackPtBin; iTrackPt++){
              
              trackPtString = Form("Track pT: %.1f-%.1f GeV", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
              compactTrackPtString = Form("_pT=%.1f-%.1f", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
              compactTrackPtString.ReplaceAll(".","v");
              
              // Draw the first uncertainty source
              uncertaintyShape[0][iJetType][JffCorrector::kTotal][iAsymmetry][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
              maxValue = uncertaintyShape[0][iJetType][JffCorrector::kTotal][iAsymmetry][iCentrality][iTrackPt]->GetMaximum();
              uncertaintyShape[0][iJetType][JffCorrector::kTotal][iAsymmetry][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0,maxValue*1.3);
              uncertaintyShape[0][iJetType][JffCorrector::kTotal][iAsymmetry][iCentrality][iTrackPt]->SetLineColor(kBlack);
              drawer->DrawHistogram(uncertaintyShape[0][iJetType][JffCorrector::kTotal][iAsymmetry][iCentrality][iTrackPt], "#Deltar", "P(#Deltar) error", " ", "");
              
              // Make a legend for the titles
              legend = new TLegend(0.12,0.77,0.4,0.92);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
              legend->AddEntry((TObject*) 0,jetShapeString.Data(),"");
              legend->AddEntry((TObject*) 0,(centralityString+asymmetryString).Data(),"");
              legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
              legend->Draw();
              
              // Create a new legend for uncertainties
              legend = new TLegend(0.5,0.6,0.9,0.92);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
              
              // Add the first uncertainty source to legend
              legend->AddEntry(uncertaintyShape[0][iJetType][JffCorrector::kTotal][iAsymmetry][iCentrality][iTrackPt], uncertaintyManager[0]->GetUncertaintyName(JffCorrector::kTotal), "l");
              
              // Loop over all the other uncertainty sources
              for(int iUncertainty = 0; iUncertainty < JffCorrector::kTotal; iUncertainty++){
                if(!drawUncertainty[iUncertainty]) continue;
                uncertaintyShape[0][iJetType][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->SetLineColor(colors[iUncertainty]);
                uncertaintyShape[0][iJetType][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->Draw("same");
                legend->AddEntry(uncertaintyShape[0][iJetType][iUncertainty][iAsymmetry][iCentrality][iTrackPt], uncertaintyManager[0]->GetUncertaintyName(iUncertainty), "l");
              }
              
              // Draw the legend to the figure
              legend->Draw();
              
              // Save the figures into a file
              if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/uncertaintySourceBreakdown%sJet%s%s%s.pdf", jetType[iJetType], compactAsymmetryString.Data(), compactCentralityString.Data(), compactTrackPtString.Data()));
              }
              
            } // Track pT loop
          } // Asymmetry loop
        } // Centrality loop
      } // Jet type loop (leading/subleading)
    } // Uncertainty source comparison
    
    // Draw single uncertainty source for all systems as a function of R in each pT bin
    if(drawUncertaintySystemComparison){
      
      for(int iJetType = 0; iJetType < 2; iJetType++){
        
        jetShapeString = Form("%s jet shape", jetType[iJetType]);
        
        for(int iUncertainty = 0; iUncertainty < JffCorrector::knUncertaintySources; iUncertainty++){
          if(!drawUncertainty[iUncertainty]) continue;
          
          for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
            
            if(iAsymmetry < nAsymmetryBins){
              asymmetryString = Form(", %.1f < x_{j} < %.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
              compactAsymmetryString = Form("_A=%.1f-%.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
              compactAsymmetryString.ReplaceAll(".","v");
            } else {
              asymmetryString = "";
              compactAsymmetryString = "";
            }
            
            for(int iTrackPt = firstDrawnTrackPtBin; iTrackPt <= lastDrawnTrackPtBin; iTrackPt++){
              
              trackPtString = Form("Track pT: %.1f-%.1f GeV", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
              compactTrackPtString = Form("_pT=%.1f-%.1f", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
              compactTrackPtString.ReplaceAll(".","v");
              
              // First, find the maximum value for uncertainty among all centrality bins
              maxValue = 0;
              for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
                uncertaintyShape[0][iJetType][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
                histogramMaximum = uncertaintyShape[0][iJetType][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->GetMaximum();
                if(histogramMaximum > maxValue) maxValue = histogramMaximum;
              }
              
              // Draw the most central bin
              uncertaintyShape[0][iJetType][iUncertainty][iAsymmetry][0][iTrackPt]->GetYaxis()->SetRangeUser(0,maxValue*1.3);
              uncertaintyShape[0][iJetType][iUncertainty][iAsymmetry][0][iTrackPt]->SetLineColor(kBlack);
              drawer->DrawHistogram(uncertaintyShape[0][iJetType][iUncertainty][iAsymmetry][0][iTrackPt], "#Deltar", "P(#Deltar) error", " ", "");
              
              // Make a legend for the titles
              legend = new TLegend(0.12,0.77,0.4,0.92);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
              legend->AddEntry((TObject*) 0,(jetShapeString+asymmetryString).Data(),"");
              legend->AddEntry((TObject*) 0,uncertaintyManager[0]->GetUncertaintyName(iUncertainty),"");
              legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
              legend->Draw();
              
              // Create a new legend systems
              legend = new TLegend(0.7,0.65,0.9,0.92);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
              
              // Add the most central bin to the legend
              legend->AddEntry(uncertaintyShape[0][iJetType][iUncertainty][iAsymmetry][0][iTrackPt], "C = 0-10 %", "l");
              
              // Loop over all the other centrality bins and pp
              for(int iCentrality = 1; iCentrality <= nCentralityBins; iCentrality++){
                
                if(iCentrality == nCentralityBins){
                  centralityString = "pp";
                } else {
                  centralityString = Form("C = %.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
                }
                
                uncertaintyShape[0][iJetType][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->SetLineColor(colors[iCentrality-1]);
                uncertaintyShape[0][iJetType][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->Draw("same");
                legend->AddEntry(uncertaintyShape[0][iJetType][iUncertainty][iAsymmetry][iCentrality][iTrackPt], centralityString.Data(), "l");
              } // Centrality bin loop
              
              // Draw the legend to the figure
              legend->Draw();
              
              // Save the figures into a file
              if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/uncertaintySystemBreakdown%sJet_%s%s%s.pdf", jetType[iJetType], uncertaintyManager[0]->GetUncertaintyName(iUncertainty).Data(), compactAsymmetryString.Data(), compactTrackPtString.Data()));
              }
              
            } // Track pT loop
          } // Asymmetry loop
        } // Uncertainty source loop
      } // Jet type loop (leading/subleading)
    } // Uncertainty system comparison
    
    // Draw comparison of leading jet uncertainties and uncertainties from the inclusive jet analysis source by source
    if(drawComparisonToInclusive){
      
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        
        // No pp available for inclusive
        if(iCentrality == nCentralityBins) continue;
        
        centralityString = Form("C = %.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
        compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
        
        for(int iTrackPt = firstDrawnTrackPtBin; iTrackPt <= lastDrawnTrackPtBin; iTrackPt++){
          
          trackPtString = Form("Track pT: %.1f-%.1f GeV", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
          compactTrackPtString = Form("_pT=%.1f-%.1f", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
          compactTrackPtString.ReplaceAll(".","v");
          
          // ======================================================
          // == Compare with inclusive on background fluctuation ==
          // ======================================================
          
          maxValue = 0;
          for(int iFile = 0; iFile < nFiles; iFile++){
            uncertaintyShape[iFile][0][JffCorrector::kBackgroundFluctuation][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
            localMaximum = uncertaintyShape[iFile][0][JffCorrector::kBackgroundFluctuation][nAsymmetryBins][iCentrality][iTrackPt]->GetMaximum();
            if(localMaximum > maxValue) maxValue = localMaximum;
          }
          
          inclusiveUncertainty[0][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
          histogramMaximum = inclusiveUncertainty[0][iCentrality][iTrackPt]->GetMaximum();
          if(maxValue < histogramMaximum) maxValue = histogramMaximum;
     
          uncertaintyShape[0][0][JffCorrector::kBackgroundFluctuation][nAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0, maxValue*1.3);
          uncertaintyShape[0][0][JffCorrector::kBackgroundFluctuation][nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(kBlack);
          if(smoothHistograms) uncertaintyShape[0][0][JffCorrector::kBackgroundFluctuation][nAsymmetryBins][iCentrality][iTrackPt]->Smooth(1,"R");
          drawer->DrawHistogram(uncertaintyShape[0][0][JffCorrector::kBackgroundFluctuation][nAsymmetryBins][iCentrality][iTrackPt], "#Deltar", "P(#Deltar) error", " ", "");
          
          for(int iFile = 1; iFile < nFiles; iFile++){
            uncertaintyShape[iFile][0][JffCorrector::kBackgroundFluctuation][nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(trueColors[iFile]);
            if(smoothHistograms) uncertaintyShape[iFile][0][JffCorrector::kBackgroundFluctuation][nAsymmetryBins][iCentrality][iTrackPt]->Smooth(1,"R");
            uncertaintyShape[iFile][0][JffCorrector::kBackgroundFluctuation][nAsymmetryBins][iCentrality][iTrackPt]->Draw("same");
          }
          
          
          inclusiveUncertainty[0][iCentrality][iTrackPt]->SetLineColor(kRed);
          inclusiveUncertainty[0][iCentrality][iTrackPt]->Draw("same");
          
          // Make a legend for the titles
          legend = new TLegend(0.12,0.77,0.4,0.92);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0,"Background fluctuation","");
          legend->AddEntry((TObject*) 0,centralityString.Data(),"");
          legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
          legend->Draw();
          
          // Create a new legend for uncertainties
          legend = new TLegend(0.6,0.85-0.05*nFiles,0.9,0.92);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
          
          // Add the first uncertainty source to legend
          for(int iFile = 0; iFile < nFiles; iFile++){
            legend->AddEntry(uncertaintyShape[iFile][0][JffCorrector::kBackgroundFluctuation][nAsymmetryBins][iCentrality][iTrackPt], legendLabels[iFile], "l");
          }
          legend->AddEntry(inclusiveUncertainty[0][iCentrality][iTrackPt], "Inclusive analysis", "l");
          
          // Draw the legend to the figure
          legend->Draw();
          
          // Save the figures into a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/uncertaintyInclusiveComparisonBackgroundFluctuation%s%s.pdf", compactCentralityString.Data(), compactTrackPtString.Data()));
          }
          
          // ======================================================
          // == Compare with inclusive on jet fragmentation bias ==
          // ======================================================
          maxValue = 0;
          for(int iFile = 0; iFile < nFiles; iFile++){
            uncertaintyShape[iFile][0][JffCorrector::kFragmentationBias][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
            localMaximum = uncertaintyShape[iFile][0][JffCorrector::kFragmentationBias][nAsymmetryBins][iCentrality][iTrackPt]->GetMaximum();
            if(localMaximum > maxValue) maxValue = localMaximum;
          }
          inclusiveUncertainty[1][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
          histogramMaximum = inclusiveUncertainty[1][iCentrality][iTrackPt]->GetMaximum();
          if(maxValue < histogramMaximum) maxValue = histogramMaximum;
          
          uncertaintyShape[0][0][JffCorrector::kFragmentationBias][nAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0,maxValue*1.3);
          uncertaintyShape[0][0][JffCorrector::kFragmentationBias][nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(kBlack);
          if(smoothHistograms){
            uncertaintyShape[0][0][JffCorrector::kFragmentationBias][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0.05,1);
            uncertaintyShape[0][0][JffCorrector::kFragmentationBias][nAsymmetryBins][iCentrality][iTrackPt]->Smooth(1,"R");
            uncertaintyShape[0][0][JffCorrector::kFragmentationBias][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
          }
          drawer->DrawHistogram(uncertaintyShape[0][0][JffCorrector::kFragmentationBias][nAsymmetryBins][iCentrality][iTrackPt], "#Deltar", "P(#Deltar) error", " ", "");
          
          for(int iFile = 1; iFile < nFiles; iFile++){
            uncertaintyShape[iFile][0][JffCorrector::kFragmentationBias][nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(trueColors[iFile]);
            if(smoothHistograms){
              uncertaintyShape[iFile][0][JffCorrector::kFragmentationBias][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0.05,1);
              uncertaintyShape[iFile][0][JffCorrector::kFragmentationBias][nAsymmetryBins][iCentrality][iTrackPt]->Smooth(1,"R");
              uncertaintyShape[iFile][0][JffCorrector::kFragmentationBias][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
            }
            uncertaintyShape[iFile][0][JffCorrector::kFragmentationBias][nAsymmetryBins][iCentrality][iTrackPt]->Draw("same");
          }
          
          inclusiveUncertainty[1][iCentrality][iTrackPt]->SetLineColor(kRed);
          inclusiveUncertainty[1][iCentrality][iTrackPt]->Draw("same");
          
          // Make a legend for the titles
          legend = new TLegend(0.12,0.77,0.4,0.92);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0,"JFF bias","");
          legend->AddEntry((TObject*) 0,centralityString.Data(),"");
          legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
          legend->Draw();
          
          // Create a new legend for uncertainties
          legend = new TLegend(0.6,0.8,0.9,0.92);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
          
          // Add the first uncertainty source to legend
          for(int iFile = 0; iFile < nFiles; iFile++){
            legend->AddEntry(uncertaintyShape[iFile][0][JffCorrector::kFragmentationBias][nAsymmetryBins][iCentrality][iTrackPt], legendLabels[iFile], "l");
          }
          legend->AddEntry(inclusiveUncertainty[1][iCentrality][iTrackPt], "Inclusive analysis", "l");
          
          // Draw the legend to the figure
          legend->Draw();
          
          // Save the figures into a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/uncertaintyInclusiveComparisonJetFragmentation%s%s.pdf", compactCentralityString.Data(), compactTrackPtString.Data()));
          }
          
          // ========================================================================
          // == Compare with inclusive on background subtraction + pair acceptance ==
          // ========================================================================
          
          // Add pair acceptance and background uncertainties in quadrature
          for(int iFile = 0; iFile < nFiles; iFile++){
            for(int iBin = 1; iBin <= uncertaintyShape[iFile][0][JffCorrector::kPairAcceptance][nAsymmetryBins][iCentrality][iTrackPt]->GetNbinsX(); iBin++){
              oldContent = uncertaintyShape[iFile][0][JffCorrector::kPairAcceptance][nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(iBin);
              addedContent = uncertaintyShape[iFile][0][JffCorrector::kBackgroundSubtraction][nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(iBin);
              newContent = TMath::Sqrt(oldContent*oldContent+addedContent*addedContent);
              uncertaintyShape[iFile][0][JffCorrector::kPairAcceptance][nAsymmetryBins][iCentrality][iTrackPt]->SetBinContent(iBin,newContent);
            } // Histogram bin loop
          } // Uncertainty file loop
          
          maxValue = 0;
          for(int iFile = 0; iFile < nFiles; iFile++){
            uncertaintyShape[iFile][0][JffCorrector::kPairAcceptance][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
            localMaximum = uncertaintyShape[iFile][0][JffCorrector::kPairAcceptance][nAsymmetryBins][iCentrality][iTrackPt]->GetMaximum();
            if(localMaximum > maxValue) maxValue = localMaximum;
          }
          
          // There is a bug in the inclusive file, the range need to be manually set to the value found from the last bin
          inclusiveUncertainty[2][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
          histogramMaximum = inclusiveUncertainty[2][iCentrality][iTrackPt]->GetBinContent(inclusiveUncertainty[2][iCentrality][iTrackPt]->FindBin(0.99));
          inclusiveUncertainty[2][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0,histogramMaximum);
          if(maxValue < histogramMaximum) maxValue = histogramMaximum;
          
          uncertaintyShape[0][0][JffCorrector::kPairAcceptance][nAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0,maxValue*1.3);
          uncertaintyShape[0][0][JffCorrector::kPairAcceptance][nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(kBlack);
          drawer->DrawHistogram(uncertaintyShape[0][0][JffCorrector::kPairAcceptance][nAsymmetryBins][iCentrality][iTrackPt], "#Deltar", "P(#Deltar) error", " ", "");
          
          for(int iFile = 1; iFile < nFiles; iFile++){
            uncertaintyShape[iFile][0][JffCorrector::kPairAcceptance][nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(trueColors[iFile]);
            uncertaintyShape[iFile][0][JffCorrector::kPairAcceptance][nAsymmetryBins][iCentrality][iTrackPt]->Draw("same");
          }
          
          inclusiveUncertainty[2][iCentrality][iTrackPt]->SetLineColor(kRed);
          inclusiveUncertainty[2][iCentrality][iTrackPt]->Draw("same");
          
          // Make a legend for the titles
          legend = new TLegend(0.12,0.77,0.4,0.92);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0,"Acceptance + background","");
          legend->AddEntry((TObject*) 0,centralityString.Data(),"");
          legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
          legend->Draw();
          
          // Create a new legend for uncertainties
          legend = new TLegend(0.6,0.8,0.9,0.92);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
          
          // Add the first uncertainty source to legend
          for(int iFile = 0; iFile < nFiles; iFile++){
            legend->AddEntry(uncertaintyShape[iFile][0][JffCorrector::kPairAcceptance][nAsymmetryBins][iCentrality][iTrackPt], legendLabels[iFile], "l");
          }
          legend->AddEntry(inclusiveUncertainty[2][iCentrality][iTrackPt], "Inclusive analysis", "l");
          
          // Draw the legend to the figure
          legend->Draw();
          
          // Save the figures into a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/uncertaintyInclusiveComparisonBackgroundAcceptance%s%s.pdf", compactCentralityString.Data(), compactTrackPtString.Data()));
          }
          
          // =================================================
          // == Compare with inclusive on the rest combined ==
          // =================================================
          
          // Add pair acceptance and background uncertainties in quadrature
          for(int iFile = 0; iFile < nFiles; iFile++){
            for(int iBin = 1; iBin <= uncertaintyShape[iFile][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt]->GetNbinsX(); iBin++){
              oldContent = uncertaintyShape[iFile][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(iBin);
              addedContent = uncertaintyShape[iFile][0][JffCorrector::kTrackingEfficiency][nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(iBin);
              histogramMaximum = uncertaintyShape[iFile][0][JffCorrector::kResidualTracking][nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(iBin);
              maxValue = uncertaintyShape[iFile][0][JffCorrector::kTrackingDeltaR][nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(iBin);
              newContent = TMath::Sqrt(oldContent*oldContent+addedContent*addedContent+histogramMaximum*histogramMaximum+maxValue*maxValue);
              uncertaintyShape[iFile][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt]->SetBinContent(iBin,newContent);
            } // Histogram bin loop
          } // Uncertainty file loop
          
          maxValue = 0;
          for(int iFile = 0; iFile < nFiles; iFile++){
            uncertaintyShape[iFile][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
            localMaximum = uncertaintyShape[iFile][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt]->GetMaximum();
            if(localMaximum > maxValue) maxValue = localMaximum;
          }
          inclusiveUncertainty[3][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
          histogramMaximum = inclusiveUncertainty[3][iCentrality][iTrackPt]->GetMaximum();
          if(maxValue < histogramMaximum) maxValue = histogramMaximum;
          
          uncertaintyShape[0][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0,maxValue*1.3);
          uncertaintyShape[0][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(kBlack);
          if(smoothHistograms){
            uncertaintyShape[0][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0.05,1);
            uncertaintyShape[0][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt]->Smooth(1,"R");
            uncertaintyShape[0][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
          }
          drawer->DrawHistogram(uncertaintyShape[0][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt], "#Deltar", "P(#Deltar) error", " ", "");
          
          for(int iFile = 1; iFile < nFiles; iFile++){
            uncertaintyShape[iFile][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(trueColors[iFile]);
            if(smoothHistograms){
              uncertaintyShape[iFile][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0.05,1);
              uncertaintyShape[iFile][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt]->Smooth(1,"R");
              uncertaintyShape[iFile][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
            }
            uncertaintyShape[iFile][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt]->Draw("same");
          }
          
          inclusiveUncertainty[3][iCentrality][iTrackPt]->SetLineColor(kRed);
          inclusiveUncertainty[3][iCentrality][iTrackPt]->Draw("same");
          
          // Make a legend for the titles
          legend = new TLegend(0.12,0.77,0.4,0.92);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0,"Tracking + JES","");
          legend->AddEntry((TObject*) 0,centralityString.Data(),"");
          legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
          legend->Draw();
          
          // Create a new legend for uncertainties
          legend = new TLegend(0.6,0.8,0.9,0.92);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
          
          // Add the first uncertainty source to legend
          for(int iFile = 0; iFile < nFiles; iFile++){
            legend->AddEntry(uncertaintyShape[iFile][0][JffCorrector::kJetEnergyScale][nAsymmetryBins][iCentrality][iTrackPt], legendLabels[iFile], "l");
          }
          legend->AddEntry(inclusiveUncertainty[3][iCentrality][iTrackPt], "Inclusive analysis", "l");
          
          // Draw the legend to the figure
          legend->Draw();
          
          // Save the figures into a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/uncertaintyInclusiveComparisonTrackingJetEnergyScale%s%s.pdf", compactCentralityString.Data(), compactTrackPtString.Data()));
          }
          
          // ==============================================
          // == Compare with total inclusive uncertainty ==
          // ==============================================
          
          maxValue = 0;
          for(int iFile = 0; iFile < nFiles; iFile++){
            uncertaintyShape[iFile][0][JffCorrector::kTotal][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
            localMaximum = uncertaintyShape[iFile][0][JffCorrector::kTotal][nAsymmetryBins][iCentrality][iTrackPt]->GetMaximum();
            if(localMaximum > maxValue) maxValue = localMaximum;
          }
          inclusiveUncertainty[4][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
          histogramMaximum = inclusiveUncertainty[4][iCentrality][iTrackPt]->GetMaximum();
          if(maxValue < histogramMaximum) maxValue = histogramMaximum;
          
          uncertaintyShape[0][0][JffCorrector::kTotal][nAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0,maxValue*1.3);
          uncertaintyShape[0][0][JffCorrector::kTotal][nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(kBlack);
          if(smoothHistograms){
            uncertaintyShape[0][0][JffCorrector::kTotal][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0.05,1);
            uncertaintyShape[0][0][JffCorrector::kTotal][nAsymmetryBins][iCentrality][iTrackPt]->Smooth(1,"R");
            uncertaintyShape[0][0][JffCorrector::kTotal][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
          }
          drawer->DrawHistogram(uncertaintyShape[0][0][JffCorrector::kTotal][nAsymmetryBins][iCentrality][iTrackPt], "#Deltar", "P(#Deltar) error", " ", "");
          
          for(int iFile = 1; iFile < nFiles; iFile++){
            uncertaintyShape[iFile][0][JffCorrector::kTotal][nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(trueColors[iFile]);
            if(smoothHistograms){
              uncertaintyShape[iFile][0][JffCorrector::kTotal][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0.05,1);
              uncertaintyShape[iFile][0][JffCorrector::kTotal][nAsymmetryBins][iCentrality][iTrackPt]->Smooth(1,"R");
              uncertaintyShape[iFile][0][JffCorrector::kTotal][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
            }
            uncertaintyShape[iFile][0][JffCorrector::kTotal][nAsymmetryBins][iCentrality][iTrackPt]->Draw("same");
          }
          
          inclusiveUncertainty[4][iCentrality][iTrackPt]->SetLineColor(kRed);
          inclusiveUncertainty[4][iCentrality][iTrackPt]->Draw("same");
          
          // Make a legend for the titles
          legend = new TLegend(0.12,0.77,0.4,0.92);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0,"Total uncertainty","");
          legend->AddEntry((TObject*) 0,centralityString.Data(),"");
          legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
          legend->Draw();
          
          // Create a new legend for uncertainties
          legend = new TLegend(0.6,0.8,0.9,0.92);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
          
          // Add the first uncertainty source to legend
          for(int iFile = 0; iFile < nFiles; iFile++){
            legend->AddEntry(uncertaintyShape[iFile][0][JffCorrector::kTotal][nAsymmetryBins][iCentrality][iTrackPt], legendLabels[iFile], "l");
          }
          legend->AddEntry(inclusiveUncertainty[4][iCentrality][iTrackPt], "Inclusive analysis", "l");
          
          // Draw the legend to the figure
          legend->Draw();
          
          // Save the figures into a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/uncertaintyInclusiveComparisonTotal%s%s.pdf", compactCentralityString.Data(), compactTrackPtString.Data()));
          }
          
        } // Track pT loop
      } // Centrality loop
      
    } // Comparison to inclusive
    
  } // Drawing the figures if
}

