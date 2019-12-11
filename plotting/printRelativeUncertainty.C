#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "JDrawer.h"

/*
 * Find the maximum value over centrality bins, when other bins are kept constants
 *
 *  TH1D* histogramArray[2][JffCorrector::knUncertaintySources][5][8] = Array of all uncertainty histograms
 *  int iJetType = Constant value for the jet type bin
 *  int iUncertainty = Constant value for the uncertainty source bin
 *  int iTrackPt = Constant value for the track pT bin
 */
double findMaximumValue(TH1D* histogramArray[2][JffCorrector::knUncertaintySources][5][8], int iJetType, int iUncertainty, int iTrackPt){
  
  double maxValue = 0;
  double histogramMaximum = 0;
  
  for(int iCentrality = 0; iCentrality < 5; iCentrality++){
    histogramArray[iJetType][iUncertainty][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
    histogramMaximum = histogramArray[iJetType][iUncertainty][iCentrality][iTrackPt]->GetMaximum();
    if(histogramMaximum > maxValue) maxValue = histogramMaximum;
  }
  
  return maxValue;
}

/*
 * Macro for printing the relative uncertainties in latex slides
 */
void printRelativeUncertainty(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString dataFileName = "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsWithShifterCentrality_trackDeltaRonlyInLowPt_processed_2019-10-17.root"; // Compare also with uncorrected data
  // data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root
  TString ppFileName = "data/ppData2017_highForest_pfJets_20EventsMixed_finalTrackCorr_xjBins_JECv4_wtaAxis_tunedSeagull_allCorrections_processed_2019-10-17.root";
  TString uncertaintyFileName = "uncertainties/systematicUncertaintyForPbPb_25eveMix_oldJES_15percentSpill10Jff_2019-10-17.root";
  TString ppUncertainryFileName = "uncertainties/systematicUncertaintyForPp_20percentSpillJff_2019-09-30.root";
  
  bool printSlides = false;  // Print slides showing the R-integrated uncertainty in each pT bin
  bool drawUncertaintySourceComparison = false; // Draw all uncertainty sources as a function of R in each pT bin
  bool drawUncertaintySystemComparison = true; // Draw single uncertainty source for all systems as a function of R in each pT bin
  bool saveFigures = true;   // Save the drawn figures to file
  
  // ==================================================================
  // ====================== Configuration done ========================
  // ==================================================================
  
  // Open the input files
  TFile *dataFile = TFile::Open(dataFileName);
  TFile *uncertaintyFile = TFile::Open(uncertaintyFileName);
  TFile *ppFile = TFile::Open(ppFileName);
  TFile *ppUncertainryFile = TFile::Open(ppUncertainryFileName);
  
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
  
  JffCorrector *uncertaintyManager = new JffCorrector();
  uncertaintyManager->ReadSystematicFile(uncertaintyFile);
  
  JffCorrector *ppUncertaintyManager = new JffCorrector();
  ppUncertaintyManager->ReadSystematicFile(ppUncertainryFile);
  
  const int nCentralityBins = dataManager->GetNCentralityBins();
  const int nTrackPtBins = dataManager->GetNTrackPtBins();
  const int nAsymmetryBins = dataManager->GetNAsymmetryBins();
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  // Select drawn bins
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins-1;
  
  int firstDrawnAsymmetryBin = 0;
  int lastDrawnAsymmetryBin = nAsymmetryBins;
  
  // If we are drawing system comparison or printing slides, we must use all centrality bins
  if(drawUncertaintySystemComparison || printSlides){
    firstDrawnCentralityBin = 0;
    lastDrawnCentralityBin = nCentralityBins;
  }
  
  // Define arrays for the jet shapes
  TH1D *dataShape[2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins+1];
  TH1D *uncertaintyShape[2][JffCorrector::knUncertaintySources][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins+1];
  
  // Define arrays for yield
  double dataYield[2][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins+1];
  double uncertaintyYield[2][JffCorrector::knUncertaintySources][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins+1];
  
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
        cout << "iCentrality: " << iCentrality << " iTrackPt: " << iTrackPt << " iAsymmetry: " << iAsymmetry << endl;
        dataYield[0][iAsymmetry][iCentrality][iTrackPt] = dataShape[0][iAsymmetry][iCentrality][iTrackPt]->Integral(1, dataShape[0][iAsymmetry][iCentrality][iTrackPt]->FindBin(0.99), "width");
        dataYield[1][iAsymmetry][iCentrality][iTrackPt] = dataShape[1][iAsymmetry][iCentrality][iTrackPt]->Integral(1, dataShape[1][iAsymmetry][iCentrality][iTrackPt]->FindBin(0.99), "width");
        
        for(int iUncertainty = 0; iUncertainty < JffCorrector::knUncertaintySources; iUncertainty++){
          
          // Different reader for pp and PbPb uncertainty files
          if(iCentrality == nCentralityBins){
            uncertaintyShape[0][iUncertainty][iAsymmetry][iCentrality][iTrackPt] = ppUncertaintyManager->GetJetShapeSystematicUncertainty(DijetHistogramManager::kPtWeightedTrackLeadingJet, 0, iTrackPt, iAsymmetry, iUncertainty);
            uncertaintyShape[1][iUncertainty][iAsymmetry][iCentrality][iTrackPt] = ppUncertaintyManager->GetJetShapeSystematicUncertainty(DijetHistogramManager::kPtWeightedTrackSubleadingJet, 0, iTrackPt, iAsymmetry, iUncertainty);
          } else {
            uncertaintyShape[0][iUncertainty][iAsymmetry][iCentrality][iTrackPt] = uncertaintyManager->GetJetShapeSystematicUncertainty(DijetHistogramManager::kPtWeightedTrackLeadingJet, iCentrality, iTrackPt, iAsymmetry, iUncertainty);
            uncertaintyShape[1][iUncertainty][iAsymmetry][iCentrality][iTrackPt] = uncertaintyManager->GetJetShapeSystematicUncertainty(DijetHistogramManager::kPtWeightedTrackSubleadingJet, iCentrality, iTrackPt, iAsymmetry, iUncertainty);
          }
          
          uncertaintyYield[0][iUncertainty][iAsymmetry][iCentrality][iTrackPt] = uncertaintyShape[0][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->Integral(1, uncertaintyShape[0][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->FindBin(0.99), "width");
          uncertaintyYield[1][iUncertainty][iAsymmetry][iCentrality][iTrackPt] = uncertaintyShape[1][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->Integral(1, uncertaintyShape[1][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->FindBin(0.99), "width");
        } // Uncertainty type loop
      } // Asymmetry loop
    } // Track pT loop
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
        cout << "  \\begin{tabular}{cccccc}" << endl;
        cout << "    \\toprule" << endl;
        cout << "    Source & C: 0-10 & C: 10-30 & C: 30-50 & C: 50-90 & pp \\\\" << endl;
        cout << "    \\midrule" << endl;
        
        // Set the correct precision for printing floating point numbers
        cout << fixed << setprecision(3);
        
        // Print one line for each track pT bin
        for(int iUncertainty = 0; iUncertainty < JffCorrector::knUncertaintySources; iUncertainty++){
          cout << uncertaintyManager->GetUncertaintyName(iUncertainty);
          for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
            cout << " & $" << uncertaintyYield[iJetType][iUncertainty][nAsymmetryBins][iCentrality][iTrackPt] / dataYield[iJetType][nAsymmetryBins][iCentrality][iTrackPt] << "$";
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
  
  if(drawUncertaintySourceComparison || drawUncertaintySystemComparison){
    
    JDrawer *drawer = new JDrawer();
    TLegend *legend;
    int colors[] = {kBlue, kRed, kMagenta, kCyan, kGreen+3, kViolet+7, kOrange+7, kSpring, kGray};
    double maxValue = 0;
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
              uncertaintyShape[iJetType][JffCorrector::kTotal][iAsymmetry][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
              maxValue = uncertaintyShape[iJetType][JffCorrector::kTotal][iAsymmetry][iCentrality][iTrackPt]->GetMaximum();
              uncertaintyShape[iJetType][JffCorrector::kTotal][iAsymmetry][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0,maxValue*1.3);
              uncertaintyShape[iJetType][JffCorrector::kTotal][iAsymmetry][iCentrality][iTrackPt]->SetLineColor(kBlack);
              drawer->DrawHistogram(uncertaintyShape[iJetType][JffCorrector::kTotal][iAsymmetry][iCentrality][iTrackPt], "#Deltar", "P(#Deltar) error", " ", "");
              
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
              legend->AddEntry(uncertaintyShape[iJetType][JffCorrector::kTotal][iAsymmetry][iCentrality][iTrackPt], uncertaintyManager->GetUncertaintyName(JffCorrector::kTotal), "l");
              
              // Loop over all the other uncertainty sources
              for(int iUncertainty = 0; iUncertainty < JffCorrector::kTotal; iUncertainty++){
                uncertaintyShape[iJetType][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->SetLineColor(colors[iUncertainty]);
                uncertaintyShape[iJetType][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->Draw("same");
                legend->AddEntry(uncertaintyShape[iJetType][iUncertainty][iAsymmetry][iCentrality][iTrackPt], uncertaintyManager->GetUncertaintyName(iUncertainty), "l");
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
                uncertaintyShape[iJetType][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
                histogramMaximum = uncertaintyShape[iJetType][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->GetMaximum();
                if(histogramMaximum > maxValue) maxValue = histogramMaximum;
              }
              
              // Draw the most central bin
              uncertaintyShape[iJetType][iUncertainty][iAsymmetry][0][iTrackPt]->GetYaxis()->SetRangeUser(0,maxValue*1.3);
              uncertaintyShape[iJetType][iUncertainty][iAsymmetry][0][iTrackPt]->SetLineColor(kBlack);
              drawer->DrawHistogram(uncertaintyShape[iJetType][iUncertainty][iAsymmetry][0][iTrackPt], "#Deltar", "P(#Deltar) error", " ", "");
              
              // Make a legend for the titles
              legend = new TLegend(0.12,0.77,0.4,0.92);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
              legend->AddEntry((TObject*) 0,(jetShapeString+asymmetryString).Data(),"");
              legend->AddEntry((TObject*) 0,uncertaintyManager->GetUncertaintyName(iUncertainty),"");
              legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
              legend->Draw();
              
              // Create a new legend systems
              legend = new TLegend(0.7,0.65,0.9,0.92);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
              
              // Add the most central bin to the legend
              legend->AddEntry(uncertaintyShape[iJetType][iUncertainty][iAsymmetry][0][iTrackPt], "C = 0-10 %", "l");
              
              // Loop over all the other centrality bins and pp
              for(int iCentrality = 1; iCentrality <= nCentralityBins; iCentrality++){
                
                if(iCentrality == nCentralityBins){
                  centralityString = "pp";
                } else {
                  centralityString = Form("C = %.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
                }
                
                uncertaintyShape[iJetType][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->SetLineColor(colors[iCentrality-1]);
                uncertaintyShape[iJetType][iUncertainty][iAsymmetry][iCentrality][iTrackPt]->Draw("same");
                legend->AddEntry(uncertaintyShape[iJetType][iUncertainty][iAsymmetry][iCentrality][iTrackPt], centralityString.Data(), "l");
              } // Centrality bin loop
              
              // Draw the legend to the figure
              legend->Draw();
              
              // Save the figures into a file
              if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/uncertaintySystemBreakdown%sJet_%s%s%s.pdf", jetType[iJetType], uncertaintyManager->GetUncertaintyName(iUncertainty).Data(), compactAsymmetryString.Data(), compactTrackPtString.Data()));
              }
              
            } // Track pT loop
          } // Asymmetry loop
        } // Uncertainty source loop
      } // Jet type loop (leading/subleading)
    } // Uncertainty system comparison
    
  } // Drawing the figures if
}

