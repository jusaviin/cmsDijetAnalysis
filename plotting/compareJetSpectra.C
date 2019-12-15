#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

/*
 * Macro for comparing jets shapes in different pT bins between pp and PbPb
 */
void compareJetSpectra(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString fileName = "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsWithCentShift_trackDeltaRonlyLowPt_processed_2019-10-16.root";
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsWithCentShift_trackDeltaRonlyLowPt_processed_2019-10-16.root
  
  // Open the comparison file
  TFile *file = TFile::Open(fileName);
  
  // Create histogram managers to read the histograms from the file
  DijetHistogramManager *histograms = new DijetHistogramManager(file);
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = true;
  
  // Get the number of asymmetry bins
  const int nAsymmetryBins = histograms->GetNAsymmetryBins();
  const int nTrackPtBins = histograms->GetNTrackPtBins();
  const int nCentralityBins = histograms->GetNCentralityBins();
  
  double centralityBinBorders[] = {0,10,30,50,90};
  double asymmetryBinBorders[] = {0,0.6,0.8,1};
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};
  
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int iJetTrack = DijetHistogramManager::kTrackLeadingJet; // DijetHistogramManager::kTrackSubleadingJet
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Load leading and subleading histograms from the file
  histograms->SetLoadLeadingJetHistograms(true);
  histograms->SetLoadSubleadingJetHistograms(true);
  histograms->LoadProcessedHistograms();

  // Histograms for jet spectra
  TH1D *jetPt[nAsymmetryBins+1][nCentralityBins];
  TH1D *jetPtRatio[nAsymmetryBins][nCentralityBins];
  
  // Read the jet shape histograms from the file and normalize them to one
  for(int iAsymmetry = nAsymmetryBins; iAsymmetry >= 0; iAsymmetry--){
    
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      
      jetPt[iAsymmetry][iCentrality] = histograms->GetHistogramJetPt(iJetTrack, iCentrality, iAsymmetry);
      jetPt[iAsymmetry][iCentrality]->Scale(1.0/jetPt[iAsymmetry][iCentrality]->Integral());
      
      if(iAsymmetry < nAsymmetryBins){
        
        jetPtRatio[iAsymmetry][iCentrality] = (TH1D*) jetPt[iAsymmetry][iCentrality]->Clone(Form("jetRatio%d%d", iCentrality, iAsymmetry));
        jetPtRatio[iAsymmetry][iCentrality]->Divide(jetPt[nAsymmetryBins][iCentrality]);
      }
    } // Centrality loop
  } // Track pT loop
  
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  
  TLegend *legend;
  
  TString centralityString;
  TString compactCentralityString;
  TString asymmetryString[] = {"0.0 < x_{j} < 0.6", "0.6 < x_{j} < 0.8", "0.8 < x_{j} < 1.0", "0.0 < x_{j} < 1.0"};
  
  TString figureName;
  int veryNiceColors[] = {kRed,kBlue,kGreen+3,kBlack};
  const char *labels[] = {"Leading","Leading","Leading","Subleading","Subleading","Subleading","Inclusive","Inclusive"};
  
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    
    centralityString = Form("%.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    
    drawer->CreateSplitCanvas();
    drawer->SetLogY(true);

    jetPt[nAsymmetryBins][iCentrality]->SetLineColor(kBlack);
    jetPt[nAsymmetryBins][iCentrality]->GetXaxis()->SetRangeUser(120,200);
    drawer->DrawHistogramToUpperPad(jetPt[nAsymmetryBins][iCentrality],"p_{T} (GeV)","#frac{dN}{dp_{T}}"," ");
    
    legend = new TLegend(0.25,0.05,0.65,0.4);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*) 0,Form("%s jet",labels[iJetTrack]),"");
    legend->AddEntry((TObject*) 0,centralityString.Data(),"");
    legend->AddEntry(jetPt[nAsymmetryBins][iCentrality],asymmetryString[nAsymmetryBins],"l");
    
    for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
      jetPt[iAsymmetry][iCentrality]->SetLineColor(veryNiceColors[iAsymmetry]);
      jetPt[iAsymmetry][iCentrality]->Draw("same");
      legend->AddEntry(jetPt[iAsymmetry][iCentrality],asymmetryString[iAsymmetry],"l");
    }
    
    legend->Draw();
    
    // Draw the ratio to the lower pad
    drawer->SetLogY(false);
    
    jetPtRatio[0][iCentrality]->GetYaxis()->SetRangeUser(0.7,1.3);
    jetPtRatio[0][iCentrality]->GetXaxis()->SetRangeUser(120,200);
    jetPtRatio[0][iCentrality]->SetLineColor(veryNiceColors[0]);
    drawer->DrawHistogramToLowerPad(jetPtRatio[0][iCentrality],"p_{T} (GeV)","Color / All", " ");
    
    for(int iAsymmetry = 1; iAsymmetry < nAsymmetryBins; iAsymmetry++){
      jetPtRatio[iAsymmetry][iCentrality]->SetLineColor(veryNiceColors[iAsymmetry]);
      jetPtRatio[iAsymmetry][iCentrality]->Draw("same");
    }
    
    // Save the figures into a file
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/jetSpectraComparison%s%s.pdf", labels[iJetTrack], compactCentralityString.Data()));
    }
    
  } // Centrality loop

}
