#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "JDrawer.h"

/*
 * Macro for producing the JFF correction from PYTHIA simulation results
 */
void compareLongRangeAsymmetry(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString fileName = "data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_adjustedBackground_processed_2019-07-05.root";
  
  bool saveFigures = true;
  
  // Open the input files
  TFile *dataFile = TFile::Open(fileName);
  
  // Create a reader for the histograms and define binning
  DijetHistogramManager *reader = new DijetHistogramManager(dataFile);
  
  const int nCentralityBins = reader->GetNCentralityBins();
  const int nTrackPtBins = reader->GetNTrackPtBins();
  const int nAsymmetryBins = reader->GetNAsymmetryBins();
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  // Define arrays for the histograms
  TH1D *background[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TF1 *backgroundFit[nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
  // Read the histograms from the input file
  reader->SetLoadTrackLeadingJetCorrelations(true);
  reader->SetAsymmetryBinRange(0,nAsymmetryBins);
  reader->LoadProcessedHistograms();
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        background[iAsymmetry][iCentrality][iTrackPt] = reader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta);
        backgroundFit[iAsymmetry][iCentrality][iTrackPt] = background[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
      }
    }
  }
  
  // Draw all asymmery bins deltaR curves to a single plot to see if the spillover is similar regardless of asymmetry bin
  JDrawer *drawer = new JDrawer();
  
  // Make subeNon0 to spillover comparison in all bins
  TLegend *legend;
  int colors[] = {kBlue,kRed,kGreen+4,kBlack};
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      
      legend = new TLegend(0.35,0.75,0.75,0.99);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->SetHeader(Form("%.1f < p_{T} < %.1f, C: %.0f-%.0f %%",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1], centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
      
      
      /*background[nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(colors[nAsymmetryBins]);
      background[nAsymmetryBins][iCentrality][iTrackPt]->Rebin(4);
      background[nAsymmetryBins][iCentrality][iTrackPt]->Scale(1/4.0);
      backgroundFit[nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(colors[nAsymmetryBins]);
      //backgroundFit[nAsymmetryBins][iCentrality][iTrackPt]->SetLineWidth(0);
      legend->AddEntry(background[nAsymmetryBins][iCentrality][iTrackPt],"All x_{j}","l");*/
      
      drawer->DrawHistogram(background[0][iCentrality][iTrackPt],"#Delta#phi","#frac{dN}{d#Delta#phi}"," ");
      
      for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
        background[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(colors[iAsymmetry]);
        background[iAsymmetry][iCentrality][iTrackPt]->Rebin(4);
        background[iAsymmetry][iCentrality][iTrackPt]->Scale(1/4.0);
        background[iAsymmetry][iCentrality][iTrackPt]->Draw("same");
        backgroundFit[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(colors[iAsymmetry]);
        //backgroundFit[iAsymmetry][iCentrality][iTrackPt]->SetLineWidth(0);
        legend->AddEntry(background[iAsymmetry][iCentrality][iTrackPt],Form("%.1f < x_{j} < %.1f",xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]),"l");
      }
      
      legend->Draw();
      
      
      if(saveFigures){
       gPad->GetCanvas()->SaveAs(Form("figures/longRangeAsymmetryComparison_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
      }
    } // Track pt Loop
  } // Centrality loop
  
}

