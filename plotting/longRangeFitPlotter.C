#include "JDrawer.h"

/*
 * Macro for studying the different components of the long range fits
 */
void longRangeFitPlotter(){

  // File from which the long range distributions are plotted
  TString directoryName = "flowGraphs/";
  TString fileName = "flowGraphs_PbPbMC2018_subeNon0_4pCentShift_caloJets_noQcut_correctedJetHadron_correctedDihadron_2021-03-04.root";
  // "flowGraphs_PbPb2018_caloJets_fixedJEC_correctedJetHadron_correctedDihadron_2021-02-26.root" "flowGraphs_PbPbMC2018_subeNon0_4pCentShift_caloJets_noQcut_correctedJetHadron_correctedDihadron_2021-03-04.root"
  
  TString systemAndEnergy = "Pythia+Hydjet 5.02 TeV";
  
  // Open the file
  TFile *longRangeFile = TFile::Open(directoryName+fileName);
  
  // Binning information
  const int nCentralityBins = 3;
  const int nTrackPtBins = 4;
  const int nAsymmetryBins = 3;
  const int nFlowComponents = 4;
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  // Configuration
  const bool drawJetHadron = true;  // Draw the jet hadron long range deltaPhi distribution
  const bool drawDihadron = false;   // Draw the dihadron long range deltaPhi distribution

  const bool drawOverallFit = false;  // Draw the overall Fourier fit to the distribution
  const bool drawFitDecomposition = false;  // Draw different fit components seperately
  
  const bool saveFigures = false;
  TString saveComment = "_mcWithFullFit";
  
  TString typeString[2] = {"Jet-hadron", "Dihadron"};
  TString saveTypeString[2] = {"JetHadron", "Dihadron"};
  
  // Fill drawing configuration arrays based on user input
  bool drawType[2];
  drawType[0] = drawJetHadron;
  drawType[1] = drawDihadron;
  
  // Long range distributions for jet-hadron and dihadron
  TH1D *longRangeDeltaPhi[2][nCentralityBins][nTrackPtBins]; // First bin: 0 = jet-hadron, 1 = dihadron
  
  // Read the long range distributions from the file
  char histogramNamer[150];
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      
      sprintf(histogramNamer,"longRangeJetHadron/longRangeJetHadron_A3C%dT%d", iCentrality, iTrackPt);
      longRangeDeltaPhi[0][iCentrality][iTrackPt] = (TH1D*) longRangeFile->Get(histogramNamer);
      
      sprintf(histogramNamer,"longRangeDihadron/longRangeDihadron_A3C%dT%d", iCentrality, iTrackPt);
      longRangeDeltaPhi[1][iCentrality][iTrackPt] = (TH1D*) longRangeFile->Get(histogramNamer);
      
    } // track pT loop
  } // centrality loop
  
  // Configure the histogram drawing class
  JDrawer *drawer = new JDrawer();
  
  double maxYscale, minYscale, yDifference;
  TLegend *legend;
    
  // Draw the distributions
  for(int iType = 0; iType < 2; iType++){
    if(!drawType[iType]) continue;  // Only draw jet-hadron and dihadron if selected
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){

        // If you do not want to draw the overall fit, remove it from long range histogram
        if(!drawOverallFit){
          TF1 *fourierFit = longRangeDeltaPhi[iType][iCentrality][iTrackPt]->GetFunction("fourier");
          fourierFit->SetLineWidth(0);
        }
        
        // Set the y-axis scaling so that there is some room for legend
        maxYscale = longRangeDeltaPhi[iType][iCentrality][iTrackPt]->GetMaximum();
        minYscale = longRangeDeltaPhi[iType][iCentrality][iTrackPt]->GetMinimum();
        yDifference = maxYscale - minYscale;
        maxYscale = maxYscale + 0.5 * yDifference;
        minYscale = minYscale - 0.12 * yDifference;
        longRangeDeltaPhi[iType][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(minYscale, maxYscale);
        
        // Draw the histogram to the canvas
        drawer->DrawHistogram(longRangeDeltaPhi[iType][iCentrality][iTrackPt], Form("%s #Delta#varphi", typeString[iType].Data()), "#frac{1}{N_{jet}} #frac{dN}{d#Delta#varphi}", " ");
        
        // Draw a legend for the histogram
        legend = new TLegend(0.17,0.7,0.37,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        legend->AddEntry((TObject*) 0, Form("C = %.0f-%.0f%%", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "");
        legend->AddEntry((TObject*) 0, Form("%.1f < p_{T} < %.1f GeV", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]), "");
        
        if(drawOverallFit){
          TF1 *fourierFit = longRangeDeltaPhi[iType][iCentrality][iTrackPt]->GetFunction("fourier");
          legend->AddEntry(fourierFit, "Fourier fit", "l");
        } // Drawing overall fit
        
        legend->Draw();
        
        // If we want to draw a decomposition of the fourier fit, do it
        if(drawFitDecomposition){
          TF1 *fourierFit = longRangeDeltaPhi[iType][iCentrality][iTrackPt]->GetFunction("fourier");
          
          TLegend *vLegend = new TLegend(0.56,0.7,0.76,0.9);
          vLegend->SetFillStyle(0);vLegend->SetBorderSize(0);vLegend->SetTextSize(0.045);vLegend->SetTextFont(62);
          
          TF1 *fitComponents[nFlowComponents];
          int componentColors[] = {kGreen+4, kBlue, kMagenta, kCyan, kViolet};
          for(int iFlow = 0; iFlow < nFlowComponents; iFlow++){
            fitComponents[iFlow] = new TF1(Form("fv%d", iFlow+1), Form("[0]+[0]*[1]*2.0*TMath::Cos(%d.0*x)", iFlow+1), -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
            fitComponents[iFlow]->SetParameter(0,fourierFit->GetParameter(0));
            fitComponents[iFlow]->SetParameter(1,fourierFit->GetParameter(iFlow+1));
            fitComponents[iFlow]->SetLineColor(componentColors[iFlow]);
            fitComponents[iFlow]->Draw("same");
            vLegend->AddEntry(fitComponents[iFlow],Form("v_{%d}", iFlow+1),"l");
          }
          
          vLegend->Draw();
          
        } // Drawing fit composition
        
        // Save the figures to file
        if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/longRangeDistribution%s%s_C=%.0f-%.0f_T=%.0f-%.0f.pdf", saveTypeString[iType].Data(), saveComment.Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
        } // Saving figures
        
      } // track pT loop
    } // centrality loop
  }
  
  
}
