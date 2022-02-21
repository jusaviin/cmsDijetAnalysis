#include "JDrawer.h"

/*
 * Macro for studying the different components of the long range fits
 */
void longRangeFitComparer(){

  // File from which the long range distributions are plotted
  TString directoryName = "flowGraphs/";
  const int nFiles = 2;
  TString fileName[] = {"flowGraphs_PbPb2018_caloJets_jetEta1v3_correctedJetHadron_correctedDihadron_2022-02-11.root", "flowGraphs_PbPb2018_akPfCsJets_correctedJetHadron_correctedDihadronSmallStats_2021-05-13.root", "flowGraphs_PbPb2018MC_caloJets_4pCentShift_subeNon0_multiplicityWeight_2022-01-26.root",  "flowGraphs_PbPbMC2018_pfCsJets_multWeight_subeNon0_jetEta1v6_correctedJetHadron_correctedDihadron_2022-02-12.root",   "flowGraphs_PbPbMC2018_subeNon0_4pCentShift_pfCsJets_noQcut_correctedJetHadron_correctedDihadronFromCalo_2021-06-04.root", "flowGraphs_PbPb2018MC_caloJets_4pCentShift_subeNon0_jetEta1v6_consistencyCheck_2022-02-08.root"};
  // "flowGraphs_PbPb2018_caloJets_fixedJEC_correctedJetHadron_correctedDihadron_2021-02-26.root"
  // "flowGraphs_PbPbMC2018_subeNon0_4pCentShift_caloJets_noQcut_correctedJetHadron_correctedDihadron_2021-03-04.root"
  //  flowGraphs_PbPb2018_caloJets_jetEta1v3_correctedJetHadron_correctedDihadron_2022-02-11.root
  
  TString legendComment[] = {"Calo dijet", "PFCS dijet"};
  
  TString systemAndEnergy = "PbPb 5.02 TeV";
  
  // Open the file
  TFile *longRangeFile[nFiles];
  for(int iFile = 0; iFile < nFiles; iFile++){
    longRangeFile[iFile] = TFile::Open(directoryName+fileName[iFile]);
  }
  
  // Binning information
  const int nCentralityBins = 3;
  const int nTrackPtBins = 4;
  const int nAsymmetryBins = 3;
  const int nFlowComponents = 4;
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  bool scaleToOne = true;  // To see more clearly the vn variation, we can scale the average of the histogram to one
  
  // Configuration
  const bool drawJetHadron = true;  // Draw the jet hadron long range deltaPhi distribution
  const bool drawDihadron = false;   // Draw the dihadron long range deltaPhi distribution

  const bool drawOverallFit = false;  // Draw the overall Fourier fit to the distribution
  const bool drawV1 = false;         // Draw different fit components seperately
  const bool drawV2 = true;
  const bool drawV3 = false;
  const bool drawV4 = false;
  const bool drawFitDecomposition = drawV1 || drawV2 || drawV3 || drawV4;
  bool drawVn[4] = {drawV1, drawV2, drawV3, drawV4};
  
  const bool drawRatio = true;
  
  const bool saveFigures = true;
  TString saveComment = "_jetCollectionData";
  TString figureFormat = "png";
  
  TString typeString[2] = {"Jet-hadron", "Dihadron"};
  TString saveTypeString[2] = {"JetHadron", "Dihadron"};
  
  // Fill drawing configuration arrays based on user input
  bool drawType[2];
  drawType[0] = drawJetHadron;
  drawType[1] = drawDihadron;
  
  // Coloring scheme
  int lineColors[] = {kBlue, kRed, kGreen+3, kCyan, kMagenta};
  
  // Long range distributions for jet-hadron and dihadron
  TH1D *longRangeDeltaPhi[nFiles][2][nCentralityBins][nTrackPtBins]; // First bin: 0 = jet-hadron, 1 = dihadron
  TH1D *longRangeDeltaPhiRatio[nFiles-1][2][nCentralityBins][nTrackPtBins];
  TF1 *fitRatio[nFiles-1][2][nCentralityBins][nTrackPtBins];
  
  // Read the long range distributions from the file
  char histogramNamer[150];
  double averageValue;
  TH1D *yieldFinder;
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        sprintf(histogramNamer,"longRangeJetHadron/longRangeJetHadron_A3C%dT%d", iCentrality, iTrackPt);
        longRangeDeltaPhi[iFile][0][iCentrality][iTrackPt] = (TH1D*) longRangeFile[iFile]->Get(histogramNamer);
        longRangeDeltaPhi[iFile][0][iCentrality][iTrackPt]->SetLineColor(lineColors[iFile]);
        
        sprintf(histogramNamer,"longRangeDihadron/longRangeDihadron_A3C%dT%d", iCentrality, iTrackPt);
        longRangeDeltaPhi[iFile][1][iCentrality][iTrackPt] = (TH1D*) longRangeFile[iFile]->Get(histogramNamer);
        longRangeDeltaPhi[iFile][1][iCentrality][iTrackPt]->SetLineColor(lineColors[iFile]);
        
        // Scale the average yields to one
        if(scaleToOne){
          for(int iType = 0; iType < 2; iType++){
            yieldFinder = (TH1D*)longRangeDeltaPhi[iFile][iType][iCentrality][iTrackPt]->Clone(Form("yieldFinder%d%d%d%d", iFile, iCentrality, iTrackPt, iType));
            yieldFinder->Fit("pol0");
            averageValue = yieldFinder->GetFunction("pol0")->GetParameter(0);
            longRangeDeltaPhi[iFile][iType][iCentrality][iTrackPt]->Scale(1/averageValue);
            
            TF1 *fourierFit = longRangeDeltaPhi[iFile][iType][iCentrality][iTrackPt]->GetFunction("fourier");
            fourierFit->SetParameter(0,1);
            
          }
        }
        
        // Calculate the ratios of the histograms
        if(iFile > 0){
          for(int iType = 0; iType < 2; iType++){
            
            longRangeDeltaPhiRatio[iFile-1][iType][iCentrality][iTrackPt] = (TH1D*) longRangeDeltaPhi[iFile][iType][iCentrality][iTrackPt]->Clone(Form("%sRatio%d%d%d", saveTypeString[iType].Data(), iFile, iCentrality, iTrackPt));
            longRangeDeltaPhiRatio[iFile-1][iType][iCentrality][iTrackPt]->Divide(longRangeDeltaPhi[0][iType][iCentrality][iTrackPt]);
            TF1 *fourierFit = longRangeDeltaPhiRatio[iFile-1][iType][iCentrality][iTrackPt]->GetFunction("fourier");
            fourierFit->SetLineWidth(0);
            
          } // Type loop
        } // Determining ratios
        
      } // track pT loop
    } // centrality loop
  } // File loop
  
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
          for(int iFile = 0; iFile < nFiles; iFile++){
            TF1 *fourierFit = longRangeDeltaPhi[iFile][iType][iCentrality][iTrackPt]->GetFunction("fourier");
            fourierFit->SetLineWidth(0);
          }
        } else {
          for(int iFile = 0; iFile < nFiles; iFile++){
            TF1 *fourierFit = longRangeDeltaPhi[iFile][iType][iCentrality][iTrackPt]->GetFunction("fourier");
            fourierFit->SetLineColor(lineColors[iFile]);
          }
        }
        
        // Set the y-axis scaling so that there is some room for legend
        maxYscale = longRangeDeltaPhi[1][iType][iCentrality][iTrackPt]->GetMaximum();
        minYscale = longRangeDeltaPhi[1][iType][iCentrality][iTrackPt]->GetMinimum();
        yDifference = maxYscale - minYscale;
        maxYscale = maxYscale + 0.5 * yDifference;
        minYscale = minYscale - 0.12 * yDifference;
        longRangeDeltaPhi[0][iType][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(minYscale, maxYscale);
        
        // Draw the histogram to the canvas
        drawer->DrawHistogram(longRangeDeltaPhi[0][iType][iCentrality][iTrackPt], Form("%s #Delta#varphi", typeString[iType].Data()), "#frac{1}{N_{jet}} #frac{dN}{d#Delta#varphi}", " ");
        
        for(int iFile = 1; iFile < nFiles; iFile++){
          longRangeDeltaPhi[iFile][iType][iCentrality][iTrackPt]->Draw("same");
        }
        
        // Draw a legend for the histogram
        legend = new TLegend(0.17,0.7,0.37,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        
        legend->AddEntry((TObject*) 0, systemAndEnergy.Data(), "");
        legend->AddEntry((TObject*) 0, Form("C = %.0f-%.0f%%", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "");
        legend->AddEntry((TObject*) 0, Form("%.1f < p_{T} < %.1f GeV", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]), "");
        legend->Draw();
        
        TLegend *vLegend = new TLegend(0.56,0.7,0.76,0.9);
        vLegend->SetFillStyle(0);vLegend->SetBorderSize(0);vLegend->SetTextSize(0.045);vLegend->SetTextFont(62);
        for(int iFile = 0; iFile < nFiles; iFile++){
          vLegend->AddEntry(longRangeDeltaPhi[iFile][iType][iCentrality][iTrackPt],legendComment[iFile],"l");
        }
        vLegend->Draw();
        
        // If we want to draw a decomposition of the fourier fit, do it
        if(drawFitDecomposition){
          
          for(int iFile = 0; iFile < nFiles; iFile++){
            TF1 *fourierFit = longRangeDeltaPhi[iFile][iType][iCentrality][iTrackPt]->GetFunction("fourier");
            TF1 *fitComponents[nFlowComponents];
            int componentColors[] = {kGreen+4, kBlue, kMagenta, kCyan, kViolet};
            for(int iFlow = 0; iFlow < nFlowComponents; iFlow++){
              if(!drawVn[iFlow]) continue;
              fitComponents[iFlow] = new TF1(Form("fv%d", iFlow+1), Form("[0]+[0]*[1]*2.0*TMath::Cos(%d.0*x)", iFlow+1), -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
              fitComponents[iFlow]->SetParameter(0,fourierFit->GetParameter(0));
              fitComponents[iFlow]->SetParameter(1,fourierFit->GetParameter(iFlow+1));
              //fitComponents[iFlow]->SetLineColor(componentColors[iFlow]);
              fitComponents[iFlow]->SetLineColor(lineColors[iFile]);
              fitComponents[iFlow]->Draw("same");
              
            }
          }
    
          
        } // Drawing fit composition
        
        // Save the figures to file
        if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/longRangeDistributionComparison%s%s_C=%.0f-%.0f_T=%.0f-%.0f.%s", saveTypeString[iType].Data(), saveComment.Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1], figureFormat.Data()));
        } // Saving figures
        
        if(drawRatio){
          
          legend = new TLegend(0.17,0.7,0.37,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0, Form("C = %.0f-%.0f%%", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "");
          legend->AddEntry((TObject*) 0, Form("%.1f < p_{T} < %.1f GeV", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]), "");
          
          for(int iFile = 1; iFile < nFiles; iFile++){
            
            if(iFile == 1){
              drawer->DrawHistogram(longRangeDeltaPhiRatio[iFile-1][iType][iCentrality][iTrackPt], Form("%s #Delta#varphi", typeString[iType].Data()), "Ratio", " ");
            } else {
              longRangeDeltaPhiRatio[iFile-1][iType][iCentrality][iTrackPt]->Draw("same");
            }
            
            legend->AddEntry(longRangeDeltaPhiRatio[iFile-1][iType][iCentrality][iTrackPt], Form("%s / %s", legendComment[iFile].Data(), legendComment[0].Data()), "l");
            
          } // Loop over files
          
          legend->Draw();
          
          // Save the figures to file
          if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/longRangeDistributionRatio%s%s_C=%.0f-%.0f_T=%.0f-%.0f.%s", saveTypeString[iType].Data(), saveComment.Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1], figureFormat.Data()));
          } // Saving figures
          
        } // Drawing ratio
        
      } // track pT loop
    } // centrality loop
  } // Jet type loop
  
  
}
