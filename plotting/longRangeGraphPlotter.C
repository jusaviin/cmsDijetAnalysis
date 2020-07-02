#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "DijetMethods.h"
#include "JDrawer.h"
#include "JffCorrector.h"

/*
* Macro for plotting long range correlation results.
* The results are fully corrected down to jet vn level.
* Systematic uncertainties can be plotted together with data points.
*
*  Drawing styles:
*  - different asymmetry bins in the same plot
*  - different flow components in the same plot
*/
void longRangeGraphPlotter(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Main files from which the long range asymmetries are obtained
  TString graphFileName = "finalGraphTestWithXj.root";
  TFile *graphFile = TFile::Open(graphFileName);
  
  const int nCentralityBins = 4;
  const int nTrackPtBins = 7;
  const int nAsymmetryBins = 3;
  const int nFlowComponents = 4;
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj

  const bool drawGraphAsymmetryComparison = false;    // Draw all selected asymmetry bins to the same graph
  const bool drawGraphVnComparison = false;           // Draw selected flow components to the same graph
  const bool drawGraphStages = true;                 // Draw all intermediate steps leading to jet vn
  
  const bool drawSystematicUncertainties = false;     // Include systematic uncertainties in the plots
  
  const bool saveFigures = false;                     // Save the figures in a file
  TString saveComment = "_initialCheck";              // String to be added to saved file names
  
  int firstDrawnAsymmetryBin = 2;
  int lastDrawnAsymmetryBin = 2;
  
  int firstDrawnVn = 2;
  int lastDrawnVn = 2;
  
  double maxTrackPt = 4.5;
  
  // Load the graphs from the input file
  TGraphErrors *flowGraphJetHadron[nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowGraphJetHadronCorrected[nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowGraphDihadron[nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowGraphHadron[nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowGraphJet[nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  
  TGraphErrors *flowSystematicsJetHadron[nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowSystematicsJetHadronCorrected[nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowSystematicsDihadron[nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowSystematicsHadron[nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowSystematicsJet[nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  
  char histogramNamer[150];
  int nPoints;
  double *xAxisPoints;
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
      for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
        
        sprintf(histogramNamer,"jetHadronVn/jetHadronV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
        flowGraphJetHadron[iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile->Get(histogramNamer);
        
        sprintf(histogramNamer,"jetHadronVn/jetHadronV%dSystematics_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
        flowSystematicsJetHadron[iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile->Get(histogramNamer);
        
        sprintf(histogramNamer,"jetHadronVnCorrected/jetHadronV%dCorrected_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
        flowGraphJetHadronCorrected[iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile->Get(histogramNamer);
        
        sprintf(histogramNamer,"jetHadronVnCorrected/jetHadronV%dCorrectedSystematics_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
        flowSystematicsJetHadronCorrected[iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile->Get(histogramNamer);
        
        sprintf(histogramNamer,"dihadronVn/dihadronV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
        flowGraphDihadron[iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile->Get(histogramNamer);
        
        sprintf(histogramNamer,"dihadronVn/dihadronV%dSystematics_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
        flowSystematicsDihadron[iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile->Get(histogramNamer);
        
        sprintf(histogramNamer,"hadronVn/hadronV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
        flowGraphHadron[iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile->Get(histogramNamer);
        
        sprintf(histogramNamer,"hadronVn/hadronV%dSystematics_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
        flowSystematicsHadron[iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile->Get(histogramNamer);
        
        sprintf(histogramNamer,"jetVn/jetV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
        flowGraphJet[iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile->Get(histogramNamer);
        
        sprintf(histogramNamer,"jetVn/jetV%dSystematics_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
        flowSystematicsJet[iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile->Get(histogramNamer);
        
        
        // Remove the points that are above the maximum pT limit from the graphs
        xAxisPoints = flowGraphJetHadron[iAsymmetry][iCentrality][iFlow]->GetX();
        nPoints = flowGraphJetHadron[iAsymmetry][iCentrality][iFlow]->GetN();
        for(int iX = nPoints-1; iX >= 0; iX--){
          if(xAxisPoints[iX] > maxTrackPt){
            flowGraphJetHadron[iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
            flowGraphJetHadronCorrected[iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
            flowGraphDihadron[iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
            flowGraphHadron[iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
            flowGraphJet[iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
            
            flowSystematicsJetHadron[iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
            flowSystematicsJetHadronCorrected[iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
            flowSystematicsDihadron[iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
            flowSystematicsHadron[iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
            flowSystematicsJet[iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
          } else {
            break;
          }
        }
        
      } // Flow component loop
    } // Asymmetry loop
  } // Centrality loop
  
  // Setup the drawer for graphs
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceGraph();
  
  
  // How to zoom vn plots    //   0-10 10-30 30-50 50-100  pp
  double vZoomTable[4][5] = {{    0.2,  0.25, 0.35,  1,   1.4},  // v1
                             {    0.12, 0.2,  0.3,  0.5,  0.6}, // v2
                             {    0.12, 0.12, 0.12, 0.2,  0.25},  // v3
                             {    0.1,  0.1,  0.1,  0.1,  0.1}}; // v4

  // Helper variables for drawing figures
  TLegend *legend;
  TLegend *vLegend;
  int markers[] = {kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross, kOpenStar};
  int fullMarkers[] = {kFullCircle, kFullSquare, kFullDiamond, kFullCross, kFullStar};
  int colors[] = {kBlue,kRed,kGreen+2,kBlack, kMagenta};
  int flowColors[] = {kBlue, kBlack, kRed, kGreen+3, kMagenta};
  TString asymmetryString[] = {" 0.0 < x_{j} < 0.6", " 0.6 < x_{j} < 0.8", " 0.8 < x_{j} < 1.0", ""};
  TString asymmetryLegend[] = {"0.0 < x_{j} < 0.6", "0.6 < x_{j} < 0.8", "0.8 < x_{j} < 1.0", "x_{j} integrated"};
  TString compactAsymmetryString[] = {"_A=0v0-0v6", "_A=0v6-0v8", "_A=0v8-1v0", ""};

  TLine *atlasV2;
  double atlasV2Number[] = {0.018, 0.03, 0.035, 0.03};
  
  TString legendString;
  char namerY[100];
  
  TLine *zeroLine = new TLine(0,0,maxTrackPt,0);
  zeroLine->SetLineStyle(2);
  
  // Draw graph with all the stages leading to final jet vn visible
  if(drawGraphStages){
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
          
          sprintf(namerY,"V_{%d}",iFlow+1);
          flowGraphJetHadron[iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(kFullCircle);
          flowGraphJetHadron[iAsymmetry][iCentrality][iFlow]->SetMarkerColor(kBlack);
          drawer->DrawGraph(flowGraphJetHadron[iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, 0, 0.35, "Track p_{T} (GeV)", namerY, " ", "p");
          
          flowGraphJetHadronCorrected[iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(kFullSquare);
          flowGraphJetHadronCorrected[iAsymmetry][iCentrality][iFlow]->SetMarkerColor(kGreen+3);
          flowGraphJetHadronCorrected[iAsymmetry][iCentrality][iFlow]->Draw("p,same");
          
          flowGraphDihadron[iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(kFullCross);
          flowGraphDihadron[iAsymmetry][iCentrality][iFlow]->SetMarkerColor(kRed);
          flowGraphDihadron[iAsymmetry][iCentrality][iFlow]->Draw("p,same");
          
          flowGraphHadron[iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(kFullDiamond);
          flowGraphHadron[iAsymmetry][iCentrality][iFlow]->SetMarkerColor(kMagenta);
          flowGraphHadron[iAsymmetry][iCentrality][iFlow]->Draw("p,same");
          
          flowGraphJet[iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(kFullStar);
          flowGraphJet[iAsymmetry][iCentrality][iFlow]->SetMarkerColor(kBlue);
          flowGraphJet[iAsymmetry][iCentrality][iFlow]->Draw("p,same");
          
          if(iFlow == 1){
            atlasV2 = new TLine(0.5,atlasV2Number[iCentrality],maxTrackPt-0.5,atlasV2Number[iCentrality]);
            atlasV2->SetLineStyle(2);
            atlasV2->SetLineColor(kBlue);
            atlasV2->Draw();
          }
          
          // Add a legend to the figure
          legend = new TLegend(0.2,0.55,0.5,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->SetHeader(Form("Cent: %.0f-%.0f%%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
          legend->AddEntry(flowGraphJetHadron[iAsymmetry][iCentrality][iFlow], Form("Jet-track V_{%d}", iFlow+1), "p");
          legend->AddEntry(flowGraphJetHadronCorrected[iAsymmetry][iCentrality][iFlow], Form("Jet-track V_{%d}, corrected", iFlow+1), "p");
          legend->AddEntry(flowGraphDihadron[iAsymmetry][iCentrality][iFlow], Form("Dihadron V_{%d}", iFlow+1), "p");
          legend->AddEntry(flowGraphHadron[iAsymmetry][iCentrality][iFlow], Form("Hadron v_{%d}", iFlow+1), "p");
          legend->AddEntry(flowGraphJet[iAsymmetry][iCentrality][iFlow], Form("Jet v_{%d}", iFlow+1), "p");
          
          if(iFlow == 1){
            legend->AddEntry(atlasV2, "ATLAS jet v_{2}", "l");
          }
          
          legend->Draw();
          
          // Save the figures to file
          if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/flowStagesV%d%s%s_C=%.0f-%.0f.pdf", iFlow+1, saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
          } // Saving figures
          
        } // Asymmetry loop
      } // Flow component loop
    } // Centrality loop
    
  } // Draw stages leading to jet vn

  // Draw all the different asymmetry bins in the same plot
  if(drawGraphAsymmetryComparison){
    
    // Draw the graphs as a function of pT
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
        
        // Setup a legend for the figure
        legend = new TLegend(0.2,0.6,0.6,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->SetHeader(Form("PbPb C: %.0f-%.0f %%", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
        
        // First, draw the systematic uncertainty bands to the canvas
        if(drawSystematicUncertainties){
          flowSystematicsJet[firstDrawnAsymmetryBin][iCentrality][iFlow]->SetFillColorAlpha(colors[firstDrawnAsymmetryBin],0.25);
          drawer->DrawGraph(flowSystematicsJet[firstDrawnAsymmetryBin][iCentrality][iFlow], 0, maxTrackPt, -0.05, 0.35, "Track p_{T} (GeV)", Form("Jet v_{%d}", iFlow+1), " ", "2,same");
          
          for(int iAsymmetry = firstDrawnAsymmetryBin+1; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
            flowSystematicsJet[iAsymmetry][iCentrality][iFlow]->SetFillColorAlpha(colors[iAsymmetry],0.25);
            flowSystematicsJet[iAsymmetry][iCentrality][iFlow]->Draw("2,same");
          }
        }
        
        // After systematic uncertainties are drawn, draw the points on top of the bands
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
          flowGraphJet[iAsymmetry][iCentrality][iFlow]->SetMarkerColor(colors[iAsymmetry]);
          flowGraphJet[iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iAsymmetry]);
          if(iAsymmetry == firstDrawnAsymmetryBin && !drawSystematicUncertainties){
            drawer->DrawGraph(flowGraphJet[iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, -0.05, 0.35, "Track p_{T} (GeV)", Form("Jet v_{%d}", iFlow+1), " ", "p");
          } else {
            flowGraphJet[iAsymmetry][iCentrality][iFlow]->Draw("psame");
          }
          legend->AddEntry(flowGraphJet[iAsymmetry][iCentrality][iFlow],asymmetryLegend[iAsymmetry],"p");
          
          zeroLine->Draw();
        }
        
        legend->Draw();
        
        // Save the figures to file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/flowAsymmetryComparison%s%s_C=%.0f-%.0f.pdf", Form("V%d", iFlow+1), saveComment.Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
          
        } // Saving figures
        
      } // Asymmetry loop
    } // Centrality loop
  } // Draw flow component comparison graphs
  
  // Draw all the different Vn components in the same plot
  if(drawGraphVnComparison){
    
    // Draw the graphs as a function of pT
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
        
        // Setup a legend for the figure
        legend = new TLegend(0.2,0.6,0.6,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->SetHeader(Form("PbPb C: %.0f-%.0f %%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
        
        // First, draw the systematic uncertainty bands to the canvas
        if(drawSystematicUncertainties){
          flowSystematicsJet[iAsymmetry][iCentrality][firstDrawnVn-1]->SetFillColorAlpha(flowColors[firstDrawnVn-1],0.25);
          drawer->DrawGraph(flowSystematicsJet[iAsymmetry][iCentrality][firstDrawnVn-1], 0, maxTrackPt, -0.05, 0.35, "Track p_{T} (GeV)", "v_{n}", " ", "2,same");
          
          for(int iFlow = firstDrawnVn; iFlow <= lastDrawnVn-1; iFlow++){
            flowSystematicsJet[iAsymmetry][iCentrality][iFlow]->SetFillColorAlpha(flowColors[iFlow],0.25);
            flowSystematicsJet[iAsymmetry][iCentrality][iFlow]->Draw("2,same");
          }
        }
        
        // After systematic uncertainties are drawn, draw the points on top of the bands
        for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
          flowGraphJet[iAsymmetry][iCentrality][iFlow]->SetMarkerColor(flowColors[iFlow]);
          flowGraphJet[iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iFlow]);
          if(iFlow == firstDrawnVn-1 && !drawSystematicUncertainties){
            drawer->DrawGraph(flowGraphJet[iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, -0.05, 0.35, "Track p_{T} (GeV)", "v_{n}", " ", "p");
          } else {
            flowGraphJet[iAsymmetry][iCentrality][iFlow]->Draw("psame");
          }
          legend->AddEntry(flowGraphJet[iAsymmetry][iCentrality][iFlow],Form("Jet v_{%d}",iFlow+1),"p");
          
          zeroLine->Draw();
        }
        
        legend->Draw();
        
        // Save the figures to file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/flowVnComparison%s%s_C=%.0f-%.0f.pdf", saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
          
        } // Saving figures
        
      } // Asymmetry loop
    } // Centrality loop
  } // Draw flow component comparison graphs
  
}

