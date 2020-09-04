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
*  - different stages of the analysis in the same plot
*  - different asymmetry bins in the same plot
*  - different flow components in the same plot
*  - chosen histogram from different files in the same plot
*/
void longRangeGraphPlotter(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Main files from which the long range asymmetries are obtained
  const int maxFiles = 6;
  TString directoryName = "flowGraphs/";
  TString graphFileName = "testDihadron_sameEvent_allNormQ.root";
  TFile *graphFile[maxFiles];
  graphFile[0] = TFile::Open(directoryName+graphFileName);
  
  // Other files whose results can be compared with the nominal file
  int nComparisonFiles = 3;
  TString comparisonFileName[] = {"testDihadron_sameEvent_lowNormQ_largeGap.root", "testDihadron_sameEvent_highNormQ_largeGap.root", "flowGraphs_PbPbData_noJetReconstructionCorrection_fullDihadronStats.root", "finalGraphTestNew.root", ""};
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    graphFile[iFile+1] = TFile::Open(directoryName+comparisonFileName[iFile]);
  }
  
  // Legend text given to each compared file
  TString fileLegend[] = {"All Q", "Low Q", "High Q", "Data", "Fifth file"};
  
  const int nCentralityBins = 3;
  const int nTrackPtBins = 7;
  const int nAsymmetryBins = 3;
  const int nFlowComponents = 4;
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj

  // Plots to be drawn from the main file
  const bool drawGraphAsymmetryComparison = false;    // Draw all selected asymmetry bins to the same graph
  const bool drawGraphVnComparison = false;           // Draw selected flow components to the same graph
  const bool drawGraphStages = false;                 // Draw all intermediate steps leading to jet vn
  
  // Plots to be compared between files
  const bool drawJetHadronVnFileComparison = false;
  const bool drawDihadronVnFileComparison = true;
  const bool drawHadronVnFileComparison = true;
  const bool drawJetVnFileComparison = false;
  const bool drawFileComparison = drawJetHadronVnFileComparison || drawDihadronVnFileComparison || drawHadronVnFileComparison || drawJetVnFileComparison;
  
  const bool drawSystematicUncertainties = false;     // Include systematic uncertainties in the plots
  
  const bool saveFigures = false;                     // Save the figures in a file
  TString saveComment = "_normalizedQvectorComparisonLargeGap";              // String to be added to saved file names
  
  int firstDrawnAsymmetryBin = nAsymmetryBins;
  int lastDrawnAsymmetryBin = nAsymmetryBins;
  
  int firstDrawnVn = 2;
  int lastDrawnVn = 2;
  
  double maxTrackPt = 4.5;
  
  // Load the graphs from the input file
  TGraphErrors *flowGraphJetHadron[maxFiles][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowGraphJetHadronCorrected[maxFiles][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowGraphDihadron[maxFiles][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowGraphHadron[maxFiles][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowGraphJet[maxFiles][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  
  TGraphErrors *flowSystematicsJetHadron[maxFiles][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowSystematicsJetHadronCorrected[maxFiles][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowSystematicsDihadron[maxFiles][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowSystematicsHadron[maxFiles][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowSystematicsJet[maxFiles][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  
  char histogramNamer[150];
  int nPoints;
  double *xAxisPoints;
  
  for(int iFile = 0; iFile < nComparisonFiles+1; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
        for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
          
          sprintf(histogramNamer,"jetHadronVn/jetHadronV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
          flowGraphJetHadron[iFile][iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          sprintf(histogramNamer,"jetHadronVn/jetHadronV%dSystematics_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
          flowSystematicsJetHadron[iFile][iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          sprintf(histogramNamer,"jetHadronVnCorrected/jetHadronV%dCorrected_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
          flowGraphJetHadronCorrected[iFile][iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          sprintf(histogramNamer,"jetHadronVnCorrected/jetHadronV%dCorrectedSystematics_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
          flowSystematicsJetHadronCorrected[iFile][iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          sprintf(histogramNamer,"dihadronVn/dihadronV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
          flowGraphDihadron[iFile][iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          sprintf(histogramNamer,"dihadronVn/dihadronV%dSystematics_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
          flowSystematicsDihadron[iFile][iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          sprintf(histogramNamer,"hadronVn/hadronV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
          flowGraphHadron[iFile][iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          sprintf(histogramNamer,"hadronVn/hadronV%dSystematics_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
          flowSystematicsHadron[iFile][iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          sprintf(histogramNamer,"jetVn/jetV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
          flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          sprintf(histogramNamer,"jetVn/jetV%dSystematics_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
          flowSystematicsJet[iFile][iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          
          // Remove the points that are above the maximum pT limit from the graphs
          xAxisPoints = flowGraphJetHadron[iFile][iAsymmetry][iCentrality][iFlow]->GetX();
          nPoints = flowGraphJetHadron[iFile][iAsymmetry][iCentrality][iFlow]->GetN();
          for(int iX = nPoints-1; iX >= 0; iX--){
            if(xAxisPoints[iX] > maxTrackPt){
              flowGraphJetHadron[iFile][iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
              flowGraphJetHadronCorrected[iFile][iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
              flowGraphDihadron[iFile][iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
              flowGraphHadron[iFile][iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
              flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
              
              flowSystematicsJetHadron[iFile][iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
              flowSystematicsJetHadronCorrected[iFile][iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
              flowSystematicsDihadron[iFile][iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
              flowSystematicsHadron[iFile][iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
              flowSystematicsJet[iFile][iAsymmetry][iCentrality][iFlow]->RemovePoint(iX);
            } else {
              break;
            }
          }
          
        } // Flow component loop
      } // Asymmetry loop
    } // Centrality loop
  } // File loop
  
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
  int fileColors[] = {kBlack, kBlue, kRed, kGreen+3, kMagenta, kCyan};
  TString asymmetryString[] = {" 0.0 < x_{j} < 0.6", " 0.6 < x_{j} < 0.8", " 0.8 < x_{j} < 1.0", ""};
  TString asymmetryLegend[] = {"0.0 < x_{j} < 0.6", "0.6 < x_{j} < 0.8", "0.8 < x_{j} < 1.0", "x_{j} integrated"};
  TString compactAsymmetryString[] = {"_A=0v0-0v6", "_A=0v6-0v8", "_A=0v8-1v0", ""};

  TLine *atlasV2;
  double atlasV2Number[] = {0.018, 0.03, 0.035, 0.03};
  
  TString legendString;
  char namerY[100];
  
  TLine *zeroLine = new TLine(0,0,maxTrackPt,0);
  zeroLine->SetLineStyle(2);
  
  double minZoom, maxZoom;
  
  // Draw graph with all the stages leading to final jet vn visible
  if(drawGraphStages){
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
          
          sprintf(namerY,"V_{%d}",iFlow+1);
          flowGraphJetHadron[0][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(kFullCircle);
          flowGraphJetHadron[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(kBlack);
          drawer->DrawGraph(flowGraphJetHadron[0][iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, 0, 0.35, "Track p_{T} (GeV)", namerY, " ", "p");
          
          flowGraphJetHadronCorrected[0][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(kFullSquare);
          flowGraphJetHadronCorrected[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(kGreen+3);
          flowGraphJetHadronCorrected[0][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
          
          flowGraphDihadron[0][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(kFullCross);
          flowGraphDihadron[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(kRed);
          flowGraphDihadron[0][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
          
          flowGraphHadron[0][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(kFullDiamond);
          flowGraphHadron[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(kMagenta);
          flowGraphHadron[0][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
          
          flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(kFullStar);
          flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(kBlue);
          flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
          
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
          legend->AddEntry(flowGraphJetHadron[0][iAsymmetry][iCentrality][iFlow], Form("Jet-track V_{%d}", iFlow+1), "p");
          legend->AddEntry(flowGraphJetHadronCorrected[0][iAsymmetry][iCentrality][iFlow], Form("Jet-track V_{%d}, corrected", iFlow+1), "p");
          legend->AddEntry(flowGraphDihadron[0][iAsymmetry][iCentrality][iFlow], Form("Dihadron V_{%d}", iFlow+1), "p");
          legend->AddEntry(flowGraphHadron[0][iAsymmetry][iCentrality][iFlow], Form("Hadron v_{%d}", iFlow+1), "p");
          legend->AddEntry(flowGraphJet[0][iAsymmetry][iCentrality][iFlow], Form("Jet v_{%d}", iFlow+1), "p");
          
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
          flowSystematicsJet[0][firstDrawnAsymmetryBin][iCentrality][iFlow]->SetFillColorAlpha(colors[firstDrawnAsymmetryBin],0.25);
          drawer->DrawGraph(flowSystematicsJet[0][firstDrawnAsymmetryBin][iCentrality][iFlow], 0, maxTrackPt, -0.05, 0.35, "Track p_{T} (GeV)", Form("Jet v_{%d}", iFlow+1), " ", "2,same");
          
          for(int iAsymmetry = firstDrawnAsymmetryBin+1; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
            flowSystematicsJet[0][iAsymmetry][iCentrality][iFlow]->SetFillColorAlpha(colors[iAsymmetry],0.25);
            flowSystematicsJet[0][iAsymmetry][iCentrality][iFlow]->Draw("2,same");
          }
        }
        
        // After systematic uncertainties are drawn, draw the points on top of the bands
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
          flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(colors[iAsymmetry]);
          flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iAsymmetry]);
          if(iAsymmetry == firstDrawnAsymmetryBin && !drawSystematicUncertainties){
            drawer->DrawGraph(flowGraphJet[0][iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, -0.05, 0.35, "Track p_{T} (GeV)", Form("Jet v_{%d}", iFlow+1), " ", "p");
          } else {
            flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->Draw("psame");
          }
          legend->AddEntry(flowGraphJet[0][iAsymmetry][iCentrality][iFlow],asymmetryLegend[iAsymmetry],"p");
          
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
          flowSystematicsJet[0][iAsymmetry][iCentrality][firstDrawnVn-1]->SetFillColorAlpha(flowColors[firstDrawnVn-1],0.25);
          drawer->DrawGraph(flowSystematicsJet[0][iAsymmetry][iCentrality][firstDrawnVn-1], 0, maxTrackPt, -0.05, 0.35, "Track p_{T} (GeV)", "v_{n}", " ", "2,same");
          
          for(int iFlow = firstDrawnVn; iFlow <= lastDrawnVn-1; iFlow++){
            flowSystematicsJet[0][iAsymmetry][iCentrality][iFlow]->SetFillColorAlpha(flowColors[iFlow],0.25);
            flowSystematicsJet[0][iAsymmetry][iCentrality][iFlow]->Draw("2,same");
          }
        }
        
        // After systematic uncertainties are drawn, draw the points on top of the bands
        for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
          flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(flowColors[iFlow]);
          flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iFlow]);
          if(iFlow == firstDrawnVn-1 && !drawSystematicUncertainties){
            drawer->DrawGraph(flowGraphJet[0][iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, -0.05, 0.35, "Track p_{T} (GeV)", "v_{n}", " ", "p");
          } else {
            flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->Draw("psame");
          }
          legend->AddEntry(flowGraphJet[0][iAsymmetry][iCentrality][iFlow],Form("Jet v_{%d}",iFlow+1),"p");
          
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
  
  // Compare graphs from different files
  if(drawFileComparison){
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
          
          if(drawJetHadronVnFileComparison){
          
            sprintf(namerY,"Jet-hadron V_{%d}",iFlow+1);
            legend = new TLegend(0.2,0.55,0.5,0.9);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->SetHeader(Form("Cent: %.0f-%.0f%%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
            
            for(int iFile = 0; iFile < nComparisonFiles+1; iFile++){
              flowGraphJetHadron[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iFile]);
              flowGraphJetHadron[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(fileColors[iFile]);
              flowGraphJetHadron[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerSize(1.3);
              if(iFile == 0){
                 drawer->DrawGraph(flowGraphJetHadron[iFile][iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, 0, 0.08, "Track p_{T} (GeV)", namerY, " ", "p");
              } else {
                flowGraphJetHadron[iFile][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
              }
              
              legend->AddEntry(flowGraphJetHadron[iFile][iAsymmetry][iCentrality][iFlow], fileLegend[iFile], "p");
              
            }
            
            legend->Draw();
            
            // Save the figures to file
            if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/jetHadronV%dComparison%s%s_C=%.0f-%.0f.pdf", iFlow+1, saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
            }
            
          } // File comparison for jet-hadron Vn
          
          if(drawDihadronVnFileComparison){
          
            sprintf(namerY,"Dihadron V_{%d}",iFlow+1);
            legend = new TLegend(0.2,0.55,0.5,0.9);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->SetHeader(Form("Cent: %.0f-%.0f%%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
            
            for(int iFile = 0; iFile < nComparisonFiles+1; iFile++){
              flowGraphDihadron[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iFile]);
              flowGraphDihadron[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(fileColors[iFile]);
              flowGraphDihadron[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerSize(1.3);
              if(iFile == 0){
                 drawer->DrawGraph(flowGraphDihadron[iFile][iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, 0, 0.08, "Track p_{T} (GeV)", namerY, " ", "p");
              } else {
                flowGraphDihadron[iFile][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
              }
              
              legend->AddEntry(flowGraphDihadron[iFile][iAsymmetry][iCentrality][iFlow], fileLegend[iFile], "p");
              
            }
            
            legend->Draw();
            
            // Save the figures to file
            if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/dihadronV%dComparison%s%s_C=%.0f-%.0f.pdf", iFlow+1, saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
            }
            
          } // File comparison for jet-hadron Vn
          
          if(drawHadronVnFileComparison){
          
            sprintf(namerY,"Hadron v_{%d}",iFlow+1);
            legend = new TLegend(0.2,0.55,0.5,0.9);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->SetHeader(Form("Cent: %.0f-%.0f%%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
            
            for(int iFile = 0; iFile < nComparisonFiles+1; iFile++){
              flowGraphHadron[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iFile]);
              flowGraphHadron[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(fileColors[iFile]);
              flowGraphHadron[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerSize(1.3);
              if(iFile == 0){
                 drawer->DrawGraph(flowGraphHadron[iFile][iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, 0, 0.3, "Track p_{T} (GeV)", namerY, " ", "p");
              } else {
                flowGraphHadron[iFile][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
              }
              
              legend->AddEntry(flowGraphHadron[iFile][iAsymmetry][iCentrality][iFlow], fileLegend[iFile], "p");
              
            }
            
            legend->Draw();
            
            // Save the figures to file
            if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/hadronV%dComparison%s%s_C=%.0f-%.0f.pdf", iFlow+1, saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
            }
            
          } // File comparison for jet-hadron Vn
          
          if(drawJetVnFileComparison){
          
            sprintf(namerY,"Jet v_{%d}",iFlow+1);
            legend = new TLegend(0.2,0.55,0.5,0.9);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->SetHeader(Form("Cent: %.0f-%.0f%%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
            
            minZoom = 0;
            maxZoom = 0.25;
            
            if(iCentrality == 0){
              minZoom = -0.05;
              maxZoom = 0.2;
            }
            
            for(int iFile = 0; iFile < nComparisonFiles+1; iFile++){
              flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iFile]);
              flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(fileColors[iFile]);
              flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerSize(1.3);
              if(iFile == 0){
                 drawer->DrawGraph(flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, minZoom, maxZoom, "Track p_{T} (GeV)", namerY, " ", "p");
              } else {
                flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
              }
              
              legend->AddEntry(flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow], fileLegend[iFile], "p");
              
            }
            
            if(iFlow == 1){
              atlasV2 = new TLine(0.5,atlasV2Number[iCentrality],maxTrackPt-0.5,atlasV2Number[iCentrality]);
              atlasV2->SetLineStyle(2);
              atlasV2->SetLineColor(kBlue);
              atlasV2->Draw();
              legend->AddEntry(atlasV2, "ATLAS jet v_{2}", "l");
            }
            
            legend->Draw();
            
            // Save the figures to file
            if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/jetV%dComparison%s%s_C=%.0f-%.0f.pdf", iFlow+1, saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
            }
            
          } // File comparison for jet vn
          
        } // Asymmetry loop
      } // Flow component loop
    } // Centrality loop
    
  } // Draw stages leading to jet vn
  
}

