#include "JDrawer.h"  R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"

/*
 * Macro for comparing jet reconstruction bias corrections for long range study
 */
void compareBiasCorrections(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Main files from which the long range asymmetries are obtained
  TString correctionFileName[] = {"corrections/jetReconstructionBiasCorrection_noShiftFitUpToV4_forTestingPurposes.txt", "corrections/jetReconstructionBiasCorrection_noShiftFitUpToV4_forTestingPurposes_old.txt", "", "", ""};
  
  TString correctionLegendName[] = {"New", "Old", "", "", ""};
  
  const int maxFiles = 5;
  int nFilesToCompare = 1;
  for(int iFile = 1; iFile < 5; iFile++){
    if(correctionFileName[iFile] == "") break;
    nFilesToCompare++;
  }

  // Binning information
  const int nCentralityBins = 4;
  const int nTrackPtBins = 7;
  const int nAsymmetryBins = 3;
  const int nFlow = 4;
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  int firstGraphFlowComponent = 2;    // First drawn flow component to comparison graphs
  int lastGraphFlowComponent = 2;     // Last drawn flow component to the comparison graphs
  
  const bool drawGraphFileComparison = true;          // Compare different files
  const bool drawGraphAsymmetryComparison = false;    // Compare xj bins from first file
  const bool drawGraphFlowComparison = false;          // Compare vn components from first file
  
  const bool saveFigures = false;
  TString saveComment = "_dihadronFit";
  
  int firstDrawnAsymmetryBin = nAsymmetryBins;
  int lastDrawnAsymmetryBin = nAsymmetryBins;
  
  // ==================================================================
  // ====================== Configuration done ========================
  // ==================================================================
  
  // Open the pp input file

  
  // Read the uncertainties from the uncertainty file
  JffCorrector *correctionProvider[maxFiles];
  for(int iFile = 0; iFile < nFilesToCompare; iFile++){
    correctionProvider[iFile] = new JffCorrector();
    correctionProvider[iFile]->ReadJetReconstructionBiasFile(correctionFileName[iFile]);
  }
  
  // Define arrays for corrections
  double masterFlowTable[maxFiles][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nFlow];
  double masterFlowError[maxFiles][nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nFlow];
    
  // Read the corrections from the file
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
        for(int iFlow = 0; iFlow < nFlow; iFlow++){
          for(int iFile = 0; iFile < nFilesToCompare; iFile++){
            masterFlowTable[iFile][iAsymmetry][iCentrality][iTrackPt][iFlow] = correctionProvider[iFile]->GetJetReconstructionBiasCorrection(iFlow, iCentrality, iTrackPt, iAsymmetry);;
            masterFlowError[iFile][iAsymmetry][iCentrality][iTrackPt][iFlow] = 0;
            
            
          } // files
        } // flow components
      } // Asymmetry loop
    } // track pT
  } // centrality
    
  // Put the extracted values into graphs
  TGraphErrors *flowGraphPt[maxFiles][nAsymmetryBins+1][nCentralityBins+1][nFlow];
  
  // Arrays for track pT for graphs
  double graphPointsX[nTrackPtBins-2];      // x-axis points in flow graphs
  double graphErrorsX[nTrackPtBins-2];      // No errors for x-axis
  double graphPointsY[nTrackPtBins-2];      // Vn values
  double graphErrorsY[nTrackPtBins-2];      // Statistical errors for Vn
  
  double defaultXpoints[] = {0.85, 1.5, 2.5, 3.5, 6, 10, 14};
  
  // Style settings for graphs
  int markers[] = {kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross};
  int fullMarkers[] = {kFullCircle, kFullSquare, kFullDiamond, kFullCross};
  int colors[] = {kBlue,kRed,kGreen+2,kBlack};
  int flowColors[] = {kBlue, kBlack, kRed, kGreen+3};
  
  for(int iFile = 0; iFile < nFilesToCompare; iFile++){
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins - 2; iTrackPt++){
        
        // Initialize the arrays for graphs
        graphPointsX[iTrackPt] = defaultXpoints[iTrackPt];
        graphErrorsX[iTrackPt] = 0;
        graphPointsY[iTrackPt] = 0;
        graphErrorsY[iTrackPt] = 0;
        
      }
      
      // Create an array for the y-axis and make a graph out of vn values
      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
        for(int iFlow = 0; iFlow < nFlow; iFlow++){
          for(int iTrackPt = 0; iTrackPt < nTrackPtBins-2; iTrackPt++){
            graphPointsY[iTrackPt] = masterFlowTable[iFile][iAsymmetry][iCentrality][iTrackPt][iFlow];
            graphErrorsY[iTrackPt] = masterFlowError[iFile][iAsymmetry][iCentrality][iTrackPt][iFlow];
          }
          
          flowGraphPt[iFile][iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsY, graphErrorsX, graphErrorsY);
          flowGraphPt[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(colors[iAsymmetry]);
          flowGraphPt[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iAsymmetry]);
          
        } // Loop over flow components
      } // Loop over asymmetry bins
    } // Centrality loop
  } // File loop
  
  
  // Configuration for plotting
  
  // How to zoom vn plots    //   0-10 10-30 30-50 50-100  pp
  double vZoomTable[4][5] = {{    0.2,  0.25, 0.35,  1,   1.4},  // v1
                             {    0.12, 0.2,  0.3,  0.5,  0.6},  // v2
                             {    0.12, 0.12, 0.12, 0.2,  0.25},  // v3
                             {    0.1,  0.1,  0.1,  0.1,  0.1}}; // v4
  double multiplier;
  double factors[] = {0.7,0.6,0.4,0.6};
  
  // Helper variables for drawing figures
  TLegend *legend;
  TLegend *vLegend;

  char namerY[100];
  TString legendString;
  
  TString asymmetryString;
  TString compactAsymmetryString;
  
  // Initialize the graph drawer
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceGraph();
  
  // Helper variables for finding track pT information
  // double yZoomForFlowLow[] = {-0.2,-0.005,-0.005,-0.005};
  double yZoomForFlowLow[] = {-0.2,-0.05,-0.005,-0.005};
  double yZoomForFlowHigh[] = {0.06,0.12,0.05,0.02};
  double legendY1, legendY2;
  
  TLine *zeroLine = new TLine(0,0,8,0);
  zeroLine->SetLineStyle(2);
  
  // Draw all selected flow components to the same graph
  if(drawGraphAsymmetryComparison){
    
    // Draw the graphs as a function of pT
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      for(int iFlow = firstGraphFlowComponent-1; iFlow <= lastGraphFlowComponent-1; iFlow++){
        
        legendY1 = 0.6; legendY2 = 0.9;
        if(iFlow == 0){
          legendY1 = 0.2; legendY2 = 0.5;
        }
        
        // Setup a legend for the figure
        legend = new TLegend(0.2,legendY1,0.6,legendY2);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        if(iCentrality == nCentralityBins){
          legend->SetHeader("pp");
        } else {
          legend->SetHeader(Form("PbPb C: %.0f-%.0f %%", centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
        }
        
        // First, draw the systematic uncertainty bands to the canvas
        sprintf(namerY,"V_{%d}",iFlow+1);
        
        // After systematic uncertainties are drawn, draw the points on top of the bands
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
          if(iAsymmetry == firstDrawnAsymmetryBin) {
            drawer->DrawGraph(flowGraphPt[0][firstDrawnAsymmetryBin][iCentrality][iFlow], 0, 8, yZoomForFlowLow[iFlow], yZoomForFlowHigh[iFlow], "Track p_{T} (GeV)", namerY, " ", "p");
          } else {
            flowGraphPt[0][iAsymmetry][iCentrality][iFlow]->Draw("psame");
          }
          if(iAsymmetry == nAsymmetryBins){
            legend->AddEntry(flowGraphPt[0][iAsymmetry][iCentrality][iFlow],"x_{j} > 0.0","p");
          } else {
            legend->AddEntry(flowGraphPt[0][iAsymmetry][iCentrality][iFlow],Form("%.1f < x_{j} < %.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]),"p");
          }
          
        }
        
        legend->Draw();
        
        if(iFlow == 0) zeroLine->Draw();
        
        // Save the figures to file
        if(saveFigures){
          if(iCentrality == nCentralityBins){
            gPad->GetCanvas()->SaveAs(Form("figures/biasCorrectionPt%s_v%d_pp.pdf", saveComment.Data(), iFlow+1));
          } else {
            gPad->GetCanvas()->SaveAs(Form("figures/biasCorrectionPt%s_v%d_C=%.0f-%.0f.pdf", saveComment.Data(), iFlow+1, centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
          }
        } // Saving figures
        
      } // Flow loop
    } // Centrality loop
  } // Draw asymmetry comparison graphs
  
  // Draw all the different asymmetry bins to the same graph
  if(drawGraphFlowComparison){

    // Draw the graphs as a function of pT
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){

        // Set the asymmetry string based on the selected asymmetry bin
        if(iAsymmetry < nAsymmetryBins){
          asymmetryString = Form(", %.1f < x_{j} < %.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
          compactAsymmetryString = Form("_A=%.1f-%.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
          compactAsymmetryString.ReplaceAll(".","v");
        } else {
          asymmetryString = ", x_{j} > 0.0";
          compactAsymmetryString = "";
        }

        legendY1 = 0.6; legendY2 = 0.9;

        // Setup a legend for the figure
        legend = new TLegend(0.2,legendY1,0.6,legendY2);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        if(iCentrality == nCentralityBins){
          legend->SetHeader(Form("pp%s", asymmetryString.Data()));
        } else {
          legend->SetHeader(Form("PbPb C: %.0f-%.0f %%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString.Data()));
        }

        // After systematic uncertainties are drawn, draw the points on top of the bands
        for(int iFlow = firstGraphFlowComponent-1; iFlow <= lastGraphFlowComponent-1; iFlow++){
          flowGraphPt[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(flowColors[iFlow]);
          if(iFlow == firstGraphFlowComponent-1){
            drawer->DrawGraph(flowGraphPt[0][iAsymmetry][iCentrality][firstGraphFlowComponent-1], 0, 8, yZoomForFlowLow[firstGraphFlowComponent-1], yZoomForFlowHigh[firstGraphFlowComponent-1], "Track p_{T} (GeV)", "V_{n}", " ", "p");
          } else {
            flowGraphPt[0][iAsymmetry][iCentrality][iFlow]->Draw("psame");
          }
          legend->AddEntry(flowGraphPt[0][iAsymmetry][iCentrality][iFlow],Form("V_{%d}",iFlow+1),"p");

        }

        legend->Draw();

        // Save the figures to file
        if(saveFigures){
          if(iCentrality == nCentralityBins){
            gPad->GetCanvas()->SaveAs(Form("figures/biasVnComparison%s%s_pp.pdf", saveComment.Data(), compactAsymmetryString.Data()));
          } else {
            gPad->GetCanvas()->SaveAs(Form("figures/biasVnComparison%s%s_C=%.0f-%.0f.pdf", saveComment.Data(), compactAsymmetryString.Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
          }
        } // Saving figures

      } // Asymmetry loop
    } // Centrality loop
  } // Draw flow component comparison
    
  
  // Draw all the different asymmetry bins to the same graph
  if(drawGraphFileComparison){
    
    // Draw the graphs as a function of pT
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){

        // Set the asymmetry string based on the selected asymmetry bin
        if(iAsymmetry < nAsymmetryBins){
          asymmetryString = Form(", %.1f < x_{j} < %.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
          compactAsymmetryString = Form("_A=%.1f-%.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
          compactAsymmetryString.ReplaceAll(".","v");
        } else {
          asymmetryString = ", x_{j} > 0.0";
          compactAsymmetryString = "";
        }

        legendY1 = 0.6; legendY2 = 0.9;

        // Setup a legend for the figure
        legend = new TLegend(0.2,legendY1,0.6,legendY2);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        if(iCentrality == nCentralityBins){
          legend->SetHeader(Form("pp%s", asymmetryString.Data()));
        } else {
          legend->SetHeader(Form("PbPb C: %.0f-%.0f %%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString.Data()));
        }

        // After systematic uncertainties are drawn, draw the points on top of the bands
        for(int iFlow = firstGraphFlowComponent-1; iFlow <= lastGraphFlowComponent-1; iFlow++){
          
          for(int iFile = 0; iFile < nFilesToCompare; iFile++){
          
            flowGraphPt[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(flowColors[iFile]);
            
            if(iFile == 0){
              drawer->DrawGraph(flowGraphPt[iFile][iAsymmetry][iCentrality][firstGraphFlowComponent-1], 0, 8, yZoomForFlowLow[firstGraphFlowComponent-1], yZoomForFlowHigh[firstGraphFlowComponent-1], "Track p_{T} (GeV)", "V_{n}", " ", "p");
            } else {
              flowGraphPt[iFile][iAsymmetry][iCentrality][iFlow]->Draw("psame"); 
            }
            
            legend->AddEntry(flowGraphPt[iFile][iAsymmetry][iCentrality][iFlow], correctionLegendName[iFile], "p");

          } // File loop
          
          legend->Draw();
          
          // Save the figures to file
          if(saveFigures){
            if(iCentrality == nCentralityBins){
              gPad->GetCanvas()->SaveAs(Form("figures/biasComparison%s_v%d%s_pp.pdf", saveComment.Data(), iFlow+1, compactAsymmetryString.Data()));
            } else {
              gPad->GetCanvas()->SaveAs(Form("figures/biasComparison%s_v%d%s_C=%.0f-%.0f.pdf", saveComment.Data(), iFlow+1, compactAsymmetryString.Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
            }
          } // Saving figures
          
        } // Flow component loop

      } // Asymmetry loop
    } // Centrality loop
  } // Draw file comparison
  
}

