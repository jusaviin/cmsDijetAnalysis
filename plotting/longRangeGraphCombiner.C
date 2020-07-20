#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "DijetMethods.h"
#include "JDrawer.h"
#include "JffCorrector.h"

/*
* Macro for plotting different presentation of the input graphs.
*
*  Drawing styles:
*  - Centrality in the x-axis, combine results from different files
*/
void longRangeGraphCombiner(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Main files from which the long range asymmetries are obtained
  const int maxFiles = 6;
  TString directoryName = "flowGraphs/";
  TString graphFileName = "exampleDihadronFromMC.root";//"flowGraphs_PbPbMC2018_RecoGen_noCentShift.root";
  TFile *graphFile[maxFiles];
  graphFile[0] = TFile::Open(directoryName+graphFileName);
  
  // Other files whose results can be compared with the nominal file
  int nComparisonFiles = 2;
  TString comparisonFileName[] = {"flowGraphs_PbPbMC2018_RecoGen_5pCentShift.root", "flowGraphs_PbPbData_noJetReconstructionCorrection_fullDihadronStats.root", "finalGraphTest_noThirdJet_jetLevel_noAsymmetry.root", "finalGraphTestNew.root", ""};
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    graphFile[iFile+1] = TFile::Open(directoryName+comparisonFileName[iFile]);
  }
  
  // Legend text given to each compared file
  TString fileLegend[] = {"Nominal", "No rebin", "Third jet cut", "Fourth file", "Fifth file"};
  
  const int nCentralityBins = 4;
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
  const int nCentralityBinsToCompare = 3;      // How many centrality bins are taken into account in the comparison
  const int ptBinForCentrality = 3;           // Which track pT bin is selected for each pT bin
  const bool drawJetHadronVnCentralityComparison = false;
  const bool drawDihadronVnCentralityComparison = true;
  const bool drawHadronVnCentralityComparison = true;
  const bool drawJetVnCentralityComparison = true;
  const bool drawCentralityComparison = drawJetHadronVnCentralityComparison || drawDihadronVnCentralityComparison || drawHadronVnCentralityComparison || drawJetVnCentralityComparison;
    
  const bool saveFigures = false;                     // Save the figures in a file
  TString saveComment = "_mcComparison";              // String to be added to saved file names
  
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
  
  char histogramNamer[150];
  int nPoints;
  double *xAxisPoints;
  
  for(int iFile = 0; iFile < nComparisonFiles+1; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
        for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
          
          sprintf(histogramNamer,"jetHadronVn/jetHadronV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
          flowGraphJetHadron[iFile][iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          sprintf(histogramNamer,"jetHadronVnCorrected/jetHadronV%dCorrected_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
          flowGraphJetHadronCorrected[iFile][iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          sprintf(histogramNamer,"dihadronVn/dihadronV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
          flowGraphDihadron[iFile][iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          sprintf(histogramNamer,"hadronVn/hadronV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
          flowGraphHadron[iFile][iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          sprintf(histogramNamer,"jetVn/jetV%d_A%dC%d", iFlow+1, iAsymmetry, iCentrality);
          flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          
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
              
            } else {
              break;
            }
          }
          
        } // Flow component loop
      } // Asymmetry loop
    } // Centrality loop
  } // File loop
  
  // Construct Vn graphs as a function of centrality. Combine different centrality bins from different files
  double xValuesCentrality[] = {1,2,3,4,5,6,7,8};
  double xErrorsCentrality[] = {0,0,0,0,0,0,0,0};
  double yValuesCentrality[] = {1,2,3,4,5,6,7,8};
  double yErrorsCentrality[] = {0,0,0,0,0,0,0,0};
  double xValuesCentralityData[] = {1,3,5,7};
  double xErrorsCentralityData[] = {0,0,0,0};
  double yValuesCentralityData[] = {1,3,5,7};
  double yErrorsCentralityData[] = {0,0,0,0};
  
  TString binLabels[] = {"0-10%","5-15%","10-30%","15-35%","30-50%","35-55%","50-90%","55-95%"};
  
  int binIndex;
  double *yPointer;
  
  TGraphErrors *flowGraphJetHadronCentrality[2][nAsymmetryBins+1][nFlowComponents];
  TGraphErrors *flowGraphJetHadronCorrectedCentrality[2][nAsymmetryBins+1][nFlowComponents];
  TGraphErrors *flowGraphDihadronCentrality[2][nAsymmetryBins+1][nFlowComponents];
  TGraphErrors *flowGraphHadronCentrality[2][nAsymmetryBins+1][nFlowComponents];
  TGraphErrors *flowGraphJetCentrality[2][nAsymmetryBins+1][nFlowComponents];
  
  for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      
      // Construct the jet-hadron Vn centrality graphs from three files, MC with nominal centrality, MC with 5 % shifted centrality and data
      for(int iCentrality = 0; iCentrality < nCentralityBinsToCompare; iCentrality++){
        yPointer = flowGraphJetHadron[0][iAsymmetry][iCentrality][iFlow]->GetY();
        yValuesCentrality[iCentrality*2] = yPointer[ptBinForCentrality];
        yErrorsCentrality[iCentrality*2] = flowGraphJetHadron[0][iAsymmetry][iCentrality][iFlow]->GetErrorY(ptBinForCentrality);
        
        yPointer = flowGraphJetHadron[1][iAsymmetry][iCentrality][iFlow]->GetY();
        yValuesCentrality[iCentrality*2+1] = yPointer[ptBinForCentrality];
        yErrorsCentrality[iCentrality*2+1] = flowGraphJetHadron[1][iAsymmetry][iCentrality][iFlow]->GetErrorY(ptBinForCentrality);
        
        yPointer = flowGraphJetHadron[2][iAsymmetry][iCentrality][iFlow]->GetY();
        yValuesCentralityData[iCentrality] = yPointer[ptBinForCentrality];
        yErrorsCentralityData[iCentrality] = flowGraphJetHadron[2][iAsymmetry][iCentrality][iFlow]->GetErrorY(ptBinForCentrality);
      }
      
      flowGraphJetHadronCentrality[0][iAsymmetry][iFlow] = new TGraphErrors(2*nCentralityBinsToCompare, xValuesCentrality, yValuesCentrality, xErrorsCentrality, yErrorsCentrality);
      flowGraphJetHadronCentrality[1][iAsymmetry][iFlow] = new TGraphErrors(nCentralityBinsToCompare, xValuesCentralityData, yValuesCentralityData, xErrorsCentralityData, yErrorsCentralityData);
      
      // Construct the corrected jet-hadron Vn centrality graphs from two files. The second one is assumed to have 5 % shifted centrality
      for(int iCentrality = 0; iCentrality < nCentralityBinsToCompare; iCentrality++){
        yPointer = flowGraphJetHadronCorrected[0][iAsymmetry][iCentrality][iFlow]->GetY();
        yValuesCentrality[iCentrality*2] = yPointer[ptBinForCentrality];
        yErrorsCentrality[iCentrality*2] = flowGraphJetHadronCorrected[0][iAsymmetry][iCentrality][iFlow]->GetErrorY(ptBinForCentrality);
        
        yPointer = flowGraphJetHadronCorrected[1][iAsymmetry][iCentrality][iFlow]->GetY();
        yValuesCentrality[iCentrality*2+1] = yPointer[ptBinForCentrality];
        yErrorsCentrality[iCentrality*2+1] = flowGraphJetHadronCorrected[1][iAsymmetry][iCentrality][iFlow]->GetErrorY(ptBinForCentrality);
        
        yPointer = flowGraphJetHadronCorrected[2][iAsymmetry][iCentrality][iFlow]->GetY();
        yValuesCentralityData[iCentrality] = yPointer[ptBinForCentrality];
        yErrorsCentralityData[iCentrality] = flowGraphJetHadronCorrected[2][iAsymmetry][iCentrality][iFlow]->GetErrorY(ptBinForCentrality);
      }
      
      flowGraphJetHadronCorrectedCentrality[0][iAsymmetry][iFlow] = new TGraphErrors(2*nCentralityBinsToCompare, xValuesCentrality, yValuesCentrality, xErrorsCentrality, yErrorsCentrality);
      flowGraphJetHadronCorrectedCentrality[1][iAsymmetry][iFlow] = new TGraphErrors(nCentralityBinsToCompare, xValuesCentralityData, yValuesCentralityData, xErrorsCentralityData, yErrorsCentralityData);
      
      // Construct the dihadron Vn centrality graphs from two files. The second one is assumed to have 5 % shifted centrality
      for(int iCentrality = 0; iCentrality < nCentralityBinsToCompare; iCentrality++){
        yPointer = flowGraphDihadron[0][iAsymmetry][iCentrality][iFlow]->GetY();
        yValuesCentrality[iCentrality*2] = yPointer[ptBinForCentrality];
        yErrorsCentrality[iCentrality*2] = flowGraphDihadron[0][iAsymmetry][iCentrality][iFlow]->GetErrorY(ptBinForCentrality);
        
        yPointer = flowGraphDihadron[1][iAsymmetry][iCentrality][iFlow]->GetY();
        yValuesCentrality[iCentrality*2+1] = yPointer[ptBinForCentrality];
        yErrorsCentrality[iCentrality*2+1] = flowGraphDihadron[1][iAsymmetry][iCentrality][iFlow]->GetErrorY(ptBinForCentrality);
        
        yPointer = flowGraphDihadron[2][iAsymmetry][iCentrality][iFlow]->GetY();
        yValuesCentralityData[iCentrality] = yPointer[ptBinForCentrality];
        yErrorsCentralityData[iCentrality] = flowGraphDihadron[2][iAsymmetry][iCentrality][iFlow]->GetErrorY(ptBinForCentrality);
      }
      
      flowGraphDihadronCentrality[0][iAsymmetry][iFlow] = new TGraphErrors(2*nCentralityBinsToCompare, xValuesCentrality, yValuesCentrality, xErrorsCentrality, yErrorsCentrality);
      flowGraphDihadronCentrality[1][iAsymmetry][iFlow] = new TGraphErrors(nCentralityBinsToCompare, xValuesCentralityData, yValuesCentralityData, xErrorsCentralityData, yErrorsCentralityData);
      
      // Construct the hadron vn centrality graphs from two files. The second one is assumed to have 5 % shifted centrality
      for(int iCentrality = 0; iCentrality < nCentralityBinsToCompare; iCentrality++){
        yPointer = flowGraphHadron[0][iAsymmetry][iCentrality][iFlow]->GetY();
        yValuesCentrality[iCentrality*2] = yPointer[ptBinForCentrality];
        yErrorsCentrality[iCentrality*2] = flowGraphHadron[0][iAsymmetry][iCentrality][iFlow]->GetErrorY(ptBinForCentrality);
        
        yPointer = flowGraphHadron[1][iAsymmetry][iCentrality][iFlow]->GetY();
        yValuesCentrality[iCentrality*2+1] = yPointer[ptBinForCentrality];
        yErrorsCentrality[iCentrality*2+1] = flowGraphHadron[1][iAsymmetry][iCentrality][iFlow]->GetErrorY(ptBinForCentrality);
        
        yPointer = flowGraphHadron[2][iAsymmetry][iCentrality][iFlow]->GetY();
        yValuesCentralityData[iCentrality] = yPointer[ptBinForCentrality];
        yErrorsCentralityData[iCentrality] = flowGraphHadron[2][iAsymmetry][iCentrality][iFlow]->GetErrorY(ptBinForCentrality);
      }
      
      flowGraphHadronCentrality[0][iAsymmetry][iFlow] = new TGraphErrors(2*nCentralityBinsToCompare, xValuesCentrality, yValuesCentrality, xErrorsCentrality, yErrorsCentrality);
      flowGraphHadronCentrality[1][iAsymmetry][iFlow] = new TGraphErrors(nCentralityBinsToCompare, xValuesCentralityData, yValuesCentralityData, xErrorsCentralityData, yErrorsCentralityData);
      
      // Construct the jet vn centrality graphs from two files. The second one is assumed to have 5 % shifted centrality
      for(int iCentrality = 0; iCentrality < nCentralityBinsToCompare; iCentrality++){
        yPointer = flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->GetY();
        yValuesCentrality[iCentrality*2] = yPointer[ptBinForCentrality];
        yErrorsCentrality[iCentrality*2] = flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->GetErrorY(ptBinForCentrality);
        
        yPointer = flowGraphJet[1][iAsymmetry][iCentrality][iFlow]->GetY();
        yValuesCentrality[iCentrality*2+1] = yPointer[ptBinForCentrality];
        yErrorsCentrality[iCentrality*2+1] = flowGraphJet[1][iAsymmetry][iCentrality][iFlow]->GetErrorY(ptBinForCentrality);
        
        yPointer = flowGraphJet[2][iAsymmetry][iCentrality][iFlow]->GetY();
        yValuesCentralityData[iCentrality] = yPointer[ptBinForCentrality];
        yErrorsCentralityData[iCentrality] = flowGraphJet[2][iAsymmetry][iCentrality][iFlow]->GetErrorY(ptBinForCentrality);
      }
      
      flowGraphJetCentrality[0][iAsymmetry][iFlow] = new TGraphErrors(2*nCentralityBinsToCompare, xValuesCentrality, yValuesCentrality, xErrorsCentrality, yErrorsCentrality);
      flowGraphJetCentrality[1][iAsymmetry][iFlow] = new TGraphErrors(nCentralityBinsToCompare, xValuesCentralityData, yValuesCentralityData, xErrorsCentralityData, yErrorsCentralityData);
      
      // Set the bin labels for x-axis
      for(int iCentrality = 0; iCentrality < 2*nCentralityBinsToCompare; iCentrality++){
        flowGraphJetHadronCentrality[0][iAsymmetry][iFlow]->GetXaxis()->ChangeLabel(iCentrality+1,330,-1,-1,-1,-1,binLabels[iCentrality]);
        flowGraphJetHadronCorrectedCentrality[0][iAsymmetry][iFlow]->GetXaxis()->ChangeLabel(iCentrality+1,330,-1,-1,-1,-1,binLabels[iCentrality]);
        flowGraphDihadronCentrality[0][iAsymmetry][iFlow]->GetXaxis()->ChangeLabel(iCentrality+1,330,-1,-1,-1,-1,binLabels[iCentrality]);
        flowGraphHadronCentrality[0][iAsymmetry][iFlow]->GetXaxis()->ChangeLabel(iCentrality+1,330,-1,-1,-1,-1,binLabels[iCentrality]);
        flowGraphJetCentrality[0][iAsymmetry][iFlow]->GetXaxis()->ChangeLabel(iCentrality+1,330,-1,-1,-1,-1,binLabels[iCentrality]);
      }
      
    } // Flow component loop
  } // Asymmetry loop
  
  // Setup the drawer for graphs
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceGraph();

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
//
//  // Draw graph with all the stages leading to final jet vn visible
//  if(drawGraphStages){
//
//    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
//      for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
//        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
//
//          sprintf(namerY,"V_{%d}",iFlow+1);
//          flowGraphJetHadron[0][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(kFullCircle);
//          flowGraphJetHadron[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(kBlack);
//          drawer->DrawGraph(flowGraphJetHadron[0][iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, 0, 0.35, "Track p_{T} (GeV)", namerY, " ", "p");
//
//          flowGraphJetHadronCorrected[0][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(kFullSquare);
//          flowGraphJetHadronCorrected[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(kGreen+3);
//          flowGraphJetHadronCorrected[0][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
//
//          flowGraphDihadron[0][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(kFullCross);
//          flowGraphDihadron[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(kRed);
//          flowGraphDihadron[0][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
//
//          flowGraphHadron[0][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(kFullDiamond);
//          flowGraphHadron[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(kMagenta);
//          flowGraphHadron[0][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
//
//          flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(kFullStar);
//          flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(kBlue);
//          flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
//
//          if(iFlow == 1){
//            atlasV2 = new TLine(0.5,atlasV2Number[iCentrality],maxTrackPt-0.5,atlasV2Number[iCentrality]);
//            atlasV2->SetLineStyle(2);
//            atlasV2->SetLineColor(kBlue);
//            atlasV2->Draw();
//          }
//
//          // Add a legend to the figure
//          legend = new TLegend(0.2,0.55,0.5,0.9);
//          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
//          legend->SetHeader(Form("Cent: %.0f-%.0f%%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
//          legend->AddEntry(flowGraphJetHadron[0][iAsymmetry][iCentrality][iFlow], Form("Jet-track V_{%d}", iFlow+1), "p");
//          legend->AddEntry(flowGraphJetHadronCorrected[0][iAsymmetry][iCentrality][iFlow], Form("Jet-track V_{%d}, corrected", iFlow+1), "p");
//          legend->AddEntry(flowGraphDihadron[0][iAsymmetry][iCentrality][iFlow], Form("Dihadron V_{%d}", iFlow+1), "p");
//          legend->AddEntry(flowGraphHadron[0][iAsymmetry][iCentrality][iFlow], Form("Hadron v_{%d}", iFlow+1), "p");
//          legend->AddEntry(flowGraphJet[0][iAsymmetry][iCentrality][iFlow], Form("Jet v_{%d}", iFlow+1), "p");
//
//          if(iFlow == 1){
//            legend->AddEntry(atlasV2, "ATLAS jet v_{2}", "l");
//          }
//
//          legend->Draw();
//
//          // Save the figures to file
//          if(saveFigures){
//              gPad->GetCanvas()->SaveAs(Form("figures/flowStagesV%d%s%s_C=%.0f-%.0f.pdf", iFlow+1, saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
//          } // Saving figures
//
//        } // Asymmetry loop
//      } // Flow component loop
//    } // Centrality loop
//
//  } // Draw stages leading to jet vn
//
//  // Draw all the different asymmetry bins in the same plot
//  if(drawGraphAsymmetryComparison){
//
//    // Draw the graphs as a function of pT
//    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
//      for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
//
//        // Setup a legend for the figure
//        legend = new TLegend(0.2,0.6,0.6,0.9);
//        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
//        legend->SetHeader(Form("PbPb C: %.0f-%.0f %%", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
//
//        // First, draw the systematic uncertainty bands to the canvas
//        if(drawSystematicUncertainties){
//          flowSystematicsJet[0][firstDrawnAsymmetryBin][iCentrality][iFlow]->SetFillColorAlpha(colors[firstDrawnAsymmetryBin],0.25);
//          drawer->DrawGraph(flowSystematicsJet[0][firstDrawnAsymmetryBin][iCentrality][iFlow], 0, maxTrackPt, -0.05, 0.35, "Track p_{T} (GeV)", Form("Jet v_{%d}", iFlow+1), " ", "2,same");
//
//          for(int iAsymmetry = firstDrawnAsymmetryBin+1; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
//            flowSystematicsJet[0][iAsymmetry][iCentrality][iFlow]->SetFillColorAlpha(colors[iAsymmetry],0.25);
//            flowSystematicsJet[0][iAsymmetry][iCentrality][iFlow]->Draw("2,same");
//          }
//        }
//
//        // After systematic uncertainties are drawn, draw the points on top of the bands
//        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
//          flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(colors[iAsymmetry]);
//          flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iAsymmetry]);
//          if(iAsymmetry == firstDrawnAsymmetryBin && !drawSystematicUncertainties){
//            drawer->DrawGraph(flowGraphJet[0][iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, -0.05, 0.35, "Track p_{T} (GeV)", Form("Jet v_{%d}", iFlow+1), " ", "p");
//          } else {
//            flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->Draw("psame");
//          }
//          legend->AddEntry(flowGraphJet[0][iAsymmetry][iCentrality][iFlow],asymmetryLegend[iAsymmetry],"p");
//
//          zeroLine->Draw();
//        }
//
//        legend->Draw();
//
//        // Save the figures to file
//        if(saveFigures){
//          gPad->GetCanvas()->SaveAs(Form("figures/flowAsymmetryComparison%s%s_C=%.0f-%.0f.pdf", Form("V%d", iFlow+1), saveComment.Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
//
//        } // Saving figures
//
//      } // Asymmetry loop
//    } // Centrality loop
//  } // Draw flow component comparison graphs
//
//  // Draw all the different Vn components in the same plot
//  if(drawGraphVnComparison){
//
//    // Draw the graphs as a function of pT
//    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
//      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
//
//        // Setup a legend for the figure
//        legend = new TLegend(0.2,0.6,0.6,0.9);
//        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
//        legend->SetHeader(Form("PbPb C: %.0f-%.0f %%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
//
//        // First, draw the systematic uncertainty bands to the canvas
//        if(drawSystematicUncertainties){
//          flowSystematicsJet[0][iAsymmetry][iCentrality][firstDrawnVn-1]->SetFillColorAlpha(flowColors[firstDrawnVn-1],0.25);
//          drawer->DrawGraph(flowSystematicsJet[0][iAsymmetry][iCentrality][firstDrawnVn-1], 0, maxTrackPt, -0.05, 0.35, "Track p_{T} (GeV)", "v_{n}", " ", "2,same");
//
//          for(int iFlow = firstDrawnVn; iFlow <= lastDrawnVn-1; iFlow++){
//            flowSystematicsJet[0][iAsymmetry][iCentrality][iFlow]->SetFillColorAlpha(flowColors[iFlow],0.25);
//            flowSystematicsJet[0][iAsymmetry][iCentrality][iFlow]->Draw("2,same");
//          }
//        }
//
//        // After systematic uncertainties are drawn, draw the points on top of the bands
//        for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
//          flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(flowColors[iFlow]);
//          flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iFlow]);
//          if(iFlow == firstDrawnVn-1 && !drawSystematicUncertainties){
//            drawer->DrawGraph(flowGraphJet[0][iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, -0.05, 0.35, "Track p_{T} (GeV)", "v_{n}", " ", "p");
//          } else {
//            flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->Draw("psame");
//          }
//          legend->AddEntry(flowGraphJet[0][iAsymmetry][iCentrality][iFlow],Form("Jet v_{%d}",iFlow+1),"p");
//
//          zeroLine->Draw();
//        }
//
//        legend->Draw();
//
//        // Save the figures to file
//        if(saveFigures){
//          gPad->GetCanvas()->SaveAs(Form("figures/flowVnComparison%s%s_C=%.0f-%.0f.pdf", saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
//
//        } // Saving figures
//
//      } // Asymmetry loop
//    } // Centrality loop
//  } // Draw flow component comparison graphs
//
  // Compare graphs from different files
  if(drawCentralityComparison){
    
    drawer->SetNDivisionsX(510);
    drawer->SetBottomMargin(0.18);
    drawer->SetTitleOffsetX(1.63);
    drawer->SetLabelOffsetX(0.04);
    
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
        
        if(drawJetHadronVnCentralityComparison){
          
          sprintf(namerY,"Jet-hadron V_{%d}",iFlow+1);
          legend = new TLegend(0.2,0.65,0.5,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->SetHeader(Form("Track p_{T}: %.1f-%.1f GeV%s", trackPtBinBorders[ptBinForCentrality], trackPtBinBorders[ptBinForCentrality+1], asymmetryString[iAsymmetry].Data()));
          
          
          flowGraphJetHadronCentrality[0][iAsymmetry][iFlow]->SetMarkerStyle(kFullSquare);
          flowGraphJetHadronCentrality[0][iAsymmetry][iFlow]->SetMarkerColor(kBlack);
          flowGraphJetHadronCentrality[0][iAsymmetry][iFlow]->SetMarkerSize(1.3);
          drawer->DrawGraphCustomAxes(flowGraphJetHadronCentrality[0][iAsymmetry][iFlow], 0, 2*nCentralityBinsToCompare+1, 0, 0.08, "Centrality", namerY, " ", "ap");
          
          flowGraphJetHadronCentrality[1][iAsymmetry][iFlow]->SetMarkerStyle(kFullCircle);
          flowGraphJetHadronCentrality[1][iAsymmetry][iFlow]->SetMarkerColor(kRed);
          flowGraphJetHadronCentrality[1][iAsymmetry][iFlow]->SetMarkerSize(1.3);
          flowGraphJetHadronCentrality[1][iAsymmetry][iFlow]->Draw("p,same");
          
          legend->AddEntry(flowGraphJetHadronCentrality[0][iAsymmetry][iFlow], "MC, RecoGen", "p");
          legend->AddEntry(flowGraphJetHadronCentrality[1][iAsymmetry][iFlow], "PbPb data", "p");
          
          legend->Draw();
          
          // Save the figures to file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/jetHadronV%dCentralityComparison%s%s_T=%.0f-%.0f.pdf", iFlow+1, saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), trackPtBinBorders[ptBinForCentrality], trackPtBinBorders[ptBinForCentrality+1]));
          }
          
        } // Centrality comparison for jet-hadron Vn
        
        if(drawDihadronVnCentralityComparison){
          
          sprintf(namerY,"Dihadron V_{%d}",iFlow+1);
          legend = new TLegend(0.2,0.65,0.5,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->SetHeader(Form("Track p_{T}: %.1f-%.1f GeV%s", trackPtBinBorders[ptBinForCentrality], trackPtBinBorders[ptBinForCentrality+1], asymmetryString[iAsymmetry].Data()));
          
          flowGraphDihadronCentrality[0][iAsymmetry][iFlow]->SetMarkerStyle(kFullSquare);
          flowGraphDihadronCentrality[0][iAsymmetry][iFlow]->SetMarkerColor(kBlack);
          flowGraphDihadronCentrality[0][iAsymmetry][iFlow]->SetMarkerSize(1.3);
          drawer->DrawGraphCustomAxes(flowGraphDihadronCentrality[0][iAsymmetry][iFlow], 0, 2*nCentralityBinsToCompare+1, 0, 0.08, "Centrality", namerY, " ", "ap");
          
          flowGraphDihadronCentrality[1][iAsymmetry][iFlow]->SetMarkerStyle(kFullCircle);
          flowGraphDihadronCentrality[1][iAsymmetry][iFlow]->SetMarkerColor(kRed);
          flowGraphDihadronCentrality[1][iAsymmetry][iFlow]->SetMarkerSize(1.3);
          flowGraphDihadronCentrality[1][iAsymmetry][iFlow]->Draw("p,same");
          
          legend->AddEntry(flowGraphDihadronCentrality[0][iAsymmetry][iFlow], "MC, gen particles", "p");
          legend->AddEntry(flowGraphDihadronCentrality[1][iAsymmetry][iFlow], "PbPb data", "p");
          
          legend->Draw();
          
          // Save the figures to file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/dihadronV%dCentralityComparison%s%s_T=%.0f-%.0f.pdf", iFlow+1, saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), trackPtBinBorders[ptBinForCentrality], trackPtBinBorders[ptBinForCentrality+1]));
          }
          
        } // Centrality comparison for dihadron Vn
        
        if(drawHadronVnCentralityComparison){
          
          sprintf(namerY,"Hadron v_{%d}",iFlow+1);
          legend = new TLegend(0.2,0.65,0.5,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->SetHeader(Form("Track p_{T}: %.1f-%.1f GeV%s", trackPtBinBorders[ptBinForCentrality], trackPtBinBorders[ptBinForCentrality+1], asymmetryString[iAsymmetry].Data()));
          
          flowGraphHadronCentrality[0][iAsymmetry][iFlow]->SetMarkerStyle(kFullSquare);
          flowGraphHadronCentrality[0][iAsymmetry][iFlow]->SetMarkerColor(kBlack);
          flowGraphHadronCentrality[0][iAsymmetry][iFlow]->SetMarkerSize(1.3);
          drawer->DrawGraphCustomAxes(flowGraphHadronCentrality[0][iAsymmetry][iFlow], 0, 2*nCentralityBinsToCompare+1, 0, 0.3, "Centrality", namerY, " ", "ap");
          
          flowGraphHadronCentrality[1][iAsymmetry][iFlow]->SetMarkerStyle(kFullCircle);
          flowGraphHadronCentrality[1][iAsymmetry][iFlow]->SetMarkerColor(kRed);
          flowGraphHadronCentrality[1][iAsymmetry][iFlow]->SetMarkerSize(1.3);
          flowGraphHadronCentrality[1][iAsymmetry][iFlow]->Draw("p,same");
          
          legend->AddEntry(flowGraphHadronCentrality[0][iAsymmetry][iFlow], "MC, gen particles", "p");
          legend->AddEntry(flowGraphHadronCentrality[1][iAsymmetry][iFlow], "PbPb data", "p");
          
          legend->Draw();
          
          // Save the figures to file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/hadronV%dCentralityComparison%s%s_T=%.0f-%.0f.pdf", iFlow+1, saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), trackPtBinBorders[ptBinForCentrality], trackPtBinBorders[ptBinForCentrality+1]));
          }
          
        } // Centrality comparison for hadron vn
        
        if(drawJetVnCentralityComparison){
          
          sprintf(namerY,"Jet v_{%d}",iFlow+1);
          legend = new TLegend(0.2,0.25,0.5,0.5);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->SetHeader(Form("Track p_{T}: %.1f-%.1f GeV%s", trackPtBinBorders[ptBinForCentrality], trackPtBinBorders[ptBinForCentrality+1], asymmetryString[iAsymmetry].Data()));
          
          flowGraphJetCentrality[0][iAsymmetry][iFlow]->SetMarkerStyle(kFullSquare);
          flowGraphJetCentrality[0][iAsymmetry][iFlow]->SetMarkerColor(kBlack);
          flowGraphJetCentrality[0][iAsymmetry][iFlow]->SetMarkerSize(1.3);
          drawer->DrawGraphCustomAxes(flowGraphJetCentrality[0][iAsymmetry][iFlow], 0, 2*nCentralityBinsToCompare+1, 0, 0.3, "Centrality", namerY, " ", "ap");
          
          flowGraphJetCentrality[1][iAsymmetry][iFlow]->SetMarkerStyle(kFullCircle);
          flowGraphJetCentrality[1][iAsymmetry][iFlow]->SetMarkerColor(kRed);
          flowGraphJetCentrality[1][iAsymmetry][iFlow]->SetMarkerSize(1.3);
          flowGraphJetCentrality[1][iAsymmetry][iFlow]->Draw("p,same");
          
          legend->AddEntry(flowGraphJetCentrality[0][iAsymmetry][iFlow], "MC, reco jets", "p");
          legend->AddEntry(flowGraphJetCentrality[1][iAsymmetry][iFlow], "PbPb data", "p");
          
          //          if(iFlow == 1){
          //            atlasV2 = new TLine(0.5,atlasV2Number[iCentrality],maxTrackPt-0.5,atlasV2Number[iCentrality]);
          //            atlasV2->SetLineStyle(2);
          //            atlasV2->SetLineColor(kBlue);
          //            atlasV2->Draw();
          //            legend->AddEntry(atlasV2, "ATLAS jet v_{2}", "l");
          //          }
          
          legend->Draw();
          
          // Save the figures to file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/jetV%dCentralityComparison%s%s_T=%.0f-%.0f.pdf", iFlow+1, saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), trackPtBinBorders[ptBinForCentrality], trackPtBinBorders[ptBinForCentrality+1]));
          }
          
        } // Centrality comparison for jet vn
        
      } // Asymmetry loop
    } // Flow component loop
    
  } // Draw stages leading to jet vn

}

