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
  TString graphFileName = "flowGraphs_PbPb2018_correctedJetHadron_correctedDihadron_2021-01-08.root";
  // flowGraphs_PbPb2018_caloJets_correctedJetHadron_correctedDihadron_xjBins_2020-12-03.root
  // flowGraphs_PbPb2018_systematicUncertainties_regularJetHadron_2020-12-04.root
  TFile *graphFile[maxFiles];
  graphFile[0] = TFile::Open(directoryName+graphFileName);
  
  // Other files whose results can be compared with the nominal file
  int nComparisonFiles = 1;
  TString comparisonFileName[] = {  "flowGraphs_PbPbMC2018_caloJets_correctedJetHadron_correctedEventDihadron_bestQvectorForEachCentrality_2021-01-14.root", "flowGraphs_PbPbMC2018_caloJets_correctedJetHadron_correctedDihadron_5pCentShift_qVectorBelow2p2OnlyForDihadron_2021-01-19.root", "flowGraphs_PbPbMC2018_caloJets_correctedJetHadron_correctedDihadron_5pCentShift_qVectorBelow2p8OnlyForDihadron_2021-01-19.root", "flowGraphs_PbPbMC2018_caloJets_correctedJetHadron_correctedDihadron_5pCentShift_qVectorBelow3p3OnlyForDihadron_2021-01-19.root", "flowGraphs_PbPbMC2018_5pCentShift_qVectorCutBelow5ForDihadron_correctedJetHadron_correctedDihadron_2021-01-11.root", "flowGraphs_PbPbMC2018_5pCentShift_qVectorCutBelow6ForDihadron_correctedJetHadron_correctedDihadron_2021-01-11.root", "flowGraphs_PbPbMC2018_5pCentShift_qVectorCutBelow7ForDihadron_correctedJetHadron_correctedDihadron_2021-01-11.root", "flowGraphs_PbPbMC2018_caloJets_5pCentShift_correctedJetHadron_sameEventDihadron_2020-11-18.root", "flowGraphs_PbPb2018_systematicUncertainties_backgroundAdjustedJetHadron_2020-12-04.root", "flowGraphs_PbPb2018_systematicUncertainties_nearEtaJetHadron_2020-12-04.root", "flowGraphs_PbPb2018_systematicUncertainties_farEtaJetHadron_2020-12-04.root", "flowGraphs_PbPbMC2018_caloJets_improvisedMixingJetHadron_sameEventDihadron_2020-11-13.root", "flowGraphs_PbPbMC2018_caloJets_5pCentShift_correctedJetHadron_sameEventDihadron_2020-11-18.root", "flowGraphs_PbPb2018_caloJets_improvisedMixingJetHadron_correctedDihadron_noJetCorrections_2020-11-05.root", "flowGraphs_PbPbMC2018_caloJets_improvisedMixingJetHadron_sameEventDihadron_2020-11-05.root", "flowGraphs_PbPb2018_caloJets_improvisedMixingJetHadron_correctedDihadron_noJetCorrections_2020-11-05.root", "qVectorStudy_manualCut6_recoJets_sameEventJetHadron_sameEventDihadron_2020-10-20.root", "qVectorStudy_manualCut7_recoJets_sameEventJetHadron_sameEventDihadron_2020-10-20.root", "qVectorStudy_noCut_correctedJetHadron_correctedDihadron.root", "flowGraphs_PbPbData_noJetReconstructionCorrection_fullDihadronStats.root", "finalGraphTestNew.root", ""};
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    graphFile[iFile+1] = TFile::Open(directoryName+comparisonFileName[iFile]);
  }
  
  // Legend text given to each compared file
  TString fileLegend[] = {"Calo jets", "This analysis", "MC+5%, Q < 2.2", "MC+5%, Q < 2.8",  "MC+5%, Q < 3.3", "MC+5%, Q < 3.33"};
  
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
  const bool drawAtlasV2 = false;                     // Draw a line showing the v2 result from ATLAS
  const bool fitJetVn = true;                         // Fit a constant line to jet vn points
  const bool doSummaryCorrection = true;              // Correct the jet vn summary plots based on fits
  
  // Plots to be compared between files
  const bool drawJetHadronVnFileComparison = false;
  const bool drawDihadronVnFileComparison = false;
  const bool drawHadronVnFileComparison = false;
  const bool drawJetVnFileComparison = true;
  const bool drawJetHadronYieldFileComparison = false;
  const bool drawDihadronYieldFileComparison = false;
  const bool drawFileComparison = drawJetHadronVnFileComparison || drawDihadronVnFileComparison || drawHadronVnFileComparison || drawJetVnFileComparison || drawJetHadronYieldFileComparison || drawDihadronYieldFileComparison;
  
  const bool drawSystematicUncertainties = false;     // Include systematic uncertainties in the plots
  
  const bool saveFigures = true;                     // Save the figures in a file
  TString saveComment = "_caloJetCurrent";              // String to be added to saved file names
  
  int firstDrawnAsymmetryBin = nAsymmetryBins;
  int lastDrawnAsymmetryBin = nAsymmetryBins;
  
  int firstDrawnVn = 2;
  int lastDrawnVn = 3;
  
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
  
  TGraphErrors *yieldGraphJetHadron[maxFiles][nAsymmetryBins+1][nCentralityBins];
  TGraphErrors *yieldGraphDihadron[maxFiles][nAsymmetryBins+1][nCentralityBins];
  
  // Jet vn summary graph in all centrality bins
  TGraphErrors *flowSummaryJet[maxFiles][nAsymmetryBins+1][nFlowComponents];
  TGraphErrors *atlasJetV2graph;
  TGraphErrors *cmsHighPtV2;
  double summaryXaxis[nCentralityBins];
  double summaryXaxisError[nCentralityBins];
  double summaryYaxis[maxFiles][nAsymmetryBins+1][nFlowComponents][nCentralityBins];
  double summaryYaxisError[maxFiles][nAsymmetryBins+1][nFlowComponents][nCentralityBins];
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    summaryXaxis[iCentrality] = iCentrality+1;
    summaryXaxisError[iCentrality] = 0;
    for(int iFile = 0; iFile < maxFiles; iFile++){
      for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins+1; iAsymmetry++){
        for(int iFlow = 0; iFlow < nFlowComponents; iFlow++){
          summaryYaxis[iFile][iAsymmetry][iFlow][iCentrality] = 0;
          summaryYaxisError[iFile][iAsymmetry][iFlow][iCentrality] = 0;
        }
      }
    }
  }
  TString binLabels[] = {"0-10%"," ","10-30%"," ","30-50%"," ","50-90%"};
  
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
        
        // If drawn, load the yield graphs
        if(drawJetHadronYieldFileComparison || drawDihadronYieldFileComparison){
          sprintf(histogramNamer,"longRangeJetHadronYield/longRangeJetHadronYield_A%dC%d", iAsymmetry, iCentrality);
          yieldGraphJetHadron[iFile][iAsymmetry][iCentrality] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          sprintf(histogramNamer,"longRangeDihadronYield/longRangeDihadronYield_A%dC%d", iAsymmetry, iCentrality);
          yieldGraphDihadron[iFile][iAsymmetry][iCentrality] = (TGraphErrors*) graphFile[iFile]->Get(histogramNamer);
          
          // Remove the points that are above the maximum pT limit from the graphs
          xAxisPoints = yieldGraphJetHadron[iFile][iAsymmetry][iCentrality]->GetX();
          nPoints = yieldGraphJetHadron[iFile][iAsymmetry][iCentrality]->GetN();
          for(int iX = nPoints-1; iX >= 0; iX--){
            if(xAxisPoints[iX] > maxTrackPt){
              yieldGraphJetHadron[iFile][iAsymmetry][iCentrality]->RemovePoint(iX);
              yieldGraphDihadron[iFile][iAsymmetry][iCentrality]->RemovePoint(iX);
            } else {
              break;
            }
          }
          
        }
        
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
  
  // How to zoom vn plots    //    v1,   v2,   v3,   v4
  double jetHadronZoomTable[4] = {0.1, 0.05, 0.005, 0.01};
  double dihadronZoomTable[4]  = {0.08, 0.08, 0.03, 0.03};

  // Helper variables for drawing figures
  TLegend *legend;
  TLegend *vLegend;
  int markers[] = {kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross, kOpenStar};
  int fullMarkers[] = {kFullSquare, kFullCircle, kFullDiamond, kFullCross, kFullStar};
  int secondMarkers[] = {kFullCircle, kFullCross, kFullSquare, kFullCircle, kFullFourTrianglesPlus};
  int colors[] = {kBlue,kRed,kGreen+2,kBlack, kMagenta};
  int flowColors[] = {kBlue, kBlack, kRed, kGreen+3, kMagenta};
  int fileColors[] = {kBlack, kBlue, kRed, kGreen+3, kMagenta, kCyan};
  TString asymmetryString[] = {" 0.0 < x_{j} < 0.6", " 0.6 < x_{j} < 0.8", " 0.8 < x_{j} < 1.0", ""};
  TString asymmetryLegend[] = {"0.0 < x_{j} < 0.6", "0.6 < x_{j} < 0.8", "0.8 < x_{j} < 1.0", "x_{j} integrated"};
  TString compactAsymmetryString[] = {"_A=0v0-0v6", "_A=0v6-0v8", "_A=0v8-1v0", ""};

  // Numbers from HP2020 conference presentation
  TLine *atlasV2;
  double atlasV2Number[] = {0.018, 0.03, 0.035, 0.03};
  
  // Numbers from fitting the data from paper arXiv:1702.00630 at large pT
  double cmsHighPtV2Number[] = {0.0220, 0.0376, 0.0431, 0.04};
  double cmsHighPtV2Error[] = {0.0019, 0.0016, 0.0027, 0.04};
  
  TString legendString;
  char namerY[100];
  
  TLine *zeroLine = new TLine(0,0,maxTrackPt,0);
  zeroLine->SetLineStyle(2);
  
  TLine *shortZeroLine = new TLine(0.8,0,3.2,0);
  shortZeroLine->SetLineStyle(2);
  
  TLine *vnLine = new TLine(0,0,maxTrackPt,0);
  vnLine->SetLineStyle(2);
  
  TH1D *vnLineError = new TH1D("vnLineError","vnLineError",1,0,maxTrackPt);
  
  double vnValue, vnError;
  
  //double minZoomJetVn[] = {0.085, 0.11, 0.075, 0.05};  // Nominal
  //double maxZoomJetVn[] = {0.135, 0.16, 0.125, 0.15};  // Nominal
  
  double minZoomJetVn[] = {0.085, 0.05, 0.03, 0};
  double maxZoomJetVn[] = {0.25, 0.2, 0.2, 0.15};
  
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
          
          if(iFlow == 1 && drawAtlasV2){
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
          
          if(iFlow == 1 && drawAtlasV2){
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
        
        // Asymmetry comparison for jet-hadron Vn
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
          flowGraphJetHadron[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(colors[iAsymmetry]);
          flowGraphJetHadron[0][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iAsymmetry]);
          if(iAsymmetry == firstDrawnAsymmetryBin && !drawSystematicUncertainties){
            drawer->DrawGraph(flowGraphJetHadron[0][iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, 0, 0.05, "Track p_{T} (GeV)", Form("Jet-hadron V_{%d}", iFlow+1), " ", "p");
          } else {
            flowGraphJetHadron[0][iAsymmetry][iCentrality][iFlow]->Draw("psame");
          }
          //legend->AddEntry(flowGraphJetHadron[0][iAsymmetry][iCentrality][iFlow],asymmetryLegend[iAsymmetry],"p");
          
          zeroLine->Draw();
        }
        
        legend->Draw();
        
        // Save the figures to file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/jetHadronAsymmetryComparison%s%s_C=%.0f-%.0f.pdf", Form("V%d", iFlow+1), saveComment.Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
        }
        
        // Asymmetry comparison for dihadron Vn
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
          flowGraphDihadron[0][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(colors[iAsymmetry]);
          flowGraphDihadron[0][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iAsymmetry]);
          if(iAsymmetry == firstDrawnAsymmetryBin && !drawSystematicUncertainties){
            drawer->DrawGraph(flowGraphDihadron[0][iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, 0, 0.1, "Track p_{T} (GeV)", Form("Dihadron V_{%d}", iFlow+1), " ", "p");
          } else {
            flowGraphDihadron[0][iAsymmetry][iCentrality][iFlow]->Draw("psame");
          }
          //legend->AddEntry(flowGraphDihadron[0][iAsymmetry][iCentrality][iFlow],asymmetryLegend[iAsymmetry],"p");
          
          zeroLine->Draw();
        }
        
        legend->Draw();
        
        // Save the figures to file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/dihadronAsymmetryComparison%s%s_C=%.0f-%.0f.pdf", Form("V%d", iFlow+1), saveComment.Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
        }
        
      } // Flow component loop
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
                 drawer->DrawGraph(flowGraphJetHadron[iFile][iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, 0, jetHadronZoomTable[iFlow], "Track p_{T} (GeV)", namerY, " ", "p");
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
                 drawer->DrawGraph(flowGraphDihadron[iFile][iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, 0, dihadronZoomTable[iFlow], "Track p_{T} (GeV)", namerY, " ", "p");
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
            
          } // File comparison for hadron vn
          
          if(drawJetVnFileComparison){
          
            sprintf(namerY,"Jet v_{%d}",iFlow+1);
            legend = new TLegend(0.2,0.55,0.5,0.9);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->SetHeader(Form("Cent: %.0f-%.0f%%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
            
            for(int iFile = 0; iFile < nComparisonFiles+1; iFile++){
              flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iFile]);
              flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(fileColors[iFile]);
              flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerSize(1.3);
              if(iFile == 0){
                 drawer->DrawGraph(flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, minZoomJetVn[iCentrality], maxZoomJetVn[iCentrality], "Track p_{T} (GeV)", namerY, " ", "p");
              } else {
                flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
              }
              
              // Fit a constant line to the jet v_{n} values
              if(fitJetVn){
                if(iFile == 1){
                  flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->Fit("pol0","0","",0,3);
                } else {
                  flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->Fit("pol0","0");
                }
                
                vnValue = flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->GetFunction("pol0")->GetParameter(0);
                vnError = flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->GetFunction("pol0")->GetParError(0);
                
                vnLineError->SetBinContent(1, vnValue);
                vnLineError->SetBinError(1, vnError);
                vnLineError->SetFillColorAlpha(fileColors[iFile], 0.25);
                vnLineError->DrawCopy("e2,same");
                
                vnLine->SetLineColor(fileColors[iFile]);
                vnLine->DrawLine(0, vnValue, maxTrackPt, vnValue);
                
                summaryYaxis[iFile][iAsymmetry][iFlow][iCentrality] = vnValue;
                summaryYaxisError[iFile][iAsymmetry][iFlow][iCentrality] = vnError;
              }
              
              legend->AddEntry(flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow], fileLegend[iFile], "p");
              
            }
            
            if(iFlow == 1 && drawAtlasV2){
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
    
    // File comparison for yields
    if(drawJetHadronYieldFileComparison || drawDihadronYieldFileComparison){
      
      double maxJetHadronYield[] = {3000, 2000, 1000, 500};
      double maxDihadronYield[] = {5500000, 2000000, 400000, 100000};
      drawer->SetTitleOffsetY(1.6);
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
          
          if(drawJetHadronYieldFileComparison){
          
            legend = new TLegend(0.5,0.55,0.8,0.9);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->SetHeader(Form("Cent: %.0f-%.0f%%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
            
            for(int iFile = 0; iFile < nComparisonFiles+1; iFile++){
              yieldGraphJetHadron[iFile][iAsymmetry][iCentrality]->SetMarkerStyle(fullMarkers[iFile]);
              yieldGraphJetHadron[iFile][iAsymmetry][iCentrality]->SetMarkerColor(fileColors[iFile]);
              yieldGraphJetHadron[iFile][iAsymmetry][iCentrality]->SetMarkerSize(1.3);
              if(iFile == 0){
                 drawer->DrawGraph(yieldGraphJetHadron[iFile][iAsymmetry][iCentrality], 0, maxTrackPt, 0, maxJetHadronYield[iCentrality], "Track p_{T} (GeV)", "Jet-hadron yield", " ", "p");
              } else {
                yieldGraphJetHadron[iFile][iAsymmetry][iCentrality]->Draw("p,same");
              }
              
              legend->AddEntry(yieldGraphJetHadron[iFile][iAsymmetry][iCentrality], fileLegend[iFile], "p");
              
            }
            
            legend->Draw();
            
            // Save the figures to file
            if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/jetHadronYieldComparison%s%s_C=%.0f-%.0f.pdf", saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
            }
            
          } // File comparison for jet-hadron yields
          
          if(drawDihadronYieldFileComparison){
          
            legend = new TLegend(0.5,0.55,0.8,0.9);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->SetHeader(Form("Cent: %.0f-%.0f%%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
            
            for(int iFile = 0; iFile < nComparisonFiles+1; iFile++){
              yieldGraphDihadron[iFile][iAsymmetry][iCentrality]->SetMarkerStyle(fullMarkers[iFile]);
              yieldGraphDihadron[iFile][iAsymmetry][iCentrality]->SetMarkerColor(fileColors[iFile]);
              yieldGraphDihadron[iFile][iAsymmetry][iCentrality]->SetMarkerSize(1.3);
              if(iFile == 0){
                 drawer->DrawGraph(yieldGraphDihadron[iFile][iAsymmetry][iCentrality], 0, maxTrackPt, 0, maxDihadronYield[iCentrality], "Track p_{T} (GeV)", "Dihadron yield", " ", "p");
              } else {
                yieldGraphDihadron[iFile][iAsymmetry][iCentrality]->Draw("p,same");
              }
              
              legend->AddEntry(yieldGraphDihadron[iFile][iAsymmetry][iCentrality], fileLegend[iFile], "p");
              
            }
            
            legend->Draw();
            
            // Save the figures to file
            if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/dihadronYieldComparison%s%s_C=%.0f-%.0f.pdf", saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
            }
            
          } // File comparison for dihadron yields
          
        } // Asymmetry loop
      } // Centrality loop
    } // File comparison for yields
    
    drawer->SetTitleOffsetY(1.1);
    
  } // Draw file comparison plots
  
  // If we draw jet vn file comparison and fit the vn plots, we can also draw summary plots
  if(drawJetVnFileComparison && fitJetVn){
    
    // For summary correction it is assumed that first the data is given and then MC.
    if(doSummaryCorrection){
      for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
          for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
            for(int iFile = 1; iFile < nComparisonFiles+1; iFile++){
              summaryYaxis[iFile][iAsymmetry][iFlow][iCentrality] = summaryYaxis[0][iAsymmetry][iFlow][iCentrality] - summaryYaxis[iFile][iAsymmetry][iFlow][iCentrality];
              summaryYaxisError[iFile][iAsymmetry][iFlow][iCentrality] = TMath::Sqrt(TMath::Power(summaryYaxisError[0][iAsymmetry][iFlow][iCentrality],2) + TMath::Power(summaryYaxisError[iFile][iAsymmetry][iFlow][iCentrality],2));
            }
          } // Centrality loop
        } // Asymmetry loop
      } // Flow component loop
      nComparisonFiles++;
    }
    
    drawer->SetNDivisionsX(510);
    drawer->SetBottomMargin(0.18);
    drawer->SetTitleOffsetX(1.63);
    drawer->SetLabelOffsetX(0.04);
    drawer->SetTitleOffsetY(1.6);
    
    // First, we need to construct the graphs based on the fit values
    atlasJetV2graph = new TGraphErrors(nCentralityBins, summaryXaxis, atlasV2Number, summaryXaxisError, summaryXaxisError);
    atlasJetV2graph->SetMarkerStyle(kFullDiamond);
    atlasJetV2graph->SetMarkerColor(kViolet-2);
    atlasJetV2graph->SetMarkerSize(1.3);
    
    cmsHighPtV2 = new TGraphErrors(nCentralityBins, summaryXaxis, cmsHighPtV2Number, summaryXaxisError, cmsHighPtV2Error);
    cmsHighPtV2->SetMarkerStyle(kFullStar);
    cmsHighPtV2->SetMarkerColor(kAzure+9);
    cmsHighPtV2->SetMarkerSize(1.3);
    
    for(int iFile = 0; iFile < nComparisonFiles+1; iFile++){
      for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
          flowSummaryJet[iFile][iAsymmetry][iFlow] = new TGraphErrors(nCentralityBins, summaryXaxis, summaryYaxis[iFile][iAsymmetry][iFlow], summaryXaxisError, summaryYaxisError[iFile][iAsymmetry][iFlow]);
          flowSummaryJet[iFile][iAsymmetry][iFlow]->SetMarkerStyle(secondMarkers[iFile]);
          flowSummaryJet[iFile][iAsymmetry][iFlow]->SetMarkerColor(fileColors[iFile]);
          flowSummaryJet[iFile][iAsymmetry][iFlow]->SetMarkerSize(1.3);
          
          // Set the bin labels for x-axis
          for(int iCentrality = 0; iCentrality < nCentralityBins*2; iCentrality++){
            flowSummaryJet[iFile][iAsymmetry][iFlow]->GetXaxis()->ChangeLabel(iCentrality+1,-1,-1,-1,-1,-1,binLabels[iCentrality]);
          } // Centrality loop
        } // Asymmetry loop
      } // Flow component loop
    } // File loop
    
    // Once the graphs are constructed, they can be plotted
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
        
        legend = new TLegend(0.2,0.7,0.5,0.9); //0.2,0.6,0.5,0.9
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        
        for(int iFile = 1; iFile < nComparisonFiles; iFile++){
          if(iFile == 1){
            sprintf(namerY,"Jet v_{%d}",iFlow+1);
            drawer->DrawGraphCustomAxes(flowSummaryJet[iFile][iAsymmetry][iFlow], 0, 4, -0.02, 0.1, "Centrality", namerY, " ", "ap"); // 0, 4, -0.05, 0.3
          } else {
            flowSummaryJet[iFile][iAsymmetry][iFlow]->Draw("p,same");
          }
          legend->AddEntry(flowSummaryJet[iFile][iAsymmetry][iFlow], fileLegend[iFile], "p");
        } // File loop
        
        if(iFlow == 1 && iAsymmetry == nAsymmetryBins){
          atlasJetV2graph->Draw("p,same");
          legend->AddEntry(atlasJetV2graph, "ATLAS v_{2}", "p"); // (#scale[0.8]{HP 2020})
          
          cmsHighPtV2->Draw("p,same");
          legend->AddEntry(cmsHighPtV2, "CMS high p_{T} v_{2}", "p"); // (#scale[0.8]{PLB 776 (2018) 195})
        }
        
        shortZeroLine->Draw();
        legend->Draw();
        
        // Save the figures to file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/jetV%dSummary%s%s.pdf", iFlow+1, saveComment.Data(), compactAsymmetryString[iAsymmetry].Data()));
        }
      } // Asymmetry loop
    } // Flow component loop
    
    // For testing porposes, just draw something
    //drawer->DrawGraph(flowSummaryJet[0][nAsymmetryBins][1], 0, 4, 0, 0.3, "Centrality", "Jet v2", " ", "p");
    //drawer->DrawGraphCustomAxes(flowSummaryJet[0][nAsymmetryBins][1], 0, 4, 0, 0.3, "Centrality", "Jet v2", " ", "ap");

  }
  
}

