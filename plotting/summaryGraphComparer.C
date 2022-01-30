#include "LongRangeSystematicOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

/*
 * Macro for plotting final jet v2 and v3 results with systematic uncertainties.
 */
void summaryGraphComparer(){

  // ============= //
  // Configuration //
  // ============= //
  
  // Input file name for data
  TString directoryName = "flowGraphs/";
  const int nInputFiles = 3;
  TString inputFileName[] =  {"summaryPlot_multiplicityScheme_2022-01-26.root", "summaryPlot_shiftToMatchV2Check_2022-01-27.root", "summaryPlot_akCaloJet_matchHadronV2ScaleYieldWith4pCentShift_2021-06-03.root",  "summaryPlot_multiplicitySchemeVnScale_2022-01-27.root", "summaryPlot_jetEventPlaneCheck_2022-01-27.root", "summaryPlot_akCaloJet_qVectorWeight_2021-08-10.root",  "summaryPlot_akCaloJet_smearedJER_lowStatDihadron_2021-08-16.root", "summaryPlot_akCaloJet_dihadronDeltaEta2to3v5_2021-08-06.root", "summaryPlot_akCaloJet_dihadronDeltaEta2v5to4_2021-08-06.root", "summaryPlot_akCaloJet_correctionWith25pMoreQuarkJets_2021-07-26.root",  "summaryPlot_akCaloJet_noTrackEfficiency_2021-07-14.root",  "summaryPlot_akCaloJet_matchHadronV2ScaleYieldWith3pCentShift_2021-06-03.root", "summaryPlot_akCaloJet_matchHadronV2ScaleYieldWith5pCentShift_2021-06-03.root", "summaryPlot_akCaloJet_matchHadronV2ScaleYieldWith4v5pCentShift_2021-06-03.root", "summaryPlot_akCaloJet_averageCorrectionWith3pCentShift_2021-06-02.root", "summaryPlot_akCaloJet_averageCorrection_2021-06-01.root",  "summaryPlot_akCaloJet_averageCorrectionWith4v5pCentShift_2021-06-02.root",  "summaryPlot_akCaloJet_averageCorrectionWith5pCentShift_2021-06-02.root", "summaryPlot_akCaloJet_matchHadronV2CorrectionWith3pCentShift_2021-06-02.root", "summaryPlot_akCaloJet_matchHadronV2CorrectionWith4pCentShift_2021-06-02.root", "summaryPlot_akCaloJet_matchHadronV2CorrectionWith4v5pCentShift_2021-06-02.root", "summaryPlot_akCaloJet_matchHadronV2CorrectionWith5pCentShift_2021-06-02.root"};
  
  TString uncertaintyFileName = "systematicUncertainties_allSources_finalCorrection_2021-08-10.root";
  
  // Text to be put into legend for the input graphs
  TString legendText[] = {"Multiplicity match", "1.5% shift + scale", "Nominal (Q-cut)", "Correction centrality shift 4.5 %", "Extrapolate to hadron v_{2} , 3 %", "Extrapolate to hadron v_{2} , 4 %", "Extrapolate to hadron v_{2} , 4.5 %", "Extrapolate to hadron v_{2} , 5 %"};
  
  // Define the bins that are drawn
  const int nCentralityBins = 3;  // Number of drawn centrality bins
  
  const int maxVn = 4;            // Maximum defined vn. Plots are made upto v4.
  const int firstDrawnVn = 3;     // First drawn flow component
  const int lastDrawnVn = 3;      // Last drawn flow component
  
  // Define if previous results should be included in the plot
  const bool drawAtlasJetV2 = false;
  const bool drawCmsHigtPtV2 = false;
  const bool drawUncertaintyBand = false;
  const bool drawRatio = true;
  
  // Save the final plots
  const bool saveFigures = false;
  TString saveComment = "_scaleCheck";
  TString figureFormat = "png";
  
  // =========== //
  // Read graphs //
  // =========== //
  
  TGraphErrors *jetVnGraph[nInputFiles][maxVn];
  TGraphErrors *jetVnRatioGraph[nInputFiles][maxVn];
  TGraphErrors *jetVnUncertainty[maxVn];
  
  // Initialize the graphs to NULL
  for(int iFile = 0; iFile < nInputFiles; iFile++){
    for(int iFlow = 0; iFlow < maxVn; iFlow++){
      jetVnGraph[iFile][iFlow] = NULL;
    }
  }
  
  // Open input file for reading and read graphs from it
  TFile *inputFile[nInputFiles];
  
  for(int iFile = 0; iFile < nInputFiles; iFile++){
    inputFile[iFile] = TFile::Open(directoryName+inputFileName[iFile]);
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      jetVnGraph[iFile][iFlow] = (TGraphErrors*) inputFile[iFile]->Get(Form("summaryV%d", iFlow+1));
      
      // If the graph we wanted to load does not exist, inform the user and end program
      if(jetVnGraph[iFile][iFlow] == NULL){
        cout << "Hey dude! The file: " << directoryName.Data() << inputFileName[iFile].Data() << " does not contain graph: " << Form("summaryV%d", iFlow+1) <<  "." << endl;
        cout << "Cannot do the plotting, mate! Be a good lad and make sure the graph is there next time, ok?" << endl;
        return;
      }
      
      // For the ratio graphs, clone the Vn graphs
      jetVnRatioGraph[iFile][iFlow] = (TGraphErrors*) jetVnGraph[iFile][iFlow]->Clone(Form("ratio%d%d",iFile,iFlow));
      
    } // Flow component loop
  } // File loop
  
  // Calculate the ratio to the first file
  double xPoint1, yPoint1, xPoint2, yPoint2, yError1, yError2, ratioValue, combinedError;
  TString binLabels[] = {"0-10%"," ","10-30%"," ","30-50%"," ","50-90%"};
  
  for(int iFile = 0; iFile < nInputFiles; iFile++){
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      
      for(int iPoint = 0; iPoint < jetVnGraph[0][iFlow]->GetN(); iPoint++){
        
        jetVnGraph[0][iFlow]->GetPoint(iPoint, xPoint1, yPoint1);
        yError1 = jetVnGraph[0][iFlow]->GetErrorY(iPoint);
        jetVnGraph[iFile][iFlow]->GetPoint(iPoint, xPoint2, yPoint2);
        yError2 = jetVnGraph[iFile][iFlow]->GetErrorY(iPoint);
        
        ratioValue = yPoint2 / yPoint1;
        combinedError = TMath::Sqrt(TMath::Power(yError2/yPoint1,2) + TMath::Power((yPoint2*yError1)/(yPoint1*yPoint1),2));
        jetVnRatioGraph[iFile][iFlow]->SetPoint(iPoint,xPoint1,ratioValue);
        jetVnRatioGraph[iFile][iFlow]->SetPointError(iPoint, 0, combinedError);
        
      } // Point loop
      
      // Set the bin labels for x-axis
      for(int iPoint = 0; iPoint < jetVnGraph[iFile][iFlow]->GetN()*2; iPoint++){
        jetVnRatioGraph[iFile][iFlow]->GetXaxis()->ChangeLabel(iPoint+1,-1,-1,-1,-1,-1,binLabels[iPoint]);
      }
    } // Flow component loop
  } // Comparison file loop
  
  // Read the systematic uncertainty graphs
  TFile *uncertaintyFile = TFile::Open(directoryName+uncertaintyFileName);
  LongRangeSystematicOrganizer *uncertaintyOrganizer = new LongRangeSystematicOrganizer(uncertaintyFile);
  
  for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
    jetVnUncertainty[iFlow] = uncertaintyOrganizer->GetLongRangeSystematicUncertainty(iFlow);
    
    // If the graph we wanted to load does not exist, inform the user and end program
    if(jetVnUncertainty[iFlow] == NULL){
      cout << "Hey dude! The file: " << directoryName.Data() << uncertaintyFileName.Data() << " does not contain uncertainties for " << Form("v%d", iFlow+1) <<  "." << endl;
      cout << "No uncertainties can be plotted." << endl;
    }
  }
  
  // ===================== //
  // Previous measurements //
  // ===================== //
  
  // Previous results that can be plotted together with data from this analysis
  double summaryXaxis[nCentralityBins];
  double summaryXaxisError[nCentralityBins];
  double summaryYaxisError[nCentralityBins];
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    summaryXaxis[iCentrality] = iCentrality+1;
    summaryXaxisError[iCentrality] = 0;
    summaryYaxisError[iCentrality] = 0;
  }
  
  // ATLAS results taken from those presented in Hard Probes 2020: ATLAS-CONF-2020-019 (https://cds.cern.ch/record/2720249?ln=en)
  const double atlasV2Number[] = {0.018, 0.03, 0.035, 0.03};
  TGraphErrors* atlasJetV2graph = new TGraphErrors(nCentralityBins, summaryXaxis, atlasV2Number, summaryXaxisError, summaryXaxisError);
  atlasJetV2graph->SetMarkerStyle(kFullDiamond);
  atlasJetV2graph->SetMarkerColor(kViolet-2);
  atlasJetV2graph->SetMarkerSize(1.8);
  
  // CMS high pT results extracted from analysis arXiv:1702.00630 (PLB 776 (2018) 195)
  const double cmsHighPtV2Number[] = {0.0220, 0.0376, 0.0431, 0.04};
  const double cmsHighPtV2Error[] = {0.0019, 0.0016, 0.0027, 0.04};
  TGraphErrors* cmsHighPtV2 = new TGraphErrors(nCentralityBins, summaryXaxis, cmsHighPtV2Number, summaryXaxisError, cmsHighPtV2Error);
  cmsHighPtV2->SetMarkerStyle(kFullStar);
  cmsHighPtV2->SetMarkerColor(kAzure+9);
  cmsHighPtV2->SetMarkerSize(1.8);
  
  // ============== //
  // Draw the plots //
  // ============== //
  
  // Setup the drawer for graphs
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceGraph();
  drawer->SetNDivisionsX(510);
  drawer->SetBottomMargin(0.18);
  drawer->SetTitleOffsetX(1.63);
  drawer->SetLabelOffsetX(0.04);
  drawer->SetTitleOffsetY(1.6);
  
  double minZoom[] = {0,0,-0.03,-0.2};
  double maxZoom[] = {0.12,0.12,0.03,0.2};
  double errorY;
  
  // Draw the graphs for selected flow components
  TLegend *legend;
  const int markers[] = {kFullSquare, 89, kFullCross, kFullFourTrianglesPlus, kFullDoubleDiamond, kFullDiamond, kFullStar};
  const int colors[] = {kBlue, kRed, kGreen+3, kMagenta, kCyan, kBlack, kViolet};
  
  const int markers2[] = {kFullSquare, kFullDiamond, kFullDoubleDiamond, kFullCross, kFullFourTrianglesPlus, kFullStar};
  const int colors2[] = {kBlack, kBlue, kRed, kGreen+3, kMagenta, kCyan, kBlack, kViolet};
  
  for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
    legend = new TLegend(0.2,0.9 - 0.05*nInputFiles - 0.05*drawAtlasJetV2 - 0.05*drawCmsHigtPtV2,0.5,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    
    // Set the style for the jet vn values
    for(int iFile = 0; iFile < nInputFiles; iFile++){
      jetVnGraph[iFile][iFlow]->SetMarkerStyle(markers2[iFile]);
      jetVnGraph[iFile][iFlow]->SetMarkerColor(colors[iFile]);
      //jetVnGraph[iFile][iFlow]->SetLineWidth(0);
      legend->AddEntry(jetVnGraph[iFile][iFlow], legendText[iFile], "p");
    }
    
    // Draw the first graph to the canvas
    if(drawUncertaintyBand){
      
      // Set the style for uncertainties
      for(int iCentrality; iCentrality < jetVnUncertainty[iFlow]->GetN(); iCentrality++){
        errorY = jetVnUncertainty[iFlow]->GetErrorY(iCentrality);
        jetVnUncertainty[iFlow]->SetPointError(iCentrality, 0.1, errorY);
      }
      jetVnUncertainty[iFlow]->SetFillColorAlpha(kBlue, 0.3);
      drawer->DrawGraphCustomAxes(jetVnUncertainty[iFlow], 0, 4, minZoom[iFlow], maxZoom[iFlow], "Centrality", Form("Jet v_{%d}", iFlow+1), " ", "a,e2");
      
      jetVnGraph[0][iFlow]->Draw("p,same");
    } else {
      drawer->DrawGraphCustomAxes(jetVnGraph[0][iFlow], 0, 4, minZoom[iFlow], maxZoom[iFlow], "Centrality", Form("Jet v_{%d}", iFlow+1), " ", "ap");
    }
    
    // Draw the rest of the graphs to the same canvas
    for(int iFile = 1; iFile < nInputFiles; iFile++){
      jetVnGraph[iFile][iFlow]->Draw("p,same");
    }
    
    
    if(iFlow == 1 && drawAtlasJetV2){
      atlasJetV2graph->Draw("p,same");
      legend->AddEntry(atlasJetV2graph, "ATLAS v_{2}", "p"); // (#scale[0.8]{HP 2020})
    }
    
    if(iFlow == 1 && drawCmsHigtPtV2){
      cmsHighPtV2->Draw("p,same");
      legend->AddEntry(cmsHighPtV2, "CMS high p_{T} v_{2}", "p"); // (#scale[0.8]{PLB 776 (2018) 195})
    }
    
    legend->Draw();
    
    // Save the figures to file
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/jetV%dComparison%s.%s", iFlow+1, saveComment.Data(), figureFormat.Data()));
    }
  }
  
  // Draw ratios of jet vn plots
  if(drawRatio){
    
    TLine *oneLine = new TLine(0.8,1,3.2,1);
    oneLine->SetLineStyle(2);
    
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      
      legend = new TLegend(0.2,0.9 - 0.05*nInputFiles,0.5,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      
      // Set the style for the jet vn values
      for(int iFile = 1; iFile < nInputFiles; iFile++){
        jetVnRatioGraph[iFile][iFlow]->SetMarkerStyle(markers2[iFile]);
        jetVnRatioGraph[iFile][iFlow]->SetMarkerColor(colors2[iFile]);
        //jetVnRatioGraph[iFile][iFlow]->SetLineWidth(0);
        legend->AddEntry(jetVnRatioGraph[iFile][iFlow], Form("%s / %s", legendText[iFile].Data(), legendText[0].Data()), "p");
      }
      
      drawer->DrawGraphCustomAxes(jetVnRatioGraph[1][iFlow], 0, 4, 0.75, 1.25, "Centrality", Form("Jet v_{%d} ratio", iFlow+1), " ", "ap");
      
      // Draw the rest of the graphs to the same canvas
      for(int iFile = 2; iFile < nInputFiles; iFile++){
        jetVnRatioGraph[iFile][iFlow]->Draw("p,same");
      }
      
      legend->Draw();
      oneLine->Draw();
      
      // Save the figures to file
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/jetV%dRatioComparison%s.%s", iFlow+1, saveComment.Data(), figureFormat.Data()));
      }
      
    }
  }
}
