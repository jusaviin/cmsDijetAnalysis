#include "LongRangeSystematicOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"
#include "xCanvas.h"

/*
 * Macro for plotting final jet v2 and v3 results with systematic uncertainties.
 */
void finalLongRangePlotter(){

  // ============= //
  // Configuration //
  // ============= //
  
  // Input file name for data
  TString directoryName = "flowGraphs/";
  TString inputFileName = "summaryPlot_multiplicityMatch_caloJets_integratedBinsUpTo3_nominalCorrection_jetEta1v3_2022-03-04.root";
  // summaryPlot_multiplicityMatch_caloJets_integratedBinsUpTo3_nominalCorrection_jetEta1v3_2022-03-04.root
  // summaryPlot_multiplicityMatch_caloJets_updatedScaling_jetEta1v3_2022-02-28.root
  // summaryPlot_multiplicityMatch_jetEta1v3_2022-02-13.root
  // summaryPlot_akCaloJet_matchHadronV2ScaleYieldWith4pCentShift_2021-06-03.root
  TString uncertaintyFileName = "systematicUncertainties_multiMatchNominal_scalingUpdates_2022-03-04.root";
  // systematicUncertainties_addMultiplicityMatch_updateCentralValue_2022-02-14.root
  // systematicUncertainties_allSources_finalCorrection_2021-08-10.root
  
  // Text to be put into legend for the input graphs
  const char* legendText = "This analysis";
  
  // Define the bins that are drawn
  const int nCentralityBins = 3;  // Number of drawn centrality bins
  
  const int maxVn = 4;            // Maximum defined vn. Plots are made upto v4.
  const int firstDrawnVn = 2;     // First drawn flow component
  const int lastDrawnVn = 4;      // Last drawn flow component
  
  // Define if previous results should be included in the plot
  const bool drawAtlasJetV2 = false;
  const bool drawCmsHigtPtV2 = true;
  
  // Drawing for preliminary tag
  const bool drawBigCanvas = true;
  const bool drawPreliminaryTag = false;
  
  // Save the final plots
  const bool saveFigures = true;
  TString saveComment = "_initialStyle";
  
  // =========== //
  // Read graphs //
  // =========== //
  
  TGraphErrors *jetVnGraph[maxVn];
  TGraphErrors *jetVnUncertainty[maxVn];
  
  // Initialize the graphs to NULL
  for(int iFlow = 0; iFlow < maxVn; iFlow++){
    jetVnGraph[iFlow] = NULL;
    jetVnUncertainty[iFlow] = NULL;
  }
  
  // Open input file for reading and read graphs from it
  TFile *inputFile = TFile::Open(directoryName+inputFileName);
  for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
    jetVnGraph[iFlow] = (TGraphErrors*) inputFile->Get(Form("summaryV%d", iFlow+1));
    
    // If the graph we wanted to load does not exist, inform the user and end program
    if(jetVnGraph[iFlow] == NULL){
      cout << "Hey dude! The file: " << directoryName.Data() << inputFileName.Data() << " does not contain graph: " << Form("summaryV%d", iFlow+1) <<  "." << endl;
      cout << "Cannot do the plotting, mate! Be a good lad and make sure the graph is there next time, ok?" << endl;
      return;
    }
  }
  
  // Read the systematic uncertainty graphs
  TFile *uncertaintyFile = TFile::Open(directoryName+uncertaintyFileName);
  LongRangeSystematicOrganizer *uncertaintyOrganizer = new LongRangeSystematicOrganizer(uncertaintyFile);
  uncertaintyOrganizer->AdjustCentralPoints(jetVnGraph);
  
  for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
    jetVnUncertainty[iFlow] = uncertaintyOrganizer->GetLongRangeSystematicUncertainty(iFlow);
    
    // If the graph we wanted to load does not exist, inform the user and end program
    if(jetVnUncertainty[iFlow] == NULL){
      cout << "Hey dude! The file: " << directoryName.Data() << uncertaintyFileName.Data() << " does not contain uncertainties for " << Form("v%d", iFlow+1) <<  "." << endl;
      cout << "No uncertainties can be plotted." << endl;
    }
    
  }
  
  // Set the bin labels for x-axis
  TString binLabels[] = {"0-10%"," ","10-30%"," ","30-50%"," ","50-90%"};
  for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
    for(int iCentrality = 0; iCentrality < jetVnUncertainty[iFlow]->GetN()*2; iCentrality++){
      jetVnUncertainty[iFlow]->GetXaxis()->ChangeLabel(iCentrality+1,-1,-1,-1,-1,-1,binLabels[iCentrality]);
    } // Centrality loop
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
  
  // Draw the graphs for selected flow components
  TLegend *legend;
  double errorY;
  
  double minZoom[] = {0,-0.045,-0.045,-0.045};
  double maxZoom[] = {0.1,0.085,0.085,0.085};
  
  TLine *zeroLine = new TLine(0.75,0,3.5,0);
  zeroLine->SetLineStyle(2);
  
  // =======================================================================
  // == Drawing style, where each histogram is drawn into it's own canvas ==
  // =======================================================================
  
  for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
    legend = new TLegend(0.2,0.7,0.5,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    
    // Set the style for the jet vn values
    jetVnGraph[iFlow]->SetMarkerStyle(kFullCircle);
    jetVnGraph[iFlow]->SetMarkerColor(kBlue);
    
    if(jetVnUncertainty[iFlow] != NULL){
      
      // Set the style for uncertainties
      for(int iCentrality; iCentrality < jetVnUncertainty[iFlow]->GetN(); iCentrality++){
        errorY = jetVnUncertainty[iFlow]->GetErrorY(iCentrality);
        jetVnUncertainty[iFlow]->SetPointError(iCentrality, 0.1, errorY);
      }
      jetVnUncertainty[iFlow]->SetFillColorAlpha(kBlue, 0.3);
      drawer->DrawGraphCustomAxes(jetVnUncertainty[iFlow], 0, 4, minZoom[iFlow], maxZoom[iFlow], "Centrality", Form("Jet v_{%d}", iFlow+1), " ", "a,e2");
      
      jetVnGraph[iFlow]->Draw("p,same");
    } else {
      drawer->DrawGraphCustomAxes(jetVnGraph[iFlow], 0, 4, minZoom[iFlow], maxZoom[iFlow], "Centrality", Form("Jet v_{%d}", iFlow+1), " ", "ap");
    }
    
    zeroLine->Draw();
    
    legend->AddEntry(jetVnGraph[iFlow], legendText, "p");
    
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
      gPad->GetCanvas()->SaveAs(Form("figures/finalJetV%d%s.pdf", iFlow+1, saveComment.Data()));
    }
  }
  
  // =============================================
  // ===== Drawing style with one big canvas =====
  // =============================================
  
  if(drawBigCanvas){
    
    // Draw all the distributions to big canvases
    auxi_canvas *bigCanvas;
    TLatex *mainTitle;
    
    // Draw a big canvas and put all the plots in it
    bigCanvas = new auxi_canvas("bigCanvas", "", 1500, 550);
    bigCanvas->SetMargin(0.07, 0.01, 0.16, 0.01);
    bigCanvas->divide(1,3);
    
    mainTitle = new TLatex();
    
    int bigCanvasColor[] = {kBlack, kBlue, kRed, kGreen+3};
    int bigCanvasMarker[] = {kFullSquare, kFullCircle, kFullDiamond, kFullCross};
    
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      bigCanvas->CD(iFlow);
      
      // Adjust the number of divisions
      jetVnUncertainty[iFlow]->GetYaxis()->SetNdivisions(510);
      jetVnUncertainty[iFlow]->GetYaxis()->SetLabelOffset(0.01);
      jetVnUncertainty[iFlow]->GetYaxis()->SetTitle("Dijet v_{n}");
      jetVnUncertainty[iFlow]->GetYaxis()->SetTitleOffset(1.5);
      jetVnUncertainty[iFlow]->GetYaxis()->SetTitleSize(0.06);
      
      // Adjust label sizes for all but the leftmost plot
      if(iFlow > firstDrawnVn-1){
        jetVnUncertainty[iFlow]->GetXaxis()->SetTitleSize(0.069);
        jetVnUncertainty[iFlow]->GetXaxis()->SetTitleOffset(1.09);
        jetVnUncertainty[iFlow]->GetXaxis()->SetLabelSize(0.058);
        jetVnUncertainty[iFlow]->GetXaxis()->SetLabelOffset(0.0135);
      } else {
        jetVnUncertainty[iFlow]->GetXaxis()->SetTitleOffset(1.25);
        jetVnUncertainty[iFlow]->GetXaxis()->SetLabelOffset(0.02);
      }
      
      // Set consistent marker style with data
      jetVnUncertainty[iFlow]->SetMarkerSize(2);
      jetVnUncertainty[iFlow]->SetMarkerStyle(bigCanvasMarker[iFlow]);
      jetVnUncertainty[iFlow]->SetLineColor(bigCanvasColor[iFlow]);
      jetVnUncertainty[iFlow]->SetFillColorAlpha(bigCanvasColor[iFlow], 0.3);
      jetVnUncertainty[iFlow]->SetMarkerColor(bigCanvasColor[iFlow]);
      jetVnGraph[iFlow]->SetMarkerStyle(bigCanvasMarker[iFlow]);
      jetVnGraph[iFlow]->SetMarkerColor(bigCanvasColor[iFlow]);
      jetVnGraph[iFlow]->SetMarkerSize(2);
      
      jetVnUncertainty[iFlow]->Draw("a,e2");
      jetVnGraph[iFlow]->Draw("same,p");
      
      // Create a legend for the plot division
      if(iFlow > 1) {
        legend = new TLegend(0.06,0.89,0.7,0.98);
      } else {
        legend = new TLegend(0.23,0.8,0.73,0.98);
      }
      
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      if(iFlow > 1) legend->SetTextSize(0.06);  // Need to increase text size for smaller divisions
      
      legend->AddEntry(jetVnUncertainty[iFlow], Form("Dijet v_{%d}", iFlow+1), "pf");
      
      if(iFlow == 1 && drawAtlasJetV2){
        atlasJetV2graph->SetMarkerSize(2);
        atlasJetV2graph->Draw("p,same");
        legend->AddEntry(atlasJetV2graph, "ATLAS v_{2}", "p"); // (#scale[0.8]{HP 2020})
      }
      
      if(iFlow == 1 && drawCmsHigtPtV2){
        cmsHighPtV2->SetMarkerSize(2);
        cmsHighPtV2->Draw("p,same");
        legend->AddEntry(cmsHighPtV2, "CMS hadron v_{2} , p_{T} > 20 GeV", "p"); // (#scale[0.8]{PLB 776 (2018) 195})
      }
      
      legend->Draw();
      
      zeroLine->Draw();
      
    } // Flow component loop
    
    // Draw all necessary CMS text to the plot
    bigCanvas->cd(0);
    
    mainTitle->SetTextFont(62);
    mainTitle->SetTextSize(0.065);
    double cmsYPosition = 0.89;
    if(drawPreliminaryTag) cmsYPosition = 0.91;
    mainTitle->DrawLatexNDC(0.9, cmsYPosition, "CMS");
    
    mainTitle->SetTextFont(42);
    mainTitle->SetTextSize(0.055);
    if(drawPreliminaryTag)mainTitle->DrawLatexNDC(0.878, 0.86, "Preliminary");
    mainTitle->DrawLatexNDC(0.72, 0.76, "PbPb #sqrt{s_{NN}} = 5.02 TeV, 1.7 nb^{-1}");
    mainTitle->DrawLatexNDC(0.54, 0.92, "anti-k_{T} R = 0.4");
    mainTitle->DrawLatexNDC(0.54, 0.84, "|#eta_{jet}| < 1.3");
    mainTitle->DrawLatexNDC(0.54, 0.76, "p_{T,1} > 120 GeV");
    mainTitle->DrawLatexNDC(0.54, 0.68, "p_{T,2} > 50 GeV");
    mainTitle->DrawLatexNDC(0.54, 0.60, "#Delta#varphi_{1,2} > #frac{5#pi}{6}");
    mainTitle->DrawLatexNDC(0.13, 0.33, "Factorization region:");
    mainTitle->DrawLatexNDC(0.13, 0.25, "0.7 < Hadron p_{T} < 3 GeV");
    
    // Save the figures to file
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/finalBigCanvasJetVn%s.pdf", saveComment.Data()));
    }
    
  }
}
