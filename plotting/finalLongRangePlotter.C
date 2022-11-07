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
  // systematicUncertainties_multiMatchNominal_addPtBinVariation_2022-04-21.root
  // systematicUncertainties_multiMatchNominal_scalingUpdates_2022-03-04.root
  
  // Define the bins that are drawn
  const int nCentralityBins = 3;  // Number of drawn centrality bins
  
  const int maxVn = 4;            // Maximum defined vn. Plots are made upto v4.
  const int firstDrawnVn = 2;     // First drawn flow component
  const int lastDrawnVn = 4;      // Last drawn flow component
  
  // Define if previous results should be included in the plot
  const bool drawAtlasJetV2 = false;
  const bool drawCmsHigtPtV2 = false;
  
  // Drawing for preliminary tag
  bool drawIndividualGraphs = true;
  const bool drawSmallCanvas = false;
  const bool drawBigCanvas = false;
  const bool drawPreliminaryTag = false;
  const bool leadSubTitles = false; // true: Use lead and sub instead of 1 and 2 in the figure legend
  
  // The correct style is set while drawing individual canvases. We must include that when doing big canvas
  if(drawBigCanvas) drawIndividualGraphs = true;
  
  // Printing of values to console
  const bool printCentralValues = false;  // Print the central vn values
  
  // Save the final plots
  const bool saveFigures = true;
  TString saveComment = "_supplementaryUpdates";
  const char* figureFormat = "pdf";
  
  // Marker colors and styles
  int bigCanvasColor[] = {kBlack, kBlue, kRed, kGreen+3};
  int bigCanvasMarker[] = {kFullSquare, kFullCircle, kFullDiamond, kFullCross};
  
  // Zooming for y-axis
  double minZoom[] = {0,-0.045,-0.045,-0.045};
  double maxZoom[] = {0.1,0.085,0.085,0.085};
  
  // Save the final results for HepData
  bool saveGraphsForHepData = false;
  
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
    
    // Set the style for the jet vn values
    jetVnGraph[iFlow]->SetMarkerStyle(bigCanvasMarker[iFlow]);
    jetVnGraph[iFlow]->SetMarkerColor(bigCanvasColor[iFlow]);
    jetVnGraph[iFlow]->SetMarkerSize(2);
    
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
    
    // Set the style for the jet vn uncertainties
    jetVnUncertainty[iFlow]->SetMarkerStyle(bigCanvasMarker[iFlow]);
    jetVnUncertainty[iFlow]->SetLineColor(bigCanvasColor[iFlow]);
    jetVnUncertainty[iFlow]->SetFillColorAlpha(bigCanvasColor[iFlow], 0.3);
    jetVnUncertainty[iFlow]->SetMarkerColor(bigCanvasColor[iFlow]);
    jetVnUncertainty[iFlow]->SetMarkerSize(2);
    
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
    summaryXaxis[iCentrality] = iCentrality+0.5;
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
  cmsHighPtV2->SetMarkerColor(kMagenta);
  cmsHighPtV2->SetMarkerSize(2.8); // 2.8
  
  // ============== //
  // Draw the plots //
  // ============== //
  
  // Setup the drawer for graphs
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceGraph();
  drawer->SetRelativeCanvasSize(2,0.9);
  drawer->SetNDivisionsX(5);
  drawer->SetNDivisionsY(510);
  drawer->SetTopMargin(0.16);
  drawer->SetBottomMargin(0.16);
  drawer->SetTitleOffsetX(1.3);
  drawer->SetLabelOffsetX(10);
  drawer->SetTitleOffsetY(1.77);
  drawer->SetLabelOffsetY(0.01);
  
  // Draw the graphs for selected flow components
  TLegend *legend;
  TLegend *anotherLegend;
  double errorY;
  double cmsYPosition, cmsXPosition;
  double currentX;
  
  TLine *zeroLine = new TLine(0.25,0,2.85,0);
  zeroLine->SetLineStyle(2);
  
  TLatex *preliminaryText = new TLatex();
  
  // =======================================================================
  // == Drawing style, where each histogram is drawn into it's own canvas ==
  // =======================================================================
  
  if(drawIndividualGraphs){
    
    // Positioning of the selection information
    //                                 v1     v2     v3     v4
    double ptChargedPositionX[] =    {0.23,  0.23,  0.23,  0.23};
    double ptChargedPositionY[] =    {0.34,  0.34,  0.77,  0.77};
    double jetRadiusPositionX[] =    {0.23,  0.23,  0.23,  0.23};
    double jetRadiusPositionY[] =    {0.275, 0.275, 0.705, 0.705};
    double etaJetPositionX[] =       {0.23,  0.23,  0.23,  0.23};
    double etaJetPositionY[] =       {0.21,  0.21,  0.64,  0.64};
    double ptLeadingPositionX[] =    {0.58,  0.58,  0.58,  0.58};
    double ptLeadingPositionY[] =    {0.34,  0.34,  0.77,  0.77};
    double ptSubleadingPositionX[] = {0.58,  0.58,  0.58,  0.58};
    double ptSubleadingPositionY[] = {0.275, 0.275, 0.705, 0.705};
    double deltaPhiPositionX[] =     {0.58,  0.58,  0.58,  0.58};
    double deltaPhiPositionY[] =     {0.21,  0.21,  0.64,  0.64};
    
    // Positioning of the legend
    //                    v1     v2     v3     v4
    double legendX1[] = {0.2,   0.2,   0.2,   0.2};
    double legendX2[] = {0.7,   0.7,   0.7,   0.7};
    double legendY1[] = {0.73,  0.73,  0.53,  0.21};
    double legendY2[] = {0.80,  0.80,  0.60,  0.28};
    
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      
      // Displace the points a bit in x-direction for nicer drawing result
      for(int iPoint = 0; iPoint < jetVnUncertainty[iFlow]->GetN(); iPoint++){
        currentX = jetVnUncertainty[iFlow]->GetPointX(iPoint);
        jetVnUncertainty[iFlow]->SetPointX(iPoint, currentX-0.5);
        jetVnGraph[iFlow]->SetPointX(iPoint, currentX-0.5);
      }
      
      // Create legends for the plot
      if(iFlow == 1 && drawCmsHigtPtV2){
        legend = new TLegend(0.68,0.75,1.1,0.82);
        cmsHighPtV2->SetMarkerSize(4);
      } else {
        legend = new TLegend(legendX1[iFlow],legendY1[iFlow],legendX2[iFlow],legendY2[iFlow]);
      }
      anotherLegend = new TLegend(0.19,0.67,0.4,0.81);
      
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.045);legend->SetTextFont(62);
      anotherLegend->SetFillStyle(0);anotherLegend->SetBorderSize(0);anotherLegend->SetTextSize(0.04);anotherLegend->SetTextFont(62);
      
      if(drawCmsHigtPtV2) legend->SetTextSize(0.04);
      
      if(jetVnUncertainty[iFlow] != NULL){
        
        // Set the style for uncertainties
        for(int iCentrality = 0; iCentrality < jetVnUncertainty[iFlow]->GetN(); iCentrality++){
          errorY = jetVnUncertainty[iFlow]->GetErrorY(iCentrality);
          jetVnUncertainty[iFlow]->SetPointError(iCentrality, 0.1, errorY);
        }
        
        //jetVnUncertainty[iFlow]->GetYaxis()->SetNdivisions(510);
        jetVnUncertainty[iFlow]->SetMarkerSize(4);
        drawer->DrawGraphCustomAxes(jetVnUncertainty[iFlow], 0, 4, minZoom[iFlow], maxZoom[iFlow], "Centrality", Form("Dijet v_{%d}", iFlow+1), " ", "a,e2");
        
        jetVnGraph[iFlow]->SetMarkerSize(4);
        jetVnGraph[iFlow]->Draw("p,same");
        
        legend->AddEntry(jetVnUncertainty[iFlow], Form("Dijet v_{%d}", iFlow+1), "pf");
        
      } else {
        drawer->DrawGraphCustomAxes(jetVnGraph[iFlow], 0, 4, minZoom[iFlow], maxZoom[iFlow], "Centrality", "Dijet v_{n}", " ", "ap");
        legend->AddEntry(jetVnGraph[iFlow], Form("Dijet v_{%d}", iFlow+1), "p");
      }
      
      zeroLine->Draw();
      
      if(iFlow == 1 && drawAtlasJetV2){
        atlasJetV2graph->Draw("p,same");
        legend->AddEntry(atlasJetV2graph, "ATLAS v_{2}", "p"); // (#scale[0.8]{HP 2020})
      }
      
      if(iFlow == 1 && drawCmsHigtPtV2){
        cmsHighPtV2->Draw("p,same");
        //legend->AddEntry(cmsHighPtV2, "CMS high p_{T} v_{2}", "p"); // (#scale[0.8]{PLB 776 (2018) 195})
        anotherLegend->AddEntry(cmsHighPtV2, "CMS charged hadron v_{2}", "p"); // (#scale[0.8]{PLB 776 (2018) 195})
        anotherLegend->AddEntry((TObject*)0, "p_{T} > 20 GeV, |#eta| < 1", "");
        anotherLegend->AddEntry((TObject*)0, "PLB 776 (2018) 195", "");
        anotherLegend->Draw();
      }
      
      legend->Draw();
      
      // Draw CMS supplementary tag
      preliminaryText->SetTextFont(62);
      preliminaryText->SetTextSize(0.06);
      preliminaryText->DrawLatexNDC(0.08, 0.941, "CMS");
      
      preliminaryText->SetTextFont(52);
      preliminaryText->SetTextSize(0.055);
      preliminaryText->DrawLatexNDC(0.215, 0.941, "Supplementary");
      
      preliminaryText->SetTextFont(42);
      preliminaryText->SetTextSize(0.05);
      preliminaryText->DrawLatexNDC(0.625, 0.941, "arXiv:2210.08325");
      
      // Draw luminosity and selection information to the figure
      preliminaryText->DrawLatexNDC(0.23, 0.87, "PbPb #sqrt{s_{NN}} = 5.02 TeV, 1.69 nb^{-1}");
      
      preliminaryText->SetTextSize(0.04);
      preliminaryText->DrawLatexNDC(ptChargedPositionX[iFlow], ptChargedPositionY[iFlow], "0.7 < p^{ch}_{T} < 3 GeV");
      preliminaryText->DrawLatexNDC(jetRadiusPositionX[iFlow], jetRadiusPositionY[iFlow], "anti-k_{T} R = 0.4");
      preliminaryText->DrawLatexNDC(etaJetPositionX[iFlow], etaJetPositionY[iFlow], "|#eta_{jet}| < 1.3");
      preliminaryText->DrawLatexNDC(ptLeadingPositionX[iFlow], ptLeadingPositionY[iFlow], "p_{T}^{lead} > 120 GeV");
      preliminaryText->DrawLatexNDC(ptSubleadingPositionX[iFlow], ptSubleadingPositionY[iFlow], "p_{T}^{sub} > 50 GeV");
      preliminaryText->DrawLatexNDC(deltaPhiPositionX[iFlow], deltaPhiPositionY[iFlow], "#Delta#varphi > #frac{5#pi}{6}");
      
      // Draw labels to centrality axis
      preliminaryText->SetTextSize(0.05);
      preliminaryText->DrawLatexNDC(0.21, 0.1, "0#minus10%");
      preliminaryText->DrawLatexNDC(0.47, 0.1, "10#minus30%");
      preliminaryText->DrawLatexNDC(0.75, 0.1, "30#minus50%");
      
      // Save the figures to file
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/finalJetV%d%s.%s", iFlow+1, saveComment.Data(), figureFormat));
      }
    } // Flow component loop
  } // Drawing individual graphs
  
  // ========================================================================
  // == Drawing style, where each histogram is drawn into one small canvas ==
  // ========================================================================
  
  if(drawSmallCanvas){
    
    double xAdder[] = {-0.2,-0.65,-0.5,-0.35};
    int flowOrder[] = {3,2,4};
    int iFlow;
    
    // Only one legend for the plot
    legend = new TLegend(0.2, 0.625, 0.7, 0.805);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.045);legend->SetTextFont(62);
    
    // Loop over flows in a funny order
    for(int flowIndex = 0; flowIndex < 3; flowIndex++){
      
      iFlow = flowOrder[flowIndex]-1;
      
      // Displace the points a bit in x-direction for nicer drawing result
      for(int iPoint = 0; iPoint < jetVnUncertainty[iFlow]->GetN(); iPoint++){
        currentX = jetVnUncertainty[iFlow]->GetPointX(iPoint);
        jetVnUncertainty[iFlow]->SetPointX(iPoint, currentX+xAdder[iFlow]);
        jetVnGraph[iFlow]->SetPointX(iPoint, currentX+xAdder[iFlow]);
      }
      
      if(jetVnUncertainty[iFlow] != NULL){
        
        // Set the style for uncertainties
        for(int iCentrality = 0; iCentrality < jetVnUncertainty[iFlow]->GetN(); iCentrality++){
          errorY = jetVnUncertainty[iFlow]->GetErrorY(iCentrality);
          jetVnUncertainty[iFlow]->SetPointError(iCentrality, 0.075, errorY);
        }
        
        //jetVnUncertainty[iFlow]->GetYaxis()->SetNdivisions(510);
        jetVnUncertainty[iFlow]->SetMarkerSize(4);
        if(flowIndex == 0){
          drawer->DrawGraphCustomAxes(jetVnUncertainty[iFlow], 0, 4, minZoom[iFlow], maxZoom[iFlow], "Centrality", "Dijet v_{n}", " ", "a,e2");
        } else {
          jetVnUncertainty[iFlow]->Draw("same,e2");
        }
        
        jetVnGraph[iFlow]->SetMarkerSize(4);
        jetVnGraph[iFlow]->Draw("p,same");
                
      } else {
        if(flowIndex == 0){
          drawer->DrawGraphCustomAxes(jetVnGraph[iFlow], 0, 4, minZoom[iFlow], maxZoom[iFlow], "Centrality", "Dijet v_{n}", " ", "ap");
        } else {
          jetVnGraph[iFlow]->Draw("same,p");
        }
      }
      
    } // Flow component loop
    
    // Add legend to entry in different order
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      legend->AddEntry(jetVnUncertainty[iFlow], Form("Dijet v_{%d}", iFlow+1), "pf");
    }
    
    // Draw illustrative line and the legend
    zeroLine->Draw();
    legend->Draw();
    
    // Draw CMS supplementary tag
    preliminaryText->SetTextFont(62);
    preliminaryText->SetTextSize(0.06);
    preliminaryText->DrawLatexNDC(0.08, 0.941, "CMS");
    
    preliminaryText->SetTextFont(52);
    preliminaryText->SetTextSize(0.055);
    preliminaryText->DrawLatexNDC(0.215, 0.941, "Supplementary");
    
    preliminaryText->SetTextFont(42);
    preliminaryText->SetTextSize(0.05);
    preliminaryText->DrawLatexNDC(0.625, 0.941, "arXiv:2210.08325");
    
    // Draw luminosity and selection information to the figure
    preliminaryText->DrawLatexNDC(0.23, 0.87, "PbPb #sqrt{s_{NN}} = 5.02 TeV, 1.69 nb^{-1}");
    
    preliminaryText->SetTextSize(0.04);
    preliminaryText->DrawLatexNDC(0.57, 0.77, "0.7 < p^{ch}_{T} < 3 GeV");
    preliminaryText->DrawLatexNDC(0.57, 0.705, "anti-k_{T} R = 0.4");
    preliminaryText->DrawLatexNDC(0.57, 0.64, "|#eta_{jet}| < 1.3");
    preliminaryText->DrawLatexNDC(0.23, 0.275, "p_{T}^{lead} > 120 GeV");
    preliminaryText->DrawLatexNDC(0.23, 0.21, "p_{T}^{sub} > 50 GeV");
    preliminaryText->DrawLatexNDC(0.57, 0.21, "#Delta#varphi > #frac{5#pi}{6}");
    
    // Draw labels to centrality axis
    preliminaryText->SetTextSize(0.05);
    preliminaryText->DrawLatexNDC(0.21, 0.1, "0#minus10%");
    preliminaryText->DrawLatexNDC(0.47, 0.1, "10#minus30%");
    preliminaryText->DrawLatexNDC(0.75, 0.1, "30#minus50%");
    
    // Save the figures to file
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/smallCanvasJetVn%s.%s", saveComment.Data(), figureFormat));
    }
    
  } // Drawing everything in one small canvas
  
  // =============================================
  // ===== Drawing style with one big canvas =====
  // =============================================
  
  if(drawBigCanvas){
    
    cmsHighPtV2->SetMarkerSize(2.8);
    
    // Draw all the distributions to big canvases
    auxi_canvas *bigCanvas;
    TLatex *mainTitle;
    
    // Draw a big canvas and put all the plots in it
    bigCanvas = new auxi_canvas("bigCanvas", "", 1500, 550);
    bigCanvas->SetMargin(0.07, 0.01, 0.16, 0.01);
    bigCanvas->divide(1,3);
    
    mainTitle = new TLatex();
    
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      bigCanvas->CD(iFlow);
      
      
      // Adjust marker sizes
      jetVnUncertainty[iFlow]->SetMarkerSize(2);
      jetVnGraph[iFlow]->SetMarkerSize(2);
      
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
        jetVnUncertainty[iFlow]->GetXaxis()->SetLabelOffset(10);
      } else {
        jetVnUncertainty[iFlow]->GetXaxis()->SetTitleOffset(1.25);
        jetVnUncertainty[iFlow]->GetXaxis()->SetLabelOffset(10);
      }
      
      // Draw the graphs
      jetVnUncertainty[iFlow]->Draw("a,e2");
      jetVnGraph[iFlow]->Draw("same,p");
      
      // Create a legend for the plot division
      if(iFlow > 1) {
        legend = new TLegend(0.06,0.89,0.7,0.98);
        anotherLegend = new TLegend(0.23,0.7,0.73,0.85);
      } else {
        legend = new TLegend(0.23,0.89,0.73,0.98);
        anotherLegend = new TLegend(0.23,0.77,0.73,0.89);
      }
      
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.054);legend->SetTextFont(62);
      anotherLegend->SetFillStyle(0);anotherLegend->SetBorderSize(0);anotherLegend->SetTextSize(0.054);anotherLegend->SetTextFont(62);
      
      if(iFlow > 1) legend->SetTextSize(0.06);  // Need to increase text size for smaller divisions
      
      legend->AddEntry(jetVnUncertainty[iFlow], Form("Dijet v_{%d}", iFlow+1), "pf");
      
      if(iFlow == 1 && drawAtlasJetV2){
        atlasJetV2graph->SetMarkerSize(2);
        atlasJetV2graph->Draw("p,same");
        legend->AddEntry(atlasJetV2graph, "ATLAS v_{2}", "p"); // (#scale[0.8]{HP 2020})
      }
      
      if(iFlow == 1 && drawCmsHigtPtV2){
        cmsHighPtV2->Draw("p,same");
        anotherLegend->AddEntry(cmsHighPtV2, "CMS charged hadron v_{2}", "p"); // (#scale[0.8]{PLB 776 (2018) 195})
        anotherLegend->AddEntry((TObject*)0, "p_{T} > 20 GeV, |#eta| < 1", "");
      }
      
      legend->Draw();
      
      zeroLine->Draw();
      
      if(iFlow == 1) anotherLegend->Draw();
      
    } // Flow component loop
    
    // Draw all necessary CMS text to the plot
    bigCanvas->cd(0);
    
    mainTitle->SetTextFont(62);
    mainTitle->SetTextSize(0.065);
    cmsYPosition = 0.89;
    if(drawPreliminaryTag) cmsYPosition = 0.91;
    mainTitle->DrawLatexNDC(0.9, cmsYPosition, "CMS");
    
    mainTitle->SetTextFont(42);
    mainTitle->SetTextSize(0.055);
    if(drawPreliminaryTag)mainTitle->DrawLatexNDC(0.878, 0.86, "Preliminary");
    mainTitle->DrawLatexNDC(0.715, 0.76, "PbPb #sqrt{s_{NN}} = 5.02 TeV, 1.69 nb^{-1}");
    mainTitle->DrawLatexNDC(0.54, 0.92, "anti-k_{T} R = 0.4");
    mainTitle->DrawLatexNDC(0.54, 0.84, "|#eta_{jet}| < 1.3");
    if(leadSubTitles){
    mainTitle->DrawLatexNDC(0.54, 0.76, "p_{T}^{lead} > 120 GeV");
    mainTitle->DrawLatexNDC(0.54, 0.68, "p_{T}^{sub} > 50 GeV");
    mainTitle->DrawLatexNDC(0.54, 0.60, "#Delta#varphi > #frac{5#pi}{6}");
    } else {
      mainTitle->DrawLatexNDC(0.54, 0.76, "p_{T,1} > 120 GeV");
      mainTitle->DrawLatexNDC(0.54, 0.68, "p_{T,2} > 50 GeV");
      mainTitle->DrawLatexNDC(0.54, 0.60, "|#Delta#varphi_{1,2}| > #frac{5#pi}{6}");
    }
    mainTitle->DrawLatexNDC(0.13, 0.33, "Factorization region:");
    mainTitle->DrawLatexNDC(0.13, 0.25, "0.7 < Hadron p_{T} < 3 GeV");
    
    // Print the centrality bins
    mainTitle->SetTextFont(42);
    mainTitle->SetTextSize(0.045);
    preliminaryText->DrawLatexNDC(0.09, 0.1, "0#minus10%");
    preliminaryText->DrawLatexNDC(0.2, 0.1, "10#minus30%");
    preliminaryText->DrawLatexNDC(0.31, 0.1, "30#minus50%");
    preliminaryText->DrawLatexNDC(0.395, 0.1, "0#minus10%");
    preliminaryText->DrawLatexNDC(0.51, 0.1, "10#minus30%");
    preliminaryText->DrawLatexNDC(0.615, 0.1, "30#minus50%");
    preliminaryText->DrawLatexNDC(0.70, 0.1, "0#minus10%");
    preliminaryText->DrawLatexNDC(0.81, 0.1, "10#minus30%");
    preliminaryText->DrawLatexNDC(0.925, 0.1, "30#minus50%");
    
    // Save the figures to file
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/finalBigCanvasJetVn%s.%s", saveComment.Data(), figureFormat));
    }
    
  }
  
  // Print the central values of the points to the console
  if(printCentralValues){
    
    double vnValue;
    int centralityBinBorders[] = {0, 10, 30, 50, 90};
    
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      
      for(int iPoint = 0; iPoint < 3; iPoint++){
        vnValue = jetVnGraph[iFlow]->GetPointY(iPoint);
        cout << "Dijet v" << iFlow+1 << " C: " << centralityBinBorders[iPoint] << "-" << centralityBinBorders[iPoint+1] << " = " << vnValue << endl;
      }
      
    }
  }
  
  // Save the histograms to a file for HepData submission
  if(saveGraphsForHepData){
    TString outputFileName = "hepdata/hepdata_dijetVnCentrality_hin-21-002.root";
    TFile *outputFile = TFile::Open(outputFileName,"UPDATE");
    
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      jetVnGraph[iFlow]->Write(Form("dijetV%dCentrality", iFlow+1), TObject::kOverwrite);
      jetVnUncertainty[iFlow]->Write(Form("dijetV%dCentralityError", iFlow+1), TObject::kOverwrite);
    }
    
    outputFile->Close();
  }
}
