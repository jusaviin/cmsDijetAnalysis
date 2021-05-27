#include "LongRangeSystematicOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

/*
 * Macro for plotting final jet v2 and v3 results with systematic uncertainties.
 */
void finalLongRangePlotter(){

  // ============= //
  // Configuration //
  // ============= //
  
  // Input file name for data
  TString directoryName = "flowGraphs/";
  TString inputFileName = "summaryPlot_akCaloJet_nominalCorrection_2021-05-20.root";
  TString uncertaintyFileName = "systematicUncertainties_justTestingTwo.root";
  
  // Text to be put into legend for the input graphs
  const char* legendText = "Nominal calo jets";
  
  // Define the bins that are drawn
  const int nCentralityBins = 3;  // Number of drawn centrality bins
  
  const int maxVn = 4;            // Maximum defined vn. Plots are made upto v4.
  const int firstDrawnVn = 2;     // First drawn flow component
  const int lastDrawnVn = 2;      // Last drawn flow component
  
  // Define if previous results should be included in the plot
  const bool drawAtlasJetV2 = true;
  const bool drawCmsHigtPtV2 = true;
  
  // Save the final plots
  const bool saveFigures = false;
  TString saveComment = "";
  
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
  
  // Draw the graphs for selected flow components
  TLegend *legend;
  double errorY;
  
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
      drawer->DrawGraphCustomAxes(jetVnUncertainty[iFlow], 0, 4, 0, 0.1, "Centrality", Form("Jet v_{%d}", iFlow+1), " ", "a,e2");
      
      jetVnGraph[iFlow]->Draw("p,same");
    } else {
      drawer->DrawGraphCustomAxes(jetVnGraph[iFlow], 0, 4, 0, 0.1, "Centrality", Form("Jet v_{%d}", iFlow+1), " ", "ap");
    }
    
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
}
