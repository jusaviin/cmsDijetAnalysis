#include "LongRangeSystematicOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

/*
 * Macro for visualizing the systematic uncertainties coming from different sources.
 */
void longRangeSystematicExplorer(){

  // ============= //
  // Configuration //
  // ============= //
  
  // Input file name for uncertainties
  TString directoryName = "flowGraphs/";
  TString uncertaintyFileName = "systematicUncertainties_justTestingWithFun_finalCorrection.root";
  
  // Define the bins that are drawn
  const int nCentralityBins = 3;  // Number of drawn centrality bins
  const double centralityBinBorders[] = {0, 10, 30, 50, 90}; // Bin borders for centrality bins
  
  const int maxVn = 4;            // Maximum defined vn. Plots are made upto v4.
  const int firstDrawnVn = 2;     // First drawn flow component
  const int lastDrawnVn = 2;      // Last drawn flow component
  
  // Choose which histograms to draw
  const bool drawRelativeUncertainties = true;
  const bool drawAbsoluteUncertainties = false;
  
  // Save the final plots
  const bool saveFigures = false;
  TString saveComment = "_initialCheck";
  
  // =========== //
  // Read graphs //
  // =========== //
  
  TGraphErrors *jetVnUncertainty[LongRangeSystematicOrganizer::knUncertaintySources][maxVn];
  
  // Initialize the graphs to NULL
  for(int iUncertainty = 0; iUncertainty < LongRangeSystematicOrganizer::knUncertaintySources; iUncertainty++){
    for(int iFlow = 0; iFlow < maxVn; iFlow++){
      jetVnUncertainty[iUncertainty][iFlow] = NULL;
    } // Flow component loop
  } // Uncertainty loop
  
  // Read the systematic uncertainty graphs
  TFile *uncertaintyFile = TFile::Open(directoryName+uncertaintyFileName);
  LongRangeSystematicOrganizer *uncertaintyOrganizer = new LongRangeSystematicOrganizer(uncertaintyFile);
  
  for(int iUncertainty = 0; iUncertainty < LongRangeSystematicOrganizer::knUncertaintySources; iUncertainty++){
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      jetVnUncertainty[iUncertainty][iFlow] = uncertaintyOrganizer->GetLongRangeSystematicUncertainty(iFlow, iUncertainty);
      
      // TODO: Maybe some error handling if some uncertainties are not included in the file
      if(jetVnUncertainty[iUncertainty][iFlow] == NULL){
        cout << "So, this is awkward..." << endl;
        cout << "It seems that I could not find the graph for v" << iFlow+1 << " " << uncertaintyOrganizer->GetLongRangeUncertaintyName(iUncertainty) << endl;
      }
      
    } // Flow component loop
  } // Uncertainty loop
  
  // ============================================================================================== //
  // Transform the uncertainty graphs into histograms illustrating different sources of uncertainty //
  // ============================================================================================== //
  
  TH1D* systematicIllustrationAbsolute[maxVn][nCentralityBins];
  TH1D* systematicIllustrationRelative[maxVn][nCentralityBins];
  TString histogramName;
  
  // Initilize the histograms
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iFlow = 0; iFlow < maxVn; iFlow++){
      
      histogramName = Form("systematicIllustrationAbsolute_v%d_C%d", iFlow+1, iCentrality);
      systematicIllustrationAbsolute[iFlow][iCentrality] = new TH1D(histogramName, histogramName, LongRangeSystematicOrganizer::knUncertaintySources, -0.5, LongRangeSystematicOrganizer::knUncertaintySources-0.5);
      systematicIllustrationAbsolute[iFlow][iCentrality]->Sumw2();
      
      histogramName = Form("systematicIllustrationRelative_v%d_C%d", iFlow+1, iCentrality);
      systematicIllustrationRelative[iFlow][iCentrality] = new TH1D(histogramName, histogramName, LongRangeSystematicOrganizer::knUncertaintySources, -0.5, LongRangeSystematicOrganizer::knUncertaintySources-0.5);
      systematicIllustrationRelative[iFlow][iCentrality]->Sumw2();
      
      for(int iUncertainty = 0; iUncertainty < LongRangeSystematicOrganizer::knUncertaintySources; iUncertainty++){
        
        // Change the bin labels to refer to different sources of uncertainty
        systematicIllustrationAbsolute[iFlow][iCentrality]->GetXaxis()->ChangeLabel(iUncertainty+1, 315, -1, -1, -1, -1, uncertaintyOrganizer->GetUncertaintyAxisName(iUncertainty));
        systematicIllustrationRelative[iFlow][iCentrality]->GetXaxis()->ChangeLabel(iUncertainty+1, 315, -1, -1, -1, -1, uncertaintyOrganizer->GetUncertaintyAxisName(iUncertainty));
        
      } // Uncertainty loop
    } // Flow component loop
  } // Centrality loop
  
  double valueY, errorY, ratioY;
  
  // Fill the histograms
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      for(int iUncertainty = 0; iUncertainty < LongRangeSystematicOrganizer::knUncertaintySources; iUncertainty++){
        
        // For the absolute uncertainty, read the number directly from the uncertainty histogram
        errorY = jetVnUncertainty[iUncertainty][iFlow]->GetErrorY(iCentrality);
        systematicIllustrationAbsolute[iFlow][iCentrality]->SetBinContent(iUncertainty+1,errorY);
        
        // For the relative uncertainty, calculate the ratio between uncertainty and central value
        jetVnUncertainty[iUncertainty][iFlow]->GetPoint(iCentrality, ratioY, valueY);
        ratioY = errorY / valueY;
        systematicIllustrationRelative[iFlow][iCentrality]->SetBinContent(iUncertainty+1,ratioY);
        
      } // Uncertainty loop
    } // Flow component loop
  } // Centrality loop
  
  // ============== //
  // Draw the plots //
  // ============== //
  
  // Setup the drawer for the uncertainty sources
  JDrawer *drawer = new JDrawer();
  drawer->SetCanvasSize(1000,600);
  drawer->SetLeftMargin(0.1);
  drawer->SetRightMargin(0.01);
  drawer->SetBottomMargin(0.3);
  drawer->SetTitleOffsetY(0.8);
  drawer->SetLabelOffsetX(0.09);
  drawer->SetTitleOffsetX(2.8);
  drawer->SetNDivisionsX(100 + LongRangeSystematicOrganizer::knUncertaintySources);
  
  TLegend *legend;
  double maxValue;
  
  // Draw the relative uncertainties
  if(drawRelativeUncertainties){
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
                
        legend = new TLegend(0.2,0.7,0.5,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        
        maxValue = systematicIllustrationRelative[iFlow][iCentrality]->GetBinContent(LongRangeSystematicOrganizer::kAll+1);
        maxValue = maxValue + 0.1 * maxValue;
        systematicIllustrationRelative[iFlow][iCentrality]->GetYaxis()->SetRangeUser(0,maxValue);
        
        drawer->DrawHistogram(systematicIllustrationRelative[iFlow][iCentrality], "Source", "Relative uncertainty", " ");
        
        legend->AddEntry((TObject*) 0, Form("Cent: %.0f-%.0f%%", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "");
        legend->Draw();
        
        // Save the figures to file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/relativeUncertainties%s_jetV%d_C=%.0f-%.0f.pdf", saveComment.Data(), iFlow+1, centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
        }
      } // Centrality loop
    } // Flow component loop
  } // Drawing relative uncertainty plots
  
  // Draw the absolute uncertainties
  if(drawAbsoluteUncertainties){
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
                
        legend = new TLegend(0.2,0.7,0.5,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        
        maxValue = systematicIllustrationAbsolute[iFlow][iCentrality]->GetBinContent(LongRangeSystematicOrganizer::kAll+1);
        maxValue = maxValue + 0.1 * maxValue;
        systematicIllustrationAbsolute[iFlow][iCentrality]->GetYaxis()->SetRangeUser(0,maxValue);
        
        drawer->DrawHistogram(systematicIllustrationAbsolute[iFlow][iCentrality], "Source", "Absolute uncertainty", " ");
        
        legend->AddEntry((TObject*) 0, Form("Cent: %.0f-%.0f%%", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "");
        legend->Draw();
        
        // Save the figures to file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/absoluteUncertainties%s_jetV%d_C=%.0f-%.0f.pdf", saveComment.Data(), iFlow+1, centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
        }
      } // Centrality loop
    } // Flow component loop
  } // Drawing relative uncertainty plots
}
