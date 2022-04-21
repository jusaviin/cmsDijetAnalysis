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
  const int nFiles = 2;
  TString uncertaintyFileName[] = { "systematicUncertainties_multiMatchNominal_scalingUpdates_2022-03-04.root", "systematicUncertainties_multiMatchNominal_addPtBinVariation_2022-04-21.root",    "systematicUncertainties_addMultiplicityMatch_noCentralityShift_finalCorrection_2022-28-01.root"};
  // systematicUncertainties_addMultiplicityMatch_noCentralityShift_finalCorrection_2022-28-01.root
  // systematicUncertainties_allSources_finalCorrection_2021-08-10.root
  
  TString legendComment[] = {"Updated scaling", "Add pT bin variation", "Test2", "Dummy"};
  
  // Define the bins that are drawn
  const int nCentralityBins = 3;  // Number of drawn centrality bins
  const double centralityBinBorders[] = {0, 10, 30, 50, 90}; // Bin borders for centrality bins
  
  const int maxVn = 4;            // Maximum defined vn. Plots are made upto v4.
  const int firstDrawnVn = 2;     // First drawn flow component
  const int lastDrawnVn = 4;      // Last drawn flow component
  
  // Choose which histograms to draw
  const bool drawRelativeUncertainties = false;
  const bool drawAbsoluteUncertainties = true;
  const bool groupUncertainties = false; // Group different uncertainty sources based on similarity
  const int groupStrategy = 0;  // Strategy used to group the uncertainty histograms. 0 = Paper. 1 = AN.
  
  // Save the final plots
  const bool saveFigures = false;
  TString saveComment = "_scalingUpdate";
  TString fileFormat = "pdf";
  
  // Print the uncertainties to a slide
  const bool printUncertainties = false;
  const int printUncertaintyType = 0; // 0 = Absolute uncertainties. 1 = Relative uncertainties
  
  // =========== //
  // Read graphs //
  // =========== //
  
  // First, read the input files and create systematic uncertainty organizers
  TFile *uncertaintyFile[nFiles];
  LongRangeSystematicOrganizer *uncertaintyOrganizer[nFiles];
  for(int iFile = 0; iFile < nFiles; iFile++){
    uncertaintyFile[iFile] = TFile::Open(directoryName+uncertaintyFileName[iFile]);
    uncertaintyOrganizer[iFile] = new LongRangeSystematicOrganizer(uncertaintyFile[iFile]);
    uncertaintyOrganizer[iFile]->SetGroupingStrategy(groupStrategy);
  }
  
  // Determine the number of uncertainty sources from the file
  const int nUncertaintySources = uncertaintyOrganizer[0]->GetNUncertaintySources(groupUncertainties);
  
  // Create graphs for each uncertainty source
  TGraphErrors *jetVnUncertainty[nFiles][nUncertaintySources][maxVn];
  
  // Initialize the graphs to NULL
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iUncertainty = 0; iUncertainty < nUncertaintySources; iUncertainty++){
      for(int iFlow = 0; iFlow < maxVn; iFlow++){
        jetVnUncertainty[iFile][iUncertainty][iFlow] = NULL;
      } // Flow component loop
    } // Uncertainty loop
  } // File loop
  
  // Read the systematic uncertainty graphs
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iUncertainty = 0; iUncertainty < nUncertaintySources; iUncertainty++){
      for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
        jetVnUncertainty[iFile][iUncertainty][iFlow] = uncertaintyOrganizer[iFile]->GetLongRangeSystematicUncertainty(iFlow, iUncertainty, -1, groupUncertainties);
        
        // TODO: Maybe some error handling if some uncertainties are not included in the file
        if(jetVnUncertainty[iFile][iUncertainty][iFlow] == NULL){
          cout << "So, this is awkward..." << endl;
          cout << "It seems that I could not find the graph for v" << iFlow+1 << " " << uncertaintyOrganizer[iFile]->GetLongRangeUncertaintyName(iUncertainty, groupUncertainties) << endl;
        }
        
      } // Flow component loop
    } // Uncertainty loop
  } // File loop
  
  // ============================================================================================== //
  // Transform the uncertainty graphs into histograms illustrating different sources of uncertainty //
  // ============================================================================================== //
  
  TH1D* systematicIllustrationAbsolute[nFiles][maxVn][nCentralityBins];
  TH1D* systematicIllustrationRelative[nFiles][maxVn][nCentralityBins];
  TString histogramName;
  
  // Initilize the histograms
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iFlow = 0; iFlow < maxVn; iFlow++){
        
        histogramName = Form("systematicIllustrationAbsolute%d_v%d_C%d", iFile, iFlow+1, iCentrality);
        systematicIllustrationAbsolute[iFile][iFlow][iCentrality] = new TH1D(histogramName, histogramName, nUncertaintySources, -0.5, nUncertaintySources-0.5);
        systematicIllustrationAbsolute[iFile][iFlow][iCentrality]->Sumw2();
        
        histogramName = Form("systematicIllustrationRelative%d_v%d_C%d", iFile, iFlow+1, iCentrality);
        systematicIllustrationRelative[iFile][iFlow][iCentrality] = new TH1D(histogramName, histogramName, nUncertaintySources, -0.5, nUncertaintySources-0.5);
        systematicIllustrationRelative[iFile][iFlow][iCentrality]->Sumw2();
        
        for(int iUncertainty = 0; iUncertainty < nUncertaintySources; iUncertainty++){
          
          // Change the bin labels to refer to different sources of uncertainty
          systematicIllustrationAbsolute[iFile][iFlow][iCentrality]->GetXaxis()->ChangeLabel(iUncertainty+1, 315, -1, -1, -1, -1, uncertaintyOrganizer[iFile]->GetUncertaintyAxisName(iUncertainty, groupUncertainties));
          systematicIllustrationRelative[iFile][iFlow][iCentrality]->GetXaxis()->ChangeLabel(iUncertainty+1, 315, -1, -1, -1, -1, uncertaintyOrganizer[iFile]->GetUncertaintyAxisName(iUncertainty, groupUncertainties));
          
        } // Uncertainty loop
      } // Flow component loop
    } // Centrality loop
  } // File loop
  
  double valueY, errorY, ratioY;
  
  // Fill the histograms
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
        for(int iUncertainty = 0; iUncertainty < nUncertaintySources; iUncertainty++){
          
          // For the absolute uncertainty, read the number directly from the uncertainty histogram
          errorY = jetVnUncertainty[iFile][iUncertainty][iFlow]->GetErrorY(iCentrality);
          systematicIllustrationAbsolute[iFile][iFlow][iCentrality]->SetBinContent(iUncertainty+1,errorY);
          
          // For the relative uncertainty, calculate the ratio between uncertainty and central value
          jetVnUncertainty[iFile][iUncertainty][iFlow]->GetPoint(iCentrality, ratioY, valueY);
          ratioY = errorY / TMath::Abs(valueY);
          systematicIllustrationRelative[iFile][iFlow][iCentrality]->SetBinContent(iUncertainty+1,ratioY);
          
        } // Uncertainty loop
      } // Flow component loop
    } // Centrality loop
  } // File loop
  
  // =========================================================================== //
  // Print slides summarizing all the different sources of systemtic uncertainty //
  // =========================================================================== //
  
  if(printUncertainties){
    
    // Print a slide with uncertainties for each source and each centrality for each V
    char namer[100];
    double nominalX, nominalY;
    const char* uncertaintyType[2] = {"absolute","relative"};
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      
      sprintf(namer,"\\frametitle{Uncertainty for jet $v_{%d}$}", iFlow+1);
      
      cout << endl;
      cout << "\\begin{frame}" << endl;
      cout << namer << endl;
      cout << "\\begin{center}" << endl;
      cout << "  \\begin{tabular}{cccccc}" << endl;
      cout << "    \\toprule" << endl;
      cout << "    Source & C: 0-10 & C: 10-30 & C: 30-50 \\\\" << endl;
      cout << "    \\midrule" << endl;
      
      // Set the correct precision for printing floating point numbers
      cout << fixed << setprecision(4);
      
      // First, print the actual values of vn:s
      cout << "$v_{" << iFlow+1 << "}$ ";
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        jetVnUncertainty[0][0][iFlow]->GetPoint(iCentrality, nominalX, nominalY);
        cout << " & $" << nominalY << "$";
      }
      cout << " \\\\" << endl;
      cout << "    \\midrule" << endl;
      
      if(printUncertaintyType > 0){
        
        // Print the relative uncertainties
        for(int iUncertainty = 0; iUncertainty < nUncertaintySources; iUncertainty++){
          cout << uncertaintyOrganizer[0]->GetLongRangeUncertaintyName(iUncertainty, groupUncertainties).Data();
          for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
            cout << " & $" << systematicIllustrationRelative[0][iFlow][iCentrality]->GetBinContent(iUncertainty+1) << "$";
          }
          cout << " \\\\" << endl;
        }
        
      } else {
        // Print the absolute uncertainties
        
        for(int iUncertainty = 0; iUncertainty < nUncertaintySources; iUncertainty++){
          cout << uncertaintyOrganizer[0]->GetLongRangeUncertaintyName(iUncertainty, groupUncertainties).Data();
          for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
            cout << " & $" << systematicIllustrationAbsolute[0][iFlow][iCentrality]->GetBinContent(iUncertainty+1) << "$";
          }
          cout << " \\\\" << endl;
        }
        
      }
      
      cout << "    \\bottomrule" << endl;
      cout << "  \\end{tabular}" << endl;
      cout << "\\end{center}" << endl;
      cout << "\\begin{itemize}" << endl;
      cout << "  \\item Top: value. Bottom: " << uncertaintyType[printUncertaintyType] << " uncertainty." << endl;
      cout << "\\end{itemize}" << endl;
      cout << "\\end{frame}" << endl;
      cout << endl;
      
    } // vn loop
    
  } // Printing uncertainties
  
  // ============== //
  // Draw the plots //
  // ============== //
  
  // Setup the drawer for the uncertainty sources
  JDrawer *drawer = new JDrawer();
  drawer->SetCanvasSize(1000,600);
  drawer->SetLeftMargin(0.1);
  drawer->SetRightMargin(0.01);
  drawer->SetBottomMargin(0.3);
  drawer->SetTitleOffsetY(0.9);
  drawer->SetLabelOffsetX(0.09);
  drawer->SetTitleOffsetX(2.8);
  drawer->SetNDivisionsX(100 + nUncertaintySources);
  
  TLegend *legend;
  double maxValue;
  double highValue;
  int lineColor[] = {kBlue, kRed, kCyan, kMagenta, kGreen+3};
  
  // Draw the relative uncertainties
  if(drawRelativeUncertainties){
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
                
        legend = new TLegend(0.11,0.7,0.41,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        
        // Find the maximum values over all drawn histograms
        maxValue = 0;
        for(int iFile = 0; iFile < nFiles; iFile++){
          systematicIllustrationRelative[iFile][iFlow][iCentrality]->SetLineColor(lineColor[iFile]);
          highValue = systematicIllustrationRelative[iFile][iFlow][iCentrality]->GetBinContent(nUncertaintySources);
          highValue = highValue + 0.1 * highValue;
          if(highValue > maxValue) maxValue = highValue;
        }
        systematicIllustrationRelative[0][iFlow][iCentrality]->GetYaxis()->SetRangeUser(0,maxValue);
        
        drawer->DrawHistogram(systematicIllustrationRelative[0][iFlow][iCentrality], "Source", "Relative uncertainty", " ");
        for(int iFile = 1; iFile < nFiles; iFile++){
          systematicIllustrationRelative[iFile][iFlow][iCentrality]->Draw("same");
        }
        
        legend->AddEntry((TObject*) 0, Form("Jet v_{%d}", iFlow+1), "");
        legend->AddEntry((TObject*) 0, Form("Cent: %.0f-%.0f%%", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "");
        if(nFiles > 1){
          for(int iFile = 0; iFile < nFiles; iFile++){
            legend->AddEntry(systematicIllustrationRelative[iFile][iFlow][iCentrality], legendComment[iFile], "l");
          }
        }
        legend->Draw();
        
        // Save the figures to file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/relativeUncertainties%s_jetV%d_C=%.0f-%.0f.%s", saveComment.Data(), iFlow+1, centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], fileFormat.Data()));
        }
      } // Centrality loop
    } // Flow component loop
  } // Drawing relative uncertainty plots
  
  // Draw the absolute uncertainties
  if(drawAbsoluteUncertainties){
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
                
        legend = new TLegend(0.11,0.7,0.41,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        
        // Find the maximum values over all drawn histograms
        maxValue = 0;
        for(int iFile = 0; iFile < nFiles; iFile++){
          systematicIllustrationAbsolute[iFile][iFlow][iCentrality]->SetLineColor(lineColor[iFile]);
          highValue = systematicIllustrationAbsolute[iFile][iFlow][iCentrality]->GetBinContent(nUncertaintySources);
          highValue = highValue + 0.1 * highValue;
          if(highValue > maxValue) maxValue = highValue;
        }
        systematicIllustrationAbsolute[0][iFlow][iCentrality]->GetYaxis()->SetRangeUser(0,maxValue);
        
        drawer->DrawHistogram(systematicIllustrationAbsolute[0][iFlow][iCentrality], "Source", "Absolute uncertainty", " ");
        for(int iFile = 1; iFile < nFiles; iFile++){
          systematicIllustrationAbsolute[iFile][iFlow][iCentrality]->Draw("same");
        }
        
        legend->AddEntry((TObject*) 0, Form("Jet v_{%d}", iFlow+1), "");
        legend->AddEntry((TObject*) 0, Form("Cent: %.0f-%.0f%%", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "");
        if(nFiles > 1){
          for(int iFile = 0; iFile < nFiles; iFile++){
            legend->AddEntry(systematicIllustrationAbsolute[iFile][iFlow][iCentrality], legendComment[iFile], "l");
          }
        }
        legend->Draw();
        
        // Save the figures to file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/absoluteUncertainties%s_jetV%d_C=%.0f-%.0f.%s", saveComment.Data(), iFlow+1, centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], fileFormat.Data()));
        }
      } // Centrality loop
    } // Flow component loop
  } // Drawing relative uncertainty plots
}
