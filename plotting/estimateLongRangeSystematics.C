// #include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

/*
 * Macro for estimating systematic uncertainties for long range correlation results
 *
 * The following sources of systematic uncertainty are estimated:
 *  - Adjustment on the gluing point of leading/subleading distributions (TODO)
 *  - Fit projections over regions -2.5 < deltaEta < -1.5 and 1.5 < deltaEta < 2.5 and compare to nominal (TODO)
 *  - Fit projections over regions 1.5 < |deltaEta| < 2 and 2 < |deltaEta| < 2.5 and compare to nominal (TODO)
 *  - Fit projections before and after mixed event correction (TODO)
 *  - Fit different v_{z} selection (like positive and negative) (TODO)
 *  - Differences coming from selecting different jet collections
 *
 *  The results will be saved to a file and also a slide with different contribution can be printed
 */
void estimateLongRangeSystematics(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Indices for currently implemented uncertainty sources
  enum enumLongRangeUncertaintySources{kJetCollection, kAll, knLongRangeUncertaintySources};
  TString uncertaintyNameString[knLongRangeUncertaintySources] = {"jetCollection","all"};
  
  // Nominal result file
  TString nominalResultFileName = "flowGraphs/summaryPlot_akCaloJet_nominalCorrection_2021-05-20.root";
  
  // Results using pfCs jets instead of calo jets
  TString pfcsResultFileName = "flowGraphs/summaryPlot_akPfCsJet_correctionOnlyQcut2p5_2021-05-20.root";
  
  // Flow component configuration
  const int maxVn = 4;      // Maximum number of flow components
  const int firstFlow = 2;  // First analyzed flow component
  const int lastFlow = 2;   // Last analyzed flow component
  
  const bool printUncertainties = true;
  
  bool plotExample = false;
  
  TString outputFileName = "flowGraphs/systematicUncertainties_justTesting.root";
  
  // ==================================================================
  // ====================== Configuration done ========================
  // ==================================================================
  
  // Open the input files
  TFile* nominalResultFile = TFile::Open(nominalResultFileName);
  TFile* pfcsResultFile = TFile::Open(pfcsResultFileName);
  
  // Result graphs
  TGraphErrors* nominalResultGraph[maxVn];
  TGraphErrors* pfcsResultGraph[maxVn];
  
  TString graphName;
  
  // Load the graphs from the files
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    
    graphName = Form("summaryV%d", iFlow+1);
    
    nominalResultGraph[iFlow] = (TGraphErrors*) nominalResultFile->Get(graphName);
    pfcsResultGraph[iFlow] = (TGraphErrors*) pfcsResultFile->Get(graphName);
    
  }
  
  const int nCentralityBins = nominalResultGraph[firstFlow-1]->GetN();
  
  // Make a big table with all relative and absolute uncertainties
  double relativeUncertaintyTable[knLongRangeUncertaintySources][maxVn][nCentralityBins];
  double absoluteUncertaintyTable[knLongRangeUncertaintySources][maxVn][nCentralityBins];
  
  // Initialize the tables to zero
  for(int iUncertainty = 0; iUncertainty < knLongRangeUncertaintySources; iUncertainty++){
    for(int iFlow = 0; iFlow < maxVn; iFlow++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        relativeUncertaintyTable[iUncertainty][iFlow][iCentrality] = 0;
        absoluteUncertaintyTable[iUncertainty][iFlow][iCentrality] = 0;
      }
    }
  }
  
  // Helper variables for reading points from graphs
  double nominalX, nominalY, nominalYerror;
  double comparisonX, comparisonY, comparisonYerror;
  
  // ================================================ //
  // Uncertainty coming from jet collection selection //
  // ================================================ //
  
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      nominalResultGraph[iFlow]->GetPoint(iCentrality, nominalX, nominalY);
      pfcsResultGraph[iFlow]->GetPoint(iCentrality, comparisonX, comparisonY);
      
      absoluteUncertaintyTable[kJetCollection][iFlow][iCentrality] = TMath::Abs(nominalY-comparisonY);
      relativeUncertaintyTable[kJetCollection][iFlow][iCentrality] = TMath::Abs(1 - comparisonY/nominalY);
      
    }
  }
  
  // ========================================= //
  // Add all uncertainty sources in quadrature //
  // ========================================= //
  
  double sumOfSquares;
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      nominalResultGraph[iFlow]->GetPoint(iCentrality, nominalX, nominalY);
      
      sumOfSquares = 0;
      for(int iUncertainty = 0; iUncertainty < kAll; iUncertainty++){
        sumOfSquares = sumOfSquares + TMath::Power(absoluteUncertaintyTable[iUncertainty][iFlow][iCentrality],2);
      }
    
      absoluteUncertaintyTable[kAll][iFlow][iCentrality] = TMath::Sqrt(sumOfSquares);
      relativeUncertaintyTable[kAll][iFlow][iCentrality] = TMath::Abs(absoluteUncertaintyTable[kAll][iFlow][iCentrality]/nominalY);
      
    }
  }
  
  // ================================================ //
  // Collect the calculated uncertainties into graphs //
  // ================================================ //
  
  TGraphErrors* systematicUncertaintyGraph[knLongRangeUncertaintySources][maxVn];
  
  for(int iUncertainty = 0; iUncertainty < knLongRangeUncertaintySources; iUncertainty++){
    for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
      graphName = Form("systematicUncertainty_v%d_%s", iFlow+1, uncertaintyNameString[iUncertainty].Data());
      systematicUncertaintyGraph[iUncertainty][iFlow] = (TGraphErrors*) nominalResultGraph[iFlow]->Clone(graphName);
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        systematicUncertaintyGraph[iUncertainty][iFlow]->SetPointError(iCentrality, 0, absoluteUncertaintyTable[iUncertainty][iFlow][iCentrality]);
      }
    }
  }
  
  // =================================================== //
  // Write the systematic uncertainty graphs into a file //
  // =================================================== //
  
  // If an output file name is given, put the total uncertainties to an output file
  if(outputFileName.EndsWith(".root")){
    TFile *outputFile = new TFile(outputFileName,"UPDATE");

    for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
      for(int iUncertainty = 0; iUncertainty < knLongRangeUncertaintySources; iUncertainty++){
        graphName = Form("systematicUncertainty_v%d_%s", iFlow+1, uncertaintyNameString[iUncertainty].Data());
        systematicUncertaintyGraph[iUncertainty][iFlow]->Write(graphName, TObject::kOverwrite);
      } // Uncertainty loop
    } // Flow component loop

    outputFile->Close();
  }
  
  // =========================================================================== //
  // Print slides summarizing all the different sources of systemtic uncertainty //
  // =========================================================================== //
  
  if(printUncertainties){
    
    // Print a slide with uncertainties for each source and each centrality for each V
    char namer[100];
    for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
      
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
        nominalResultGraph[iFlow]->GetPoint(iCentrality, nominalX, nominalY);
        cout << " & $" << nominalY << "$";
      }
      cout << " \\\\" << endl;
      cout << "    \\midrule" << endl;
      
      // Print the relative uncertainties
      
      for(int iUncertainty = 0; iUncertainty < knLongRangeUncertaintySources; iUncertainty++){
        cout << uncertaintyNameString[iUncertainty].Data();
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          cout << " & $" << relativeUncertaintyTable[iUncertainty][iFlow][iCentrality] << "$";
        }
        cout << " \\\\" << endl;
      }
      
      cout << "    \\midrule" << endl;
      
      // Print the absolute uncertainties
      
      for(int iUncertainty = 0; iUncertainty < knLongRangeUncertaintySources; iUncertainty++){
        cout << uncertaintyNameString[iUncertainty].Data();
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          cout << " & $" << absoluteUncertaintyTable[iUncertainty][iFlow][iCentrality] << "$";
        }
        cout << " \\\\" << endl;
      }
      
      
      cout << "    \\bottomrule" << endl;
      cout << "  \\end{tabular}" << endl;
      cout << "\\end{center}" << endl;
      cout << "\\begin{itemize}" << endl;
      cout << "  \\item Top: value. Middle: relative uncertainty. Bottom: absolute uncertainty." << endl;
      cout << "\\end{itemize}" << endl;
      cout << "\\end{frame}" << endl;
      cout << endl;
      
    } // vn loop
    
  } // Printing uncertainties
}
