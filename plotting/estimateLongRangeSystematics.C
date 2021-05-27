#include "LongRangeSystematicOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

/*
 * Macro for checking if two values are within statistical uncertainties
 *
 *  double nominalValue = Nominal jet vn value to be checked agains
 *  double comparisonValue = Jet vn value to be compared with the nominal
 *  double nominalError = Statistical uncertainty of the nominal vn value
 *  double comparisonError = Statistical uncertainty of the jet vn value compared with the nominal
 *
 *  return: True if the two values are within statistical uncertainties, false if not
 */
bool isWithinStatisticalErrors(double nominalValue, double comparisonValue, double nominalError, double comparisonError){
  return TMath::Abs(nominalValue-comparisonValue) - (nominalError + comparisonError) < 0;
}

/*
 * Macro for estimating systematic uncertainties for long range correlation results
 *
 * The following sources of systematic uncertainty are estimated:
 *  - Adjustment on the gluing point of leading/subleading distributions (TODO)
 *  - Fit projections over regions -2.5 < deltaEta < -1.5 and 1.5 < deltaEta < 2.5 and compare to nominal
 *  - Fit projections over regions 1.5 < |deltaEta| < 2 and 2 < |deltaEta| < 2.5 and compare to nominal
 *  - Fit different v_{z} selection (like positive and negative) (TODO)
 *  - Differences coming from selecting different jet collections
 *  - Jet axis smearing
 *
 *  The results will be saved to a file and also a slide with different contribution can be printed
 */
void estimateLongRangeSystematics(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Nominal result file
  TString nominalResultFileName = "flowGraphs/summaryPlot_akCaloJet_nominalCorrection_2021-05-20.root";
  
  // Results using pfCs jets instead of calo jets
  TString pfcsResultFileName = "flowGraphs/summaryPlot_akPfCsJet_correctionOnlyQcut2p5_2021-05-20.root";
  
  // Results with smeared jet axis
  TString angleSmearResultFileName = "flowGraphs/summaryPlot_akCaloJet_angleSmear_nominalCorrection_2021-05-24.root";
  
  // Results with modified deltaEta region for jet-hadron long range definition
  TString narrowDeltaEtaJetHadronResultFileName = "flowGraphs/summaryPlot_akCaloJet_narrowDeltaEtaJetHadron_nominalCorrection_2021-05-26.root";
  TString wideDeltaEtaJetHadronResultFileName = "flowGraphs/summaryPlot_akCaloJet_wideDeltaEtaJetHadron_nominalCorrection_2021-05-26.root";
  TString negativeDeltaEtaJetHadronResultFileName = "flowGraphs/summaryPlot_akCaloJet_negativeDeltaEtaJetHadron_nominalCorrection_2021-05-26.root";
  TString positiveDeltaEtaJetHadronResultFileName = "flowGraphs/summaryPlot_akCaloJet_positiveDeltaEtaJetHadron_nominalCorrection_2021-05-26.root";
  
  // Results with minimum bias dihadron
  TString minBiasDihadronResultFileName = "flowGraphs/summaryPlot_akCaloJet_minBiasDihadron_nominalCorrection_2021-05-26.root";
  
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
  TFile* angleSmearResultFile = TFile::Open(angleSmearResultFileName);
  TFile* narrowDeltaEtaJetHadronResultFile = TFile::Open(narrowDeltaEtaJetHadronResultFileName);
  TFile* wideDeltaEtaJetHadronResultFile = TFile::Open(wideDeltaEtaJetHadronResultFileName);
  TFile* negativeDeltaEtaJetHadronResultFile = TFile::Open(negativeDeltaEtaJetHadronResultFileName);
  TFile* positiveDeltaEtaJetHadronResultFile = TFile::Open(positiveDeltaEtaJetHadronResultFileName);
  TFile* minBiasDihadronResultFile = TFile::Open(minBiasDihadronResultFileName);
  
  // Result graphs
  TGraphErrors* nominalResultGraph[maxVn];
  TGraphErrors* pfcsResultGraph[maxVn];
  TGraphErrors* angleSmearResultGraph[maxVn];
  TGraphErrors* deltaEtaRegionJetHadronGraph[2][maxVn];
  TGraphErrors* deltaEtaSideJetHadronGraph[2][maxVn];
  TGraphErrors* minBiasDihadronResultGraph[maxVn];
  
  TString graphName;
  
  // Load the graphs from the files
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    
    graphName = Form("summaryV%d", iFlow+1);
    
    nominalResultGraph[iFlow] = (TGraphErrors*) nominalResultFile->Get(graphName);
    pfcsResultGraph[iFlow] = (TGraphErrors*) pfcsResultFile->Get(graphName);
    angleSmearResultGraph[iFlow] = (TGraphErrors*) angleSmearResultFile->Get(graphName);
    deltaEtaRegionJetHadronGraph[0][iFlow] = (TGraphErrors*) narrowDeltaEtaJetHadronResultFile->Get(graphName);
    deltaEtaRegionJetHadronGraph[1][iFlow] = (TGraphErrors*) wideDeltaEtaJetHadronResultFile->Get(graphName);
    deltaEtaSideJetHadronGraph[0][iFlow] = (TGraphErrors*) negativeDeltaEtaJetHadronResultFile->Get(graphName);
    deltaEtaSideJetHadronGraph[1][iFlow] = (TGraphErrors*) positiveDeltaEtaJetHadronResultFile->Get(graphName);
    minBiasDihadronResultGraph[iFlow] = (TGraphErrors*) minBiasDihadronResultFile->Get(graphName);
    
  }
  
  const int nCentralityBins = nominalResultGraph[firstFlow-1]->GetN();
  
  // Histogram showing if difference between nominal and uncertainty is less than the statistical uncertainties
  TH1D* statisticallyInsignificant[maxVn][nCentralityBins];
  
  // Make a big table with all relative and absolute uncertainties
  double relativeUncertaintyTable[LongRangeSystematicOrganizer::knUncertaintySources][maxVn][nCentralityBins];
  double absoluteUncertaintyTable[LongRangeSystematicOrganizer::knUncertaintySources][maxVn][nCentralityBins];
  
  // Initialize the tables to zero
  for(int iUncertainty = 0; iUncertainty < LongRangeSystematicOrganizer::knUncertaintySources; iUncertainty++){
    for(int iFlow = 0; iFlow < maxVn; iFlow++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        relativeUncertaintyTable[iUncertainty][iFlow][iCentrality] = 0;
        absoluteUncertaintyTable[iUncertainty][iFlow][iCentrality] = 0;
      }
    }
  }
  
  // Initialize the significance histograms to zero
  for(int iFlow = 0; iFlow < maxVn; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      graphName = Form("statisticallyInsignificant_v%d_C%d", iFlow+1, iCentrality);
      statisticallyInsignificant[iFlow][iCentrality] = new TH1D(graphName, graphName, LongRangeSystematicOrganizer::knUncertaintySources, -0.5, LongRangeSystematicOrganizer::knUncertaintySources-0.5);
    }
  }
  
  // Helper variables for reading points from graphs
  double nominalX, nominalY, nominalYerror;
  double comparisonX, comparisonY, comparisonYerror;
  double absoluteUncertainty;
  int maxErrorIndex;
  bool isInsignificant;
  
  // ========================================= //
  // Uncertainty coming from jet axis smearing //
  // ========================================= //
  
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      nominalResultGraph[iFlow]->GetPoint(iCentrality, nominalX, nominalY);
      nominalYerror = nominalResultGraph[iFlow]->GetErrorY(iCentrality);
      angleSmearResultGraph[iFlow]->GetPoint(iCentrality, comparisonX, comparisonY);
      comparisonYerror = angleSmearResultGraph[iFlow]->GetErrorY(iCentrality);
      
      absoluteUncertaintyTable[LongRangeSystematicOrganizer::kAngleSmear][iFlow][iCentrality] = TMath::Abs(nominalY-comparisonY);
      relativeUncertaintyTable[LongRangeSystematicOrganizer::kAngleSmear][iFlow][iCentrality] = TMath::Abs(1 - comparisonY/nominalY);
      
      // Check if nominal and smeared v2 are within statistical errors
      isInsignificant = isWithinStatisticalErrors(nominalY, comparisonY, nominalYerror, comparisonYerror);
      if(isInsignificant) statisticallyInsignificant[iFlow][iCentrality]->Fill(LongRangeSystematicOrganizer::kAngleSmear);
    }
  }
  
  // ================================================================================ //
  // Uncertainty coming from different deltaEta regions in jet-hadron correlation fit //
  // ================================================================================ //
  
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iDeltaEtaRegion = 0; iDeltaEtaRegion < 2; iDeltaEtaRegion++){
        
        nominalResultGraph[iFlow]->GetPoint(iCentrality, nominalX, nominalY);
        nominalYerror = nominalResultGraph[iFlow]->GetErrorY(iCentrality);
        deltaEtaRegionJetHadronGraph[iDeltaEtaRegion][iFlow]->GetPoint(iCentrality, comparisonX, comparisonY);
        comparisonYerror = deltaEtaRegionJetHadronGraph[iDeltaEtaRegion][iFlow]->GetErrorY(iCentrality);
        
        absoluteUncertainty = TMath::Abs(nominalY-comparisonY);
        
        if(absoluteUncertainty > absoluteUncertaintyTable[LongRangeSystematicOrganizer::kDeltaEtaRegion][iFlow][iCentrality]){
          absoluteUncertaintyTable[LongRangeSystematicOrganizer::kDeltaEtaRegion][iFlow][iCentrality] = absoluteUncertainty;
          relativeUncertaintyTable[LongRangeSystematicOrganizer::kDeltaEtaRegion][iFlow][iCentrality] = TMath::Abs(1 - comparisonY/nominalY);
          maxErrorIndex = iDeltaEtaRegion;
        }
        
      } // DeltaEta region loop
      
      deltaEtaRegionJetHadronGraph[maxErrorIndex][iFlow]->GetPoint(iCentrality, comparisonX, comparisonY);
      comparisonYerror = deltaEtaRegionJetHadronGraph[maxErrorIndex][iFlow]->GetErrorY(iCentrality);
      
      // Check if nominal and smaller deltaEta region v2 are within statistical errors
      isInsignificant = isWithinStatisticalErrors(nominalY, comparisonY, nominalYerror, comparisonYerror);
      if(isInsignificant) statisticallyInsignificant[iFlow][iCentrality]->Fill(LongRangeSystematicOrganizer::kDeltaEtaRegion);
      
    } // Centrality loop
  } // Flow component loop
  
  // ============================================================================== //
  // Uncertainty coming from different deltaEta sides in jet-hadron correlation fit //
  // ============================================================================== //
  
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iDeltaEtaRegion = 0; iDeltaEtaRegion < 2; iDeltaEtaRegion++){
        
        nominalResultGraph[iFlow]->GetPoint(iCentrality, nominalX, nominalY);
        nominalYerror = nominalResultGraph[iFlow]->GetErrorY(iCentrality);
        deltaEtaSideJetHadronGraph[iDeltaEtaRegion][iFlow]->GetPoint(iCentrality, comparisonX, comparisonY);
        comparisonYerror = deltaEtaSideJetHadronGraph[iDeltaEtaRegion][iFlow]->GetErrorY(iCentrality);
        
        absoluteUncertainty = TMath::Abs(nominalY-comparisonY);
        
        if(absoluteUncertainty > absoluteUncertaintyTable[LongRangeSystematicOrganizer::kDeltaEtaSide][iFlow][iCentrality]){
          absoluteUncertaintyTable[LongRangeSystematicOrganizer::kDeltaEtaSide][iFlow][iCentrality] = absoluteUncertainty;
          relativeUncertaintyTable[LongRangeSystematicOrganizer::kDeltaEtaSide][iFlow][iCentrality] = TMath::Abs(1 - comparisonY/nominalY);
          maxErrorIndex = iDeltaEtaRegion;
        }
        
      } // DeltaEta region loop
      
      deltaEtaSideJetHadronGraph[maxErrorIndex][iFlow]->GetPoint(iCentrality, comparisonX, comparisonY);
      comparisonYerror = deltaEtaSideJetHadronGraph[maxErrorIndex][iFlow]->GetErrorY(iCentrality);
      
      // Check if nominal and one deltaEta sided v2 are within statistical errors
      isInsignificant = isWithinStatisticalErrors(nominalY, comparisonY, nominalYerror, comparisonYerror);
      if(isInsignificant) statisticallyInsignificant[iFlow][iCentrality]->Fill(LongRangeSystematicOrganizer::kDeltaEtaSide);
      
    } // Centrality loop
  } // Flow component loop
  
  // ================================================ //
  // Uncertainty coming from jet collection selection //
  // ================================================ //
  
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      nominalResultGraph[iFlow]->GetPoint(iCentrality, nominalX, nominalY);
      nominalYerror = nominalResultGraph[iFlow]->GetErrorY(iCentrality);
      pfcsResultGraph[iFlow]->GetPoint(iCentrality, comparisonX, comparisonY);
      comparisonYerror = pfcsResultGraph[iFlow]->GetErrorY(iCentrality);
      
      absoluteUncertaintyTable[LongRangeSystematicOrganizer::kJetCollection][iFlow][iCentrality] = TMath::Abs(nominalY-comparisonY);
      relativeUncertaintyTable[LongRangeSystematicOrganizer::kJetCollection][iFlow][iCentrality] = TMath::Abs(1 - comparisonY/nominalY);
      
      // Check if nominal and pfcs v2 are within statistical errors
      isInsignificant = isWithinStatisticalErrors(nominalY, comparisonY, nominalYerror, comparisonYerror);
      if(isInsignificant) statisticallyInsignificant[iFlow][iCentrality]->Fill(LongRangeSystematicOrganizer::kJetCollection);
      
    } // Centrality loop
  } // Flow component loop
  
  // ======================================================================================= //
  // Uncertainty coming from dihadron correlation selection (dijet events / min bias events) //
  // ======================================================================================= //
  
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      nominalResultGraph[iFlow]->GetPoint(iCentrality, nominalX, nominalY);
      nominalYerror = nominalResultGraph[iFlow]->GetErrorY(iCentrality);
      minBiasDihadronResultGraph[iFlow]->GetPoint(iCentrality, comparisonX, comparisonY);
      comparisonYerror = minBiasDihadronResultGraph[iFlow]->GetErrorY(iCentrality);
      
      absoluteUncertaintyTable[LongRangeSystematicOrganizer::kMinBias][iFlow][iCentrality] = TMath::Abs(nominalY-comparisonY);
      relativeUncertaintyTable[LongRangeSystematicOrganizer::kMinBias][iFlow][iCentrality] = TMath::Abs(1 - comparisonY/nominalY);
      
      // Check if nominal and pfcs v2 are within statistical errors
      isInsignificant = isWithinStatisticalErrors(nominalY, comparisonY, nominalYerror, comparisonYerror);
      if(isInsignificant) statisticallyInsignificant[iFlow][iCentrality]->Fill(LongRangeSystematicOrganizer::kMinBias);
      
    } // Centrality loop
  } // Flow component loop
  
  // ========================================= //
  // Add all uncertainty sources in quadrature //
  // ========================================= //
  
  double sumOfSquares;
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      nominalResultGraph[iFlow]->GetPoint(iCentrality, nominalX, nominalY);
      
      sumOfSquares = 0;
      for(int iUncertainty = 0; iUncertainty < LongRangeSystematicOrganizer::kAll; iUncertainty++){
        sumOfSquares = sumOfSquares + TMath::Power(absoluteUncertaintyTable[iUncertainty][iFlow][iCentrality],2);
      }
      
      absoluteUncertaintyTable[LongRangeSystematicOrganizer::kAll][iFlow][iCentrality] = TMath::Sqrt(sumOfSquares);
      relativeUncertaintyTable[LongRangeSystematicOrganizer::kAll][iFlow][iCentrality] = TMath::Abs(absoluteUncertaintyTable[LongRangeSystematicOrganizer::kAll][iFlow][iCentrality]/nominalY);
      
    } // Centrality loop
  } // Flow component loop
  
  // ================================================ //
  // Collect the calculated uncertainties into graphs //
  // ================================================ //
  
  TGraphErrors* systematicUncertaintyGraph[LongRangeSystematicOrganizer::knUncertaintySources][maxVn];
  TString uncertaintyNameString;
  LongRangeSystematicOrganizer *uncertaintyNamer = new LongRangeSystematicOrganizer();
  
  for(int iUncertainty = 0; iUncertainty < LongRangeSystematicOrganizer::knUncertaintySources; iUncertainty++){
    uncertaintyNameString = uncertaintyNamer->GetLongRangeUncertaintyName(iUncertainty);
    for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
      graphName = Form("systematicUncertainty_v%d_%s", iFlow+1, uncertaintyNameString.Data());
      systematicUncertaintyGraph[iUncertainty][iFlow] = (TGraphErrors*) nominalResultGraph[iFlow]->Clone(graphName);
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        systematicUncertaintyGraph[iUncertainty][iFlow]->SetPointError(iCentrality, 0, absoluteUncertaintyTable[iUncertainty][iFlow][iCentrality]);
      } // Centrality loop
    } // Flow component loop
  } // Uncertainty source loop
  
  // =================================================== //
  // Write the systematic uncertainty graphs into a file //
  // =================================================== //
  
  // If an output file name is given, put the total uncertainties to an output file
  if(outputFileName.EndsWith(".root")){
    TFile *outputFile = new TFile(outputFileName,"UPDATE");

    for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
      for(int iUncertainty = 0; iUncertainty < LongRangeSystematicOrganizer::knUncertaintySources; iUncertainty++){
        uncertaintyNameString = uncertaintyNamer->GetLongRangeUncertaintyName(iUncertainty);
        graphName = Form("systematicUncertainty_v%d_%s", iFlow+1, uncertaintyNameString.Data());
        systematicUncertaintyGraph[iUncertainty][iFlow]->Write(graphName, TObject::kOverwrite);
      } // Uncertainty loop
    } // Flow component loop
    
    for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        graphName = Form("statisticallyInsignificant_v%d_C%d", iFlow+1, iCentrality);
        statisticallyInsignificant[iFlow][iCentrality]->Write(graphName, TObject::kOverwrite);
      } // Centrality loop
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
      
      for(int iUncertainty = 0; iUncertainty < LongRangeSystematicOrganizer::knUncertaintySources; iUncertainty++){
        cout << uncertaintyNamer->GetLongRangeUncertaintyName(iUncertainty).Data();
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          cout << " & $" << relativeUncertaintyTable[iUncertainty][iFlow][iCentrality] << "$";
        }
        cout << " \\\\" << endl;
      }
      
      cout << "    \\midrule" << endl;
      
      // Print the absolute uncertainties
      
      for(int iUncertainty = 0; iUncertainty < LongRangeSystematicOrganizer::knUncertaintySources; iUncertainty++){
        cout << uncertaintyNamer->GetLongRangeUncertaintyName(iUncertainty).Data();
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
