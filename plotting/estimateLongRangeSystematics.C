#include "LongRangeSystematicOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

// Function definitions. Implementations after the main macro.
std::tuple<double,double,bool> findTheDifference(TGraphErrors *nominalGraph, TGraphErrors *comparisonGraph[], const int nComparisonGraphs, const int iPoint);
std::tuple<double,double,bool> findTheDifference(TGraphErrors *nominalGraph, TGraphErrors *comparisonGraph, const int iPoint);
void drawIllustratingPlots(JDrawer *drawer, TGraphErrors *nominalGraph, TGraphErrors *comparisonGraph[], const int nComparisonGraphs, TString comparisonLegend[], const int iFlow, TString plotName);
void drawIllustratingPlots(JDrawer *drawer, TGraphErrors *nominalGraph, TGraphErrors *comparisonGraph, TString comparisonLegend, const int iFlow, TString plotName);

/*
 * Macro for estimating systematic uncertainties for long range correlation results
 *
 * The following sources of systematic uncertainty are estimated:
 *  - Adjustment on the gluing point of leading/subleading distributions
 *  - Fit projections over regions -2.5 < deltaEta < -1.5 and 1.5 < deltaEta < 2.5 and compare to nominal
 *  - Fit projections over regions 1.5 < |deltaEta| < 2 and 2 < |deltaEta| < 2.5 and compare to nominal
 *  - Fit different v_{z} selection (like positive and negative)
 *  - Differences coming from selecting different jet collections
 *  - Jet axis smearing
 *  - Comparing dihadron correlations from dijet and minimum bias events
 *
 *  The results will be saved to a file and also slides with different contributions can be printed
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
  
  // Results from different vz regions
  TString negativeVzResultFileName = "flowGraphs/summaryPlot_akCaloJet_negativeVzJetHadron_nominalCorrection_2021-05-27.root";
  TString positiveVzResultFileName = "flowGraphs/summaryPlot_akCaloJet_positiveVzJetHadron_nominalCorrection_2021-05-27.root";
  
  // Results with minimum bias dihadron
  TString minBiasDihadronResultFileName = "flowGraphs/summaryPlot_akCaloJet_minBiasDihadron_nominalCorrection_2021-05-26.root";
  
  // Results with final correction TODO: Check this is the actual nominal correction.
  TString finalResultFileName = "flowGraphs/summaryPlot_akCaloJet_matchHadronV2ScaleYieldWith4pCentShift_2021-06-03.root";
  TString shiftedResultFileName1 = "flowGraphs/summaryPlot_akCaloJet_matchHadronV2ScaleYieldWith3pCentShift_2021-06-03.root";
  TString shiftedResultFileName2 = "flowGraphs/summaryPlot_akCaloJet_matchHadronV2ScaleYieldWith5pCentShift_2021-06-03.root";
  
  // Results without tracking efficiency correction
  TString noTrackEfficiencyFileName = "flowGraphs/summaryPlot_akCaloJet_noTrackEfficiency_2021-07-14.root";
  
  // Results with varied quark/gluon fraction
  TString quarkGluonFractionFileName = "flowGraphs/summaryPlot_akCaloJet_correctionWith25pMoreQuarkJets_2021-07-26.root";
  
  // Results where jet energy correction is altered by its uncertainties
  TString jecMinusFileName = "flowGraphs/summaryPlot_akCaloJet_JECminus_2021-08-05.root";
  TString jecPlusFileName = "flowGraphs/summaryPlot_akCaloJet_JECplus_2021-08-05.root";
  
  // Flow component configuration
  const int maxVn = 4;      // Maximum number of flow components
  const int firstFlow = 2;  // First analyzed flow component
  const int lastFlow = 3;   // Last analyzed flow component
  
  const bool printUncertainties = false;
  
  bool plotExample = false;
  
  TString outputFileName = "flowGraphs/systematicUncertainties_includeJEC_finalCorrection_2021-08-05.root";
  
  // ==================================================================
  // ====================== Configuration done ========================
  // ==================================================================
  
  // Long range systematic organizer can easily provide names for all your naming purposes
  LongRangeSystematicOrganizer *nameGiver = new LongRangeSystematicOrganizer();
  
  // Open the input files
  TFile* nominalResultFile = TFile::Open(nominalResultFileName);
  TFile* pfcsResultFile = TFile::Open(pfcsResultFileName);
  TFile* angleSmearResultFile = TFile::Open(angleSmearResultFileName);
  TFile* narrowDeltaEtaJetHadronResultFile = TFile::Open(narrowDeltaEtaJetHadronResultFileName);
  TFile* wideDeltaEtaJetHadronResultFile = TFile::Open(wideDeltaEtaJetHadronResultFileName);
  TFile* negativeDeltaEtaJetHadronResultFile = TFile::Open(negativeDeltaEtaJetHadronResultFileName);
  TFile* positiveDeltaEtaJetHadronResultFile = TFile::Open(positiveDeltaEtaJetHadronResultFileName);
  TFile* negativeVzResultFile = TFile::Open(negativeVzResultFileName);
  TFile* positiveVzResultFile = TFile::Open(positiveVzResultFileName);
  TFile* minBiasDihadronResultFile = TFile::Open(minBiasDihadronResultFileName);
  TFile* finalResultFile = TFile::Open(finalResultFileName);
  TFile* shiftedResultFile1 = TFile::Open(shiftedResultFileName1);
  TFile* shiftedResultFile2 = TFile::Open(shiftedResultFileName2);
  TFile* noTrackEfficiencyFile = TFile::Open(noTrackEfficiencyFileName);
  TFile* quarkGluonFractionFile = TFile::Open(quarkGluonFractionFileName);
  TFile* jecFile1 = TFile::Open(jecMinusFileName);
  TFile* jecFile2 = TFile::Open(jecPlusFileName);
  
  // Result graphs
  TGraphErrors* nominalResultGraph[maxVn];
  TGraphErrors* pfcsResultGraph[maxVn];
  TGraphErrors* angleSmearResultGraph[maxVn];
  TGraphErrors* deltaEtaRegionJetHadronGraph[maxVn][2];
  TGraphErrors* deltaEtaSideJetHadronGraph[maxVn][2];
  TGraphErrors* vzRegionJetHadronGraph[maxVn][2];
  TGraphErrors* minBiasDihadronResultGraph[maxVn];
  TGraphErrors* finalResultGraph[maxVn];
  TGraphErrors* shiftedResultGraph[maxVn][2];
  TGraphErrors* noTrackEfficiencyGraph[maxVn];
  TGraphErrors* quarkGluonFractionGraph[maxVn];
  TGraphErrors* jecGraph[maxVn][2];
  
  TString graphName;
  
  // Load the graphs from the files
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    
    graphName = Form("summaryV%d", iFlow+1);
    
    nominalResultGraph[iFlow] = (TGraphErrors*) nominalResultFile->Get(graphName);
    pfcsResultGraph[iFlow] = (TGraphErrors*) pfcsResultFile->Get(graphName);
    angleSmearResultGraph[iFlow] = (TGraphErrors*) angleSmearResultFile->Get(graphName);
    deltaEtaRegionJetHadronGraph[iFlow][0] = (TGraphErrors*) narrowDeltaEtaJetHadronResultFile->Get(graphName);
    deltaEtaRegionJetHadronGraph[iFlow][1] = (TGraphErrors*) wideDeltaEtaJetHadronResultFile->Get(graphName);
    deltaEtaSideJetHadronGraph[iFlow][0] = (TGraphErrors*) negativeDeltaEtaJetHadronResultFile->Get(graphName);
    deltaEtaSideJetHadronGraph[iFlow][1] = (TGraphErrors*) positiveDeltaEtaJetHadronResultFile->Get(graphName);
    vzRegionJetHadronGraph[iFlow][0] = (TGraphErrors*) negativeVzResultFile->Get(graphName);
    vzRegionJetHadronGraph[iFlow][1] = (TGraphErrors*) positiveVzResultFile->Get(graphName);
    minBiasDihadronResultGraph[iFlow] = (TGraphErrors*) minBiasDihadronResultFile->Get(graphName);
    finalResultGraph[iFlow] = (TGraphErrors*) finalResultFile->Get(graphName);
    shiftedResultGraph[iFlow][0] = (TGraphErrors*) shiftedResultFile1->Get(graphName);
    shiftedResultGraph[iFlow][1] = (TGraphErrors*) shiftedResultFile2->Get(graphName);
    noTrackEfficiencyGraph[iFlow] = (TGraphErrors*) noTrackEfficiencyFile->Get(graphName);
    quarkGluonFractionGraph[iFlow] = (TGraphErrors*) quarkGluonFractionFile->Get(graphName);
    jecGraph[iFlow][0] = (TGraphErrors*) jecFile1->Get(graphName);
    jecGraph[iFlow][1] = (TGraphErrors*) jecFile2->Get(graphName);
    
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
  
  // Helper variables for reading points from the graphs
  double nominalX, nominalY;
  
  // Helper variables for uncertainties
  double absoluteUncertainty, relativeUncertainty;
  bool isInsignificant;
  
  // Drawer for drawing graphs
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceGraph();
  drawer->SetNDivisionsX(510);
  drawer->SetBottomMargin(0.18);
  drawer->SetTitleOffsetX(1.63);
  drawer->SetLabelOffsetX(0.04);
  drawer->SetTitleOffsetY(1.6);
  
  // Legends given to drawns graphs
  TString legendNames[2];
  
  // ========================================= //
  // Uncertainty coming from jet axis smearing //
  // ========================================= //
  
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      std::tie(absoluteUncertainty, relativeUncertainty, isInsignificant) = findTheDifference(nominalResultGraph[iFlow], angleSmearResultGraph[iFlow], iCentrality);
      
      absoluteUncertaintyTable[LongRangeSystematicOrganizer::kAngleSmear][iFlow][iCentrality] = absoluteUncertainty;
      relativeUncertaintyTable[LongRangeSystematicOrganizer::kAngleSmear][iFlow][iCentrality] = relativeUncertainty;
      if(isInsignificant) statisticallyInsignificant[iFlow][iCentrality]->Fill(LongRangeSystematicOrganizer::kAngleSmear);
      
    } // Centrality loop
    
    // Draw example plots on how the uncertainty is obtained
    if(plotExample){
      legendNames[0] = "Smeared axis";
      drawIllustratingPlots(drawer, nominalResultGraph[iFlow], angleSmearResultGraph[iFlow], legendNames[0], iFlow, nameGiver->GetLongRangeUncertaintyName(LongRangeSystematicOrganizer::kAngleSmear));
    }
    
  } // Flow component loop
  
  // ================================================================================ //
  // Uncertainty coming from different deltaEta regions in jet-hadron correlation fit //
  // ================================================================================ //
  
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      std::tie(absoluteUncertainty, relativeUncertainty, isInsignificant) = findTheDifference(nominalResultGraph[iFlow], deltaEtaRegionJetHadronGraph[iFlow], 2, iCentrality);
      
      absoluteUncertaintyTable[LongRangeSystematicOrganizer::kDeltaEtaRegion][iFlow][iCentrality] = absoluteUncertainty;
      relativeUncertaintyTable[LongRangeSystematicOrganizer::kDeltaEtaRegion][iFlow][iCentrality] = relativeUncertainty;
      if(isInsignificant) statisticallyInsignificant[iFlow][iCentrality]->Fill(LongRangeSystematicOrganizer::kDeltaEtaRegion);
      
    } // Centrality loop
    
    // Draw example plots on how the uncertainty is obtained
    if(plotExample){
      legendNames[0] = "1.5 < |#Delta#eta| < 2.0";
      legendNames[1] = "2.0 < |#Delta#eta| < 2.5";
      drawIllustratingPlots(drawer, nominalResultGraph[iFlow], deltaEtaRegionJetHadronGraph[iFlow], 2, legendNames, iFlow, nameGiver->GetLongRangeUncertaintyName(LongRangeSystematicOrganizer::kDeltaEtaRegion));
    }
    
  } // Flow component loop
  
  // ============================================================================== //
  // Uncertainty coming from different deltaEta sides in jet-hadron correlation fit //
  // ============================================================================== //
  
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      std::tie(absoluteUncertainty, relativeUncertainty, isInsignificant) = findTheDifference(nominalResultGraph[iFlow], deltaEtaSideJetHadronGraph[iFlow], 2, iCentrality);
      
      absoluteUncertaintyTable[LongRangeSystematicOrganizer::kDeltaEtaSide][iFlow][iCentrality] = absoluteUncertainty;
      relativeUncertaintyTable[LongRangeSystematicOrganizer::kDeltaEtaSide][iFlow][iCentrality] = relativeUncertainty;
      if(isInsignificant) statisticallyInsignificant[iFlow][iCentrality]->Fill(LongRangeSystematicOrganizer::kDeltaEtaSide);

      
    } // Centrality loop
    
    // Draw example plots on how the uncertainty is obtained
    if(plotExample){
      legendNames[0] = "-2.5 < #Delta#eta < -1.5";
      legendNames[1] = "1.5 < #Delta#eta < 2.5";
      drawIllustratingPlots(drawer, nominalResultGraph[iFlow], deltaEtaSideJetHadronGraph[iFlow], 2, legendNames, iFlow, nameGiver->GetLongRangeUncertaintyName(LongRangeSystematicOrganizer::kDeltaEtaSide));
    }
    
  } // Flow component loop
  
  // ======================================================================== //
  // Uncertainty coming from different vz sides in jet-hadron correlation fit //
  // ======================================================================== //
  
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      std::tie(absoluteUncertainty, relativeUncertainty, isInsignificant) = findTheDifference(nominalResultGraph[iFlow], vzRegionJetHadronGraph[iFlow], 2, iCentrality);
      
      absoluteUncertaintyTable[LongRangeSystematicOrganizer::kVz][iFlow][iCentrality] = absoluteUncertainty;
      relativeUncertaintyTable[LongRangeSystematicOrganizer::kVz][iFlow][iCentrality] = relativeUncertainty;
      if(isInsignificant) statisticallyInsignificant[iFlow][iCentrality]->Fill(LongRangeSystematicOrganizer::kVz);
      
    } // Centrality loop
    
    // Draw example plots on how the uncertainty is obtained
    if(plotExample){
      legendNames[0] = "v_{z} < 0";
      legendNames[1] = "v_{z} > 0";
      drawIllustratingPlots(drawer, nominalResultGraph[iFlow], vzRegionJetHadronGraph[iFlow], 2, legendNames, iFlow, nameGiver->GetLongRangeUncertaintyName(LongRangeSystematicOrganizer::kVz));
    }
    
  } // Flow component loop
  
  // ================================================ //
  // Uncertainty coming from jet collection selection //
  // ================================================ //
  
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      std::tie(absoluteUncertainty, relativeUncertainty, isInsignificant) = findTheDifference(nominalResultGraph[iFlow], pfcsResultGraph[iFlow], iCentrality);
      
      absoluteUncertaintyTable[LongRangeSystematicOrganizer::kJetCollection][iFlow][iCentrality] = absoluteUncertainty;
      relativeUncertaintyTable[LongRangeSystematicOrganizer::kJetCollection][iFlow][iCentrality] = relativeUncertainty;
      if(isInsignificant) statisticallyInsignificant[iFlow][iCentrality]->Fill(LongRangeSystematicOrganizer::kJetCollection);
      
    } // Centrality loop
    
    // Draw example plots on how the uncertainty is obtained
    if(plotExample){
      legendNames[0] = "akCs4PFJet";
      drawIllustratingPlots(drawer, nominalResultGraph[iFlow], pfcsResultGraph[iFlow], legendNames[0], iFlow, nameGiver->GetLongRangeUncertaintyName(LongRangeSystematicOrganizer::kJetCollection));
    }
    
  } // Flow component loop
  
  // ========================================== //
  // Uncertainty coming from Monte Carlo tuning //
  // ========================================== //
  
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      std::tie(absoluteUncertainty, relativeUncertainty, isInsignificant) = findTheDifference(finalResultGraph[iFlow], shiftedResultGraph[iFlow], 2, iCentrality);
      
      absoluteUncertaintyTable[LongRangeSystematicOrganizer::kMCTuning][iFlow][iCentrality] = absoluteUncertainty;
      relativeUncertaintyTable[LongRangeSystematicOrganizer::kMCTuning][iFlow][iCentrality] = relativeUncertainty;
      if(isInsignificant) statisticallyInsignificant[iFlow][iCentrality]->Fill(LongRangeSystematicOrganizer::kMCTuning);
      
    } // Centrality loop
    
    // Draw example plots on how the uncertainty is obtained
    if(plotExample){
      legendNames[0] = "MC + 3%";
      legendNames[1] = "MC + 5%";
      drawIllustratingPlots(drawer, finalResultGraph[iFlow], shiftedResultGraph[iFlow], 2, legendNames, iFlow, nameGiver->GetLongRangeUncertaintyName(LongRangeSystematicOrganizer::kMCTuning));
    }
    
  } // Flow component loop
  
  // ======================================================================================= //
  // Uncertainty coming from dihadron correlation selection (dijet events / min bias events) //
  // ======================================================================================= //
  
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      std::tie(absoluteUncertainty, relativeUncertainty, isInsignificant) = findTheDifference(nominalResultGraph[iFlow], minBiasDihadronResultGraph[iFlow], iCentrality);
      
      absoluteUncertaintyTable[LongRangeSystematicOrganizer::kMinBias][iFlow][iCentrality] = absoluteUncertainty;
      relativeUncertaintyTable[LongRangeSystematicOrganizer::kMinBias][iFlow][iCentrality] = relativeUncertainty;
      if(isInsignificant) statisticallyInsignificant[iFlow][iCentrality]->Fill(LongRangeSystematicOrganizer::kMinBias);
      
    } // Centrality loop
    
    // Draw example plots on how the uncertainty is obtained
    if(plotExample){
      legendNames[0] = "Minimum bias";
      drawIllustratingPlots(drawer, nominalResultGraph[iFlow], minBiasDihadronResultGraph[iFlow], legendNames[0], iFlow, nameGiver->GetLongRangeUncertaintyName(LongRangeSystematicOrganizer::kMinBias));
    }
    
  } // Flow component loop
  
  // =========================================== //
  // Uncertainty coming tracking related sources //
  // =========================================== //
  
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      std::tie(absoluteUncertainty, relativeUncertainty, isInsignificant) = findTheDifference(finalResultGraph[iFlow], noTrackEfficiencyGraph[iFlow], iCentrality);
      
      absoluteUncertaintyTable[LongRangeSystematicOrganizer::kTracking][iFlow][iCentrality] = absoluteUncertainty;
      relativeUncertaintyTable[LongRangeSystematicOrganizer::kTracking][iFlow][iCentrality] = relativeUncertainty;
      if(isInsignificant) statisticallyInsignificant[iFlow][iCentrality]->Fill(LongRangeSystematicOrganizer::kTracking);
      
    } // Centrality loop
    
    // Draw example plots on how the uncertainty is obtained
    if(plotExample){
      legendNames[0] = "No track efficiency";
      drawIllustratingPlots(drawer, finalResultGraph[iFlow], noTrackEfficiencyGraph[iFlow], legendNames[0], iFlow, nameGiver->GetLongRangeUncertaintyName(LongRangeSystematicOrganizer::kTracking));
    }
    
  } // Flow component loop
  
  // ===================================================================== //
  // Uncertainty coming uncertainty on quark/gluon fraction in Monte Carlo //
  // ===================================================================== //
  
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      std::tie(absoluteUncertainty, relativeUncertainty, isInsignificant) = findTheDifference(finalResultGraph[iFlow], quarkGluonFractionGraph[iFlow], iCentrality);
      
      absoluteUncertaintyTable[LongRangeSystematicOrganizer::kQuarkGluonFraction][iFlow][iCentrality] = absoluteUncertainty;
      relativeUncertaintyTable[LongRangeSystematicOrganizer::kQuarkGluonFraction][iFlow][iCentrality] = relativeUncertainty;
      if(isInsignificant) statisticallyInsignificant[iFlow][iCentrality]->Fill(LongRangeSystematicOrganizer::kQuarkGluonFraction);
      
    } // Centrality loop
    
    // Draw example plots on how the uncertainty is obtained
    if(plotExample){
      legendNames[0] = "Adjusted q/g fraction";
      drawIllustratingPlots(drawer, finalResultGraph[iFlow], quarkGluonFractionGraph[iFlow], legendNames[0], iFlow, nameGiver->GetLongRangeUncertaintyName(LongRangeSystematicOrganizer::kQuarkGluonFraction));
    }
    
  } // Flow component loop
  
  // ============================================= //
  // Uncertainty coming from jet energy correction //
  // ============================================= //
  
  for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      std::tie(absoluteUncertainty, relativeUncertainty, isInsignificant) = findTheDifference(finalResultGraph[iFlow], jecGraph[iFlow], 2, iCentrality);
      
      absoluteUncertaintyTable[LongRangeSystematicOrganizer::kJEC][iFlow][iCentrality] = absoluteUncertainty;
      relativeUncertaintyTable[LongRangeSystematicOrganizer::kJEC][iFlow][iCentrality] = relativeUncertainty;
      if(isInsignificant) statisticallyInsignificant[iFlow][iCentrality]->Fill(LongRangeSystematicOrganizer::kJEC);
      
    } // Centrality loop
    
    // Draw example plots on how the uncertainty is obtained
    if(plotExample){
      legendNames[0] = "JEC minus";
      legendNames[1] = "JEC plus";
      drawIllustratingPlots(drawer, finalResultGraph[iFlow], jecGraph[iFlow], 2, legendNames, iFlow, nameGiver->GetLongRangeUncertaintyName(LongRangeSystematicOrganizer::kJEC));
    }
    
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
      //systematicUncertaintyGraph[iUncertainty][iFlow] = (TGraphErrors*) nominalResultGraph[iFlow]->Clone(graphName);
      systematicUncertaintyGraph[iUncertainty][iFlow] = (TGraphErrors*) finalResultGraph[iFlow]->Clone(graphName);
      
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

/*
 * Function for checking if two values are within statistical uncertainties
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
 * Function for finding the realtive and absolute difference of points in two graphs
 *
 *  TGraphErrors *nominalGraph = Graph containing nominal values
 *  TGraphErrors *comparisonGraph[] = Array of graphs containing comparison values
 *  const int nComparisonGraphs = Number of comparison graphs in the array
 *  const int iPoint = Index of the point that is compared
 *
 *  return: Tuple containing the absolute uncertainty, the relative uncertainty and information if the difference in within statistical uncertainites
 *
 */
std::tuple<double,double,bool> findTheDifference(TGraphErrors *nominalGraph, TGraphErrors *comparisonGraph[], const int nComparisonGraphs, const int iPoint){
  
  // Helper variables for reading points from graphs
  double nominalX, nominalY, nominalYerror;
  double comparisonX, comparisonY, comparisonYerror;
  double currentUncertainty, relativeUncertainty;
  double absoluteUncertainty = 0;
  int maxErrorIndex;
  bool isInsignificant;
  
  for(int iComparison = 0; iComparison < nComparisonGraphs; iComparison++){
    
    nominalGraph->GetPoint(iPoint, nominalX, nominalY);
    nominalYerror = nominalGraph->GetErrorY(iPoint);
    comparisonGraph[iComparison]->GetPoint(iPoint, comparisonX, comparisonY);
    comparisonYerror = comparisonGraph[iComparison]->GetErrorY(iPoint);
    
    currentUncertainty = TMath::Abs(nominalY-comparisonY);
    
    if(currentUncertainty > absoluteUncertainty){
      absoluteUncertainty = currentUncertainty;
      relativeUncertainty = TMath::Abs(1 - comparisonY/nominalY);
      maxErrorIndex = iComparison;
    }
    
  } // Comparison graph loop
  
  comparisonGraph[maxErrorIndex]->GetPoint(iPoint, comparisonX, comparisonY);
  comparisonYerror = comparisonGraph[maxErrorIndex]->GetErrorY(iPoint);
  
  // Check if nominal and comparisonv2 are within statistical errors
  isInsignificant = isWithinStatisticalErrors(nominalY, comparisonY, nominalYerror, comparisonYerror);
  
  // Return the uncertainties and the information whether the numbers are within statistical uncertainties
  return std::make_tuple(absoluteUncertainty,relativeUncertainty,isInsignificant);
  
}

/*
 * Function for finding the realtive and absolute difference of points in two graphs
 *
 *  TGraphErrors *nominalGraph = Graph containing nominal values
 *  TGraphErrors *comparisonGraph = Graph containing comparison values
 *  const int iPoint = Index of the point that is compared
 *
 *  return: Tuple containing the absolute uncertainty, the relative uncertainty and information if the difference in within statistical uncertainites
 *
 */
std::tuple<double,double,bool> findTheDifference(TGraphErrors *nominalGraph, TGraphErrors *comparisonGraph, const int iPoint){
  
  // Enclose the single graph to an array and use the difference finder for graph array
  TGraphErrors *comparisonArray[1] = {comparisonGraph};
  return findTheDifference(nominalGraph, comparisonArray, 1, iPoint);
  
}

/*
 * Draw plots illustrating how systematic uncertainties are estimated form a certain source
 *
 *  JDrawer *drawer = JDrawer doing the dirty work in drawing
 *  TGraphErrors *nominalGraph = Graph containing nominal results
 *  TGraphErrors *comparisonGraph[] = Array of graphs containing the graphs compared with the nominal
 *  const int nComparisonGraphs = Number of comparison graphs in the array
 *  TString comparisonLegend[] = An array of strings describing each comparison graphs
 *  const int iFlow = Index of the plotted flow component
 *  TString plotName = String added to saved plots. If left empty, the plots are not saved into files.
 */
void drawIllustratingPlots(JDrawer *drawer, TGraphErrors *nominalGraph, TGraphErrors *comparisonGraph[], const int nComparisonGraphs, TString comparisonLegend[], const int iFlow, TString plotName){
  
  TLegend *legend = new TLegend(0.23,0.84-nComparisonGraphs*0.06,0.53,0.9);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  const int markers[] = {kFullDiamond, kFullDoubleDiamond, kFullCross, kFullFourTrianglesPlus, kFullSquare, kFullStar};
  const int colors[] = { kBlue, kRed, kGreen+3, kMagenta, kCyan, kBlack, kViolet};
  
  // Zooming ranges for different flow components
  double minZoom[] = {0,0,-0.03,0};
  double maxZoom[] = {0.1,0.1,0.03,0.1};
  
  TLine *zeroLine = new TLine(0.75,0,3.25,0);
  zeroLine->SetLineStyle(2);
  
  // Draw the nominal graph to the canvas
  nominalGraph->SetMarkerStyle(kFullCircle);
  nominalGraph->SetMarkerSize(1.3);
  nominalGraph->SetMarkerColor(kBlack);
  drawer->DrawGraphCustomAxes(nominalGraph, 0, 4, minZoom[iFlow], maxZoom[iFlow], "Centrality", Form("Jet v_{%d}", iFlow+1), " ", "ap");
  legend->AddEntry(nominalGraph, "Nominal result", "p");
  
  // Draw the comparison graphs to the same
  for(int iComparison = 0; iComparison < nComparisonGraphs; iComparison++){
    comparisonGraph[iComparison]->SetMarkerStyle(markers[iComparison]);
    comparisonGraph[iComparison]->SetMarkerSize(1.3);
    comparisonGraph[iComparison]->SetMarkerColor(colors[iComparison]);
    comparisonGraph[iComparison]->Draw("p,same");
    legend->AddEntry(comparisonGraph[iComparison], comparisonLegend[iComparison], "p");
  }
  
  if(iFlow == 2) zeroLine->Draw();
  legend->Draw();
  
  // If a plot name is given, save the plot in a file
  if(plotName != ""){
    gPad->GetCanvas()->SaveAs(Form("figures/systematicUncertaintyJetV%d_%s.pdf", iFlow+1, plotName.Data()));
  }
  
  // Calculate ratios between nominal and comparison graphs
  TGraphErrors *ratioGraph[nComparisonGraphs];
  double xPoint1, yPoint1, xPoint2, yPoint2, yError1, yError2, ratioValue, combinedError;
  TString binLabels[] = {"0-10%"," ","10-30%"," ","30-50%"," ","50-90%"};
  
  for(int iComparison = 0; iComparison < nComparisonGraphs; iComparison++){
    ratioGraph[iComparison] = (TGraphErrors*) comparisonGraph[iComparison]->Clone(Form("ratio_%s", comparisonLegend[iComparison].Data()));
    
    for(int iPoint = 0; iPoint < nominalGraph->GetN(); iPoint++){
      
      nominalGraph->GetPoint(iPoint, xPoint1, yPoint1);
      yError1 = nominalGraph->GetErrorY(iPoint);
      comparisonGraph[iComparison]->GetPoint(iPoint, xPoint2, yPoint2);
      yError2 = comparisonGraph[iComparison]->GetErrorY(iPoint);
      
      ratioValue = yPoint2 / yPoint1;
      combinedError = TMath::Sqrt(TMath::Power(yError2/yPoint1,2) + TMath::Power((yPoint2*yError1)/(yPoint1*yPoint1),2));
      ratioGraph[iComparison]->SetPoint(iPoint,xPoint1,ratioValue);
      ratioGraph[iComparison]->SetPointError(iPoint, 0, combinedError);
      
    } // Point loop
    
    // Set the bin labels for x-axis
    for(int iPoint = 0; iPoint < nominalGraph->GetN()*2; iPoint++){
      ratioGraph[iComparison]->GetXaxis()->ChangeLabel(iPoint+1,-1,-1,-1,-1,-1,binLabels[iPoint]);
    }
    
  } // Comparison graph loop
    
  // Draw the ratio plots
  TLine *oneLine = new TLine(0.75,1,3.25,1);
  oneLine->SetLineStyle(2);
  
  legend = new TLegend(0.25,0.9-nComparisonGraphs*0.06,0.55,0.9);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  
  // Zooming options for different flow components
  double minZoomRatio[] = {0,0.65,0,0};
  double maxZoomRatio[] = {2,1.35,2,2};
  
  for(int iComparison = 0; iComparison < nComparisonGraphs; iComparison++){
    if(iComparison == 0){
      drawer->DrawGraphCustomAxes(ratioGraph[iComparison], 0, 4, minZoomRatio[iFlow], maxZoomRatio[iFlow], "Centrality", Form("Jet v_{%d} ratio", iFlow+1), " ", "ap");
    } else {
      ratioGraph[iComparison]->Draw("p,same");
    }
    legend->AddEntry(ratioGraph[iComparison], Form("%s / Nominal",comparisonLegend[iComparison].Data()), "p");
  }
  
  oneLine->Draw();
  legend->Draw();
  
  // If a plot name is given, save the plot in a file
  if(plotName != ""){
    gPad->GetCanvas()->SaveAs(Form("figures/systematicUncertaintyJetV%dRatio_%s.pdf", iFlow+1, plotName.Data()));
  }
}

/*
 * Draw plots illustrating how systematic uncertainties are estimated form a certain source
 *
 *  JDrawer *drawer = JDrawer doing the dirty work in drawing
 *  TGraphErrors *nominalGraph = Graph containing nominal results
 *  TGraphErrors *comparisonGraph = The graphs compared with the nominal one
 *  const int nComparisonGraphs = Number of comparison graphs in the array
 *  TString comparisonLegend = A strings describing the comparison graph
 *  const int iFlow = Index of the plotted flow component
 *  TString plotName = String added to saved plots. If left empty, the plots are not saved into files.
 */
void drawIllustratingPlots(JDrawer *drawer, TGraphErrors *nominalGraph, TGraphErrors *comparisonGraph, TString comparisonLegend, const int iFlow, TString plotName){
  TGraphErrors *comparisonGraphArray[1] = {comparisonGraph};
  TString comparisonLegendArray[1] = {comparisonLegend};
  drawIllustratingPlots(drawer, nominalGraph, comparisonGraphArray, 1, comparisonLegendArray, iFlow, plotName);
}
