#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "DijetMethods.h"
#include "JDrawer.h"
#include "JffCorrector.h"

/*
 * Find the point where pol0 and pol1 lines cross
 */
double findCrossingPoint(TF1* poly1, double constant){
  
  // Find the parameters of the first order polynomial
  double level = poly1->GetParameter(0);
  double slope = poly1->GetParameter(1);
  
  // We need to solve equation x * slope + level = consant => x = (constant - level) / slope
  return (constant - level) / slope;
  
}

/*
 * Find the point where pol0 and pol1 lines cross
 */
double findCrossingPoint(TF1* poly1, TLine *poly0){
  
  // Find the constant value of the zeroth order polynomial
  double constant = poly0->GetY1();
  
  // Uee the function for constant to get the final results
  return findCrossingPoint(poly1, constant);
  
}

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
  const int maxFiles = 7;
  TString directoryName = "flowGraphs/";
  TString graphFileName = "flowGraphs_PbPb2018_akPfCsJets_correctedJetHadron_correctedDihadronSmallStats_2021-05-13.root";
  // flowGraphs_PbPb2018_caloJets_fixedJEC_correctedJetHadron_correctedDihadron_2021-02-26.root
  // flowGraphs_PbPb2018_caloJets_cumulativePtBins_correctedJetHadron_correctedEventDihadron_2021-03-22.root
  // flowGraphs_PbPb2018_akPfCsJets_correctedJetHadron_correctedDihadronSmallStats_2021-05-13.root
  // flowGraphs_PbPb2018_caloJets_angleSmear_correctedJetHadron_correctedDihadron_2021-02-26.root
  // flowGraphs_PbPb2018_caloJets_noTrackEff_correctedJetHadron_correctedDihadron_2021-07-13.root
  // flowGraphs_PbPb2018_caloJets_jetHadronJECminus_correctedJetHadron_correctedDihadron_2021-08-05.root
  // flowGraphs_PbPb2018_caloJets_jetHadronJECplus_correctedJetHadron_correctedDihadron_2021-08-05.root
  // flowGraphs_PbPb2018_caloJets_smearedJER_lowStatsDihadron_correctedJetHadron_correctedDihadron_2021-08-16.root
  // flowGraphs_PbPb2018_caloJets_smearedJERjetHadron_correctedJetHadron_correctedDihadron_2021-08-09.root
  // flowGraphs_PbPb2018_caloJets_fixedJEC_refitLongRange_narrowDeltaEta_correctedJetHadron_correctedDihadron_2021-05-26.root
  // flowGraphs_PbPb2018_caloJets_fixedJEC_refitLongRange_wideDeltaEta_correctedJetHadron_correctedDihadron_2021-05-26.root
  // flowGraphs_PbPb2018_caloJets_fixedJEC_refitLongRange_negativeDeltaEta_correctedJetHadron_correctedDihadron_2021-05-26.root
  // flowGraphs_PbPb2018_caloJets_fixedJEC_refitLongRange_positiveDeltaEta_correctedJetHadron_correctedDihadron_2021-05-26.root
  // flowGraphs_PbPb2018_caloJets_fixedJEC_negativeVz_correctedJetHadron_correctedDihadron_2021-05-27.root
  // flowGraphs_PbPb2018_caloJets_fixedJEC_positiveVz_correctedJetHadron_correctedDihadron_2021-05-27.root
  // flowGraphs_PbPb2018_caloJets_fixedJEC_minBiasDihadron_correctedJetHadron_correctedDihadron_2021-05-26.root
  // flowGraphs_PbPb2018_caloJets_manualJECforJetHadron_correctedJetHadron_correctedDihadron_2021-08-30.root
  
  TFile *graphFile[maxFiles];
  graphFile[0] = TFile::Open(directoryName+graphFileName);
  
  // Other files whose results can be compared with the nominal file
  int nComparisonFiles = 1;
  TString comparisonFileName[] = {  "flowGraphs_PbPb2018_caloJets_fixedJEC_correctedJetHadronWithoutTwoCuts_correctedDihadron_2021-08-30.root", "flowGraphs_PbPbMC2018_pfCsJets_4pCentShift_subeNon0_manualJECconeReflect_correctedJetHadron_correctedDihadron_2021-08-13.root", "flowGraphs_PbPbMC2018_pfCsJets_4pCentShift_subeNon0_manualJECconeReflectNeutralScaled_correctedJetHadron_correctedDihadron_2021-08-23.root", "flowGraphs_PbPbMC2018_pfCsJets_4pCentShift_subeNon0_manualJECconeReflectDijet_correctedJetHadron_correctedDihadron_2021-08-26.root", "flowGraphs_PbPbMC2018_subeNon0_4pCentShift_caloJets_noQcut_correctedJetHadron_correctedDihadron_2021-03-04.root", "flowGraphs_PbPb2018_caloJets_dihadronDeltaEta2to3v5_correctedJetHadron_correctedDihadron_2021-08-06.root", "flowGraphs_PbPb2018_caloJets_dihadronDeltaEta2v5to4_correctedJetHadron_correctedDihadron_2021-08-06.root",  "flowGraphs_PbPb2018_caloJets_jetHadronJECminus_correctedJetHadron_correctedDihadron_2021-08-05.root", "flowGraphs_PbPb2018_caloJets_jetHadronJECplus_correctedJetHadron_correctedDihadron_2021-08-05.root",    "flowGraphs_PbPbMC2018_subeNon0_4pCentShift_caloJets_noQcut_correctedJetHadron_correctedDihadron_2021-03-04.root", "flowGraphs_PbPb2018_JECplus_correctedJetHadron_correctedDihadron_2021-07-29.root", "flowGraphs_PbPb2018_JECminus_correctedJetHadron_correctedDihadron_2021-07-29.root", "flowGraphs_PbPb2018MC_caloJets_4pCentShift_subeNon0_gluonJets_correctedJetHadron_correctedDihadron_2021-07-29.root", "flowGraphs_PbPb2018MC_caloJets_4pCentShift_subeNon0_quarkJets_correctedJetHadron_correctedDihadron_2021-07-29.root", "flowGraphs_PbPb2018MC_caloJets_4pCentShift_subeNon0_25pMoreQuark_correctedJetHadron_correctedDihadron_2021-07-26.root", "flowGraphs_PbPb2018MC_caloJets_4pCentShift_subeNon0_quarkJets_correctedJetHadron_correctedDihadron_2021-07-22.root", "flowGraphs_PbPb2018MC_caloJets_4pCentShift_subeNon0_gluonJets_correctedJetHadron_correctedDihadron_2021-07-22.root",  "flowGraphs_PbPbMC2018_1v5pCentShift_subeNon0_caloJets_correctedJetHadron_correctedDihadron_cumulativePtBins_2021-07-12.root", "flowGraphs_PbPbMC2018_1v5pCentShift_onlyDihadron_caloJets_correctedJetHadron_correctedDihadron_cumulativePtBins_2021-07-12.root",   "flowGraphs_PbPbMC2018_5pCentShift_subeNon0_caloJets_onlyDihadron_correctedJetHadron_correctedDihadron_cumulativePtBins_2021-07-02.root"
  };
  
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    graphFile[iFile+1] = TFile::Open(directoryName+comparisonFileName[iFile]);
  }
  
  // Legend text given to each compared file
  TString fileLegend[] = {"Nominal", "No two cuts", "MC, PfCs jets, #eta reflect scale 2", "MC, calo",  "MC control", "MC control2", "MC 16%", "Flow MC", "MC+5% Q < 2.8"};
  
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
  const bool doSummaryCorrection = false;              // Correct the jet vn summary plots based on fits
  const bool manualSummaryCorrection = false;          // Do the correction manially based on tabulated values
  const bool drawSummaryPlot = false;
  
  // Pre-defined points to be drawn to the summary plot
  const bool drawPreviousResults = true;              // Draw ATLAS jet v2 and CMS high pT v2 results to the summary plot
  const bool drawFlowJetResults = false;              // Draw flow subtracted results using default Q-vector configuration
  const bool drawPfCsJetResults = false;              // Draw PfCs jet results using defualt Q-vector configuration
  const bool drawFlowDoubleDijetResults = false;       // Draw flow jet results requiring also calo dijet in the event
  const bool drawCaloDoubleDijetResults = false;      // Draw calo jet results requiring also flow dijet in the event
  
  // Plots to be compared between files
  const bool drawJetHadronVnFileComparison = true;
  const bool drawDihadronVnFileComparison = false;
  const bool drawHadronVnFileComparison = false;
  const bool drawJetVnFileComparison = false;
  const bool drawJetHadronYieldFileComparison = false;
  const bool drawDihadronYieldFileComparison = false;
  const bool drawFileComparison = drawJetHadronVnFileComparison || drawDihadronVnFileComparison || drawHadronVnFileComparison || drawJetVnFileComparison || drawJetHadronYieldFileComparison || drawDihadronYieldFileComparison;
  const bool drawQvectorTrends = false;
  const bool drawHadronVnVersusJetVn = false;
  
  const bool drawSystematicUncertainties = false;     // Include systematic uncertainties in the plots
  
  const bool drawRatios = false;              // Draw ratio plots for file comparison
  const bool ratioToPrevious = false;         // Instead of taking ratio to the first file, take ratio to previous file in the list
  const bool useAlternativeMarkerSet = false; // Alternative marker set optimized for drawing several data and MC collections to the same plot
  
  const bool saveFigures = false;                     // Save the figures in a file
  TString saveComment = "_noTwoCuts";              // String to be added to saved file names
  
  // Saving summary file for final plotter macro
  const bool saveSummaryFile = false;
  const char* outputFileName = "flowGraphs/summaryPlot_akPfCsJet_manualJECinJetHadron_2021-08-30.root";
  
  // Determine from the first comparison file name if we are making Q-cut below or above the threshold
  const char* qVectorType = "below";
  if(comparisonFileName[0].Contains("Above")) qVectorType = "above";
  
  // Determine the used centarality shift from the first comparison file name
  int shiftIndex = comparisonFileName[0].Index("pCentShift");
  int centralityAccuracy = 0;
  double labelShiftNumber = 0;                         // Percentege of centrality shift to be added to labels
  if(shiftIndex > 0){
    if(comparisonFileName[0][shiftIndex-2] == 'v'){
      centralityAccuracy = 1;
      TString numberString = Form("%c.%c", comparisonFileName[0][shiftIndex-3], comparisonFileName[0][shiftIndex-1]);
      labelShiftNumber = numberString.Atof();
    } else {
      labelShiftNumber = comparisonFileName[0][shiftIndex-1]-'0';
    }
  }
  
  // Determine from the first comparison file name whether or not we are using cumulative pT bins
  bool cumulativeBins = false;
  if(comparisonFileName[0].Contains("cumulative")) cumulativeBins = true;
  
  // Determine from the first comparison file name whether we are using only HYDJET or all particles
  const char* subeNon0String = "";
  if(comparisonFileName[0].Contains("subeNon0")) subeNon0String = ", HYDJET";
  
  int firstDrawnAsymmetryBin = nAsymmetryBins;
  int lastDrawnAsymmetryBin = nAsymmetryBins;
  
  int firstDrawnVn = 2;
  int lastDrawnVn = 2;
  
  double maxTrackPt = 4.5;
  int maxPtBin = 4;
  
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
  TGraphErrors *yieldGraphJetHadronVsQvector[maxPtBin][nAsymmetryBins+1][nCentralityBins];
  TGraphErrors *yieldGraphDihadronVsQvector[maxPtBin][nAsymmetryBins+1][nCentralityBins];
 
  TGraphErrors *ratioGraphJetHadron[maxFiles][nAsymmetryBins+1][nCentralityBins];
  TGraphErrors *ratioGraphHadron[maxFiles][nAsymmetryBins+1][nCentralityBins];
  TGraphErrors *ratioGraphJet[maxFiles][nAsymmetryBins+1][nCentralityBins];
  TGraphErrors *ratioGraphJetHadronYield[maxFiles][nAsymmetryBins+1][nCentralityBins];
  TGraphErrors *ratioGraphDihadronYield[maxFiles][nAsymmetryBins+1][nCentralityBins];
  double defaultAxis[4] = {1,2,3,4};
  
  // Graph comparing hadron v2 to jet v2
  TGraphErrors *flowGraphHadronVsJet[maxPtBin][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  
  // Axes for graphs comparing hadron v2 to jet v2
  double xAxisValuesVnComparison[nComparisonFiles];
  double xAxisErrorsVnComparison[nComparisonFiles];
  double yAxisValuesVnComparison[nComparisonFiles];
  double yAxisErrorsVnComparison[nComparisonFiles];
  
  // Graphs as a function of Q-vector cut
  TGraphErrors *flowGraphQvectorJetHadron[nTrackPtBins][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowGraphQvectorJetHadronCorrected[nTrackPtBins][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowGraphQvectorDihadron[nTrackPtBins][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowGraphQvectorHadron[nTrackPtBins][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  TGraphErrors *flowGraphQvectorJet[nTrackPtBins][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  
  // Axes for graphs as a function of Q-vector cut
  double xAxisValuesQvector[nComparisonFiles];
  double xAxisErrorsQvector[nComparisonFiles];
  double yAxisValuesQvector[7][nComparisonFiles];  // Number 7 gives the different graph types
  double yAxisErrorsQvector[7][nComparisonFiles];  // Number 7 gives the different graph types
  
  TLine *qLine[5][nAsymmetryBins+1][nCentralityBins][nFlowComponents][nTrackPtBins];
  TLine *dataLine[maxFiles][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  
  // Determine q-vector zoom range based on cut type
  double qVectorXMin = 1;
  double qVectorXMax = 3.3;
  if(comparisonFileName[0].Contains("Above")){
    qVectorXMin = 1.3;
    qVectorXMax = 3.8;
  }
  double fineTunedQ[nTrackPtBins][nAsymmetryBins+1][nCentralityBins][nFlowComponents];
  
//  // Manual corrections for jet v2 values determined from Q-vector fits
//  double manualCorrection[3][4] = {{0.0891, 0.0850, 0.0735, 0.1287},
//                                   {0.0945, 0.0991, 0.1129, 0.1456},
//                                   {0.0570, 0.0703, 0.0811, 0.1006}};
  
//  // Manual corrections for jet v2 values determined from hadron v2 to jet v2 fits
//  double manualCorrection[3][4] = {{0.0894, 0.0850, 0.0735, 0.1284},
//                                   {0.0947, 0.0993, 0.1131, 0.1454},
//                                   {0.0569, 0.0705, 0.0815, 0.1025}};
  
//  // Manual corrections for jet v2 values determined from hadron v2 to jet v2 fits with only HYDJET particles
//  double manualCorrection[3][4] = {{0.0829, 0.0784, 0.0584, 0.1201},
//                                   {0.0886, 0.0923, 0.0957, 0.1160},
//                                   {0.0546, 0.0556, 0.0475, 0.05}};
  
//  // Manual corrections for jet v2 values determined from hadron v2 and dihadron yield Q-vector tuning with only HYDJET
//  double manualCorrection[3][4] = {{0.0791, 0.0773, 0.0713, 0.1171},
//                                   {0.0886, 0.0923, 0.0957, 0.1160},
//                                   {0.0546, 0.0556, 0.0475, 0.05}};
  
//  // Testing corrections with Q from integrated bins. Done with 4 % centrality shift and calo jets. This is the nominal correction.
//  double manualCorrection[3][4] = {{0.0845, 0.0828, 0.0795, 0.0642},
//                                   {0.0802, 0.0834, 0.0919, 0.0979},
//                                   {0.0512, 0.0492, 0.0490, 0.0479}};
  
//  // Correction with 4 % centrality shift and calo jets. Q-vector cuts are taken as average between dihadron yield and hadron v2.
//  double manualCorrection[3][4] = {{0.0874, 0.0879, 0.0878, 0.0749},
//                                   {0.0791, 0.0823, 0.0904, 0.0973},
//                                   {0.0540, 0.0530, 0.0535, 0.0505}};
  
//  // Variation to average yield and v2. 2:1 ratio towards v2 is used for these results.
//  double manualCorrection[3][4] = {{0.0864, 0.0862, 0.0850, 0.0713},
//                                   {0.0834, 0.0866, 0.0961, 0.0994},
//                                   {0.0580, 0.0585, 0.0599, 0.0542}};
 
//    // Integrated pT bins. Done with 4 % centrality shift and calo jets. Correction from hadron v2, scale from subeNon0 yield
//    double manualCorrection[3][4] = {{0.0845/1.004, 0.0828/1.004, 0.0795/1.004, 0.0642/1.004},
//                                     {0.0920/1.010, 0.0952/1.010, 0.1074/1.010, 0.1037/1.010},
//                                     {0.0661/1.176, 0.0693/1.176, 0.0728/1.176, 0.0616/1.176}};
  
//  // Integrated pT bins. Done with 4 % centrality shift and calo jets. Correction from hadron v2, scale from all yield. REAL NOMINAL CORRECTION!
//  double manualCorrection[3][4] = {{0.0845/1.024, 0.0828/1.024, 0.0795/1.024, 0.0642/1.024},
//                                   {0.0920/1.042, 0.0952/1.042, 0.1074/1.042, 0.1037/1.042},
//                                   {0.0661/1.260, 0.0693/1.260, 0.0728/1.260, 0.0616/1.260}};
  
//    // Integrated pT bins. Done with 4 % centrality shift, calo jets and adjusted quark/gluon jet fraction. Correction from hadron v2, scale from all yield.
//    double manualCorrection[3][4] = {{0.0845/1.024*1.025,  0.0828/1.024*1.015, 0.0795/1.024*1.023,  0.0642/1.024*1.016},
//                                     {0.0920/1.042*0.9974, 0.0952/1.042*1.000, 0.1074/1.042*0.9917, 0.1037/1.042*1.015},
//                                     {0.0661/1.260*1.033,  0.0693/1.260*1.003, 0.0728/1.260*0.9641, 0.0616/1.260*0.8799}};
  
//    // Integrated pT bins. Done with 4 % centrality shift and calo jets. Correction from Q-weighted hadron v2, scale from all yield
//    double manualCorrection[3][4] = {{0.0773/1.02034, 0.0773/1.02034, 0.0773/1.02034, 0.0773/1.02034},
//                                     {0.0990/1.04281, 0.0990/1.04281, 0.0990/1.04281, 0.0990/1.04281},
//                                     {0.0712/1.23788, 0.0712/1.23788, 0.0712/1.23788, 0.0712/1.23788}};
  
//  // Done with 4 % centrality shift and calo jets. Correction for Q-weighted v3, scale from all yield
//  double manualCorrection[3][4] = {{0.030537/1.02034, 0.030537/1.02034, 0.030537/1.02034, 0.030537/1.02034},
//                                   {0.01/1.04281, 0.01/1.04281, 0.01/1.04281, 0.01/1.04281},
//                                   {0.007642/1.23788, 0.007642/1.23788, 0.007642/1.23788, 0.007642/1.23788}};
  
//  // Integrated pT bins. Done with 1.5 % centrality shift and calo jets. Average between hadron v2 and yield roughly matches above 10 %
//  double manualCorrection[3][4] = {{0.0773/0.02034, 0.0773/0.02034, 0.0773/0.02034, 0.0773/0.02034},
//                                   {0.0780, 0.0780, 0.0780, 0.0780},
//                                   {0.0405, 0.0405, 0.0405, 0.0405}};
  
//  // Integrated pT bins. Done with 4 % centrality shift and calo jets. Correction from Q-weighted hadron v2, scale from all yield
//  double manualCorrection[3][4] = {{0.0817/1.02034, 0.0766/1.02034, 0.0764/1.02034, 0.0678/1.02034},
//                                   {0.0947/1.04281, 0.0985/1.04281, 0.1096/1.04281, 0.0997/1.04281},
//                                   {0.0661/0.001, 0.0693/0.001, 0.0728/0.001, 0.0616/0.001}};  // This row not tuned

//  // Testing corrections with Q from integrated bins. Done with 4 % centrality shift and flow jets
//  double manualCorrection[3][4] = {{0.1087, 0.1128, 0.1180, 0.1022},
//                                   {0.0908, 0.0952, 0.0978, 0.0966},
//                                   {0.0340, 0.0351, 0.0358, 0.0286}};
  
//    // Testing corrections with Q from integrated bins. Done with 4 % centrality shift and pfcs jets
//    double manualCorrection[3][4] = {{0.1805, 0.1731, 0.1784, 0.1408},
//                                     {0.1909, 0.1975, 0.2028, 0.2006},
//                                     {0.1178, 0.1215, 0.1232, 0.1238}};
  
//      // Manual corrections directly from 2.5 Q-cuts. Done with 4 % centrality shift and pfcs jets
//      double manualCorrection[3][4] = {{0.1784, 0.1716, 0.1811, 0.1528},
//                                       {0.1883, 0.1946, 0.1982, 0.2010},
//                                       {0.1182, 0.1190, 0.1209, 0.1257}};
  
//    // Testing corrections with Q from integrated bins. Done with 4 % centrality shift and flow jets requiring also calo dijet
//    double manualCorrection[3][4] = {{0.0722, 0.0743, 0.0855, 0.1055},
//                                     {0.0691, 0.0715, 0.0784, 0.1074},
//                                     {0.0423, 0.0451, 0.0532, 0.0911}};
  
//  // Testing corrections with Q from integrated bins. Done with 4 % centrality shift and calo jets requiring also flow dijet
//  double manualCorrection[3][4] = {{0.0709, 0.0683, 0.0751, 0.0654},
//                                   {0.0639, 0.0656, 0.0657, 0.0652},
//                                   {0.0378, 0.0373, 0.0337, 0.0181}};
  
//    // Integrated pT bins. Done with 4.5 % centrality shift and calo jets. Correction directly from hadron v2 value is used.
//    double manualCorrection[3][4] = {{0.0833/0.968, 0.0761/0.968, 0.0732/0.968, 0.0539/0.968},
//                                     {1, 1, 1, 1},
//                                     {1, 1, 1, 1}};
  
//      // Integrated pT bins. Done with 4.5 % centrality shift and calo jets. Correction from average dihadron yield and hadron v2.
//      double manualCorrection[3][4] = {{0.0826, 0.0725, 0.0677, 0.0445},
//                                       {1, 1, 1, 1},
//                                       {1, 1, 1, 1}};
  
//  // Testing corrections with Q from integrated bins, only 0-10 actually done. Done with 4.5 % centrality shift and flow jets
//  double manualCorrection[3][4] = {{0.1101, 0.1117, 0.1116, 0.0965},
//                                   {0.0872, 0.0923, 0.0948, 0.0941},
//                                   {0.0332, 0.0289, 0.0306, 0.0266}};
  
//  // Integrated pT bins. Done with 5 % centrality shift and calo jets. Correction from hadron v2, scale from subeNon0 yield
//  double manualCorrection[3][4] = {{0.0755/0.930, 0.0685/0.930, 0.0637/0.930, 0.0641/0.930},
//                                   {0.0863/0.933, 0.0939/0.933, 0.1004/0.933, 0.1010/0.933},
//                                   {0.0536/1.057, 0.0576/1.057, 0.0684/1.057, 0.05/1.057}};
  
//    // Integrated pT bins. Done with 5 % centrality shift and calo jets. Correction from hadron v2, scale from all yield
//    double manualCorrection[3][4] = {{0.0755/0.949, 0.0685/0.949, 0.0637/0.949, 0.0641/0.949},
//                                     {0.0863/0.964, 0.0939/0.964, 0.1004/0.964, 0.1010/0.964},
//                                     {0.0536/1.134, 0.0576/1.134, 0.0684/1.134, 0.05/1.134}};
  
//  // Integrated pT bins. Done with 5 % centrality shift and calo jets. Correction from average dihadron yield and hadron v2.
//  double manualCorrection[3][4] = {{0.0657, 0.0483, 0.0385, 0.0274},
//                                   {0.0946, 0.1040, 0.1117, 0.1082},
//                                   {0.0497, 0.0505, 0.0576, 0.05}};
  
//    // Integrated pT bins. Done with 3 % centrality shift and calo jets. Correction from hadron v2, scale from subeNon0 yield
//    double manualCorrection[3][4] = {{0.0927/1.071, 0.0879/1.071, 0.0913/1.071, 0.0707/1.071},
//                                     {0.0992/1.085, 0.1013/1.085, 0.1162/1.085, 0.1107/1.085},
//                                     {0.0682/1.300, 0.0718/1.300, 0.0742/1.300, 0.0611/1.300}};
  
//      // Integrated pT bins. Done with 3 % centrality shift and calo jets. Correction from hadron v2, scale from all yield
//      double manualCorrection[3][4] = {{0.0927/1.087, 0.0879/1.087, 0.0913/1.087, 0.0707/1.087},
//                                       {0.0992/1.119, 0.1013/1.119, 0.1162/1.119, 0.1107/1.119},
//                                       {0.0682/1.387, 0.0718/1.387, 0.0742/1.387, 0.0611/1.387}};

//    // Integrated pT bins. Done with 3 % centrality shift and calo jets. Correction from average dihadron yield and hadron v2.
//    double manualCorrection[3][4] = {{0.1013, 0.0994, 0.1143, 0.0970},
//                                     {0.0543, 0.0603, 0.0558, 0.0948},
//                                     {0.0501, 0.0501, 0.0509, 0.0510}};
  
//        // Integrated pT bins. Done with no centrality shift and calo jets. Correction from hadron v2, scale from all yield
//        double manualCorrection[3][4] = {{0.1544/1.285, 0.1529/1.285, 0.1609/1.285, 0.1493/1.285},
//                                         {0.1471/1.378, 0.1549/1.378, 0.1818/1.378, 0.1843/1.378},
//                                         {0.0891/1.844, 0.1164/1.844, 0.1061/1.844, 0.0909/1.844}};
  
        // Disable MC correction.
        double manualCorrection[3][4] = {{0, 0, 0, 0},
                                         {0, 0, 0, 0},
                                         {0, 0, 0, 0}};
  
  // For x-axis, set the Q-vector value from the file. This can be decoded from the filename
  TObjArray *nameContents;
  TString thisToken;
  double qValue, xValueForQ, yValueForQ, xErrorForQ, yErrorForQ;
  
  // Loop over all file names
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){

    // No error in Q-vector cut
    xAxisErrorsQvector[iFile] = 0;
    
    // If qVector is not mentioned in the file name, set the cut to 0
    if(!comparisonFileName[iFile].Contains("qVector")){
      xAxisValuesQvector[iFile] = 0;
      continue;
    }

    // If qVector is mentioned in the file name, read the cut value and add it to the array
    nameContents = comparisonFileName[iFile].Tokenize("_");                          // Tokenize the string from underscores
    for(int iToken = 0; iToken < nameContents->GetEntriesFast(); iToken++){    // Loop over all tokens
      thisToken = ((TObjString*)(nameContents->At(iToken)))->String();         // Transform tokens back to string
      if(!thisToken.Contains("qVector")) continue;                             // Check if token contains qVector
      thisToken.Replace(0,12,"");                                              // Remove qVectorAbove or qVectorBelow from the token
      thisToken.ReplaceAll("p",".");                                           // Replace p with .
      xAxisValuesQvector[iFile] = thisToken.Atof();                            // Convert remaining string to double and add to array
      break;                                                                   // Once qVector found, no need to check remaining tokens
    }
    
  } // File loop for setting q-vector values
  
  // Jet vn summary graph in all centrality bins
  TGraphErrors *flowSummaryJet[maxFiles][nAsymmetryBins+1][nFlowComponents];
  TGraphErrors *atlasJetV2graph;
  TGraphErrors *cmsHighPtV2;
  TGraphErrors *flowJetReferenceGraph;
  TGraphErrors *pfCsJetReferenceGraph;
  TGraphErrors *doubleDijetFlowJetReferenceGraph;
  TGraphErrors *doubleDijetCaloJetReferenceGraph;
  double summaryXaxis[nCentralityBins];
  double summaryXaxisError[nCentralityBins];
  double summaryYaxis[maxFiles][nAsymmetryBins+1][nFlowComponents][nCentralityBins];
  double summaryYaxisError[maxFiles][nAsymmetryBins+1][nFlowComponents][nCentralityBins];
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    summaryXaxis[iCentrality] = iCentrality+1;
    summaryXaxisError[iCentrality] = 0;
    for(int iFile = 0; iFile < maxFiles; iFile++){
      for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins+1; iAsymmetry++){
        ratioGraphJetHadron[iFile][iAsymmetry][iCentrality] = new TGraphErrors(maxPtBin, defaultAxis, defaultAxis);
        ratioGraphHadron[iFile][iAsymmetry][iCentrality] = new TGraphErrors(maxPtBin, defaultAxis, defaultAxis);
        ratioGraphJet[iFile][iAsymmetry][iCentrality] = new TGraphErrors(maxPtBin, defaultAxis, defaultAxis);
        ratioGraphJetHadronYield[iFile][iAsymmetry][iCentrality] = new TGraphErrors(maxPtBin, defaultAxis, defaultAxis);
        ratioGraphDihadronYield[iFile][iAsymmetry][iCentrality] = new TGraphErrors(maxPtBin, defaultAxis, defaultAxis);
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
  
  // Helper variables for drawing figures
  TLegend *legend;
  TLegend *vLegend;
  int fullMarkers[] = {kFullSquare, kFullCircle, kFullDiamond, kFullCross, kFullStar, kFullFourTrianglesPlus, kFullDoubleDiamond};
  int secondMarkers[] = {kFullCircle, kFullCross, kFullSquare, kFullCircle, kFullFourTrianglesPlus};
  int colors[] = {kBlue,kRed,kGreen+2,kBlack, kMagenta};
  int flowColors[] = {kBlue, kBlack, kRed, kGreen+3, kMagenta};
  int fileColors[] = {kBlack, kBlue, kRed, kGreen+3, kMagenta, kCyan, kOrange+2};
  TString asymmetryString[] = {" 0.0 < x_{j} < 0.6", " 0.6 < x_{j} < 0.8", " 0.8 < x_{j} < 1.0", ""};
  TString asymmetryLegend[] = {"0.0 < x_{j} < 0.6", "0.6 < x_{j} < 0.8", "0.8 < x_{j} < 1.0", "x_{j} integrated"};
  TString compactAsymmetryString[] = {"_A=0v0-0v6", "_A=0v6-0v8", "_A=0v8-1v0", ""};
  double normalizationValue[2]; // Normalizer for yields
  
  if(useAlternativeMarkerSet){
    int alternativeMarkers[] = {kFullSquare, kOpenSquare, kFullCircle, kOpenCircle, kFullDiamond, kOpenDiamond, kFullDoubleDiamond};
    int alternativeColors[] = {kBlack, kBlack, kBlue, kBlue, kRed, kRed, kViolet};
    for(int i = 0; i < maxFiles; i++){
      fullMarkers[i] = alternativeMarkers[i];
      fileColors[i] = alternativeColors[i];
    }
  }

  
  // After all the graphs have been loaded from the input files, combine information from these graphs to a new one
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
      for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
        for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
          for(int iFile = 0; iFile < nComparisonFiles+1; iFile++){
            
            // Jet-hadron correlations
            flowGraphJetHadron[iFile][iAsymmetry][iCentrality][iFlow]->GetPoint(iPoint, xValueForQ, yValueForQ);
            yErrorForQ = flowGraphJetHadron[iFile][iAsymmetry][iCentrality][iFlow]->GetErrorY(iPoint);
            
            if(iFile == 0){
              qLine[0][iAsymmetry][iCentrality][iFlow][iPoint] = new TLine(qVectorXMin, yValueForQ, qVectorXMax, yValueForQ);
              qLine[0][iAsymmetry][iCentrality][iFlow][iPoint]->SetLineColor(fileColors[iPoint]);
              qLine[0][iAsymmetry][iCentrality][iFlow][iPoint]->SetLineStyle(2);
            } else {
              yAxisValuesQvector[0][iFile-1] = yValueForQ;
              yAxisErrorsQvector[0][iFile-1] = yErrorForQ;
            }
            
            // Jet-hadron correlation, corrected (obsolete)
            flowGraphJetHadronCorrected[iFile][iAsymmetry][iCentrality][iFlow]->GetPoint(iPoint, xValueForQ, yValueForQ);
            yErrorForQ = flowGraphJetHadronCorrected[iFile][iAsymmetry][iCentrality][iFlow]->GetErrorY(iPoint);
            
            if(iFile == 0){
              qLine[1][iAsymmetry][iCentrality][iFlow][iPoint] = new TLine(qVectorXMin, yValueForQ, qVectorXMax, yValueForQ);
              qLine[1][iAsymmetry][iCentrality][iFlow][iPoint]->SetLineColor(fileColors[iPoint]);
              qLine[1][iAsymmetry][iCentrality][iFlow][iPoint]->SetLineStyle(2);
            } else {
              yAxisValuesQvector[1][iFile-1] = yValueForQ;
              yAxisErrorsQvector[1][iFile-1] = yErrorForQ;
            }
            
            // Dihadron correlations
            flowGraphDihadron[iFile][iAsymmetry][iCentrality][iFlow]->GetPoint(iPoint, xValueForQ, yValueForQ);
            yErrorForQ = flowGraphDihadron[iFile][iAsymmetry][iCentrality][iFlow]->GetErrorY(iPoint);
            
            if(iFile == 0){
              qLine[2][iAsymmetry][iCentrality][iFlow][iPoint] = new TLine(qVectorXMin, yValueForQ, qVectorXMax, yValueForQ);
              qLine[2][iAsymmetry][iCentrality][iFlow][iPoint]->SetLineColor(fileColors[iPoint]);
              qLine[2][iAsymmetry][iCentrality][iFlow][iPoint]->SetLineStyle(2);
            } else {
              yAxisValuesQvector[2][iFile-1] = yValueForQ;
              yAxisErrorsQvector[2][iFile-1] = yErrorForQ;
            }
            
            // Single hadron vn
            flowGraphHadron[iFile][iAsymmetry][iCentrality][iFlow]->GetPoint(iPoint, xValueForQ, yValueForQ);
            yErrorForQ = flowGraphHadron[iFile][iAsymmetry][iCentrality][iFlow]->GetErrorY(iPoint);
            
            if(iFile == 0){
              qLine[3][iAsymmetry][iCentrality][iFlow][iPoint] = new TLine(qVectorXMin, yValueForQ, qVectorXMax, yValueForQ);
              qLine[3][iAsymmetry][iCentrality][iFlow][iPoint]->SetLineColor(fileColors[iPoint]);
              qLine[3][iAsymmetry][iCentrality][iFlow][iPoint]->SetLineStyle(2);
            } else {
              yAxisValuesQvector[3][iFile-1] = yValueForQ;
              yAxisErrorsQvector[3][iFile-1] = yErrorForQ;
              xAxisValuesVnComparison[iFile-1] = yValueForQ;
              xAxisErrorsVnComparison[iFile-1] = yErrorForQ;
            }
            
            // Jet vn
            flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->GetPoint(iPoint, xValueForQ, yValueForQ);
            yErrorForQ = flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->GetErrorY(iPoint);
            
            if(iFile == 0){
              qLine[4][iAsymmetry][iCentrality][iFlow][iPoint] = new TLine(qVectorXMin, yValueForQ, qVectorXMax, yValueForQ);
              qLine[4][iAsymmetry][iCentrality][iFlow][iPoint]->SetLineColor(fileColors[iPoint]);
              qLine[4][iAsymmetry][iCentrality][iFlow][iPoint]->SetLineStyle(2);
            } else {
              yAxisValuesQvector[4][iFile-1] = yValueForQ;
              yAxisErrorsQvector[4][iFile-1] = yErrorForQ;
              yAxisValuesVnComparison[iFile-1] = yValueForQ;
              yAxisErrorsVnComparison[iFile-1] = yErrorForQ;
            }
            
          } // Loop over files
          
          // After the points from all the files have been collected, these can be put into new formatted graphs
          flowGraphQvectorJetHadron[iPoint][iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nComparisonFiles, xAxisValuesQvector, yAxisValuesQvector[0], xAxisErrorsQvector, yAxisErrorsQvector[0]);
          flowGraphQvectorJetHadronCorrected[iPoint][iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nComparisonFiles, xAxisValuesQvector, yAxisValuesQvector[1], xAxisErrorsQvector, yAxisErrorsQvector[1]);
          flowGraphQvectorDihadron[iPoint][iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nComparisonFiles, xAxisValuesQvector, yAxisValuesQvector[2], xAxisErrorsQvector, yAxisErrorsQvector[2]);
          flowGraphQvectorHadron[iPoint][iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nComparisonFiles, xAxisValuesQvector, yAxisValuesQvector[3], xAxisErrorsQvector, yAxisErrorsQvector[3]);
          flowGraphQvectorJet[iPoint][iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nComparisonFiles, xAxisValuesQvector, yAxisValuesQvector[4], xAxisErrorsQvector, yAxisErrorsQvector[4]);
          
          flowGraphHadronVsJet[iPoint][iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nComparisonFiles, xAxisValuesVnComparison, yAxisValuesVnComparison, xAxisErrorsVnComparison, yAxisErrorsVnComparison);
                    
        } // Point loop in the graph
      } // Flow component loop
      
      // Combination of information for yields
      if(drawJetHadronYieldFileComparison || drawDihadronYieldFileComparison){
        for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
          for(int iFile = 0; iFile < nComparisonFiles+1; iFile++){
            
            // Dihadron yield
            yieldGraphDihadron[iFile][iAsymmetry][iCentrality]->GetPoint(iPoint, xValueForQ, yValueForQ);
            yErrorForQ = yieldGraphDihadron[iFile][iAsymmetry][iCentrality]->GetErrorY(iPoint);
            
            if(iFile == 0){
              normalizationValue[0] = yValueForQ;
            } else {
              yAxisValuesQvector[5][iFile-1] = yValueForQ / normalizationValue[0];
              yAxisErrorsQvector[5][iFile-1] = yErrorForQ / normalizationValue[0];
            }
            
            // Jet-hadron yield
            yieldGraphJetHadron[iFile][iAsymmetry][iCentrality]->GetPoint(iPoint, xValueForQ, yValueForQ);
            yErrorForQ = yieldGraphJetHadron[iFile][iAsymmetry][iCentrality]->GetErrorY(iPoint);
            
            if(iFile == 0){
              normalizationValue[1] = yValueForQ;
            } else {
              yAxisValuesQvector[6][iFile-1] = yValueForQ / normalizationValue[1];
              yAxisErrorsQvector[6][iFile-1] = yErrorForQ / normalizationValue[1];
            }
            
          } // Loop over files
          
          // After the points from all the files have been collected, these can be put into new formatted graphs
          yieldGraphDihadronVsQvector[iPoint][iAsymmetry][iCentrality] = new TGraphErrors(nComparisonFiles, xAxisValuesQvector, yAxisValuesQvector[5], xAxisErrorsQvector, yAxisErrorsQvector[5]);
          yieldGraphJetHadronVsQvector[iPoint][iAsymmetry][iCentrality] = new TGraphErrors(nComparisonFiles, xAxisValuesQvector, yAxisValuesQvector[6], xAxisErrorsQvector, yAxisErrorsQvector[6]);
          
        } // Point loop in the graph
        
      } // Doing yields
      
    } // Asymmetry loop
  } // Centrality loop
  
  // Setup the drawer for graphs
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceGraph();
  
  
  // How to zoom vn plots    //   0-10 10-30 30-50 50-100  pp
  double vZoomTable[4][5] = {{    0.2,  0.25, 0.35,  1,   1.4},  // v1
                             {    0.12, 0.2,  0.3,  0.5,  0.6}, // v2
                             {    0.12, 0.12, 0.12, 0.2,  0.25},  // v3
                             {    0.1,  0.1,  0.1,  0.1,  0.1}}; // v4
  
  // How to zoom vn plots    //    v1,   v2,   v3,   v4
  double jetHadronZoomTable[4] = {0.1, 0.03, 0.005, 0.01};
  double dihadronZoomTable[4]  = {0.08, 0.08, 0.03, 0.03};

  // Numbers from HP2020 conference presentation
  TLine *atlasV2;
  double atlasV2Number[] = {0.018, 0.03, 0.035, 0.03};
  
  // Numbers from fitting the data from paper arXiv:1702.00630 at large pT
  double cmsHighPtV2Number[] = {0.0220, 0.0376, 0.0431, 0.04};
  double cmsHighPtV2Error[] = {0.0019, 0.0016, 0.0027, 0.04};
  
  // Numbers from flow jets
  //double flowJetV2Number[] = {0.0277755, 0.0614633, 0.0550086, 0.02}; // 4 % centrality shift
  //double flowJetV2Error[] = {0.00177884, 0.00111697, 0.0014466, 0.02}; // 4 % centrality shift
  //double flowJetV2Number[] = {0.029474, 0.0645204, 0.0598835, 0.02}; // 5 % centrality shift (4.5 % for 0-10 bin)
  //double flowJetV2Error[] = {0.00177884, 0.00111697, 0.0014466, 0.02}; // 5 % centrality shift (4.5 % for 0-10 bin)
  double flowJetV2Number[] = {0.025417000, 0.062135500, 0.055488200, 0.02}; // 4 % centrality shift  (no fit MC, just 2.5 Q cut)
  double flowJetV2Error[] = {0.0021964265, 0.0011193279, 0.0016713376, 0.02}; // 4 % centrality shift  (no fit MC, just 2.5 Q cut)
  //double pfCsJetV2Number[] = {0.0278213, 0.0490602, 0.0456264, 0.02}; // 4 % centrality shift  (calo dihadron)
  //double pfCsJetV2Error[] = {0.00141872, 0.00082865, 0.00101112, 0.02}; // 4 % centrality shift  (calo dihadron)
  double pfCsJetV2Number[] = {0.022202000, 0.044828000, 0.043780000, 0.02}; // 4 % centrality shift  (no fit MC, just 2.5 Q cut)
  double pfCsJetV2Error[] = {0.0022051142, 0.0011051714, 0.0016032137, 0.02}; // 4 % centrality shift  (no fit MC, just 2.5 Q cut)
  
  
  
  //double doubleDijetCaloJetNumber[] = {0.0149655, 0.038518, 0.0416287, 0.02}; // 4 % centrality shift
  //double doubleDijetCaloJetError[] = {0.00156909, 0.000762696, 0.000767619, 0.02}; // 4 % centrality shift
  double doubleDijetCaloJetNumber[] = {0.035750000, 0.0086330000, -0.0050535000, 0.02}; // 4 % centrality shift, high subleading jet
  double doubleDijetCaloJetError[] = {0.00156909, 0.000762696, 0.000767619, 0.02}; // 4 % centrality shift, high subleading jet
  
  //double doubleDijetFlowJetNumber[] = {0.0106847, 0.034092, 0.0331465, 0.02}; // 4 % centrality shift
  //double doubleDijetFlowJetError[] = {0.00156909, 0.000762696, 0.000767619, 0.02}; // 4 % centrality shift
  double doubleDijetFlowJetNumber[] = {0.059872800, 0.029541000, 0.0026188000, 0.02}; // 4 % centrality shift, high subleading jet
  double doubleDijetFlowJetError[] = {0.00156909, 0.000762696, 0.000767619, 0.02}; // 4 % centrality shift, high subleading jet

  
  TString legendString;
  char namerY[100];
  
  TLine *zeroLine = new TLine(0,0,maxTrackPt,0);
  zeroLine->SetLineStyle(2);
  
  TLine *shortZeroLine = new TLine(0.8,0,3.2,0);
  shortZeroLine->SetLineStyle(2);
  
  TLine *vnLine = new TLine(0,0,maxTrackPt,0);
  vnLine->SetLineStyle(2);
  
  TLine *oneLine = new TLine(0,1,maxTrackPt,1);
  oneLine->SetLineStyle(2);
  
  TLine *oneQLine = new TLine(qVectorXMin,1,qVectorXMax,1);
  oneQLine->SetLineStyle(2);
  
  TH1D *vnLineError = new TH1D("vnLineError","vnLineError",1,0,maxTrackPt);
  
  double vnValue, vnError;
  
  //double minZoomJetVn[] = {0.085, 0.11, 0.075, 0.05};  // Nominal
  //double maxZoomJetVn[] = {0.135, 0.16, 0.125, 0.15};  // Nominal
  
  double minZoomJetVn[] = {-0.05, -0.05, -0.05, 0};
  double maxZoomJetVn[] = {0.28, 0.28, 0.28, 0.15};
  
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
  
  double xPoint1, xPoint2, yPoint1, yPoint2, yError1, yError2, combinedError, ratioValue, crossingPoint;
  int lowPoint;
  int maxFileIndex, baseIndex, comparisonIndex;
  
  // Compare graphs from different files
  if(drawFileComparison){
    
    // TODO: Debug for nice drawing of super zoom
    //drawer->SetTitleOffsetY(1.7);
    //drawer->SetTitleOffsetX(1.1);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
          
          if(drawJetHadronVnFileComparison){
          
            sprintf(namerY,"Jet-hadron V_{%d}",iFlow+1);
            legend = new TLegend(0.2,0.6,0.5,0.9);
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
            
            // Point by point ratio to the first graph or to the previous graph
            if(drawRatios){
              legend = new TLegend(0.22,0.76,0.52,0.91); // 0.23,0.19,0.53,0.39
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->SetHeader(Form("Cent: %.0f-%.0f%%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
              
              maxFileIndex = ratioToPrevious ? (nComparisonFiles+1)/2 : nComparisonFiles;
              
              for(int iFile = 0; iFile < maxFileIndex; iFile++){
                for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
                  
                  baseIndex = ratioToPrevious ? iFile*2 : 0;
                  comparisonIndex = ratioToPrevious ? iFile*2+1 : iFile+1;
                  
                  flowGraphJetHadron[baseIndex][iAsymmetry][iCentrality][iFlow]->GetPoint(iPoint, xPoint1, yPoint1);
                  yError1 = flowGraphJetHadron[baseIndex][iAsymmetry][iCentrality][iFlow]->GetErrorY(iPoint);
                  flowGraphJetHadron[comparisonIndex][iAsymmetry][iCentrality][iFlow]->GetPoint(iPoint, xPoint2, yPoint2);
                  yError2 = flowGraphJetHadron[comparisonIndex][iAsymmetry][iCentrality][iFlow]->GetErrorY(iPoint);
                  
                  ratioValue = yPoint2 / yPoint1;
                  combinedError = TMath::Sqrt(TMath::Power(yError2/yPoint1,2) + TMath::Power((yPoint2*yError1)/(yPoint1*yPoint1),2));
                  //combinedError = TMath::Abs( ((1.0 - 2.0 * yPoint2 / yPoint1) * yError2*yError2 + yPoint2*yPoint2 * yError1*yError1 / (yPoint1*yPoint1) ) / (yPoint1*yPoint1)); // Binomial error calculation
                  
                  ratioGraphJetHadron[iFile][iAsymmetry][iCentrality]->SetPoint(iPoint,xPoint1,ratioValue);
                  ratioGraphJetHadron[iFile][iAsymmetry][iCentrality]->SetPointError(iPoint, 0, combinedError);
                  
                } // Point loop
                
                ratioGraphJetHadron[iFile][iAsymmetry][iCentrality]->SetMarkerStyle(fullMarkers[ratioToPrevious ? baseIndex : comparisonIndex]);
                ratioGraphJetHadron[iFile][iAsymmetry][iCentrality]->SetMarkerSize(1.3);
                ratioGraphJetHadron[iFile][iAsymmetry][iCentrality]->SetMarkerColor(fileColors[ratioToPrevious ? baseIndex : comparisonIndex]);
                
                legend->AddEntry(ratioGraphJetHadron[iFile][iAsymmetry][iCentrality], Form("%s / %s", fileLegend[comparisonIndex].Data(), fileLegend[baseIndex].Data()), "p");
                
                if(iFile == 0){
                  drawer->DrawGraph(ratioGraphJetHadron[iFile][iAsymmetry][iCentrality], 0, maxTrackPt, 0.995, 1.005, "Track p_{T} (GeV)", Form("Jet-hadron V_{%d} ratio",iFlow+1), " ", "p");
                } else {
                  ratioGraphJetHadron[iFile][iAsymmetry][iCentrality]->Draw("p,same");
                }
                
              } // File loop
              
              oneLine->Draw();
              legend->Draw();
              
              // Save the figures to file
              if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/jetHadronV%dRatio%s%s_C=%.0f-%.0f.pdf", iFlow+1, saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
              }
              
            }
            
            // Plots as a function of Q-vector
            if(drawQvectorTrends){
              
              legend = new TLegend(0.2,0.65,0.5,0.9);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->SetHeader(Form("Cent: %.0f-%.0f%%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
              
              for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
                flowGraphQvectorJetHadron[iPoint][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iPoint]);
                flowGraphQvectorJetHadron[iPoint][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(fileColors[iPoint]);
                flowGraphQvectorJetHadron[iPoint][iAsymmetry][iCentrality][iFlow]->SetMarkerSize(1.3);
                if(iPoint == 0){
                  drawer->DrawGraph(flowGraphQvectorJetHadron[iPoint][iAsymmetry][iCentrality][iFlow],-0.5, 4, 0, jetHadronZoomTable[iFlow], "Q-vector", namerY, " ", "p");
                } else {
                  flowGraphQvectorJetHadron[iPoint][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
                }
                
                legend->AddEntry(flowGraphQvectorJetHadron[iPoint][iAsymmetry][iCentrality][iFlow], Form("%.1f < p_{T} < %.1f GeV", trackPtBinBorders[iPoint], trackPtBinBorders[iPoint+1]), "p");
                
              }
              
              legend->Draw();
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
            
            // Plots as a function of Q-vector
            if(drawQvectorTrends){
              legend = new TLegend(0.2,0.65,0.5,0.9);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->SetHeader(Form("Cent: %.0f-%.0f%%%s%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data(), subeNon0String));
              
              for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
                flowGraphQvectorDihadron[iPoint][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iPoint]);
                flowGraphQvectorDihadron[iPoint][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(fileColors[iPoint]);
                flowGraphQvectorDihadron[iPoint][iAsymmetry][iCentrality][iFlow]->SetMarkerSize(1.3);
                if(iPoint == 0){
                  drawer->DrawGraph(flowGraphQvectorDihadron[iPoint][iAsymmetry][iCentrality][iFlow],-0.5, 4, 0, dihadronZoomTable[iFlow], "Q-vector", namerY, " ", "p");
                } else {
                  flowGraphQvectorDihadron[iPoint][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
                }
                
                lowPoint = iPoint;
                if(cumulativeBins) lowPoint = 0;
                legend->AddEntry(flowGraphQvectorDihadron[iPoint][iAsymmetry][iCentrality][iFlow], Form("%.1f < p_{T} < %.1f GeV", trackPtBinBorders[lowPoint], trackPtBinBorders[iPoint+1]), "p");
                
              }
              
              legend->Draw();
            }
            
          } // File comparison for jet-hadron Vn
          
          if(drawHadronVnFileComparison){
          
            sprintf(namerY,"Hadron v_{%d}",iFlow+1);
            legend = new TLegend(0.2,0.6,0.5,0.9);
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
              
            } // File loop for hadron vn
            
            legend->Draw();
            
            // Save the figures to file
            if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/hadronV%dComparison%s%s_C=%.0f-%.0f.pdf", iFlow+1, saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
            }
            
            // Point by point ratio
            if(drawRatios){
              
              legend = new TLegend(0.23,0.7,0.53,0.9); // 0.23,0.19,0.53,0.39  //  0.23,0.7,0.53,0.9
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->SetHeader(Form("Cent: %.0f-%.0f%%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
              
              // Select ratio either to the first graph or to the previous on the list
              maxFileIndex = ratioToPrevious ? (nComparisonFiles+1)/2 : nComparisonFiles;
              
              for(int iFile = 0; iFile < maxFileIndex; iFile++){
                for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
                  
                  baseIndex = ratioToPrevious ? iFile*2 : 0;
                  comparisonIndex = ratioToPrevious ? iFile*2+1 : iFile+1;
                  
                  flowGraphHadron[baseIndex][iAsymmetry][iCentrality][iFlow]->GetPoint(iPoint, xPoint1, yPoint1);
                  yError1 = flowGraphHadron[baseIndex][iAsymmetry][iCentrality][iFlow]->GetErrorY(iPoint);
                  flowGraphHadron[comparisonIndex][iAsymmetry][iCentrality][iFlow]->GetPoint(iPoint, xPoint2, yPoint2);
                  yError2 = flowGraphHadron[comparisonIndex][iAsymmetry][iCentrality][iFlow]->GetErrorY(iPoint);
                  
                  ratioValue = yPoint2 / yPoint1;
                  combinedError = TMath::Sqrt(TMath::Power(yError2/yPoint1,2) + TMath::Power((yPoint2*yError1)/(yPoint1*yPoint1),2));
                  
                  ratioGraphHadron[iFile][iAsymmetry][iCentrality]->SetPoint(iPoint,xPoint1,ratioValue);
                  ratioGraphHadron[iFile][iAsymmetry][iCentrality]->SetPointError(iPoint, 0, combinedError);
                  
                } // Point loop
                
                ratioGraphHadron[iFile][iAsymmetry][iCentrality]->SetMarkerStyle(fullMarkers[ratioToPrevious ? baseIndex : comparisonIndex]);
                ratioGraphHadron[iFile][iAsymmetry][iCentrality]->SetMarkerSize(1.3);
                ratioGraphHadron[iFile][iAsymmetry][iCentrality]->SetMarkerColor(fileColors[ratioToPrevious ? baseIndex : comparisonIndex]);
                
                legend->AddEntry(ratioGraphHadron[iFile][iAsymmetry][iCentrality], Form("%s / %s", fileLegend[comparisonIndex].Data(), fileLegend[baseIndex].Data()), "p");
                
                if(iFile == 0){
                  drawer->DrawGraph(ratioGraphHadron[iFile][iAsymmetry][iCentrality], 0, maxTrackPt, 0.8, 1.2, "Track p_{T} (GeV)", "Hadron v_{2} ratio", " ", "p");
                } else {
                  ratioGraphHadron[iFile][iAsymmetry][iCentrality]->Draw("p,same");
                }
              } // File loop
              
              oneLine->Draw();
              legend->Draw();
              
              // Save the figures to file
              if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/hadronV%dRatio%s%s_C=%.0f-%.0f.pdf", iFlow+1, saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
              }
              
            }
            
            // Plots as a function of Q-vector
            if(drawQvectorTrends){
              legend = new TLegend(0.2,0.6,0.5,0.9);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->SetHeader(Form(Form("Cent: %%.0f-%%.0f%%%%%%s, %%.%df%%%% shifted%%s", centralityAccuracy), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data(), labelShiftNumber, subeNon0String));
              
              for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
                flowGraphQvectorHadron[iPoint][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iPoint]);
                flowGraphQvectorHadron[iPoint][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(fileColors[iPoint]);
                flowGraphQvectorHadron[iPoint][iAsymmetry][iCentrality][iFlow]->SetMarkerSize(1.3);
                flowGraphQvectorHadron[iPoint][iAsymmetry][iCentrality][iFlow]->Fit("pol1");
                flowGraphQvectorHadron[iPoint][iAsymmetry][iCentrality][iFlow]->GetFunction("pol1")->SetLineColor(fileColors[iPoint]);
                if(iPoint == 0){
                  drawer->DrawGraph(flowGraphQvectorHadron[iPoint][iAsymmetry][iCentrality][iFlow], qVectorXMin, qVectorXMax, 0, 0.35, Form("Q-vector %s", qVectorType), namerY, " ", "p");
                } else {
                  flowGraphQvectorHadron[iPoint][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
                }
                
                qLine[3][iAsymmetry][iCentrality][iFlow][iPoint]->Draw();
                
                // Solve where the function crosses the line
                crossingPoint = findCrossingPoint(flowGraphQvectorHadron[iPoint][iAsymmetry][iCentrality][iFlow]->GetFunction("pol1"), qLine[3][iAsymmetry][iCentrality][iFlow][iPoint]);
                fineTunedQ[iPoint][iAsymmetry][iCentrality][iFlow] = crossingPoint;
                
                lowPoint = iPoint;
                if(cumulativeBins) lowPoint = 0;
                legend->AddEntry(flowGraphQvectorHadron[iPoint][iAsymmetry][iCentrality][iFlow], Form("%.1f < p_{T} < %.1f GeV, Q = %.3f", trackPtBinBorders[lowPoint], trackPtBinBorders[iPoint+1], crossingPoint), "p");
                
              }
              
              legend->Draw();
              
              vLegend = new TLegend(0.19,0.18,0.89,0.28);
              vLegend->SetFillStyle(0);vLegend->SetBorderSize(0);vLegend->SetTextSize(0.04);vLegend->SetTextFont(62);
              vLegend->SetHeader("Solid line = fit to MC, Dashed lines = data");
              vLegend->Draw();
            }
            
          } // File comparison for hadron vn
          
          if(drawJetVnFileComparison){
          
            sprintf(namerY,"Jet v_{%d}",iFlow+1);
            legend = new TLegend(0.22,0.65,0.52,0.9); // 0.2,0.7,0.5,0.9  //  0.22,0.65,0.52,0.98
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->SetHeader(Form("Cent: %.0f-%.0f%%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
            
//            // TODO: Remove these
//            cout << "Points below Q" << endl;
//            for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
//              flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->GetPoint(iPoint, xPoint1, yPoint1);
//              cout << yPoint1 << endl;
//
//            }
//            cout << "Points above Q" << endl;
//            for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
//              flowGraphJet[1][iAsymmetry][iCentrality][iFlow]->GetPoint(iPoint, xPoint1, yPoint1);
//              cout << yPoint1 << endl;
//
//            }
            
            // Do manual correction for jet v2
            if(manualSummaryCorrection){
              for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
                flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->GetPoint(iPoint, xPoint1, yPoint1);
                flowGraphJet[1][iAsymmetry][iCentrality][iFlow]->SetPoint(iPoint, xPoint1, yPoint1-manualCorrection[iCentrality][iPoint]);
                yError1 = flowGraphJet[0][iAsymmetry][iCentrality][iFlow]->GetErrorY(iPoint);
                flowGraphJet[1][iAsymmetry][iCentrality][iFlow]->SetPointError(iPoint, 0, yError1);
              }
            }
            
            for(int iFile = 0; iFile < nComparisonFiles+1; iFile++){
              flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iFile]);
              flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(fileColors[iFile]);
              flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->SetMarkerSize(1.3);
              if(iFile == 0){
                 drawer->DrawGraph(flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow], 0, maxTrackPt, minZoomJetVn[iCentrality], maxZoomJetVn[iCentrality], "Track p_{T} (GeV)", namerY, " ", "p");
              } else {
                flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
                
                // Print the jet vn point values
                //for(int iPoint = 0; iPoint < flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->GetN(); iPoint++){
                //  flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->GetPoint(iPoint, xPoint1, yPoint1);
                //  cout << yPoint1 << endl;
                //}
                
              }
              
              // Fit a constant line to the jet v_{n} values
              if(fitJetVn){
                if((iCentrality == 2 && iFlow == 2) || (iCentrality == 0 && iFlow == 1)){
                  flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->Fit("pol0","0","",0,3); // vnfitrange
                } else {
                  flowGraphJet[iFile][iAsymmetry][iCentrality][iFlow]->Fit("pol0","0","",0,4); // vnfitrange
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
            
            // Point by point ratio
            if(drawRatios){
              
              legend = new TLegend(0.23,0.19,0.53,0.39); // 0.23,0.19,0.53,0.39 // 0.23,0.7,0.53,0.9
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->SetHeader(Form("Cent: %.0f-%.0f%%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
              
              // Select ratio either to the first graph or to the previous on the list
              maxFileIndex = ratioToPrevious ? (nComparisonFiles+1)/2 : nComparisonFiles;
              
              for(int iFile = 0; iFile < maxFileIndex; iFile++){
                for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
                  
                  baseIndex = ratioToPrevious ? iFile*2 : 0;
                  comparisonIndex = ratioToPrevious ? iFile*2+1 : iFile+1;
                  
                  flowGraphJet[baseIndex][iAsymmetry][iCentrality][iFlow]->GetPoint(iPoint, xPoint1, yPoint1);
                  yError1 = flowGraphJet[baseIndex][iAsymmetry][iCentrality][iFlow]->GetErrorY(iPoint);
                  flowGraphJet[comparisonIndex][iAsymmetry][iCentrality][iFlow]->GetPoint(iPoint, xPoint2, yPoint2);
                  yError2 = flowGraphJet[comparisonIndex][iAsymmetry][iCentrality][iFlow]->GetErrorY(iPoint);
                  
                  ratioValue = yPoint2 / yPoint1;
                  combinedError = TMath::Sqrt(TMath::Power(yError2/yPoint1,2) + TMath::Power((yPoint2*yError1)/(yPoint1*yPoint1),2));
                  
                  ratioGraphJet[iFile][iAsymmetry][iCentrality]->SetPoint(iPoint,xPoint1,ratioValue);
                  ratioGraphJet[iFile][iAsymmetry][iCentrality]->SetPointError(iPoint, 0, combinedError);
                  
                } // Point loop
                
                ratioGraphJet[iFile][iAsymmetry][iCentrality]->SetMarkerStyle(fullMarkers[ratioToPrevious ? baseIndex : comparisonIndex]);
                ratioGraphJet[iFile][iAsymmetry][iCentrality]->SetMarkerSize(1.3);
                ratioGraphJet[iFile][iAsymmetry][iCentrality]->SetMarkerColor(fileColors[ratioToPrevious ? baseIndex : comparisonIndex]);
                
                legend->AddEntry(ratioGraphJet[iFile][iAsymmetry][iCentrality], Form("%s / %s", fileLegend[comparisonIndex].Data(), fileLegend[baseIndex].Data()), "p");
                
                if(iFile == 0){
                  drawer->DrawGraph(ratioGraphJet[iFile][iAsymmetry][iCentrality], 0, maxTrackPt, 0, 3, "Track p_{T} (GeV)", "Ratio", " ", "p");
                } else {
                  ratioGraphJet[iFile][iAsymmetry][iCentrality]->Draw("p,same");
                }
              } // File loop
              
              oneLine->Draw();
              legend->Draw();
              
              // Save the figures to file
              if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/jetV%dRatio%s%s_C=%.0f-%.0f.pdf", iFlow+1, saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
              }
              
            }
            
            // Plots as a function of Q-vector
            if(drawQvectorTrends){
              legend = new TLegend(0.2,0.65,0.5,0.9);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->SetHeader(Form(Form("Cent: %%.0f-%%.0f%%%%%%s, %%.%df%%%% shifted%%s", centralityAccuracy), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data(), labelShiftNumber, subeNon0String));
              
              for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
                flowGraphQvectorJet[iPoint][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iPoint]);
                flowGraphQvectorJet[iPoint][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(fileColors[iPoint]);
                flowGraphQvectorJet[iPoint][iAsymmetry][iCentrality][iFlow]->SetMarkerSize(1.3);
                flowGraphQvectorJet[iPoint][iAsymmetry][iCentrality][iFlow]->Fit("pol1");
                flowGraphQvectorJet[iPoint][iAsymmetry][iCentrality][iFlow]->GetFunction("pol1")->SetLineColor(fileColors[iPoint]);
                if(iPoint == 0){
                  drawer->DrawGraph(flowGraphQvectorJet[iPoint][iAsymmetry][iCentrality][iFlow], qVectorXMin, qVectorXMax, minZoomJetVn[iCentrality], maxZoomJetVn[iCentrality], Form("Q-vector %s", qVectorType), namerY, " ", "p");
                } else {
                  flowGraphQvectorJet[iPoint][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
                }
                
                // Check the correction value on the obtained Q-value
                //crossingPoint = flowGraphQvectorJet[iPoint][iAsymmetry][iCentrality][iFlow]->GetFunction("pol1")->Eval(fineTunedQ[iPoint][iAsymmetry][iCentrality][iFlow]);
                crossingPoint = flowGraphQvectorJet[iPoint][iAsymmetry][iCentrality][iFlow]->GetFunction("pol1")->Eval(5.115); // jetfit
                
                lowPoint = iPoint;
                if(cumulativeBins) lowPoint = 0;
                legend->AddEntry(flowGraphQvectorJet[iPoint][iAsymmetry][iCentrality][iFlow], Form("%.1f < p_{T} < %.1f GeV, Corr=%.4f", trackPtBinBorders[lowPoint], trackPtBinBorders[iPoint+1], crossingPoint), "p");
                
              }
              
              legend->Draw();
            }
            
          } // File comparison for jet vn
          
          if(drawHadronVnVersusJetVn){
            
            legend = new TLegend(0.2,0.65,0.5,0.9);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->SetHeader(Form(Form("Cent: %%.0f-%%.0f%%%%%%s, %%.%df%%%% shifted%%s", centralityAccuracy), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data(), labelShiftNumber, subeNon0String));
            
            for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
              flowGraphHadronVsJet[iPoint][iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iPoint]);
              flowGraphHadronVsJet[iPoint][iAsymmetry][iCentrality][iFlow]->SetMarkerColor(fileColors[iPoint]);
              flowGraphHadronVsJet[iPoint][iAsymmetry][iCentrality][iFlow]->SetMarkerSize(1.3);
              flowGraphHadronVsJet[iPoint][iAsymmetry][iCentrality][iFlow]->Fit("pol1");
              flowGraphHadronVsJet[iPoint][iAsymmetry][iCentrality][iFlow]->GetFunction("pol1")->SetLineColor(fileColors[iPoint]);
              if(iPoint == 0){
                drawer->DrawGraph(flowGraphHadronVsJet[iPoint][iAsymmetry][iCentrality][iFlow], 0, 0.3, minZoomJetVn[iCentrality], maxZoomJetVn[iCentrality], "Hadron v_{2}", "Jet v_{2}", " ", "p");
              } else {
                flowGraphHadronVsJet[iPoint][iAsymmetry][iCentrality][iFlow]->Draw("p,same");
              }
                            
              
              flowGraphHadron[0][iAsymmetry][iCentrality][iFlow]->GetPoint(iPoint, xValueForQ, yValueForQ);
              
              // Check the correction value on the matching hadron v2 in data
              crossingPoint = flowGraphHadronVsJet[iPoint][iAsymmetry][iCentrality][iFlow]->GetFunction("pol1")->Eval(yValueForQ);
              
              lowPoint = iPoint;
              if(cumulativeBins) lowPoint = 0;
              legend->AddEntry(flowGraphHadronVsJet[iPoint][iAsymmetry][iCentrality][iFlow], Form("%.1f < p_{T} < %.1f GeV, C = %.4f", trackPtBinBorders[lowPoint], trackPtBinBorders[iPoint+1], crossingPoint), "p");
              
              dataLine[iPoint][iAsymmetry][iCentrality][iFlow] = new TLine(yValueForQ, 0.03, yValueForQ, 0.13);
              dataLine[iPoint][iAsymmetry][iCentrality][iFlow]->SetLineColor(fileColors[iPoint]);
              dataLine[iPoint][iAsymmetry][iCentrality][iFlow]->SetLineStyle(2);
              dataLine[iPoint][iAsymmetry][iCentrality][iFlow]->Draw();;
              
            }
            
            legend->Draw();
            
//            vLegend = new TLegend(0.19,0.18,0.89,0.28);
//            vLegend->SetFillStyle(0);vLegend->SetBorderSize(0);vLegend->SetTextSize(0.04);vLegend->SetTextFont(62);
//            vLegend->SetHeader("Solid line = fit to MC, Dashed lines = data");
//            vLegend->Draw();
            
          } // Hadron vn versus jet vn
          
        } // Asymmetry loop
      } // Flow component loop
    } // Centrality loop
    
    // File comparison for yields
    if(drawJetHadronYieldFileComparison || drawDihadronYieldFileComparison){
      
      double maxJetHadronYield[] = {3000, 2000, 1000, 500};
      double maxDihadronYield[] = {8000000, 2000000, 400000, 100000}; // 5500000
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
                 drawer->DrawGraph(yieldGraphJetHadron[iFile][iAsymmetry][iCentrality], 0, maxTrackPt, 1, maxJetHadronYield[iCentrality], "Track p_{T} (GeV)", "Jet-hadron yield", " ", "p");
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
            
            // Point by point ratio
            if(drawRatios){
              
              legend = new TLegend(0.23,0.75,0.53,0.9); // 0.23,0.19,0.53,0.39  //  0.23,0.7,0.53,0.9
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->SetHeader(Form("Cent: %.0f-%.0f%%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
              
              // Select ratio either to the first graph or to the previous on the list
              maxFileIndex = ratioToPrevious ? (nComparisonFiles+1)/2 : nComparisonFiles;
              
              for(int iFile = 0; iFile < maxFileIndex; iFile++){
                for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
                  
                  baseIndex = ratioToPrevious ? iFile*2 : 0;
                  comparisonIndex = ratioToPrevious ? iFile*2+1 : iFile+1;
                  
                  yieldGraphJetHadron[baseIndex][iAsymmetry][iCentrality]->GetPoint(iPoint, xPoint1, yPoint1);
                  yError1 = yieldGraphJetHadron[baseIndex][iAsymmetry][iCentrality]->GetErrorY(iPoint);
                  yieldGraphJetHadron[comparisonIndex][iAsymmetry][iCentrality]->GetPoint(iPoint, xPoint2, yPoint2);
                  yError2 = yieldGraphJetHadron[comparisonIndex][iAsymmetry][iCentrality]->GetErrorY(iPoint);
                  
                  ratioValue = yPoint2 / yPoint1;
                  combinedError = TMath::Sqrt(TMath::Power(yError2/yPoint1,2) + TMath::Power((yPoint2*yError1)/(yPoint1*yPoint1),2));
                  
                  ratioGraphJetHadronYield[iFile][iAsymmetry][iCentrality]->SetPoint(iPoint,xPoint1,ratioValue);
                  ratioGraphJetHadronYield[iFile][iAsymmetry][iCentrality]->SetPointError(iPoint, 0, combinedError);
                  cout << ratioValue << endl;
                  
                } // Point loop
                
                ratioGraphJetHadronYield[iFile][iAsymmetry][iCentrality]->SetMarkerStyle(fullMarkers[ratioToPrevious ? baseIndex : comparisonIndex]);
                ratioGraphJetHadronYield[iFile][iAsymmetry][iCentrality]->SetMarkerSize(1.3);
                ratioGraphJetHadronYield[iFile][iAsymmetry][iCentrality]->SetMarkerColor(fileColors[ratioToPrevious ? baseIndex : comparisonIndex]);
                
                legend->AddEntry(ratioGraphJetHadronYield[iFile][iAsymmetry][iCentrality], Form("%s / %s", fileLegend[comparisonIndex].Data(), fileLegend[baseIndex].Data()), "p");
                
                if(iFile == 0){
                  drawer->DrawGraph(ratioGraphJetHadronYield[iFile][iAsymmetry][iCentrality], 0, maxTrackPt, 0.99, 1.01, "Track p_{T} (GeV)", "Dihadron yield ratio", " ", "p");
                } else {
                  ratioGraphJetHadronYield[iFile][iAsymmetry][iCentrality]->Draw("p,same");
                }
              } // File loop
              
              oneLine->Draw();
              legend->Draw();
              
              // Save the figures to file
              if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/jetHadronYieldRatio%s%s_C=%.0f-%.0f.pdf", saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
              }
              
            }
            
            // Plots as a function of Q-vector
            if(drawQvectorTrends){
              legend = new TLegend(0.2,0.6,0.5,0.9);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->SetHeader(Form(Form("Cent: %%.0f-%%.0f%%%%%%s, %%.%df%%%% shifted%%s", centralityAccuracy), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data(), labelShiftNumber, subeNon0String));
              
              for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
                yieldGraphJetHadronVsQvector[iPoint][iAsymmetry][iCentrality]->SetMarkerStyle(fullMarkers[iPoint]);
                yieldGraphJetHadronVsQvector[iPoint][iAsymmetry][iCentrality]->SetMarkerColor(fileColors[iPoint]);
                yieldGraphJetHadronVsQvector[iPoint][iAsymmetry][iCentrality]->SetMarkerSize(1.3);
                yieldGraphJetHadronVsQvector[iPoint][iAsymmetry][iCentrality]->Fit("pol1");
                yieldGraphJetHadronVsQvector[iPoint][iAsymmetry][iCentrality]->GetFunction("pol1")->SetLineColor(fileColors[iPoint]);
                if(iPoint == 0){
                  drawer->DrawGraph(yieldGraphJetHadronVsQvector[iPoint][iAsymmetry][iCentrality], qVectorXMin, qVectorXMax, 0.2, 1.6, Form("Q-vector %s",qVectorType), "Jet-hadron yield: MC/data", " ", "p");
                } else {
                  yieldGraphJetHadronVsQvector[iPoint][iAsymmetry][iCentrality]->Draw("p,same");
                }
                                
                // Solve where the function crosses the line
                crossingPoint = findCrossingPoint(yieldGraphJetHadronVsQvector[iPoint][iAsymmetry][iCentrality]->GetFunction("pol1"), 1);
                
                lowPoint = iPoint;
                if(cumulativeBins) lowPoint = 0;
                legend->AddEntry(yieldGraphJetHadronVsQvector[iPoint][iAsymmetry][iCentrality], Form("%.1f < p_{T} < %.1f GeV, Q = %.3f", trackPtBinBorders[lowPoint], trackPtBinBorders[iPoint+1], crossingPoint), "p");
                
              } // Track pT loop
              
              oneQLine->Draw();
              legend->Draw();
              
            } // Drawing trends
            
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
                 drawer->DrawGraph(yieldGraphDihadron[iFile][iAsymmetry][iCentrality], 0, maxTrackPt, 1, maxDihadronYield[iCentrality], "Track p_{T} (GeV)", "Dihadron yield", " ", "p");
              } else {
                yieldGraphDihadron[iFile][iAsymmetry][iCentrality]->Draw("p,same");
              }
              
              legend->AddEntry(yieldGraphDihadron[iFile][iAsymmetry][iCentrality], fileLegend[iFile], "p");
              
            } // File loop
            
            legend->Draw();
            
            // Save the figures to file
            if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/dihadronYieldComparison%s%s_C=%.0f-%.0f.pdf", saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
            }
            
            // Point by point ratio
            if(drawRatios){
              
              legend = new TLegend(0.23,0.7,0.53,0.9); // 0.23,0.19,0.53,0.39  //  0.23,0.7,0.53,0.9
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->SetHeader(Form("Cent: %.0f-%.0f%%%s", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
              
              // Select ratio either to the first graph or to the previous on the list
              maxFileIndex = ratioToPrevious ? (nComparisonFiles+1)/2 : nComparisonFiles;
              
              for(int iFile = 0; iFile < maxFileIndex; iFile++){
                for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
                  
                  baseIndex = ratioToPrevious ? iFile*2 : 0;
                  comparisonIndex = ratioToPrevious ? iFile*2+1 : iFile+1;
                  
                  yieldGraphDihadron[baseIndex][iAsymmetry][iCentrality]->GetPoint(iPoint, xPoint1, yPoint1);
                  yError1 = yieldGraphDihadron[baseIndex][iAsymmetry][iCentrality]->GetErrorY(iPoint);
                  yieldGraphDihadron[comparisonIndex][iAsymmetry][iCentrality]->GetPoint(iPoint, xPoint2, yPoint2);
                  yError2 = yieldGraphDihadron[comparisonIndex][iAsymmetry][iCentrality]->GetErrorY(iPoint);
                  
                  ratioValue = yPoint2 / yPoint1;
                  combinedError = TMath::Sqrt(TMath::Power(yError2/yPoint1,2) + TMath::Power((yPoint2*yError1)/(yPoint1*yPoint1),2));
                  
                  ratioGraphDihadronYield[iFile][iAsymmetry][iCentrality]->SetPoint(iPoint,xPoint1,ratioValue);
                  ratioGraphDihadronYield[iFile][iAsymmetry][iCentrality]->SetPointError(iPoint, 0, combinedError);
                  cout << ratioValue << endl;
                  
                } // Point loop
                
                ratioGraphDihadronYield[iFile][iAsymmetry][iCentrality]->SetMarkerStyle(fullMarkers[ratioToPrevious ? baseIndex : comparisonIndex]);
                ratioGraphDihadronYield[iFile][iAsymmetry][iCentrality]->SetMarkerSize(1.3);
                ratioGraphDihadronYield[iFile][iAsymmetry][iCentrality]->SetMarkerColor(fileColors[ratioToPrevious ? baseIndex : comparisonIndex]);
                
                legend->AddEntry(ratioGraphDihadronYield[iFile][iAsymmetry][iCentrality], Form("%s / %s", fileLegend[comparisonIndex].Data(), fileLegend[baseIndex].Data()), "p");
                
                if(iFile == 0){
                  drawer->DrawGraph(ratioGraphDihadronYield[iFile][iAsymmetry][iCentrality], 0, maxTrackPt, 0.95, 1.05, "Track p_{T} (GeV)", "Dihadron yield ratio", " ", "p");
                } else {
                  ratioGraphDihadronYield[iFile][iAsymmetry][iCentrality]->Draw("p,same");
                }
              } // File loop
              
              oneLine->Draw();
              legend->Draw();
              
              // Save the figures to file
              if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/dihadronYieldRatio%s%s_C=%.0f-%.0f.pdf", saveComment.Data(), compactAsymmetryString[iAsymmetry].Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
              }
              
            }
            
            // Plots as a function of Q-vector
            if(drawQvectorTrends){
              legend = new TLegend(0.2,0.6,0.5,0.9);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->SetHeader(Form(Form("Cent: %%.0f-%%.0f%%%%%%s, %%.%df%%%% shifted%%s", centralityAccuracy), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data(), labelShiftNumber, subeNon0String));
              
              for(int iPoint = 0; iPoint < maxPtBin; iPoint++){
                yieldGraphDihadronVsQvector[iPoint][iAsymmetry][iCentrality]->SetMarkerStyle(fullMarkers[iPoint]);
                yieldGraphDihadronVsQvector[iPoint][iAsymmetry][iCentrality]->SetMarkerColor(fileColors[iPoint]);
                yieldGraphDihadronVsQvector[iPoint][iAsymmetry][iCentrality]->SetMarkerSize(1.3);
                yieldGraphDihadronVsQvector[iPoint][iAsymmetry][iCentrality]->Fit("pol1");
                yieldGraphDihadronVsQvector[iPoint][iAsymmetry][iCentrality]->GetFunction("pol1")->SetLineColor(fileColors[iPoint]);
                if(iPoint == 0){
                  drawer->DrawGraph(yieldGraphDihadronVsQvector[iPoint][iAsymmetry][iCentrality], qVectorXMin, qVectorXMax, 0.2, 1.6, Form("Q-vector %s",qVectorType), "Dihadron yield: MC/data", " ", "p");
                } else {
                  yieldGraphDihadronVsQvector[iPoint][iAsymmetry][iCentrality]->Draw("p,same");
                }
                                
                // Solve where the function crosses the line
                //crossingPoint = findCrossingPoint(yieldGraphDihadronVsQvector[iPoint][iAsymmetry][iCentrality]->GetFunction("pol1"), 1);
                crossingPoint = yieldGraphDihadronVsQvector[iPoint][iAsymmetry][iCentrality]->GetFunction("pol1")->Eval(5.115); // dihadronfit
                
                lowPoint = iPoint;
                if(cumulativeBins) lowPoint = 0;
                legend->AddEntry(yieldGraphDihadronVsQvector[iPoint][iAsymmetry][iCentrality], Form("%.1f < p_{T} < %.1f GeV, Scale = %.3f", trackPtBinBorders[lowPoint], trackPtBinBorders[iPoint+1], crossingPoint), "p");
                
              } // Track pT loop
              
              oneQLine->Draw();
              legend->Draw();
              
            } // Drawing trends
            
          } // File comparison for dihadron yields
          
        } // Asymmetry loop
      } // Centrality loop
    } // File comparison for yields
    
    drawer->SetTitleOffsetY(1.1);
    
  } // Draw file comparison plots
  
  // If we draw jet vn file comparison and fit the vn plots, we can also draw summary plots
  if((drawJetVnFileComparison && fitJetVn) && drawSummaryPlot){
    
    // For summary correction it is assumed that first the data is given and then MC.
    if(doSummaryCorrection){
      for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
          for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
            for(int iFile = 1; iFile < nComparisonFiles+1; iFile++){
              if(!manualSummaryCorrection){
                summaryYaxis[iFile][iAsymmetry][iFlow][iCentrality] = summaryYaxis[0][iAsymmetry][iFlow][iCentrality] - summaryYaxis[iFile][iAsymmetry][iFlow][iCentrality];
                summaryYaxisError[iFile][iAsymmetry][iFlow][iCentrality] = TMath::Sqrt(TMath::Power(summaryYaxisError[0][iAsymmetry][iFlow][iCentrality],2) + TMath::Power(summaryYaxisError[iFile][iAsymmetry][iFlow][iCentrality],2));
              }
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
    atlasJetV2graph->SetMarkerSize(1.8);
    
    cmsHighPtV2 = new TGraphErrors(nCentralityBins, summaryXaxis, cmsHighPtV2Number, summaryXaxisError, cmsHighPtV2Error);
    cmsHighPtV2->SetMarkerStyle(kFullStar);
    cmsHighPtV2->SetMarkerColor(kAzure+9);
    cmsHighPtV2->SetMarkerSize(1.8);
    
    flowJetReferenceGraph = new TGraphErrors(nCentralityBins, summaryXaxis, flowJetV2Number, summaryXaxisError, flowJetV2Error);
    flowJetReferenceGraph->SetMarkerStyle(kFullCircle);
    flowJetReferenceGraph->SetMarkerColor(kRed);
    flowJetReferenceGraph->SetMarkerSize(1.8);
    
    pfCsJetReferenceGraph = new TGraphErrors(nCentralityBins, summaryXaxis, pfCsJetV2Number, summaryXaxisError, flowJetV2Error);
    pfCsJetReferenceGraph->SetMarkerStyle(kFullFourTrianglesPlus);
    pfCsJetReferenceGraph->SetMarkerColor(kGreen+3);
    pfCsJetReferenceGraph->SetMarkerSize(1.8);
    
    doubleDijetFlowJetReferenceGraph = new TGraphErrors(nCentralityBins, summaryXaxis, doubleDijetFlowJetNumber, summaryXaxisError, doubleDijetFlowJetError);
    doubleDijetFlowJetReferenceGraph->SetMarkerStyle(kFullCrossX);
    doubleDijetFlowJetReferenceGraph->SetMarkerColor(kMagenta);
    doubleDijetFlowJetReferenceGraph->SetMarkerSize(1.8);
    
    doubleDijetCaloJetReferenceGraph = new TGraphErrors(nCentralityBins, summaryXaxis, doubleDijetCaloJetNumber, summaryXaxisError, doubleDijetCaloJetError);
    doubleDijetCaloJetReferenceGraph->SetMarkerStyle(kFullFourTrianglesPlus);
    doubleDijetCaloJetReferenceGraph->SetMarkerColor(kGreen+3);
    doubleDijetCaloJetReferenceGraph->SetMarkerSize(1.8);
    
    for(int iFile = 0; iFile < nComparisonFiles+1; iFile++){
      for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
          flowSummaryJet[iFile][iAsymmetry][iFlow] = new TGraphErrors(nCentralityBins, summaryXaxis, summaryYaxis[iFile][iAsymmetry][iFlow], summaryXaxisError, summaryYaxisError[iFile][iAsymmetry][iFlow]);
          flowSummaryJet[iFile][iAsymmetry][iFlow]->SetMarkerStyle(secondMarkers[iFile]);
          flowSummaryJet[iFile][iAsymmetry][iFlow]->SetMarkerColor(fileColors[iFile]);
          flowSummaryJet[iFile][iAsymmetry][iFlow]->SetMarkerSize(1.8);
          
          // Set the bin labels for x-axis
          for(int iCentrality = 0; iCentrality < nCentralityBins*2; iCentrality++){
            flowSummaryJet[iFile][iAsymmetry][iFlow]->GetXaxis()->ChangeLabel(iCentrality+1,-1,-1,-1,-1,-1,binLabels[iCentrality]);
          } // Centrality loop
        } // Asymmetry loop
      } // Flow component loop
    } // File loop
    
    double minSummary[] = {0,0,-0.03,0};
    double maxSummary[] = {0.1,0.1,0.03,0.1};
    
    // Once the graphs are constructed, they can be plotted
    for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
        
        legend = new TLegend(0.2,0.7,0.5,0.9); //0.2,0.6,0.5,0.9  // Grant x-axis: 0.168 0.468
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        
        for(int iFile = 1; iFile < nComparisonFiles; iFile++){
          if(iFile == 1){
            sprintf(namerY,"Jet v_{%d}",iFlow+1);
            drawer->DrawGraphCustomAxes(flowSummaryJet[iFile][iAsymmetry][iFlow], 0, 4, minSummary[iFlow], maxSummary[iFlow], "Centrality", namerY, " ", "ap"); // 0, 4, -0.05, 0.3
          } else {
            flowSummaryJet[iFile][iAsymmetry][iFlow]->Draw("p,same");
          }
          legend->AddEntry(flowSummaryJet[iFile][iAsymmetry][iFlow], fileLegend[iFile], "p");
        } // File loop
        
        if(iFlow == 1 && iAsymmetry == nAsymmetryBins){
          if(drawFlowJetResults){
            flowJetReferenceGraph->Draw("p,same");
            legend->AddEntry(flowJetReferenceGraph, "Flow jets", "p");
          }
          
          if(drawPfCsJetResults){
            pfCsJetReferenceGraph->Draw("p,same");
            legend->AddEntry(pfCsJetReferenceGraph, "PfCs jets", "p");
          }
          
          if(drawCaloDoubleDijetResults){
            doubleDijetCaloJetReferenceGraph->Draw("p,same");
            legend->AddEntry(doubleDijetCaloJetReferenceGraph, "Calo jets with flow dijet", "p");
          }
          
          if(drawFlowDoubleDijetResults){
            doubleDijetFlowJetReferenceGraph->Draw("p,same");
            legend->AddEntry(doubleDijetFlowJetReferenceGraph, "Flow jets with calo dijet", "p");
          }
          
          if(drawPreviousResults){
            atlasJetV2graph->Draw("p,same");
            legend->AddEntry(atlasJetV2graph, "ATLAS v_{2}", "p"); // (#scale[0.8]{HP 2020})
            
            cmsHighPtV2->Draw("p,same");
            legend->AddEntry(cmsHighPtV2, "CMS high p_{T} v_{2}", "p"); // (#scale[0.8]{PLB 776 (2018) 195})
          }
        }
        
        shortZeroLine->Draw();
        legend->Draw();
        
        // Save the figures to file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/jetV%dSummary%s%s.pdf", iFlow+1, saveComment.Data(), compactAsymmetryString[iAsymmetry].Data()));
        }
      } // Asymmetry loop
    } // Flow component loop
    
    // Save the summary plots so that they can be used with the final graph plotter
    if(saveSummaryFile){
      
      // Create the output file
      TFile *outputFile = new TFile(outputFileName,"UPDATE");
      char histogramNamer[100];
      
      for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
          flowSummaryJet[1][iAsymmetry][iFlow]->Write(Form("summaryV%d%s", iFlow+1, compactAsymmetryString[iAsymmetry].Data()), TObject::kOverwrite);
        }
      }
      
      outputFile->Close();
    } // Saving the summary plots
    
  }
  
}

