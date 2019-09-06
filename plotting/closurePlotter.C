#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

//Define types for histogram arrays
enum enumCollisionSystem {kPp, kPbPb, knCollisionSystems};
enum enumDataType {kData, kMC, knDataTypes};
enum enumMonteCarloType{kRecoReco, kRecoGen, kGenReco, kGenGen, knMonteCarloTypes};

/*
 * Our stylist will give your histogram a fresh new look!
 *
 *  TH1* histogram = Histogram in need of a stylist
 *  int markerStyle = Marker style to be set for the histogram
 *  int markerColor = Marker color to be set for the histogram
 *  int rebin = Rebin to be applied for the histogram
 */
void stylistForHistogram(TH1* histogram, int markerStyle, int markerColor, int rebin){
  
  // Set the merker style and color
  histogram->SetMarkerStyle(markerStyle);
  histogram->SetMarkerColor(markerColor);
  histogram->SetLineColor(markerColor);
  
  // Do the rebin an scaling to retain the normalization
  if(rebin > 1){
    histogram->Rebin(rebin);
    histogram->Scale(1.0/rebin);
  }
}

/*
 * Normalize a given histogram to one and set a nice style for it
 *
 *  TH1 *histogram = Histogram in need of a stylist
 *  const int iDataType = Data type index, used by the stylist
 */
void setStyleAndNormalize(TH1 *histogram, int iDataType){
  
  // Transform the given data type to allowed interval
  if(iDataType < 0) {
    iDataType = 0;
    cout << "Negative data type not allowed! Setting it to 0 for styling purposes." << endl;
  }
  if(iDataType > 2){
    iDataType = 2;
    cout << "Our stylist cannot handle data types larger than 2! Using 2 as data type for styling purposes." << endl;
  }
  
  // Common style settings for all the figures
  int drawColor[3] = {kBlack,kRed,kBlue};
  int goodStyle[3] = {21,20,20};
  
  // Normalize the histogram to 1
  double integral = histogram->Integral();
  histogram->Scale(1.0/integral);
  
  // Put the histogram in the bench of a stylist
  histogram->SetMarkerStyle(goodStyle[iDataType]);
  histogram->SetMarkerColor(drawColor[iDataType]);
}

/*
 * Plotter for track closure histograms
 *
 *  JDrawer *drawer = Drawer doing all the dirty work for drawing the histograms
 *  TH1 *genHistogram = Histogram with generator level tracks
 *  TH1 *recoHistogram = Histogram with corrected reconstructed tracks
 *  TH1 *uncorrectedHistogram = Histogram with uncorrected reconstructed tracks
 *  const char *xTitle = Title given for the x-axis of the histograms
 *  const char *yTitle = Title given for the y-axis of the histograms
 *  bool logScale = Use logarithmic scale for the main histogram
 *  int rebin = Rebin for the histograms
 *  double yZoom = Range for the y-axis in the main histogram
 *  bool bigZoom = Zoom to very close range to show exactly the deviation from one
 *  legendX = Left side x-position of the title
 *  legendY = Bottom side y-position of the title
 *  const char *header = Header for the legend
 *  const char *centralityString = String for the centrality
 *  const char *trackPtString = String for track pT
 */
void plotTrackClosure(JDrawer *drawer, TH1 *genHistogram, TH1 *recoHistogram, TH1 *uncorrectedHistogram, const char *xTitle, const char *yTitle, bool logScale, int rebin, double yZoom, bool bigZoom, double legendX, double legendY, const char *header, const char *centralityString, const char *trackPtString){

  // Helper variable to name the histograms
  char namer[200];
  
  // Style settings for the markers
  int markerStyleGen = 21;     // Full color square marker
  int markerColorGen = 38;     // Grey-blue
  int markerStyleReco = 20;    // Full color round masker
  int markerColorCorr = kRed;  // Red
  int markerColorReco = kBlue; // Blue
 
  stylistForHistogram(genHistogram,markerStyleGen,markerColorGen,rebin);
  genHistogram->SetMarkerSize(1.2);
  
  drawer->SetLogY(logScale);
  drawer->CreateSplitCanvas();
  if(yZoom > 0) genHistogram->GetYaxis()->SetRangeUser(0,yZoom);
  drawer->DrawHistogramToUpperPad(genHistogram,xTitle,yTitle," ");

  double ySpace = 0;
  if(strncmp(centralityString,"",2) != 0) ySpace += 0.025;
  if(strncmp(trackPtString,"",2) != 0) ySpace += 0.025;
  
  stylistForHistogram(recoHistogram,markerStyleReco,markerColorCorr,rebin);
  recoHistogram->Draw("same");
  stylistForHistogram(uncorrectedHistogram,markerStyleReco,markerColorReco,rebin);
  if(!bigZoom)  uncorrectedHistogram->Draw("same");

  TLegend *legend = new TLegend(legendX,legendY-ySpace,legendX+0.3,legendY+0.25+ySpace);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

  legend->SetHeader(header);
  if(strncmp(centralityString,"",2) != 0) legend->AddEntry((TObject*) 0,centralityString,"");
  if(strncmp(trackPtString,"",2) != 0) legend->AddEntry((TObject*) 0,trackPtString,"");
  legend->AddEntry(genHistogram,"Gen. particles","pl");
  legend->AddEntry(recoHistogram,"Corr. tracks","pl");
  if(!bigZoom) legend->AddEntry(uncorrectedHistogram,"Reco. tracks","pl");
  legend->Draw();

  sprintf(namer,"%sRatio",recoHistogram->GetName());
  TH1D *recoRatio = (TH1D*) recoHistogram->Clone(namer);
  recoRatio->Divide(genHistogram);
  
  sprintf(namer,"%sRatio",uncorrectedHistogram->GetName());
  TH1D *uncorrectedRatio = (TH1D*) uncorrectedHistogram->Clone(namer);
  uncorrectedRatio->Divide(genHistogram);
  
  drawer->SetLogY(false);
  if(bigZoom){
    recoRatio->GetYaxis()->SetRangeUser(0.9,1.1);
  } else {
    recoRatio->GetYaxis()->SetRangeUser(0,2);
  }
  drawer->DrawHistogramToLowerPad(recoRatio,xTitle,"Reco/Gen"," ");
  if(!bigZoom) uncorrectedRatio->Draw("same");

}

/*
 * Plotter for jet kinematics comparison between data and simulations
 *
 *  JDrawer *drawer = Drawer doing all the dirty work for drawing the histograms
 *  TH1 *dataHistogram = Histogram with jet kinematics from real data
 *  TH1 *mcHistogram = Histogram with jet kinematics from simulation
 *  const char *xTitle = Title given for the x-axis of the histograms
 *  const char *yTitle = Title given for the y-axis of the histograms
 *  bool logScale = Use logarithmic scale for the main histogram
 *  int rebin = Rebin for the histograms
 *  double xRangeMin = Start for the x-range in the main histogram
 *  double yZoom = Range for the y-axis in the main histogram
 *  legendX = Left side x-position of the title
 *  legendY = Bottom side y-position of the title
 *  const char *header = Header for the legend
 *  const char *centralityString = String for the centrality
 *  const char *dataString = String describing data in the legend
 *  const char *mcString = String decribing the simulation in the legend
 *  int markerStyleData = Marker style to be used with the data histogram
 *  int markerColorData = Marker color to be used with the data histogram
 *  int markerStyleMC = Marker style to be used with the simulation histogram
 *  int markerColorMC = Marker color to be used with the simulation histogram
 */
void plotJetKinematics(JDrawer *drawer, TH1 *dataHistogram, TH1 *mcHistogram, const char *xTitle, const char *yTitle, bool logScale, int rebin, double xRangeMin, double yZoom, double legendX, double legendY, const char *header, const char *centralityString, const char *dataString, const char *mcString, int markerStyleData, int markerColorData, int markerStyleMC, int markerColorMC){
  
  // Helper variable for naming
  char namer[100];
  
  // Set up the histogram style
  stylistForHistogram(dataHistogram,markerStyleData,markerColorData,rebin);
  stylistForHistogram(mcHistogram,markerStyleMC,markerColorMC,rebin);
  
  drawer->SetLogY(logScale);
  double minY = logScale ? 8e-5 : 0;
  if(xRangeMin > -100) mcHistogram->GetXaxis()->SetRangeUser(xRangeMin,300);
  if(yZoom > 0) mcHistogram->GetYaxis()->SetRangeUser(minY,yZoom);
  
  // Special case for inclusive pT plots to force the x-axis range from 50 to 300
  if(logScale && (xRangeMin > -100) && (yZoom > 0)){
    drawer->CreateCanvas(xRangeMin,300,minY,yZoom,xTitle,yTitle," ");
    mcHistogram->Draw("same");
  } else {
    drawer->DrawHistogram(mcHistogram,xTitle,yTitle," ");
  }
  
  double ySpace = 0;
  if(strncmp(centralityString,"",2) != 0) ySpace += 0.025;
  
  dataHistogram->Draw("same");
  
  TLegend *legend = new TLegend(legendX,legendY-ySpace,legendX+0.3,legendY+0.15+ySpace);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  
  sprintf(namer,"%s jet",header);
  legend->SetHeader(namer);
  if(strncmp(centralityString,"",2) != 0) legend->AddEntry((TObject*) 0,centralityString,"");
  legend->AddEntry(dataHistogram,dataString,"pl");
  legend->AddEntry(mcHistogram,mcString,"pl");
  legend->Draw();

  drawer->SetLogY(false);
  
}

/*
 * Plotter for jet kinematics histograms in asymmetry bins
 *
 *  JDrawer *drawer = Drawer doing all the dirty work for drawing the histograms
 *  TH1 *allLeadingHistogram = Histogram with all leading jets
 *  TH1 *allSubleadingHistogram = Histogram with all subleading jets
 *  TH1 *leadingAsymmetryHistogram = Histograms with leading jets in a selected asymmetry bin
 *  TH1* subleadingAsymmetryHistogram = Histograms with subleading jets in a selected asymmetry bin
 *  const char *xTitle = Title given for the x-axis of the histograms
 *  const char *yTitle = Title given for the y-axis of the histograms
 *  bool logScale = Use logarithmic scale for the main histogram
 *  int rebin = Rebin for the histograms
 *  double xRangeMin = Start for the x-range in the main histogram
 *  double yZoom = Range for the y-axis in the main histogram
 *  bool drawLegend = Draw a legend to the figure
 *  legendX = Left side x-position of the title
 *  legendY = Bottom side y-position of the title
 *  const char *header = System to be put to legend header
 *  const char *centralityString = String for the centrality
 *  const char *asymmetryString = String for the asymmetry
 */
void plotJetKinematicsAsymmetry(JDrawer *drawer, TH1 *allLeadingHistogram, TH1 *allSubleadingHistogram, TH1 *leadingAsymmetryHistogram, TH1 *subleadingAsymmetryHistogram, const char *xTitle, const char *yTitle, bool logScale, int rebin, double xRangeMin, double yZoom, bool drawLegend, double legendX, double legendY, const char *header, const char *centralityString, const char *asymmetryString){
  
  // Helper variable to name the histograms
  char namer[200];
  
  // Set up the histogram style
  int markerStyle = 20;
  stylistForHistogram(allLeadingHistogram,markerStyle,kRed,rebin);
  stylistForHistogram(allSubleadingHistogram,markerStyle,kGreen+3,rebin);
  stylistForHistogram(leadingAsymmetryHistogram,markerStyle,kViolet-6,rebin);
  stylistForHistogram(subleadingAsymmetryHistogram,markerStyle,kBlue,rebin);
  
  // Set drawing styles and draw the first histogram
  drawer->SetLogY(logScale);
  double minY = logScale ? 8e-5 : 0;
  if(xRangeMin > -100) allLeadingHistogram->GetXaxis()->SetRangeUser(xRangeMin,300);
  if(yZoom > 0) allLeadingHistogram->GetYaxis()->SetRangeUser(minY,yZoom);
  drawer->DrawHistogram(allLeadingHistogram,xTitle,yTitle," ");
  
  // Draw the rest of the histograms to the same figure
  allSubleadingHistogram->Draw("same");
  leadingAsymmetryHistogram->Draw("same");
  subleadingAsymmetryHistogram->Draw("same");
  
  // Draw a legend if specified
  if(drawLegend){
    
    TLegend *legend = new TLegend(legendX,legendY,legendX+0.3,legendY+0.3);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    
    sprintf(namer,"%s",header);
    if(strncmp(centralityString,"",2) != 0) sprintf(namer,"%s %s",header,centralityString);
    legend->SetHeader(namer);
    legend->AddEntry(allLeadingHistogram,"All leading jets","pl");
    legend->AddEntry(allSubleadingHistogram,"All subleading jets","pl");
    sprintf(namer,"Leading %s",asymmetryString);
    legend->AddEntry(leadingAsymmetryHistogram,namer,"pl");
    sprintf(namer,"Subleading %s",asymmetryString);
    legend->AddEntry(subleadingAsymmetryHistogram,namer,"pl");
    legend->Draw();
    
  }
  
  drawer->SetLogY(false);
  
}

/*
 * Plotter for closure histograms to the analysis note
 */
void closurePlotter(){
  
  ///////////////////
  // Configuration //
  ///////////////////
  
  bool saveFigures = true;          // Save the figures to a file
  
  bool drawCentrality = false;        // Draw the QA plots for spillover correction
  bool drawVz = false;                // Draw the QA plots for seagull correction
  bool drawTrackClosure = false;       // Draw the tracking closures
  
  bool drawJetKinematicsMcComparison = true;     // Draw the jet kinematics figures comparing data and simulation
  bool drawJetKinematicsAsymmerty = false;        // Draw the jet kinematics figures in different asymmetry bins
  bool drawJetKinematics = (drawJetKinematicsAsymmerty || drawJetKinematicsMcComparison);
  
  int ptRebin = 10;                  // Rebin for track pT closure histograms (there are 500 bins)
  int trackAngleRebin = 2;          // Rebin for track eta and phi histograms
  
  // Read the number of bins from histogram manager
  const int nCentralityBins = 4; // Expected number of centrality bins
  const int nTrackPtBins = 7; // Expected number of ttack pT bins
  const int nAsymmetryBins = 4; // Maximum number of asymmetry bins that can be in a file
  double centralityBinBorders[] = {0,10,30,50,90};       // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};       // Bin borders for track pT
  
  // Normalize to high pT
  bool normalizeToHighPt = false;  // If the trigger is inefficient, we can try to normalize only to high pT to look at the shape there
  
  // Select which bins to plot
  int firstCentralityBin = 0;
  int lastCentralityBin = nCentralityBins-1;
  
  int firstTrackPtBin = 0;
  int lastTrackPtBin = nTrackPtBins-1;
  
  // Zoom very close to one in closure plots
  bool bigClosureZoom = false;
  
  /////////////////
  // Config done //
  /////////////////
  
  const char *legendNames[knCollisionSystems][knDataTypes+1] = {{"pp","Raw Pythia", "Weighted Pythia"},{"PbPb","Raw P+H", "Weighted P+H"}};
  const char *monteCarloName[knCollisionSystems] = {"Pythia8","Pythia+Hydjet"};
  const char *inclusiveLegend[2] = {"dijet","inclusive"};
  TString systemString[knCollisionSystems] = {"Pp","PbPb"};
  
  // Open files containing the QA histograms
  DijetHistogramManager *inputManager[knCollisionSystems][knDataTypes];
  TFile *inputFile[knCollisionSystems][knDataTypes];
  inputFile[kPp][kData] = TFile::Open("data/ppData2017_highForest_pfJets_noMixing_onlyJets_JECv3_wtaAxis_processed_2019-08-30.root");
  // data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_wtaAxis_allCorrections_noSmoothing_tightSideBand_processed_2019-08-13.root
  // data/ppData2017_highForest_pfJets_onlyJets_wtaAxis_processed_2019-08-05.root
  // data/dijet_pp_highForest_pfJets_noUncOrInc_allCorrections_wtaAxis_processed_2019-07-13.root
  inputFile[kPbPb][kData] = TFile::Open("data/dijetPbPb2018_highForest_akFlowPu4CsPFJets_JECv5b_onlyJets_processed_2019-08-28.root");
  // data/dijetPbPb2018_highForest_akFlowPu4CsPFJets_JECv5b_onlyJets_processed_2019-08-28.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_onlyJets_onlyL2RelV4_wtaAxis_processed_2019-08-13.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_onlyJets_rawPt_wtaAxis_processed_2019-08-05.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_onlyJets_wtaAxis_processed_2019-08-13.root
  // data/dijetPbPb2018_highForest_flowPuCs4PfJets_noUncorr_wtaAxis_onlyJets_processed_2019-08-02.root
  // data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root
  // data/dijetPbPb_pfCsJets_xj_noCorrelations_jetCountWTA_processed_2019-07-01.root
  // data/dijetPbPb_skims_pfJets_noUncorrected_10mixedEvents_smoothedMixing_noCorrections_processed_2019-01-07.root
  inputFile[kPp][kMC] = TFile::Open("data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_onlyJets_JECv3_processed_2019-08-30.root");
  // data/dijet_ppMC_RecoReco_Pythia8_pfJets_wtaAxis_tracksAndJets_processed_2019-08-12.root
  // data/dijet_ppMC_RecoGen_Pythia6_pfCsJets_xjBins_wtaAxis_onlySeagull_processed_2019-07-13.root
  inputFile[kPbPb][kMC] = TFile::Open("data/PbPbMC_RecoReco_akFlowPuCs4PfJets_onlyJets_JECv5b_processed_2019-08-29.root");
  // data/PbPbMC_RecoReco_akFlowPuCs4PfJets_onlyJets_JECv5b_processed_2019-08-28.root
  // data/PbPbMC_RecoReco_akFlowPuCsPfJets_noUncorr_improvisedMixing_JECv4_noCorrections_processed_2019-08-09.root
  // data/PbPbMC_RecoReco_pfCsJets_xjBins_noUncOrInc_improvisedMixing_onlySeagull_wtaAxis_processed_2019-07-12.root
  // data/PbPbMC_RecoGen_pfCsJets_noUncorr_matchedCaloJets_subeNon0_improvisedMixing_onlySeagull_processed_2019-07-03.root
  // data/PbPbMC_RecoReco_skims_pfJets_noMixing_processed_2019-01-04.root
  
  // Load the necessary histograms to histogram managers
  for(int iSystem = 0; iSystem < knCollisionSystems; iSystem++){
    for(int iDataType = 0; iDataType < knDataTypes; iDataType++){
      inputManager[iSystem][iDataType] = new DijetHistogramManager(inputFile[iSystem][iDataType]);
      inputManager[iSystem][iDataType]->SetLoadEventInformation(true);
      inputManager[iSystem][iDataType]->SetLoadAllJets(true,true,true,true);
      inputManager[iSystem][iDataType]->LoadProcessedHistograms();
    }
  }
  
  // Open files for the closure tests
  DijetHistogramManager *closureManager[knCollisionSystems][knMonteCarloTypes];
  TFile *closureFile[knCollisionSystems][knMonteCarloTypes];
  closureFile[kPp][kRecoReco] = TFile::Open("data/dijet_ppMC_RecoReco_Pythia8_pfJets_wtaAxis_tracksAndJets_processed_2019-08-12.root");
  // data/dijet_ppMC_RecoReco_Pythia8_pfJets_wtaAxis_tracksAndJets_processed_2019-08-12.root
  // data/dijet_ppMC_RecoReco_Pythia6_pfCsJets_noUncOrInc_wtaAxis_onlySeagull_processed_2019-07-13.root
  closureFile[kPp][kRecoGen] = TFile::Open("data/dijet_ppMC_RecoGen_Pythia8_pfJets_wtaAxis_tracksAndJets_processed_2019-08-12.root");
  // data/dijet_ppMC_RecoGen_Pythia8_pfJets_wtaAxis_tracksAndJets_processed_2019-08-12.root
  // data/dijet_ppMC_RecoGen_Pythia6_pfCsJets_xjBins_wtaAxis_onlySeagull_processed_2019-07-13.root
  closureFile[kPp][kGenReco] = TFile::Open("data/dijet_ppMC_GenReco_Pythia6_pfCsJets_noUncOrInc_wtaAxis_onlySeagull_processed_2019-07-13.root");
  closureFile[kPp][kGenGen] = TFile::Open("data/dijet_ppMC_GenGen_Pythia6_pfCsJets_xjBins_wtaAxis_onlySeagull_processed_2019-07-13.root");
  closureFile[kPbPb][kRecoReco] = TFile::Open("data/PbPbMC_RecoReco_akFlowPuCsPfJets_noUncorr_improvisedMixing_JECv4_preprocessed_2019-08-09.root");
  // data/PbPbMC_RecoReco_akFlowPuCsPfJets_noUncorr_improvisedMixing_JECv4_preprocessed_2019-08-09.root
  // data/PbPbMC_RecoReco_pfCsJets_xjBins_noUncOrInc_improvisedMixing_onlySeagull_wtaAxis_processed_2019-07-12.root
  // data/PbPbMC_RecoReco_pfCsJets_noUncorr_5eveStrictMix_xjBins_seagullCheck_processed_2019-06-16.root
  // data/PbPbMC_RecoReco_skims_pfJets_noMixing_noJetLimit_noCorrelations_processed_2019-01-10.root
  // data/PbPbMC_RecoReco_skims_pfJets_noMixing_processed_2019-01-04.root
  closureFile[kPbPb][kRecoGen] = TFile::Open("data/PbPbMC_RecoGen_akFlowPuCsPfJets_noUncorr_improvisedMixing_noCorrections_JECv4_processed_2019-08-09.root");
  // data/PbPbMC_RecoGen_akFlowPuCsPfJets_noUncorr_improvisedMixing_noCorrections_JECv4_processed_2019-08-09.root
  // data/PbPbMC_RecoGen_pfCsJets_noUncOrInc_xjBins_improvisedMixing_onlySeagull_wtaAxis_processed_2019-07-12.root
  // data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_xjBins_2019-06-12_onlySeagull_processed.root
  // data/PbPbMC_RecoGen_skims_pfJets_noMixing_noJetLimit_noCorrelations_processed_2019-01-10.root
  // data/PbPbMC_RecoGen_skims_pfJets_noMixing_noJetLimit_noCorrelations_processed_2019-01-10.root
  closureFile[kPbPb][kGenReco] = TFile::Open("data/PbPbMC_GenReco_pfCsJets_noUncOrInc_xjBins_improvisedMixing_onlySeagull_wtaAxis_processed_2019-07-12.root");
  closureFile[kPbPb][kGenGen] = TFile::Open("data/PbPbMC_GenGen_pfCsJets_noUncOrInc_xjBins_improvisedMixing_onlySeagull_wtaAxis_processed_2019-07-12.root");
  
  // Load the necessary histograms to histogram managers
  for(int iSystem = 0; iSystem < knCollisionSystems; iSystem++){
    for(int iDataType = 0; iDataType < knMonteCarloTypes; iDataType++){
      closureManager[iSystem][iDataType] = new DijetHistogramManager(closureFile[iSystem][iDataType]);
      closureManager[iSystem][iDataType]->SetLoadAllTracks(true,true);
      closureManager[iSystem][iDataType]->SetLoadAllInclusiveTracks(true,true);
      closureManager[iSystem][iDataType]->LoadProcessedHistograms();
    }
  }
  
  // Define needed histograms
  TH1D *hVz[knCollisionSystems][knDataTypes+1];
  TH1D *hCentrality[knDataTypes+1];
  TH1D *trackPt[2][2][knCollisionSystems][knMonteCarloTypes][nCentralityBins]; // First bin = Dijet/inclusive
  TH1D *trackEta[2][2][knCollisionSystems][knMonteCarloTypes][nCentralityBins][nTrackPtBins+1]; // Second bin = Corrected/Uncorrected
  TH1D *trackPhi[2][2][knCollisionSystems][knMonteCarloTypes][nCentralityBins][nTrackPtBins+1]; // First bin = Dijet/inclusive
  TH1D *trackPtRatio[2][2][knCollisionSystems][nCentralityBins]; // Second bin = Corrected/Uncorrected
  TH1D *trackEtaRatio[2][2][knCollisionSystems][nCentralityBins][nTrackPtBins+1]; // First bin = Dijet/inclusive
  TH1D *trackPhiRatio[2][2][knCollisionSystems][nCentralityBins][nTrackPtBins+1]; // Second bin = Corrected/Uncorrected
  TH1D *jetPtNormalizer[2][knCollisionSystems][knMonteCarloTypes][nCentralityBins]; // First bin = Leading/inclusive
  TH1D *jetPt[3][knCollisionSystems][knDataTypes][nCentralityBins];  // First bin = Leading/inclusive/subleading
  TH1D *jetEta[3][knCollisionSystems][knDataTypes][nCentralityBins];  // First bin = Leading/inclusive/subleading
  TH1D *jetPhi[3][knCollisionSystems][knDataTypes][nCentralityBins];  // First bin = Leading/inclusive/subleading
  TH1D *jetPtDijet[2][nAsymmetryBins][knCollisionSystems][nCentralityBins];  // First bin = Leading/subleading
  TH1D *jetEtaDijet[2][nAsymmetryBins][knCollisionSystems][nCentralityBins]; // First bin = Leading/subleading
  TH1D *jetPhiDijet[2][nAsymmetryBins][knCollisionSystems][nCentralityBins]; // First bin = Leading/subleading
  
  // String for finding inclusive histograms
  const char *inclusiveString[2] = {"","Inclusive"};   // 0 = Tracks in dijet events, 1 = Inclusive tracks
  const char *inclusiveJetString[3] = {"leading","any","subleading"}; // 0 = Jets in dijet events, 1 = Inclusive jets, 2 = Subleading jets
  const char *leadingJetString[2] = {"leading","subleading"}; // 0 = Leading jets, 1 = Subleading jets
  const char *correctionString[2] = {"","Uncorrected"}; // 0 = Track correction included, 1 = No tracking corrections
  const char *inclusiveJetHeader[3] = {"Leading","Inclusive","Subleading"};
  const char *cutFinder[3] = {"","A0",""}; // A0 selects inclusive jets with pT cut of 120.
  double normalizationFactor;
  int jetTypeIndex[3] = {DijetHistogramManager::kLeadingJet, DijetHistogramManager::kAnyJet, DijetHistogramManager::kSubleadingJet};
  int asymmetryBinIndex[3] = {DijetHistogramManager::kMaxAsymmetryBins, 0, DijetHistogramManager::kMaxAsymmetryBins}; // Selecting asymmetry bin 0 for inclusive jet gives only jets above 120 GeV.
  
  // ******************************************
  // **    Read the histograms from files    **
  // ******************************************
  
  for(int iSystem = 0; iSystem < knCollisionSystems; iSystem++){
    for(int iDataType = 0; iDataType < knDataTypes; iDataType++){
      
      // Vz plots for both pp and PbPb
      hVz[iSystem][iDataType] = inputManager[iSystem][iDataType]->GetHistogramVertexZ();
      setStyleAndNormalize(hVz[iSystem][iDataType],iDataType);
      if(iDataType == kMC) {
        hVz[iSystem][iDataType+1] = inputManager[iSystem][iDataType]->GetHistogramVertexZWeighted();
        setStyleAndNormalize(hVz[iSystem][iDataType+1],iDataType+1);
      }
      
      // Centrality plots only for PbPb
      if(iSystem == kPbPb){
        hCentrality[iDataType] = inputManager[iSystem][iDataType]->GetHistogramCentrality();
        setStyleAndNormalize(hCentrality[iDataType],iDataType);
        if(iDataType == kMC){
          hCentrality[iDataType+1] = inputManager[iSystem][iDataType]->GetHistogramCentralityWeighted();
          setStyleAndNormalize(hCentrality[iDataType+1],iDataType+1);
        }
      }
    } // Data type loop
    
    // Read the jet kinematics histograms for kinematics figures
    for(int iJetType = 0; iJetType < 3; iJetType++){
      for(int iCentrality = firstCentralityBin; iCentrality <= lastCentralityBin; iCentrality++){
        
        // No centrality binning for pp
        if(iSystem == kPp && iCentrality > 0) continue;
        
        // We compare the jet kinematics data with Monte Carlo
        for(int iDataType = 0; iDataType < knDataTypes; iDataType++){

          // Read jet pT histograms
          jetPt[iJetType][iSystem][iDataType][iCentrality] = inputManager[iSystem][iDataType]->GetHistogramJetPt(jetTypeIndex[iJetType], iCentrality, asymmetryBinIndex[iJetType]);
          if(normalizeToHighPt){
            normalizationFactor = 1.0/jetPt[iJetType][iSystem][iDataType][iCentrality]->Integral(jetPt[iJetType][iSystem][iDataType][iCentrality]->FindBin(200), jetPt[iJetType][iSystem][iDataType][iCentrality]->GetNbinsX(), "width");
          } else {
            normalizationFactor = 1.0/jetPt[iJetType][iSystem][iDataType][iCentrality]->Integral("width");
          }
          jetPt[iJetType][iSystem][iDataType][iCentrality]->Scale(normalizationFactor);
          
          // Read jet eta histograms and scale them with the number of jets
          jetEta[iJetType][iSystem][iDataType][iCentrality] = inputManager[iSystem][iDataType]->GetHistogramJetEta(jetTypeIndex[iJetType], iCentrality, asymmetryBinIndex[iJetType]);
          jetEta[iJetType][iSystem][iDataType][iCentrality]->Scale(normalizationFactor);
          
          // Read jet phi histograms and scale them with the number of jets
          jetPhi[iJetType][iSystem][iDataType][iCentrality] = inputManager[iSystem][iDataType]->GetHistogramJetPhi(jetTypeIndex[iJetType], iCentrality, asymmetryBinIndex[iJetType]);
          jetPhi[iJetType][iSystem][iDataType][iCentrality]->Scale(normalizationFactor);
          
        } // Data type loop (real data / simulation)
         
        if(iJetType == 2) continue;  // Asymmetry binning only for leading and subleading jets (not for inclusive)
        
        // For the asymmetry binning, we only make the figures for real data
        for(int iAsymmetry = 0; iAsymmetry < inputManager[iSystem][kData]->GetNAsymmetryBins(); iAsymmetry++){
          
          // Read jet pT histograms in dijet asymmetry bins
          jetPtDijet[iJetType][iAsymmetry][iSystem][iCentrality] = inputManager[iSystem][kData]->GetHistogramJetPt(iJetType, iCentrality, iAsymmetry);
          if(normalizeToHighPt){
            normalizationFactor = 1.0/jetPtDijet[iJetType][iAsymmetry][iSystem][iCentrality]->Integral(jetPtDijet[iJetType][iAsymmetry][iSystem][iCentrality]->FindBin(200), jetPtDijet[iJetType][iAsymmetry][iSystem][iCentrality]->GetNbinsX(), "width");
          } else {
            normalizationFactor = 1.0/jetPtDijet[iJetType][iAsymmetry][iSystem][iCentrality]->Integral("width");
          }
          jetPtDijet[iJetType][iAsymmetry][iSystem][iCentrality]->Scale(normalizationFactor);
          
          // Read jet eta histograms in dijet asymmetry bins and scale them with the number of jets
          jetEtaDijet[iJetType][iAsymmetry][iSystem][iCentrality] = inputManager[iSystem][kData]->GetHistogramJetEta(iJetType, iCentrality, iAsymmetry);
          jetEtaDijet[iJetType][iAsymmetry][iSystem][iCentrality]->Scale(normalizationFactor);
          
          // Read jet phi histograms in dijet asymmetry bins and scale them with the number of jets
          jetPhiDijet[iJetType][iAsymmetry][iSystem][iCentrality] = inputManager[iSystem][kData]->GetHistogramJetPhi(iJetType, iCentrality, iAsymmetry);
          jetPhiDijet[iJetType][iAsymmetry][iSystem][iCentrality]->Scale(normalizationFactor);
          
        } // Asymmetry bin loop
      } // Centrality loop
    } // Jet type loop (leading/subleading)
    
    if(drawTrackClosure){
      // Read all the histograms needed for tracking closure plots
      for(int iInclusive = 0; iInclusive < 2; iInclusive++){
        
        for(int iMonteCarloType = 0; iMonteCarloType < knMonteCarloTypes; iMonteCarloType++){
          for(int iCentrality = firstCentralityBin; iCentrality <= lastCentralityBin; iCentrality++){
            
            // No centrality binning for pp
            if(iSystem == kPp && iCentrality > 0) continue;
            
            // Read jet pT histograms for normalizing the track histograms per jet
            jetPtNormalizer[iInclusive][iSystem][iMonteCarloType][iCentrality] =  closureManager[iSystem][iMonteCarloType]->GetHistogramJetPt(jetTypeIndex[iInclusive], iCentrality, DijetHistogramManager::kMaxAsymmetryBins);
            normalizationFactor = 1.0/jetPtNormalizer[iInclusive][iSystem][iMonteCarloType][iCentrality]->Integral("width");
            jetPtNormalizer[iInclusive][iSystem][iMonteCarloType][iCentrality]->Scale(normalizationFactor);
            
            for(int iCorrection = 0; iCorrection < 2; iCorrection++){
              
              // Read track pT histograms
              trackPt[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality] = closureManager[iSystem][iMonteCarloType]->GetHistogramTrackPt(iInclusive*2 + iCorrection, 0, iCentrality);
              trackPt[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality]->Scale(normalizationFactor);
              
              // Read track eta histograms without pT cut
              trackEta[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][nTrackPtBins] = closureManager[iSystem][iMonteCarloType]->GetHistogramTrackEta(iInclusive*2 + iCorrection, 0, iCentrality, closureManager[iSystem][iMonteCarloType]->GetNTrackPtBins());
              trackEta[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][nTrackPtBins]->Scale(normalizationFactor);
              
              // Read track phi histograms without pT cut
              trackPhi[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][nTrackPtBins] = closureManager[iSystem][iMonteCarloType]->GetHistogramTrackPhi(iInclusive*2 + iCorrection, 0, iCentrality, closureManager[iSystem][iMonteCarloType]->GetNTrackPtBins());
              trackPhi[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][nTrackPtBins]->Scale(normalizationFactor);
              
              for(int iTrackPt = firstTrackPtBin; iTrackPt <= lastTrackPtBin; iTrackPt++){
                
                // Read track eta histograms in pT bins
                trackEta[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][iTrackPt] = closureManager[iSystem][iMonteCarloType]->GetHistogramTrackEta(iInclusive*2 + iCorrection, 0, iCentrality, iTrackPt);
                trackEta[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][iTrackPt]->Scale(normalizationFactor);
                
                // Read track phi histograms in pT bins
                trackPhi[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][iTrackPt] = closureManager[iSystem][iMonteCarloType]->GetHistogramTrackPhi(iInclusive*2 + iCorrection, 0, iCentrality, iTrackPt);
                trackPhi[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][iTrackPt]->Scale(normalizationFactor);
                
              } // Track pT loop
            } // Corrected - Uncorrected loop
          } // Centrality loop
        } // Loop over Monte Carlo types
        
      } // Dijet - inclusive loop
    } // Collision system loop
  }
  
  // Drawing class for drawing your favorite root style
  JDrawer *drawer = new JDrawer();
  drawer->SetCanvasSize(700,700);
  TLegend *legend;
  
  // **************************************
  // **         Drawing vz plots         **
  // **************************************
  
  if(drawVz){
    
    double vzZoom[knCollisionSystems] = {0.042,0.052};
    
    for(int iSystem = 0; iSystem < knCollisionSystems; iSystem++){
      
      // Draw the first histogram and set up legend
      hVz[iSystem][0]->GetYaxis()->SetRangeUser(0,vzZoom[iSystem]);
      drawer->DrawHistogram(hVz[iSystem][0],"v_{z} (cm)","AU"," ");
      legend = new TLegend(0.155,0.78-0.08*iSystem,0.455,0.93-0.03*iSystem);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
      legend->AddEntry(hVz[iSystem][0],legendNames[iSystem][0],"p");
      
      // Draw the other histograms to the same figure
      for(int iDataType = 1; iDataType < knDataTypes+1; iDataType++){
        hVz[iSystem][iDataType]->Draw("same");
        legend->AddEntry(hVz[iSystem][iDataType],legendNames[iSystem][iDataType],"p");
      }
      
      // Draw also the legend
      legend->Draw();
      
      // Save the figures
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/vzComparison%s.pdf",systemString[iSystem].Data()));
      }
    }
  }

  // **************************************
  // **     Drawing centrality plots     **
  // **************************************
  
  if(drawCentrality){
    // Draw the first histogram and set up legend
    drawer->DrawHistogram(hCentrality[0],"centrality (%)","AU"," ");
    legend = new TLegend(0.6,0.7,0.9,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
    legend->AddEntry(hCentrality[0],legendNames[kPbPb][0],"p");
    
    // Draw the other histograms to the same figure
    for(int iDataType = 1; iDataType < knDataTypes+1; iDataType++){
      hCentrality[iDataType]->Draw("same");
      legend->AddEntry(hCentrality[iDataType],legendNames[kPbPb][iDataType],"p");
    }
    
    // Draw also the legend
    legend->Draw();
    
    // Save the figures
    if(saveFigures){
      gPad->GetCanvas()->SaveAs("figures/centralityComparison.pdf");
    }
  }
  
  // ***************************************
  // **    Drawing track closure plots    **
  // ***************************************
  
  if(drawTrackClosure){
    
    char header[100];
    char trackPtString[100];
    char trackPtName[100];
    char centralityString[100];
    char centralityName[100];
    char figureName[100];
    const char* zoomerName;
    drawer->SetDefaultAppearanceSplitCanvas();
    
    if(bigClosureZoom){
      zoomerName = "bigZoom_";
    } else {
      zoomerName = "";
    }
    
    for(int iInclusive = 0; iInclusive < 2; iInclusive++){
      for(int iSystem = 0; iSystem < knCollisionSystems; iSystem++){
        for(int iCentrality = firstCentralityBin; iCentrality <= lastCentralityBin; iCentrality++){
          
          if(iSystem == kPp && iCentrality > 0) continue;
          
          // Closure plots for track pT
          sprintf(header,"%s, %s", monteCarloName[iSystem],inclusiveLegend[iInclusive]);
          if(iSystem == kPp){
            sprintf(centralityString,"");
            sprintf(centralityName,"");
          } else {
            sprintf(centralityString,"Cent: %.0f-%.0f%%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
            sprintf(centralityName,"_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
          }
          
          sprintf(trackPtString,"p_{T} inclusive");
          
          /*
           * Parameters for plotter function
           *
           *  JDrawer *drawer = Drawer doing all the dirty work for drawing the histograms
           *  TH1 *genHistogram = Histogram with generator level tracks
           *  TH1 *recoHistogram = Histogram with corrected reconstructed tracks
           *  TH1 *uncorrectedHistogram = Histogram with uncorrected reconstructed tracks
           *  const char *xTitle = Title given for the x-axis of the histograms
           *  const char *yTitle = Title given for the y-axis of the histograms
           *  bool logScale = Use logarithmic scale for the main histogram
           *  int rebin = Rebin for the histograms
           *  double yZoom = Range for the y-axis in the main histogram
           *  bool bigZoom = Zoom to very close range to show exactly the deviation from one
           *  legendX = Left side x-position of the title
           *  legendY = Bottom side y-position of the title
           *  const char *header = Header for the legend
           *  const char *centralityString = String for the centrality
           *  const char *trackPtString = String for track pT
           */
          
          // Plot the closure for track pT
          plotTrackClosure(drawer, trackPt[iInclusive][0][iSystem][kRecoGen][iCentrality], trackPt[iInclusive][0][iSystem][kRecoReco][iCentrality], trackPt[iInclusive][1][iSystem][kRecoReco][iCentrality], "p_{T} (GeV)", "#frac{dN}{dp_{T}} (1/GeV)", true, ptRebin, -1, bigClosureZoom, 0.5, 0.5, header, centralityString, "");
          
          // Save the figures for track pT closure
          if(saveFigures){
            sprintf(figureName,"figures/trackPtClosure_%s%s%s%s.pdf", zoomerName, inclusiveJetHeader[iInclusive], systemString[iSystem].Data(), centralityName);
            gPad->GetCanvas()->SaveAs(figureName);
          }
          
          // Plot the closure for pT inclusive track eta
          plotTrackClosure(drawer, trackEta[iInclusive][0][iSystem][kRecoGen][iCentrality][nTrackPtBins], trackEta[iInclusive][0][iSystem][kRecoReco][iCentrality][nTrackPtBins], trackEta[iInclusive][1][iSystem][kRecoReco][iCentrality][nTrackPtBins], "#eta", "#frac{dN}{d#eta}", false, trackAngleRebin, trackEta[iInclusive][0][iSystem][kRecoReco][iCentrality][nTrackPtBins]->GetMaximum()*2.6, bigClosureZoom, 0.5, 0.5, header, centralityString, trackPtString);
          
          // Save the figures for track eta closure
          if(saveFigures){
            sprintf(figureName,"figures/trackEtaClosure_%s%s%s%s.pdf",zoomerName,inclusiveJetHeader[iInclusive],systemString[iSystem].Data(),centralityName);
            gPad->GetCanvas()->SaveAs(figureName);
          }
          
          // Plot the closure for pT inclusive track phi
          plotTrackClosure(drawer, trackPhi[iInclusive][0][iSystem][kRecoGen][iCentrality][nTrackPtBins], trackPhi[iInclusive][0][iSystem][kRecoReco][iCentrality][nTrackPtBins], trackPhi[iInclusive][1][iSystem][kRecoReco][iCentrality][nTrackPtBins], "#phi", "#frac{dN}{d#phi}", false, trackAngleRebin, trackPhi[iInclusive][0][iSystem][kRecoReco][iCentrality][nTrackPtBins]->GetMaximum()*2.6, bigClosureZoom, 0.5, 0.5, header, centralityString, trackPtString);
          
          // Save the figures for track phi closure
          if(saveFigures){
            sprintf(figureName,"figures/trackPhiClosure_%s%s%s%s.pdf",zoomerName,inclusiveJetHeader[iInclusive],systemString[iSystem].Data(),centralityName);
            gPad->GetCanvas()->SaveAs(figureName);
          }
          
          for(int iTrackPt = firstTrackPtBin; iTrackPt <= lastTrackPtBin; iTrackPt++){
            
            sprintf(trackPtString,"%.1f < p_{T} < %.1f",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
            sprintf(trackPtName,"_pT=%.0f-%.0f",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
            
            // Closure plots for track eta
            plotTrackClosure(drawer, trackEta[iInclusive][0][iSystem][kRecoGen][iCentrality][iTrackPt], trackEta[iInclusive][0][iSystem][kRecoReco][iCentrality][iTrackPt], trackEta[iInclusive][1][iSystem][kRecoReco][iCentrality][iTrackPt], "#eta", "#frac{dN}{d#eta}", false, trackAngleRebin, trackEta[iInclusive][0][iSystem][kRecoReco][iCentrality][iTrackPt]->GetMaximum()*2.6, bigClosureZoom, 0.5, 0.5, header, centralityString, trackPtString);
            
            // Save the figures for track eta closure
            if(saveFigures){
              sprintf(figureName,"figures/trackEtaClosure_%s%s%s%s%s.pdf",zoomerName,inclusiveJetHeader[iInclusive],systemString[iSystem].Data(),centralityName,trackPtName);
              gPad->GetCanvas()->SaveAs(figureName);
            }
            
            // Closure plots for track phi
            plotTrackClosure(drawer, trackPhi[iInclusive][0][iSystem][kRecoGen][iCentrality][iTrackPt], trackPhi[iInclusive][0][iSystem][kRecoReco][iCentrality][iTrackPt], trackPhi[iInclusive][1][iSystem][kRecoReco][iCentrality][iTrackPt], "#phi", "#frac{dN}{d#phi}", false, trackAngleRebin,  trackPhi[iInclusive][0][iSystem][kRecoReco][iCentrality][iTrackPt]->GetMaximum()*2.6, bigClosureZoom, 0.5, 0.5, header, centralityString, trackPtString);
            
            // Save the figures for track phi closure
            if(saveFigures){
              sprintf(figureName,"figures/trackPhiClosure_%s%s%s%s%s.pdf",zoomerName,inclusiveJetHeader[iInclusive],systemString[iSystem].Data(),centralityName,trackPtName);
              gPad->GetCanvas()->SaveAs(figureName);
            }
          } // Track pT loop
          
        } // Centrality loop
      } // Collision system loop
    } // Dijet-inclusive loop
  } // Plotting track closure
  
  // ****************************************
  // **    Drawing jet kinematics plots    **
  // ****************************************
  
  if(drawJetKinematics){
    
    char centralityString[100];
    char centralityName[100];
    char asymmetryString[100];
    char compactAsymmetryString[100];
    char figureName[200];
    int markerStyleData = 20;
    int markerStyleMC = 34;
    int markerColorData[3] = {kRed,kBlue,kBlack};
    int markerColorMC[3] = {kOrange,kCyan,kGreen};
    
    // Set nice drawing options for the drawer class
    drawer->Reset();
    drawer->SetCanvasSize(700,700);
    drawer->SetLeftMargin(0.22);
    drawer->SetTitleOffsetY(1.55);
    drawer->SetLabelOffsetY(0.01);
    
    // System ans centrality loops are common for data-simulation comparison and kinematics in asymmetry bins
    for(int iSystem = 0; iSystem < knCollisionSystems; iSystem++){
      for(int iCentrality = firstCentralityBin; iCentrality <= lastCentralityBin; iCentrality++){
        
        // No centrality binning for pp
        if(iSystem == kPp && iCentrality > 0) continue;
        
        if(iSystem == kPp){
          sprintf(centralityString,"");
          sprintf(centralityName,"");
        } else {
          sprintf(centralityString,"Cent: %.0f-%.0f%%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
          sprintf(centralityName,"_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
        }
        
        // Draw the comparison plots between data and simulation
        if(drawJetKinematicsMcComparison){
          
          /*
           * Parameters for jet kinematics plotter
           *
           *  JDrawer *drawer = Drawer doing all the dirty work for drawing the histograms
           *  TH1 *dataHistogram = Histogram with jet kinematics from real data
           *  TH1 *mcHistogram = Histogram with jet kinematics from simulation
           *  const char *xTitle = Title given for the x-axis of the histograms
           *  const char *yTitle = Title given for the y-axis of the histograms
           *  bool logScale = Use logarithmic scale for the main histogram
           *  int rebin = Rebin for the histograms
           *  double xRangeMin = Start for the x-range in the main histogram
           *  double yZoom = Range for the y-axis in the main histogram
           *  legendX = Left side x-position of the title
           *  legendY = Bottom side y-position of the title
           *  const char *header = System to be put to legend header
           *  const char *centralityString = String for the centrality
           *  const char *dataString = String describing data in the legend
           *  const char *mcString = String decribing the simulation in the legend
           *  int markerStyleData = Marker style to be used with the data histogram
           *  int markerColorData = Marker color to be used with the data histogram
           *  int markerStyleMC = Marker style to be used with the simulation histogram
           *  int markerColorMC = Marker color to be used with the simulation histogram
           */
          
          // Jet type loop here, because it is not needed for asymmetry plots
          for(int iJetType = 0; iJetType < 3; iJetType++){
            
            // Plot the kinematics figures for jet pT
            plotJetKinematics(drawer,jetPt[iJetType][iSystem][kData][iCentrality],jetPt[iJetType][iSystem][kMC][iCentrality],"p_{T} (GeV)","#frac{1}{N_{dijet}} #frac{dN}{dp_{T}} (1/GeV)",true,1,50,0.08,0.25,0.25,inclusiveJetHeader[iJetType],centralityString,legendNames[iSystem][0],monteCarloName[iSystem],markerStyleData,markerColorData[iJetType],markerStyleMC,markerColorMC[iJetType]);
            
            // Save the figures for jet pT
            if(saveFigures){
              sprintf(figureName,"figures/jetPtKinematics_%s%s%s.pdf",inclusiveJetHeader[iJetType],systemString[iSystem].Data(),centralityName);
              gPad->GetCanvas()->SaveAs(figureName);
            }
            
            // Plot the kinematics figures for jet eta
            plotJetKinematics(drawer,jetEta[iJetType][iSystem][kData][iCentrality],jetEta[iJetType][iSystem][kMC][iCentrality],"#eta","#frac{1}{N_{dijet}} #frac{dN}{d#eta}",false,1,-1000,0.45,0.4,0.25,inclusiveJetHeader[iJetType],centralityString,legendNames[iSystem][0],monteCarloName[iSystem],markerStyleData,markerColorData[iJetType],markerStyleMC,markerColorMC[iJetType]);
            
            // Save the figures for jet eta
            if(saveFigures){
              sprintf(figureName,"figures/jetEtaKinematics_%s%s%s.pdf",inclusiveJetHeader[iJetType],systemString[iSystem].Data(),centralityName);
              gPad->GetCanvas()->SaveAs(figureName);
            }
            
            // Plot the kinematics figures for jet phi
            plotJetKinematics(drawer,jetPhi[iJetType][iSystem][kData][iCentrality],jetPhi[iJetType][iSystem][kMC][iCentrality],"#phi","#frac{1}{N_{dijet}} #frac{dN}{d#phi}",false,1,-1000,0.22,0.4,0.25,inclusiveJetHeader[iJetType],centralityString,legendNames[iSystem][0],monteCarloName[iSystem],markerStyleData,markerColorData[iJetType],markerStyleMC,markerColorMC[iJetType]);
            
            // Save the figures for jet phi
            if(saveFigures){
              sprintf(figureName,"figures/jetPhiKinematics_%s%s%s.pdf",inclusiveJetHeader[iJetType],systemString[iSystem].Data(),centralityName);
              gPad->GetCanvas()->SaveAs(figureName);
            }
            
            // Compare PbPb pT spectrum to pp pT spectrum
            if(iSystem == kPbPb){
              plotJetKinematics(drawer, jetPt[iJetType][kPbPb][kData][iCentrality], jetPt[iJetType][kPp][kData][0], "p_{T} (GeV)", "#frac{1}{N_{dijet}} #frac{dN}{dp_{T}} (1/GeV)", true, 1, 50, 0.08, 0.25, 0.25, inclusiveJetHeader[iJetType], centralityString, legendNames[iSystem][0], "pp", markerStyleData, markerColorData[iJetType], markerStyleMC, markerColorMC[iJetType]);
              
              // Save pT spectrum comparison between pp and PbPb
              if(saveFigures){
                sprintf(figureName, "figures/jetPtSystemComparison_%s%s.pdf", inclusiveJetHeader[iJetType], centralityName);
                gPad->GetCanvas()->SaveAs(figureName);
              }
            }
          } // Jet type loop (leading/inclusive/subleading)
          
        } // Jet kinematics MC comparison if
        
        // Draw the jet kinematics plots in asymmetry bins
        if(drawJetKinematicsAsymmerty){
          
          /*
           * Parameters for the jet kinematics plotter in asymmetry bins
           *
           *  JDrawer *drawer = Drawer doing all the dirty work for drawing the histograms
           *  TH1 *allLeadingHistogram = Histogram with all leading jets
           *  TH1 *allSubleadingHistogram = Histogram with all subleading jets
           *  TH1 *leadingAsymmetryHistogram = Histograms with leading jets in a selected asymmetry bin
           *  TH1* subleadingAsymmetryHistogram = Histograms with subleading jets in a selected asymmetry bin
           *  const char *xTitle = Title given for the x-axis of the histograms
           *  const char *yTitle = Title given for the y-axis of the histograms
           *  bool logScale = Use logarithmic scale for the main histogram
           *  int rebin = Rebin for the histograms
           *  double xRangeMin = Start for the x-range in the main histogram
           *  double yZoom = Range for the y-axis in the main histogram
           *  bool drawLegend = Draw a legend to the figure
           *  legendX = Left side x-position of the title
           *  legendY = Bottom side y-position of the title
           *  const char *header = Header for the legend
           *  const char *centralityString = String for the centrality
           *  const char *asymmetryString = String for the asymmetry
           */
          
          // Asymmetry binning only for special asymmetry plots
          for(int iAsymmetry = 0; iAsymmetry < inputManager[iSystem][kData]->GetNAsymmetryBins() ; iAsymmetry++){
            
            // Define asymmetry strings
            sprintf(asymmetryString,"%.1f < %s < %.1f",inputManager[iSystem][kData]->GetCard()->GetLowBinBorderAsymmetry(iAsymmetry) ,inputManager[iSystem][kData]->GetCard()->GetAsymmetryBinType(), inputManager[iSystem][kData]->GetCard()->GetHighBinBorderAsymmetry(iAsymmetry));
            sprintf(compactAsymmetryString,"A%d",iAsymmetry);
            
            // Plot the kinematics figures for jet pT
            // Note: For jetPt 0 = leading and 2 = subleading and for jetPtDijet 0 = leading and 1 = subleading
            plotJetKinematicsAsymmetry(drawer, jetPt[0][iSystem][kData][iCentrality], jetPt[2][iSystem][kData][iCentrality], jetPtDijet[0][iAsymmetry][iSystem][iCentrality], jetPtDijet[1][iAsymmetry][iSystem][iCentrality], "p_{T} (GeV)", "#frac{1}{N_{dijet}} #frac{dN}{dp_{T}} (1/GeV)", true, 1, 50, 0.06, false, 0.2, 0.2, legendNames[iSystem][0], centralityString, asymmetryString);
            
            // Save the figures for jet pT
            if(saveFigures){
              sprintf(figureName,"figures/jetPtAsymmetryKinematics_%s%s_%s.pdf", legendNames[iSystem][0], centralityName, compactAsymmetryString);
              gPad->GetCanvas()->SaveAs(figureName);
            }
            
            // Plot the kinematics figures for jet eta
            // Note: For jetEta 0 = leading and 2 = subleading and for jetEtaDijet 0 = leading and 1 = subleading
            plotJetKinematicsAsymmetry(drawer, jetEta[0][iSystem][kData][iCentrality], jetEta[2][iSystem][kData][iCentrality], jetEtaDijet[0][iAsymmetry][iSystem][iCentrality], jetEtaDijet[1][iAsymmetry][iSystem][iCentrality], "#eta", "#frac{1}{N_{dijet}} #frac{dN}{d#eta}", false, 1, -1000, 1, true, 0.26, 0.6, legendNames[iSystem][0], centralityString, asymmetryString);
            
            // Save the figures for jet eta
            if(saveFigures){
              sprintf(figureName, "figures/jetEtaAsymmetryKinematics_%s%s_%s.pdf", legendNames[iSystem][0], centralityName, compactAsymmetryString);
              gPad->GetCanvas()->SaveAs(figureName);
            }
            
            // Plot the kinematics figures for jet phi
            // Note: For jetPhi 0 = leading and 2 = subleading and for jetPhiDijet 0 = leading and 1 = subleading
            plotJetKinematicsAsymmetry(drawer, jetPhi[0][iSystem][kData][iCentrality], jetPhi[2][iSystem][kData][iCentrality], jetPhiDijet[0][iAsymmetry][iSystem][iCentrality], jetPhiDijet[1][iAsymmetry][iSystem][iCentrality], "#phi","#frac{1}{N_{dijet}} #frac{dN}{d#phi}", false, 1, -1000, 0.3, false, 0.3, 0.2, legendNames[iSystem][0], centralityString, asymmetryString);
            
            // Save the figures for jet phi
            if(saveFigures){
              sprintf(figureName, "figures/jetPhiAsymmetryKinematics_%s%s_%s.pdf", legendNames[iSystem][0], centralityName, compactAsymmetryString);
              gPad->GetCanvas()->SaveAs(figureName);
            }
            
          } // Dijet asymmetry loop
          
        } // Jet kinematics asymmetry if
        
      } // Centrality loop
    } // System loop (pp/PbPb)
    
    
  } // Drawing jet kinematics if
}
