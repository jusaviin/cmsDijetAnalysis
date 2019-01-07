#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

//Define types for histogram arrays
enum enumCollisionSystem {kPp, kPbPb, knCollisionSystems};
enum enumDataType {kData, kMC, knDataTypes};
enum enumMonteCarloType{kRecoReco, kRecoGen, kGenReco, kGenGen, knMonteCarloTypes};

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
  int goodStyle = 20;
  
  // Normalize the histogram to 1
  double integral = histogram->Integral();
  histogram->Scale(1.0/integral);
  
  // Put the histogram in the bench of a stylist
  histogram->SetMarkerStyle(goodStyle);
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
 *  legendX = Left side x-position of the title
 *  legendY = Bottom side y-position of the title
 *  const char *header = Header for the legend
 *  const char *centralityString = String for the centrality
 *  const char *trackPtString = String for track pT
 */
void plotTrackClosure(JDrawer *drawer, TH1 *genHistogram, TH1 *recoHistogram, TH1 *uncorrectedHistogram, const char *xTitle, const char *yTitle, bool logScale, int rebin, double yZoom, double legendX, double legendY, const char *header, const char *centralityString, const char *trackPtString){

  // Helper variable to name the histograms
  char namer[200];
  
  drawer->SetLogY(logScale);
  drawer->CreateSplitCanvas();
  genHistogram->Rebin(rebin);
  if(yZoom > 0) genHistogram->GetYaxis()->SetRangeUser(0,yZoom);
  drawer->DrawHistogramToUpperPad(genHistogram,xTitle,yTitle," ");

  double ySpace = 0;
  if(strncmp(centralityString,"",2) != 0) ySpace += 0.025;
  if(strncmp(trackPtString,"",2) != 0) ySpace += 0.025;
  
  recoHistogram->SetLineColor(kBlue);
  recoHistogram->Rebin(rebin);
  recoHistogram->Draw("same");
  uncorrectedHistogram->SetLineColor(kRed);
  uncorrectedHistogram->Rebin(rebin);
  uncorrectedHistogram->Draw("same");

  TLegend *legend = new TLegend(legendX,legendY-ySpace,legendX+0.3,legendY+0.25+ySpace);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

  legend->SetHeader(header);
  if(strncmp(centralityString,"",2) != 0) legend->AddEntry((TObject*) 0,centralityString,"");
  if(strncmp(trackPtString,"",2) != 0) legend->AddEntry((TObject*) 0,trackPtString,"");
  legend->AddEntry(genHistogram,"Gen. particles","l");
  legend->AddEntry(recoHistogram,"Corr. tracks","l");
  legend->AddEntry(uncorrectedHistogram,"Reco. tracks","l");
  legend->Draw();

  sprintf(namer,"%sRatio",recoHistogram->GetName());
  TH1D *recoRatio = (TH1D*) recoHistogram->Clone(namer);
  recoRatio->Divide(genHistogram);

  sprintf(namer,"%sRatio",uncorrectedHistogram->GetName());
  TH1D *uncorrectedRatio = (TH1D*) uncorrectedHistogram->Clone(namer);
  uncorrectedRatio->Divide(genHistogram);

  drawer->SetLogY(false);
  recoRatio->GetYaxis()->SetRangeUser(0,2);
  drawer->DrawHistogramToLowerPad(recoRatio,xTitle,"Reco/Gen"," ");
  uncorrectedRatio->Draw("same");

}

/*
 * Plotter for closure histograms to the analysis note
 */
void closurePlotter(){
  
  ///////////////////
  // Configuration //
  ///////////////////
  
  bool saveFigures = false;          // Save the figures to a file
  
  bool drawCentrality = false;        // Draw the QA plots for spillover correction
  bool drawVz = false;                // Draw the QA plots for seagull correction
  bool drawTrackClosure = false;      // Draw the tracking closures
  bool drawJetKinematics = false;     // Draw the jet kinematics figures
  
  int ptRebin = 4;                  // Rebin for track pT closure histograms (there are 500 bins)
  
  // Read the number of bins from histogram manager
  DijetHistogramManager *dummyManager = new DijetHistogramManager();
  const int nCentralityBins = dummyManager->GetNCentralityBins();
  const int nTrackPtBins = dummyManager->GetNTrackPtBins();
  const int nAsymmetryBins = 4;
  double centralityBinBorders[] = {0,10,30,50,100};       // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,300};       // Bin borders for track pT
  double asymmetryBinBorders[] = {0,0.11,0.22,0.33,0.75}; // Bin borders for dijet asymmetry
  
  // Select which bins to plot
  int firstCentralityBin = 0;
  int lastCentralityBin = nCentralityBins-1;
  
  int firstTrackPtBin = 5;
  int lastTrackPtBin = 5;
  
  /////////////////
  // Config done //
  /////////////////
  
  const char *legendNames[knCollisionSystems][knDataTypes+1] = {{"pp","Raw Pythia", "Weighted Pythia"},{"PbPb","Raw P+H", "Weighted P+H"}};
  const char *monteCarloName[knCollisionSystems] = {"Pythia6","Pythia+Hydjet"};
  const char *inclusiveLegend[2] = {"dijet","inclusive"};
  TString systemString[knCollisionSystems] = {"Pp","PbPb"};
  
  // Open files containing the QA histograms
  TFile *inputFile[knCollisionSystems][knDataTypes];
  inputFile[kPp][kData] = TFile::Open("data/dijet_pp_highForest_pfJets_processed_2018-09-14.root");
  inputFile[kPbPb][kData] = TFile::Open("data/dijetPbPb_skims_pfJets_pfCandAxis_noInclusiveOrUncorrected_10mixedEvents_smoothedMixing_processed_2018-11-19.root");
  inputFile[kPp][kMC] = TFile::Open("data/dijet_ppMC_RecoReco_mergedSkims_Pythia6_pfJets_processed_2018-09-15.root");
  inputFile[kPbPb][kMC] = TFile::Open("data/PbPbMC_GenGen_skims_pfJets_pfCandAxis_noInclOrUncorr_10eveMixed_sube0_smoothedMixing_processed_2018-11-19.root");
  
  // Open files for the closure tests
  TFile *closureFile[knCollisionSystems][knMonteCarloTypes];
  closureFile[kPp][kRecoReco] = TFile::Open("data/dijet_ppMC_RecoReco_mergedSkims_Pythia6_pfJets_processed_2018-09-15.root");
  closureFile[kPp][kRecoGen] = TFile::Open("data/dijet_ppMC_RecoGen_mergedSkims_Pythia6_pfJets_processed_2018-09-15.root");
  closureFile[kPp][kGenReco] = TFile::Open("data/dijet_ppMC_GenReco_mergedSkims_Pythia6_processed_2018-08-16.root");
  closureFile[kPp][kGenGen] = TFile::Open("data/dijet_ppMC_GenGen_mergedSkims_Pythia6_pfJets_processed_2018-09-15.root");
  closureFile[kPbPb][kRecoReco] = TFile::Open("data/PbPbMC_RecoReco_skims_pfJets_pfCandAxis_noMixing_processed_2018-12-06.root");
  closureFile[kPbPb][kRecoGen] = TFile::Open("data/PbPbMC_RecoGen_skims_pfJets_pfCandAxis_noMixing_processed_2018-12-06.root");
  closureFile[kPbPb][kGenReco] = TFile::Open("data/PbPbMC_GenReco_skims_pfJets_pfCandAxis_noMixing_processed_2018-12-06.root");
  closureFile[kPbPb][kGenGen] = TFile::Open("data/PbPbMC_GenGen_skims_pfJets_pfCandAxis_noMixing_processed_2018-12-06.root");
  
  // Define needed histograms
  TH1D *hVz[knCollisionSystems][knDataTypes+1];
  TH1D *hCentrality[knDataTypes+1];
  TH1D *trackPt[2][2][knCollisionSystems][knMonteCarloTypes][nCentralityBins]; // First bin = Dijet/inclusive
  TH1D *trackEta[2][2][knCollisionSystems][knMonteCarloTypes][nCentralityBins][nTrackPtBins+1]; // Second bin = Corrected/Uncorrected
  TH1D *trackPhi[2][2][knCollisionSystems][knMonteCarloTypes][nCentralityBins][nTrackPtBins+1]; // First bin = Dijet/inclusive
  TH1D *trackPtRatio[2][2][knCollisionSystems][nCentralityBins]; // Second bin = Corrected/Uncorrected
  TH1D *trackEtaRatio[2][2][knCollisionSystems][nCentralityBins][nTrackPtBins+1]; // First bin = Dijet/inclusive
  TH1D *trackPhiRatio[2][2][knCollisionSystems][nCentralityBins][nTrackPtBins+1]; // Second bin = Corrected/Uncorrected
  TH1D *jetPt[3][knCollisionSystems][knMonteCarloTypes][nCentralityBins];  // First bin = Leading/inclusive/subleading
  TH1D *jetEta[3][knCollisionSystems][knMonteCarloTypes][nCentralityBins];  // First bin = Leading/inclusive/subleading
  TH1D *jetPhi[3][knCollisionSystems][knMonteCarloTypes][nCentralityBins];  // First bin = Leading/inclusive/subleading
  TH1D *jetPtDijet[2][nAsymmetryBins][knCollisionSystems][knMonteCarloTypes][nCentralityBins];  // First bin = Leading/subleading
  TH1D *jetEtaDijet[2][nAsymmetryBins][knCollisionSystems][knMonteCarloTypes][nCentralityBins]; // First bin = Leading/subleading
  TH1D *jetPhiDijet[2][nAsymmetryBins][knCollisionSystems][knMonteCarloTypes][nCentralityBins]; // First bin = Leading/subleading
  
  // String for finding inclusive histograms
  const char *inclusiveString[2] = {"","Inclusive"};   // 0 = Tracks in dijet events, 1 = Inclusive tracks
  const char *inclusiveJetString[3] = {"leading","any","subleading"}; // 0 = Jets in dijet events, 1 = Inclusive jets, 2 = Subleading jets
  const char *leadingJetString[2] = {"leading","subleading"}; // 0 = Leading jets, 1 = Subleading jets
  const char *correctionString[2] = {"","Uncorrected"}; // 0 = Track correction included, 1 = No tracking corrections
  char namer[200];
  double normalizationFactor;
  
  // Read the histograms from files
  char histogramNamer[200];
  for(int iSystem = 0; iSystem < knCollisionSystems; iSystem++){
    for(int iDataType = 0; iDataType < knDataTypes; iDataType++){
      
      // Vz plots for both pp and PbPb
      hVz[iSystem][iDataType] = (TH1D*) inputFile[iSystem][iDataType]->Get("vertexZ");
      setStyleAndNormalize(hVz[iSystem][iDataType],iDataType);
      if(iDataType == kMC) {
        hVz[iSystem][iDataType+1] = (TH1D*) inputFile[iSystem][iDataType]->Get("vertexZweighted");
        setStyleAndNormalize(hVz[iSystem][iDataType+1],iDataType+1);
      }
      
      // Centrality plots only for PbPb
      if(iSystem == kPbPb){
        hCentrality[iDataType] = (TH1D*) inputFile[iSystem][iDataType]->Get("centrality");
        setStyleAndNormalize(hCentrality[iDataType],iDataType);
        if(iDataType == kMC){
          hCentrality[iDataType+1] = (TH1D*) inputFile[iSystem][iDataType]->Get("centralityWeighted");
          setStyleAndNormalize(hCentrality[iDataType+1],iDataType+1);
        }
      }
    } // Data type loop
    
    for(int iInclusive = 0; iInclusive < 3; iInclusive++){
      
      for(int iMonteCarloType = 0; iMonteCarloType < knMonteCarloTypes; iMonteCarloType++){
        for(int iCentrality = firstCentralityBin; iCentrality <= lastCentralityBin; iCentrality++){
          
          // No centrality binning for pp
          if(iSystem == kPp && iCentrality > 0) continue;
          
          // Dijet asymmetry binning is implemented only for leading and subleading jet histograms
          if(iInclusive < 2){
            for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
              
              // Read jet pT histograms in dijet asymmetry bins
              sprintf(namer,"%sJet/%sJetPt_C%dA%d",leadingJetString[iInclusive],leadingJetString[iInclusive],iCentrality,iAsymmetry);
              jetPtDijet[iInclusive][iAsymmetry][iSystem][iMonteCarloType][iCentrality] = (TH1D*) closureFile[iSystem][iMonteCarloType]->Get(namer);
              normalizationFactor = 1.0/jetPtDijet[iInclusive][iAsymmetry][iSystem][iMonteCarloType][iCentrality]->Integral("width");
              jetPtDijet[iInclusive][iAsymmetry][iSystem][iMonteCarloType][iCentrality]->Scale(normalizationFactor);
              
              // Read jet eta histograms in dijet asymmetry bins and scale them with the number of jets
              sprintf(namer,"%sJet/%sJetEta_C%dA%d",leadingJetString[iInclusive],leadingJetString[iInclusive],iCentrality,iAsymmetry);
              jetEtaDijet[iInclusive][iAsymmetry][iSystem][iMonteCarloType][iCentrality] = (TH1D*) closureFile[iSystem][iMonteCarloType]->Get(namer);
              jetEtaDijet[iInclusive][iAsymmetry][iSystem][iMonteCarloType][iCentrality]->Scale(normalizationFactor);
              
              // Read jet phi histograms in dijet asymmetry bins and scale them with the number of jets
              sprintf(namer,"%sJet/%sJetPhi_C%dA%d",leadingJetString[iInclusive],leadingJetString[iInclusive],iCentrality,iAsymmetry);
              jetPhiDijet[iInclusive][iAsymmetry][iSystem][iMonteCarloType][iCentrality] = (TH1D*) closureFile[iSystem][iMonteCarloType]->Get(namer);
              jetPhiDijet[iInclusive][iAsymmetry][iSystem][iMonteCarloType][iCentrality]->Scale(normalizationFactor);
              
            }
          }
          
          // Read jet pT histograms
          sprintf(namer,"%sJet/%sJetPt_C%d",inclusiveJetString[iInclusive],inclusiveJetString[iInclusive],iCentrality);
          jetPt[iInclusive][iSystem][iMonteCarloType][iCentrality] = (TH1D*) closureFile[iSystem][iMonteCarloType]->Get(namer);
          normalizationFactor = 1.0/jetPt[iInclusive][iSystem][iMonteCarloType][iCentrality]->Integral("width");
          jetPt[iInclusive][iSystem][iMonteCarloType][iCentrality]->Scale(normalizationFactor);
          
          // Read jet eta histograms and scale them with the number of jets
          sprintf(namer,"%sJet/%sJetEta_C%d",inclusiveJetString[iInclusive],inclusiveJetString[iInclusive],iCentrality);
          jetEta[iInclusive][iSystem][iMonteCarloType][iCentrality] = (TH1D*) closureFile[iSystem][iMonteCarloType]->Get(namer);
          jetEta[iInclusive][iSystem][iMonteCarloType][iCentrality]->Scale(normalizationFactor);
          
          // Read jet phi histograms and scale them with the number of jets
          sprintf(namer,"%sJet/%sJetPhi_C%d",inclusiveJetString[iInclusive],inclusiveJetString[iInclusive],iCentrality);
          jetPhi[iInclusive][iSystem][iMonteCarloType][iCentrality] = (TH1D*) closureFile[iSystem][iMonteCarloType]->Get(namer);
          jetPhi[iInclusive][iSystem][iMonteCarloType][iCentrality]->Scale(normalizationFactor);
          
          if(iInclusive == 2) continue; // Track histograms only for dijet events and all events
          
          for(int iCorrection = 0; iCorrection < 2; iCorrection++){
            
            // Read track pT histograms
            sprintf(namer,"track%s%s/track%s%sPt_SameEvent_C%d",inclusiveString[iInclusive],correctionString[iCorrection],inclusiveString[iInclusive],correctionString[iCorrection],iCentrality);
            trackPt[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality] = (TH1D*) closureFile[iSystem][iMonteCarloType]->Get(namer);
            trackPt[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality]->Scale(normalizationFactor);
            
            // Read track eta histograms without pT cut
            sprintf(namer,"track%s%s/track%s%sEta_SameEvent_C%dT6",inclusiveString[iInclusive],correctionString[iCorrection],inclusiveString[iInclusive],correctionString[iCorrection],iCentrality);
            trackEta[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][nTrackPtBins] = (TH1D*) closureFile[iSystem][iMonteCarloType]->Get(namer);
            trackEta[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][nTrackPtBins]->Scale(normalizationFactor);
            
            // Read track phi histograms without pT cut
            sprintf(namer,"track%s%s/track%s%sPhi_SameEvent_C%dT6",inclusiveString[iInclusive],correctionString[iCorrection],inclusiveString[iInclusive],correctionString[iCorrection],iCentrality);
            trackPhi[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][nTrackPtBins] = (TH1D*) closureFile[iSystem][iMonteCarloType]->Get(namer);
            trackPhi[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][nTrackPtBins]->Scale(normalizationFactor);
            
            for(int iTrackPt = firstTrackPtBin; iTrackPt <= lastTrackPtBin; iTrackPt++){
              
              // Read track eta histograms in pT bins
              sprintf(namer,"track%s%s/track%s%sEta_SameEvent_C%dT%d",inclusiveString[iInclusive],correctionString[iCorrection],inclusiveString[iInclusive],correctionString[iCorrection],iCentrality,iTrackPt);
              trackEta[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][iTrackPt] = (TH1D*) closureFile[iSystem][iMonteCarloType]->Get(namer);
              trackEta[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][iTrackPt]->Scale(normalizationFactor);
              
              // Read track phi histograms in pT bins
              sprintf(namer,"track%s%s/track%s%sPhi_SameEvent_C%dT%d",inclusiveString[iInclusive],correctionString[iCorrection],inclusiveString[iInclusive],correctionString[iCorrection],iCentrality,iTrackPt);
              trackPhi[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][iTrackPt] = (TH1D*) closureFile[iSystem][iMonteCarloType]->Get(namer);
              trackPhi[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][iTrackPt]->Scale(normalizationFactor);
              
            } // Track pT loop
          } // Corrected - Uncorrected loop
        } // Centrality loop
      } // Loop over Monte Carlo types
      
    } // Dijet - inclusive loop
  } // Collision system loop
  
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
    char centralityString[100];
    drawer->SetDefaultAppearanceSplitCanvas();
    
    // Bins: Dijet/inclusive, pp/PbPb, centrality, track pT
    double phiZoom[2][2][4][7] =
    // Track pT = 0.7-1  1-2   2-3   3-4   4-8   8-300   Inclusive
    // Dijet Pythia
              {{{{  5,    6,   2.2,  1.4,  2.2,   1.6,        20   },
                {0,0,0,0,0,0,0},    {0,0,0,0,0,0,0},    {0,0,0,0,0,0,0} },
    // Inclusive Pythia
                {{ 2.2,   3,    1,   0.5,   1,    0.6,        10   },
                {0,0,0,0,0,0,0},    {0,0,0,0,0,0,0},    {0,0,0,0,0,0,0}}},
    // Dijet Pythia+Hydjet Centrality 0-10 %
                {{{350,  380,   60,   12,   6,     2,       800   },
    // Dijet Pythia+Hydjet Centrality 10-30 %
                  {240,  240,   40,   8,   4.5,    2,       500   },
    // Dijet Pythia+Hydjet Centrality 30-50 %
                  { 90,   90,   18,   4,   3.4,    2,       220   },
    // Dijet Pythia+Hydjet Centrality 50-100 %
                  { 40,   40,    7,   2,   2.6,    2,        78   }},
    // Inclusive Pythia+Hydjet Centrality 0-10 %
                 {{ 40,   44,    7,  1.4,  0.6,   0.15,       85   },
    // Inclusive Pythia+Hydjet Centrality 10-30 %
                  { 40,   44,    8,  1.6,  0.8,   0.25,        90   },
    // Inclusive Pythia+Hydjet Centrality 30-50 %
                  { 34,   34,    6,  1.5,   1,    0.45,        75   },
    // Inclusive Pythia+Hydjet Centrality 50-100 %
                  { 15,   18,    3,   1,    1,    0.55,        40   }}}};
    
    // Bins: Dijet/inclusive, pp/PbPb, centrality, track pT
    double etaZoom[2][2][4][7] =
    // Track pT = 0.7-1  1-2   2-3   3-4   4-8   8-300   Inclusive
    // Dijet Pythia
             {{{{   6,   10,    4,   2.5,  4.5,   3.2,        32   },
          {0,0,0,0,0,0,0},    {0,0,0,0,0,0,0},    {0,0,0,0,0,0,0} },
      // Inclusive Pythia
               {{   3,    5,   1.6,   1,   1.6,    1,        14   },
         {0,0,0,0,0,0,0},    {0,0,0,0,0,0,0},    {0,0,0,0,0,0,0}}},
      // Dijet Pythia+Hydjet Centrality 0-10 %
              {{{ 550,  550,   100,   20,  10,     4,       1200   },
        // Dijet Pythia+Hydjet Centrality 10-30 %
                { 300,  340,    70,   15,   9,     4,       750   },
        // Dijet Pythia+Hydjet Centrality 30-50 %
                { 130,  130,    30,    7,   6,     4,       330   },
        // Dijet Pythia+Hydjet Centrality 50-100 %
                {  50,   60,    12,    4,   5,     4,       120   }},
        // Inclusive Pythia+Hydjet Centrality 0-10 %
                {{ 60,   60,    12,   2.4,  1,    0.3,       140   },
          // Inclusive Pythia+Hydjet Centrality 10-30 %
                 { 60,   60,    12,   2.6, 1.4,   0.4,       135   },
          // Inclusive Pythia+Hydjet Centrality 30-50 %
                 { 50,   50,     8,   2.4, 1.8,   0.8,        120   },
          // Inclusive Pythia+Hydjet Centrality 50-100 %
                 { 25,   25,     5,   1.6, 1.6,    1,        60   }}}};
    
    for(int iInclusive = 0; iInclusive < 2; iInclusive++){
      for(int iSystem = 0; iSystem < knCollisionSystems; iSystem++){
        for(int iCentrality = firstCentralityBin; iCentrality <= lastCentralityBin; iCentrality++){
          
          if(iSystem == kPp && iCentrality > 0) continue;
          
          // Closure plots for track pT
          sprintf(header,"%s, %s", monteCarloName[iSystem],inclusiveLegend[iInclusive]);
          if(iSystem == kPp){
            sprintf(centralityString,"");
          } else {
            sprintf(centralityString,"Cent: %.0f-%.0f%%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
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
           *  legendX = Left side x-position of the title
           *  legendY = Bottom side y-position of the title
           *  const char *header = Header for the legend
           *  const char *centralityString = String for the centrality
           *  const char *trackPtString = String for track pT
           */
          
          // Plot the closure for track pT
          plotTrackClosure(drawer, trackPt[iInclusive][0][iSystem][kRecoGen][iCentrality], trackPt[iInclusive][0][iSystem][kRecoReco][iCentrality], trackPt[iInclusive][1][iSystem][kRecoReco][iCentrality], "p_{T} (GeV)", "#frac{dN}{dp_{T}} (1/GeV)", true, ptRebin, -1, 0.5, 0.5, header, centralityString, "");
          
          // Plot the closure for pT inclusive track eta
          plotTrackClosure(drawer, trackEta[iInclusive][0][iSystem][kRecoGen][iCentrality][nTrackPtBins], trackEta[iInclusive][0][iSystem][kRecoReco][iCentrality][nTrackPtBins], trackEta[iInclusive][1][iSystem][kRecoReco][iCentrality][nTrackPtBins], "#eta", "#frac{dN}{d#eta}", false, 1, etaZoom[iSystem][iInclusive][iCentrality][nTrackPtBins], 0.5, 0.5, header, centralityString, trackPtString);
          
          // Plot the closure for pT inclusive track phi
          plotTrackClosure(drawer, trackPhi[iInclusive][0][iSystem][kRecoGen][iCentrality][nTrackPtBins], trackPhi[iInclusive][0][iSystem][kRecoReco][iCentrality][nTrackPtBins], trackPhi[iInclusive][1][iSystem][kRecoReco][iCentrality][nTrackPtBins], "#phi", "#frac{dN}{d#phi}", false, 1, phiZoom[iSystem][iInclusive][iCentrality][nTrackPtBins], 0.5, 0.5, header, centralityString, trackPtString);
          
          for(int iTrackPt = firstTrackPtBin; iTrackPt <= lastTrackPtBin; iTrackPt++){
            
            sprintf(trackPtString,"%.1f < p_{T} < %.1f",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
            
            // Closure plots for track eta
            plotTrackClosure(drawer, trackEta[iInclusive][0][iSystem][kRecoGen][iCentrality][iTrackPt], trackEta[iInclusive][0][iSystem][kRecoReco][iCentrality][iTrackPt], trackEta[iInclusive][1][iSystem][kRecoReco][iCentrality][iTrackPt], "#eta", "#frac{dN}{d#eta}", false, 1, etaZoom[iSystem][iInclusive][iCentrality][iTrackPt], 0.5, 0.5, header, centralityString, trackPtString);
            
            // Closure plots for track phi
            plotTrackClosure(drawer, trackPhi[iInclusive][0][iSystem][kRecoGen][iCentrality][iTrackPt], trackPhi[iInclusive][0][iSystem][kRecoReco][iCentrality][iTrackPt], trackPhi[iInclusive][1][iSystem][kRecoReco][iCentrality][iTrackPt], "#phi", "#frac{dN}{d#phi}", false, 1, phiZoom[iSystem][iInclusive][iCentrality][iTrackPt], 0.5, 0.5, header, centralityString, trackPtString);
          } // Track pT loop
          
        } // Centrality loop
      } // Collision system loop
    } // Dijet- inclusive loop
  }
}
