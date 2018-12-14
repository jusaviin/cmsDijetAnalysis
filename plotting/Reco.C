#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

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
 * Plotter for closure histograms to the analysis note
 */
void closurePlotter(){
  
  ///////////////////
  // Configuration //
  ///////////////////
  
  bool saveFigures = false;          // Save the figures to a file
  
  bool drawCentrality = false;        // Draw the QA plots for spillover correction
  bool drawVz = false;                // Draw the QA plots for seagull correction
  bool drawTrackClosure = true;      // Draw the tracking closures
  
  /////////////////
  // Config done //
  /////////////////
  
  //Define types for histogram arrays
  enum enumCollisionSystem {kPp, kPbPb, knCollisionSystems};
  enum enumDataType {kData, kMC, knDataTypes};
  enum enumMonteCarloType{kRecoReco, kRecoGen, kGenReco, kGenGen, knMonteCarloTypes};
  
  const char *legendNames[knCollisionSystems][knDataTypes+1] = {{"pp","Raw Pythia", "Weighted Pythia"},{"PbPb","Raw P+H", "Weighted P+H"}};
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
  
  // Read the number of bins from histogram manager
  DijetHistogramManager *dummyManager = new DijetHistogramManager();
  const int nCentralityBins = dummyManager->GetNCentralityBins();
  const int nTrackPtBins = dummyManager->GetNTrackPtBins();
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,300};  // Bin borders for track pT
  
  // Define needed histograms
  TH1D *hVz[knCollisionSystems][knDataTypes+1];
  TH1D *hCentrality[knDataTypes+1];
  TH1D *trackPt[2][2][knCollisionSystems][knMonteCarloTypes][nCentralityBins]; // First bin = Dijet/inclusive
  TH1D *trackEta[2][2][knCollisionSystems][knMonteCarloTypes][nCentralityBins][nTrackPtBins+1]; // Second bin = Corrected/Uncorrected
  TH1D *trackPhi[2][2][knCollisionSystems][knMonteCarloTypes][nCentralityBins][nTrackPtBins+1];
  
  // String for finding inclusive histograms
  const char *inclusiveString[2] = {"","Inclusive"};   // 0 = Tracks in dijet events, 1 = Inclusive tracks
  const char *correctionString[2] = {"","Uncorrected"}; // 0 = Track correction included, 1 = No tracking corrections
  char namer[200];
  
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
    
    for(int iInclusive = 0; iInclusive < 2; iInclusive++){
      for(int iCorrection = 0; iCorrection < 2; iCorrection++){
        for(int iMonteCarloType = 0; iMonteCarloType < knMonteCarloTypes; iMonteCarloType++){
          for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
            
            // Read track pT histograms
            sprintf(namer,"track%s%s/track%s%sPt_SameEvent_C%d",inclusiveString[iInclusive],correctionString[iCorrection],inclusiveString[iInclusive],correctionString[iCorrection],iCentrality);
            trackPt[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality] = (TH1D*) closureFile[iSystem][iMonteCarloType]->Get(namer);
            
            // Read track eta histograms without pT cut
            sprintf(namer,"track%s%s/track%s%sEta_SameEvent_C%dT6",inclusiveString[iInclusive],correctionString[iCorrection],inclusiveString[iInclusive],correctionString[iCorrection],iCentrality);
            trackEta[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][nTrackPtBins] = (TH1D*) closureFile[iSystem][iMonteCarloType]->Get(namer);
            
            // Read track phi histograms without pT cut
            sprintf(namer,"track%s%s/track%s%sPhi_SameEvent_C%dT6",inclusiveString[iInclusive],correctionString[iCorrection],inclusiveString[iInclusive],correctionString[iCorrection],iCentrality);
            trackPhi[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][nTrackPtBins] = (TH1D*) closureFile[iSystem][iMonteCarloType]->Get(namer);
            
            for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
              
              // Read track eta histograms in pT bins
              sprintf(namer,"track%s%s/track%s%sEta_SameEvent_C%dT%d",inclusiveString[iInclusive],correctionString[iCorrection],inclusiveString[iInclusive],correctionString[iCorrection],iCentrality,iTrackPt);
              trackEta[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][nTrackPtBins] = (TH1D*) closureFile[iSystem][iMonteCarloType]->Get(namer);
              
              // Read track phi histograms in pT bins
              sprintf(namer,"track%s%s/track%s%sPhi_SameEvent_C%dT%d",inclusiveString[iInclusive],correctionString[iCorrection],inclusiveString[iInclusive],correctionString[iCorrection],iCentrality,iTrackPt);
              trackPhi[iInclusive][iCorrection][iSystem][iMonteCarloType][iCentrality][nTrackPtBins] = (TH1D*) closureFile[iSystem][iMonteCarloType]->Get(namer);
              
            } // Track pT loop
          } // Centrality loop
        } // Loop over Monte Carlo types
      } // Corrected - Uncorrected loop
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
    
    drawer->SetDefaultAppearanceSplitCanvas();
    drawer->SetLogY(logAxis);
    drawer->DrawHistogramToUpperPad(fMainHistogram,xTitle,yTitle," ");
    
    for(int iInclusive = 0; iInclusive < 2; iInclusive++){
      for(int iSystem = 0; iSystem < knCollisionSystems; iSystem++){
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          
          // No centrality binning for pp
          if(iSystem == kPp && iCentrality > 0) continue;
          
          drawer->CreateSplitCanvas();
          TH1D *trackPt[iInclusive][0][iSystem][knMonteCarloTypes][nCentralityBins];
          
        }
      } // Collision system loop
    } // Dijet- inclusive loop
  }
}
