#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

/*
 * Read the flow fits from background histograms and put the results to graph
 *
 */
void graphFlow(){
  
  ///////////////////
  // Configuration //
  ///////////////////
  
  bool saveFigures = true;          // Save the figures to a file
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = false;   // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = false;    // Produce the correction for inclusive jet-track correlations
  bool drawSubleading = false;       // Draw the subleading side also
  
  bool drawHydjet = true;
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,(regularJetTrack&&drawSubleading),(uncorrectedJetTrack&&drawSubleading),(ptWeightedJetTrack&&drawSubleading),inclusiveJetTrack,inclusiveJetTrack};
  const char *titleAddition[DijetHistogramManager::knJetTrackCorrelations] = {"",", uncorrected",", $p_{\\mathrm{T}}$ weighted","",", uncorrected",", $p_{\\mathrm{T}}$ weighted","",", $p_{\\mathrm{T}}$ weighted"};
  
  /////////////////
  // Config done //
  /////////////////
  
  // Open files containing the background histograms with the fit
  TFile *backgroundFile = TFile::Open("data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_modifiedSeagull_noErrorJff_averagePeakMixing_processed_2019-08-13_fiveJobsMissing.root");
  
  TFile *hydjetFile = TFile::Open("data/PbPbMC_RecoGen_akFlowPuCsPfJets_noUncorr_xjBins_improvisedMixing_subeNon0_JECv4_noCorrections_processed_2019-08-09.root");
  
  // Read the number of bins from histogram manager
  DijetHistogramManager *dummyManager = new DijetHistogramManager(backgroundFile);
  
  // Histogram manager reads the number of track pT bins from the DijetCard, so they are automatically
  // updated based on how the backgroundFile is created
  const int nCentralityBins = 4;
  const int nTrackPtBins = 6;//dummyManager->GetNTrackPtBins();
  const int nFlowComponents = dummyManager->knFittedFlowComponents;
  
  // Todo: Read also bin borders from the Card
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,300};  // Bin borders for track pT
  double graphPoints[] = {0.85,1.5,2.5,3.5,6,10};    // x-axis points in flow graphs
  double graphErrors[] = {0,0,0,0,0,0};              // No errors for x-axis
  double meanPt[nCentralityBins][nTrackPtBins];      // Mean pT within each pT bin
  
  // Define needed histograms
  TH1D *backgroundDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *hydjetDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *trackPt[nCentralityBins];
  
  // Create a double array to which all the vn:s are read
  double masterFlowTable[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nFlowComponents][nTrackPtBins];
  double masterFlowError[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nFlowComponents][nTrackPtBins];
  double hydjetFlowTable[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nFlowComponents][nTrackPtBins];
  double hydjetFlowError[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nFlowComponents][nTrackPtBins];
  double jetFlowTable[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nFlowComponents][nTrackPtBins];
  double jetFlowError[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nFlowComponents][nTrackPtBins];
  TF1 *fourierFit;
  TF1 *hydjetFit;
  
  // Rough flow coefficients from paper arXiv:1201.3158. This does not have the same track pT or centrality bins
  // and values are read by eye, but they should be good enough to give an order of magnitude
  if(nCentralityBins != 4 || nFlowComponents != 4 || nTrackPtBins != 6){
    cout << "Error! wrong length for reference array! Check the code!" << endl;
    return;
  }
  double roughHadronicFlow[4][4][6] = {
     // pT 0.7-1,  1-2    2-3    3-4    4-8    8-300 GeV
     //              Centrality 0-10 %
    {{       1,     1,     1,     1,     1,     1    },   // v1
     {     0.045, 0.065, 0.08,  0.09,  0.08,  0.04   },   // v2
     {       1,     1,     1,     1,     1,     1    },   // v3
     {       1,     1,     1,     1,     1,     1    }},  // v4
     //              Centrality 10-30 %
    {{       1,     1,     1,     1,     1,     1    },   // v1
     {      0.1,  0.14,  0.17,  0.19,  0.155, 0.075  },   // v2
     {       1,     1,     1,     1,     1,     1    },   // v3
     {       1,     1,     1,     1,     1,     1    }},  // v4
     //              Centrality 30-50 %
    {{       1,     1,     1,     1,     1,     1    },   // v1
     {     0.14,  0.175, 0.21,  0.225, 0.18,   0.1   },   // v2
     {       1,     1,     1,     1,     1,     1    },   // v3
     {       1,     1,     1,     1,     1,     1    }},  // v4
     //              Centrality 50-100 %
    {{       1,     1,     1,     1,     1,     1    },   // v1
     {     0.125, 0.15,  0.19,  0.18,  0.16,  0.15   },   // v2
     {       1,     1,     1,     1,     1,     1    },   // v3
     {       1,     1,     1,     1,     1,     1    }}}; // v4
  
  // Read the vn:s from the background histograms that are in the background file
  char histogramNamer[200];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        sprintf(histogramNamer,"%s/%sDeltaPhi_Background_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        backgroundDeltaPhi[iJetTrack][iCentrality][iTrackPt] = (TH1D*) backgroundFile->Get(histogramNamer);
        
        fourierFit = backgroundDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetFunction("fourier");
        
        sprintf(histogramNamer,"%s/%sDeltaPhi_Background_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        hydjetDeltaPhi[iJetTrack][iCentrality][iTrackPt] = (TH1D*) hydjetFile->Get(histogramNamer);
        
        hydjetFit = hydjetDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetFunction("fourier");
        
        for(int iFlow = 0; iFlow < nFlowComponents; iFlow++){
          masterFlowTable[iJetTrack][iCentrality][iFlow][iTrackPt] = fourierFit->GetParameter(iFlow+1);
          masterFlowError[iJetTrack][iCentrality][iFlow][iTrackPt] = fourierFit->GetParError(iFlow+1);
          hydjetFlowTable[iJetTrack][iCentrality][iFlow][iTrackPt] = hydjetFit->GetParameter(iFlow+1);
          hydjetFlowError[iJetTrack][iCentrality][iFlow][iTrackPt] = hydjetFit->GetParError(iFlow+1);
          jetFlowTable[iJetTrack][iCentrality][iFlow][iTrackPt] = masterFlowTable[iJetTrack][iCentrality][iFlow][iTrackPt] / roughHadronicFlow[iCentrality][iFlow][iTrackPt];
          jetFlowError[iJetTrack][iCentrality][iFlow][iTrackPt] = masterFlowError[iJetTrack][iCentrality][iFlow][iTrackPt] / roughHadronicFlow[iCentrality][iFlow][iTrackPt];
        } // flow components
      } // track pT
    } // centrality
  } // jet track
  
  // Read track pT to determine a good place inside at pT bin to put the graph point
  int lowBin, highBin;
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    sprintf(histogramNamer,"track/trackPt_SameEvent_C%d",iCentrality);
    trackPt[iCentrality] = (TH1D*) backgroundFile->Get(histogramNamer);
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      lowBin = trackPt[iCentrality]->GetXaxis()->FindBin(trackPtBinBorders[iTrackPt]+0.001);
      if(trackPtBinBorders[iTrackPt+1] > 20){
        highBin = trackPt[iCentrality]->GetNbinsX();
      } else {
        highBin = trackPt[iCentrality]->GetXaxis()->FindBin(trackPtBinBorders[iTrackPt+1]-0.001);
      }
      trackPt[iCentrality]->GetXaxis()->SetRange(lowBin,highBin);
      meanPt[iCentrality][iTrackPt] = trackPt[iCentrality]->GetMean();
    }
  }
  
  // Define graphs that will take inside them the flow values
  TGraphErrors *flowGraph[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nFlowComponents];
  TGraphErrors *hydjetGraph[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nFlowComponents];
  TGraphErrors *jetFlowGraph[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nFlowComponents];
  TGraph *crudeGraph[nCentralityBins][nFlowComponents];
  TLegend *legend;
  char namer[100];
  char namerY[20];
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceGraph();
  
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iFlow = 1; iFlow < 2; iFlow++){
        flowGraph[iJetTrack][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins,meanPt[iCentrality],masterFlowTable[iJetTrack][iCentrality][iFlow],graphErrors,masterFlowError[iJetTrack][iCentrality][iFlow]);
        flowGraph[iJetTrack][iCentrality][iFlow]->SetMarkerStyle(21);
        flowGraph[iJetTrack][iCentrality][iFlow]->SetMarkerColor(kBlue);
        
        sprintf(namerY,"v_{%d}",iFlow+1);
        sprintf(namer,"Centrality: %.0f-%.0f%%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
        drawer->CreateCanvas(0,7,0,0.1,"p_{T} (GeV)",namerY,namer);
        flowGraph[iJetTrack][iCentrality][iFlow]->Draw("psame");
        
        hydjetGraph[iJetTrack][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins,meanPt[iCentrality],hydjetFlowTable[iJetTrack][iCentrality][iFlow],graphErrors,hydjetFlowError[iJetTrack][iCentrality][iFlow]);
        hydjetGraph[iJetTrack][iCentrality][iFlow]->SetMarkerStyle(21);
        hydjetGraph[iJetTrack][iCentrality][iFlow]->SetMarkerColor(kMagenta);
        if(drawHydjet) hydjetGraph[iJetTrack][iCentrality][iFlow]->Draw("psame");
        
        jetFlowGraph[iJetTrack][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins,meanPt[iCentrality],jetFlowTable[iJetTrack][iCentrality][iFlow],graphErrors,jetFlowError[iJetTrack][iCentrality][iFlow]);
        jetFlowGraph[iJetTrack][iCentrality][iFlow]->SetMarkerStyle(21);
        jetFlowGraph[iJetTrack][iCentrality][iFlow]->SetMarkerColor(kGreen+4);
        //jetFlowGraph[iJetTrack][iCentrality][iFlow]->Draw("psame");
        
        crudeGraph[iCentrality][iFlow] = new TGraph(nTrackPtBins,meanPt[iCentrality],roughHadronicFlow[iCentrality][iFlow]);
        crudeGraph[iCentrality][iFlow]->SetMarkerStyle(21);
        crudeGraph[iCentrality][iFlow]->SetMarkerColor(kRed);
        //crudeGraph[iCentrality][iFlow]->Draw("psame");
        
        legend = new TLegend(0.2,0.7,0.5,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
        sprintf(namer,"Track-jet V_{%d}",iFlow+1);
        legend->AddEntry(flowGraph[iJetTrack][iCentrality][iFlow],namer,"p");
        sprintf(namer,"Hydjet V_{%d}",iFlow+1);
        if(drawHydjet) legend->AddEntry(hydjetGraph[iJetTrack][iCentrality][iFlow],namer,"p");
        sprintf(namer,"Calculated v_{%d}^{jet}",iFlow+1);
        //legend->AddEntry(jetFlowGraph[iJetTrack][iCentrality][iFlow],namer,"p");
        sprintf(namer,"Rough v_{%d}^{track}",iFlow+1);
        //legend->AddEntry(crudeGraph[iCentrality][iFlow],namer,"p");
        legend->Draw();
        
      } // flow components
    } // centrality
  } // jet track
  
}
