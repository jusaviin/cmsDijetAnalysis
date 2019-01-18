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
  
  bool saveFigures = false;          // Save the figures to a file
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = false;   // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = false;    // Produce the correction for inclusive jet-track correlations
  bool drawSubleading = false;       // Draw the subleading side also
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,(regularJetTrack&&drawSubleading),(uncorrectedJetTrack&&drawSubleading),(ptWeightedJetTrack&&drawSubleading),inclusiveJetTrack,inclusiveJetTrack};
  const char *titleAddition[DijetHistogramManager::knJetTrackCorrelations] = {"",", uncorrected",", $p_{\\mathrm{T}}$ weighted","",", uncorrected",", $p_{\\mathrm{T}}$ weighted","",", $p_{\\mathrm{T}}$ weighted"};
  
  /////////////////
  // Config done //
  /////////////////
  
  // Open files containing the background histograms with the fit
  TFile *backgroundFile = TFile::Open("data/dijetPbPb_skims_pfJets_noUncorr_improvedPoolMixing_noJetLimit_noCorrections_processed_2019-01-09.root");
  
  // Read the number of bins from histogram manager
  DijetHistogramManager *dummyManager = new DijetHistogramManager();
  
  // In principle it would be better to read these from card, since in that case the number of bins can be changed
  // in HistogramManager without breaking background compatibility. But this work also for fixed number of bins.
  const int nCentralityBins = dummyManager->knCentralityBins;
  const int nTrackPtBins = dummyManager->knTrackPtBins;
  const int nFlowComponents = dummyManager->knFittedFlowComponents;
  
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,300};  // Bin borders for track pT
  double graphPoints[] = {0.85,1.5,2.5,3.5,6,10};    // x-axis points in flow graphs
  double graphErrors[] = {0,0,0,0,0,0};              // No errors for x-axis
  double meanPt[nCentralityBins][nTrackPtBins];      // Mean pT within each pT bin
  
  // Define needed histograms
  TH1D *backgroundDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *trackPt[nCentralityBins];
  
  // Create a double array to which all the vn:s are read
  double masterFlowTable[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nFlowComponents][nTrackPtBins];
  double masterFlowError[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nFlowComponents][nTrackPtBins];
  double jetFlowTable[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nFlowComponents][nTrackPtBins];
  double jetFlowError[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nFlowComponents][nTrackPtBins];
  TF1 *fourierFit;
  
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
     {     0.045, 0.065, 0.08,  0.09,  0.06,  0.04   },   // v2
     {       1,     1,     1,     1,     1,     1    },   // v3
     {       1,     1,     1,     1,     1,     1    }},  // v4
     //              Centrality 10-30 %
    {{       1,     1,     1,     1,     1,     1    },   // v1
     {      0.1,  0.14,  0.17,  0.19,  0.12,  0.075  },   // v2
     {       1,     1,     1,     1,     1,     1    },   // v3
     {       1,     1,     1,     1,     1,     1    }},  // v4
     //              Centrality 30-50 %
    {{       1,     1,     1,     1,     1,     1    },   // v1
     {     0.14,  0.175, 0.21,  0.225, 0.155,  0.1   },   // v2
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
        
        for(int iFlow = 0; iFlow < nFlowComponents; iFlow++){
          masterFlowTable[iJetTrack][iCentrality][iFlow][iTrackPt] = fourierFit->GetParameter(iFlow+1);
          masterFlowError[iJetTrack][iCentrality][iFlow][iTrackPt] = fourierFit->GetParError(iFlow+1);
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
        drawer->CreateCanvas(0,7,0,0.6,"p_{T} (GeV)",namerY,namer);
        flowGraph[iJetTrack][iCentrality][iFlow]->Draw("psame");
        
        jetFlowGraph[iJetTrack][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins,meanPt[iCentrality],jetFlowTable[iJetTrack][iCentrality][iFlow],graphErrors,jetFlowError[iJetTrack][iCentrality][iFlow]);
        jetFlowGraph[iJetTrack][iCentrality][iFlow]->SetMarkerStyle(21);
        jetFlowGraph[iJetTrack][iCentrality][iFlow]->SetMarkerColor(kGreen+4);
        jetFlowGraph[iJetTrack][iCentrality][iFlow]->Draw("psame");
        
        crudeGraph[iCentrality][iFlow] = new TGraph(nTrackPtBins,meanPt[iCentrality],roughHadronicFlow[iCentrality][iFlow]);
        crudeGraph[iCentrality][iFlow]->SetMarkerStyle(21);
        crudeGraph[iCentrality][iFlow]->SetMarkerColor(kRed);
        crudeGraph[iCentrality][iFlow]->Draw("psame");
        
        legend = new TLegend(0.2,0.7,0.5,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
        sprintf(namer,"V_{%d} from fit",iFlow+1);
        legend->AddEntry(flowGraph[iJetTrack][iCentrality][iFlow],namer,"p");
        sprintf(namer,"Calculated v_{%d}^{jet}",iFlow+1);
        legend->AddEntry(jetFlowGraph[iJetTrack][iCentrality][iFlow],namer,"p");
        sprintf(namer,"Rough v_{%d}^{track}",iFlow+1);
        legend->AddEntry(crudeGraph[iCentrality][iFlow],namer,"p");
        legend->Draw();
        
      } // flow components
    } // centrality
  } // jet track
  
}
