#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"
#include "DijetMethods.h"

/*
 * Read the flow fits from background histograms and put the results to graph
 *
 */
void hadronFlow(){
  
  ///////////////////
  // Configuration //
  ///////////////////
  
  const bool saveFigures = false;     // Save the figures to a file
  
  const bool inclusiveDihadron = false;  // Dihadron correlations from all events
  const bool dijetDihadron = true;       // Dihadron correlations from dijet events
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {dijetDihadron,false,false,false,false,false,inclusiveDihadron,false};
  const char *titleAddition[DijetHistogramManager::knJetTrackCorrelations] = {"",", uncorrected",", $p_{\\mathrm{T}}$ weighted","",", uncorrected",", $p_{\\mathrm{T}}$ weighted","",", $p_{\\mathrm{T}}$ weighted"};
  
  const bool drawFits = true;
  
  /////////////////
  // Config done //
  /////////////////
  
  // Open files containing the background histograms with the fit
  TString backgroundFileName = "data/dihadronPbPb2018_sameTriggerAssoc_5eventMixed_noCorrections_processed_2020-06-18_smallStats.root";
  // dihadronPbPb2018_sameTriggerAssoc_5eventMixed_noCorrections_processed_2020-06-18_smallStats.root
  // data/dihadronPbPb2018_sameTriggerAssoc_xjBins_improvisedMixing_preprocessed_2020-06-18.root
  // data/dihadronPbPb2018_trigger2-3_processed_2020-03-23.root
  // data/dihadronPbPb2018_trigger3-4_processed_2020-03-23.root
  // data/dihadronPbPb2018_trigger4-8_processed_2020-03-23.root
  // "data/dijetPbPb2018_highForest_minBias_noMixing_dihadron_oldBinning_combinedCent_processed_2019-08-26.root";
  // data/dijetPbPb2018_highForest_minBias_improvisedMixing_dihadron_oldBinning_largeEtaGap_processed_2019-08-26.root
  // data/dijetPbPb2018_highForest_minBias_improvisedMixing_dihadron_oldBinning_processed_2019-08-26.root
  // data/dijetPbPb2018_highForest_minBias_improvisedMixing_dihadron_processed_2019-08-23.root
  // data/dijetPbPb2018_highForest_minBias_improvisedMixing_dihadron_processed_2019-08-22.root
  TFile *backgroundFile = TFile::Open(backgroundFileName);
  
  // Read the number of bins from histogram manager
  DijetHistogramManager *dummyManager = new DijetHistogramManager(backgroundFile);
  
  // Histogram manager reads the number of track pT bins from the DijetCard, so they are automatically
  // updated based on how the backgroundFile is created
  const int nCentralityBins = 3;
  const int nTrackPtBins = 4;//dummyManager->GetNTrackPtBins();
  const int nFlowComponents = dummyManager->knFittedFlowComponents;
  const int triggerBin = -1;   // Bin used for normalizing dihadron V2 to hadron v2. For -1, each bin is normalized by the square root of that bin
  
  const bool useSameEvent = false;  // true = Use same event long range, false = Use mixed event corrected long range
  
  // Todo: Read also bin borders from the Card
  //double centralityBinBorders[] = {0,5,10,20,90,100};  // Bin borders for centrality
  //double trackPtBinBorders[] = {1,1.5,2,2.5,3,3.5,4,4.5,5,8,300};  // Bin borders for track pT
  double centralityBinBorders[] = {0, 10, 30, 50, 90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double graphPoints[] = {0.85,1.5,2.5,6,10,14};    // x-axis points in flow graphs
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
  /*if(nCentralityBins != 4 || nFlowComponents != 4 || nTrackPtBins != 5){
    cout << "Error! wrong length for reference array! Check the code!" << endl;
    return;
  }*/
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
  
  // v2 from the v2{SP} method from the paper arXiv:1702.00630
  double flowCMS5TeVCentral[6][11] = {
// pT  1.1    1.4    1.7    2.2    2.7    3.2    3.7    4.4    5.4    6.4    7.4
    {0.0381,0.0441,0.0512,0.0603,0.0674,0.0712,0.0710,0.0656,0.0513,0.0381,0.0353},   // Centrality = 0-5 %
    {0.0645,0.0748,0.0876,0.1036,0.1160,0.1234,0.1242,0.1156,0.0914,0.0690,0.0543},   // Centrality = 5-10 %
    {0.0945,0.1095,0.1278,0.1503,0.1669,0.1748,0.1742,0.1609,0.1282,0.1033,0.0831},   // Centrality = 10-20 %
    {0.1250,0.1445,0.1677,0.1949,0.2125,0.2187,0.2149,0.1957,0.1571,0.1254,0.1052},   // Centrality = 20-30 %
    {0.1454,0.1675,0.1929,0.2205,0.2357,0.2383,0.2295,0.2069,0.1667,0.1357,0.1140},   // Centrality = 30-40 %
    {0.1570,0.1799,0.2049,0.2294,0.2402,0.2376,0.2251,0.2017,0.1650,0.1367,0.1190}    // Centrality = 40-50 %
  };
  
  double errorCMS5TeVCentral[6][11] = {
// pT 1.1    1.4    1.7    2.2    2.7    3.2    3.7    4.4    5.4    6.4    7.4
    {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.002, 0.002},  // Centrality = 0-5 %
    {0.002, 0.002, 0.002, 0.002, 0.002, 0.003, 0.003, 0.003, 0.003, 0.003, 0.003},  // Centrality = 5-10 %
    {0.002, 0.002, 0.003, 0.003, 0.003, 0.004, 0.004, 0.003, 0.003, 0.003, 0.003},  // Centrality = 10-20 %
    {0.003, 0.003, 0.003, 0.004, 0.004, 0.004, 0.004 ,0.004, 0.004, 0.004, 0.004},  // Centrality = 20-30 %
    {0.003, 0.003, 0.004, 0.004, 0.005, 0.005, 0.005, 0.005, 0.005, 0.004, 0.004},  // Centrality = 30-40 %
    {0.003, 0.004, 0.004, 0.005, 0.005, 0.005, 0.005, 0.005, 0.004, 0.004, 0.004}   // Centrality = 40-50 %
  };
  
  double averagedFlowCMS5TeVCentral[3][11];
  double averagedErrorCMS5TeVCentral[3][11];
  
  // Average the centrality bins to match the binning in this analysis
  for(int iCentrality = 0; iCentrality < 6; iCentrality = iCentrality+2){
    for(int iTrackPt = 0; iTrackPt < 11; iTrackPt++){
      averagedFlowCMS5TeVCentral[iCentrality/2][iTrackPt] = (flowCMS5TeVCentral[iCentrality][iTrackPt] + flowCMS5TeVCentral[iCentrality+1][iTrackPt]) / 2.0;
      
      averagedErrorCMS5TeVCentral[iCentrality/2][iTrackPt] = (errorCMS5TeVCentral[iCentrality][iTrackPt] + errorCMS5TeVCentral[iCentrality+1][iTrackPt]) / 2.0;
    }
  }
  
  
  double xAxisCMS5TeVCentral[] = {1.1,1.4,1.7,2.2,2.7,3.2,3.7,4.4,5.4,6.4,7.4};
  double xAxisErrorCMS5TeVCentral[] = {0,0,0,0,0,0,0,0,0,0,0};
  
  // Read the vn:s from the background histograms that are in the background file
  char histogramNamer[200];
  TH2D *helperHistogram;
  DijetMethods *fitter = new DijetMethods();
  int nBins;
  double binError, errorScale;
  
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Use same event long range for getting the hadron v2
        if(useSameEvent){
          sprintf(histogramNamer,"%s/%sDeltaEtaDeltaPhi_SameEvent_C%dT%d", dummyManager->GetJetTrackHistogramName(iJetTrack), dummyManager->GetJetTrackHistogramName(iJetTrack), iCentrality, iTrackPt);
          helperHistogram = (TH2D*) backgroundFile->Get(histogramNamer);
          
          fitter->SubtractBackground(helperHistogram, helperHistogram, 4, false);
          helperHistogram = fitter->GetBackground();
          
          nBins = helperHistogram->GetNbinsY();
          backgroundDeltaPhi[iJetTrack][iCentrality][iTrackPt] = helperHistogram->ProjectionX(Form("sameLong%d%d%d", iJetTrack, iCentrality, iTrackPt), 1, nBins);  // Exclude underflow and overflow bins by specifying range
          
          backgroundDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Scale(helperHistogram->GetYaxis()->GetBinWidth(1) );  // For correct normalization, need to divide out deltaEta bin width of the two-dimensional histogram
          
          /*
           * Do error scaling and for the background deltaPhi distribution
           *
           * Error scaling is needed because information only from 1.5 < |deltaEta| < 2.5 is used to determine the background.
           * When doing projection, root by default scaled the histogram with square root of bins projected over
           * The whole range of the analysis is |deltaEta| < 4, which means that four times the bins used to determine
           * the background are used in normalixing the error. We can fix this be scaling the error up by sqrt(4) = 2.
           * This number is provided by the DijetMethods, since if we change the limits described above, the number is
           * automatically adjusted for the new limits inside DijetMethods.
           */
          
          for(int iBin = 1; iBin < backgroundDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetNbinsX(); iBin++){
            binError = backgroundDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetBinError(iBin);
            errorScale = fitter->GetBackgroundErrorScalingFactor();
            backgroundDeltaPhi[iJetTrack][iCentrality][iTrackPt]->SetBinError(iBin, binError*errorScale);
          }
          
          fitter->FourierFit(backgroundDeltaPhi[iJetTrack][iCentrality][iTrackPt], 4);
          fourierFit = backgroundDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetFunction("fourier");
          
        // Use mixed event corrected long range for getting the hadron v2
        } else {
          sprintf(histogramNamer,"%s/%sDeltaPhi_Background_C%dT%d", dummyManager->GetJetTrackHistogramName(iJetTrack), dummyManager->GetJetTrackHistogramName(iJetTrack), iCentrality, iTrackPt);
          backgroundDeltaPhi[iJetTrack][iCentrality][iTrackPt] = (TH1D*) backgroundFile->Get(histogramNamer);
          cout << "iJetTrack = " << iJetTrack << " iCentrality = " << iCentrality << " iTrackPt = " << iTrackPt << endl;
          fourierFit = backgroundDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetFunction("fourier");
        }
        
        for(int iFlow = 0; iFlow < nFlowComponents; iFlow++){
          masterFlowTable[iJetTrack][iCentrality][iFlow][iTrackPt] = fourierFit->GetParameter(iFlow+1);
          masterFlowError[iJetTrack][iCentrality][iFlow][iTrackPt] = fourierFit->GetParError(iFlow+1);
        } // flow components
      } // track pT
    } // centrality
  } // jet track
  
  // Calculate the hadronic v2
  double flowNormalizer = 0;
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iFlow = 0; iFlow < nFlowComponents; iFlow++){
        if(triggerBin > -1) flowNormalizer = TMath::Sqrt(masterFlowTable[iJetTrack][iCentrality][iFlow][triggerBin]);
        if(iFlow == 1){
          cout << "Flow normalizer " << iCentrality << " is: " << flowNormalizer << endl;
        }
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          if(triggerBin < 0) flowNormalizer = TMath::Sqrt(masterFlowTable[iJetTrack][iCentrality][iFlow][iTrackPt]);
          jetFlowTable[iJetTrack][iCentrality][iFlow][iTrackPt] = masterFlowTable[iJetTrack][iCentrality][iFlow][iTrackPt] / flowNormalizer;
          jetFlowError[iJetTrack][iCentrality][iFlow][iTrackPt] = masterFlowError[iJetTrack][iCentrality][iFlow][iTrackPt] / flowNormalizer;
          if(iFlow == 1){
            cout << "iCentrality: " << iCentrality << " iTrackPt: " << " v2 = " << jetFlowTable[iJetTrack][iCentrality][iFlow][iTrackPt] << endl;
          }
        } // track pT
      } // flow components
    } // centrality
  } // jet track
  
  // Read track pT to determine a good place inside at pT bin to put the graph point
  int lowBin, highBin;
  double defaultPt[] = {0.835, 1.353, 2.36, 3.375, 6, 10, 14};
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    sprintf(histogramNamer,"track/trackPt_SameEvent_C%d",iCentrality);
    trackPt[iCentrality] = (TH1D*) backgroundFile->Get(histogramNamer);
    if(trackPt[iCentrality] == NULL){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        meanPt[iCentrality][iTrackPt] = defaultPt[iTrackPt];
      }
    } else {
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
  }
  
  // Define graphs that will take inside them the flow values
  TGraphErrors *flowGraph[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nFlowComponents];
  TGraphErrors *jetFlowGraph[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nFlowComponents];
  TGraphErrors *cms5TeVGraph[nCentralityBins][nFlowComponents];
  TGraphErrors *crudeGraph[nCentralityBins][nFlowComponents];
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
        drawer->CreateCanvas(0,7,0,0.35,"p_{T} (GeV)",namerY,namer);
        flowGraph[iJetTrack][iCentrality][iFlow]->Draw("psame");
        
        jetFlowGraph[iJetTrack][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins,meanPt[iCentrality],jetFlowTable[iJetTrack][iCentrality][iFlow],graphErrors,jetFlowError[iJetTrack][iCentrality][iFlow]);
        jetFlowGraph[iJetTrack][iCentrality][iFlow]->SetMarkerStyle(21);
        jetFlowGraph[iJetTrack][iCentrality][iFlow]->SetMarkerColor(kGreen+4);
        jetFlowGraph[iJetTrack][iCentrality][iFlow]->Draw("psame");
        
        cms5TeVGraph[iCentrality][iFlow] = new TGraphErrors(11,xAxisCMS5TeVCentral,averagedFlowCMS5TeVCentral[iCentrality],xAxisErrorCMS5TeVCentral,averagedErrorCMS5TeVCentral[iCentrality]);
        cms5TeVGraph[iCentrality][iFlow]->SetMarkerStyle(21);
        cms5TeVGraph[iCentrality][iFlow]->SetMarkerColor(kRed);
        cms5TeVGraph[iCentrality][iFlow]->Draw("psame");
        
        crudeGraph[iCentrality][iFlow] = new TGraphErrors(nTrackPtBins,meanPt[iCentrality],roughHadronicFlow[iCentrality][1],graphErrors,jetFlowError[iJetTrack][iCentrality][iFlow]);
        crudeGraph[iCentrality][iFlow]->SetMarkerStyle(21);
        crudeGraph[iCentrality][iFlow]->SetMarkerColor(kMagenta);
        crudeGraph[iCentrality][iFlow]->Draw("psame");
        
        legend = new TLegend(0.2,0.7,0.5,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
        sprintf(namer,"Dihadron V_{%d} , dijet",iFlow+1);
        legend->AddEntry(flowGraph[iJetTrack][iCentrality][iFlow],namer,"p");
        sprintf(namer,"Calculated v_{%d}^{track}",iFlow+1);
        legend->AddEntry(jetFlowGraph[iJetTrack][iCentrality][iFlow],namer,"p");
        sprintf(namer,"CMS 1702.00630 v_{%d}^{track}{SP}",iFlow+1);
        legend->AddEntry(cms5TeVGraph[iCentrality][iFlow],namer,"p");
        sprintf(namer,"CMS 1201.3158 v_{%d}^{track}",iFlow+1);
        legend->AddEntry(crudeGraph[iCentrality][iFlow],namer,"p");
        legend->Draw();
        
      } // flow components
    } // centrality
  } // jet track
  
  // Draw the fits in each bin to see that nothing crazy is happening
  if(drawFits){
    
    drawer->Reset();
    
    for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
      if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        
        // Collect the y-axis information to arrays
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          drawer->DrawHistogram(backgroundDeltaPhi[iJetTrack][iCentrality][iTrackPt], "#Delta#varphi", "#frac{dN}{d#Delta#varphi}", " ");
          
          legend = new TLegend(0.2,0.7,0.5,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
          legend->SetHeader(Form("%s", dummyManager->GetJetTrackHistogramName(iJetTrack)));
          legend->AddEntry((TObject*) 0, Form("C = %.0f-%.0f", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "");
          legend->AddEntry((TObject*) 0, Form("%.1f < pT < %.1f GeV", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]), "");
          
          legend->Draw();
          
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
    
  }
  
}
