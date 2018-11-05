#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)

/*
 * Set a nice drawing style for histogram
 *
 *  TH1 * histogram = Histogram needing a nice style
 *  double maxY = Maximum drawing range for Y-axis. Not set if negative.
 *  TString titleString = Title for the histogram
 *  const char *xTitle = Title for the x-axis
 *  bool disableFit = Disable the fit from spillover histograms
 */
void setSpilloverHistogramStyle(TH1 *histogram, double maxY, TString titleString, const char *xTitle, bool disableFit = false){
  histogram->GetXaxis()->SetRangeUser(-1.5,1.5);
  if(maxY > 0) histogram->GetYaxis()->SetRangeUser(0,maxY);
  histogram->SetStats(kFALSE);
  histogram->SetTitle(titleString);
  histogram->SetLabelSize(0.09,"xy");
  histogram->SetXTitle(xTitle);
  histogram->SetTitleSize(0.09,"x");
  
  if(disableFit){
    TF1 *gaussFit = histogram->GetFunction("gaussFit");
    gaussFit->SetLineWidth(0);
  }
}

/*
 * Plotter for QA related histograms
 *
 *  Currently only plotting QA for spillover correction. More QA plots to be added.
 */
void qaPlotter(){
  
  ///////////////////
  // Configuration //
  ///////////////////
  
  bool saveFigures = false;          // Save the figures to a file
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = true;    // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = false;     // Produce the correction for inclusive jet-track correlations
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,inclusiveJetTrack,inclusiveJetTrack};
  
  /////////////////
  // Config done //
  /////////////////
  
  // Open files containing the QA histograms
  TFile *spilloverFile = TFile::Open("data/spilloverCorrection_PbPbMC_skims_pfJets_noUncorrected_3eventsMixed_subeNon0_smoothedMixing_2018-11-05_QA.root");
  
  // Read the number of bins from histogram manager
  DijetHistogramManager *dummyManager = new DijetHistogramManager();
  const int nCentralityBins = dummyManager->GetNCentralityBins();
  const int nTrackPtBins = dummyManager->GetNTrackPtBins();
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,300};  // Bin borders for track pT
  
  // Define needed histograms
  TH1D *spilloverDeltaEtaProjection[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaPhiProjection[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  
  const char *histogramNames[2] = {"spilloverQA","dataReplica"};
  
  // Read the histograms from the QA file
  char histogramNamer[200];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        for(int iDataType = 0; iDataType < 2; iDataType++){
        
          sprintf(histogramNamer,"%sDeltaEtaProjection/%s_%sDeltaEtaProjection_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),histogramNames[iDataType],dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
          spilloverDeltaEtaProjection[iDataType][iJetTrack][iCentrality][iTrackPt] = (TH1D*) spilloverFile->Get(histogramNamer);
        
          sprintf(histogramNamer,"%sDeltaPhiProjection/%s_%sDeltaPhiProjection_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),histogramNames[iDataType],dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
          spilloverDeltaPhiProjection[iDataType][iJetTrack][iCentrality][iTrackPt] = (TH1D*) spilloverFile->Get(histogramNamer);
        } // data types
      } // track pT
    } // centrality
  } // jet track
  
  // Create canvases for different qa plots
  TCanvas *deltaEtaCanvas[DijetHistogramManager::knJetTrackCorrelations];
  TCanvas *deltaEtaComparisonCanvas[DijetHistogramManager::knJetTrackCorrelations];
  TCanvas *deltaPhiCanvas[DijetHistogramManager::knJetTrackCorrelations];
  TCanvas *deltaPhiComparisonCanvas[DijetHistogramManager::knJetTrackCorrelations];
  TString titleString;
  TPad *currentPad;
  char padNamer[100];
  
  // Draw spillover corrections with the fit
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    
    // Create one big canvas with a pad for each centrality and track pT bin
    sprintf(histogramNamer,"spilloverDeltaEta%d",iJetTrack);
    sprintf(padNamer,"Spillover deltaEta %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
    deltaEtaCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1000,1800);
    deltaEtaCanvas[iJetTrack]->Divide(nCentralityBins,nTrackPtBins);
    
    sprintf(histogramNamer,"comparisonDeltaEta%d",iJetTrack);
    sprintf(padNamer,"Comparison deltaEta %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
    deltaEtaComparisonCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1000,1800);
    deltaEtaComparisonCanvas[iJetTrack]->Divide(nCentralityBins,nTrackPtBins);
    
    sprintf(histogramNamer,"spilloverDeltaPhi%d",iJetTrack);
    sprintf(padNamer,"Spillover deltaPhi %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
    deltaPhiCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1000,1800);
    deltaPhiCanvas[iJetTrack]->Divide(nCentralityBins,nTrackPtBins);
    
    sprintf(histogramNamer,"comparisonDeltaPhi%d",iJetTrack);
    sprintf(padNamer,"Comparison deltaPhi %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
    deltaPhiComparisonCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1000,1800);
    deltaPhiComparisonCanvas[iJetTrack]->Divide(nCentralityBins,nTrackPtBins);
    
    gStyle->SetTitleSize(0.09,"t");
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        titleString = Form("Cent: %.0f-%.0f%% - Track pT: %.1f-%.1f GeV",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1],trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
        
        // Find the correct pad inside the canvas
        deltaEtaCanvas[iJetTrack]->cd(iCentrality+nCentralityBins*iTrackPt+1);
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.2);
        
        // Draw the histogram to canvas
        setSpilloverHistogramStyle(spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt],0.02,titleString,"#Delta#eta");
        spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->DrawCopy();
        
        // Change to comparison canvas
        deltaEtaComparisonCanvas[iJetTrack]->cd(iCentrality+nCentralityBins*iTrackPt+1);
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.2);
        
        // Draw the histograms without fits to the canvas
        setSpilloverHistogramStyle(spilloverDeltaEtaProjection[1][iJetTrack][iCentrality][iTrackPt],-1,titleString,"#Delta#eta",true);
        setSpilloverHistogramStyle(spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt],-1,titleString,"#Delta#eta",true);
        spilloverDeltaEtaProjection[1][iJetTrack][iCentrality][iTrackPt]->SetLineColor(kRed);
        spilloverDeltaEtaProjection[1][iJetTrack][iCentrality][iTrackPt]->Draw();
        spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->Draw("same");
        
        // Find the correct pad inside the canvas
        deltaPhiCanvas[iJetTrack]->cd(iCentrality+nCentralityBins*iTrackPt+1);
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.2);
        
        // Draw the histogram to canvas
        setSpilloverHistogramStyle(spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt],0.02,titleString,"#Delta#eta");
        spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->DrawCopy();
        
        // Change to comparison canvas
        deltaPhiComparisonCanvas[iJetTrack]->cd(iCentrality+nCentralityBins*iTrackPt+1);
        gPad->SetTopMargin(0.1);
        gPad->SetBottomMargin(0.2);
        
        // Draw the histograms without fits to the canvas
        setSpilloverHistogramStyle(spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt],-1,titleString,"#Delta#eta",true);
        setSpilloverHistogramStyle(spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt],-1,titleString,"#Delta#eta",true);
        spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt]->SetLineColor(kRed);
        spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt]->Draw();
        spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->Draw("same");
        
      } // track pT
    } // centrality
    
    // After all the canvases are filled, save them
    if(saveFigures) {
      sprintf(histogramNamer,"figures/%s_spilloverDeltaEta.pdf",dummyManager->GetJetTrackHistogramName(iJetTrack));
      deltaEtaCanvas[iJetTrack]->SaveAs(histogramNamer);
      
      sprintf(histogramNamer,"figures/%s_spilloverDeltaEtaComparison.pdf",dummyManager->GetJetTrackHistogramName(iJetTrack));
      deltaEtaComparisonCanvas[iJetTrack]->SaveAs(histogramNamer);
      
      sprintf(histogramNamer,"figures/%s_spilloverDeltaPhi.pdf",dummyManager->GetJetTrackHistogramName(iJetTrack));
      deltaPhiCanvas[iJetTrack]->SaveAs(histogramNamer);
      
      sprintf(histogramNamer,"figures/%s_spilloverDeltaPhiComparison.pdf",dummyManager->GetJetTrackHistogramName(iJetTrack));
      deltaEtaComparisonCanvas[iJetTrack]->SaveAs(histogramNamer);
      
    }
  } // jet track
  
}
