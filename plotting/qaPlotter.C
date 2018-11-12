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
void setHistogramStyle(TH1 *histogram, double rangeX, double maxY, TString titleString, const char *xTitle, bool disableFit = false){
  histogram->GetXaxis()->SetRangeUser(-rangeX,rangeX);
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
 * Setup the ratio number to print to a canvas
 *
 *  TLegend* legend = Pointer to legend
 *  TH1D* dataHistogram = Histogram containing PbPb data
 *  TH1D* mcHistogram = Histogram containing hydjet simulation
 */
void setupRatioLegend(TLegend* legend, TH1* dataHistogram, TH1* mcHistogram){
  
  // Find the bins for signal region
  int lowXbin = dataHistogram->GetXaxis()->FindBin(-1.5);
  int highXbin = dataHistogram->GetXaxis()->FindBin(1.5);
  
  // Get the yields by integration and calculate ratio
  double dataYield = dataHistogram->Integral(lowXbin,highXbin,"width");
  double mcYield = mcHistogram->Integral(lowXbin,highXbin,"width");
  double yieldRatio = mcYield/dataYield;
  
  char namer[100];
  sprintf(namer,"Ratio: %.3f",yieldRatio);
  
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.08);legend->SetTextFont(62);
  legend->AddEntry((TObject*) 0,namer,"");
  legend->AddEntry(dataHistogram,"Data","l");
  legend->AddEntry(mcHistogram,"Hydjet","l");
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
  
  bool drawSpillover = true;        // Draw the QA plots for spillover correction
  bool drawSeagull = false;           // Draw the QA plots for seagull correction
  
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
  TFile *seagullFile = TFile::Open("data/dijetPbPb_pfJets_noInclusiveOrUncorrected_noCorrections_smoothedMixing_processed_2018-11-06_QA.root");
  
  // Read the number of bins from histogram manager
  DijetHistogramManager *dummyManager = new DijetHistogramManager();
  const int nCentralityBins = dummyManager->GetNCentralityBins();
  const int nTrackPtBins = dummyManager->GetNTrackPtBins();
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,300};  // Bin borders for track pT
  
  // Define needed histograms
  TH1D *spilloverDeltaEtaProjection[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaPhiProjection[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *seagullDeltaEtaWings[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  
  const char *histogramNames[2] = {"spilloverQA","dataReplica"};
  
  // Read the histograms from the QA file
  char histogramNamer[200];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        sprintf(histogramNamer,"seagullDeltaEta_%s/seagullDeltaEta_%s_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        seagullDeltaEtaWings[iJetTrack][iCentrality][iTrackPt] = (TH1D*) seagullFile->Get(histogramNamer);
        
        for(int iDataType = 0; iDataType < 2; iDataType++){
        
          sprintf(histogramNamer,"%sDeltaEtaProjection/%s_%sDeltaEtaProjection_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),histogramNames[iDataType],dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
          spilloverDeltaEtaProjection[iDataType][iJetTrack][iCentrality][iTrackPt] = (TH1D*) spilloverFile->Get(histogramNamer);
        
          sprintf(histogramNamer,"%sDeltaPhiProjection/%s_%sDeltaPhiProjection_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),histogramNames[iDataType],dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
          spilloverDeltaPhiProjection[iDataType][iJetTrack][iCentrality][iTrackPt] = (TH1D*) spilloverFile->Get(histogramNamer);
        } // data types
      } // track pT
    } // centrality
  } // jet track
  
  // ******************************************
  // **      Drawing spillover plots         **
  // ******************************************
  
  // Variables that can be reused for all QA plots
  TPad *currentPad;
  TLegend *legend;
  char padNamer[100];
  TString titleString;
  
  if(drawSpillover){
    
    // Create canvases for different qa plots
    TCanvas *deltaEtaCanvas[DijetHistogramManager::knJetTrackCorrelations];
    TCanvas *deltaEtaComparisonCanvas[DijetHistogramManager::knJetTrackCorrelations];
    TCanvas *deltaPhiCanvas[DijetHistogramManager::knJetTrackCorrelations];
    TCanvas *deltaPhiComparisonCanvas[DijetHistogramManager::knJetTrackCorrelations];
    
    // Set a good title size for big canvases
    gStyle->SetTitleSize(0.09,"t");
    
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
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          titleString = Form("Cent: %.0f-%.0f%% - Track pT: %.1f-%.1f GeV",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1],trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
          
          // Find the correct pad inside the canvas
          deltaEtaCanvas[iJetTrack]->cd(iCentrality+nCentralityBins*iTrackPt+1);
          gPad->SetTopMargin(0.1);
          gPad->SetBottomMargin(0.2);
          
          // Draw the histogram to canvas
          setHistogramStyle(spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt],1.5,0.02,titleString,"#Delta#eta");
          spilloverDeltaEtaProjection[1][iJetTrack][iCentrality][iTrackPt]->SetLineColor(kRed);
          spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->DrawCopy();
          
          legend = new TLegend(0.61,0.6,0.9,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.08);legend->SetTextFont(62);
          legend->AddEntry(spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt],"Hydjet","l");
          legend->AddEntry(spilloverDeltaEtaProjection[1][iJetTrack][iCentrality][iTrackPt],"Fit","l");
          legend->Draw();
          
          // Change to comparison canvas
          deltaEtaComparisonCanvas[iJetTrack]->cd(iCentrality+nCentralityBins*iTrackPt+1);
          gPad->SetTopMargin(0.1);
          gPad->SetBottomMargin(0.2);
          
          // Draw the histograms without fits to the canvas
          setHistogramStyle(spilloverDeltaEtaProjection[1][iJetTrack][iCentrality][iTrackPt],1.5,-1,titleString,"#Delta#eta",true);
          setHistogramStyle(spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt],1.5,-1,titleString,"#Delta#eta",true);
          spilloverDeltaEtaProjection[1][iJetTrack][iCentrality][iTrackPt]->Draw();
          spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->Draw("same");
          
          legend = new TLegend(0.61,0.6,0.9,0.9);
          setupRatioLegend(legend,spilloverDeltaEtaProjection[1][iJetTrack][iCentrality][iTrackPt],spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]);
          legend->Draw();
          
          // Find the correct pad inside the canvas
          deltaPhiCanvas[iJetTrack]->cd(iCentrality+nCentralityBins*iTrackPt+1);
          gPad->SetTopMargin(0.1);
          gPad->SetBottomMargin(0.2);
          
          // Draw the histogram to canvas
          setHistogramStyle(spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt],1.5,0.02,titleString,"#Delta#eta");
          spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt]->SetLineColor(kRed);
          spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->DrawCopy();
          
          legend = new TLegend(0.61,0.6,0.9,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.08);legend->SetTextFont(62);
          legend->AddEntry(spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt],"Hydjet","l");
          legend->AddEntry(spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt],"Fit","l");
          legend->Draw();
          
          // Change to comparison canvas
          deltaPhiComparisonCanvas[iJetTrack]->cd(iCentrality+nCentralityBins*iTrackPt+1);
          gPad->SetTopMargin(0.1);
          gPad->SetBottomMargin(0.2);
          
          // Draw the histograms without fits to the canvas
          setHistogramStyle(spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt],1.5,-1,titleString,"#Delta#eta",true);
          setHistogramStyle(spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt],1.5,-1,titleString,"#Delta#eta",true);
          spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt]->Draw();
          spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->Draw("same");
          
          legend = new TLegend(0.61,0.6,0.9,0.9);
          setupRatioLegend(legend,spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt],spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]);
          legend->Draw();
          
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
        
      } // saving figures
    } // jet track
    
  } // Drawing spillover correction
  
  // ******************************************
  // **       Drawing seagull plots          **
  // ******************************************
  
  if(drawSeagull){
    
    // Define canvases for the seagull correction
    TCanvas *seagullCanvas[DijetHistogramManager::knJetTrackCorrelations];
    
    // Set a good title size for big canvases
    gStyle->SetTitleSize(0.09,"t");
    
    // Draw spillover corrections with the fit
    for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
      if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
      
      // Create one big canvas with a pad for each centrality and track pT bin
      sprintf(histogramNamer,"seagullDeltaEta%d",iJetTrack);
      sprintf(padNamer,"Seagull deltaEta %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
      seagullCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1000,1800);
      seagullCanvas[iJetTrack]->Divide(nCentralityBins,nTrackPtBins);
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          titleString = Form("Cent: %.0f-%.0f%% - Track pT: %.1f-%.1f GeV",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1],trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
          
          // Find the correct pad inside the canvas
          seagullCanvas[iJetTrack]->cd(iCentrality+nCentralityBins*iTrackPt+1);
          gPad->SetTopMargin(0.1);
          gPad->SetBottomMargin(0.2);
          
          // Draw the histogram to canvas
          setHistogramStyle(seagullDeltaEtaWings[iJetTrack][iCentrality][iTrackPt],3,-1,titleString,"#Delta#eta");
          seagullDeltaEtaWings[iJetTrack][iCentrality][iTrackPt]->DrawCopy();
          
        } // track pT loop
      } // centrality loop
      
      // After all the canvases are filled, save them
      if(saveFigures) {
        sprintf(histogramNamer,"figures/%s_seagullDeltaEta.pdf",dummyManager->GetJetTrackHistogramName(iJetTrack));
        seagullCanvas[iJetTrack]->SaveAs(histogramNamer);
      } // saving figures
    } // jet-track loop
    
  } // Drawing the seagull correction
}
