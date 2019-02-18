#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"
#include "JDrawer.h"
#include <tuple>

/*
 * Get maximum and minimum drawing range from a histogram
 *
 *  TH1* histogram = Histogram from which drawing range in searched
 *
 *  return: Minimum and maximum value for y-axis
 */
std::tuple<double,double> getDrawMinMax(TH1* histogram){
  
  // Helper variables
  double minValue, maxValue, minMaxError;
  int minBin, maxBin;
  
  // Find a good drawing range for the distributions
  maxBin = histogram->GetMaximumBin();
  minBin = histogram->GetMinimumBin();
  maxValue = histogram->GetBinContent(maxBin);
  minMaxError = histogram->GetBinError(maxBin);
  maxValue = maxValue + minMaxError*2;
  minValue = histogram->GetBinContent(minBin);
  minMaxError = histogram->GetBinError(minBin);
  minValue = minValue - minMaxError*2;
  
  return std::make_tuple(minValue,maxValue);
}

/*
 * Set titles for the histogram
 *
 *  TH1 * histogram = Histogram needing a titles
 *  TString titleString = Title for the histogram
 *  const char *xTitle = Title for the x-axis
 */
void setHistogramTitles(TH1 *histogram, TString titleString, const char *xTitle){
  
  histogram->SetStats(kFALSE);
  histogram->SetTitle(titleString);
  histogram->SetLabelSize(0.09,"xy");
  histogram->SetXTitle(xTitle);
  histogram->SetTitleSize(0.09,"x");
  
}

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
  if(rangeX > 0) histogram->GetXaxis()->SetRangeUser(-rangeX,rangeX);
  if(maxY > 0) histogram->GetYaxis()->SetRangeUser(0,maxY);
  
  setHistogramTitles(histogram,titleString,xTitle);
  
  if(disableFit){
    TF1 *gaussFit = histogram->GetFunction("gaussFit");
    gaussFit->SetLineWidth(0);
  } else {
    TF1 *gaussFit = histogram->GetFunction("gaussFit");
    gaussFit->SetRange(-1.5,1.5);
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
 * Print a slide with background overlap information to console
 *
 *
 */
void printBackroundSumSlide(double leadingValues[4][6], double leadingErrors[4][6], double subleadingValues[4][6], double subleadingErrors[4][6]){
  
  // Create a histogram manager to facilitate binning info exchange
  const int nTrackPtBins = 6;
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,300};  // Bin borders for track pT
  
  char namer[100];
  
  cout << endl;
  cout << "\\begin{frame}" << endl;
  cout << "\\frametitle{Leading background/subleading background}" << endl;
  cout << "\\begin{center}" << endl;
  cout << "  \\begin{tabular}{ccccc}" << endl;
  cout << "    \\toprule" << endl;
  cout << "    $p_{\\mathrm{T}} (GeV)$ & C: 0-10 \\% & C: 10-30 \\% & C: 30-50 \\% & C: 50-100 \\% \\\\" << endl;
  cout << "    \\midrule" << endl;
  
  // Set the correct precision for printing floating point numbers
  cout << fixed << setprecision(3);
  
  // Print one line for each track pT bin
  for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
    sprintf(namer,"    %.1f-%.1f",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
    cout << namer;
    for(int iCentrality = 0; iCentrality < 4; iCentrality++){
      cout << " & $" << leadingValues[iCentrality][iTrackPt]/subleadingValues[iCentrality][iTrackPt] << "\\pm" << TMath::Sqrt(TMath::Power(leadingErrors[iCentrality][iTrackPt]/subleadingValues[iCentrality][iTrackPt],2)+TMath::Power(leadingValues[iCentrality][iTrackPt]*subleadingErrors[iCentrality][iTrackPt]/TMath::Power(subleadingValues[iCentrality][iTrackPt],2),2)) << "$";
    }
    cout << " \\\\" << endl;
  }
  
  cout << "    \\bottomrule" << endl;
  cout << "  \\end{tabular}" << endl;
  cout << "\\end{center}" << endl;
  cout << "\\end{frame}" << endl;
}

/*
 * Print a slide with background overlap information to console
 *
 *
 */
void printBackgroundOverlapSlide(bool subleadingOverlap, double binValues[6][6], double binErrors[6][6], double averageValues[6], double averageErrors[6], TString centrality){
  
  // Create a histogram manager to facilitate binning info exchange
  const int nTrackPtBins = 6;
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,300};  // Bin borders for track pT
  
  // Set up correct naming based on background type
  int binNumbers[3] = {101,102,103};
  char namer[100];
  sprintf(namer,"\\frametitle{Subleading bg/leading overlap %s}",centrality.Data());
  int adder = 3;
  
  if(subleadingOverlap){
    sprintf(namer,"\\frametitle{Leading bg/subleading overlap %s}",centrality.Data());
    binNumbers[0] = 98; binNumbers[1] = 99; binNumbers[2] = 100;
    adder = 0;
  }
  
  
  cout << endl;
  cout << "\\begin{frame}" << endl;
  cout << namer << endl;
  cout << "\\begin{center}" << endl;
  cout << "  \\begin{tabular}{ccccc}" << endl;
  cout << "    \\toprule" << endl;
  cout << "    $p_{\\mathrm{T}} (GeV)$ & Bin "<< binNumbers[0] << " & Bin " << binNumbers[1] << " & Bin " << binNumbers[2] << " & Average \\\\" << endl;
  cout << "    \\midrule" << endl;

  // Set the correct precision for printing floating point numbers
  cout << fixed << setprecision(3);
  
  // Print one line for each track pT bin
  for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
    sprintf(namer,"    %.1f-%.1f",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
    cout << namer;
    for(int iBin = 0; iBin < 3; iBin++){
      cout << " & $" << binValues[iTrackPt][iBin+adder] << "\\pm" << binErrors[iTrackPt][iBin+adder] << "$";
    }
    cout << " & $" << averageValues[iTrackPt] << "\\pm" << averageErrors[iTrackPt] << "$";
    cout << " \\\\" << endl;
  }
  
  cout << "    \\bottomrule" << endl;
  cout << "  \\end{tabular}" << endl;
  cout << "\\end{center}" << endl;
  cout << "\\end{frame}" << endl;
}

/*
 * Plotter for QA related histograms
 *
 *  Implemented QA checks:
 *    - spillover correction
 *    - seagull correction
 *    - background level check
 *    - jet shape corrections
 */
void qaPlotter(){
  
  ///////////////////
  // Configuration //
  ///////////////////
  
  bool saveFigures = true;          // Save the figures to a file
  
  bool drawSpillover = true;              // Draw the QA plots for spillover correction
  bool drawSeagull = false;                // Draw the QA plots for seagull correction
  bool calculateBackgroundOverlap = false; // Check difference in background overlap region of leading and subleading jets
  bool drawJetShapeCorrections = false;    // Draw the jet shape corrections as a function or R
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = true;    // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = false;     // Produce the correction for inclusive jet-track correlations
  
  bool jetShapeCorrectionComparison = false; // Draw the comparison plots between JFF and spillover corrections
  bool jetShapeCorrectionBigCanvas = false;   // Draw JFF and spillover corrections in all centrality on pT bins to big canvas
  bool constantBigCanvasScale = false;        // Use same scale for all bins in big canvas
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,(regularJetTrack && !calculateBackgroundOverlap),uncorrectedJetTrack,(ptWeightedJetTrack && !calculateBackgroundOverlap),inclusiveJetTrack,inclusiveJetTrack};
  const char *titleAddition[DijetHistogramManager::knJetTrackCorrelations] = {"",", uncorrected",", $p_{\\mathrm{T}}$ weighted","",", uncorrected",", $p_{\\mathrm{T}}$ weighted","",", $p_{\\mathrm{T}}$ weighted"};
  
  /////////////////
  // Config done //
  /////////////////
  
  // Open files containing the QA histograms
  TFile *spilloverQaFile = TFile::Open("data/spilloverCorrection_PbPbMC_skims_pfJets_noInclOrUncorr_10eventsMixed_subeNon0_smoothedMixing_revisedFit_2019-02-18_QA.root");
  // "data/spilloverCorrection_PbPbMC_skims_pfJets_noInclOrUncorr_10eventsMixed_subeNon0_smoothedMixing_revisedFit_2019-02-18_QA.root"
  // "data/spilloverCorrection_PbPbMC_skims_pfJets_noInclusiveOrUncorrected_10eventsMixed_subeNon0_smoothedMixing_2019-02-14_QA.root"
  // "data/spilloverCorrection_PbPbMC_skims_pfJets_noInclusiveOrUncorrected_10eventsMixed_subeNon0_smoothedMixing_2018-11-27_QA.root"
  TFile *seagullFile = TFile::Open("data/dijetPbPb_skims_pfJets_noUncorr_smoothedMixingAvoidPeakLowPt_noCorrections_2019-02-06_QA.root");
  // "data/dijetPbPb_skims_pfJets_noUncorr_smoothedMixingAvoidPeakLowPt_noCorrections_2019-02-06_QA.root"
  // "data/dijetPbPb_skims_pfJets_noUncorr_improvedPoolMixing_noJetLimit_noCorrections_processed_2019-01-09_QA.root"
  TFile *seagullPpFile = TFile::Open("data/dijet_pp_highForest_pfJets_noUncorr_noJetLimit_noCorrections_processed_2019-01-14_QA.root");
  TFile *backgroundFile = TFile::Open("data/dijetPbPb_skims_pfJets_noUncorr_improvedPoolMixing_noJetLimit_noCorrections_processed_2019-01-09.root");
  // "data/dijet_pp_highForest_pfJets_noUncorr_noJetLimit_noCorrections_processed_2019-01-14.root"
  // "data/dijet_pp_highForest_pfJets_noUncorr_noJetLimit_noCorrections_adjustedBackground_processed_2019-01-14.root"
  // "data/dijetPbPb_skims_pfJets_noUncorr_improvedPoolMixing_noJetLimit_noCorrections_processed_2019-01-09.root"
  // "data/dijetPbPb_skims_pfJets_noUncorr_improvedPoolMixing_noJetLimit_noCorrectionsOrAdjust_processed_2019-01-09.root"
  TFile *jffFile = TFile::Open("data/jffCorrection_PbPbMC_noInclOrUncorr_10eveMixed_sube0_smoothedMixing_adjustedBackground_2018-11-27.root");
  // "data/jffCorrection_PbPbMC_noInclOrUncorr_10eveMixed_sube0_smoothedMixing_adjustedBackground_2018-11-27.root"
  TFile *spilloverFile = TFile::Open("data/spilloverCorrection_PbPbMC_skims_pfJets_noInclOrUncorr_10eventsMixed_subeNon0_smoothedMixing_revisedFit_2019-02-18.root");
  // "data/spilloverCorrection_PbPbMC_skims_pfJets_noInclOrUncorr_10eventsMixed_subeNon0_smoothedMixing_revisedFit_2019-02-18.root"
  // "data/spilloverCorrection_PbPbMC_skims_pfJets_noInclusiveOrUncorrected_10eventsMixed_subeNon0_smoothedMixing_2019-02-14.root"
  
  // Read the number of bins from histogram manager
  DijetHistogramManager *dummyManager = new DijetHistogramManager();
  const int nCentralityBins = dummyManager->GetNCentralityBins();
  const int nTrackPtBins = dummyManager->GetNTrackPtBins();
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,300};  // Bin borders for track pT
  
  // Centrality range to be considered (only implemented to background check at the moment)
  int lastCentralityBin = 4;
  
  // Define needed histograms
  TH1D *spilloverDeltaEtaProjection[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaPhiProjection[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH2D *spilloverDeltaEtaDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *seagullDeltaEtaWings[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins+1][nTrackPtBins];
  TH1D *backgroundDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *backgroundDeltaPhiOverlap[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *correctionJetShape[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins]; // 0 = JFF, 1 = spillover
  
  const char *histogramNames[2] = {"spilloverQA","dataReplica"};
  
  // Read the histograms from the QA file
  char histogramNamer[200];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      if(iCentrality > lastCentralityBin) continue; // No centrality selection for pp
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Seagull histograms for PbPb
        sprintf(histogramNamer,"seagullDeltaEta_%s/seagullDeltaEta_%s_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        seagullDeltaEtaWings[iJetTrack][iCentrality][iTrackPt] = (TH1D*) seagullFile->Get(histogramNamer);
        
        // Seagull histograms for pp
        if(iCentrality == 0){
          sprintf(histogramNamer,"seagullDeltaEta_%s/seagullDeltaEta_%s_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
          seagullDeltaEtaWings[iJetTrack][nCentralityBins][iTrackPt] = (TH1D*) seagullPpFile->Get(histogramNamer);
        }
        
        // Jet shape histograms for JFF correction
        if(drawJetShapeCorrections){
          sprintf(histogramNamer,"JetShape_%s/jffCorrection_JetShape_%s_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
          correctionJetShape[0][iJetTrack][iCentrality][iTrackPt] = (TH1D*) jffFile->Get(histogramNamer);
          
          sprintf(histogramNamer,"%sDeltaEtaDeltaPhi/spilloverCorrection_%sDeltaEtaDeltaPhi_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
          spilloverDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt] = (TH2D*) spilloverFile->Get(histogramNamer);
        }
        
        // Background histograms for background adjustment check
        sprintf(histogramNamer,"%s/%sDeltaPhi_Background_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        backgroundDeltaPhi[iJetTrack][iCentrality][iTrackPt] = (TH1D*) backgroundFile->Get(histogramNamer);
        
        sprintf(histogramNamer,"%s/%sDeltaPhi_BackgroundOverlap_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        backgroundDeltaPhiOverlap[iJetTrack][iCentrality][iTrackPt] = (TH1D*) backgroundFile->Get(histogramNamer);
        
        // Projections for spillover
        for(int iDataType = 0; iDataType < 2; iDataType++){
        
          sprintf(histogramNamer,"%sDeltaEtaProjection/%s_%sDeltaEtaProjection_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),histogramNames[iDataType],dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
          spilloverDeltaEtaProjection[iDataType][iJetTrack][iCentrality][iTrackPt] = (TH1D*) spilloverQaFile->Get(histogramNamer);
        
          sprintf(histogramNamer,"%sDeltaPhiProjection/%s_%sDeltaPhiProjection_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),histogramNames[iDataType],dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
          spilloverDeltaPhiProjection[iDataType][iJetTrack][iCentrality][iTrackPt] = (TH1D*) spilloverQaFile->Get(histogramNamer);
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
          deltaEtaCanvas[iJetTrack]->cd(nCentralityBins-1-iCentrality+nCentralityBins*iTrackPt+1);
          gPad->SetTopMargin(0.1);
          gPad->SetBottomMargin(0.2);
          
          // Draw the histogram to canvas
          setHistogramStyle(spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt],1.5,2,titleString,"#Delta#eta");
          spilloverDeltaEtaProjection[1][iJetTrack][iCentrality][iTrackPt]->SetLineColor(kRed);
          spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->DrawCopy();
          
          legend = new TLegend(0.61,0.6,0.9,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.08);legend->SetTextFont(62);
          legend->AddEntry(spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt],"Hydjet","l");
          legend->AddEntry(spilloverDeltaEtaProjection[1][iJetTrack][iCentrality][iTrackPt],"Fit","l");
          legend->Draw();
          
          // Change to comparison canvas
          deltaEtaComparisonCanvas[iJetTrack]->cd(nCentralityBins-1-iCentrality+nCentralityBins*iTrackPt+1);
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
          deltaPhiCanvas[iJetTrack]->cd(nCentralityBins-1-iCentrality+nCentralityBins*iTrackPt+1);
          gPad->SetTopMargin(0.1);
          gPad->SetBottomMargin(0.2);
          
          // Draw the histogram to canvas
          setHistogramStyle(spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt],1.5,3,titleString,"#Delta#phi");
          spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt]->SetLineColor(kRed);
          spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->DrawCopy();
          
          legend = new TLegend(0.61,0.6,0.9,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.08);legend->SetTextFont(62);
          legend->AddEntry(spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt],"Hydjet","l");
          legend->AddEntry(spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt],"Fit","l");
          legend->Draw();
          
          // Change to comparison canvas
          deltaPhiComparisonCanvas[iJetTrack]->cd(nCentralityBins-1-iCentrality+nCentralityBins*iTrackPt+1);
          gPad->SetTopMargin(0.1);
          gPad->SetBottomMargin(0.2);
          
          // Draw the histograms without fits to the canvas
          setHistogramStyle(spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt],1.5,-1,titleString,"#Delta#phi",true);
          setHistogramStyle(spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt],1.5,-1,titleString,"#Delta#phi",true);
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
        deltaPhiComparisonCanvas[iJetTrack]->SaveAs(histogramNamer);
        
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
      seagullCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1250,1800);
      seagullCanvas[iJetTrack]->Divide(nCentralityBins+1,nTrackPtBins);
      
      for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          if(iCentrality == nCentralityBins){
            titleString = Form("pp - Track pT: %.1f-%.1f GeV",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
          } else {
            titleString = Form("Cent: %.0f-%.0f%% - Track pT: %.1f-%.1f GeV",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1],trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
          }
          
          // Find the correct pad inside the canvas
          seagullCanvas[iJetTrack]->cd(nCentralityBins-iCentrality+(nCentralityBins+1)*iTrackPt+1);
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
  
  // ******************************************
  // **    Calculating background overlap    **
  // ******************************************
  
  // In the overlap histogram the filled bins are: 98-103
  if(calculateBackgroundOverlap){
    
    // Turn off subleading side, since it is same as leading but mirrored
    for(int i = 3; i < 6; i++){
      correlationSelector[i] = false;
    }
    
    // Define helper variables to calculate ratios and errors
    double backgroundValue, overlapValue, backgroundError, overlapError;
    double ratioValue[nTrackPtBins][6], ratioError[nTrackPtBins][6];
    double leadingSum[nCentralityBins][6], subleadingSum[nCentralityBins][6];
    double leadingSumError[nCentralityBins][6], subleadingSumError[nCentralityBins][6];
    double averageValue[nTrackPtBins], averageError[nTrackPtBins];
    double averageSubleadingOverlap[nTrackPtBins], averageSubleadingOverlapError[nTrackPtBins];
    double averageLeadingOverlap[nTrackPtBins], averageLeadingOverlapError[nTrackPtBins];
    
    // Initialize all the helper variables to zero
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      averageValue[iTrackPt] = 0;
      averageError[iTrackPt] = 0;
      averageSubleadingOverlap[iTrackPt] = 0;
      averageSubleadingOverlapError[iTrackPt] = 0;
      averageLeadingOverlap[iTrackPt] = 0;
      averageLeadingOverlapError[iTrackPt] = 0;
      for(int iBin = 0; iBin < 6; iBin ++){
        ratioValue[iTrackPt][iBin] = 0;
        ratioError[iTrackPt][iBin] = 0;
      }
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        leadingSum[iCentrality][iTrackPt] = 0;
        subleadingSum[iCentrality][iTrackPt] = 0;
        leadingSumError[iCentrality][iTrackPt] = 0;
        subleadingSumError[iCentrality][iTrackPt] = 0;
      }
    }
    
    // Loop over all the bins
    for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
      if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        titleString = Form("Cent: %.0f-%.0f\\%%%s",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1],titleAddition[iJetTrack]);
        
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          // Initialize averages to 0
          averageValue[iTrackPt] = 0;
          averageError[iTrackPt] = 0;
          averageLeadingOverlap[iTrackPt] = 0;
          averageSubleadingOverlap[iTrackPt] = 0;
          averageLeadingOverlapError[iTrackPt] = 0;
          averageSubleadingOverlapError[iTrackPt] = 0;
          leadingSum[iCentrality][iTrackPt] = 0;
          subleadingSum[iCentrality][iTrackPt] = 0;
          leadingSumError[iCentrality][iTrackPt] = 0;
          subleadingSumError[iCentrality][iTrackPt] = 0;
          
          // For pp, no centrality binning
          if(iCentrality > lastCentralityBin) continue;
          
          // Loop over the bins which have content in the background overlap histogram
          for(int iBin = 98; iBin <= 103; iBin++){
            backgroundValue = backgroundDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetBinContent(iBin);
            overlapValue = backgroundDeltaPhiOverlap[iJetTrack][iCentrality][iTrackPt]->GetBinContent(iBin);
            backgroundError = backgroundDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetBinError(iBin);
            overlapError = backgroundDeltaPhiOverlap[iJetTrack][iCentrality][iTrackPt]->GetBinError(iBin);
            
            // Sum together all leading and subleading
            if(iBin < 101){  // If background = leading and overlap = subleading
              leadingSum[iCentrality][iTrackPt] += backgroundValue;
              leadingSumError[iCentrality][iTrackPt] += backgroundError;
              subleadingSum[iCentrality][iTrackPt] += overlapValue;
              subleadingSumError[iCentrality][iTrackPt] += overlapError;
            } else { // If background = subleading and overlap = leading
              leadingSum[iCentrality][iTrackPt] += overlapValue;
              leadingSumError[iCentrality][iTrackPt] += overlapError;
              subleadingSum[iCentrality][iTrackPt] += backgroundValue;
              subleadingSumError[iCentrality][iTrackPt] += backgroundError;
            }
            
            // Calcaulate the ratio and the error for the ratio for overlapping points
            ratioValue[iTrackPt][iBin-98] = backgroundValue/overlapValue;
            ratioError[iTrackPt][iBin-98] = TMath::Sqrt(TMath::Power(backgroundError/overlapValue,2)+TMath::Power((backgroundValue*overlapError)/(overlapValue*overlapValue),2));
            
            // Calculate average of the ratios in both overlapping areas
            averageValue[iTrackPt] += ratioValue[iTrackPt][iBin-98]/6.0;
            averageError[iTrackPt] += TMath::Power(ratioError[iTrackPt][iBin-98]/6.0,2);
            if(iBin <= 100){
              averageSubleadingOverlap[iTrackPt] += ratioValue[iTrackPt][iBin-98]/3.0;
              averageSubleadingOverlapError[iTrackPt] += TMath::Power(ratioError[iTrackPt][iBin-98]/3.0,2);
            } else {
              averageLeadingOverlap[iTrackPt] += ratioValue[iTrackPt][iBin-98]/3.0;
              averageLeadingOverlapError[iTrackPt] += TMath::Power(ratioError[iTrackPt][iBin-98]/3.0,2);
            }
          } // Loop over bins that are filled in the overlap histogram
          averageError[iTrackPt] = TMath::Sqrt(averageError[iTrackPt]);
          averageLeadingOverlapError[iTrackPt] = TMath::Sqrt(averageLeadingOverlapError[iTrackPt]);
          averageSubleadingOverlapError[iTrackPt] = TMath::Sqrt(averageSubleadingOverlapError[iTrackPt]);
          
        } // Track pT loop
        
        // Print the obtained results to LaTeX slides
        printBackgroundOverlapSlide(true,ratioValue,ratioError,averageSubleadingOverlap,averageSubleadingOverlapError,titleString);
        
        printBackgroundOverlapSlide(false,ratioValue,ratioError,averageLeadingOverlap,averageLeadingOverlapError,titleString);
        
      } // Centrality loop
      
      printBackroundSumSlide(leadingSum,leadingSumError,subleadingSum,subleadingSumError);
      
    } // Jet-track type loop
  } // Background overlap numbers
  
  // ***********************************************
  // **   Drawing JFF correction jet shape plots  **
  // ***********************************************
  
  if(drawJetShapeCorrections){
    
    // Create canvases for different qa plots
    TCanvas *jffCorrectionCanvas[DijetHistogramManager::knJetTrackCorrelations];
    TCanvas *spilloverCorrectionCanvas[DijetHistogramManager::knJetTrackCorrelations];
    
    DijetMethods *jetShapeCalculator = new DijetMethods();
    JDrawer *drawer = new JDrawer();
    double max1, max2, min1, min2, theMax, theMin;
    double legendX1, legendY1;
    
    for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
      if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
      
      // Create one big canvas with a pad for each centrality and track pT bin
      sprintf(histogramNamer,"jffCorrection%d",iJetTrack);
      sprintf(padNamer,"JFF correction %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
      jffCorrectionCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1000,1800);
      jffCorrectionCanvas[iJetTrack]->Divide(nCentralityBins,nTrackPtBins);
      
      sprintf(histogramNamer,"spilloverCorrection%d",iJetTrack);
      sprintf(padNamer,"Spillover correction %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
      spilloverCorrectionCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1000,1800);
      spilloverCorrectionCanvas[iJetTrack]->Divide(nCentralityBins,nTrackPtBins);
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        if(iCentrality > lastCentralityBin) continue; // No centrality selection for pp
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          // Calculate jet shape from the two-dimensional spillover distribution
          correctionJetShape[1][iJetTrack][iCentrality][iTrackPt] = jetShapeCalculator->GetJetShape(spilloverDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]);
          
          /////////////////////////////////////////////
          // First, fill the plots with big canvases
          /////////////////////////////////////////////
          
          if(jetShapeCorrectionBigCanvas){
            
            // Determine title for pads
            if(iCentrality == nCentralityBins){
              titleString = Form("pp - Track pT: %.1f-%.1f GeV",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
            } else {
              titleString = Form("Cent: %.0f-%.0f%% - Track pT: %.1f-%.1f GeV",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1],trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
            }
            
            // Change to canvas for JFF correction
            jffCorrectionCanvas[iJetTrack]->cd(nCentralityBins-1-iCentrality+nCentralityBins*iTrackPt+1);
            gPad->SetTopMargin(0.1);
            gPad->SetBottomMargin(0.2);
            
            // Draw the JFF correction histogram to the canvas
            gStyle->SetTitleSize(0.09,"t");
            setHistogramTitles(correctionJetShape[0][iJetTrack][iCentrality][iTrackPt],titleString,"#Deltar");
            if(constantBigCanvasScale){
              if(iJetTrack == DijetHistogramManager::kTrackSubleadingJet){
                correctionJetShape[0][iJetTrack][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(-0.9,0.4);
              } else if(iJetTrack == DijetHistogramManager::kPtWeightedTrackLeadingJet){
                correctionJetShape[0][iJetTrack][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(-2.5,0.7);
              } else if(iJetTrack == DijetHistogramManager::kPtWeightedTrackSubleadingJet){
                correctionJetShape[0][iJetTrack][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(-1.2,0.6);
              } else {
                correctionJetShape[0][iJetTrack][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(-0.9,0.9);
              }
            }
            correctionJetShape[0][iJetTrack][iCentrality][iTrackPt]->DrawCopy();
            
            // Change to canvas for spillover correction
            spilloverCorrectionCanvas[iJetTrack]->cd(nCentralityBins-1-iCentrality+nCentralityBins*iTrackPt+1);
            gPad->SetTopMargin(0.1);
            gPad->SetBottomMargin(0.2);
            
            // Draw the JFF correction histogram to the canvas
            gStyle->SetTitleSize(0.09,"t");
            setHistogramTitles(correctionJetShape[1][iJetTrack][iCentrality][iTrackPt],titleString,"#Deltar");
            if(constantBigCanvasScale){
              if(iJetTrack == DijetHistogramManager::kPtWeightedTrackLeadingJet){
                correctionJetShape[1][iJetTrack][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0,3.4);
              } else {
                correctionJetShape[1][iJetTrack][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0,2.4);
              }
            }
            correctionJetShape[1][iJetTrack][iCentrality][iTrackPt]->DrawCopy();
          
          } // Big canvas if
          
          //////////////////////////////////////////////////
          // Next, draw comparison plots in small canvases
          //////////////////////////////////////////////////
          
          if(jetShapeCorrectionComparison){
            
            // Gether information of the current bin
            if(lastCentralityBin == 0){
              titleString = Form("%s - pp",dummyManager->GetJetTrackAxisName(iJetTrack));
            } else {
              titleString = Form("%s - Cent: %.0f-%.0f%%",dummyManager->GetJetTrackAxisName(iJetTrack),centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
            }
            
            sprintf(padNamer,"Track pT: %.1f-%.1f GeV",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
            
            // Find a good drawing range for the distributions
            std::tie(min1,max1) = getDrawMinMax(correctionJetShape[0][iJetTrack][iCentrality][iTrackPt]);
            std::tie(min2,max2) = getDrawMinMax(correctionJetShape[1][iJetTrack][iCentrality][iTrackPt]);
            theMin = (min1 < min2) ? min1 : min2;
            theMax = (max1 > max2) ? max1 : max2;
            
            // Draw the histogram to canvas together with the one from file
            correctionJetShape[0][iJetTrack][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(theMin,theMax);
            drawer->DrawHistogram(correctionJetShape[0][iJetTrack][iCentrality][iTrackPt], "#Deltar","P(#Deltar)", " ");
            correctionJetShape[1][iJetTrack][iCentrality][iTrackPt]->SetLineColor(kRed);
            correctionJetShape[1][iJetTrack][iCentrality][iTrackPt]->Draw("same");
            
            // Draw a legend to the plot
            legendX1 = 0.47; legendY1 = 0.61;
            
            // Very confusing way to put the legend in a good position :p. Good luck changing if needed.
            if((iJetTrack == 0 || iJetTrack == 2) && iTrackPt == 4 && iCentrality < 2) legendY1 = 0.18;
            if(((iJetTrack == 0 || iJetTrack == 2) && iTrackPt < 2 && iCentrality > 1) && !(iTrackPt == 1 && iCentrality == 2)) legendY1 = 0.18;
            if((iJetTrack == 3 || iJetTrack == 5) && !(iTrackPt == 1 && iCentrality < 2)) legendY1 = 0.18;
            if((iJetTrack == 3 || iJetTrack == 5) && ((iCentrality > 1 && iTrackPt == 5) || (iCentrality == 2 && iTrackPt == 0))) legendY1 = 0.61;
            
            legend = new TLegend(legendX1,legendY1,legendX1+0.43,legendY1+0.3);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->SetHeader(titleString.Data());
            legend->AddEntry((TObject*) 0,padNamer,"");
            legend->AddEntry(correctionJetShape[0][iJetTrack][iCentrality][iTrackPt],"JFF correction","l");
            legend->AddEntry(correctionJetShape[1][iJetTrack][iCentrality][iTrackPt],"Spillover correction","l");
            legend->Draw();
            
            // Save the figure to file
            if(saveFigures) {
              gPad->GetCanvas()->SaveAs(Form("figures/jetShapeCorrectionComparison_%s_C%dT%d.pdf",dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt));
            }
          } // Correction comparison if
          
        } // Track pT loop
      } // Centrality loop
      
      // After all the canvases are filled, save them
      if(saveFigures && jetShapeCorrectionBigCanvas){
        
        // Save the big JFF correction canvas
        if(constantBigCanvasScale){
          sprintf(histogramNamer,"figures/%s_jffCorrectionQA_constantScale.pdf",dummyManager->GetJetTrackHistogramName(iJetTrack));
        } else {
          sprintf(histogramNamer,"figures/%s_jffCorrectionQA.pdf",dummyManager->GetJetTrackHistogramName(iJetTrack));
        }
        
        jffCorrectionCanvas[iJetTrack]->SaveAs(histogramNamer);
        
        // Save the big spillover correction canvas
        if(constantBigCanvasScale){
          sprintf(histogramNamer,"figures/%s_spilloverCorrectionQA_constantScale.pdf",dummyManager->GetJetTrackHistogramName(iJetTrack));
        } else {
          sprintf(histogramNamer,"figures/%s_spilloverCorrectionQA.pdf",dummyManager->GetJetTrackHistogramName(iJetTrack));
        }
        
        spilloverCorrectionCanvas[iJetTrack]->SaveAs(histogramNamer);
        
      } // saving figures
      
    } // Jet-track category loop
    
  } // Jet shape plots from JFF correction
}
