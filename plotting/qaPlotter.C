#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"
#include "JDrawer.h"
#include "DijetCard.h"
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
 * Draw a legend to a big canvas plot
 */
void drawBigCanvasLegend(TLegend* legend, TH1* histogram, const char* legendString){
  TH1D *legender = (TH1D*) histogram->Clone(Form("legender%s",histogram->GetName()));
  legend = new TLegend(0.15,0.3,0.9,0.5);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.11);legend->SetTextFont(62);
  legend->AddEntry(legender,legendString,"p");
  legend->Draw();
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
 * Draw the graphs for spillover correction
 *
 *  JDrawer *drawer = JDrawer preparing the canvas for drawing
 *  TGraphErrors *spilloverGraph[4] = Array of graphs in centrality bins
 *  int iJetTrack = Jet-track correlation type index
 *  double yZoom = Zoom for y-axis
 *  const char* yTitle = Title for the y-axis
 *  bool logY = Flag for logarithmic y-axis
 *  bool saveFigures = Flag for saving the figures to a file
 *  const char* saveName = Name given to the figure saved to the file
 *  TGraphErrors *averageGraph = Average of the other graphs
 *  const char* fitFunction = Name of the fit function
 */
void drawSpilloverGraph(JDrawer *drawer, TGraphErrors *spilloverGraph[4], int iJetTrack, double yZoom, const char* yTitle, bool logY, bool saveFigures, const char* saveName, TGraphErrors *averageGraph, const char* fitFunction){
  
  // Create a histogram manager for naming purposes
  DijetHistogramManager *dummyManager = new DijetHistogramManager();
  char namer[70];
  
  // Binning information for centrality
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  const int nCentralityBins = 4;
  
  // Setup a legend for the delta eta width
  TLegend *legend = new TLegend(0.5,0.5,0.9,0.9);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  legend->SetHeader(dummyManager->GetJetTrackAxisName(iJetTrack));
  
  // Graph styling information
  int graphMarkerStyle[] = {21,20,34,47};
  int graphMarkerColor[] = {kBlack,kRed,kBlue,kGreen+3};
  
  // Helper variables
  TString titleString;
  TF1* thisFit;
  
  // Draw the for graphs for delta eta width
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    titleString = Form("Cent: %.0f-%.0f%%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    
    // Set a nice style for the graph
    spilloverGraph[iCentrality]->SetMarkerStyle(graphMarkerStyle[iCentrality]);
    spilloverGraph[iCentrality]->SetMarkerColor(graphMarkerColor[iCentrality]);
    
    // Draw the fits to each centrality bin only if average graph is not provided
    thisFit = spilloverGraph[iCentrality]->GetFunction(fitFunction);
    if(thisFit){
      if(averageGraph == NULL){
        thisFit->SetLineColor(graphMarkerColor[iCentrality]);
        thisFit->SetRange(0.5,8);
      } else {
        thisFit->SetLineWidth(0);
      }
    }
    
    drawer->SetLogY(logY);
    // Draw the graph with a new canvas for the first centrality bin and the other centralities to the same figure
    if(iCentrality == 0){
      // For some weird reason the fit line is not fit if this histogram is not cloned and the fit is disabled later
      // for the drawing of averaged fit.
      sprintf(namer,"graphFit%s%s",saveName,dummyManager->GetJetTrackHistogramName(iJetTrack));
      drawer->DrawGraph((TGraphErrors*)spilloverGraph[iCentrality]->Clone(namer),0,12,0.00001,yZoom,"Track p_{T} (GeV)",yTitle,"","psame");
    } else {
      spilloverGraph[iCentrality]->Draw("psame");
    }
    
    // Add the histogram to legend
    legend->AddEntry(spilloverGraph[iCentrality],titleString,"p");
    
  }
  
  // Draw the average graph on top of other points
  if(averageGraph != NULL){
    
    // Set a nice style for the graph
    averageGraph->SetMarkerStyle(kFullStar);
    averageGraph->SetMarkerColor(kMagenta);
    averageGraph->SetMarkerSize(1.5);
    
    // Change the function line color to match the points
    thisFit = (TF1*) averageGraph->GetListOfFunctions()->First();
    if(thisFit) thisFit->SetLineColor(kMagenta);
    
    // Draw the average graph to the same canvas
    averageGraph->Draw("psame");
    
    // Add the graph to the legend
    legend->AddEntry(averageGraph,"Average","p");
  }
  
  // Draw the legend to the figure
  legend->Draw();
  
  // Save the figure to a file
  if(saveFigures) {
    gPad->GetCanvas()->SaveAs(Form("figures/%s_%s.pdf",saveName,dummyManager->GetJetTrackHistogramName(iJetTrack)));
  }
 
  // Delete the dummy histogram manager
  delete dummyManager;
}

/*
 * Draw the graphs for spillover correction
 *
 *  JDrawer *drawer = JDrawer preparing the canvas for drawing
 *  TGraphErrors *etaGraph = DeltaEta graph for comparison
 *  TGraphErrors *phiGraph = DeltaPhi graph for comparison
 *  int iJetTrack = Jet-track correlation type index
 *  double yZoom = Zoom for y-axis
 *  const char* yTitle = Title for the y-axis
 *  bool logY = Flag for logarithmic y-axis
 *  bool saveFigures = Flag for saving the figures to a file
 *  const char* saveName = Name given to the figure saved to the file
 *  TGraphErrors *combinedGraph = Graph combining deltaEta and deltaPhi
 *  TGraphErrors *binCountGraph = Graph obtained from bin counting in deltaEta
 *  const char* fitFunction = Name of the function fitted to the graphs
 */
void drawSpilloverGraphComparison(JDrawer *drawer, TGraphErrors *etaGraph, TGraphErrors *phiGraph, int iJetTrack, int iCentrality, double yZoom, const char* yTitle, bool logY, bool saveFigures, const char* saveName, TGraphErrors *combinedGraph = NULL, TGraphErrors *binCountGraph = NULL, const char* fitFunction = "expo"){
  
  // Create a histogram manager for naming purposes
  DijetHistogramManager *dummyManager = new DijetHistogramManager();
  TF1* thisFit;
  char namer[50];
  
  // Binning information for centrality
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  
  // Setup a legend for the delta eta width
  TLegend *legend = new TLegend(0.5,0.65,0.9,0.9);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  legend->SetHeader(Form("%s, %.0f-%.0f%%",dummyManager->GetJetTrackAxisName(iJetTrack),centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
  
  // Graph styling information
  int graphMarkerStyle[] = {21,20,34,47};
  int graphMarkerColor[] = {kBlack,kRed,kBlue,kGreen+3};
  
  // Set the style and draw the first graph
  drawer->SetLogY(logY);
  etaGraph->SetMarkerStyle(graphMarkerStyle[0]);
  etaGraph->SetMarkerColor(graphMarkerColor[0]);
  drawer->DrawGraph(etaGraph,0,12,0.001,yZoom,"Track p_{T} (GeV)",yTitle,"","psame");
  legend->AddEntry(etaGraph,"#Delta#eta fit","p");
  
  // Do not draw the fit to deltaEta if combined histogram is provided
  thisFit = etaGraph->GetFunction(fitFunction);
  if(thisFit){
    if(combinedGraph != NULL){
      thisFit->SetLineWidth(0);
    } else {
      thisFit->SetLineColor(graphMarkerColor[0]);
      thisFit->SetRange(0.5,8);
    }
  }
  
  // Set the style and draw the first graph
  phiGraph->SetMarkerStyle(graphMarkerStyle[1]);
  phiGraph->SetMarkerColor(graphMarkerColor[1]);
  phiGraph->Draw("psame");
  legend->AddEntry(phiGraph,"#Delta#phi fit","p");
  
  // Do not draw the fit to deltaPhi if combined histogram is provided
  thisFit = phiGraph->GetFunction(fitFunction);
  if(thisFit){
    if(combinedGraph != NULL){
      thisFit->SetLineWidth(0);
    } else {
      thisFit->SetLineColor(graphMarkerColor[1]);
      thisFit->SetRange(0.5,8);
    }
  }
  
  // Draw the combined graph if it is provided
  if(combinedGraph != NULL){
    combinedGraph->SetMarkerStyle(kFullStar);
    combinedGraph->SetMarkerColor(kMagenta);
    combinedGraph->SetMarkerSize(1.5);
    combinedGraph->Draw("psame");
    legend->AddEntry(combinedGraph,"Average","p");
    
    // Change the function line color to match the color of points
    thisFit = (TF1*) combinedGraph->GetListOfFunctions()->First();
    if(thisFit){
      thisFit->SetLineColor(kMagenta);
      thisFit->SetRange(0.5,8);
    }
  }
  
  // Draw the bin count brapg if it is provided
  if(binCountGraph != NULL){
    binCountGraph->SetMarkerStyle(graphMarkerStyle[2]);
    binCountGraph->SetMarkerColor(graphMarkerColor[2]);
    binCountGraph->Draw("psame");
    legend->AddEntry(binCountGraph,"Bin count","p");
    
    // Change the function line color to match the color of points
    thisFit = (TF1*) binCountGraph->GetListOfFunctions()->First();
    if(thisFit){
      thisFit->SetLineColor(graphMarkerColor[2]);
      thisFit->SetRange(0.5,8);
    }
  }
  
  // Draw the legend to the figure
  legend->Draw();
  
  // Save the figure to a file
  if(saveFigures) {
    gPad->GetCanvas()->SaveAs(Form("figures/%s_%s_C=%.0f-%.0f.pdf",saveName,dummyManager->GetJetTrackHistogramName(iJetTrack),centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
  }
  
  // Delete the dummy histogram manager
  delete dummyManager;
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
  
  bool saveFigures = false;          // Save the figures to a file
  
  bool drawSpillover = true;              // Draw the QA plots for spillover correction
  bool drawSeagull = false;                // Draw the QA plots for seagull correction
  bool calculateBackgroundOverlap = false; // Check difference in background overlap region of leading and subleading jets
  bool drawJetShapeCorrections = false;    // Draw the jet shape corrections as a function or R
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = false;    // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = false;     // Produce the correction for inclusive jet-track correlations
  
  bool jetShapeCorrectionComparison = false; // Draw the comparison plots between JFF and spillover corrections
  bool jetShapeCorrectionBigCanvas = true;   // Draw JFF and spillover corrections in all centrality on pT bins to big canvas
  bool constantBigCanvasScale = true;        // Use same scale for all bins in big canvas
  bool extraSpilloverComparison = false;      // Add constant fit deltaEta to spillover fit parameter comparison plots
  bool constantSpilloverScale = false;        // True = Draw spillover plots of one type all in the same scale, False = Zoom to fits
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,(regularJetTrack && !calculateBackgroundOverlap),uncorrectedJetTrack,(ptWeightedJetTrack && !calculateBackgroundOverlap),inclusiveJetTrack,inclusiveJetTrack};
  const char *titleAddition[DijetHistogramManager::knJetTrackCorrelations] = {"",", uncorrected",", $p_{\\mathrm{T}}$ weighted","",", uncorrected",", $p_{\\mathrm{T}}$ weighted","",", $p_{\\mathrm{T}}$ weighted"};
  
  /////////////////
  // Config done //
  /////////////////
  
  // Open files containing the QA histograms
  TFile *spilloverQaFile;
  spilloverQaFile = TFile::Open("newSpilloverTest_QA.root");
  // spillingOverTesting_QA.root
  // "data/spilloverCorrection_PbPbMC_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_smoothedMixing_refitParameters_2019-03-18_QA.root"
  TFile *seagullFile = TFile::Open("data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_subeNon0_xj_2019-06-06_onlySeagull_processed_QA.root");
  // "data/dijetPbPb_skims_pfJets_noUncorr_smoothedMixingAvoidPeakLowPt_noCorrections_2019-02-06_QA.root"
  // "data/dijetPbPb_skims_pfJets_noUncorr_improvedPoolMixing_noJetLimit_noCorrections_processed_2019-01-09_QA.root"
  // "data/PbPbMC_RecoGen_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_fixedCentality_processed_2019-02-15_QA.root"
  TFile *seagullPpFile = TFile::Open("data/dijet_ppMC_RecoGen_mergedSkims_Pythia6_pfJets_fixedJetPt_matchedJets_adjustedBackground_processed_2019-02-25_QA.root");
  // data/dijet_pp_highForest_pfJets_noUncorr_noJetLimit_noCorrections_processed_2019-01-14_QA.root
  TFile *backgroundFile = TFile::Open("data/dijetPbPb_skims_pfJets_noUncorr_improvedPoolMixing_noJetLimit_noCorrections_processed_2019-01-09.root");
  // "data/dijet_pp_highForest_pfJets_noUncorr_noJetLimit_noCorrections_processed_2019-01-14.root"
  // "data/dijet_pp_highForest_pfJets_noUncorr_noJetLimit_noCorrections_adjustedBackground_processed_2019-01-14.root"
  // "data/dijetPbPb_skims_pfJets_noUncorr_improvedPoolMixing_noJetLimit_noCorrections_processed_2019-01-09.root"
  // "data/dijetPbPb_skims_pfJets_noUncorr_improvedPoolMixing_noJetLimit_noCorrectionsOrAdjust_processed_2019-01-09.root"
  TFile *jffPbPbFile = TFile::Open("newPbPbTest.root");
  // "data/jffCorrection_PbPbMC_noInclOrUncorr_10eveMixed_sube0_smoothedMixing_adjustedBackground_2018-11-27.root"
  TFile *jffPpFile = TFile::Open("newPpTest2.root");
  // "data/jffCorrection_ppMC_mergedSkims_Pythia6_pfJets_noJetLimit_smoothedMixing_adjustedBackground_2019-01-15.root"
  TFile *spilloverFile = TFile::Open("newSpilloverTest.root");
  // data/spilloverCorrection_PbPbMC_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_smoothedMixing_refitParameters_2019-03-18.root
  
  // Read the number of bins from histogram manager
  DijetHistogramManager *dummyManager = new DijetHistogramManager();
  const int nCentralityBins = 4; // TODO: See if these can be read in a clever way
  const int nTrackPtBins = 6;    // TODO: Different files may have different binning!
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,300};  // Bin borders for track pT
  const char *saveNameComment = "";
  
  // Creata a dijet card for the JFF correction and setup asymmetry bins
  DijetCard* jffCard = new DijetCard(jffPbPbFile);
  const int nJffAsymmetryBins = jffCard->GetNAsymmetryBins();
  TString jffAsymmetryName[nJffAsymmetryBins+1];
  TString jffAsymmetryLegend[nJffAsymmetryBins+1];
  for(int iAsymmetry = 0; iAsymmetry < nJffAsymmetryBins; iAsymmetry++) {
    jffAsymmetryName[iAsymmetry] = Form("A%d",iAsymmetry);
    jffAsymmetryLegend[iAsymmetry] = Form("%.2f < %s < %.2f",jffCard->GetLowBinBorderAsymmetry(iAsymmetry),jffCard->GetAsymmetryBinType("latex"),jffCard->GetHighBinBorderAsymmetry(iAsymmetry));
  }
  jffAsymmetryName[nJffAsymmetryBins] = "";
  jffAsymmetryLegend[nJffAsymmetryBins] = Form("0.0 < %s < 1.0",jffCard->GetAsymmetryBinType("latex"));
  
  // Create a new DijetMethods so that several tricks can be easily made to distribution
  DijetMethods *methods = new DijetMethods();
  
  // Centrality range to be considered (only implemented to background check at the moment)
  int lastCentralityBin = nCentralityBins-1;
  
  // Define needed histograms
  TH1D *spilloverDeltaEtaProjection[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaPhiProjection[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH2D *spilloverDeltaEtaDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *seagullDeltaEtaWings[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins+1][nTrackPtBins];
  TH1D *backgroundDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *backgroundDeltaPhiOverlap[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *correctionJetShape[2][DijetHistogramManager::knJetTrackCorrelations][nJffAsymmetryBins+1][nCentralityBins+1][nTrackPtBins]; // 0 = JFF, 1 = spillover
  TF1 *spilloverDeltaEtaFit[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TF1 *spilloverDeltaPhiFit[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *trackPtForGraphs[nCentralityBins];  // Track pT histograms to find proper place to put pT point in graphs
  TH1D *deltaEtaProjectionFromCorrection;
  TH1D *deltaPhiProjectionFromCorrection;
  
  // Define graphs to illustrate the spillover calculation
  TGraphErrors *combinedGraphYield[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins];
  TGraphErrors *binCountGraphYield[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins];
  TGraphErrors *combinedGraphDeltaEtaWidth[DijetHistogramManager::knJetTrackCorrelations];
  TGraphErrors *combinedGraphDeltaPhiWidth[DijetHistogramManager::knJetTrackCorrelations];
  
  const char *histogramNames[2] = {"spilloverQA","dataReplica"};
  
  // Read the histograms from the QA file
  char histogramNamer[200];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    
    // Read the combined width graphs for spillover
    sprintf(histogramNamer,"%sCombinedGraph/combinedGraph_%sDeltaEta",dummyManager->GetJetTrackHistogramName(iJetTrack),dummyManager->GetJetTrackHistogramName(iJetTrack));
    combinedGraphDeltaEtaWidth[iJetTrack] = (TGraphErrors*) spilloverQaFile->Get(histogramNamer);
    
    sprintf(histogramNamer,"%sCombinedGraph/combinedGraph_%sDeltaPhi",dummyManager->GetJetTrackHistogramName(iJetTrack),dummyManager->GetJetTrackHistogramName(iJetTrack));
    combinedGraphDeltaPhiWidth[iJetTrack] = (TGraphErrors*) spilloverQaFile->Get(histogramNamer);
    
    for(int iCentrality = 0; iCentrality <= lastCentralityBin; iCentrality++){
      
      // Read the combined yield graphs for spillover
      sprintf(histogramNamer,"%sCombinedGraph/combinedGraph_%sYield_C%d",dummyManager->GetJetTrackHistogramName(iJetTrack),dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality);
      combinedGraphYield[iJetTrack][iCentrality] = (TGraphErrors*) spilloverQaFile->Get(histogramNamer);
      
      // Read the bin count yield for the spillover
      sprintf(histogramNamer,"%sCombinedGraph/binCountGraph_%sYield_C%d",dummyManager->GetJetTrackHistogramName(iJetTrack),dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality);
      binCountGraphYield[iJetTrack][iCentrality] = (TGraphErrors*) spilloverQaFile->Get(histogramNamer);
      
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
          
          // Asymmetry bins for the JFF correction
          for(int iAsymmetry = 0; iAsymmetry <= nJffAsymmetryBins; iAsymmetry++){
            
            sprintf(histogramNamer,"JetShape_%s/jffCorrection_JetShape_%s_%sC%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),dummyManager->GetJetTrackHistogramName(iJetTrack),jffAsymmetryName[iAsymmetry].Data(),iCentrality,iTrackPt);
            correctionJetShape[0][iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH1D*) jffPbPbFile->Get(histogramNamer);
            
            // Jet shape histograms for pp JFF correction
            if(iCentrality == 0){
              sprintf(histogramNamer,"JetShape_%s/jffCorrection_JetShape_%s_%sC%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),dummyManager->GetJetTrackHistogramName(iJetTrack),jffAsymmetryName[iAsymmetry].Data(),iCentrality,iTrackPt);
              correctionJetShape[0][iJetTrack][iAsymmetry][nCentralityBins][iTrackPt] = (TH1D*) jffPpFile->Get(histogramNamer);
            }
          }
        }
        
        sprintf(histogramNamer,"%sDeltaEtaDeltaPhi/fittedSpilloverCorrection_%sDeltaEtaDeltaPhi_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt] = (TH2D*) spilloverFile->Get(histogramNamer);
        
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
        
        // Fitted functions for spillover
        sprintf(histogramNamer,"%sDeltaEtaFit/%s_%sDeltaEtaFit_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),histogramNames[0],dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaEtaFit[iJetTrack][iCentrality][iTrackPt] = (TF1*) spilloverQaFile->Get(histogramNamer);
        
        sprintf(histogramNamer,"%sDeltaPhiFit/%s_%sDeltaPhiFit_C%dT%d",dummyManager->GetJetTrackHistogramName(iJetTrack),histogramNames[0],dummyManager->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaPhiFit[iJetTrack][iCentrality][iTrackPt] = (TF1*) spilloverQaFile->Get(histogramNamer);
      } // track pT
    } // centrality
  } // jet track
  
  // Read track pT histograms to find good places to put points in track pT axis for graphs
  for(int iCentrality = 0; iCentrality <= lastCentralityBin; iCentrality++){
    
    // Track pT histograms
    sprintf(histogramNamer,"track/trackPt_SameEvent_C%d",iCentrality);
    trackPtForGraphs[iCentrality] = (TH1D*) backgroundFile->Get(histogramNamer);
  }
  
  
  // ******************************************
  // **      Drawing spillover plots         **
  // ******************************************
  
  // Variables that can be reused for all QA plots
  TPad *currentPad;
  TLegend *legend;
  char padNamer[100];
  TString titleString;
  JDrawer *drawer = new JDrawer();
  double spilloverEtaIntegral, spilloverPhiIntegral;
  double spilloverScaleEta = -1;
  double spilloverScalePhi = -1;
  int deltaEtaRebinSpillover = 4;
  int deltaPhiRebinSpillover = 2;
  
  if(drawSpillover){
    
    // Create canvases for different QA plots
    TCanvas *deltaEtaCanvas[DijetHistogramManager::knJetTrackCorrelations];
    TCanvas *deltaEtaComparisonCanvas[DijetHistogramManager::knJetTrackCorrelations];
    TCanvas *deltaEtaValidationCanvas[DijetHistogramManager::knJetTrackCorrelations];
    TCanvas *deltaPhiCanvas[DijetHistogramManager::knJetTrackCorrelations];
    TCanvas *deltaPhiComparisonCanvas[DijetHistogramManager::knJetTrackCorrelations];
    TCanvas *deltaPhiValidationCanvas[DijetHistogramManager::knJetTrackCorrelations];
    
    // Create arrays to collect the information from the fits
    double deltaEtaFitWidth[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins]; // First bin: 0 = data 1 =error
    double deltaEtaFitYield[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins]; // First bin: 0 = data 1 =error
    double deltaPhiFitWidth[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins]; // First bin: 0 = data 1 =error
    double deltaPhiFitYield[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins]; // First bin: 0 = data 1 =error
    
    // Create graphs to draw the fit parameters
    TGraphErrors *graphDeltaEtaWidth[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins];
    TGraphErrors *graphDeltaEtaYield[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins];
    TGraphErrors *graphDeltaPhiWidth[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins];
    TGraphErrors *graphDeltaPhiYield[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins];
    
    // x-axis information for the graphs
    double graphPointsX[] = {0.85,1.5,2.5,3.5,6,10};    // x-axis points in flow graphs
    double graphErrorsX[] = {0,0,0,0,0,0};              // No errors for x-axis
    double yZoom = 2;
    int lowPtBin, highPtBin;
    double lowFitPt;
    char saveName[100];
    
    // Scales for different types of spillover plots
    double defaultSpilloverScalesEta[] = {2,2,2.8,1.4,1.4,2,2,2.5};
    double defaultSpilloverScalesPhi[] = {3,3,3.5,1.4,1.4,2,2.5,3.5};
    
    // Set nice drawing style for graphs
    drawer->SetDefaultAppearanceGraph();
    
    // Set a good title size for big canvases
    gStyle->SetTitleSize(0.09,"t");
    
    // Draw spillover corrections with the fit
    for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
      if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
      
      if(constantSpilloverScale){
        spilloverScaleEta = defaultSpilloverScalesEta[iJetTrack];
        spilloverScalePhi = defaultSpilloverScalesPhi[iJetTrack];
      }
      
      // Create one big canvas with a pad for each centrality and track pT bin
      sprintf(histogramNamer,"spilloverDeltaEta%d",iJetTrack);
      sprintf(padNamer,"Spillover deltaEta %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
      deltaEtaCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1000,1800);
      deltaEtaCanvas[iJetTrack]->Divide(nCentralityBins,nTrackPtBins);
      
      sprintf(histogramNamer,"comparisonDeltaEta%d",iJetTrack);
      sprintf(padNamer,"Comparison deltaEta %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
      deltaEtaComparisonCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1000,1800);
      deltaEtaComparisonCanvas[iJetTrack]->Divide(nCentralityBins,nTrackPtBins);
      
      sprintf(histogramNamer,"validationDeltaEta%d",iJetTrack);
      sprintf(padNamer,"Validation deltaEta %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
      deltaEtaValidationCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1000,1800);
      deltaEtaValidationCanvas[iJetTrack]->Divide(nCentralityBins,nTrackPtBins);
      
      sprintf(histogramNamer,"spilloverDeltaPhi%d",iJetTrack);
      sprintf(padNamer,"Spillover deltaPhi %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
      deltaPhiCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1000,1800);
      deltaPhiCanvas[iJetTrack]->Divide(nCentralityBins,nTrackPtBins);
      
      sprintf(histogramNamer,"comparisonDeltaPhi%d",iJetTrack);
      sprintf(padNamer,"Comparison deltaPhi %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
      deltaPhiComparisonCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1000,1800);
      deltaPhiComparisonCanvas[iJetTrack]->Divide(nCentralityBins,nTrackPtBins);
      
      sprintf(histogramNamer,"validationDeltaPhi%d",iJetTrack);
      sprintf(padNamer,"Validation deltaPhi %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
      deltaPhiValidationCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1000,1800);
      deltaPhiValidationCanvas[iJetTrack]->Divide(nCentralityBins,nTrackPtBins);
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          titleString = Form("Cent: %.0f-%.0f%% - Track pT: %.1f-%.1f GeV",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1],trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
          
          // Find the correct pad inside the canvas
          deltaEtaCanvas[iJetTrack]->cd(nCentralityBins-1-iCentrality+nCentralityBins*iTrackPt+1);
          gPad->SetTopMargin(0.1);
          gPad->SetBottomMargin(0.2);
          
          // Rebin the deltaEta histograms
          spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->Rebin(deltaEtaRebinSpillover);
          spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->Scale(1.0/deltaEtaRebinSpillover);
          spilloverDeltaEtaProjection[1][iJetTrack][iCentrality][iTrackPt]->Rebin(deltaEtaRebinSpillover);
          spilloverDeltaEtaProjection[1][iJetTrack][iCentrality][iTrackPt]->Scale(1.0/deltaEtaRebinSpillover);
          
          // Draw the histogram to canvas
          setHistogramStyle(spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt],2.5,spilloverScaleEta,titleString,"#Delta#eta");
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
          spilloverEtaIntegral = spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->Integral(spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->FindBin(-1.5+0.001),spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->FindBin(1.5-0.001),"width");
          
          legend = new TLegend(0.61,0.6,0.9,0.9);
          setupRatioLegend(legend,spilloverDeltaEtaProjection[1][iJetTrack][iCentrality][iTrackPt],spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]);
          legend->Draw();
          
          // Change to the validation canvas
          deltaEtaValidationCanvas[iJetTrack]->cd(nCentralityBins-1-iCentrality+nCentralityBins*iTrackPt+1);
          gPad->SetTopMargin(0.1);
          gPad->SetBottomMargin(0.2);
          
          // Find the eta projection from the final correction histogram
          sprintf(histogramNamer,"validationEtaProjection%d%d%d",iJetTrack,iCentrality,iTrackPt);
          deltaEtaProjectionFromCorrection = methods->ProjectRegionDeltaEta(spilloverDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt],-1.5,1.5,histogramNamer);
          deltaEtaProjectionFromCorrection->Scale(methods->GetNBinsProjectedOver());
          
          // Draw the hydjet and correction histograms to the same figure
          spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->Draw();
          deltaEtaProjectionFromCorrection->SetLineColor(kRed);
          deltaEtaProjectionFromCorrection->SetLineWidth(3);
          deltaEtaProjectionFromCorrection->Rebin(deltaEtaRebinSpillover);
          deltaEtaProjectionFromCorrection->Scale(1.0/deltaEtaRebinSpillover);
          deltaEtaProjectionFromCorrection->Draw("same");
          
          // Find the correct pad inside the canvas
          deltaPhiCanvas[iJetTrack]->cd(nCentralityBins-1-iCentrality+nCentralityBins*iTrackPt+1);
          gPad->SetTopMargin(0.1);
          gPad->SetBottomMargin(0.2);
          
          // Rebin the deltaPhi histograms
          spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->Rebin(deltaPhiRebinSpillover);
          spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->Scale(1.0/deltaPhiRebinSpillover);
          spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt]->Rebin(deltaPhiRebinSpillover);
          spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt]->Scale(1.0/deltaPhiRebinSpillover);
          
          // Draw the histogram to canvas
          setHistogramStyle(spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt],1.5,spilloverScalePhi,titleString,"#Delta#phi");
          spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt]->SetLineColor(kRed);
          spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->DrawCopy();
          spilloverPhiIntegral = spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->Integral(spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->FindBin(-1.5+0.001),spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->FindBin(1.5-0.001),"width");
          
          legend = new TLegend(0.61,0.6,0.9,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.08);legend->SetTextFont(62);
          legend->AddEntry(spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt],"Hydjet","l");
          legend->AddEntry(spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt],"Fit","l");
          legend->Draw();
          
          // Print out a debug message:
          //cout << "Debugging announcement for " << dummyManager->GetJetTrackHistogramName(iJetTrack) << " spillover bin " << titleString.Data() << endl;
          //cout << "DeltaEta integral: " << spilloverEtaIntegral << " and deltaPhi integral: " << spilloverPhiIntegral << endl;
          //cout << endl;
          
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
          
          // Change to the validation canvas
          deltaPhiValidationCanvas[iJetTrack]->cd(nCentralityBins-1-iCentrality+nCentralityBins*iTrackPt+1);
          gPad->SetTopMargin(0.1);
          gPad->SetBottomMargin(0.2);
          
          // Find the eta projection from the final correction histogram
          sprintf(histogramNamer,"validationPhiProjection%d%d%d",iJetTrack,iCentrality,iTrackPt);
          deltaPhiProjectionFromCorrection = methods->ProjectRegionDeltaPhi(spilloverDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt],0,1.5,histogramNamer);
          deltaPhiProjectionFromCorrection->Scale(methods->GetNBinsProjectedOver());
          
          // Draw the hydjet and correction histograms to the same figure
          spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->Draw();
          deltaPhiProjectionFromCorrection->SetLineColor(kRed);
          deltaPhiProjectionFromCorrection->SetLineWidth(3);
          deltaPhiProjectionFromCorrection->Rebin(deltaPhiRebinSpillover);
          deltaPhiProjectionFromCorrection->Scale(1.0/deltaPhiRebinSpillover);
          deltaPhiProjectionFromCorrection->Draw("same");
          
          // Read yields and widths from the fits
          deltaEtaFitWidth[0][iJetTrack][iCentrality][iTrackPt] = spilloverDeltaEtaFit[iJetTrack][iCentrality][iTrackPt]->GetParameter(1);
          deltaEtaFitWidth[1][iJetTrack][iCentrality][iTrackPt] = spilloverDeltaEtaFit[iJetTrack][iCentrality][iTrackPt]->GetParError(1);
          //deltaEtaFitYield[0][iJetTrack][iCentrality][iTrackPt] = spilloverDeltaEtaFit[iJetTrack][iCentrality][iTrackPt]->GetParameter(0); // Directly the Gauss yield parameter
          deltaEtaFitYield[0][iJetTrack][iCentrality][iTrackPt] =  spilloverDeltaEtaFit[iJetTrack][iCentrality][iTrackPt]->Integral(-1.5,1.5)-3*spilloverDeltaEtaFit[iJetTrack][iCentrality][iTrackPt]->GetParameter(2); // Integral of the fit - constant background
          deltaEtaFitYield[1][iJetTrack][iCentrality][iTrackPt] = spilloverDeltaEtaFit[iJetTrack][iCentrality][iTrackPt]->GetParError(0);
          deltaPhiFitWidth[0][iJetTrack][iCentrality][iTrackPt] = spilloverDeltaPhiFit[iJetTrack][iCentrality][iTrackPt]->GetParameter(1);
          deltaPhiFitWidth[1][iJetTrack][iCentrality][iTrackPt] = spilloverDeltaPhiFit[iJetTrack][iCentrality][iTrackPt]->GetParError(1);
          //deltaPhiFitYield[0][iJetTrack][iCentrality][iTrackPt] = spilloverDeltaPhiFit[iJetTrack][iCentrality][iTrackPt]->GetParameter(0);
          deltaPhiFitYield[0][iJetTrack][iCentrality][iTrackPt] = spilloverDeltaPhiFit[iJetTrack][iCentrality][iTrackPt]->Integral(-1.5,1.5)-3*spilloverDeltaPhiFit[iJetTrack][iCentrality][iTrackPt]->GetParameter(2); // Integral of the fit - constant background
          deltaPhiFitYield[1][iJetTrack][iCentrality][iTrackPt] = spilloverDeltaPhiFit[iJetTrack][iCentrality][iTrackPt]->GetParError(0);
          
          // XXXXX Some debuggery about the eta yield
          //cout << "JetTrack: " << iJetTrack << " Centrality: " << iCentrality << " TrackpT: " << iTrackPt << endl;
          //cout << "Delta eta yield from fit: " << deltaEtaFitYield[0][iJetTrack][iCentrality][iTrackPt] << endl;
          //cout << "Delta eta yield from bin counting: " << spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->Integral(spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->FindBin(-1.5+0.001),spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->FindBin(1.5-0.001),"width")-3*spilloverDeltaEtaFit[iJetTrack][iCentrality][iTrackPt]->GetParameter(2) << endl;
          
          // Normalize the yields to pT bin width
          deltaEtaFitYield[0][iJetTrack][iCentrality][iTrackPt] = deltaEtaFitYield[0][iJetTrack][iCentrality][iTrackPt]/(trackPtBinBorders[iTrackPt+1]-trackPtBinBorders[iTrackPt]);
          deltaEtaFitYield[1][iJetTrack][iCentrality][iTrackPt] = deltaEtaFitYield[1][iJetTrack][iCentrality][iTrackPt]/(trackPtBinBorders[iTrackPt+1]-trackPtBinBorders[iTrackPt]);
          deltaPhiFitYield[0][iJetTrack][iCentrality][iTrackPt] = deltaPhiFitYield[0][iJetTrack][iCentrality][iTrackPt]/(trackPtBinBorders[iTrackPt+1]-trackPtBinBorders[iTrackPt]);
          deltaPhiFitYield[1][iJetTrack][iCentrality][iTrackPt] = deltaPhiFitYield[1][iJetTrack][iCentrality][iTrackPt]/(trackPtBinBorders[iTrackPt+1]-trackPtBinBorders[iTrackPt]);
          
          // Find a good place to put the track pT points for the graphs
          lowPtBin = trackPtForGraphs[iCentrality]->FindBin(trackPtBinBorders[iTrackPt]);
          highPtBin = trackPtForGraphs[iCentrality]->FindBin(trackPtBinBorders[iTrackPt+1]);
          trackPtForGraphs[iCentrality]->GetXaxis()->SetRange(lowPtBin,highPtBin);
          graphPointsX[iTrackPt] = trackPtForGraphs[iCentrality]->GetMean();
          
        } // track pT
        
        // Do not do the spillover graphs for subleading jets, as there is no peak on the subleading side in Hydjet
        if(iJetTrack < DijetHistogramManager::kTrackSubleadingJet || iJetTrack > DijetHistogramManager::kPtWeightedTrackSubleadingJet){
          
          double originalValueX, originalValueY, originalErrorX, originalErrorY;
          // Exclude some bins from the fit
          lowFitPt = 0.5;
          if(iCentrality > 0) lowFitPt = 1;
          
          // Put the obtained fit parameters to graphs
          graphDeltaEtaWidth[iJetTrack][iCentrality] = new TGraphErrors(nTrackPtBins,graphPointsX,deltaEtaFitWidth[0][iJetTrack][iCentrality],graphErrorsX,deltaEtaFitWidth[1][iJetTrack][iCentrality]);
          graphDeltaEtaWidth[iJetTrack][iCentrality]->Fit("pol1","","",lowFitPt,8);
          graphDeltaEtaYield[iJetTrack][iCentrality] = new TGraphErrors(nTrackPtBins,graphPointsX,deltaEtaFitYield[0][iJetTrack][iCentrality],graphErrorsX,deltaEtaFitYield[1][iJetTrack][iCentrality]);
          
          // Remove one point from fit for 0-10 centrality bin
          if(iCentrality == 0){
            graphDeltaEtaYield[iJetTrack][iCentrality]->GetPoint(3,originalValueX,originalValueY);
            originalErrorY = graphDeltaEtaYield[iJetTrack][iCentrality]->GetErrorY(3);
            originalErrorX = graphDeltaEtaYield[iJetTrack][iCentrality]->GetErrorX(3);
            graphDeltaEtaYield[iJetTrack][iCentrality]->RemovePoint(3);
          }
          
          graphDeltaEtaYield[iJetTrack][iCentrality]->Fit("expo","","",lowFitPt,8);
          
          // Return the point to the graph after fit
          if(iCentrality == 0) {
            graphDeltaEtaYield[iJetTrack][iCentrality]->SetPoint(5,originalValueX,originalValueY);
            graphDeltaEtaYield[iJetTrack][iCentrality]->SetPointError(5,originalErrorX,originalErrorY);
          }
          
          graphDeltaPhiWidth[iJetTrack][iCentrality] = new TGraphErrors(nTrackPtBins,graphPointsX,deltaPhiFitWidth[0][iJetTrack][iCentrality],graphErrorsX,deltaPhiFitWidth[1][iJetTrack][iCentrality]);
          graphDeltaPhiWidth[iJetTrack][iCentrality]->Fit("pol1","","",lowFitPt,8);
          graphDeltaPhiYield[iJetTrack][iCentrality] = new TGraphErrors(nTrackPtBins,graphPointsX,deltaPhiFitYield[0][iJetTrack][iCentrality],graphErrorsX,deltaPhiFitYield[1][iJetTrack][iCentrality]);
          
          // Remove one point from fit for 0-10 centrality bin
          if(iCentrality == 0){
            graphDeltaPhiYield[iJetTrack][iCentrality]->GetPoint(3,originalValueX,originalValueY);
            originalErrorY = graphDeltaPhiYield[iJetTrack][iCentrality]->GetErrorY(3);
            originalErrorX = graphDeltaPhiYield[iJetTrack][iCentrality]->GetErrorX(3);
            graphDeltaPhiYield[iJetTrack][iCentrality]->RemovePoint(3);
          }
          
          graphDeltaPhiYield[iJetTrack][iCentrality]->Fit("expo","","",lowFitPt,8);
          
          // Return the point to the graph after fit
          if(iCentrality == 0) {
            graphDeltaPhiYield[iJetTrack][iCentrality]->SetPoint(5,originalValueX,originalValueY);
            graphDeltaPhiYield[iJetTrack][iCentrality]->SetPointError(5,originalErrorX,originalErrorY);
          }
          
        } // Not a subleading jet if
        
        
      } // centrality
      
      // After all the canvases are filled, save them
      if(saveFigures) {
        sprintf(histogramNamer,"figures/%s_spilloverDeltaEta%s.pdf",dummyManager->GetJetTrackHistogramName(iJetTrack),saveNameComment);
        deltaEtaCanvas[iJetTrack]->SaveAs(histogramNamer);
        
        sprintf(histogramNamer,"figures/%s_spilloverDeltaEtaComparison%s.pdf",dummyManager->GetJetTrackHistogramName(iJetTrack),saveNameComment);
        deltaEtaComparisonCanvas[iJetTrack]->SaveAs(histogramNamer);
        
        sprintf(histogramNamer,"figures/%s_spilloverDeltaEtaValidation%s.pdf",dummyManager->GetJetTrackHistogramName(iJetTrack),saveNameComment);
        deltaEtaValidationCanvas[iJetTrack]->SaveAs(histogramNamer);
        
        sprintf(histogramNamer,"figures/%s_spilloverDeltaPhi%s.pdf",dummyManager->GetJetTrackHistogramName(iJetTrack),saveNameComment);
        deltaPhiCanvas[iJetTrack]->SaveAs(histogramNamer);
        
        sprintf(histogramNamer,"figures/%s_spilloverDeltaPhiComparison%s.pdf",dummyManager->GetJetTrackHistogramName(iJetTrack),saveNameComment);
        deltaPhiComparisonCanvas[iJetTrack]->SaveAs(histogramNamer);
        
        sprintf(histogramNamer,"figures/%s_spilloverDeltaPhiValidation%s.pdf",dummyManager->GetJetTrackHistogramName(iJetTrack),saveNameComment);
        deltaPhiValidationCanvas[iJetTrack]->SaveAs(histogramNamer);
        
      } // saving figures
      
      // Draw the graphs containing the spillover fit parameters if we are not looking at subleading jets
      if(iJetTrack < DijetHistogramManager::kTrackSubleadingJet || iJetTrack > DijetHistogramManager::kPtWeightedTrackSubleadingJet){
        
        if(iJetTrack == DijetHistogramManager::kTrackLeadingJet || iJetTrack == DijetHistogramManager::kTrackInclusiveJet) yZoom = 2;
        if(iJetTrack == DijetHistogramManager::kPtWeightedTrackLeadingJet || iJetTrack == DijetHistogramManager::kPtWeightedTrackInclusiveJet) yZoom = 2.5;
        
        sprintf(saveName,"spilloverEtaFitWidthCentralityComparison%s",saveNameComment);
        drawSpilloverGraph(drawer,graphDeltaEtaWidth[iJetTrack],iJetTrack,yZoom,"Spillover #Delta#eta width",false,saveFigures,saveName,NULL,"pol1");
        sprintf(saveName,"spilloverEtaFitWidthCentralityAverage%s",saveNameComment);
        drawSpilloverGraph(drawer,graphDeltaEtaWidth[iJetTrack],iJetTrack,yZoom,"Spillover #Delta#eta width",false,saveFigures,saveName,combinedGraphDeltaEtaWidth[iJetTrack],"pol1");
        sprintf(saveName,"spilloverEtaFitYieldCentralityComparison%s",saveNameComment);
        drawSpilloverGraph(drawer,graphDeltaEtaYield[iJetTrack],iJetTrack,yZoom+4,"Spillover #Delta#eta yield",true,saveFigures,saveName,NULL,"expo");
        sprintf(saveName,"spilloverPhiFitWidthCentralityComparison%s",saveNameComment);
        drawSpilloverGraph(drawer,graphDeltaPhiWidth[iJetTrack],iJetTrack,yZoom,"Spillover #Delta#phi width",false,saveFigures,saveName,NULL,"pol1");
        sprintf(saveName,"spilloverPhiFitWidthCentralityAverage%s",saveNameComment);
        drawSpilloverGraph(drawer,graphDeltaPhiWidth[iJetTrack],iJetTrack,yZoom,"Spillover #Delta#phi width",false,saveFigures,saveName,combinedGraphDeltaPhiWidth[iJetTrack],"pol1");
        sprintf(saveName,"spilloverPhiFitYieldCentralityComparison%s",saveNameComment);
        drawSpilloverGraph(drawer,graphDeltaPhiYield[iJetTrack],iJetTrack,yZoom+4,"Spillover #Delta#phi yield",true,saveFigures,saveName,NULL,"expo");
        
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          sprintf(saveName,"spilloverFitWidthEtaPhiComparison%s",saveNameComment);
          drawSpilloverGraphComparison(drawer,graphDeltaEtaWidth[iJetTrack][iCentrality],graphDeltaPhiWidth[iJetTrack][iCentrality],iJetTrack,iCentrality,yZoom,"Spillover fit width",false,saveFigures,saveName,NULL,NULL,"pol1");
          sprintf(saveName,"spilloverFitYieldEtaPhiComparison%s",saveNameComment);
          drawSpilloverGraphComparison(drawer,graphDeltaEtaYield[iJetTrack][iCentrality],graphDeltaPhiYield[iJetTrack][iCentrality],iJetTrack,iCentrality,yZoom+4,"Spillover fit yield",true,saveFigures,saveName,NULL,binCountGraphYield[iJetTrack][iCentrality],"expo");
          sprintf(saveName,"spilloverFitYieldEtaPhiAverage%s",saveNameComment);
          drawSpilloverGraphComparison(drawer,graphDeltaEtaYield[iJetTrack][iCentrality],graphDeltaPhiYield[iJetTrack][iCentrality],iJetTrack,iCentrality,yZoom+4,"Spillover fit yield",true,saveFigures,saveName,combinedGraphYield[iJetTrack][iCentrality],binCountGraphYield[iJetTrack][iCentrality],"expo");
        } // centrality loop
        
      } // not subleading jet if
      
    } // jet track loop
    
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
          if(iCentrality != nCentralityBins){
            seagullDeltaEtaWings[iJetTrack][iCentrality][iTrackPt]->Rebin(2);
            seagullDeltaEtaWings[iJetTrack][iCentrality][iTrackPt]->Scale(1.0/2.0);
          }
          seagullDeltaEtaWings[iJetTrack][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-3,3);
          setHistogramTitles(seagullDeltaEtaWings[iJetTrack][iCentrality][iTrackPt],titleString,"#Delta#eta");
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
    TCanvas *jffAsymmetryComparisonCanvas[DijetHistogramManager::knJetTrackCorrelations];
    TCanvas *jffComparisonCanvas[2]; // One for regular and one for pT-weighted
    TCanvas *spilloverCorrectionCanvas[DijetHistogramManager::knJetTrackCorrelations];
    
    DijetMethods *jetShapeCalculator = new DijetMethods();
    drawer->Reset();
    double max1, max2, min1, min2, theMax, theMin;
    double asymmetryMin[10], asymmetryMax[10];
    double legendX1, legendY1;
    
                         // Track pT    0.7-1   1-2    2-3    3-4    4-8  8-300
    double defaultLowJffScale[8][6] = {{-0.4,  -1.25, -1.05, -0.95, -2.3, -7.8},  // Track-leading jet
                                       {-0.4,  -1.25, -1.05, -0.95, -2.3, -7.8},  // Uncorrected track-leading jet
                                       {-0.36, -1.9,  -2.6,  -3.4,   -13, -170},  // pT weighted track-leading jet
                                       {-0.12, -0.26, -0.3,  -0.36, -1.7, -7.5},  // Track-subleading jet
                                       {-0.08, -0.24, -0.3,  -0.36, -1.7, -7.5},  // Uncorrected track-subleading jet
                                       {-0.14, -0.44, -1.2,  -1.4,   -11, -190},  // pT weighted track-subleading jet
                                       {-0.36, -1.1,   -1,   -0.82, -2.1,  -7},   // Track-inclusive jet
                                       {-0.3,  -1.6,  -2.4,  -2.8,  -12.5,-150}}; // pT weighted track-inclusive jet
    
                          // Track pT    0.7-1   1-2    2-3    3-4    4-8  8-300
    double defaultHighJffScale[8][6] = {{ 0.07,  0.12,  0.1,   0.1,   0.2,  3.8},  // Track-leading jet
                                        { 0.07,  0.12,  0.1,   0.1,   0.2,  3.8},  // Uncorrected track-leading jet
                                        { 0.05,  0.2,   0.3,   0.6,    2,   120},  // pT weighted track-leading jet
                                        { 0.16,  0.4,   0.38,  0.28,  0.9,  3.5},  // Track-subleading jet
                                        { 0.14,  0.4,   0.38,  0.28,  0.9,  3.5},  // Uncorrected track-subleading jet
                                        { 0.14,  0.6,   1.0,   1.0,   5.5,   80},  // pT weighted track-subleading jet
                                        { 0.04,  0.08,  0.1,   0.1,   0.3,  3.8},  // Track-inclusive jet
                                        { 0.03,  0.12,  0.2,   0.3,   1.4,  110}}; // pT weighted track-inclusive jet
    
    const char* jffLegendName[] = {"Leading jet JFF","Leading jet JFF","p_{T}w leading jet JFF","Subleading jet JFF","Subleading jet JFF","p_{T}w subleading jet JFF","Inclusive jet JFF","p_{T}w inclusive jet JFF"};
    
    int asymmetryStyle[] = {29,33,47,34,20,21,43,45,41};
    int asymmetryColor[] = {kGreen+3,kBlue,kRed,kBlack,kMagenta,kCyan,kOrange+4,kBlue-7,kSpring-3};
    
    for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
      if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
      
      // Create one big canvas with a pad for each centrality and track pT bin
      sprintf(histogramNamer,"jffCorrection%d",iJetTrack);
      sprintf(padNamer,"JFF correction %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
      jffCorrectionCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1250,1800);
      jffCorrectionCanvas[iJetTrack]->Divide(nCentralityBins+1,nTrackPtBins);
      
      // Create comparison canvases only for leading jet correlations
      if(iJetTrack == DijetHistogramManager::kTrackLeadingJet || iJetTrack == DijetHistogramManager::kPtWeightedTrackLeadingJet){
        sprintf(histogramNamer,"jffComparison%d",iJetTrack);
        sprintf(padNamer,"JFF comparison %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
        jffComparisonCanvas[iJetTrack/2] = new TCanvas(histogramNamer,padNamer,1250,1800);
        jffComparisonCanvas[iJetTrack/2]->Divide(nCentralityBins+1,nTrackPtBins);
      }
      
      // Create asymmetry comparison canvases for leading and subleading jets, but not for inclusive jets
      if(iJetTrack < DijetHistogramManager::kTrackInclusiveJet){
        sprintf(histogramNamer,"jffAsymmetryComparison%d",iJetTrack);
        sprintf(padNamer,"JFF asymmetry comparison %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
        jffAsymmetryComparisonCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1250,1800);
        jffAsymmetryComparisonCanvas[iJetTrack]->Divide(nCentralityBins+1,nTrackPtBins);
      }
      
      sprintf(histogramNamer,"spilloverCorrection%d",iJetTrack);
      sprintf(padNamer,"Spillover correction %s",dummyManager->GetJetTrackHistogramName(iJetTrack));
      spilloverCorrectionCanvas[iJetTrack] = new TCanvas(histogramNamer,padNamer,1000,1800);
      spilloverCorrectionCanvas[iJetTrack]->Divide(nCentralityBins,nTrackPtBins);
      
      for(int iCentrality = 0; iCentrality < nCentralityBins+1; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          // Calculate jet shape from the two-dimensional spillover distribution
          if(iCentrality < nCentralityBins){ // No spillover correction for pp
            correctionJetShape[1][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt] = jetShapeCalculator->GetJetShape(spilloverDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]);
           }
          
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
            jffCorrectionCanvas[iJetTrack]->cd(nCentralityBins-iCentrality+(nCentralityBins+1)*iTrackPt+1);
            gPad->SetTopMargin(0.1);
            gPad->SetBottomMargin(0.2);
            
            // Draw the JFF correction histogram to the canvas
            gStyle->SetTitleSize(0.09,"t");
            setHistogramTitles(correctionJetShape[0][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt],titleString,"#Deltar");
            correctionJetShape[0][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt]->SetMarkerStyle(34);
            correctionJetShape[0][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt]->SetMarkerColor(kBlack);
            if(constantBigCanvasScale){
              correctionJetShape[0][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(defaultLowJffScale[iJetTrack][iTrackPt],defaultHighJffScale[iJetTrack][iTrackPt]);
            }
            correctionJetShape[0][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt]->DrawCopy();
            
            // Draw a legend only to the last canvas
            if(nCentralityBins-iCentrality+(nCentralityBins+1)*iTrackPt+1 == 30){
              legend = new TLegend(0.15,0.3,0.9,0.5);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.11);legend->SetTextFont(62);
              legend->AddEntry(correctionJetShape[0][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt],jffLegendName[iJetTrack],"p");
              legend->Draw();
            }
            
            // Draw to comparison canvases only for leading jet correlations
            if(iJetTrack == DijetHistogramManager::kTrackLeadingJet || iJetTrack == DijetHistogramManager::kPtWeightedTrackLeadingJet){
              
              // First, find the drawing scale for histograms
              theMin = defaultLowJffScale[iJetTrack][iTrackPt];
              theMax = defaultHighJffScale[iJetTrack][iTrackPt];
              if(defaultLowJffScale[iJetTrack+3][iTrackPt] < theMin) theMin = defaultLowJffScale[iJetTrack+3][iTrackPt];
              if(defaultHighJffScale[iJetTrack+3][iTrackPt] > theMax) theMax = defaultHighJffScale[iJetTrack+3][iTrackPt];
              correctionJetShape[0][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(theMin,theMax);
              
              // Do some styling for the subleading histogram
              correctionJetShape[0][iJetTrack+3][nJffAsymmetryBins][iCentrality][iTrackPt]->SetMarkerStyle(47);
              correctionJetShape[0][iJetTrack+3][nJffAsymmetryBins][iCentrality][iTrackPt]->SetMarkerColor(kRed);
              
              // Find the correct canvas to draw the histograms
              jffComparisonCanvas[iJetTrack/2]->cd(nCentralityBins-iCentrality+(nCentralityBins+1)*iTrackPt+1);
              gPad->SetTopMargin(0.1);
              gPad->SetBottomMargin(0.2);
              
              // Draw the leading and subleading JFF corrections to the same canvas
              gStyle->SetTitleSize(0.09,"t");
              correctionJetShape[0][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt]->DrawCopy();
              correctionJetShape[0][iJetTrack+3][nJffAsymmetryBins][iCentrality][iTrackPt]->DrawCopy("same");
              
              // Add legends to specific canvases
              if(nCentralityBins-iCentrality+(nCentralityBins+1)*iTrackPt+1 == 29){
                drawBigCanvasLegend(legend,correctionJetShape[0][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt],jffLegendName[iJetTrack]);
              } else if(nCentralityBins-iCentrality+(nCentralityBins+1)*iTrackPt+1 == 30){
                drawBigCanvasLegend(legend,correctionJetShape[0][iJetTrack+3][nJffAsymmetryBins][iCentrality][iTrackPt],jffLegendName[iJetTrack+3]);
              }
              
            } // If for leading-subleading comparison
            
            // Do not draw the asymmetry canvases for inclusive jet correlations
            if(iJetTrack < DijetHistogramManager::kTrackInclusiveJet){
              
              // Set a nice drawing style for the histograms
              if(iCentrality != nCentralityBins){ // XXXX TEMP TEMP TODO Add xj bins also to pp file
                // Set the style
                for(int iAsymmetry = 0; iAsymmetry <= nJffAsymmetryBins; iAsymmetry++){
                  correctionJetShape[0][iJetTrack][iAsymmetry][iCentrality][iTrackPt]->SetMarkerStyle(asymmetryStyle[iAsymmetry]);
                  correctionJetShape[0][iJetTrack][iAsymmetry][iCentrality][iTrackPt]->SetMarkerColor(asymmetryColor[iAsymmetry]);
                }
              }
              
              // Determine the drawing range for all bins from the first centrality bin
              if(iCentrality == 0){
                min1 = 1000;
                max1 = -1000;
                for(int iAsymmetry = 0; iAsymmetry <= nJffAsymmetryBins; iAsymmetry++){
                  
                  // Find the minimum and maximum from all asymmetry bins
                  for(int iBin = 1; iBin <= correctionJetShape[0][iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetNbinsX(); iBin++){
                    min2 = correctionJetShape[0][iJetTrack][iAsymmetry][iCentrality][iTrackPt]->GetBinContent(iBin);
                    if(min1 > min2){
                      min1 = min2;
                    }
                    if(max1 < min2){
                      max1 = min2;
                    }
                  } // Bins in histogram loop
                } // Asymmetry bin loop
                
                // After the minimum and maximum values over asymmetry bins are determined, leave 10 % marginal for drawing
                asymmetryMin[iTrackPt] = min1 - 0.1*(max1-min1);
                asymmetryMax[iTrackPt] = max1 + 0.1*(max1-min1);
              } // If for the first centrality bin
              
              // Find the correct canvas to draw the histograms
              jffAsymmetryComparisonCanvas[iJetTrack]->cd(nCentralityBins-iCentrality+(nCentralityBins+1)*iTrackPt+1);
              gPad->SetTopMargin(0.1);
              gPad->SetBottomMargin(0.2);
              
              // First draw the asymmetry inclusive distribution to the canvas
              gStyle->SetTitleSize(0.09,"t");
              correctionJetShape[0][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(asymmetryMin[iTrackPt],asymmetryMax[iTrackPt]);
              correctionJetShape[0][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt]->DrawCopy();
              
              // Then, draw all the asymmetry bins to the same canvas
              if(iCentrality != nCentralityBins){ // XXXX TEMP TEMP TODO Add xj bins also to pp file
                for(int iAsymmetry = 0; iAsymmetry < nJffAsymmetryBins; iAsymmetry++){
                  correctionJetShape[0][iJetTrack][iAsymmetry][iCentrality][iTrackPt]->DrawCopy("same");
                }
              }
              
              // Draw different asymmetry legends to different canvases
              for(int iAsymmetry = 0; iAsymmetry <= nJffAsymmetryBins; iAsymmetry++){
                if(nCentralityBins-iCentrality+(nCentralityBins+1)*iTrackPt+1 == 30-iAsymmetry){
                  drawBigCanvasLegend(legend,correctionJetShape[0][iJetTrack][iAsymmetry][iCentrality][iTrackPt],jffAsymmetryLegend[iAsymmetry].Data());
                }
              }
              
              // Draw the system legend to the pad before the asymmetry legends
              if(nCentralityBins-iCentrality+(nCentralityBins+1)*iTrackPt+1 == 30-nJffAsymmetryBins-1){
                legend = new TLegend(0.15,0.3,0.9,0.5);
                legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.11);legend->SetTextFont(62);
                legend->AddEntry((TObject*)0,jffLegendName[iJetTrack],"");
                legend->Draw();
              }
              
            } // Asymmetry correction canvas if
            
            if(iCentrality == nCentralityBins) continue; // No spillover correction for pp
            
            // Change to canvas for spillover correction
            spilloverCorrectionCanvas[iJetTrack]->cd(nCentralityBins-1-iCentrality+nCentralityBins*iTrackPt+1);
            gPad->SetTopMargin(0.1);
            gPad->SetBottomMargin(0.2);
            
            // Draw the JFF correction histogram to the canvas
            gStyle->SetTitleSize(0.09,"t");
            setHistogramTitles(correctionJetShape[1][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt],titleString,"#Deltar");
            if(constantBigCanvasScale){
              if(iJetTrack == DijetHistogramManager::kPtWeightedTrackLeadingJet){
                correctionJetShape[1][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0,3.4);
              } else {
                correctionJetShape[1][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0,2.4);
              }
            }
            correctionJetShape[1][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt]->DrawCopy();
          
          } // Big canvas if
          
          //////////////////////////////////////////////////
          // Next, draw comparison plots in small canvases
          //////////////////////////////////////////////////
          
          if(jetShapeCorrectionComparison){
            
            if(iCentrality == nCentralityBins) continue; // No spillover correction for pp
            
            // Gether information of the current bin
            if(lastCentralityBin == 0){
              titleString = Form("%s - pp",dummyManager->GetJetTrackAxisName(iJetTrack));
            } else {
              titleString = Form("%s - Cent: %.0f-%.0f%%",dummyManager->GetJetTrackAxisName(iJetTrack),centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
            }
            
            sprintf(padNamer,"Track pT: %.1f-%.1f GeV",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
            
            // Find a good drawing range for the distributions
            std::tie(min1,max1) = getDrawMinMax(correctionJetShape[0][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt]);
            std::tie(min2,max2) = getDrawMinMax(correctionJetShape[1][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt]);
            theMin = (min1 < min2) ? min1 : min2;
            theMax = (max1 > max2) ? max1 : max2;
            
            // Draw the histogram to canvas together with the one from file
            correctionJetShape[0][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(theMin,theMax);
            drawer->DrawHistogram(correctionJetShape[0][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt], "#Deltar","P(#Deltar)", " ");
            correctionJetShape[1][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(kRed);
            correctionJetShape[1][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt]->Draw("same");
            
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
            legend->AddEntry(correctionJetShape[0][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt],"JFF correction","l");
            legend->AddEntry(correctionJetShape[1][iJetTrack][nJffAsymmetryBins][iCentrality][iTrackPt],"Spillover correction","l");
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
        
        if(iJetTrack == DijetHistogramManager::kTrackLeadingJet){
          jffComparisonCanvas[0]->SaveAs("figures/leadingSubleadingJffComparison.pdf");
        } else if(iJetTrack == DijetHistogramManager::kPtWeightedTrackLeadingJet){
          jffComparisonCanvas[1]->SaveAs("figures/leadingSubleadingPtWeightedJffComparison.pdf");
        }
        
      } // saving figures
      
    } // Jet-track category loop
    
  } // Jet shape plots from JFF correction
}
