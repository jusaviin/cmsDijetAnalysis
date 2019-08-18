#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "DijetMethods.h"
#include "JDrawer.h"
#include "JffCorrector.h"

/*
 * Print a latex slide with given numerical information in xj bins to the console
 *
 *  Arguments:
 *   double (*numberArray)[5][7] = Array of numbers with correct asymmetry, centrality and pT binning to be printed to console
 *   int nAsymmetryBins = Number of asymmetry bins for the data
 *   const char* title = Title string to be printed on the slides
 */
void printNumberSlide(double (*numberArray)[5][7], int nAsymmetryBins, const char* title){
  
  // Create a histogram manager to facilitate binning info exchange
  const int nTrackPtBins = 6;
  const int nCentralityBins = 4;
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12};  // Bin borders for track pT
  double asymmetryBinBorders[] = {0,0.6,0.8,1};     // Bin borders for xj
  
  char namer[100];
  
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
  
    if(iAsymmetry == nAsymmetryBins){
      sprintf(namer,"\\frametitle{%s}",title);
    } else {
      sprintf(namer,"\\frametitle{%s, $%.1f < x_{j} < %.1f$}", title, asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1]);
    }
    
    cout << endl;
    cout << "\\begin{frame}" << endl;
    cout << namer << endl;
    cout << "\\begin{center}" << endl;
    cout << "  \\begin{tabular}{cccccc}" << endl;
    cout << "    \\toprule" << endl;
    cout << "    $p_{\\mathrm{T}} (GeV)$ & C: 0-10 & C: 10-30 & C: 30-50 & C: 50-100 & pp \\\\" << endl;
    cout << "    \\midrule" << endl;
    
    // Set the correct precision for printing floating point numbers
    cout << fixed << setprecision(3);
    
    // Print one line for each track pT bin
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      sprintf(namer,"    %.1f-%.1f",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
      cout << namer;
      for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
        cout << " & $" << numberArray[iAsymmetry][iCentrality][iTrackPt] << "$";
      }
      cout << " \\\\" << endl;
    }
    
    cout << "    \\bottomrule" << endl;
    cout << "  \\end{tabular}" << endl;
    cout << "\\end{center}" << endl;
    cout << "\\end{frame}" << endl;
    cout << endl;
    
  }
}

/*
 * Macro for plotting long range correlation results.
 *
 * Depending on the configuration of the macro, the following things can be plotted/printed:
 *  - Refit the long range correlation distribution with any order Fourier fit
 *  - Plot the total Fourier fit to different asymmetry bins
 *  - Plot only one component of the fit different asymmetry bins
 *  - Print a slide with chi2/NFD information of the Fourier fit
 *  - Print a slide with systematic uncertainty estimation coming from the background adjustment
 *    The systematic uncertainty is estimated as half the amount of shift needed to match leading and subleading sides
 */
void compareLongRangeAsymmetry(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Main files from which the long range asymmetries are obtained
  TString pbpbFileName = "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_modifiedSeagull_noErrorJff_averagePeakMixing_adjustedBackground_processed_2019-08-13_fiveJobsMissing.root";
  //"data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_adjustedBackground_processed_2019-07-05.root";
  TString ppFileName = "data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_adjustedBackground_processed_2019-08-13.root";
  //"data/dijet_pp_highForest_pfJets_noUncOrInc_allCorrections_adjustedBackground_wtaAxis_processed_2019-07-13.root";
  
  // For systematic uncertainty estimation, need files in which the background is not adjusted between leading and subleading side
  TString pbpbUnadjustedFileName = "data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_modifiedSeagull_noErrorJff_averagePeakMixing_processed_2019-08-13_fiveJobsMissing.root";
  //"data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root";
  TString ppUnadjustedFileName = "data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_processed_2019-08-13.root";
  //"data/dijet_pp_highForest_pfJets_noUncOrInc_allCorrections_wtaAxis_processed_2019-07-13.root";
  
  // Name for the file from which systematic uncertainties are read
  const char *uncertaintyFile = "data/vnUncertaintyPreliminaryNoVz2018.txt";
  
  // Open the PbPb input file and read bin numbers from it
  TFile *pbpbDataFile = TFile::Open(pbpbFileName);
  DijetHistogramManager *pbpbReader = new DijetHistogramManager(pbpbDataFile);
  const int nCentralityBins = pbpbReader->GetNCentralityBins();
  const int nTrackPtBins = pbpbReader->GetNTrackPtBins();
  const int nAsymmetryBins = pbpbReader->GetNAsymmetryBins();
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  const bool drawFourierFit = false;
  const bool drawFourierGraph = true;
  
  const bool printChi2 = false;
  const bool printBackgroundAdjustmentUncertainty = false;
  
  const bool saveFigures = true;
  TString saveComment = "";
  
  int firstDrawnAsymmetryBin = 0;
  int lastDrawnAsymmetryBin = nAsymmetryBins-1;
  
  const int fourierV = 0;  // Select which vn component to draw. 0 = All, 1...4 = v1...v4
  TString vString = "";
  TString vHeader = "V_{1} - V_{4}";
  if(fourierV > 0) {
    vString = Form("_v%d",fourierV);
    vHeader = Form("V_{%d}",fourierV);
  }
  
  const int nRefit = 4; // Number of vn:s included in the refit
  bool refitBackground = true; // Refit the background distribution
  int backgroundRebin = 4; // Rebin applied to the background distributions
  
  if(refitBackground && fourierV == 0) vHeader = Form("V_{1} - V_{%d}",nRefit);
  
  // ==================================================================
  // ====================== Configuration done ========================
  // ==================================================================
  
  // Open the pp input file
  TFile *ppDataFile = TFile::Open(ppFileName);
  
  // Open the files needed for systematic uncertainty estimation
  TFile *pbpbUnadjustedFile = TFile::Open(pbpbUnadjustedFileName);
  TFile *ppUnadjustedFile = TFile::Open(ppUnadjustedFileName);
  
  // Create readers for the histograms and define binning
  DijetHistogramManager *ppReader = new DijetHistogramManager(ppDataFile);
  DijetHistogramManager *pbpbUnadjustedReader = new DijetHistogramManager(pbpbUnadjustedFile);
  DijetHistogramManager *ppUnadjustedReader = new DijetHistogramManager(ppUnadjustedFile);
  
  // Read the uncertainties from the uncertainty file
  JffCorrector *uncertaintyProvider = new JffCorrector();
  uncertaintyProvider->ReadLongRangeSystematicFile(uncertaintyFile);
  
  // Define a fitter to test fitting the background distribution with different amount of vn:s
  DijetMethods *refitter = new DijetMethods();
  
  // Define arrays for the histograms
  TH1D *background[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  TF1 *backgroundFit[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  
  // Define arrays needed for systematic uncertainty estimation
  TH1D *leadingUnadjustedBackground[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  TH1D *leadingUnadjustedBackgroundOverlap[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  
  // Define histograms to draw the systematic uncertainties
  TH1D *adjustmentHistogram[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  
  // We will fit also the unadjusted histogram to estimate systematic uncertainties
  TF1 *unadjustedFit[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  TF1 *uncertaintyFunctionLow, *uncertaintyFunctionHigh;
  
  // Define track histograms that are needed to find correct place to put point for the graphs
  TH1D *tracksForGraph[nCentralityBins+1];
  
  // Define arrays for extracted vn numbers
  double masterFlowTable[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double masterFlowError[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // Read the numbers also for unadjusted fits to estimate systematic uncertainties
  double unadjustedFlowTable[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  double unadjustedFlowError[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins][nRefit];
  
  // To pass the array as an argument, the numbers will need to be given by hand.
  if(nAsymmetryBins != 3 || nCentralityBins != 4 || nTrackPtBins != 7){
    cout << "Wrong binning for chi2/ndf and background uncertainty histograms!" << endl;
    cout << "nAsymmetry = " << nAsymmetryBins << " expected = 3" << endl;
    cout << "nCentralityBins = " << nCentralityBins << " expected = 4" << endl;
    cout << "nTrackPtBins = " << nTrackPtBins << " expected = 7" << endl;
    cout << "Fix the binning and the code will run again!" << endl;
    return;
  }
  double chi2PerNdf[4][5][7];
  double backgroundAdjustmentUncertainty[4][5][7];
  double trackEfficiencyUncertainty = 0.05;  // The uncertainty for the tracking efficiency is 5 %
  
  // Read the histograms from the input files
  pbpbReader->SetLoadTracks(true);
  pbpbReader->SetLoadTrackLeadingJetCorrelations(true);
  pbpbReader->SetAsymmetryBinRange(0,nAsymmetryBins);
  pbpbReader->LoadProcessedHistograms();
  
  ppReader->SetLoadTracks(true);
  ppReader->SetLoadTrackLeadingJetCorrelations(true);
  ppReader->SetAsymmetryBinRange(0,nAsymmetryBins);
  ppReader->LoadProcessedHistograms();
  
  pbpbUnadjustedReader->SetLoadTrackLeadingJetCorrelations(true);
  pbpbUnadjustedReader->SetAsymmetryBinRange(0,nAsymmetryBins);
  pbpbUnadjustedReader->LoadProcessedHistograms();
  
  ppUnadjustedReader->SetLoadTrackLeadingJetCorrelations(true);
  ppUnadjustedReader->SetAsymmetryBinRange(0,nAsymmetryBins);
  ppUnadjustedReader->LoadProcessedHistograms();
  
  double chi2;  // Chi2 for fit quality estimation
  int ndf;      // Number of degrees of freedom for fit quality estimation
  double leadingAdjustmentFactor, subleadingAdjustmentFactor;  // Background adjustment factors for systematic unceratinty estimation
  double backgroundLevel;
  double currentV, currentUncertainty;
  double valueLow, valueHigh;
  double binCenter;
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    
    // Track histogram for pp (needed for graph binning)
    if(iCentrality == nCentralityBins){
      
      tracksForGraph[nCentralityBins] = ppReader->GetHistogramTrackPt(DijetHistogramManager::kTrack, DijetHistogramManager::kSameEvent, 0);
      
    // Track histograms for PbPb (needed for graph binning)
    } else {
      
      tracksForGraph[iCentrality] = pbpbReader->GetHistogramTrackPt(DijetHistogramManager::kTrack, DijetHistogramManager::kSameEvent, iCentrality);
      
    }
    
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        
        // Read the histograms for pp
        if(iCentrality == nCentralityBins){
          
          // Regular background histogram
          background[iAsymmetry][nCentralityBins][iTrackPt] = ppReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, iAsymmetry, 0, iTrackPt, DijetHistogramManager::kWholeEta);
          
          // Track-leading jet unadjusted background histogram for systematic uncertainty estimation
          leadingUnadjustedBackground[iAsymmetry][nCentralityBins][iTrackPt] = ppUnadjustedReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, iAsymmetry, 0, iTrackPt, DijetHistogramManager::kWholeEta);
          
          // Track-subleading jet unadjusted background histogram for systematic uncertainty estimation
          leadingUnadjustedBackgroundOverlap[iAsymmetry][nCentralityBins][iTrackPt] = ppUnadjustedReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackgroundOverlap, iAsymmetry, 0, iTrackPt, DijetHistogramManager::kWholeEta);
          
        // Read the histograms for PbPb
        } else {
          
          // Regular background histogram
          background[iAsymmetry][iCentrality][iTrackPt] = pbpbReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta);
          
          // Track-leading jet unadjusted background histogram for systematic uncertainty estimation
          leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt] = pbpbUnadjustedReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta);
          
          // Track-subleading jet unadjusted background histogram for systematic uncertainty estimation
          leadingUnadjustedBackgroundOverlap[iAsymmetry][iCentrality][iTrackPt] = pbpbUnadjustedReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackgroundOverlap, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta);
          
        }
        
        // Get the background adjustment factors from the unadjusted distributions
        std::tie(leadingAdjustmentFactor, subleadingAdjustmentFactor) = refitter->GetBackgroundAdjustmentFactors(leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt], leadingUnadjustedBackgroundOverlap[iAsymmetry][iCentrality][iTrackPt]);
        
        // Use the larger adjustment as the systematic uncertainty
        backgroundAdjustmentUncertainty[iAsymmetry][iCentrality][iTrackPt] = TMath::Max(TMath::Abs(1-leadingAdjustmentFactor), TMath::Abs(1-subleadingAdjustmentFactor));
        
        // Find the fourier fit function from the histogram
        backgroundFit[iAsymmetry][iCentrality][iTrackPt] = background[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        unadjustedFit[iAsymmetry][iCentrality][iTrackPt] = leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        
        // Option to refit the background distribution
        if(refitBackground){
          background[iAsymmetry][iCentrality][iTrackPt]->RecursiveRemove(backgroundFit[iAsymmetry][iCentrality][iTrackPt]);
          background[iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
          background[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
          
          
          //adjustmentHistogram[iAsymmetry][iCentrality][iTrackPt] = (TH1D*) background[iAsymmetry][iCentrality][iTrackPt]->Clone(Form("adjustmentHistogram%d%d%d", iAsymmetry, iCentrality, iTrackPt)); // Adjusted default
          
          refitter->FourierFit(background[iAsymmetry][iCentrality][iTrackPt], nRefit);
          backgroundFit[iAsymmetry][iCentrality][iTrackPt] = background[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
          
          // Refit also the unadjusted background to have same exactly the settings as for adjusted background
          leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt]->RecursiveRemove(unadjustedFit[iAsymmetry][iCentrality][iTrackPt]);
          leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt]->Rebin(backgroundRebin);
          leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt]->Scale(1.0/backgroundRebin);
          refitter->FourierFit(leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt], nRefit);
          unadjustedFit[iAsymmetry][iCentrality][iTrackPt] = leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
          
          // Clone the background histogram for uncertainty drawing
          adjustmentHistogram[iAsymmetry][iCentrality][iTrackPt] = (TH1D*) leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt]->Clone(Form("adjustmentHistogram%d%d%d", iAsymmetry, iCentrality, iTrackPt)); // Unadjusted default
          
        }
        
        // Calculate chi2/ndf from the fourier fit to the background distribution
        chi2 = backgroundFit[iAsymmetry][iCentrality][iTrackPt]->GetChisquare();
        ndf = backgroundFit[iAsymmetry][iCentrality][iTrackPt]->GetNDF();
        chi2PerNdf[iAsymmetry][iCentrality][iTrackPt] = chi2/ndf;
        
        // Set the errors for the uncertainty histogram from background adjustment based on the calculated values
        //backgroundLevel = backgroundFit[iAsymmetry][iCentrality][iTrackPt]->GetParameter(0); // Use adjusted as default value
        backgroundLevel = unadjustedFit[iAsymmetry][iCentrality][iTrackPt]->GetParameter(0); // Use unadjusted as default value

        for(int iBin = 1; iBin <= adjustmentHistogram[iAsymmetry][iCentrality][iTrackPt]->GetNbinsX(); iBin++){
          adjustmentHistogram[iAsymmetry][iCentrality][iTrackPt]->SetBinError(iBin, backgroundLevel*backgroundAdjustmentUncertainty[iAsymmetry][iCentrality][iTrackPt]);
        }
        
      } // asymmetry loop
    } // track pT loop
  } // centrality loop
  
  // If selected, print some information to the console
  if(printChi2) printNumberSlide(chi2PerNdf,nAsymmetryBins,"$\\chi^{2}$ / NDF for long range");
  if(printBackgroundAdjustmentUncertainty) printNumberSlide(backgroundAdjustmentUncertainty,nAsymmetryBins,"Background adjustment uncertainty");
  
  // Extract the vn parameters from the fits
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = 0; iAsymmetry <=nAsymmetryBins; iAsymmetry++){
        for(int iFlow = 0; iFlow < nRefit; iFlow++){
          masterFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = backgroundFit[iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1);
          masterFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = backgroundFit[iAsymmetry][iCentrality][iTrackPt]->GetParError(iFlow+1);
          unadjustedFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow] = unadjustedFit[iAsymmetry][iCentrality][iTrackPt]->GetParameter(iFlow+1);
          unadjustedFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow] = unadjustedFit[iAsymmetry][iCentrality][iTrackPt]->GetParError(iFlow+1);
        } // flow components
      } // Asymmetry loop
    } // track pT
  } // centrality
  
  // Draw all asymmery bins deltaR curves to a single plot to see if the spillover is similar regardless of asymmetry bin
  JDrawer *drawer = new JDrawer();
  
  // How to zoom vn plots    //   0-10 10-30 30-50 50-100  pp
  double vZoomTable[4][5] = {{    0.2,  0.25, 0.35,  1,   1.4},  // v1
                             {    0.12, 0.2,  0.3,  0.5,  0.6}, // v2
                             {    0.12, 0.12, 0.12, 0.2,  0.25},  // v3
                             {    0.1,  0.1,  0.1,  0.1,  0.1}}; // v4
  double multiplier;
  double factors[] = {0.7,0.6,0.4,0.6};
  
  // Helper variables for drawing figures
  TLegend *legend;
  TLegend *vLegend;
  int markers[] = {kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross};
  int fullMarkers[] = {kFullCircle, kFullSquare, kFullDiamond, kFullCross};
  int colors[] = {kBlue,kRed,kGreen+2,kBlack};
  double vx1 = 0.2; double vx2 = 0.25;  // x position of vn information legend
  if(fourierV == 0){
    vx1 = 0.183; vx2 = 0.233;
  }
  char namerY[100];
  TH1D *uncertaintyHistogram = new TH1D("uncertainties","uncertainties",200,-TMath::Pi()/2.0,3*TMath::Pi()/2.0);
  TString legendString;
  
  if(drawFourierFit){
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Setup a legend for the figure
        legend = new TLegend(0.39,0.87-(lastDrawnAsymmetryBin-firstDrawnAsymmetryBin)*0.06,0.75,0.99);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        if(iCentrality == nCentralityBins){
          legend->SetHeader(Form("%.1f < p_{T} < %.1f GeV, pp",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]));
        } else {
          legend->SetHeader(Form("%.1f < p_{T} < %.1f GeV, C: %.0f-%.0f %%",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1], centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
        }
        
        // If we are drawing only the fit curves, set the scale close to one as these will be normalized to the percentege
        // of the average background level
        if(fourierV > 0){
          multiplier = firstDrawnAsymmetryBin == nAsymmetryBins ? factors[fourierV-1] : 1;
          drawer->CreateCanvas(-TMath::Pi()/2.0, TMath::Pi()*3.0/2.0, 1-vZoomTable[fourierV-1][iCentrality]*multiplier, 1+vZoomTable[fourierV-1][iCentrality]*multiplier, "#Delta#phi", "#frac{dN}{d#Delta#phi}", " ");
          
        // If we are drawing the whole fits and the histograms, set the scale using systematic uncertainty histogram
        // if the background was adjusted, otherside use the asymmetric histogram
        } else {
          if(backgroundAdjustmentUncertainty[firstDrawnAsymmetryBin][iCentrality][iTrackPt] > 0){
            adjustmentHistogram[firstDrawnAsymmetryBin][iCentrality][iTrackPt]->SetLineColor(0);
            drawer->DrawHistogram(adjustmentHistogram[firstDrawnAsymmetryBin][iCentrality][iTrackPt], "#Delta#phi", "#frac{dN}{d#Delta#phi}", " ");
          } else {
            //drawer->DrawHistogram(background[firstDrawnAsymmetryBin][iCentrality][iTrackPt], "#Delta#phi", "#frac{dN}{d#Delta#phi}", " "); // Adjusted default
            drawer->DrawHistogram(leadingUnadjustedBackground[firstDrawnAsymmetryBin][iCentrality][iTrackPt], "#Delta#phi", "#frac{dN}{d#Delta#phi}", " "); // Unadjusted default
          }
        }
        
        // If we draw the distribution, draw first the uncertainty bands from background adjustment
        if(fourierV == 0){
          for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
            if(backgroundAdjustmentUncertainty[iAsymmetry][iCentrality][iTrackPt] > 0){
              adjustmentHistogram[iAsymmetry][iCentrality][iTrackPt]->SetFillColorAlpha(colors[iAsymmetry],0.25);
              adjustmentHistogram[iAsymmetry][iCentrality][iTrackPt]->Draw("same,E2");
            }
          } // Asymmetry loop
          
        // If we draw only the fits, draw the uncertainty bands from background adjustment around the fits
        } else {
          for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
            //backgroundFit[iAsymmetry][iCentrality][iTrackPt]->SetParameter(0,1); // Normalize the fits to 1 for adjusted default
            unadjustedFit[iAsymmetry][iCentrality][iTrackPt]->SetParameter(0,1); // Normalize the fits to 1 for unadjusted default
            for(int iFitComponent = 1; iFitComponent < backgroundFit[iAsymmetry][iCentrality][iTrackPt]->GetNpar(); iFitComponent++){
              if(iFitComponent != fourierV) {
                //backgroundFit[iAsymmetry][iCentrality][iTrackPt]->SetParameter(iFitComponent,0); // Adjusted default
                unadjustedFit[iAsymmetry][iCentrality][iTrackPt]->SetParameter(iFitComponent,0); // Unadjusted default
              }
            }
            
            // uncertaintyHistogram->Eval(backgroundFit[iAsymmetry][iCentrality][iTrackPt]); // Adjusted default
            uncertaintyHistogram->Eval(unadjustedFit[iAsymmetry][iCentrality][iTrackPt]); // Unadjusted default
            uncertaintyHistogram->SetFillColorAlpha(colors[iAsymmetry],0.25);
            
            // Clone the fit function to prepare for systemtic uncertainty drawing
            uncertaintyFunctionLow = (TF1*) unadjustedFit[iAsymmetry][iCentrality][iTrackPt]->Clone(Form("lowUncertainty%d%d%d", iAsymmetry, iCentrality, iTrackPt));
            uncertaintyFunctionHigh = (TF1*) unadjustedFit[iAsymmetry][iCentrality][iTrackPt]->Clone(Form("highUncertainty%d%d%d", iAsymmetry, iCentrality, iTrackPt));
            
            // Modify the value for vn based on its uncertainty to up and down
            currentV = uncertaintyFunctionLow->GetParameter(fourierV);
            currentUncertainty = uncertaintyProvider->GetLongRangeSystematicUncertainty(fourierV-1, iCentrality, iTrackPt, iAsymmetry);
            uncertaintyFunctionLow->SetParameter(fourierV,currentV-currentUncertainty);
            uncertaintyFunctionHigh->SetParameter(fourierV,currentV+currentUncertainty);
            
            // The errors are directly fractional from the background level, so they can be used directly
            // when the curves are normalized to 1
            for(int iBin = 1; iBin <= uncertaintyHistogram->GetNbinsX(); iBin++){
              
              // In each bin, choose the larger deviation as the systematic uncertainty
              binCenter = uncertaintyHistogram->GetBinCenter(iBin);
              valueLow = TMath::Abs(uncertaintyHistogram->GetBinContent(iBin) - uncertaintyFunctionLow->Eval(binCenter));
              valueHigh = TMath::Abs(uncertaintyHistogram->GetBinContent(iBin) - uncertaintyFunctionHigh->Eval(binCenter));
              uncertaintyHistogram->SetBinError(iBin, TMath::Max(valueLow,valueHigh));
            }
            
            uncertaintyHistogram->DrawCopy("same,E4");
          }
        }
        
        // Draw all the asymmetry bins to the same figure
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawnAsymmetryBin; iAsymmetry++){
          /*// Adjusted default
          background[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(colors[iAsymmetry]);
          background[iAsymmetry][iCentrality][iTrackPt]->SetMarkerColor(colors[iAsymmetry]);
          background[iAsymmetry][iCentrality][iTrackPt]->SetMarkerStyle(markers[iAsymmetry]);
          backgroundFit[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(colors[iAsymmetry]);
           */
          
          // Unadjusted default
          leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(colors[iAsymmetry]);
          leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt]->SetMarkerColor(colors[iAsymmetry]);
          leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt]->SetMarkerStyle(markers[iAsymmetry]);
          unadjustedFit[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(colors[iAsymmetry]);
          
          // Select the wanted vn component and scale it around one
          if(fourierV > 0){
            //backgroundFit[iAsymmetry][iCentrality][iTrackPt]->Draw("same"); // Adjusted default
            unadjustedFit[iAsymmetry][iCentrality][iTrackPt]->Draw("same"); // Unadjusted default
          } else {
            //background[iAsymmetry][iCentrality][iTrackPt]->Draw("same"); // Adjusted default
            leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt]->Draw("same"); // Unadjusted default
          }
          
          // Add the histogram to legend
          legendString = Form("%.1f < x_{j} < %.1f",xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]);
          if(iAsymmetry == nAsymmetryBins) legendString = "x_{j} > 0";
          //legend->AddEntry(background[iAsymmetry][iCentrality][iTrackPt], legendString.Data(), "lp"); // Adjusted default
          legend->AddEntry(leadingUnadjustedBackground[iAsymmetry][iCentrality][iTrackPt], legendString.Data(), "lp"); // Unadjusted default
        } // asymmetry loop
        
        // Draw the legend
        legend->Draw();
        
        // Make a legend for the fourier component
        vLegend = new TLegend(vx1,0.83,vx2,0.86);
        vLegend->SetFillStyle(0);vLegend->SetBorderSize(0);vLegend->SetTextSize(0.06);vLegend->SetTextFont(62);
        vLegend->SetHeader(vHeader);
        vLegend->Draw();
        
        // Save the figures to file
        if(saveFigures){
          if(iCentrality == nCentralityBins){
            gPad->GetCanvas()->SaveAs(Form("figures/longRangeAsymmetryComparison%s%s_pp_pT=%.0f-%.0f.pdf", saveComment.Data(), vString.Data(), trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
          } else {
            gPad->GetCanvas()->SaveAs(Form("figures/longRangeAsymmetryComparison%s%s_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", saveComment.Data(), vString.Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
          }
        }
      } // Track pt Loop
    } // Centrality loop
  } // Drawing Fourier fits
  
  // Drawing graphs of extracted Fourier parameters
  if(drawFourierGraph){
    
    drawer->SetDefaultAppearanceGraph();
    
    // Arrays for track pT for graphs
    double graphPointsX[nTrackPtBins-2];      // x-axis points in flow graphs
    double graphErrorsX[nTrackPtBins-2];      // No errors for x-axis
    double graphPointsY[nTrackPtBins-2];      // Vn values
    double graphErrorsY[nTrackPtBins-2];      // Statistical errors for Vn
    double graphSystematicsY[nTrackPtBins-2]; // Systematic errors for Vn
    double graphSystematicsX[nTrackPtBins-2]; // Width of the systematic band in x axis
    
    // Helper variables for finding track pT information
    int lowPtBin, highPtBin;
    double yZoomForFlowLow[] = {-0.8,-0.005,-0.005,0.005};
    double yZoomForFlowHigh[] = {0.1,0.2,0.05,0.05};
    double legendY1, legendY2;
    
    TGraphErrors *flowGraphPt[nAsymmetryBins+1][nCentralityBins+1][nRefit];
    TGraphErrors *flowUncertaintyPt[nAsymmetryBins+1][nCentralityBins+1][nRefit];
    
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins - 2; iTrackPt++){
        
        // Find a good place to put the track pT points for the graphs
        lowPtBin = tracksForGraph[iCentrality]->FindBin(trackPtBinBorders[iTrackPt]);
        highPtBin = tracksForGraph[iCentrality]->FindBin(trackPtBinBorders[iTrackPt+1]);
        tracksForGraph[iCentrality]->GetXaxis()->SetRange(lowPtBin,highPtBin);
        graphPointsX[iTrackPt] = tracksForGraph[iCentrality]->GetMean();
        
        // Initialize other arrays to zero
        graphErrorsX[iTrackPt] = 0;
        graphPointsY[iTrackPt] = 0;
        graphErrorsY[iTrackPt] = 0;
        
        // For the systematic errors, put some width in x
        graphSystematicsX[iTrackPt] = 0.1;
      }
      
      // Create an array for the y-axis and make a graph out of vn values
      for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        for(int iFlow = 0; iFlow < nRefit; iFlow++){
          for(int iTrackPt = 0; iTrackPt < nTrackPtBins-2; iTrackPt++){
            graphPointsY[iTrackPt] = masterFlowTable[iAsymmetry][iCentrality][iTrackPt][iFlow];
            graphErrorsY[iTrackPt] = masterFlowError[iAsymmetry][iCentrality][iTrackPt][iFlow];
            graphSystematicsY[iTrackPt] = uncertaintyProvider->GetLongRangeSystematicUncertainty(iFlow, iCentrality, iTrackPt, iAsymmetry);
          }
          flowGraphPt[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsY, graphErrorsX, graphErrorsY);
          flowGraphPt[iAsymmetry][iCentrality][iFlow]->SetMarkerColor(colors[iAsymmetry]);
          flowGraphPt[iAsymmetry][iCentrality][iFlow]->SetMarkerStyle(fullMarkers[iAsymmetry]);
          
          flowUncertaintyPt[iAsymmetry][iCentrality][iFlow] = new TGraphErrors(nTrackPtBins-2, graphPointsX, graphPointsY, graphSystematicsX, graphSystematicsY);
          flowUncertaintyPt[iAsymmetry][iCentrality][iFlow]->SetFillColorAlpha(colors[iAsymmetry],0.25);
          flowUncertaintyPt[iAsymmetry][iCentrality][iFlow]->SetFillStyle(1001);
          
          // Remove two last points from the peripheral bin
          if(iCentrality == 3){
            flowGraphPt[iAsymmetry][iCentrality][iFlow]->RemovePoint(4);
            flowGraphPt[iAsymmetry][iCentrality][iFlow]->RemovePoint(3);
            flowUncertaintyPt[iAsymmetry][iCentrality][iFlow]->RemovePoint(4);
            flowUncertaintyPt[iAsymmetry][iCentrality][iFlow]->RemovePoint(3);
          }
          
        } // Loop over flow components
      } // Loop over asymmetry bins
    } // Centrality loop
    
    // Draw the graphs as a function of pT
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      for(int iFlow = 0; iFlow < 3; iFlow++){
        
        legendY1 = 0.6; legendY2 = 0.9;
        if(iFlow == 0){
          legendY1 = 0.2; legendY2 = 0.5;
        }
        
        // Setup a legend for the figure
        legend = new TLegend(0.2,legendY1,0.6,legendY2);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        if(iCentrality == nCentralityBins){
          legend->SetHeader("pp");
        } else {
          legend->SetHeader(Form("PbPb C: %.0f-%.0f %%", centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
        }
        
        // First, draw the systematic uncertainty bands to the canvas
        sprintf(namerY,"V_{%d}",iFlow+1);
        drawer->DrawGraph(flowUncertaintyPt[0][iCentrality][iFlow], 0, 8, yZoomForFlowLow[iFlow], yZoomForFlowHigh[iFlow], "Track p_{T} (GeV)", namerY, " ", "2,same");
        
        for(int iAsymmetry = 1; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
          flowUncertaintyPt[iAsymmetry][iCentrality][iFlow]->Draw("2,same");
        }
        
        // After systematic uncertainties are drawn, draw the points on top of the bands
        for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
          flowGraphPt[iAsymmetry][iCentrality][iFlow]->Draw("psame");
          if(iAsymmetry == nAsymmetryBins){
            legend->AddEntry(flowGraphPt[iAsymmetry][iCentrality][iFlow],"x_{j} > 0.0","p");
          } else {
            legend->AddEntry(flowGraphPt[iAsymmetry][iCentrality][iFlow],Form("%.1f < x_{j} < %.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]),"p");
          }
        }
        
        legend->Draw();

        // Save the figures to file
        if(saveFigures){
          if(iCentrality == nCentralityBins){
            gPad->GetCanvas()->SaveAs(Form("figures/fourierGraphPt%s_v%d_pp.pdf", saveComment.Data(), iFlow+1));
          } else {
            gPad->GetCanvas()->SaveAs(Form("figures/fourierGraphPt%s_v%d_C=%.0f-%.0f.pdf", saveComment.Data(), iFlow+1, centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
          }
        } // Saving figures
        
      } // Flow loop
    } // Centrality loop
    
  } // Draw Fourier graph
  
}

