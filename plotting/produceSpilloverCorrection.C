#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"
#include "DijetCard.h"
#include "JDrawer.h"

/*
 * Macro for producing spillover correction for the analysis
 */
void produceSpilloverCorrection(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  bool yieldQA = false;     // Print out relative yields between uncorrected data and spillover distribution
  
  TString recoGenFileName = "data/PbPbMC_RecoGen_pfCsJets_xjBins_noUncOrInc_improvisedMixing_onlySeagull_subeNon0_wtaAxis_processed_2019-07-12.root";  // File from which the RecoGen histograms are read for the correction
  // data/PbPbMC_RecoGen_pfCsJets_xjBins_noUncOrInc_improvisedMixing_onlySeagull_subeNon0_wtaAxis_processed_2019-07-12.root
  // data/PbPbMC_RecoGen_pfCsJets_noUncorr_matchedDijets_subeNon0_improvisedMixing_onlySeagull_processed_2019-07-02.root
  // "data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_subeNon0_xjBinsIncluded_2019-06-06_onlySeagullCorrection_processed.root"
  // PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_subeNon0_xj_2019-06-06_noCorrections_processed.root
  // "data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_subeNon0_xj_2019-06-06_onlyImprovedSeagull_processed.root"
  // "data/PbPbMC_RecoGen_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_largeDeltaPhiBackground_processed_2019-02-15.root"
  // "data/PbPbMC_RecoGen_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_fixedCentality_processed_2019-02-15.root"
  // "data/PbPbMC_RecoGen_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_processed_2019-02-15.root"
  // "data/PbPbMC_RecoGen_skims_pfJets_noInclOrUncorr_10eveMixed_subeNon0_smoothedMixing_processed_2018-11-27.root"
  TString outputFileName = "spilloverTest.root"; // File name for the output file
  // data/spilloverCorrection_PbPbMC_pfCsJets_xjBins_noUncOrInc_improvisedMixing_wtaAxis_2019-07-15.root
  // data/spilloverCorrection_PbPbMC_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_smoothedMixing_refitParameters_2019-03-18.root
  // data/spilloverCorrection_PbPbMC_skims_pfJets_noInclOrUncorr_10eventsMixed_subeNon0_smoothedMixing_revisedFit_2019-02-18.root";
  TString uncorrectedDataFileName = "data/dijetPbPb_pfCsJets_xj_noUncorr_improvisedMixing_onlySeagull_processed_2019-07-05.root"; // Data file to compare yields with spillover file
  // data/PbPbMC_RecoGen_skims_pfJets_noInclUncorPtw_3eveMix_improvedMix_noJetLimit_noCorrections_processed_2019-02-09.root
  // "data/dijetPbPb_skims_pfJets_noUncorr_improvedPoolMixing_noJetLimit_noCorrections_processed_2019-01-09.root"
  
  bool regularJetTrack = true;       // Produce the correction for reguler jet-track correlations
  bool uncorrectedJetTrack = false;  // Produce the correction for uncorrected jet-track correlations
  bool ptWeightedJetTrack = true;    // Produce the correction for pT weighted jet-track correlations
  bool inclusiveJetTrack = false;     // Produce the correction for inclusive jet-track correlatio
  
  bool processAsymmetryBins = false; // Select if you want to make the correction in asymmetry bins
  
  bool correlationSelector[DijetHistogramManager::knJetTrackCorrelations] = {regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack,inclusiveJetTrack,inclusiveJetTrack};
  
  // Note: For fitMethod = 2, JffCorrector should read spilloverCorrection instead of fittedSpilloverCorrection
  int fitMethod = 2; // Method used for fitting deltaEta and deltaPhi distributions. 0 = Gauss, 1 = Gauss+constant, 2 = Gauss+Gauss
  
  int fixGaussYield = 0;  // 0 = Do not fix the Gauss yield in spillover fit. 1 = Fix the Gauss yield to deltaEta distribution without background subtraction by estimating background a a flat line fitted to the tails of the distributin. 2 = Fix the Gauss yield to the deltaEta distribution bin count of the deltaEta distribution after background subtraction.
  bool fixGaussWidth = false;  // Fix the gaussian fit width to the centrality averaged result before fitting for yield
  bool noRefitWidth = false;
  
  bool useIntegralYield = false; // True: For yield, use integral over range [-1.5,1.5], False: Use directly yield parameter from Gauss fit
  
  double minConstantFitRange[] = {1.5,1.5,1.5,1.5,1.5,1.5,1.5,1,1,1};  // Minimum range to which the constant is fitted
  double maxConstantFitRange[] = {2.5,2.5,2.5,2.5,2.5,2.5,2.5,2,2,2};  // Maximum range to which the constant is fitted
  
  // DeltaR cut applied for symmetrized spillover histograms
  //   pT bin border           0.7  1   2   3   4   8  12  300 GeV
  double spilloverRcut[4][7] = {{1.5,1.5,1.0,0.8,0.8,0.4,0.4},  // 0-10 % centrality
                                {1.5,1.5,1.5,0.8,0.8,0.4,0.4},  // 10-30 % centrality
                                {1.0,1.0,1.0,0.8,0.6,0.4,0.4},  // 30-50 % centrality
                                {0.5,0.5,0.5,0.5,0.5,0.4,0.4}}; // 50-100 % centrality
  
  // Open the input files
  TFile *recoGenFile = TFile::Open(recoGenFileName);
  TFile *yieldQAfile = TFile::Open(uncorrectedDataFileName);
  
  // Read the DijetCard from RecoGen file
  DijetCard *card = new DijetCard(recoGenFile);
  
  // Add the information about spillover parameters to the card
  card->AddOneDimensionalVector(DijetCard::kSpilloverFixGaussYield,fixGaussYield);
  card->AddOneDimensionalVector(DijetCard::kSpilloverFixGaussWidth,fixGaussWidth);
  card->AddOneDimensionalVector(DijetCard::kSpilloverFitMethod,fitMethod);
  card->AddOneDimensionalVector(DijetCard::kSpilloverUseIntegralYield,useIntegralYield);
  
  // Make an array of input files for easier initialization of histogram readers
  TFile *inputFiles[2] = {recoGenFile,yieldQAfile};
  DijetHistogramManager *histograms[2];
  
  // Create histogram managers to provide the histograms for the correction
  for(int iInputFile = 0; iInputFile < 2; iInputFile++){
    histograms[iInputFile] = new DijetHistogramManager(inputFiles[iInputFile]);
    histograms[iInputFile]->SetLoadTracks(true);  // Track needed for pT graphs
    histograms[iInputFile]->SetLoadAllTrackLeadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
    histograms[iInputFile]->SetLoadAllTrackSubleadingJetCorrelations(regularJetTrack,uncorrectedJetTrack,ptWeightedJetTrack);
    histograms[iInputFile]->SetLoadAllTrackInclusiveJetCorrelations(inclusiveJetTrack,inclusiveJetTrack);
    histograms[iInputFile]->SetLoad2DHistograms(true);  // Two-dimensional histograms are needed for deltaEta-deltaPhi correction
    if(iInputFile == 0) histograms[iInputFile]->SetAsymmetryBinRange(0,3);  // Load histograms in asymmetry bins
    histograms[iInputFile]->LoadProcessedHistograms();
  }
  
  // Find the correct number of centrality and track pT bins
  const int nCentralityBins = histograms[0]->GetNCentralityBins();
  const int nTrackPtBins = histograms[0]->GetNTrackPtBins();
  const int nAsymmetryBins = histograms[0]->GetNAsymmetryBins();
  
  // Initialize correction histograms and helper histograms
  TH2D *spilloverDeltaEtaDeltaPhi[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH2D *spilloverHelperDeltaEtaDeltaPhi[2][DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH2D *symmetrizedSpilloverDeltaEtaDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaEtaProjection[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaPhiProjection[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TF1 *spilloverDeltaEtaFit[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TF1 *spilloverDeltaPhiFit[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];

  // Create arrays to collect the information from a combination of fits.  First bin: 0 = data 1 =error
  double combinedFitYield[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins]; // Combine deltaEta and deltaPhi for one yield
  double deltaEtaCombinedWidth[2][DijetHistogramManager::knJetTrackCorrelations][nTrackPtBins]; // Combine all centrality bins for one width
  double deltaPhiCombinedWidth[2][DijetHistogramManager::knJetTrackCorrelations][nTrackPtBins]; // Combine all centrality bins for one width
  double spilloverYieldDeltaEtaBinCount[2][DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins]; // Find the yield also from bin counting
  
  // Create graphs to fit the fit parameters
  TGraphErrors *graphYield[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins];
  TGraphErrors *graphYieldBinCount[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins];
  TGraphErrors *graphDeltaEtaWidth[DijetHistogramManager::knJetTrackCorrelations];
  TGraphErrors *graphDeltaPhiWidth[DijetHistogramManager::knJetTrackCorrelations];
  TF1 *commonYieldFit[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins];
  TF1 *combinedDeltaEtaWidthFit[DijetHistogramManager::knJetTrackCorrelations];
  TF1 *combinedDeltaPhiWidthFit[DijetHistogramManager::knJetTrackCorrelations];
  
  // Final correction histogram and fits extracted from fitted graphs
  TH2D *spilloverCorrectionDeltaEtaDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TF1 *spilloverCorrectionDeltaEtaFit[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TF1 *spilloverCorrectionDeltaPhiFit[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  TF2 *spilloverCorrectionFiller = new TF2("gauss2D", "[0]/(2*TMath::Pi()*[1]*[2])*TMath::Exp(-0.5*TMath::Power(x/[1],2))*TMath::Exp(-0.5*TMath::Power(y/[2],2))",-TMath::Pi()/2,3*TMath::Pi()/2,-5,5);
  
  // Helper variables
  char namer[100];
  double commonYield;
  double deltaEtaWidth;
  double deltaPhiWidth;
  
  // x-axis information for the graphs
  TH1D *trackPtForGraphs[nCentralityBins];              // Track pT histograms to find proper place to put pT point in graphs
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double graphPointsX[nCentralityBins][nTrackPtBins];   // x-axis points in flow graphs
  double graphErrorsX[] = {0,0,0,0,0,0,0};              // No errors for x-axis
  int lowPtBin, highPtBin;                              // Helper variables to find the track pT bins
  
  // Initialize the arrays to NULL
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    
    graphDeltaEtaWidth[iJetTrack] = NULL;
    graphDeltaPhiWidth[iJetTrack] = NULL;
    
    // The widths will be fitted by first order polynomials
    sprintf(namer,"combinedDeltaEtaWidth_%s",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    combinedDeltaEtaWidthFit[iJetTrack] = new TF1(namer,"[0]+x*[1]",0,8);
    
    sprintf(namer,"combinedDeltaPhiWidth_%s",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    combinedDeltaPhiWidthFit[iJetTrack] = new TF1(namer,"[0]+x*[1]",0,8);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      graphYield[iJetTrack][iCentrality] = NULL;
      graphYieldBinCount[iJetTrack][iCentrality] = NULL;
      
      // The yields will be fitted by an exponential function
      sprintf(namer,"commonYield_%s_C%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality);
      commonYieldFit[iJetTrack][iCentrality] = new TF1(namer,"exp([0]+x*[1])",0,8);
      
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        spilloverCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt] = NULL;
        graphPointsX[iCentrality][iTrackPt] = 0;
        
        // Create blank Gaussians for the spillover deltaEta and deltaPhi fits
        sprintf(namer,"spilloverFitDeltaEta_%s_C%dT%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverCorrectionDeltaEtaFit[iJetTrack][iCentrality][iTrackPt] = new TF1(namer,"[0]/(TMath::Sqrt(2*TMath::Pi())*[1])*TMath::Exp(-0.5*TMath::Power(x/[1],2))",-1.5,1.5);
        
        sprintf(namer,"spilloverFitDeltaPhi_%s_C%dT%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverCorrectionDeltaPhiFit[iJetTrack][iCentrality][iTrackPt] = new TF1(namer,"[0]/(TMath::Sqrt(2*TMath::Pi())*[1])*TMath::Exp(-0.5*TMath::Power(x/[1],2))",-1.5,1.5);
        
        for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
          symmetrizedSpilloverDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
        }
        
        for(int iDataType = 0; iDataType < 2; iDataType++){
          
          for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
            spilloverHelperDeltaEtaDeltaPhi[iDataType][iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          }
          
          spilloverDeltaEtaDeltaPhi[iDataType][iJetTrack][iCentrality][iTrackPt] = NULL;
          spilloverDeltaEtaProjection[iDataType][iJetTrack][iCentrality][iTrackPt] = NULL;
          spilloverDeltaPhiProjection[iDataType][iJetTrack][iCentrality][iTrackPt] = NULL;
          spilloverDeltaEtaFit[iDataType][iJetTrack][iCentrality][iTrackPt] = NULL;
          spilloverDeltaPhiFit[iDataType][iJetTrack][iCentrality][iTrackPt] = NULL;
          
          combinedFitYield[iDataType][iJetTrack][iCentrality][iTrackPt] = 0;
          deltaEtaCombinedWidth[iDataType][iJetTrack][iTrackPt] = 0;
          deltaPhiCombinedWidth[iDataType][iJetTrack][iTrackPt] = 0;
          spilloverYieldDeltaEtaBinCount[iDataType][iJetTrack][iCentrality][iTrackPt] = 0;
        } // Data type loop
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation type loop
  
  // Variable definitions for yield QA
  double spilloverYield[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  double spilloverYieldDeltaEta[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  double spilloverYieldDeltaPhi[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  double spilloverYieldDeltaEtaFit[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  double spilloverYieldDeltaPhiFit[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  double dataYield[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  double yieldRatio[DijetHistogramManager::knJetTrackCorrelations][nCentralityBins][nTrackPtBins];
  int lowXbin, lowYbin, highXbin, highYbin;
  double spilloverEtaFitRange, spilloverPhiFitRange;
  double fixedSpilloverYield = 0;
  double fixedSpilloverYieldError = 0;
  int widthBins;
  
  // Create DijetMethods in which the correction procedure is implemented
  DijetMethods *corrector = new DijetMethods();
  
  // Read track pT histograms to find good places to put points in track pT axis for graphs
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    // Track pT histograms
    trackPtForGraphs[iCentrality] = histograms[1]->GetHistogramTrackPt(DijetHistogramManager::kTrack,DijetHistogramManager::kSameEvent,iCentrality);
  }
  
  // Calculate the correction
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // To obtain best fit to the spillover distribution, different bins use different fit ranges
        /*if(iTrackPt < 2){
          // Use wider fit range for the low pT bins
          spilloverEtaFitRange = 0.5;
          spilloverPhiFitRange = 0.5;
        } else {
          // Only fit close to 0 for the high pT bins
          spilloverEtaFitRange = 0.4;
          spilloverPhiFitRange = 0.4;
        }*/
        
        // If background is not subtracted, use a bit larger region to fit to get the constant background level
        // It is possible that the fit range can be further optimized, but 0.5 seems to be giving decent results
        spilloverEtaFitRange = 1.5;
        spilloverPhiFitRange = 1.5;
        
        for(int iDataType = 0; iDataType < 2; iDataType++){ // 0 = RecoGen, 1 = Uncorrected PbPb (for QA purposes)
          
          if(iDataType == 1 && iTrackPt == 6) continue; // TODO: Remove this line after this histogram is added to uncorrected file
          
          // For the RecoGen, read the bin count yield
          spilloverHelperDeltaEtaDeltaPhi[iDataType][iJetTrack][nAsymmetryBins][iCentrality][iTrackPt] = histograms[iDataType]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack,DijetHistogramManager::kCorrected, DijetHistogramManager::kMaxAsymmetryBins,iCentrality,iTrackPt);
          if(iDataType == 0){
            spilloverYieldDeltaEtaBinCount[0][iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverYield(spilloverHelperDeltaEtaDeltaPhi[iDataType][iJetTrack][nAsymmetryBins][iCentrality][iTrackPt],1.5,2.5);
            spilloverYieldDeltaEtaBinCount[1][iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverYieldError();
            
            // Scale the bin count yield yo the pT bin size
            spilloverYieldDeltaEtaBinCount[0][iJetTrack][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1]-trackPtBinBorders[iTrackPt]);
            spilloverYieldDeltaEtaBinCount[1][iJetTrack][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1]-trackPtBinBorders[iTrackPt]);
            
          }
          
          // If we are fixing the yield, find the number for that
          if(fixGaussYield == 1){
            fixedSpilloverYield = corrector->GetSpilloverYield(spilloverHelperDeltaEtaDeltaPhi[iDataType][iJetTrack][nAsymmetryBins][iCentrality][iTrackPt],minConstantFitRange[iTrackPt],maxConstantFitRange[iTrackPt]);
            fixedSpilloverYieldError = corrector->GetSpilloverYieldError();
          }
          
          // Get the signal histogram and extract the correction from it
          if(iDataType == 0 && iJetTrack < DijetHistogramManager::kTrackInclusiveJet && processAsymmetryBins){ // This can be modified if the data file also has asymmetry binning
            for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
              
              spilloverHelperDeltaEtaDeltaPhi[iDataType][iJetTrack][iAsymmetry][iCentrality][iTrackPt] = histograms[iDataType]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack,DijetHistogramManager::kBackgroundSubtracted, iAsymmetry,iCentrality,iTrackPt);
              
              // Symmetrize the subeNon0 distribution
              symmetrizedSpilloverDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = corrector->SymmetrizeHistogram(spilloverHelperDeltaEtaDeltaPhi[iDataType][iJetTrack][iAsymmetry][iCentrality][iTrackPt], spilloverRcut[iCentrality][iTrackPt]);
            }
          } else {
            spilloverHelperDeltaEtaDeltaPhi[iDataType][iJetTrack][nAsymmetryBins][iCentrality][iTrackPt] = histograms[iDataType]->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack,DijetHistogramManager::kBackgroundSubtracted, DijetHistogramManager::kMaxAsymmetryBins,iCentrality,iTrackPt);
            
            if(iDataType == 0){
              // Symmetrize the subeNon0 distribution
              symmetrizedSpilloverDeltaEtaDeltaPhi[iJetTrack][nAsymmetryBins][iCentrality][iTrackPt] = corrector->SymmetrizeHistogram(spilloverHelperDeltaEtaDeltaPhi[iDataType][iJetTrack][nAsymmetryBins][iCentrality][iTrackPt], spilloverRcut[iCentrality][iTrackPt]);
            }
          }
          
          if(fixGaussYield == 2){
            fixedSpilloverYield = corrector->GetSpilloverYield(spilloverHelperDeltaEtaDeltaPhi[iDataType][iJetTrack][nAsymmetryBins][iCentrality][iTrackPt],0,0);
            fixedSpilloverYieldError = corrector->GetSpilloverYieldError();
          }
          
          spilloverDeltaEtaDeltaPhi[iDataType][iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverCorrection(spilloverHelperDeltaEtaDeltaPhi[iDataType][iJetTrack][nAsymmetryBins][iCentrality][iTrackPt], fitMethod,  spilloverEtaFitRange, spilloverPhiFitRange, minConstantFitRange[iTrackPt], maxConstantFitRange[iTrackPt], fixedSpilloverYield);
          
          // Get the QA histograms and functions
          // Need to change name, because corrector gives the same name in both loops, which causes problems with root
          spilloverDeltaEtaProjection[iDataType][iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverDeltaEta();
          spilloverDeltaEtaProjection[iDataType][iJetTrack][iCentrality][iTrackPt]->SetName(Form("%s%d",spilloverDeltaEtaProjection[iDataType][iJetTrack][iCentrality][iTrackPt]->GetName(),iDataType));
          spilloverDeltaPhiProjection[iDataType][iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverDeltaPhi();
          spilloverDeltaPhiProjection[iDataType][iJetTrack][iCentrality][iTrackPt]->SetName(Form("%s%d",spilloverDeltaPhiProjection[iDataType][iJetTrack][iCentrality][iTrackPt]->GetName(),iDataType));
          spilloverDeltaEtaFit[iDataType][iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverDeltaEtaFit();
          spilloverDeltaEtaFit[iDataType][iJetTrack][iCentrality][iTrackPt]->SetName(Form("%s%d",spilloverDeltaEtaFit[iDataType][iJetTrack][iCentrality][iTrackPt]->GetName(),iDataType));
          if(fixGaussYield) spilloverDeltaEtaFit[iDataType][iJetTrack][iCentrality][iTrackPt]->SetParError(0,fixedSpilloverYieldError);
          
          spilloverDeltaPhiFit[iDataType][iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverDeltaPhiFit();
          spilloverDeltaPhiFit[iDataType][iJetTrack][iCentrality][iTrackPt]->SetName(Form("%s%d",spilloverDeltaPhiFit[iDataType][iJetTrack][iCentrality][iTrackPt]->GetName(),iDataType));
          if(fixGaussYield) spilloverDeltaPhiFit[iDataType][iJetTrack][iCentrality][iTrackPt]->SetParError(0,fixedSpilloverYieldError);
        } // Data type loop

        // Print out the QA numbers for yield
        if(yieldQA){
          
          // Find the bins for signal region
          lowXbin = spilloverDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(-1.5+0.001);
          highXbin = spilloverDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(1.5-0.001);
          lowYbin = spilloverHelperDeltaEtaDeltaPhi[0][iJetTrack][nAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->FindBin(-1.5+0.001);
          highYbin = spilloverHelperDeltaEtaDeltaPhi[0][iJetTrack][nAsymmetryBins][iCentrality][iTrackPt]->GetYaxis()->FindBin(1.5-0.001);
          
          // Get the yields by integration and calculate ratio
          spilloverYield[iJetTrack][iCentrality][iTrackPt] = spilloverDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->Integral(lowXbin,highXbin,lowYbin,highYbin,"width");
          dataYield[iJetTrack][iCentrality][iTrackPt] = spilloverHelperDeltaEtaDeltaPhi[0][iJetTrack][nAsymmetryBins][iCentrality][iTrackPt]->Integral(lowXbin,highXbin,lowYbin,highYbin,"width"); // Changed 1 to 0 in first index for checking purposes
          yieldRatio[iJetTrack][iCentrality][iTrackPt] = spilloverYield[iJetTrack][iCentrality][iTrackPt]/dataYield[iJetTrack][iCentrality][iTrackPt];
          
          // Find the bins for signal region in the projected histograms (might be some rebinning happening)
          lowXbin = spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(-1.5+0.001);
          highXbin = spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(1.5-0.001);
          lowYbin = spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(-1.5+0.001);
          highYbin = spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->GetXaxis()->FindBin(1.5-0.001);
          
          // For consistency, check also the integrals over deltaEta and deltaPhi projections
          spilloverYieldDeltaPhi[iJetTrack][iCentrality][iTrackPt] = spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->Integral(lowXbin,highXbin,"width");
          spilloverYieldDeltaEta[iJetTrack][iCentrality][iTrackPt] = spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->Integral(lowYbin,highYbin,"width");
          spilloverYieldDeltaPhiFit[iJetTrack][iCentrality][iTrackPt] = spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(0);
          spilloverYieldDeltaEtaFit[iJetTrack][iCentrality][iTrackPt] = spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(0);

        } // Yield QA if
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation type loop
  
  // In the end, print out some yield info is doing QA
  if(yieldQA){
    JDrawer *qaDrawer = new JDrawer();
    for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
      if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          cout << "Type: " << histograms[1]->GetJetTrackHistogramName(iJetTrack) << " Centrality: " << iCentrality << " pT " << iTrackPt <<" Spillover yield: " << spilloverYield[iJetTrack][iCentrality][iTrackPt] << " Yield without fit: " << dataYield[iJetTrack][iCentrality][iTrackPt] << " Ratio: " << yieldRatio[iJetTrack][iCentrality][iTrackPt] << endl;
          cout << "Yield from deltaPhi: " << spilloverYieldDeltaPhi[iJetTrack][iCentrality][iTrackPt] << "  Yield from deltaEta: " << spilloverYieldDeltaEta[iJetTrack][iCentrality][iTrackPt] << endl;
          cout << "Yield from deltaPhi fit: " << spilloverYieldDeltaPhiFit[iJetTrack][iCentrality][iTrackPt] << "  Yield from deltaEta fit: " << spilloverYieldDeltaEtaFit[iJetTrack][iCentrality][iTrackPt] << endl;
          cout << "Yield from deltaPhi fit integral: " << spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->Integral(-1.5,1.5) - 3*spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(2) << " Yield from deltaEta fit integral: " << spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->Integral(-1.5,1.5) - 3*spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(2) << endl;
          cout << "Background from deltaPhi fit: " << 3*spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(2) << " background from deltaEta fit: " << 3*spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(2) << endl;
          cout << endl;
          
         spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-3,3), qaDrawer->DrawHistogram(spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt],"#Delta#eta","dN/d#Delta#eta",Form("JetTrack: %d C: %d pT: %d",iJetTrack,iCentrality,iTrackPt));
          
          spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-3,3), qaDrawer->DrawHistogram(spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt],"#Delta#phi","dN/d#Delta#phi",Form("JetTrack: %d C: %d pT: %d",iJetTrack,iCentrality,iTrackPt));
          
        } // Track pT loop
      } // Centrality loop
    } // Jet-track correlation type loop
    
    // Do not run the actual spillover correction finding if we are only doing QA for the yields
    return;
    
  } // Yield QA if
  
  // Use drawer for a quick sanity check for fits
  JDrawer *drawer = new JDrawer();
  
  // After the distributions are fit, combine fit parameters and put these combined parameters to graphs for fitting
  // TODO: Use error weighted mean instead of actual mean?
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only do the correction for selected types
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // If the used data file is data/PbPbMC_RecoGen_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_fixedCentality_processed_2019-02-15.root
        // there are some bins where the fit convergence for yield is not good and they should be removed. These bins are
        // 3-4 pT bin for the centrality 0-10 %
        // 0.7-1 pT bin for centralities 10-100 %
        // Do not include these points to the averages, so that the fit to the resulting distribution will be better.
        
        // Note: This is dangerous. The values is left to 0 and may screw up the fits. Better to fill all and remove bad bins form graph.
        //if((iCentrality == 0 && iTrackPt == 3) || (iCentrality > 0 && iTrackPt == 0)){
        //  cout << "Fit convergence for yield is deemed bad in bin iCentrality: " << iCentrality << " iTrackPt: " << iTrackPt << endl;
        //  cout << "This bin will be excluded from the total yield calculation" << endl;
        //  cout << "If in some new version of the data this is not the case, please change the code accordingly!" << endl;
        //  cout << endl;
        //} else {
        
        // For the combined yield, calculate the average of points from deltaEta and deltaPhi fits
        if(useIntegralYield){
          combinedFitYield[0][iJetTrack][iCentrality][iTrackPt] = spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->Integral(-1.5,1.5)-3*spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(2); // Integral of the fit - constant background
          combinedFitYield[0][iJetTrack][iCentrality][iTrackPt] += (spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->Integral(-1.5,1.5)-3*spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(2)); // Integral of the fit - constant background
          
        } else {
          combinedFitYield[0][iJetTrack][iCentrality][iTrackPt] = spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(0);
          combinedFitYield[0][iJetTrack][iCentrality][iTrackPt] += spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(0);
        }
        
        combinedFitYield[0][iJetTrack][iCentrality][iTrackPt] = combinedFitYield[0][iJetTrack][iCentrality][iTrackPt] / 2.0;  // Divide the sum by 2 to get the average
        combinedFitYield[0][iJetTrack][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1]-trackPtBinBorders[iTrackPt]);  // Divide the sum by bin width to get yield
        
        // For the error, use the error of the Gauss yield parameter for now TODO: Check if this should be improved
        combinedFitYield[1][iJetTrack][iCentrality][iTrackPt] = TMath::Power(spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->GetParError(0),2);
        combinedFitYield[1][iJetTrack][iCentrality][iTrackPt] += TMath::Power(spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->GetParError(0),2); // Add the errors in quadrature
        combinedFitYield[1][iJetTrack][iCentrality][iTrackPt] = TMath::Sqrt(combinedFitYield[1][iJetTrack][iCentrality][iTrackPt])/2.0;
        combinedFitYield[1][iJetTrack][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1]-trackPtBinBorders[iTrackPt]);  // Divide the sum by bin width to get yield
        //}
        
        // If the used data file is data/PbPbMC_RecoGen_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_fixedCentality_processed_2019-02-15.root
        // there are some bins where the fit convergence for width is not good and they should be removed. These bins are
        // 0.7-1 pT bin for centralities 30-100 %
        // Do not include these points to the averages, so that the fit to the resulting distribution will be better.
        
        if(iCentrality > 1 && iTrackPt == 0){
          cout << "Fit convergence for width is deemed bad in bin iCentrality: " << iCentrality << " iTrackPt: " << iTrackPt << endl;
          cout << "This bin will be excluded from the total width calculation" << endl;
          cout << "If in some new version of the data this is not the case, please change the code accordingly!" << endl;
          cout << endl;
        } else {
          
          // Calculate the average width over all centrality bin and add errors in quadrature
          deltaEtaCombinedWidth[0][iJetTrack][iTrackPt] += spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(1);
          deltaEtaCombinedWidth[1][iJetTrack][iTrackPt] += TMath::Power(spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->GetParError(1),2);
          deltaPhiCombinedWidth[0][iJetTrack][iTrackPt] += spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(1);
          deltaPhiCombinedWidth[1][iJetTrack][iTrackPt] += TMath::Power(spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->GetParError(1),2);
        }
        
        // Find a good place to put the track pT points for the graphs
        lowPtBin = trackPtForGraphs[iCentrality]->FindBin(trackPtBinBorders[iTrackPt]+0.001);
        highPtBin = trackPtForGraphs[iCentrality]->FindBin(trackPtBinBorders[iTrackPt+1]-0.001);
        trackPtForGraphs[iCentrality]->GetXaxis()->SetRange(lowPtBin,highPtBin);
        graphPointsX[iCentrality][iTrackPt] = trackPtForGraphs[iCentrality]->GetMean();
        
      } // Track pT loop
    } // Centrality loop
    
    // After all the centrality bins for widths been have been added up, calculate the average for each track pT
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      widthBins = (iTrackPt == 0) ? 2 : nCentralityBins; // We have removed to centrality bins from the lowest track pT bin
      deltaEtaCombinedWidth[0][iJetTrack][iTrackPt] /= widthBins;
      deltaEtaCombinedWidth[1][iJetTrack][iTrackPt] = TMath::Sqrt(deltaEtaCombinedWidth[1][iJetTrack][iTrackPt])/widthBins;
      deltaPhiCombinedWidth[0][iJetTrack][iTrackPt] /= widthBins;
      deltaPhiCombinedWidth[1][iJetTrack][iTrackPt] = TMath::Sqrt(deltaPhiCombinedWidth[1][iJetTrack][iTrackPt])/widthBins;
    }
    
    // Variables needed in bin count fitting
    double originalValueX, originalValueY, originalErrorX, originalErrorY;
    double lowFitPt;
    
    // Now that we have the information from different fits combined, put this combined information into graphs for fitting
    graphDeltaEtaWidth[iJetTrack] = new TGraphErrors(nTrackPtBins,graphPointsX[0],deltaEtaCombinedWidth[0][iJetTrack],graphErrorsX,deltaEtaCombinedWidth[1][iJetTrack]);
    graphDeltaPhiWidth[iJetTrack] = new TGraphErrors(nTrackPtBins,graphPointsX[0],deltaPhiCombinedWidth[0][iJetTrack],graphErrorsX,deltaPhiCombinedWidth[1][iJetTrack]);
    
    // Fit the widths with a first order polynomial
    graphDeltaEtaWidth[iJetTrack]->Fit(combinedDeltaEtaWidthFit[iJetTrack],"","",0.5,8);
    graphDeltaEtaWidth[iJetTrack]->SetMarkerStyle(34);
    drawer->DrawGraph(graphDeltaEtaWidth[iJetTrack],0,8,0,1,"Track p_{T} (GeV)","Combined #Delta#eta width",histograms[0]->GetJetTrackAxisName(iJetTrack),"psame");
    graphDeltaPhiWidth[iJetTrack]->Fit(combinedDeltaPhiWidthFit[iJetTrack],"","",0.5,8);
    graphDeltaPhiWidth[iJetTrack]->SetMarkerStyle(34);
    drawer->DrawGraph(graphDeltaPhiWidth[iJetTrack],0,8,0,1,"Track p_{T} (GeV)","Combined #Delta#phi width",histograms[0]->GetJetTrackAxisName(iJetTrack),"psame");
    
    // If we want to fix width, read the width from the graph, refit with fixed width, and recombine yields for eta and phi
    if(fixGaussWidth){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          // Refit deltaEta and deltaPhi
          corrector->GetSpilloverCorrection(spilloverHelperDeltaEtaDeltaPhi[0][iJetTrack][nAsymmetryBins][iCentrality][iTrackPt], fitMethod, spilloverEtaFitRange, spilloverPhiFitRange, minConstantFitRange[iTrackPt], maxConstantFitRange[iTrackPt], 0, combinedDeltaEtaWidthFit[iJetTrack]->Eval(graphPointsX[iCentrality][iTrackPt]), combinedDeltaPhiWidthFit[iJetTrack]->Eval(graphPointsX[iCentrality][iTrackPt]));
          spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverDeltaEtaFit();
          spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt] = corrector->GetSpilloverDeltaPhiFit();
          
          // Note: this selection commented out, as we want to get directly the average for fixed width results and not do
          // fitting for yield parameters anymore. Thus the average must be calculated in all points.
          
          // If the used data file is data/PbPbMC_RecoGen_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_fixedCentality_processed_2019-02-15.root
          // there are some bins where the fit convergence for yield is not good and they should be removed. These bins are
          // 3-4 pT bin for the centrality 0-10 %
          // 0.7-1 pT bin for centralities 10-100 %
          // Do not include these points to the averages, so that the fit to the resulting distribution will be better.
          
          /*if((iCentrality == 0 && iTrackPt == 3) || (iCentrality > 0 && iTrackPt == 0)){
           cout << "Fit convergence for yield is deemed bad in bin iCentrality: " << iCentrality << " iTrackPt: " << iTrackPt << endl;
           cout << "This bin will be excluded from the total yield calculation" << endl;
           cout << "If in some new version of the data this is not the case, please change the code accordingly!" << endl;
           cout << endl;
           } else {*/
          
          // For the combined yield, calculate the average of points from deltaEta and deltaPhi fits
          if(useIntegralYield){
            combinedFitYield[0][iJetTrack][iCentrality][iTrackPt] = spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->Integral(-1.5,1.5)-3*spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(2); // Integral of the fit - constant background
            combinedFitYield[0][iJetTrack][iCentrality][iTrackPt] += (spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->Integral(-1.5,1.5)-3*spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(2)); // Integral of the fit - constant background
          } else {
            if(iCentrality == 3){ // DeltaPhi yield is more reliable in the peripheral bin, so use that for yield determination there
              combinedFitYield[0][iJetTrack][iCentrality][iTrackPt] = spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(0);
            } else {
              combinedFitYield[0][iJetTrack][iCentrality][iTrackPt] = spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(0);
            }
            combinedFitYield[0][iJetTrack][iCentrality][iTrackPt] += spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->GetParameter(0);
          }
          
          combinedFitYield[0][iJetTrack][iCentrality][iTrackPt] = combinedFitYield[0][iJetTrack][iCentrality][iTrackPt] / 2.0;  // Divide the sum by 2 to get the average
          combinedFitYield[0][iJetTrack][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1]-trackPtBinBorders[iTrackPt]);  // Divide the sum by bin width to get yield
          
          // For the error, use the error of the Gauss yield parameter for now TODO: Check if this should be improved
          combinedFitYield[1][iJetTrack][iCentrality][iTrackPt] = TMath::Power(spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->GetParError(0),2);
          combinedFitYield[1][iJetTrack][iCentrality][iTrackPt] += TMath::Power(spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->GetParError(0),2); // Add the errors in quadrature
          combinedFitYield[1][iJetTrack][iCentrality][iTrackPt] = TMath::Sqrt(combinedFitYield[1][iJetTrack][iCentrality][iTrackPt])/2.0;
          combinedFitYield[1][iJetTrack][iCentrality][iTrackPt] /= (trackPtBinBorders[iTrackPt+1]-trackPtBinBorders[iTrackPt]);  // Divide the sum by bin width to get yield
          //}
          
        } // track pT loop
      } // centrality loop
    } // Fix Gauss width
    
    // Yield graphs come in centrality bins
    drawer->SetLogY(true);
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      graphYield[iJetTrack][iCentrality] = new TGraphErrors(nTrackPtBins,graphPointsX[iCentrality],combinedFitYield[0][iJetTrack][iCentrality],graphErrorsX,combinedFitYield[1][iJetTrack][iCentrality]);
      
      // Remove the last point from each graph
      graphYield[iJetTrack][iCentrality]->RemovePoint(nTrackPtBins-1);
      
      // Fit also this graph, excluding some points from the fit
      lowFitPt = 0.5;
      if(iCentrality == 1 || iCentrality == 2) lowFitPt = 1;
      
      // Remove one point from fit for 0-10 centrality bin
      if(iCentrality == 0){
        graphYield[iJetTrack][iCentrality]->GetPoint(3,originalValueX,originalValueY);
        originalErrorY = graphYield[iJetTrack][iCentrality]->GetErrorY(3);
        originalErrorX = graphYield[iJetTrack][iCentrality]->GetErrorX(3);
        graphYield[iJetTrack][iCentrality]->RemovePoint(3);
      }
      
      // Fit the yields with an exponentinal function
      graphYield[iJetTrack][iCentrality]->Fit(commonYieldFit[iJetTrack][iCentrality],"","",lowFitPt,8);
      
      // Return the point to the graph after fit
      if(iCentrality == 0) {
        graphYield[iJetTrack][iCentrality]->SetPoint(nTrackPtBins-1,originalValueX,originalValueY);
        graphYield[iJetTrack][iCentrality]->SetPointError(nTrackPtBins-1,originalErrorX,originalErrorY);
      }
      
      graphYield[iJetTrack][iCentrality]->SetMarkerStyle(34);
      sprintf(namer,"%s C = %d",histograms[0]->GetJetTrackAxisName(iJetTrack),iCentrality);
      drawer->DrawGraph(graphYield[iJetTrack][iCentrality],0,8,0.001,10,"Track p_{T} (GeV)","Combined yield",namer,"psame");
      
      // Fill another graph with the values extracted from unsubtracred deltaEta histogram using bin counting
      graphYieldBinCount[iJetTrack][iCentrality] = new TGraphErrors(nTrackPtBins,graphPointsX[iCentrality],spilloverYieldDeltaEtaBinCount[0][iJetTrack][iCentrality],graphErrorsX,spilloverYieldDeltaEtaBinCount[1][iJetTrack][iCentrality]);
      
      // Remove one point from fit for 0-10 centrality bin
      if(iCentrality == 0){
        graphYieldBinCount[iJetTrack][iCentrality]->GetPoint(3,originalValueX,originalValueY);
        originalErrorY = graphYieldBinCount[iJetTrack][iCentrality]->GetErrorY(3);
        originalErrorX = graphYieldBinCount[iJetTrack][iCentrality]->GetErrorX(3);
        graphYieldBinCount[iJetTrack][iCentrality]->RemovePoint(3);
      }
      
      graphYieldBinCount[iJetTrack][iCentrality]->Fit("expo","","",lowFitPt,8);
      
      // Return the point to the graph after fit
      if(iCentrality == 0) {
        graphYieldBinCount[iJetTrack][iCentrality]->SetPoint(5,originalValueX,originalValueY);
        graphYieldBinCount[iJetTrack][iCentrality]->SetPointError(5,originalErrorX,originalErrorY);
      }
      
    }
    drawer->SetLogY(false);
    
    // After the graphs are fit, get the final corrections in each bin by reading the values for the fit parameters from graphs
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Find the widths for deltaEta and deltaPhi by reading them from the fits to graphs
        if(noRefitWidth){
          deltaEtaWidth = deltaEtaCombinedWidth[0][iJetTrack][iTrackPt];
          deltaPhiWidth = deltaPhiCombinedWidth[0][iJetTrack][iTrackPt];
        } else {
          deltaEtaWidth = combinedDeltaEtaWidthFit[iJetTrack]->Eval(graphPointsX[iCentrality][iTrackPt]);
          deltaPhiWidth = combinedDeltaPhiWidthFit[iJetTrack]->Eval(graphPointsX[iCentrality][iTrackPt]);
        }
        
        // Find the common yield for deltaEta and deltaPhi by evaluating it from the fit to the graph
        if(fixGaussWidth || iTrackPt > 4){
          // If width is fixed, do not fit the yield anymore
          // Also the yield does not work over pT = 8 GeV, so just use the average yield on those bins
          commonYield = combinedFitYield[0][iJetTrack][iCentrality][iTrackPt];
        } else {
          commonYield = commonYieldFit[iJetTrack][iCentrality]->Eval(graphPointsX[iCentrality][iTrackPt]);
        }
        commonYield *= (trackPtBinBorders[iTrackPt+1]-trackPtBinBorders[iTrackPt]);  // Remove the track pT bin size normalization
        
        // Now that we have yields and widths, construct the individual Gauss fits and the total spillover correction
        spilloverCorrectionDeltaEtaFit[iJetTrack][iCentrality][iTrackPt]->SetParameters(commonYield,deltaEtaWidth);
        spilloverCorrectionDeltaPhiFit[iJetTrack][iCentrality][iTrackPt]->SetParameters(commonYield,deltaPhiWidth);
        
        // Read the parameters for the two-dimensional fit
        spilloverCorrectionFiller->SetParameters(commonYield,deltaPhiWidth,deltaEtaWidth);
        
        // To get correct binning for the histogram, clone it from the previous two-dimensional histogram
        sprintf(namer,"spilloverCorrection_%s_C%dT%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt] = (TH2D*) spilloverDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->Clone(namer);
        
        // Set all the bins to zero before filling the histogram with function contents
        for(int iPhiBin = 0; iPhiBin < spilloverCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetNbinsX(); iPhiBin++){
          for(int iEtaBin = 0; iEtaBin < spilloverCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetNbinsY(); iEtaBin++){
            spilloverCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->SetBinContent(iPhiBin,iEtaBin,0);
            spilloverCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->SetBinError(iPhiBin,iEtaBin,0);
          }
        }
        
        // Fill the spillover correction histogram using the two dimensional function
        spilloverCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Eval(spilloverCorrectionFiller,"A");
        
        // Set the errors to zero after filling the histogram from the two-dimensional function
        // The systematic error estimation should be done separately, errors would be overestimated otherwise
        // By default, root just assigns the square root of the bin content as an error for each bin
        for(int iPhiBin = 1; iPhiBin <= spilloverCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetNbinsX(); iPhiBin++){
          for(int iEtaBin = 1; iEtaBin <= spilloverCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->GetNbinsY(); iEtaBin++){
            spilloverCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->SetBinError(iPhiBin,iEtaBin,0);
          }
        }
        
      } // Track pT loop
    } // Centrality loop
    
    
  } // Jet-track correlation type loop
  
  // Create the output file
  TFile *outputFile = new TFile(outputFileName,"RECREATE");
  
   // Save the obtained correction to the output file
  char histogramNamer[150];
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only write the corrections that are calculated
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaEtaDeltaPhi",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the histogram and write it into file
        sprintf(histogramNamer, "spilloverCorrection_%sDeltaEtaDeltaPhi_C%dT%d", histograms[0]->GetJetTrackHistogramName(iJetTrack), iCentrality, iTrackPt);
        spilloverDeltaEtaDeltaPhi[0][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
        // Create a new name to the fitted histogram and write it into file
        sprintf(histogramNamer,"fittedSpilloverCorrection_%sDeltaEtaDeltaPhi_C%dT%d", histograms[0]->GetJetTrackHistogramName(iJetTrack), iCentrality, iTrackPt);
        spilloverCorrectionDeltaEtaDeltaPhi[iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
        // Create a new name to the symmetrized subeNon0 histogram and write it into file
        sprintf(histogramNamer,"nofitSpilloverCorrection_%sDeltaEtaDeltaPhi_C%dT%d", histograms[0]->GetJetTrackHistogramName(iJetTrack), iCentrality, iTrackPt);
        //spilloverHelperDeltaEtaDeltaPhi[0][iJetTrack][nAsymmetryBins][iCentrality][iTrackPt]->Write(histogramNamer);
        symmetrizedSpilloverDeltaEtaDeltaPhi[iJetTrack][nAsymmetryBins][iCentrality][iTrackPt]->Write(histogramNamer);
        
        // No asymmetry binning for inclusive jets
        if(iJetTrack < DijetHistogramManager::kTrackInclusiveJet && processAsymmetryBins){
          for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
            sprintf(histogramNamer,"nofitSpilloverCorrection_%sDeltaEtaDeltaPhi_A%dC%dT%d", histograms[0]->GetJetTrackHistogramName(iJetTrack), iAsymmetry, iCentrality, iTrackPt);
            symmetrizedSpilloverDeltaEtaDeltaPhi[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Write(histogramNamer);
          }
        }
        
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Write also the card information to the correction file
    card->Write(outputFile);
    
  }
  
  // Close the output file and create QA file
  outputFile->Close();
  TString qaFileName = outputFileName.ReplaceAll(".root","_QA.root");
  TFile *qaFile = new TFile(qaFileName,"RECREATE");
 
  // Save the QA histograms and functions to file
  for(int iJetTrack = 0; iJetTrack < DijetHistogramManager::knJetTrackCorrelations; iJetTrack++){
    if(!correlationSelector[iJetTrack]) continue;  // Only write the corrections that are calculated
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaEtaProjection",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the RecoGen histogram and write it into file
        sprintf(histogramNamer,"spilloverQA_%sDeltaEtaProjection_C%dT%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaEtaProjection[0][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
        // Create a new name to the PbPb histogram and write it into file
        if(iTrackPt == 6) continue; // TODO: Remove this after this is added to reference data file
        sprintf(histogramNamer,"dataReplica_%sDeltaEtaProjection_C%dT%d",histograms[1]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaEtaProjection[1][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaPhiProjection",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the RecoGen histogram and write it into file
        sprintf(histogramNamer,"spilloverQA_%sDeltaPhiProjection_C%dT%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaPhiProjection[0][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
        // Create a new name to the PbPb histogram and write it into file
        if(iTrackPt == 6) continue; // TODO: Remove this after this is added to reference data file
        sprintf(histogramNamer,"dataReplica_%sDeltaPhiProjection_C%dT%d",histograms[1]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaPhiProjection[1][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
      }
    }
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaEtaFit",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the RecoGen histogram and write it into file
        sprintf(histogramNamer,"spilloverQA_%sDeltaEtaFit_C%dT%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaEtaFit[0][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
        // Create a new name to the spillover correction fit
        sprintf(histogramNamer,"spilloverCorrection_%sDeltaEtaFit_C%dT%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverCorrectionDeltaEtaFit[iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
        if(iTrackPt == 6) continue; // TODO: Remove this after this is added to reference data file
        // Create a new name to the PbPb histogram and write it into file
        sprintf(histogramNamer,"dataReplica_%sDeltaEtaFit_C%dT%d",histograms[1]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaEtaFit[1][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
      } // Track pT loop
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"%sDeltaPhiFit",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Create a new name to the RecoGen histogram and write it into file
        sprintf(histogramNamer,"spilloverQA_%sDeltaPhiFit_C%dT%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaPhiFit[0][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
        // Create a new name to the spillover correction fit
        sprintf(histogramNamer,"spilloverCorrection_%sDeltaPhiFit_C%dT%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverCorrectionDeltaPhiFit[iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
        if(iTrackPt == 6) continue; // TODO: Remove this after this is added to reference data file
        // Create a new name to the PbPb histogram and write it into file
        sprintf(histogramNamer,"dataReplica_%sDeltaPhiFit_C%dT%d",histograms[1]->GetJetTrackHistogramName(iJetTrack),iCentrality,iTrackPt);
        spilloverDeltaPhiFit[1][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
        
      } // Track pT loop
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the combined graphs used for fitting if it does not already exist
    sprintf(histogramNamer,"%sCombinedGraph",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    // Create a new name to centrality combined deltaEta width graph and save it to file
    sprintf(histogramNamer,"combinedGraph_%sDeltaEta",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    graphDeltaEtaWidth[iJetTrack]->Write(histogramNamer);
    
    // Create a new name to centrality combined deltaPhi width graph and save it to file
    sprintf(histogramNamer,"combinedGraph_%sDeltaPhi",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    graphDeltaPhiWidth[iJetTrack]->Write(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      // Create a new name to the deltaEta and deltaPhi combined yield graph and save it to file
      sprintf(histogramNamer,"combinedGraph_%sYield_C%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality);
      graphYield[iJetTrack][iCentrality]->Write(histogramNamer);
      
      // Create a new name to the yield graph obtained using bin counting to mixed event corrected distribution
      sprintf(histogramNamer,"binCountGraph_%sYield_C%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality);
      graphYieldBinCount[iJetTrack][iCentrality]->Write(histogramNamer);
      
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the fits to combined graphs used for fitting if it does not already exist
    sprintf(histogramNamer,"%sCombinedFit",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    // Create a new name to centrality combined deltaEta width graph and save it to file
    sprintf(histogramNamer,"combinedFit_%sDeltaEta",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    combinedDeltaEtaWidthFit[iJetTrack]->Write(histogramNamer);
    
    // Create a new name to centrality combined deltaPhi width graph and save it to file
    sprintf(histogramNamer,"combinedFit_%sDeltaPhi",histograms[0]->GetJetTrackHistogramName(iJetTrack));
    combinedDeltaPhiWidthFit[iJetTrack]->Write(histogramNamer);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      // Create a new name to the deltaEta and deltaPhi combined yield graph and save it to file
      sprintf(histogramNamer,"combinedFit_%sYield_C%d",histograms[0]->GetJetTrackHistogramName(iJetTrack),iCentrality);
      commonYieldFit[iJetTrack][iCentrality]->Write(histogramNamer);
      
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // In the very end, write also the card information to the correction QA file
    card->Write(qaFile);
    
  } // Jet-track loop
}
