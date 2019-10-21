#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "JDrawer.h"

/*
 * Macro for printing the relative uncertainties in latex slides
 */
void printRelativeUncertainty(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString dataFileName = "data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsWithShifterCentrality_trackDeltaRonlyInLowPt_processed_2019-10-17.root"; // Compare also with uncorrected data
  // data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root
  TString ppFileName = "data/ppData2017_highForest_pfJets_20EventsMixed_finalTrackCorr_xjBins_JECv4_wtaAxis_tunedSeagull_allCorrections_processed_2019-10-17.root";
  TString uncertaintyFileName = "uncertainties/systematicUncertaintyForPbPb_25eveMix_oldJES_15percentSpill10Jff_2019-10-17.root";
  TString ppUncertainryFileName = "uncertainties/systematicUncertaintyForPp_20percentSpillJff_2019-09-30.root";
  
  // Open the input files
  TFile *dataFile = TFile::Open(dataFileName);
  TFile *uncertaintyFile = TFile::Open(uncertaintyFileName);
  TFile *ppFile = TFile::Open(ppFileName);
  TFile *ppUncertainryFile = TFile::Open(ppUncertainryFileName);
  
  // Make readers to read the files
  DijetHistogramManager *dataManager = new DijetHistogramManager(dataFile);
  dataManager->SetLoadTrackLeadingJetCorrelationsPtWeighted(true);
  dataManager->SetLoadTrackSubleadingJetCorrelationsPtWeighted(true);
  dataManager->LoadProcessedHistograms();
  
  DijetHistogramManager *ppManager = new DijetHistogramManager(ppFile);
  ppManager->SetLoadTrackLeadingJetCorrelationsPtWeighted(true);
  ppManager->SetLoadTrackSubleadingJetCorrelationsPtWeighted(true);
  ppManager->LoadProcessedHistograms();
  
  JffCorrector *uncertaintyManager = new JffCorrector();
  uncertaintyManager->ReadSystematicFile(uncertaintyFile);
  
  JffCorrector *ppUncertaintyManager = new JffCorrector();
  ppUncertaintyManager->ReadSystematicFile(ppUncertainryFile);
  
  const int nCentralityBins = dataManager->GetNCentralityBins();
  const int nTrackPtBins = dataManager->GetNTrackPtBins();
  const int nAsymmetryBins = dataManager->GetNAsymmetryBins();
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  // Define arrays for the jet shapes
  TH1D *dataShape[2][nCentralityBins+1][nTrackPtBins+1];
  TH1D *uncertaintyShape[2][JffCorrector::knUncertaintySources][nCentralityBins+1][nTrackPtBins+1];
  
  // Define arrays for yield
  double dataYield[2][nCentralityBins+1][nTrackPtBins+1];
  double uncertaintyYield[2][JffCorrector::knUncertaintySources][nCentralityBins+1][nTrackPtBins+1];
  
  
  // Read PbPb and pp histograms from files
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt <= nTrackPtBins; iTrackPt++){
      
      // Different data reader for pp and PbPb files
      if(iCentrality == nCentralityBins){
        dataShape[0][iCentrality][iTrackPt] = ppManager->GetHistogramJetShape(DijetHistogramManager::kJetShape, DijetHistogramManager::kPtWeightedTrackLeadingJet, nAsymmetryBins, 0, iTrackPt);
        dataShape[1][iCentrality][iTrackPt] = ppManager->GetHistogramJetShape(DijetHistogramManager::kJetShape, DijetHistogramManager::kPtWeightedTrackSubleadingJet, nAsymmetryBins, 0, iTrackPt);
      } else {
        dataShape[0][iCentrality][iTrackPt] = dataManager->GetHistogramJetShape(DijetHistogramManager::kJetShape, DijetHistogramManager::kPtWeightedTrackLeadingJet, nAsymmetryBins, iCentrality, iTrackPt);
        dataShape[1][iCentrality][iTrackPt] = dataManager->GetHistogramJetShape(DijetHistogramManager::kJetShape, DijetHistogramManager::kPtWeightedTrackSubleadingJet, nAsymmetryBins, iCentrality, iTrackPt);
      }
      dataYield[0][iCentrality][iTrackPt] = dataShape[0][iCentrality][iTrackPt]->Integral(1, dataShape[0][iCentrality][iTrackPt]->FindBin(0.99), "width");
      dataYield[1][iCentrality][iTrackPt] = dataShape[1][iCentrality][iTrackPt]->Integral(1, dataShape[1][iCentrality][iTrackPt]->FindBin(0.99), "width");
      
      for(int iUncertainty = 0; iUncertainty < JffCorrector::knUncertaintySources; iUncertainty++){
        
        // Different reader for pp and PbPb uncertainty files
        if(iCentrality == nCentralityBins){
          uncertaintyShape[0][iUncertainty][iCentrality][iTrackPt] = ppUncertaintyManager->GetJetShapeSystematicUncertainty(DijetHistogramManager::kPtWeightedTrackLeadingJet, 0, iTrackPt, nAsymmetryBins, iUncertainty);
          uncertaintyShape[1][iUncertainty][iCentrality][iTrackPt] = ppUncertaintyManager->GetJetShapeSystematicUncertainty(DijetHistogramManager::kPtWeightedTrackSubleadingJet, 0, iTrackPt, nAsymmetryBins, iUncertainty);
        } else {
          uncertaintyShape[0][iUncertainty][iCentrality][iTrackPt] = uncertaintyManager->GetJetShapeSystematicUncertainty(DijetHistogramManager::kPtWeightedTrackLeadingJet, iCentrality, iTrackPt, nAsymmetryBins, iUncertainty);
          uncertaintyShape[1][iUncertainty][iCentrality][iTrackPt] = uncertaintyManager->GetJetShapeSystematicUncertainty(DijetHistogramManager::kPtWeightedTrackSubleadingJet, iCentrality, iTrackPt, nAsymmetryBins, iUncertainty);
        }
        
        uncertaintyYield[0][iUncertainty][iCentrality][iTrackPt] = uncertaintyShape[0][iUncertainty][iCentrality][iTrackPt]->Integral(1, uncertaintyShape[0][iUncertainty][iCentrality][iTrackPt]->FindBin(0.99), "width");
        uncertaintyYield[1][iUncertainty][iCentrality][iTrackPt] = uncertaintyShape[1][iUncertainty][iCentrality][iTrackPt]->Integral(1, uncertaintyShape[1][iUncertainty][iCentrality][iTrackPt]->FindBin(0.99), "width");
      } // Uncertainty type loop
    } // Track pT loop
  } // Centrality loop
  
  // Print a slide with uncertainties for each source and each centrality
  char namer[100];
  const char* jetType[2] = {"Leading","Subleading"};
  for(int iJetType = 0; iJetType < 2; iJetType++){
    for(int iTrackPt = 0; iTrackPt <= nTrackPtBins; iTrackPt++){
      
      if(iTrackPt == nTrackPtBins){
        sprintf(namer,"\\frametitle{%s jet shape error, $p_{\\mathrm{T}}$} integrated}",jetType[iJetType]);
      } else {
        sprintf(namer,"\\frametitle{%s jet shape error, $%.1f < p_{\\mathrm{T}} < %.1f$}",jetType[iJetType], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
      }
      
      cout << endl;
      cout << "\\begin{frame}" << endl;
      cout << namer << endl;
      cout << "\\begin{center}" << endl;
      cout << "  \\begin{tabular}{cccccc}" << endl;
      cout << "    \\toprule" << endl;
      cout << "    Source & C: 0-10 & C: 10-30 & C: 30-50 & C: 50-90 & pp \\\\" << endl;
      cout << "    \\midrule" << endl;
      
      // Set the correct precision for printing floating point numbers
      cout << fixed << setprecision(3);
      
      // Print one line for each track pT bin
      for(int iUncertainty = 0; iUncertainty < JffCorrector::knUncertaintySources; iUncertainty++){
        cout << uncertaintyManager->GetUncertaintyName(iUncertainty);
        for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
          cout << " & $" << uncertaintyYield[iJetType][iUncertainty][iCentrality][iTrackPt] / dataYield[iJetType][iCentrality][iTrackPt] << "$";
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
}

