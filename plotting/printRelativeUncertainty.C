#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "JDrawer.h"

/*
 * Macro for producing the JFF correction from PYTHIA simulation results
 */
void printRelativeUncertainty(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString dataFileName = "data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root"; // Compare also with uncorrected data
  TString uncertaintyFileName = "systematicTest.root";
  
  // Open the input files
  TFile *dataFile = TFile::Open(dataFileName);
  TFile *uncertaintyFile = TFile::Open(uncertaintyFileName);
  
  // Make readers to read the files
  DijetHistogramManager *dataManager = new DijetHistogramManager(dataFile);
  dataManager->SetLoadTrackLeadingJetCorrelationsPtWeighted(true);
  dataManager->SetLoadTrackSubleadingJetCorrelationsPtWeighted(true);
  dataManager->LoadProcessedHistograms();
  
  JffCorrector *uncertaintyManager = new JffCorrector();
  uncertaintyManager->ReadSystematicFile(uncertaintyFile);
  
  const int nCentralityBins = dataManager->GetNCentralityBins();
  const int nTrackPtBins = dataManager->GetNTrackPtBins();
  const int nAsymmetryBins = dataManager->GetNAsymmetryBins();
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  // Define arrays for the jet shapes
  TH1D *dataShape[2][nCentralityBins][nTrackPtBins];
  TH1D *uncertaintyShape[2][JffCorrector::knUncertaintySources][nCentralityBins][nTrackPtBins];
  
  // Define arrays for yield
  double dataYield[2][nCentralityBins][nTrackPtBins];
  double uncertaintyYield[2][JffCorrector::knUncertaintySources][nCentralityBins][nTrackPtBins];
  
  
  // Read the histograms from spillover file
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      dataShape[0][iCentrality][iTrackPt] = dataManager->GetHistogramJetShape(DijetHistogramManager::kJetShape, DijetHistogramManager::kPtWeightedTrackLeadingJet, nAsymmetryBins, iCentrality, iTrackPt);
      dataYield[0][iCentrality][iTrackPt] = dataShape[0][iCentrality][iTrackPt]->Integral(1, dataShape[0][iCentrality][iTrackPt]->FindBin(0.99), "width");
      dataShape[1][iCentrality][iTrackPt] = dataManager->GetHistogramJetShape(DijetHistogramManager::kJetShape, DijetHistogramManager::kPtWeightedTrackSubleadingJet, nAsymmetryBins, iCentrality, iTrackPt);
      dataYield[1][iCentrality][iTrackPt] = dataShape[1][iCentrality][iTrackPt]->Integral(1, dataShape[1][iCentrality][iTrackPt]->FindBin(0.99), "width");
      for(int iUncertainty = 0; iUncertainty < JffCorrector::knUncertaintySources; iUncertainty++){
        uncertaintyShape[0][iUncertainty][iCentrality][iTrackPt] = uncertaintyManager->GetJetShapeSystematicUncertainty(DijetHistogramManager::kPtWeightedTrackLeadingJet, iCentrality, iTrackPt, nAsymmetryBins, iUncertainty);
        uncertaintyYield[0][iUncertainty][iCentrality][iTrackPt] = uncertaintyShape[0][iUncertainty][iCentrality][iTrackPt]->Integral(1, uncertaintyShape[0][iUncertainty][iCentrality][iTrackPt]->FindBin(0.99), "width");
        uncertaintyShape[1][iUncertainty][iCentrality][iTrackPt] = uncertaintyManager->GetJetShapeSystematicUncertainty(DijetHistogramManager::kPtWeightedTrackSubleadingJet, iCentrality, iTrackPt, nAsymmetryBins, iUncertainty);
        uncertaintyYield[1][iUncertainty][iCentrality][iTrackPt] = uncertaintyShape[1][iUncertainty][iCentrality][iTrackPt]->Integral(1, uncertaintyShape[1][iUncertainty][iCentrality][iTrackPt]->FindBin(0.99), "width");
      }
    }
  }
  
  // Print a slide with uncertainties for each source and each centrality
  char namer[100];
  for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
    
    sprintf(namer,"\\frametitle{Subleading jet shape error, $%.1f < p_{\\mathrm{T}} < %.1f$}", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
    
    cout << endl;
    cout << "\\begin{frame}" << endl;
    cout << namer << endl;
    cout << "\\begin{center}" << endl;
    cout << "  \\begin{tabular}{ccccc}" << endl;
    cout << "    \\toprule" << endl;
    cout << "    Source & C: 0-10 & C: 10-30 & C: 30-50 & C: 50-100 \\\\" << endl;
    cout << "    \\midrule" << endl;
    
    // Set the correct precision for printing floating point numbers
    cout << fixed << setprecision(3);
    
    // Print one line for each track pT bin
    for(int iUncertainty = 0; iUncertainty < JffCorrector::knUncertaintySources; iUncertainty++){
      cout << uncertaintyManager->GetUncertaintyName(iUncertainty);
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        cout << " & $" << uncertaintyYield[1][iUncertainty][iCentrality][iTrackPt] / dataYield[1][iCentrality][iTrackPt] << "$";
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

