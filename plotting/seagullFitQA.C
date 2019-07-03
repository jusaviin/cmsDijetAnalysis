#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"
#include "DijetCard.h"
#include "JffCorrector.h"
#include "JDrawer.h"

/*
 * Combination of zeroth and first order polynomial for seagull correction
 * Functional form: f(x) = c, if |x| < 0.5  <=> f(x) = a*x^2+b*x+c if |x| > 0.5
 */
double seagullPoly2(double *x, double *par){
  if(x[0] < 0.5 && x[0] > -0.5) return par[0];
  return par[0]+par[1]*x[0]*x[0]+par[2]*x[0];
}

/*
 * Combination of zeroth and first order polynomial for seagull correction
 * Functional form: f(x) = a + b*e^(c*x)
 */
double seagullExp(double *x, double *par){
  return par[0]+par[1]*TMath::Exp(x[0]*par[2]);
  //return par[0]+x[0]*par[1]+x[0]*x[0]*par[2]+x[0]*x[0]*x[3]*par[3];
}


/*
 * Macro for some QA to determine automatically when seagull fit should be applied and when not.
 */
void seagullFitQA(){
  
  // Define the file used for fit testing
  TString recoGenFileName = "data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_xj_2019-06-12_noCorrections_processed.root";
  // "data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_subeNon0_xj_2019-06-06_noCorrections_processed.root"
  TFile *recoGenFile = TFile::Open(recoGenFileName);
  
  // Load the histograms from the file
  DijetHistogramManager *histograms;
  histograms = new DijetHistogramManager(recoGenFile);
  histograms->SetLoadAllTrackLeadingJetCorrelations(true,false,false);
  histograms->SetLoad2DHistograms(true);  // Two-dimensional histograms are needed for deltaEta-deltaPhi correction
  histograms->LoadProcessedHistograms();
  
  // Find the correct number of centrality and track pT bins
  const int nCentralityBins = histograms->GetNCentralityBins();
  const int nTrackPtBins = histograms->GetNTrackPtBins();
  
  // Rebin before fitting
  int rebinBeforeFit = 4;
  
  // Initialize correction histograms and helper histograms
  TH1D *seagullDeltaEta[nCentralityBins][nTrackPtBins];
  TF1 *seagullDeltaEtaConstantFit[nCentralityBins][nTrackPtBins];
  TF1 *seagullDeltaEtaPolyFit[nCentralityBins][nTrackPtBins];
  
  // Create DijetMethods in which the correction procedure is implemented
  DijetMethods *corrector = new DijetMethods();
  double initialLevel;
  
  // Read the histograms from the file
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      
      // Read the deltaEta distribution to which seagull fit is made form the file
      seagullDeltaEta[iCentrality][iTrackPt] = histograms->GetHistogramJetTrackDeltaEta(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kCorrected, DijetHistogramManager::kMaxAsymmetryBins, iCentrality, iTrackPt, DijetHistogramManager::kBetweenPeaks);
      seagullDeltaEta[iCentrality][iTrackPt]->Rebin(rebinBeforeFit);
      seagullDeltaEta[iCentrality][iTrackPt]->Scale(1.0/rebinBeforeFit);
      
      seagullDeltaEtaPolyFit[iCentrality][iTrackPt] = new TF1(Form("seagullPolyFit%d%d",iCentrality,iTrackPt),seagullPoly2,-3,3,3);
      initialLevel = seagullDeltaEta[iCentrality][iTrackPt]->GetBinContent(seagullDeltaEta[iCentrality][iTrackPt]->FindBin(0));
      seagullDeltaEtaPolyFit[iCentrality][iTrackPt]->SetParameters(initialLevel,0,0);
      
      // Fit the distribution with a constant
      seagullDeltaEta[iCentrality][iTrackPt]->Fit("pol0","","",-3,3);
      
      // Extort the fit function out from the histogram
      seagullDeltaEtaConstantFit[iCentrality][iTrackPt] = seagullDeltaEta[iCentrality][iTrackPt]->GetFunction("pol0");
      
      // Rename the fit function to avoit root problems
      seagullDeltaEtaConstantFit[iCentrality][iTrackPt]->SetName(Form("constantFit%d%d",iCentrality,iTrackPt));
      
      // Remove the function form the list of functions of the histogram
      seagullDeltaEta[iCentrality][iTrackPt]->RecursiveRemove(seagullDeltaEtaConstantFit[iCentrality][iTrackPt]);
      
      // Fit the function again with second order polynomial
      seagullDeltaEta[iCentrality][iTrackPt]->Fit(seagullDeltaEtaPolyFit[iCentrality][iTrackPt],"","",-3,3);
      
    } // Track pT loop
  } // Centrality loop
  
  // Define a JDrawer
  JDrawer *drawer = new JDrawer();

  // Draw the histograms from the file
  double constantChi2, polyChi2;
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    cout << endl;
    cout << "Results for centrality " << iCentrality << endl;
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      
      seagullDeltaEta[iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-3,3);
      drawer->DrawHistogram(seagullDeltaEta[iCentrality][iTrackPt],"#Delta#eta","counts",Form("C = %d, pT = %d",iCentrality,iTrackPt));

      // Print the chi2 / NDF for each bin
      constantChi2 = seagullDeltaEtaConstantFit[iCentrality][iTrackPt]->GetChisquare()/seagullDeltaEtaConstantFit[iCentrality][iTrackPt]->GetNDF();
      polyChi2 = seagullDeltaEtaPolyFit[iCentrality][iTrackPt]->GetChisquare()/seagullDeltaEtaPolyFit[iCentrality][iTrackPt]->GetNDF();
      cout << "Track pT " << iTrackPt << ": constant chi2/ndf = " << constantChi2 << " poly2 chi2/ndf = " << polyChi2 << " ratio: " << constantChi2/polyChi2 <<  endl;
      
    } // Track pT loop
  } // Centrality loop
  
}
