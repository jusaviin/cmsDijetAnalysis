#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"
#include "DijetCard.h"
#include "JffCorrector.h"
#include "JDrawer.h"

/*
 * Macro for testing different fits to spillover distribution
 */
void fitTest(){
  
  // Configuration
  bool doRebin = false;  // Rebin the histograms when checking if things are fine
  bool drawFitComparison = false; // Draw comparison between double Gauss and single Gauss fits
  bool drawTightProjection = true;  // Draw deltaEta and deltaPhi projection only from the top of the jet peak
  
  // Define the file used for fit testing
  TString recoGenFileName = "data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_subeNon0_xj_2019-06-06_onlyImprovedSeagull_processed.root";
  TFile *recoGenFile = TFile::Open(recoGenFileName);
  TFile *spilloverFile = TFile::Open("newSpilloverTest_doubleGauss.root");
  
  // Load the histograms from the file
  DijetHistogramManager *histograms;
  histograms = new DijetHistogramManager(recoGenFile);
  histograms->SetLoadAllTrackLeadingJetCorrelations(true,false,false);
  histograms->SetLoad2DHistograms(true);  // Two-dimensional histograms are needed for deltaEta-deltaPhi correction
  histograms->LoadProcessedHistograms();
  
  JffCorrector *spilloverProvider = new JffCorrector();
  spilloverProvider->ReadSpilloverFile(spilloverFile);
  
  // Find the correct number of centrality and track pT bins
  const int nCentralityBins = histograms->GetNCentralityBins();
  const int nTrackPtBins = histograms->GetNTrackPtBins();
  
  // Initialize correction histograms and helper histograms
  TH2D *spilloverCorrectionDeltaEtaDeltaPhi[nCentralityBins][nTrackPtBins];
  TH2D *spilloverDeltaEtaDeltaPhi[nCentralityBins][nTrackPtBins];
  TH2D *spilloverHelperDeltaEtaDeltaPhi[nCentralityBins][nTrackPtBins];
  TH2D *spilloverDeltaEtaDeltaPhiDoubleGauss[nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaEtaProjection[nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaPhiProjection[nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaEtaProjectionDoubleGauss[nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaPhiProjectionDoubleGauss[nCentralityBins][nTrackPtBins];
  TH1D *tightProjectionDeltaEta[nCentralityBins][nTrackPtBins];
  TH1D *tightProjectionDeltaPhi[nCentralityBins][nTrackPtBins];
  TH1D *tightCorrectionProjectionDeltaEta[nCentralityBins][nTrackPtBins];
  TH1D *tightCorrectionProjectionDeltaPhi[nCentralityBins][nTrackPtBins];
  TF1 *spilloverDeltaEtaFit[nCentralityBins][nTrackPtBins];
  TF1 *spilloverDeltaPhiFit[nCentralityBins][nTrackPtBins];
  TF1 *spilloverDeltaEtaFitDoubleGauss[nCentralityBins][nTrackPtBins];
  TF1 *spilloverDeltaPhiFitDoubleGauss[nCentralityBins][nTrackPtBins];
  
  // Create DijetMethods in which the correction procedure is implemented
  DijetMethods *corrector = new DijetMethods();
  
  // Define fit ranges
  double spilloverEtaFitRange = 1.5;
  double spilloverPhiFitRange = 1.5;
  
  // Read the histograms from the file
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins-1; iTrackPt++){
      
      // Read the two dimensional subeNon0 distribution
      spilloverHelperDeltaEtaDeltaPhi[iCentrality][iTrackPt] = histograms->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackgroundSubtracted, DijetHistogramManager::kMaxAsymmetryBins, iCentrality, iTrackPt);
      
      // Project very central deltaEta and deltaPhi distributions
      tightProjectionDeltaEta[iCentrality][iTrackPt] = corrector->ProjectRegionDeltaEta(spilloverHelperDeltaEtaDeltaPhi[iCentrality][iTrackPt], -0.4, 0.4, Form("tightDeltaEtaProjection%d%d",iCentrality,iTrackPt));
      tightProjectionDeltaPhi[iCentrality][iTrackPt] = corrector->ProjectRegionDeltaPhi(spilloverHelperDeltaEtaDeltaPhi[iCentrality][iTrackPt], -0.4, 0.4, Form("tightDeltaPhiProjection%d%d",iCentrality,iTrackPt));
      
      // Read the actual spillover correction
      spilloverCorrectionDeltaEtaDeltaPhi[iCentrality][iTrackPt] = spilloverProvider->GetDeltaEtaDeltaPhiSpilloverCorrection(DijetHistogramManager::kTrackLeadingJet, iCentrality, iTrackPt);
      
      // Project very central deltaEta and deltaPhi distributions
      tightCorrectionProjectionDeltaEta[iCentrality][iTrackPt] = corrector->ProjectRegionDeltaEta(spilloverCorrectionDeltaEtaDeltaPhi[iCentrality][iTrackPt], -0.4, 0.4, Form("tightDeltaEtaCorrectionProjection%d%d",iCentrality,iTrackPt));
      tightCorrectionProjectionDeltaPhi[iCentrality][iTrackPt] = corrector->ProjectRegionDeltaPhi(spilloverCorrectionDeltaEtaDeltaPhi[iCentrality][iTrackPt], -0.4, 0.4, Form("tightDeltaPhiCorrectionProjection%d%d",iCentrality,iTrackPt));
      
      // Fit projections with Gaus
      //tightProjectionDeltaEta[iCentrality][iTrackPt]->Fit("gaus","","",-0.4,0.4);
      //tightProjectionDeltaPhi[iCentrality][iTrackPt]->Fit("gaus","","",-0.4,0.4);
      tightCorrectionProjectionDeltaEta[iCentrality][iTrackPt]->SetLineColor(kRed);
      tightCorrectionProjectionDeltaPhi[iCentrality][iTrackPt]->SetLineColor(kRed);
      
      // Do the fitting for deltaEta and deltaPhi using double Gauss fit
      spilloverDeltaEtaDeltaPhiDoubleGauss[iCentrality][iTrackPt] = corrector->GetSpilloverCorrection(spilloverHelperDeltaEtaDeltaPhi[iCentrality][iTrackPt], 2, spilloverEtaFitRange, spilloverPhiFitRange);
      
      // Get the fitted histograms
      spilloverDeltaEtaProjectionDoubleGauss[iCentrality][iTrackPt] = corrector->GetSpilloverDeltaEta();
      spilloverDeltaEtaProjectionDoubleGauss[iCentrality][iTrackPt]->SetName(Form("deltaEtaProjectionDoubleGauss%d%d",iCentrality,iTrackPt));
      spilloverDeltaPhiProjectionDoubleGauss[iCentrality][iTrackPt] = corrector->GetSpilloverDeltaPhi();
      spilloverDeltaPhiProjectionDoubleGauss[iCentrality][iTrackPt]->SetName(Form("deltaPhiProjectionDoubleGauss%d%d",iCentrality,iTrackPt));
      
      // Get the fit functions from double Gauss fit
      spilloverDeltaEtaFitDoubleGauss[iCentrality][iTrackPt] = corrector->GetSpilloverDeltaEtaFit();
      spilloverDeltaEtaFitDoubleGauss[iCentrality][iTrackPt]->SetName(Form("deltaEtaFitDouble%d%d",iCentrality,iTrackPt));
      spilloverDeltaPhiFitDoubleGauss[iCentrality][iTrackPt] = corrector->GetSpilloverDeltaPhiFit();
      spilloverDeltaPhiFitDoubleGauss[iCentrality][iTrackPt]->SetName(Form("deltaPhiFitDouble%d%d",iCentrality,iTrackPt));
      
      // Do the fitting for deltaEta and deltaPhi using single Gauss fit
      spilloverDeltaEtaDeltaPhi[iCentrality][iTrackPt] = corrector->GetSpilloverCorrection(spilloverHelperDeltaEtaDeltaPhi[iCentrality][iTrackPt], 0, spilloverEtaFitRange, spilloverPhiFitRange);
      
      // Get the fit functions from single Gauss fit
      spilloverDeltaEtaFit[iCentrality][iTrackPt] = corrector->GetSpilloverDeltaEtaFit();
      spilloverDeltaEtaFit[iCentrality][iTrackPt]->SetName(Form("deltaEtaFitIsThis%d%d",iCentrality,iTrackPt));
      spilloverDeltaPhiFit[iCentrality][iTrackPt] = corrector->GetSpilloverDeltaPhiFit();
      spilloverDeltaPhiFit[iCentrality][iTrackPt]->SetName(Form("deltaPhiFitIsThis%d%d",iCentrality,iTrackPt));
      
    } // Track pT loop
  } // Centrality loop
  
  // Define a JDrawer
  JDrawer *drawer = new JDrawer();
  double yieldParameter;
  
  TH2D *symmetrizedHistogram = corrector->SymmetrizeHistogram(spilloverHelperDeltaEtaDeltaPhi[0][1],1.5,1.5);
  drawer->SetRightMargin(0.12);
  spilloverHelperDeltaEtaDeltaPhi[0][1]->GetXaxis()->SetRangeUser(-1.5,1.5);
  spilloverHelperDeltaEtaDeltaPhi[0][1]->GetYaxis()->SetRangeUser(-1.5,1.5);
  symmetrizedHistogram->GetXaxis()->SetRangeUser(-1.5,1.5);
  symmetrizedHistogram->GetYaxis()->SetRangeUser(-1.5,1.5);
  drawer->DrawHistogram(spilloverHelperDeltaEtaDeltaPhi[0][1],"#Delta#phi","#Delta#eta","subeNon0, 0 < C < 10 %, 1 < p_{T} < 2 GeV","colz");
  drawer->DrawHistogram(symmetrizedHistogram,"#Delta#phi","#Delta#eta","Symmetrized, 0 < C < 10 %, 1 < p_{T} < 2 GeV","colz");
  
  return;
  
  // Read the histograms from the file
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins-1; iTrackPt++){
      
      // Draw tight projections
      if(drawTightProjection){
        if(doRebin){
          tightProjectionDeltaEta[iCentrality][iTrackPt]->Rebin(4);
          tightProjectionDeltaEta[iCentrality][iTrackPt]->Scale(1.0/4);
        }
        tightProjectionDeltaEta[iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-1.5,1.5);
        drawer->DrawHistogram(tightProjectionDeltaEta[iCentrality][iTrackPt],"#Delta#eta","counts",Form("C = %d, pT = %d",iCentrality,iTrackPt));
        tightCorrectionProjectionDeltaEta[iCentrality][iTrackPt]->Draw("same");
        
        if(doRebin){
          tightProjectionDeltaPhi[iCentrality][iTrackPt]->Rebin(2);
          tightProjectionDeltaPhi[iCentrality][iTrackPt]->Scale(1.0/2);
        }
        tightProjectionDeltaPhi[iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-1.5,1.5);
        drawer->DrawHistogram(tightProjectionDeltaPhi[iCentrality][iTrackPt],"#Delta#phi","counts",Form("C = %d, pT = %d",iCentrality,iTrackPt));
        tightCorrectionProjectionDeltaPhi[iCentrality][iTrackPt]->Draw("same");
      }
      
      if(drawFitComparison){
        
        // Draw the fitted distributions to canvas
        if(doRebin){
          spilloverDeltaEtaProjectionDoubleGauss[iCentrality][iTrackPt]->Rebin(4);
          spilloverDeltaEtaProjectionDoubleGauss[iCentrality][iTrackPt]->Scale(1.0/4);
        }
        spilloverDeltaEtaProjectionDoubleGauss[iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-1.5,1.5);
        drawer->DrawHistogram(spilloverDeltaEtaProjectionDoubleGauss[iCentrality][iTrackPt],"#Delta#eta","counts",Form("C = %d, pT = %d",iCentrality,iTrackPt));
        spilloverDeltaEtaFit[iCentrality][iTrackPt]->SetLineColor(kBlue);
        spilloverDeltaEtaFit[iCentrality][iTrackPt]->Draw("same");
        
        // Draw the different components from double fit
        /*yieldParameter = spilloverDeltaEtaFitDoubleGauss[iCentrality][iTrackPt]->GetParameter(2);
         spilloverDeltaEtaFitDoubleGauss[iCentrality][iTrackPt]->SetParameter(2,0);
         spilloverDeltaEtaFitDoubleGauss[iCentrality][iTrackPt]->SetLineColor(kMagenta);
         spilloverDeltaEtaFitDoubleGauss[iCentrality][iTrackPt]->DrawCopy("same");
         spilloverDeltaEtaFitDoubleGauss[iCentrality][iTrackPt]->SetParameter(2,yieldParameter);
         yieldParameter = spilloverDeltaEtaFitDoubleGauss[iCentrality][iTrackPt]->GetParameter(0);
         spilloverDeltaEtaFitDoubleGauss[iCentrality][iTrackPt]->SetParameter(0,0);
         spilloverDeltaEtaFitDoubleGauss[iCentrality][iTrackPt]->SetLineColor(kCyan);
         spilloverDeltaEtaFitDoubleGauss[iCentrality][iTrackPt]->DrawCopy("same");
         spilloverDeltaEtaFitDoubleGauss[iCentrality][iTrackPt]->SetParameter(0,yieldParameter);*/
        
        if(doRebin){
          spilloverDeltaPhiProjectionDoubleGauss[iCentrality][iTrackPt]->Rebin(2);
          spilloverDeltaPhiProjectionDoubleGauss[iCentrality][iTrackPt]->Scale(1.0/2);
        }
        spilloverDeltaPhiProjectionDoubleGauss[iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-1.5,1.5);
        drawer->DrawHistogram(spilloverDeltaPhiProjectionDoubleGauss[iCentrality][iTrackPt],"#Delta#phi","counts",Form("C = %d, pT = %d",iCentrality,iTrackPt));
        spilloverDeltaPhiFit[iCentrality][iTrackPt]->SetLineColor(kBlue);
        spilloverDeltaPhiFit[iCentrality][iTrackPt]->Draw("same");
        
        // Draw the different components from double fit
        /*yieldParameter = spilloverDeltaPhiFitDoubleGauss[iCentrality][iTrackPt]->GetParameter(2);
         spilloverDeltaPhiFitDoubleGauss[iCentrality][iTrackPt]->SetParameter(2,0);
         spilloverDeltaPhiFitDoubleGauss[iCentrality][iTrackPt]->SetLineColor(kMagenta);
         spilloverDeltaPhiFitDoubleGauss[iCentrality][iTrackPt]->DrawCopy("same");
         spilloverDeltaPhiFitDoubleGauss[iCentrality][iTrackPt]->SetParameter(2,yieldParameter);
         yieldParameter = spilloverDeltaPhiFitDoubleGauss[iCentrality][iTrackPt]->GetParameter(0);
         spilloverDeltaPhiFitDoubleGauss[iCentrality][iTrackPt]->SetParameter(0,0);
         spilloverDeltaPhiFitDoubleGauss[iCentrality][iTrackPt]->SetLineColor(kCyan);
         spilloverDeltaPhiFitDoubleGauss[iCentrality][iTrackPt]->DrawCopy("same");
         spilloverDeltaPhiFitDoubleGauss[iCentrality][iTrackPt]->SetParameter(0,yieldParameter);*/
        
      }
    } // Track pT loop
  } // Centrality loop
  
}
