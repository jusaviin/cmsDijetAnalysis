#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "DijetMethods.h"
#include "JDrawer.h"

/*
 * Macro for producing the JFF correction from PYTHIA simulation results
 */
void checkSube(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString sube0FileName = "data/PbPbMC_RecoGen_pfCsJets_onlyLeading_subeNon0_matchLeading_improvisedMixing_noDijet_wtaAxis_noSeagull_processed_2019-07-16.root";  // File from which the RecoGen histograms are read for the correction
  // "data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_sube0_xj_2019-06-10_onlyNecessarySeagull_processed.root"
  // data/PbPbMC_RecoGen_skims_pfJets_sube0_noUncorr_matchedJets_improvisedMixing_xj_processed_2019-03-18.root
  // data/PbPbMC_RecoGen_skims_pfJets_noUncorr_xj_sube0_improvisedMixing_processed_2019-03-28.root
  TString subeNon0FileName = "data/PbPbMC_RecoGen_pfCsJets_onlyLeading_subeNon0_antimatchLeading_improvisedMixing_noDijet_wtaAxis_noSeagull_processed_2019-07-16.root";   // File from which the GenGen histograms are read for the correction
  // "data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_subeNon0_xj_2019-06-06_onlyOccasionalSeagull_processed.root"
  TString anySubeFileName = "data/PbPbMC_RecoGen_pfCsJets_onlyLeading_subeNon0_matchJets_improvisedMixing_noDijet_wtaAxis_noSeagull_processed_2019-07-16.root";   // File name for the output file
  // "data/PbPbMC_RecoGen_pfCsJets_noUncorr_5eveStrictMix_xj_2019-06-12_onlyImprovedSeagull_processed.root"
  TString spilloverFileName = "newSpilloverTest_symmetrizedDistribution.root";

  // Possible distribution to compare: kSameEvent, kMixedEvent, kMixedEventNormalized, kCorrected, kBackgroundSubtracted, kBackground, kBackgroundOverlap
  int comparedDistribution = DijetHistogramManager::kCorrected;
  bool saveFigures = false;
  
  // Config done

  // Open the input files
  TFile *sube0File = TFile::Open(sube0FileName);
  TFile *subeNon0File = TFile::Open(subeNon0FileName);
  TFile *anySubeFile = TFile::Open(anySubeFileName);
  TFile *spilloverFile = TFile::Open(spilloverFileName);
  
  // Create histogram managers to provide the histograms for the correction
  DijetHistogramManager *sube0Histograms = new DijetHistogramManager(sube0File);
  sube0Histograms->SetLoadAllTrackLeadingJetCorrelations(true,false,false);
  sube0Histograms->SetLoadTrackInclusiveJetCorrelations(true);
  sube0Histograms->SetLoad2DHistograms(true);              // Two-dimensional histograms are needed for deltaEta-deltaPhi correction
  sube0Histograms->LoadProcessedHistograms();
  
  DijetHistogramManager *subeNon0Histograms = new DijetHistogramManager(subeNon0File);
  subeNon0Histograms->SetLoadAllTrackLeadingJetCorrelations(true,false,false);
  subeNon0Histograms->SetLoadTrackInclusiveJetCorrelations(true);
  subeNon0Histograms->SetLoad2DHistograms(true);              // Two-dimensional histograms are needed for deltaEta-deltaPhi correction
  subeNon0Histograms->LoadProcessedHistograms();
  
  DijetHistogramManager *anySubeHistograms = new DijetHistogramManager(anySubeFile);
  anySubeHistograms->SetLoadAllTrackLeadingJetCorrelations(true,false,false);
  anySubeHistograms->SetLoadTrackInclusiveJetCorrelations(true);
  anySubeHistograms->SetLoad2DHistograms(true);              // Two-dimensional histograms are needed for deltaEta-deltaPhi correction
  anySubeHistograms->LoadProcessedHistograms();
  
  // Find the correct number of centrality, track pT and asymmetry bins
  const int nCentralityBins = sube0Histograms->GetNCentralityBins();
  const int nTrackPtBins = sube0Histograms->GetNTrackPtBins();
  
  int iCentralityBin = 0;
  int iTrackPtBin = 1;
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12};  // Bin borders for track pT
  
  // Define arrays for the histograms
  TH2D *sube0TestHistogram[nCentralityBins][nTrackPtBins];
  TH2D *subeNon0TestHistogram[nCentralityBins][nTrackPtBins];
  TH2D *anySubeTestHistogram[nCentralityBins][nTrackPtBins];
  TH2D *spilloverHistogram[nCentralityBins][nTrackPtBins];
  TF1 *zeroLine = new TF1("Nolla","[0]",0,1.5);
  zeroLine->SetParameter(0,0);
  zeroLine->SetLineColor(kBlack);
  
  // For a quick check, read the same histogram from all three files and check that sube0 + subeNon0 = anySube
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      sube0TestHistogram[iCentrality][iTrackPt] = sube0Histograms->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackInclusiveJet,comparedDistribution,DijetHistogramManager::kMaxAsymmetryBins,iCentrality,iTrackPt);
      subeNon0TestHistogram[iCentrality][iTrackPt] = subeNon0Histograms->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackInclusiveJet,comparedDistribution,DijetHistogramManager::kMaxAsymmetryBins,iCentrality,iTrackPt);
      anySubeTestHistogram[iCentrality][iTrackPt] = anySubeHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(DijetHistogramManager::kTrackInclusiveJet,comparedDistribution,DijetHistogramManager::kMaxAsymmetryBins,iCentrality,iTrackPt);
      spilloverHistogram[iCentrality][iTrackPt] = (TH2D*) spilloverFile->Get(Form("trackLeadingJetDeltaEtaDeltaPhi/nofitSpilloverCorrection_trackLeadingJetDeltaEtaDeltaPhi_C%dT%d",iCentrality,iTrackPt));
    }
  }
    
  // Add the different subes
  TH2D *oneHistogram = (TH2D*) sube0TestHistogram[iCentralityBin][iTrackPtBin]->Clone("oneHistogram");
  double nJetsSube0 = sube0Histograms->GetInclusiveJetPtIntegral(iCentralityBin);
  oneHistogram->Scale(nJetsSube0);
  cout << "MAthced: " << nJetsSube0 << endl;
  TH2D *prutski = (TH2D*) subeNon0TestHistogram[iCentralityBin][iTrackPtBin]->Clone("twoHistogram");
  double nJetsSubeNon0 = subeNon0Histograms->GetInclusiveJetPtIntegral(iCentralityBin);
  cout << "Anti: " << nJetsSubeNon0 << endl;
  cout << "Ratio: " << nJetsSubeNon0/nJetsSube0 << endl;
  prutski->Scale(nJetsSubeNon0);
  oneHistogram->Add(prutski);
  oneHistogram->Scale(1.0/(nJetsSube0+nJetsSubeNon0));
  
  // Take the ratio to any sube
  oneHistogram->Divide(anySubeTestHistogram[iCentralityBin][iTrackPtBin]);
  
  // Look also at the ratio between subeNon0 and spillover correction
  TH2D *spilloverEvaluation = (TH2D*) subeNon0TestHistogram[iCentralityBin][iTrackPtBin]->Clone("spilloverEvaluation");
  spilloverEvaluation->Divide(spilloverHistogram[iCentralityBin][iTrackPtBin]);
  
  JDrawer *drawer = new JDrawer();
  
  // Draw the resulting histogram for inspection
  oneHistogram->GetXaxis()->SetRangeUser(-1.5,1.5);
  oneHistogram->GetYaxis()->SetRangeUser(-1.5,1.5);
  drawer->SetLeftMargin(0.15);
  drawer->SetRightMargin(0.12);
  drawer->DrawHistogram(oneHistogram,"#Delta#varphi","#Delta#eta","(Match + Antimatch)/All","colz");
  
  return;
  
  // Draw the spillover evaluation
  spilloverEvaluation->GetXaxis()->SetRangeUser(-1.5,1.5);
  spilloverEvaluation->GetYaxis()->SetRangeUser(-1.5,1.5);
  drawer->DrawHistogram(spilloverEvaluation,"#Delta#varphi","#Delta#eta","SubeNon0/Spillover","colz");
  
  // Draw the subeNon0 histogram
  subeNon0TestHistogram[iCentralityBin][iTrackPtBin]->Rebin2D(2,2);
  subeNon0TestHistogram[iCentralityBin][iTrackPtBin]->Scale(1.0/4.0);
  subeNon0TestHistogram[iCentralityBin][iTrackPtBin]->GetXaxis()->SetRangeUser(-1.5,1.5);
  subeNon0TestHistogram[iCentralityBin][iTrackPtBin]->GetYaxis()->SetRangeUser(-1.5,1.5);
  drawer->DrawHistogram(subeNon0TestHistogram[iCentralityBin][iTrackPtBin],"#Delta#varphi","#Delta#eta","SubeNon0","colz");
  
  // Draw the spillover histogram close to zero
  spilloverHistogram[iCentralityBin][iTrackPtBin]->GetXaxis()->SetRangeUser(-1.5,1.5);
  spilloverHistogram[iCentralityBin][iTrackPtBin]->GetYaxis()->SetRangeUser(-1.5,1.5);
  drawer->DrawHistogram(spilloverHistogram[iCentralityBin][iTrackPtBin],"#Delta#varphi","#Delta#eta","Spillover","colz");
  drawer->Reset();
  
  // Check also deltaR calculated from the histograms
  DijetMethods *calculator = new DijetMethods;
  TH1D *sube0DeltaR = calculator->GetJetShape(sube0TestHistogram[iCentralityBin][iTrackPtBin]);
  sube0DeltaR->SetName("subo0DeltaR");
  drawer->DrawHistogram(sube0DeltaR,"#DeltaR","P(#DeltaR)","sube0");
  TH1D *anySubeDeltaR = calculator->GetJetShape(anySubeTestHistogram[iCentralityBin][iTrackPtBin]);
  anySubeDeltaR->SetName("anySubeDeltaR");
  drawer->DrawHistogram(anySubeDeltaR,"#DeltaR","P(#DeltaR)","Any sube");
  TH1D *subeNon0DeltaRForSum = calculator->GetJetShape(subeNon0TestHistogram[iCentralityBin][iTrackPtBin]);
  subeNon0DeltaRForSum->SetName("suboNon0ForSum");
  TH1D *deltaRsum = (TH1D*) sube0DeltaR->Clone("deltaRsum");
  deltaRsum->Add(subeNon0DeltaRForSum);
  deltaRsum->SetLineColor(kRed);
  deltaRsum->Draw("same");
  
  // Make subeNon0 to spillover comparison in all bins
  TH1D *subeNon0DeltaR[nCentralityBins][nTrackPtBins];
  TH1D *spilloverDeltaR[nCentralityBins][nTrackPtBins];
  TH1D *ratioOfDeltaR[nCentralityBins][nTrackPtBins];
  TLegend *legend;
  drawer->SetDefaultAppearanceSplitCanvas();
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      drawer->CreateSplitCanvas();
      subeNon0DeltaR[iCentrality][iTrackPt] = calculator->GetJetShape(subeNon0TestHistogram[iCentrality][iTrackPt]);
      drawer->DrawHistogramToUpperPad(subeNon0DeltaR[iCentrality][iTrackPt],"#DeltaR","P(#DeltaR)","subeNon0");
      spilloverDeltaR[iCentrality][iTrackPt] = calculator->GetJetShape(spilloverHistogram[iCentrality][iTrackPt]);
      spilloverDeltaR[iCentrality][iTrackPt]->SetLineColor(kRed);
      spilloverDeltaR[iCentrality][iTrackPt]->Draw("same");
      zeroLine->DrawCopy("same");
      
      legend = new TLegend(0.5,0.55,0.9,0.8);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->SetHeader(Form("%.1f < p_{T} < %.1f, C: %.0f-%.0f %%",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1], centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
      legend->AddEntry(subeNon0DeltaR[iCentrality][iTrackPt],"SubeNon0","l");
      legend->AddEntry(spilloverDeltaR[iCentrality][iTrackPt],"Spillover","l");
      legend->Draw();
      
      ratioOfDeltaR[iCentrality][iTrackPt] = (TH1D*) spilloverDeltaR[iCentrality][iTrackPt]->Clone(Form("ratioOfDeltaR%d%d",iCentrality,iTrackPt));
      ratioOfDeltaR[iCentrality][iTrackPt]->Divide(subeNon0DeltaR[iCentrality][iTrackPt]);
      ratioOfDeltaR[iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(0.5,1.5);
      drawer->DrawHistogramToLowerPad(ratioOfDeltaR[iCentrality][iTrackPt],"#DeltaR","Corr/SubeNon0");
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/spilloverDeltaRComparison_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
      }
    }
  }
  
  // Print the configuration for each data set to the console
  DijetCard *sube0Card = new DijetCard(sube0File);
  cout << endl;
  cout << "Configuration for sube0:" << endl;
  sube0Card->Print();
  cout << endl;
  
  DijetCard *subeNon0Card = new DijetCard(subeNon0File);
  cout << "Configuration for subeNon0:" << endl;
  subeNon0Card->Print();
  cout << endl;
  
  DijetCard *anySubeCard = new DijetCard(anySubeFile);
  cout << "Configuration for anySube:" << endl;
  anySubeCard->Print();
  cout << endl;
  
}

