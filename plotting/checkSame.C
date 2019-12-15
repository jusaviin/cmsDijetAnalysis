#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "DijetMethods.h"
#include "JDrawer.h"

/*
 * Macro for producing the JFF correction from PYTHIA simulation results
 */
void checkSame(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString fileName1 = "data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_xjBins_sube0_wtaAxis_onlyJffCorrection_jffExplore_processed_JECv6_2019-11-26.root";  // File from which the RecoGen histograms are read for the correction
  // PbPbMC_RecoGen_skims_pfJets_noIncOrUnc_5eveStrictMix_matchedDijets_subeNon0_noCorrections_processed_2019-05-08.root
  // PbPbMC_RecoGen_skims_pfJets_noUncOrInc_5eveImprovedMix_subeNon0_2019-02-15_noCorrections_processed.root
  // data/PbPbMC_RecoGen_skims_pfJets_noUncorr_xj_sube0_improvisedMixing_processed_2019-03-28.root
  TString fileName2 = "data/PbPbMC_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_sube0_wtaAxis_jet100trigger_JECv6_processed_2019-09-26.root";   // File from which the GenGen histograms are read for the correction

  // Possible distribution to compare: kSameEvent, kMixedEvent, kMixedEventNormalized, kCorrected, kBackgroundSubtracted, kBackground, kBackgroundOverlap
  int comparedDistribution = DijetHistogramManager::kBackgroundSubtracted;
  bool saveFigures = false;
  
  // Config done

  // Open the input files
  TFile *file1 = TFile::Open(fileName1);
  TFile *file2 = TFile::Open(fileName2);
  
  // Create histogram managers to provide the histograms for the correction
  DijetHistogramManager *firstHistograms = new DijetHistogramManager(file1);
  firstHistograms->SetLoadAllTrackLeadingJetCorrelations(true,false,false);
  firstHistograms->SetLoad2DHistograms(true);              // Two-dimensional histograms are needed for deltaEta-deltaPhi correction
  firstHistograms->LoadProcessedHistograms();
  
  DijetHistogramManager *secondHistograms = new DijetHistogramManager(file2);
  secondHistograms->SetLoadAllTrackLeadingJetCorrelations(true,false,false);
  secondHistograms->SetLoad2DHistograms(true);              // Two-dimensional histograms are needed for deltaEta-deltaPhi correction
  secondHistograms->LoadProcessedHistograms();
  
  // Find the correct number of centrality, track pT and asymmetry bins
  const int nCentralityBins = firstHistograms->GetNCentralityBins();
  const int nTrackPtBins = firstHistograms->GetNTrackPtBins();
  
  int iCentralityBin = 0;
  int iTrackPtBin = 1;
  double value250 = 0;
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,300};  // Bin borders for track pT
  
  // Define arrays for the histograms
  TH2D *firstTestHistogram[nCentralityBins][nTrackPtBins];
  TH2D *secondTestHistogram[nCentralityBins][nTrackPtBins];
  TH2D *testRatioHistogram[nCentralityBins][nTrackPtBins];
  TH1D *firstTestProjection[nCentralityBins][nTrackPtBins];
  TH1D *secondTestProjection[nCentralityBins][nTrackPtBins];
  TH1D *testRatioProjection[nCentralityBins][nTrackPtBins];
  TH1D *manualProjection1[nCentralityBins][nTrackPtBins];
  TH1D *manualProjection2[nCentralityBins][nTrackPtBins];
  TH1D *manualTestRatio1[nCentralityBins][nTrackPtBins];
  TH1D *manualTestRatio2[nCentralityBins][nTrackPtBins];
  
  JDrawer *drawer = new JDrawer();
  
  int lowBin;
  int highBin;
  int nProjectedBins;
  
  // For a quick check, read the same histogram from the two files and check 1/2 = 1
  for(int iCentrality = 1; iCentrality < 2; iCentrality++){
    for(int iTrackPt = 5; iTrackPt < 6; iTrackPt++){
      firstTestHistogram[iCentrality][iTrackPt] = firstHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(0,comparedDistribution,DijetHistogramManager::kMaxAsymmetryBins,iCentrality,iTrackPt);
      secondTestHistogram[iCentrality][iTrackPt] = secondHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(0,comparedDistribution,DijetHistogramManager::kMaxAsymmetryBins,iCentrality,iTrackPt);
      
      // Divide the histograms, the ratio should be one if the histograms are the same
      testRatioHistogram[iCentrality][iTrackPt] = (TH2D*) firstTestHistogram[iCentrality][iTrackPt]->Clone(Form("testRatio%d%d",iCentrality,iTrackPt));
      testRatioHistogram[iCentrality][iTrackPt]->Divide(secondTestHistogram[iCentrality][iTrackPt]);
      
      drawer->DrawHistogram(testRatioHistogram[iCentrality][iTrackPt],"#phi","#eta",Form("C:%d T:%d",iCentrality,iTrackPt),"colz");
      
      // Find the bins to do the projection over a background region
      lowBin = firstTestHistogram[iCentrality][iTrackPt]->GetXaxis()->FindBin(1.5);
      highBin = firstTestHistogram[iCentrality][iTrackPt]->GetXaxis()->FindBin(3*TMath::Pi()/2-0.001);
      nProjectedBins = highBin - lowBin + 1;
      
      value250 = 0;
      for(int i = lowBin; i <= highBin; i++){
        for(int j = 250; j <= 250; j++){
          value250 += firstTestHistogram[iCentrality][iTrackPt]->GetBinContent(i,j);
        }
      }
      value250 = value250*firstTestHistogram[iCentrality][iTrackPt]->GetXaxis()->GetBinWidth(1)/(nProjectedBins*1.0);
      cout << "Calculated value in bin 250: " << value250 << endl;
      
      firstTestProjection[iCentrality][iTrackPt] = firstHistograms->GetHistogramJetTrackDeltaEta(0,comparedDistribution,DijetHistogramManager::kMaxAsymmetryBins,iCentrality,iTrackPt,DijetHistogramManager::kBetweenPeaks);
      secondTestProjection[iCentrality][iTrackPt] = secondHistograms->GetHistogramJetTrackDeltaEta(0,comparedDistribution,DijetHistogramManager::kMaxAsymmetryBins,iCentrality,iTrackPt,DijetHistogramManager::kBetweenPeaks);
      
      drawer->DrawHistogram(firstTestProjection[iCentrality][iTrackPt],"eta","dn/dEta",Form("C:%d T:%d",iCentrality,iTrackPt));
      
      // Divide the histograms, the ratio should be one if the histograms are the same
      testRatioProjection[iCentrality][iTrackPt] = (TH1D*) firstTestProjection[iCentrality][iTrackPt]->Clone(Form("testRatioProj%d%d",iCentrality,iTrackPt));
      testRatioProjection[iCentrality][iTrackPt]->Divide(secondTestProjection[iCentrality][iTrackPt]);
      
      drawer->DrawHistogram(testRatioProjection[iCentrality][iTrackPt],"eta","ratio",Form("C:%d T:%d",iCentrality,iTrackPt));
      
      cout << "Value in bin 250 from the file: " << firstTestProjection[iCentrality][iTrackPt]->GetBinContent(250) << endl;
      
      // Project the deltaEta from the two dimensional histograms
      manualProjection1[iCentrality][iTrackPt] = firstTestHistogram[iCentrality][iTrackPt]->ProjectionY(Form("manualProjection1%d%d",iCentrality,iTrackPt),lowBin,highBin);
      manualProjection2[iCentrality][iTrackPt] = secondTestHistogram[iCentrality][iTrackPt]->ProjectionY(Form("manualProjection2%d%d",iCentrality,iTrackPt),lowBin,highBin);
      
      manualProjection1[iCentrality][iTrackPt]->Scale(firstTestHistogram[iCentrality][iTrackPt]->GetXaxis()->GetBinWidth(1)/(nProjectedBins*1.0));
      manualProjection2[iCentrality][iTrackPt]->Scale(firstTestHistogram[iCentrality][iTrackPt]->GetXaxis()->GetBinWidth(1)/(nProjectedBins*1.0));
      
      cout << "Value in bin 250 from manual calculation: " << manualProjection1[iCentrality][iTrackPt]->GetBinContent(250) << endl;
      
      // Check the ratio with projection, they should be the same
      manualTestRatio1[iCentrality][iTrackPt] = (TH1D*) manualProjection1[iCentrality][iTrackPt]->Clone(Form("manualTestRatio1%d%d",iCentrality,iTrackPt));
      manualTestRatio2[iCentrality][iTrackPt] = (TH1D*) manualProjection2[iCentrality][iTrackPt]->Clone(Form("manualTestRatio2%d%d",iCentrality,iTrackPt));
      
      //drawer->DrawHistogram(manualProjection1[iCentrality][iTrackPt],"#Delta#eta","#frac{dN}{d#Delta#eta}"," ");
      
      
      manualTestRatio1[iCentrality][iTrackPt]->Divide(firstTestProjection[iCentrality][iTrackPt]);
      manualTestRatio2[iCentrality][iTrackPt]->Divide(secondTestProjection[iCentrality][iTrackPt]);
      
      drawer->DrawHistogram(manualTestRatio1[iCentrality][iTrackPt],"eta","ratio",Form("TestRatioPreprocess C:%d T:%d",iCentrality,iTrackPt));
      drawer->DrawHistogram(manualTestRatio2[iCentrality][iTrackPt],"eta","ratio",Form("TestRatioNopreprocess C:%d T:%d",iCentrality,iTrackPt));
    }
  }
  
}

