#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"
#include "DijetMethods.h"
#include "JffCorrector.h"

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void compareDeltaEta(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  const int nFilesToCompare = 1;
  
  // Open data files for pp and PbPb data
  TFile *ppFile = TFile::Open("data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_processed_2019-08-13.root");
  // data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_processed_2019-08-13.root
  // data/dijet_pp_highForest_pfJets_noUncOrInc_allCorrections_wtaAxis_processed_2019-07-13.root
  TFile *pbpbFile[nFilesToCompare];
  
  pbpbFile[0] = TFile::Open("data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsWithCentShift_trackDeltaRonlyLowPt_processed_2019-10-16.root");
  //pbpbFile[1] = TFile::Open("data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_allCorrections_noStatErrorFromCorrections_lowPtResidualTrack_processed_2019-10-07_fiveJobsMissing.root");
  //pbpbFile[2] = TFile::Open("data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_modifiedSeagull_noErrorJff_averagePeakMixing_processed_2019-08-13_fiveJobsMissing.root");
  //pbpbFile[3] = TFile::Open("data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_allCorrections_newTryOnSeagull_JECv4_processed_2019-08-13_fiveJobsMissing.root");
  // dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_JECv6_finalTrack_eschemeAxis_allCorrections_allPtTrackDeltaR_processed_2019-10-16.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trigger_JECv6_finalTrack_wtaAxis_allCorrectionsWithCentralityShift_modifiedDeltaR_processed_2019-10-18.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trigger_JECv6_finalTrack_wtaAxis_allCorrectionsWithCentralityShift_processed_2019-10-18.root
  
  // data/dijetPbPb2018_akPu4CaloJets_jet100trig_eschemeAxis_improvisedMixing_centFix_oldJetAndTrackCorr_allCorrectionsWithJEC2015_processed_2019-09-22.root
  // data/dijetPbPb2018_akPu4CaloJets_jet100trig_eschemeAxis_improvisedMixing_centFix_oldJetAndTrackCorr_adhocScaling_allCorrectionsWithJEC2015_processed_2019-09-22.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PFJets_jet80Trigger_5eveMixingFromWTA_allCorrections_eschemeAxis_JECv5b_processed_2019-09-10.root
  // data/dijetPbPb_akPu4CaloJets_eschemeAxis_allCorrectionsWith2018WtaMC_noUncorr_improvisedMixing_processed_2019-09-14.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_25eveMix_allCorrections_calo80Trigger_wtaAxis_JECv5b_processed_2019-09-10.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_calo80Trigger_weightedJffCorrection_wtaAxis_JECv5b_processed_2019-09-05.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_calo80Trigger_xjBins_wtaAxis_allCorrections_JECv5b_processed_2019-09-05.root
  // data/dijetPbPb2018_highForest_akPu4CaloJets_jet80Trigger_5eveMix_xjBins_wtaAxis_allCorrections_JECv5b_processed_2019-08-30.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_calo80Trigger_wtaAxis_allCorrections_JECv5b_processed_2019-09-05_part5missing.root
  //  data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_modifiedSeagull_noErrorJff_averagePeakMixing_processed_2019-08-13_fiveJobsMissing.root"
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_allCorrections_newTryOnSeagull_JECv4_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_testNoJffCorr_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_allCorrections_modifiedSeagull_wtaAxis_JECv4_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root
  
  // Read systemtatic uncertainty file
  TFile *uncertaintyFile = TFile::Open("uncertainties/systematicUncertaintyForPbPb_15percentSpill5Jff_2019-10-14.root");
  // "uncertainties/systematicUncertaintyForPbPb_15percentSpill5Jff_2019-10-14.root"
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_oldJES_15percentSpill10Jff_2019-10-17.root
  JffCorrector *uncertaintyProvider = new JffCorrector();
  uncertaintyProvider->ReadSystematicFile(uncertaintyFile);
  
  // Create histogram managers for pp and PbPb
  DijetHistogramManager *ppHistograms = new DijetHistogramManager(ppFile);
  DijetHistogramManager *pbpbHistograms[nFilesToCompare];
  for(int iFile = 0; iFile < nFilesToCompare; iFile++){
    pbpbHistograms[iFile] = new DijetHistogramManager(pbpbFile[iFile]);
  }
  // Choose which figure sets to draw
  bool drawJetShapePptoPbPbRatio = false;
  bool drawJetShapeSymmetricAsymmetricRatio = true;
  
  // Choose if you want to write the figures to pdf file
  const int nHistogramTypesToCompare = 1;
  bool saveFigures = true;
  
  // Get the number of asymmetry bins
  const int nAsymmetryBins = pbpbHistograms[0]->GetNAsymmetryBins();
  const int nTrackPtBins = pbpbHistograms[0]->GetNTrackPtBins();
  const int nCentralityBins = pbpbHistograms[0]->GetNCentralityBins();
  
  double centralityBinBorders[] = {0,10,30,50,90};
  double asymmetryBinBorders[] = {0,0.6,0.8,1};
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};
 
  // Logarithmic scales for figures
  bool logPt = true;          // pT distributions
  bool logCorrelation = true; // track-jet deltaPhi-deltaEta distributions
  bool logJetShape = true;    // Jet shapes
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Load leading and inclusive correlation histograms from the files
  ppHistograms->SetLoadTrackLeadingJetCorrelations(true);
  ppHistograms->SetLoadTrackInclusiveJetCorrelations(true);
  ppHistograms->SetLoad2DHistograms(true);
  ppHistograms->SetAsymmetryBinRange(0,nAsymmetryBins);
  ppHistograms->LoadProcessedHistograms();
  
  for(int iFile = 0; iFile < nFilesToCompare; iFile++){
    pbpbHistograms[iFile]->SetLoadTrackLeadingJetCorrelations(true);
    pbpbHistograms[iFile]->SetLoadTrackInclusiveJetCorrelations(true);
    pbpbHistograms[iFile]->SetLoad2DHistograms(true);
    pbpbHistograms[iFile]->SetAsymmetryBinRange(0,nAsymmetryBins);
    pbpbHistograms[iFile]->LoadProcessedHistograms();
  }
  
  // Open a file to include results from inclusive jet analysis HIN-16-020
  TFile *inclusiveResultFile = TFile::Open("data/publishedResults/officialHist_py_deta_16_020.root");
  
  // Methods for doing stuff
  DijetMethods *projector = new DijetMethods();
  double projectionRegion = 1;  // Region in deltaPhi which is used to project the deltaEta peak
  const int nDeltaEtaBinsRebin = 21;
  const double deltaEtaBinBordersRebin[nDeltaEtaBinsRebin+1] = {-4,-3,-2,-1.5,-1,-0.8,-0.6,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.6,0.8,1,1.5,2,3,4};
  
  // Histograms needed to calculate the stacked deltaEta distributions
  int jetTrackTypes[] = {DijetHistogramManager::kTrackLeadingJet,DijetHistogramManager::kTrackInclusiveJet};
  TH1D *deltaEtaArray[nFilesToCompare][nHistogramTypesToCompare][nCentralityBins+1][nTrackPtBins];
  TH1D *deltaEtaUncertainty[nFilesToCompare][nHistogramTypesToCompare][nCentralityBins+1][nTrackPtBins];
  TH2D *helperHistogram;
  TH1D *addedHistogram;
  TH1D *rebinnedHistogram;
  
  // Histgrams for comparison with inclusive results
  TH1D *inclusiveDeltaEtaArray[nCentralityBins+1][nTrackPtBins];
  
  // Read two dimensional deltaEta-deltaPhi histograms and project the deltaEta yields out from them
  for(int iFilledHistograms = 0; iFilledHistograms < nHistogramTypesToCompare; iFilledHistograms++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iFile = 0; iFile < nFilesToCompare; iFile++){
          
          helperHistogram = pbpbHistograms[iFile]->GetHistogramJetTrackDeltaEtaDeltaPhi(jetTrackTypes[iFilledHistograms], DijetHistogramManager::kBackgroundSubtracted, nAsymmetryBins, iCentrality, iTrackPt);
          
          addedHistogram = projector->ProjectRegionDeltaEta(helperHistogram, -projectionRegion, projectionRegion, Form("DeltaEtaStack%d%d%d%d%d",iFilledHistograms,nAsymmetryBins,iCentrality,iTrackPt,iFile));
          
          // Since we want to plot yield, we do not want to normalize over the number of bins projected over
          // but simply look at the yield in certain region
          addedHistogram->Scale(projector->GetNBinsProjectedOver());
          
          // The different pT bins can have different size, so we need to normalize the yield with the pT bin width
          addedHistogram->Scale(1/(pbpbHistograms[iFile]->GetTrackPtBinBorder(iTrackPt+1) - pbpbHistograms[iFile]->GetTrackPtBinBorder(iTrackPt)));
          
          // Rebin the histogram to match the binning in the inclusive jet shape paper
          deltaEtaArray[iFile][iFilledHistograms][iCentrality][iTrackPt] = (TH1D*) projector->RebinAsymmetric(addedHistogram, nDeltaEtaBinsRebin, deltaEtaBinBordersRebin)->Clone(Form("deltaEtaArray%d%d%d%d", iFilledHistograms, iCentrality, iTrackPt, iFile));
          
          // Construct also uncertainty histogram
          deltaEtaUncertainty[iFile][iFilledHistograms][iCentrality][iTrackPt] = (TH1D*) deltaEtaArray[iFile][iFilledHistograms][iCentrality][iTrackPt]->Clone(Form("deltaEtaUncertainty%d%d%d%d", iFilledHistograms, iCentrality, iTrackPt, iFile));
          
          addedHistogram = uncertaintyProvider->GetDeltaEtaSystematicUncertainty(0, iCentrality, iTrackPt, nAsymmetryBins);
          
          for(int iBin = 1; iBin <= deltaEtaUncertainty[iFile][iFilledHistograms][iCentrality][iTrackPt]->GetNbinsX(); iBin++){
            
            deltaEtaUncertainty[iFile][iFilledHistograms][iCentrality][iTrackPt]->SetBinError(iBin, addedHistogram->GetBinContent(iBin));
            
          } // bin loop
          
        } // File loop
      } // Centrality loop
      
      helperHistogram = ppHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(jetTrackTypes[iFilledHistograms], DijetHistogramManager::kBackgroundSubtracted, nAsymmetryBins, 0, iTrackPt);
      
      addedHistogram = projector->ProjectRegionDeltaEta(helperHistogram, -projectionRegion, projectionRegion, Form("DeltaEtaStack%d%d%d",iFilledHistograms,nCentralityBins,iTrackPt));
      
      // Since we want to plot yield, we do not want to normalize over the number of bins projected over
      // but simply look at the yield in certain region
      addedHistogram->Scale(projector->GetNBinsProjectedOver());
      
      // The different pT bins can have different size, so we need to normalize the yield with the pT bin width
      addedHistogram->Scale(1/(ppHistograms->GetTrackPtBinBorder(iTrackPt+1) - ppHistograms->GetTrackPtBinBorder(iTrackPt)));
      
      // Rebin the histogram to match the binning in the inclusive jet shape paper
      deltaEtaArray[0][iFilledHistograms][nCentralityBins][iTrackPt] = (TH1D*) projector->RebinAsymmetric(addedHistogram, nDeltaEtaBinsRebin, deltaEtaBinBordersRebin)->Clone(Form("deltaEtaArray%d%d%d", iFilledHistograms, nCentralityBins, iTrackPt));
    } // Track pT loop
  } // Filled histogram loop

  
  // Read the inclusive histograms to make a comparison to the published results
  int inclusiveCentralityBins[] = {0,10,30,50,100};
  const char *inclusiveTrackPtBins[] = {"07","1","2","3","4","8","12"};
  for(int iTrackPt = 0; iTrackPt < 6; iTrackPt++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      inclusiveDeltaEtaArray[iCentrality][iTrackPt] = (TH1D*) inclusiveResultFile->Get(Form("Proj_dEta_PbPb_Cent%d_Cent%d_TrkPt%s_TrkPt%s_Rebin_rebin", inclusiveCentralityBins[iCentrality], inclusiveCentralityBins[iCentrality+1], inclusiveTrackPtBins[iTrackPt], inclusiveTrackPtBins[iTrackPt+1]));
      
    } // Centrality loop for published inclusive results
    
    inclusiveDeltaEtaArray[nCentralityBins][iTrackPt] = (TH1D*) inclusiveResultFile->Get(Form("Proj_dEta_pp_TrkPt%s_TrkPt%s_Rebin_rebin", inclusiveTrackPtBins[iTrackPt], inclusiveTrackPtBins[iTrackPt+1]));
    
  } // Track pT loop for published inclusive results
  
  
  JDrawer *drawer = new JDrawer();
  TLegend *legend;
  
  TString centralityString;
  TString compactCentralityString;
  TString trackString;
  TString compactTrackString;
  
  TString figureName;
  int veryNiceColors[] = {kBlue,kRed,kRed,kGreen+3};
  const char* labels[] = {"Inclusive","No spillover shift","Calo100","Calo100 TRG spill"};
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    if(iCentrality == nCentralityBins){
      centralityString = "pp";
      compactCentralityString = "_pp";
    } else {
      centralityString = Form("C = %.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
      compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    }
    
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins-1; iTrackPt++){
      
      trackString = Form("%.1f < p_{T} < %.1f GeV",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
      compactTrackString = Form("_T=%.0f-%.0f",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
      
      //inclusiveDeltaEtaArray[iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-2,2);
      //drawer->DrawHistogram(inclusiveDeltaEtaArray[iCentrality][iTrackPt],"#Delta#eta","#frac{dN}{d#Delta#eta}"," ");
      
      deltaEtaUncertainty[0][0][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-2,2);
      deltaEtaUncertainty[0][0][iCentrality][iTrackPt]->SetFillColorAlpha(veryNiceColors[0], 0.2);
      drawer->DrawHistogram(deltaEtaUncertainty[0][0][iCentrality][iTrackPt],"#Delta#eta","#frac{dN}{d#Delta#eta}"," ","e2");
      
      inclusiveDeltaEtaArray[iCentrality][iTrackPt]->Draw("same");
      
      for(int iFilledHistograms = 0; iFilledHistograms < nHistogramTypesToCompare; iFilledHistograms++){
        for(int iFile = 0; iFile < nFilesToCompare; iFile++){
          deltaEtaArray[iFile][iFilledHistograms][iCentrality][iTrackPt]->SetLineColor(veryNiceColors[iFile]);
          deltaEtaArray[iFile][iFilledHistograms][iCentrality][iTrackPt]->Draw("same");
          
          
        }
      }
      
      legend = new TLegend(0.6,0.6,0.9,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0,centralityString.Data(),"");
      legend->AddEntry((TObject*) 0,trackString.Data(),"");
      legend->AddEntry(inclusiveDeltaEtaArray[iCentrality][iTrackPt],"Published inclusive","l");
      for(int iFilledHistograms = 0; iFilledHistograms < nHistogramTypesToCompare; iFilledHistograms++){
        for(int iFile = 0; iFile < nFilesToCompare; iFile++){
          legend->AddEntry(deltaEtaArray[iFile][iFilledHistograms][iCentrality][iTrackPt],labels[iFile],"l");
        }
      }
      legend->Draw();
      
      // Save the figures into a file
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/deltaEtaInclusiveComparison%s%s.pdf", compactCentralityString.Data(), compactTrackString.Data()));
      }
      
    } // Track pT loop
  } // Centrality loop
  
}
