#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"
#include "DijetMethods.h"

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void compareDeltaEta(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Open data files for pp and PbPb data
  TFile *ppFile = TFile::Open("data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_processed_2019-08-13.root");
  // data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_processed_2019-08-13.root
  // data/dijet_pp_highForest_pfJets_noUncOrInc_allCorrections_wtaAxis_processed_2019-07-13.root
  TFile *pbpbFile = TFile::Open("data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_calo80Trigger_xjBins_wtaAxis_allCorrections_JECv5b_processed_2019-09-05.root");
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_calo80Trigger_xjBins_wtaAxis_allCorrections_JECv5b_processed_2019-09-05.root
  // data/dijetPbPb2018_highForest_akPu4CaloJets_jet80Trigger_5eveMix_xjBins_wtaAxis_allCorrections_JECv5b_processed_2019-08-30.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_calo80Trigger_wtaAxis_allCorrections_JECv5b_processed_2019-09-05_part5missing.root
  //  data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_modifiedSeagull_noErrorJff_averagePeakMixing_processed_2019-08-13_fiveJobsMissing.root"
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_allCorrections_newTryOnSeagull_JECv4_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_testNoJffCorr_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_allCorrections_modifiedSeagull_wtaAxis_JECv4_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root
  
  // Create histogram managers for pp and PbPb
  DijetHistogramManager *ppHistograms = new DijetHistogramManager(ppFile);
  DijetHistogramManager *pbpbHistograms = new DijetHistogramManager(pbpbFile);
  // Choose which figure sets to draw
  bool drawJetShapePptoPbPbRatio = false;
  bool drawJetShapeSymmetricAsymmetricRatio = true;
  
  // Choose if you want to write the figures to pdf file
  const int nHistogramTypesToCompare = 2;
  bool saveFigures = false;
  
  // Get the number of asymmetry bins
  const int nAsymmetryBins = pbpbHistograms->GetNAsymmetryBins();
  const int nTrackPtBins = pbpbHistograms->GetNTrackPtBins();
  const int nCentralityBins = pbpbHistograms->GetNCentralityBins();
  
  double centralityBinBorders[] = {0,10,30,50,90};
  double asymmetryBinBorders[] = {0,0.6,0.8,1};
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};
 
  // Logarithmic scales for figures
  bool logPt = true;          // pT distributions
  bool logCorrelation = true; // track-jet deltaPhi-deltaEta distributions
  bool logJetShape = true;    // Jet shapes
  
  int iJetTrack = DijetHistogramManager::kTrackLeadingJet;  // DijetHistogramManager::kTrackSubleadingJet
  //int iJetTrack = DijetHistogramManager::kTrackInclusiveJet;
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Load only jet shape histograms from these
  ppHistograms->SetLoadTrackLeadingJetCorrelations(true);
  ppHistograms->SetLoadTrackInclusiveJetCorrelations(true);
  ppHistograms->SetLoad2DHistograms(true);
  ppHistograms->SetAsymmetryBinRange(0,nAsymmetryBins);
  ppHistograms->LoadProcessedHistograms();
  
  pbpbHistograms->SetLoadTrackLeadingJetCorrelations(true);
  pbpbHistograms->SetLoadTrackInclusiveJetCorrelations(true);
  pbpbHistograms->SetLoad2DHistograms(true);
  pbpbHistograms->SetAsymmetryBinRange(0,nAsymmetryBins);
  pbpbHistograms->LoadProcessedHistograms();
  
  // Open a file to include results from inclusive jet analysis HIN-16-020
  TFile *inclusiveResultFile = TFile::Open("data/publishedResults/officialHist_py_deta_16_020.root");
  
  // Methods for doing stuff
  DijetMethods *projector = new DijetMethods();
  double projectionRegion = 1;  // Region in deltaPhi which is used to project the deltaEta peak
  const int nDeltaEtaBinsRebin = 21;
  const double deltaEtaBinBordersRebin[nDeltaEtaBinsRebin+1] = {-4,-3,-2,-1.5,-1,-0.8,-0.6,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.6,0.8,1,1.5,2,3,4};
  
  // Histograms needed to calculate the stacked deltaEta distributions
  int jetTrackTypes[] = {DijetHistogramManager::kTrackLeadingJet,DijetHistogramManager::kTrackInclusiveJet};
  TH1D *deltaEtaArray[nHistogramTypesToCompare][nCentralityBins+1][nTrackPtBins];
  TH2D *helperHistogram;
  TH1D *addedHistogram;
  TH1D *rebinnedHistogram;
  
  // Histgrams for comparison with inclusive results
  TH1D *inclusiveDeltaEtaArray[nCentralityBins+1][nTrackPtBins];
  
  // Read two dimensional deltaEta-deltaPhi histograms and project the deltaEta yields out from them
  for(int iFilledHistograms = 0; iFilledHistograms < nHistogramTypesToCompare; iFilledHistograms++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        
        helperHistogram = pbpbHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(jetTrackTypes[iFilledHistograms], DijetHistogramManager::kBackgroundSubtracted, nAsymmetryBins, iCentrality, iTrackPt);
        
        addedHistogram = projector->ProjectRegionDeltaEta(helperHistogram, -projectionRegion, projectionRegion, Form("DeltaEtaStack%d%d%d%d",iFilledHistograms,nAsymmetryBins,iCentrality,iTrackPt));
        
        // Since we want to plot yield, we do not want to normalize over the number of bins projected over
        // but simply look at the yield in certain region
        addedHistogram->Scale(projector->GetNBinsProjectedOver());
        
        // The different pT bins can have different size, so we need to normalize the yield with the pT bin width
        addedHistogram->Scale(1/(pbpbHistograms->GetTrackPtBinBorder(iTrackPt+1) - pbpbHistograms->GetTrackPtBinBorder(iTrackPt)));
        
        // Rebin the histogram to match the binning in the inclusive jet shape paper
        deltaEtaArray[iFilledHistograms][iCentrality][iTrackPt] = (TH1D*) projector->RebinAsymmetric(addedHistogram, nDeltaEtaBinsRebin, deltaEtaBinBordersRebin)->Clone(Form("deltaEtaArray%d%d%d", iFilledHistograms, iCentrality, iTrackPt));
      } // Centrality loop
      
      helperHistogram = ppHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(jetTrackTypes[iFilledHistograms], DijetHistogramManager::kBackgroundSubtracted, nAsymmetryBins, 0, iTrackPt);
      
      addedHistogram = projector->ProjectRegionDeltaEta(helperHistogram, -projectionRegion, projectionRegion, Form("DeltaEtaStack%d%d%d",iFilledHistograms,nCentralityBins,iTrackPt));
      
      // Since we want to plot yield, we do not want to normalize over the number of bins projected over
      // but simply look at the yield in certain region
      addedHistogram->Scale(projector->GetNBinsProjectedOver());
      
      // The different pT bins can have different size, so we need to normalize the yield with the pT bin width
      addedHistogram->Scale(1/(ppHistograms->GetTrackPtBinBorder(iTrackPt+1) - ppHistograms->GetTrackPtBinBorder(iTrackPt)));
      
      // Rebin the histogram to match the binning in the inclusive jet shape paper
      deltaEtaArray[iFilledHistograms][nCentralityBins][iTrackPt] = (TH1D*) projector->RebinAsymmetric(addedHistogram, nDeltaEtaBinsRebin, deltaEtaBinBordersRebin)->Clone(Form("deltaEtaArray%d%d%d", iFilledHistograms, nCentralityBins, iTrackPt));
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
  int veryNiceColors[] = {kBlue,kRed};
  const char* labels[] = {"Leading jets","Inclusive jets"};
  
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    
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
      
      inclusiveDeltaEtaArray[iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-2,2);
      drawer->DrawHistogram(inclusiveDeltaEtaArray[iCentrality][iTrackPt],"#Delta#eta","#frac{dN}{d#Delta#eta}"," ");
      for(int iFilledHistograms = 0; iFilledHistograms < nHistogramTypesToCompare; iFilledHistograms++){
        deltaEtaArray[iFilledHistograms][iCentrality][iTrackPt]->SetLineColor(veryNiceColors[iFilledHistograms]);
        deltaEtaArray[iFilledHistograms][iCentrality][iTrackPt]->Draw("same");
      }
      
      legend = new TLegend(0.6,0.6,0.9,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0,centralityString.Data(),"");
      legend->AddEntry((TObject*) 0,trackString.Data(),"");
      legend->AddEntry(inclusiveDeltaEtaArray[iCentrality][iTrackPt],"Published inclusive","l");
      for(int iFilledHistograms = 0; iFilledHistograms < nHistogramTypesToCompare; iFilledHistograms++){
        legend->AddEntry(deltaEtaArray[iFilledHistograms][iCentrality][iTrackPt],labels[iFilledHistograms],"l");
      }
      legend->Draw();
      
      // Save the figures into a file
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/deltaEtaInclusiveComparison%s%s.pdf", compactCentralityString.Data(), compactTrackString.Data()));
      }
      
    } // Track pT loop
  } // Centrality loop
  
}
