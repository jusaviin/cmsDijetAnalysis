#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "stackHist.h"
#include "xCanvas.h"
#include "JDrawer.h"
#include "DijetMethods.h"

void plotDeltaEtaXiao(DijetHistogramManager *ppHistograms, DijetHistogramManager *pbpbHistograms, JffCorrector *ppUncertaintyProvider, JffCorrector *pbpbUncertaintyProvider, int iJetTrack){

  // Determine the number of bins from the histogram manager
  const int nTrackPtBins = pbpbHistograms->GetNTrackPtBins();
  const int nCentralityBins = pbpbHistograms->GetNCentralityBins();
  const int nAsymmetryBins = pbpbHistograms->GetNAsymmetryBins();
  
  // Define labels and plotting layout
  const char* deltaEtaTitle[] = {"Particle yield vs. leading #Delta#eta","Particle yield vs. subleading #Delta#eta","Particle yield vs. #Delta#eta"};
  const char* deltaEtaSaveName[] = {"trackLeadingJet","trackSubleadingJet","trackInclusiveJet"};
  const char* xjString[] = {"0.0 < x_{j} < 0.6","0.6 < x_{j} < 0.8","0.8 < x_{j} < 1.0","x_{j} inclusive"};
  const char* asymmetrySaveName[] = {"_A=0v0-0v6","_A=0v6-0v8","_A=0v8-1v0",""};
  double deltaEtaZoom[] = {30,40,30};
  double subtractionZoom[] = {10,18,10};
  
  double centralityBinBorders[] = {0,10,30,50,90};
  double asymmetryBinBorders[] = {0,0.6,0.8,1};
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};
  
  // Open a file to include results from inclusive jet analysis HIN-16-020
  TFile *inclusiveResultFile = TFile::Open("data/publishedResults/officialHist_py_deta_16_020.root");
  
  // Methods for doing stuff
  DijetMethods *projector = new DijetMethods();
  double projectionRegion = 1;  // Region in deltaPhi which is used to project the deltaEta peak
  const int nDeltaEtaBinsRebin = 21;
  const double deltaEtaBinBordersRebin[nDeltaEtaBinsRebin+1] = {-4,-3,-2,-1.5,-1,-0.8,-0.6,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.6,0.8,1,1.5,2,3,4};
  
  // Histograms needed to calculate the stacked deltaEta distributions
  TH1D *sumHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *sumHistogramPp[nAsymmetryBins+1];
  TH1D *deltaEtaArray[nCentralityBins+1][nTrackPtBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPp[nAsymmetryBins+1];
  TH2D *helperHistogram;
  TH1D *addedHistogram;
  TH1D *rebinnedHistogram;
  double shapeIntegralPbPb[nCentralityBins][nAsymmetryBins+1];
  double shapeIntegralPp[nAsymmetryBins+1];
  
  // Histgrams for comparison with inclusive results
  TH1D *sumHistogramInclusive[nCentralityBins+1];
  TH1D *inclusiveDeltaEtaArray[nCentralityBins+1][nTrackPtBins];
  
  // Read two dimensional deltaEta-deltaPhi histograms and project the deltaEta yields out from them
  for(int iAsymmetry = nAsymmetryBins; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      helperHistogram = pbpbHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, 0);
      
      addedHistogram = projector->ProjectRegionDeltaEta(helperHistogram, -projectionRegion, projectionRegion, Form("DeltaEtaStack%d%d%d%d",iJetTrack,iAsymmetry,iCentrality,0));
      
      // Since we want to plot yield, we do not want to normalize over the number of bins projected over
      // but simply look at the yield in certain region
      addedHistogram->Scale(projector->GetNBinsProjectedOver());
      
      // The different pT bins can have different size, so we need to normalize the yield with the pT bin width
      addedHistogram->Scale(1/(pbpbHistograms->GetTrackPtBinBorder(1) - pbpbHistograms->GetTrackPtBinBorder(0)));
      
      // Rebin the histogram to match the binning in the inclusive jet shape paper
      sumHistogramPbPb[iCentrality][iAsymmetry] = (TH1D*) projector->RebinAsymmetric(addedHistogram, nDeltaEtaBinsRebin, deltaEtaBinBordersRebin)->Clone(Form("sumPbPb%d%d",iCentrality,iAsymmetry));

      deltaEtaArray[iCentrality][0][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("deltaEtaArray%d%d%d", iCentrality, iAsymmetry, 0));
    }
    
    helperHistogram = ppHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, 0, 0);
    
    addedHistogram = projector->ProjectRegionDeltaEta(helperHistogram, -projectionRegion, projectionRegion, Form("DeltaEtaStack%d%d%d%d",iJetTrack,iAsymmetry,nCentralityBins,0));
    
    // Since we want to plot yield, we do not want to normalize over the number of bins projected over
    // but simply look at the yield in certain region
    addedHistogram->Scale(projector->GetNBinsProjectedOver());
    
    // The different pT bins can have different size, so we need to normalize the yield with the pT bin width
    addedHistogram->Scale(1/(ppHistograms->GetTrackPtBinBorder(1) - ppHistograms->GetTrackPtBinBorder(0)));
    
    // Rebin the histogram to match the binning in the inclusive jet shape paper
    sumHistogramPp[iAsymmetry] = (TH1D*) projector->RebinAsymmetric(addedHistogram, nDeltaEtaBinsRebin, deltaEtaBinBordersRebin)->Clone(Form("sumPp%d",iAsymmetry));
    
    deltaEtaArray[nCentralityBins][0][iAsymmetry] = (TH1D*) sumHistogramPp[iAsymmetry]->Clone(Form("deltaEtaArray%d%d%d", nCentralityBins, iAsymmetry, 0));
    
  }
  
  // Sum the pT:s
  for(int iAsymmetry = nAsymmetryBins; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iTrackPt = 1; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        
        helperHistogram = pbpbHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, iTrackPt);
        
        addedHistogram = projector->ProjectRegionDeltaEta(helperHistogram, -projectionRegion, projectionRegion, Form("DeltaEtaStack%d%d%d%d",iJetTrack,iAsymmetry,iCentrality,iTrackPt));
        
        // Since we want to plot yield, we do not want to normalize over the number of bins projected over
        // but simply look at the yield in certain region
        addedHistogram->Scale(projector->GetNBinsProjectedOver());
        
        // The different pT bins can have different size, so we need to normalize the yield with the pT bin width
        addedHistogram->Scale(1/(pbpbHistograms->GetTrackPtBinBorder(iTrackPt+1) - pbpbHistograms->GetTrackPtBinBorder(iTrackPt)));
        
        // Rebin the histogram to match the binning in the inclusive jet shape paper
        rebinnedHistogram = projector->RebinAsymmetric(addedHistogram,nDeltaEtaBinsRebin,deltaEtaBinBordersRebin);
        
        sumHistogramPbPb[iCentrality][iAsymmetry]->Add(rebinnedHistogram);
        
        deltaEtaArray[iCentrality][iTrackPt][iAsymmetry] = (TH1D*) rebinnedHistogram->Clone(Form("deltaEtaArray%d%d%d", iCentrality, iAsymmetry, iTrackPt));
      }
      
      helperHistogram = ppHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, 0, iTrackPt);
      
      addedHistogram = projector->ProjectRegionDeltaEta(helperHistogram, -projectionRegion, projectionRegion, Form("DeltaEtaStack%d%d%d%d",iJetTrack,iAsymmetry,nCentralityBins,iTrackPt));
      
      // Since we want to plot yield, we do not want to normalize over the number of bins projected over
      // but simply look at the yield in certain region
      addedHistogram->Scale(projector->GetNBinsProjectedOver());
      
      // The different pT bins can have different size, so we need to normalize the yield with the pT bin width
      addedHistogram->Scale(1/(ppHistograms->GetTrackPtBinBorder(iTrackPt+1) - ppHistograms->GetTrackPtBinBorder(iTrackPt)));
      
      // Rebin the histogram to match the binning in the inclusive jet shape paper
      rebinnedHistogram = projector->RebinAsymmetric(addedHistogram,nDeltaEtaBinsRebin,deltaEtaBinBordersRebin);
      
      sumHistogramPp[iAsymmetry]->Add(rebinnedHistogram);
      
      deltaEtaArray[nCentralityBins][iTrackPt][iAsymmetry] = (TH1D*) rebinnedHistogram->Clone(Form("deltaEtaArray%d%d%d", nCentralityBins, iAsymmetry, iTrackPt));
    } // track pT
  }
  
  // Read the inclusive histograms to make a comparison to the published results
  int inclusiveCentralityBins[] = {0,10,30,50,100};
  int inclusiveTrackPtBins[] = {0,1,2,3,4,8,12};
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    sumHistogramInclusive[iCentrality] = (TH1D*) inclusiveResultFile->Get(Form("Proj_dEta_PbPb_Cent%d_Cent%d_TrkPt07_TrkPt1_Rebin_rebin", inclusiveCentralityBins[iCentrality], inclusiveCentralityBins[iCentrality+1]))->Clone(Form("inclusiveResults%d", iCentrality));
    
    inclusiveDeltaEtaArray[iCentrality][0] = (TH1D*) inclusiveResultFile->Get(Form("Proj_dEta_PbPb_Cent%d_Cent%d_TrkPt07_TrkPt1_Rebin_rebin", inclusiveCentralityBins[iCentrality], inclusiveCentralityBins[iCentrality+1]))->Clone(Form("inclusiveArray%d0", iCentrality));
    
    for(int iTrackPt = 1; iTrackPt < 6; iTrackPt++){
      
      inclusiveDeltaEtaArray[iCentrality][iTrackPt] = (TH1D*) inclusiveResultFile->Get(Form("Proj_dEta_PbPb_Cent%d_Cent%d_TrkPt%d_TrkPt%d_Rebin_rebin", inclusiveCentralityBins[iCentrality], inclusiveCentralityBins[iCentrality+1], inclusiveTrackPtBins[iTrackPt], inclusiveTrackPtBins[iTrackPt+1]));
      
      sumHistogramInclusive[iCentrality]->Add(inclusiveDeltaEtaArray[iCentrality][iTrackPt]);
      
    } // Track pT loop for published inclusive results
    
  } // Centrality loop for published inclusive results
  
  sumHistogramInclusive[nCentralityBins] = (TH1D*) inclusiveResultFile->Get("Proj_dEta_pp_TrkPt07_TrkPt1_Rebin_rebin")->Clone("inclusiveResultsPp");
  
  inclusiveDeltaEtaArray[nCentralityBins][0] = (TH1D*) inclusiveResultFile->Get("Proj_dEta_pp_TrkPt07_TrkPt1_Rebin_rebin")->Clone("inclusiveArrayPp");
  
  for(int iTrackPt = 1; iTrackPt < 6; iTrackPt++){
    
    inclusiveDeltaEtaArray[nCentralityBins][iTrackPt] = (TH1D*) inclusiveResultFile->Get(Form("Proj_dEta_pp_TrkPt%d_TrkPt%d_Rebin_rebin", inclusiveTrackPtBins[iTrackPt], inclusiveTrackPtBins[iTrackPt+1]));
    
    sumHistogramInclusive[nCentralityBins]->Add(inclusiveDeltaEtaArray[nCentralityBins][iTrackPt]);
    
  } // Track pT loop for published inclusive results
  
  
  JDrawer *drawer = new JDrawer();
  TLegend *legend;
  //drawer->SetDefaultAppearanceSplitCanvas();
  //TH1D *deltaEtaRatio[nCentralityBins][nTrackPtBins]
  
  TString centralityString;
  TString trackString;
  
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    
    if(iCentrality == nCentralityBins){
      centralityString = "pp";
    } else {
      centralityString = Form("C = %.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    }
    
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins-1; iTrackPt++){
      
      trackString = Form("%.1f < p_{T} < %.1f GeV",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
      
      //drawer->CreateSplitCanvas();
      //drawer->DrawHistogramToUpperPad(deltaEtaArray[iCentrality][iTrackPt][nAsymmetryBins],"#Delta#eta","##frac{dN}{d#Delta#eta}"," ");
      inclusiveDeltaEtaArray[iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-2,2);
      drawer->DrawHistogram(inclusiveDeltaEtaArray[iCentrality][iTrackPt],"#Delta#eta","#frac{dN}{d#Delta#eta}"," ");
      deltaEtaArray[iCentrality][iTrackPt][nAsymmetryBins]->SetLineColor(kRed);
      deltaEtaArray[iCentrality][iTrackPt][nAsymmetryBins]->Draw("same");
      
      legend = new TLegend(0.6,0.6,0.9,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0,centralityString.Data(),"");
      legend->AddEntry((TObject*) 0,trackString.Data(),"");
      legend->AddEntry(inclusiveDeltaEtaArray[iCentrality][iTrackPt],"Published inclusive","l");
      legend->AddEntry(deltaEtaArray[iCentrality][iTrackPt][nAsymmetryBins],"This analysis","l");
      legend->Draw();
    }
  }
}

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
  TFile *pbpbFile = TFile::Open("data/dijetPbPb2018_highForest_akPu4CaloJets_jet80Trigger_5eveMix_xjBins_wtaAxis_allCorrections_JECv5b_processed_2019-08-30.root");
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_calo80Trigger_xjBins_wtaAxis_allCorrections_JECv5b_processed_2019-09-05.root
  // data/dijetPbPb2018_highForest_akPu4CaloJets_jet80Trigger_5eveMix_xjBins_wtaAxis_allCorrections_JECv5b_processed_2019-08-30.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_calo80Trigger_wtaAxis_allCorrections_JECv5b_processed_2019-09-05_part5missing.root
  //  data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_modifiedSeagull_noErrorJff_averagePeakMixing_processed_2019-08-13_fiveJobsMissing.root"
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_allCorrections_newTryOnSeagull_JECv4_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_testNoJffCorr_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_allCorrections_modifiedSeagull_wtaAxis_JECv4_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root
  
  TFile *ppUncertaintyFile = TFile::Open("systematicTestPpUnmix.root");
  TFile *pbpbUncertaintyFile = TFile::Open("systematicTestPbPbUnmix.root");
  
  // Create histogram managers for pp and PbPb
  DijetHistogramManager *ppHistograms = new DijetHistogramManager(ppFile);
  DijetHistogramManager *pbpbHistograms = new DijetHistogramManager(pbpbFile);
  JffCorrector *ppUncertaintyProvider = new JffCorrector();
  JffCorrector *pbpbUncertaintyProvider = new JffCorrector();
  
  // Choose which figure sets to draw
  bool drawJetShapePptoPbPbRatio = false;
  bool drawJetShapeSymmetricAsymmetricRatio = true;
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = true;
  const char* figureFormat = "pdf";
  TString saveName;
  
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
  
  //int iAsymmetry = 0;
  //int iJetTrack = DijetHistogramManager::kTrackLeadingJet;  // DijetHistogramManager::kTrackSubleadingJet
  int iJetTrack = DijetHistogramManager::kTrackInclusiveJet;
  int lowestTrackPtBin = 0;
  
  TString asymmetryString[nAsymmetryBins+1];
  TString asymmetrySaveString[nAsymmetryBins+1];
  for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
    asymmetryString[iAsymmetry] = Form("   %.1f < x_{j} < %.1f",asymmetryBinBorders[iAsymmetry],asymmetryBinBorders[iAsymmetry+1]);
    asymmetrySaveString[iAsymmetry] = Form("_A=%.1f-%.1f",asymmetryBinBorders[iAsymmetry],asymmetryBinBorders[iAsymmetry+1]);
  }
  asymmetryString[nAsymmetryBins] = "";
  
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
  
  // Load the systematic uncertainties from the files
  ppUncertaintyProvider->ReadSystematicFile(ppUncertaintyFile);
  pbpbUncertaintyProvider->ReadSystematicFile(pbpbUncertaintyFile);
  
  // Plot the figures using Xiao's plotting macro
  plotDeltaEtaXiao(ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, iJetTrack);
  
}
