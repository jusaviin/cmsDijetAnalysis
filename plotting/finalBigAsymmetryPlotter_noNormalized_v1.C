#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "stackHistXiao.h"
#include "xCanvas.h"
#include "JDrawer.h"

void plotJetShapeBigAsymmetry(DijetHistogramManager *ppHistograms, DijetHistogramManager *pbpbHistograms, JffCorrector *ppUncertaintyProvider, JffCorrector *pbpbUncertaintyProvider, int iJetTrack){
  
  const int nTrackPtBins = pbpbHistograms->GetNTrackPtBins();
  const int nCentralityBins = pbpbHistograms->GetNCentralityBins();
  const int nAsymmetryBins = pbpbHistograms->GetNAsymmetryBins();
  
  TString dbtype[] = {"leadingJet", "subleadingJet"};
  const char* jetShapeTitle[] = {"Leading jet shape","Subleading jet shape","Jet shape"};
  const char* jetShapeSaveName[] = {"trackLeadingJet","trackSubleadingJet","trackInclusiveJet"};
  const char* xjString[] = {"0.0 < x_{j} < 0.6","0.6 < x_{j} < 0.8","0.8 < x_{j} < 1.0","x_{j} inclusive"};
  const char* asymmetrySaveName[] = {"_A=0v0-0v6","_A=0v6-0v8","_A=0v8-1v0",""};
  const char* saveString = {"JetShape"};
  
  bool drawUncertainties = true;
  bool monteCarloLabels = false;
  bool normalizeJetShape = false;           // True: draw rho. False: draw P.
  
  // Change the titles if the jet shape is not normalized to one
  if(!normalizeJetShape){
    jetShapeTitle[0] = "Leading jet radial momentum distribution";
    jetShapeTitle[1] = "Subleading jet radial momentum distribution";
    jetShapeTitle[2] = "Jet radial momentum distribution";
    saveString = "JetRadialMomentum";
  }
  
  // Temporary: Get the ratio between pp and PbPb jet shape
  TH1D *sumHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *sumHistogramPp[nAsymmetryBins+1];
  TH1D *jetShapeArray[nCentralityBins+1][nTrackPtBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPp[nAsymmetryBins+1];
  TH1D *sumUncertainty[nCentralityBins+1][nAsymmetryBins+1];
  TH1D *helperHistogram;
  double shapeIntegralPbPb[nCentralityBins][nAsymmetryBins+1];
  double shapeIntegralPp[nAsymmetryBins+1];
  
  // Read the pT summed jet shape histograms
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      sumHistogramPbPb[iCentrality][iAsymmetry] = (TH1D*) pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, 0)->Clone(Form("sumPbPb%d",iCentrality));
      jetShapeArray[iCentrality][0][iAsymmetry] = pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, 0);
    }
    sumHistogramPp[iAsymmetry] = (TH1D*) ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, 0)->Clone("sumPp");
    jetShapeArray[nCentralityBins][0][iAsymmetry] = ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, 0);
  }
  
  // Sum the pT:s
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iTrackPt = 1; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        helperHistogram = (TH1D*)pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,iCentrality,iTrackPt)->Clone();
        sumHistogramPbPb[iCentrality][iAsymmetry]->Add(helperHistogram);
        jetShapeArray[iCentrality][iTrackPt][iAsymmetry] = pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,iCentrality,iTrackPt);
      }
      helperHistogram = (TH1D*)ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,0,iTrackPt)->Clone();
      sumHistogramPp[iAsymmetry]->Add(helperHistogram);
      jetShapeArray[nCentralityBins][iTrackPt][iAsymmetry] = ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, iTrackPt);
    } // track pT
  }
  
  // Normalize the jet shape histograms
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      shapeIntegralPbPb[iCentrality][iAsymmetry] = sumHistogramPbPb[iCentrality][iAsymmetry]->Integral(1,sumHistogramPbPb[iCentrality][iAsymmetry]->FindBin(0.99),"width");
      if(normalizeJetShape) sumHistogramPbPb[iCentrality][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iCentrality][iAsymmetry]);
      
      // To draw the systematic uncertainties, clone the sum histograms to uncertainty histograms
      sumUncertainty[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("sumUncertaintyPbPb%d%d", iCentrality, iAsymmetry));
      
      if(normalizeJetShape){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          jetShapeArray[iCentrality][iTrackPt][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iCentrality][iAsymmetry]);
        }
      }
      
    }
    
    shapeIntegralPp[iAsymmetry] = sumHistogramPp[iAsymmetry]->Integral(1,sumHistogramPp[iAsymmetry]->FindBin(0.99),"width");
    if(normalizeJetShape) sumHistogramPp[iAsymmetry]->Scale(1.0/shapeIntegralPp[iAsymmetry]);
    
    // To draw the systematic uncertainties, clone the sum histograms to uncertainty histograms
    sumUncertainty[nCentralityBins][iAsymmetry] = (TH1D*) sumHistogramPp[iAsymmetry]->Clone(Form("sumUncertaintyPp%d", iAsymmetry));
    
    if(normalizeJetShape){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        jetShapeArray[nCentralityBins][iTrackPt][iAsymmetry]->Scale(1.0/shapeIntegralPp[iAsymmetry]);
      }
    }
    
    // Read the uncertainties for the pT summed jet shapes and scale them for jet shapes
    uncertaintyHistogramPp[iAsymmetry] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry);
    if(normalizeJetShape) uncertaintyHistogramPp[iAsymmetry]->Scale(1.0/shapeIntegralPp[iAsymmetry]);
    
    // Set the bin errors in the sumUncertainty histograms to match the uncertainties read from the file
    for(int iBin = 1; iBin <= sumUncertainty[nAsymmetryBins][iAsymmetry]->GetNbinsX(); iBin++){
      sumUncertainty[nCentralityBins][iAsymmetry]->SetBinError(iBin, uncertaintyHistogramPp[iAsymmetry]->GetBinContent(iBin));
    }
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      uncertaintyHistogramPbPb[iCentrality][iAsymmetry] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, nTrackPtBins, iAsymmetry);
      if(normalizeJetShape) uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iCentrality][iAsymmetry]);
      
      // Set the bin errors in the sumUncertainty histograms to match the uncertainties read from the file
      for(int iBin = 1; iBin <= sumUncertainty[iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
        
        sumUncertainty[iCentrality][iAsymmetry]->SetBinError(iBin, uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->GetBinContent(iBin));
      }
      
    }
  }
  
  TH1D* ratioHistogram[nCentralityBins][nAsymmetryBins+1];
  
  // Create uncertainty histograms for the ratio using standard jet shape binning
  TH1D* ratioUncertainty[nCentralityBins][nAsymmetryBins+1];
  double ppValue, pbpbValue, ratioValue;
  double ppUncertainty, pbpbUncertainty, ratioUncertaintyValue;
  
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      //ratioUncertainty[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("ratioUncertainty%d%d", iCentrality, iAsymmetry));
      ratioHistogram[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("ratio%d%d", iCentrality, iAsymmetry));
      ratioHistogram[iCentrality][iAsymmetry]->Divide(sumHistogramPp[iAsymmetry]);
      
      ratioUncertainty[iCentrality][iAsymmetry] = (TH1D*) sumUncertainty[iCentrality][iAsymmetry]->Clone(Form("uncertaintyOfRatio%d%d", iCentrality, iAsymmetry));
      ratioUncertainty[iCentrality][iAsymmetry]->Divide(sumUncertainty[nCentralityBins][iAsymmetry]);
      
    }
  }
  
  TString cent_lab[4] = {"0-10%", "10-30%", "30-50%", "50-90%"};
  stackHist *jetShapeStack[nCentralityBins+1][nAsymmetryBins+1];
  
  // Create the stacked jet shape histograms and set a good drawing style for the uncertainty histgorams
  for(int iCentrality = 0; iCentrality < nCentralityBins+1; iCentrality++){
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      
      sumUncertainty[iCentrality][iAsymmetry]->SetFillStyle(1001);
      sumUncertainty[iCentrality][iAsymmetry]->SetFillColorAlpha(kGray+3,0.4);
      sumUncertainty[iCentrality][iAsymmetry]->SetMarkerStyle(20);
      sumUncertainty[iCentrality][iAsymmetry]->SetMarkerSize(1.6);
      sumUncertainty[iCentrality][iAsymmetry]->SetMarkerColor(kBlack);
      sumUncertainty[iCentrality][iAsymmetry]->SetLineColor(kBlack);
      
      jetShapeStack[iCentrality][iAsymmetry] = new stackHist(Form("st_%d",iCentrality));
      jetShapeStack[iCentrality][iAsymmetry]->setRange(0.0, 0.99, "x");
      if(normalizeJetShape){
        jetShapeStack[iCentrality][iAsymmetry]->setRange(0.005, 30, "y");
      } else {
        jetShapeStack[iCentrality][iAsymmetry]->setRange(0.5, 3000, "y");
      }
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins ; iTrackPt++){
        jetShapeStack[iCentrality][iAsymmetry]->addHist((TH1*) jetShapeArray[iCentrality][iTrackPt][iAsymmetry]);
      }
      //js_dr_err_all[iCentrality]->Scale(1.0/fac);
      //jetShapeSum[iCentrality]->SetMarkerColor(0);
    } // Asymmetry loop
  } // Centrality loop
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
      ratioUncertainty[iCentrality][iAsymmetry]->SetFillStyle(1001);
      ratioUncertainty[iCentrality][iAsymmetry]->SetFillColorAlpha(kGray+3, 0.4);
      if(monteCarloLabels){
        ratioUncertainty[iCentrality][iAsymmetry]->SetMarkerColor(kBlack);
      } else {
        ratioUncertainty[iCentrality][iAsymmetry]->SetMarkerColor(0);
      }
      ratioUncertainty[iCentrality][iAsymmetry]->SetMarkerStyle(21);
      ratioUncertainty[iCentrality][iAsymmetry]->SetLineColor(kBlack);
    }
    if(monteCarloLabels){
    }
  }
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins+1; iAsymmetry++){
      cout<<iCentrality<<" : "<<iAsymmetry<<endl;
      ratioUncertainty[iCentrality][iAsymmetry]->SetName(dbtype[iJetTrack/3]+Form("_ratio_systUncert_%d_%d", iCentrality, iAsymmetry));
      ratioHistogram[iCentrality][iAsymmetry]->SetName(dbtype[iJetTrack/3]+Form("_ratio_hist_%d_%d", iCentrality, iAsymmetry));
      ratioUncertainty[iCentrality][iAsymmetry]->Write();
      ratioHistogram[iCentrality][iAsymmetry]->Write();
    }
  }
  for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins+1; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins+1; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        jetShapeArray[iCentrality][iTrackPt][iAsymmetry]->SetName(dbtype[iJetTrack/3]+Form("_js_hist_%d_%d_%d", iCentrality, iTrackPt, iAsymmetry));
        jetShapeArray[iCentrality][iTrackPt][iAsymmetry]->Write();
      }
      sumUncertainty[iCentrality][iAsymmetry]->SetName(dbtype[iJetTrack/3]+Form("_sum_systUncert_%d_%d", iCentrality, iAsymmetry));
      sumUncertainty[iCentrality][iAsymmetry]->Write();
      if(iCentrality == nCentralityBins) {
        sumHistogramPp[iAsymmetry]->SetName(dbtype[iJetTrack/3]+Form("_sum_hist_%d_%d", iCentrality, iAsymmetry));
        cout<<sumHistogramPp[iAsymmetry]->GetName()<<endl;
        sumHistogramPp[iAsymmetry]->Write();
      }else {
        sumHistogramPbPb[iCentrality][iAsymmetry]->SetName(dbtype[iJetTrack/3]+Form("_sum_hist_%d_%d", iCentrality, iAsymmetry));
        sumHistogramPbPb[iCentrality][iAsymmetry]->Write();
      }
    }
  }
  
  
  
}

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void finalBigAsymmetryPlotter_noNormalized_v1(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString dir_str = "../data";
  // Open data files for pp and PbPb data
  TFile *ppFile = TFile::Open("data/ppData2017_highForest_pfJets_20EveMixed_xjBins_wtaAxis_allCorrections_processed_2020-02-04.root");
  TFile *pbpbFile = TFile::Open("data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_allCorrections_finalTuning_onlyJetShapes_processed_2020-02-17.root");
  
  TFile *ppUncertaintyFile = TFile::Open("uncertainties/systematicUncertaintyForPp_20eveMix_xjBins_fixJES_2020-02-03.root");
  TFile *pbpbUncertaintyFile = TFile::Open("uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_addNewSources_smoothedPairBackground_2020-05-18.root");
  DijetHistogramManager *ppHistograms = new DijetHistogramManager(ppFile);
  DijetHistogramManager *pbpbHistograms = new DijetHistogramManager(pbpbFile);
  JffCorrector *ppUncertaintyProvider = new JffCorrector();
  JffCorrector *pbpbUncertaintyProvider = new JffCorrector();
  
  // Choose which figure sets to draw
  bool drawJetShapePptoPbPbRatio = true;
  bool drawJetShapeSymmetricAsymmetricRatio = false;
  
  // Choose if you want to write the figures to pdf file
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
  int jetTrackIndex[2] = {DijetHistogramManager::kPtWeightedTrackLeadingJet, DijetHistogramManager::kPtWeightedTrackSubleadingJet};
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
  ppHistograms->SetLoadTrackLeadingJetCorrelationsPtWeighted(true);
  ppHistograms->SetAsymmetryBinRange(0,nAsymmetryBins);
  ppHistograms->LoadProcessedHistograms();
  
  pbpbHistograms->SetLoadTrackLeadingJetCorrelationsPtWeighted(true);
  pbpbHistograms->SetAsymmetryBinRange(0,nAsymmetryBins);
  pbpbHistograms->LoadProcessedHistograms();
  
  // Load the systematic uncertainties from the files
  ppUncertaintyProvider->ReadSystematicFile(ppUncertaintyFile);
  pbpbUncertaintyProvider->ReadSystematicFile(pbpbUncertaintyFile);
  
  // Plot the figures using Xiao's plotting macro
  auto wf = TFile::Open("hin-19-013_jetShapes_errorUpdate.root", "recreate");
  plotJetShapeBigAsymmetry(ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kPtWeightedTrackLeadingJet);
  plotJetShapeBigAsymmetry(ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kPtWeightedTrackSubleadingJet);
  wf->Close();
}
