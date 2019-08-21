#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "stackHistXiao.h"
#include "xCanvas.h"
#include "JDrawer.h"

void plotJetShapeXiao(DijetHistogramManager *ppHistograms, DijetHistogramManager *pbpbHistograms){

  const int nTrackPtBins = pbpbHistograms->GetNTrackPtBins();
  const int nCentralityBins = pbpbHistograms->GetNCentralityBins();
  
  // Read the jet shape histograms to an array
  TH1D *jetShapeArray[nTrackPtBins][nCentralityBins+1];
  
  for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
    jetShapeArray[iTrackPt][nCentralityBins] = ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, DijetHistogramManager::kPtWeightedTrackLeadingJet, DijetHistogramManager::kMaxAsymmetryBins, 0, iTrackPt);
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      jetShapeArray[iTrackPt][iCentrality] = pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, DijetHistogramManager::kPtWeightedTrackLeadingJet, DijetHistogramManager::kMaxAsymmetryBins, iCentrality, iTrackPt);
    }
  }
  TString cent_lab[4] = {"0-10%", "10-30%", "30-50%", "50-100%"};
  stackHist *st[nCentralityBins+1];
  //labels_PAS* lb = new labels_PAS();
  TBox *box = new TBox();
  auto tl = new TLatex();
  
  TH1D* ratio[nCentralityBins];
  //TH1D* ratio_auxi[4];
  TLine *line = new TLine();
  
  TH1D* jetShapeSum[nCentralityBins+1];
  char namer[100];
  for(int iCentrality = 0; iCentrality < nCentralityBins+1; iCentrality++){
    sprintf(namer,"jetShapeSum%d",iCentrality);
    jetShapeSum[iCentrality] = (TH1D*) jetShapeArray[iCentrality][0]->Clone(namer);
    for(int iTrackPt = 1; iTrackPt < nTrackPtBins; iTrackPt++){
      jetShapeSum[iCentrality]->Add(jetShapeArray[iCentrality][iTrackPt]);
    }
  }
  
  for(int i = 0; i < nCentralityBins+1; i++){
    //js_dr_err_all[i]->SetFillStyle(1001);
    //js_dr_err_all[i]->SetFillColorAlpha(kGray+3,.4);
    //js_dr_err_all[i]->SetMarkerStyle(20);
    //js_dr_err_all[i]->SetMarkerSize(1.6);
    //js_dr_err_all[i]->SetMarkerColor(kBlack);
    //js_dr_err_all[i]->SetLineColor(kBlack);
    
    st[i]= new stackHist(Form("st_%d",i));
    st[i]->setRange(0., 0.99, "x");
    st[i]->setRange(0.005, 30, "y");
    int firstBin = 1;
    int lastBin = jetShapeSum[i]->GetXaxis()->FindBin(0.99);
    double fac = jetShapeSum[i]->Integral(firstBin,lastBin,"width");
    cout<<fac<<endl;
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins ; iTrackPt++){
      jetShapeArray[iTrackPt][i]->Scale(1.0/fac);
      st[i]->addHist((TH1*) jetShapeArray[iTrackPt][i]);
    }
    //js_dr_err_all[i]->Scale(1.0/fac);
    jetShapeSum[i]->Scale(1.0/fac);
    jetShapeSum[i]->SetMarkerColor(0);
  }
  
  for(int i=0 ;i<nCentralityBins; ++i){
    ratio[i]=(TH1D*)jetShapeSum[i]->Clone(Form("js_ratio_%d",i));
    ratio[i]->Divide(jetShapeSum[nCentralityBins]);
    //ratio_auxi[i]=(TH1D*) js_dr_err_all[i]->Clone(Form("js_syst_ratio_%d",i));
    //ratio_auxi[i]->Divide(js_dr_err_all[4]);
    //ratio_auxi[i]->SetFillStyle(1001);
    //ratio_auxi[i]->SetFillColorAlpha(kGray+3, .4);
    //ratio_auxi[i]->SetMarkerColor(0);
    //    ratio_auxi[i]->SetLineColor(kBlack);
  }
  auto c = new auxi_canvas("c", "", 2500, 2000);
  c->SetMargin(0.06, 0.01, 0.08, 0.02);
  c->divide(3,4);
  
  //cout<<ratio_auxi[0]->GetName()<<endl;
  //cout<<"err = "<<ratio_auxi[0]->GetBinError(ratio_auxi[0]->FindBin(0.95))<<endl;;
  
  for(int i=0; i<nCentralityBins; ++i){
    c->CD(8-i);
    gPad->SetLogy();
    st[i]->drawStack();
    st[i]->hst->GetYaxis()->SetNdivisions(505);
    st[i]->hst->GetYaxis()->SetLabelSize(0.08);
    st[i]->hst->GetYaxis()->SetTitleOffset(.9);
    st[i]->hst->GetYaxis()->SetTitleSize(0.1);
    st[i]->hst->GetYaxis()->SetTitle("#rho(#Deltar)");
    st[i]->hst->Draw();
    //js_dr_err_all[i]->Draw("same e2");
    if(i==3 ){
      tl->SetTextFont(22);
      tl->SetTextSize(.085);
      tl->DrawLatexNDC(.35, 0.9, "PbPb");
      tl->DrawLatexNDC(.35, .82, cent_lab[i]);
    }
    else {
      tl->SetTextFont(22);
      tl->SetTextSize(.09);
      tl->DrawLatexNDC(.19, 0.9, "PbPb");
      tl->DrawLatexNDC(.19, .82, cent_lab[i]);
    }
    c->CD(12-i);
    //ratio[i]->GetYaxis()->SetNdivisions(505);
    ratio[i]->GetXaxis()->SetTitleOffset(1.1);
    ratio[i]->GetXaxis()->SetTitle("#Deltar");
    ratio[i]->SetAxisRange(0., 6.2, "Y");
    ratio[i]->SetAxisRange(0, .99, "X");
    if( i<3 )  {
      ratio[i]->GetXaxis()->SetTitleOffset(0.7);
      ratio[i]->GetXaxis()->SetTitleSize(0.11);
      ratio[i]->GetXaxis()->SetNdivisions(505);
      ratio[i]->GetXaxis()->SetLabelSize(0.08);
      ratio[i]->GetYaxis()->SetNdivisions(505);
    }
    if(i==3 ){
      ratio[i]->GetXaxis()->SetTitleOffset(0.92);
      ratio[i]->GetXaxis()->SetTitleSize(0.086);
      ratio[i]->GetXaxis()->SetNdivisions(505);
      ratio[i]->GetXaxis()->SetLabelOffset(0.016);
      ratio[i]->GetXaxis()->SetLabelSize(0.064);
      ratio[i]->GetYaxis()->SetNdivisions(505);
      ratio[i]->GetYaxis()->SetLabelOffset(0.02);
      ratio[i]->GetYaxis()->SetLabelSize(0.07);
      ratio[i]->GetYaxis()->SetTitleOffset(0.9);
      ratio[i]->GetYaxis()->SetTitleSize(0.08);
      ratio[i]->GetYaxis()->SetTitle("#rho(#Deltar)_{PbPb}/#rho(#Deltar)_{pp}");
      tl->SetTextSize(0.073);
      //tl->DrawLatexNDC(.25, .92, "PbPb - pp");
      //tl->DrawLatexNDC(.25, .84, cent_lab[i]);
    }
    else{
      //tl->SetTextSize(0.09);
      //tl->DrawLatexNDC(.05, .92, "PbPb - pp");
      //tl->DrawLatexNDC(.05, .84, cent_lab[i]);
    }
    ratio[i]->GetYaxis()->CenterTitle();
    ratio[i]->GetXaxis()->CenterTitle();
    ratio[i]->SetStats(0);
    ratio[i]->Draw("same");
    //ratio_auxi[i]->Draw("same e2");
    line->SetLineStyle(2);
    line->DrawLine(0, 1, 1, 1);
  }
  c->CD(1);
  gPad->SetLogy();
  st[4]->drawStack();
  st[4]->hst->GetYaxis()->SetNdivisions(505);
  st[4]->hst->GetYaxis()->SetLabelSize(0.08);
  st[4]->hst->GetYaxis()->SetTitleOffset(.9);
  st[4]->hst->GetYaxis()->SetTitleSize(0.1);
  st[4]->hst->GetYaxis()->SetTitle("#rho(#Deltar)");
  tl->SetTextSize(.085);
  tl->DrawLatexNDC(.35, .88, "pp reference");
  //js_dr_err_all[4]->Draw("same e2");
  
  TLegend* lt1 = new TLegend(0.01,0.1,1.,0.5);
  TLegend* lt2 = new TLegend(0.0,0.1,1.,0.5);
  TLegend* lt3 = new TLegend(0.0,0.1,1 ,0.5);
  TLegend* lt5 = new TLegend(0.01,0.5,1 ,.68);
  //TLegend* lt4 = new TLegend(0.25,0.85,.85,0.95);
  lt1->SetTextSize(0.07);
  lt1->SetLineColor(kWhite);
  lt1->SetFillColor(kWhite);
  lt2->SetTextSize(0.07);
  lt2->SetLineColor(kWhite);
  lt2->SetFillColor(kWhite);
  //lt3->SetTextSize(0.07);
  //lt3->SetLineColor(kWhite);
  //lt3->SetFillColor(kWhite);
  //lt4->SetTextSize(0.06);
  //lt4->SetLineColor(kWhite);
  //lt4->SetFillColor(kWhite);
  lt5->SetTextSize(0.07);
  lt5->SetLineColor(kWhite);
  lt5->SetFillColor(kWhite);
  
  
  lt1->AddEntry(st[4]->hist_trunk.at(0), "0.7 < p_{T}^{trk}< 1 GeV","f");
  lt1->AddEntry(st[4]->hist_trunk.at(1), "1 < p_{T}^{trk}< 2 GeV","f");
  lt1->AddEntry(st[4]->hist_trunk.at(2), "2 < p_{T}^{trk}< 3 GeV","f");
  
  lt2->AddEntry(st[4]->hist_trunk.at(3), "3 < p_{T}^{trk}< 4 GeV","f");
  lt2->AddEntry(st[4]->hist_trunk.at(4), "4 < p_{T}^{trk}< 8 GeV","f");
  lt2->AddEntry(st[4]->hist_trunk.at(5), "8 < p_{T}^{trk}< 12 GeV","f");
  
  
  lt3->AddEntry(st[4]->hist_trunk.at(6), "12 < p_{T}^{trk}< 300 GeV","f");
  //lt3->AddEntry(st[4]->hist_trunk.at(7), "16 < p_{T}^{trk}< 20 GeV","f");
  //lt3->AddEntry(st[4]->hist_trunk.at(8), "20 < p_{T}^{trk}< 300 GeV","f");
  
  //lt4->AddEntry(ratio_auxi[0], "#rho(#Deltar)_{PbPb}/#rho(#Deltar)_{pp}");
  //lt4->AddEntry(ratio_auxi[0], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
  //c->CD(9);
  //lt4->Draw();
  
  //c->CD(1);
  //TLegend* lt6 = new TLegend(0.3,0.7,0.98,0.85);
  //lt6->SetTextSize(0.07);
  //lt6->SetLineColor(kWhite);
  //lt6->SetFillColor(kWhite);
  //lt6->AddEntry(js_dr_err_all[0], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
  //lt6->Draw();
  
  c->CD(2);
  //  lt5->AddEntry(js_dr_err_all[0], "0.7 < p_{T}^{track}< 20 GeV","lpfe");
  //  lt5->Draw();
  line->SetLineStyle(1);
  line->DrawLineNDC(0, 0, 0, 1);
  lt1->Draw();
  c->CD(3);
  lt2->Draw();
  //c->CD(4);
  //lt3->Draw();
  
  c->cd(0);
  tl->SetTextFont(62);
  tl->SetTextSize(0.035);
  tl->DrawLatexNDC(0.5, 0.94, "CMS");
  
  tl->SetTextFont(42);
  tl->SetTextSize(0.035);
  tl->DrawLatexNDC(0.6, 0.94, "Jet shapes");
  
  tl->SetTextSize(0.03);
  tl->DrawLatexNDC(.4, .9, "pp 27.4 pb^{-1} (5.02 TeV)  PbPb 404 #mub^{-1} (5.02 TeV)");
  tl->SetTextSize(0.025);
  tl->DrawLatexNDC(.45, .86, "anti-k_{T} R=0.4 jets, p_{T}> 120 GeV, |#eta_{jet}|<1.6");
  //  lb->drawText("(p_{T}> 120 GeV, |#eta_{jet}|<1.6)", 0.2, 0.25, 4);
  box->SetFillColor(kWhite);
  c->cd(0);
  box->DrawBox(0.285,.047, 0.3, 0.072);
  tl->SetTextSize(.025);
  tl->DrawLatex(0.29, 0.055, "0");
  box->DrawBox(0.518,.047, 0.533, 0.072);
  box->DrawBox(0.75,.047, 0.765, 0.072);
  
  tl->DrawLatex(0.523, 0.055, "0");
  tl->DrawLatex(0.755, 0.055, "0");
  tl->DrawLatex(0.985, 0.055, "1");
  
  //c->SaveAs("js_dr_normal_new.eps");
  //c->SaveAs("js_dr_normal_new.pdf");
  //c->SaveAs("js_dr_normal_v3.eps");
  //c->SaveAs("js_dr_normal_v3.pdf");
}

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void finalResultPlotter(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Open data files for pp and PbPb data
  TFile *ppFile = TFile::Open("data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_processed_2019-08-13.root");
  // data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_processed_2019-08-13.root
  // data/dijet_pp_highForest_pfJets_noUncOrInc_allCorrections_wtaAxis_processed_2019-07-13.root
  TFile *pbpbFile = TFile::Open("data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_modifiedSeagull_noErrorJff_averagePeakMixing_processed_2019-08-13_fiveJobsMissing.root");
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
  bool saveFigures = false;
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
  int iJetTrack = DijetHistogramManager::kPtWeightedTrackLeadingJet;  // DijetHistogramManager::kPtWeightedTrackSubleadingJet
  int lowestTrackPtBin = 0;
  
  TString asymmetryString[nAsymmetryBins+1];
  for(int i = 0; i < nAsymmetryBins; i++){
    asymmetryString[i] = Form("   %.1f < x_{j} < %.1f",asymmetryBinBorders[i],asymmetryBinBorders[i+1]);
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
  //plotJetShapeXiao(ppHistograms,pbpbHistograms);
  
  // Temporary: Get the ratio between pp and PbPb jet shape
  TH1D *sumHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *sumHistogramPp[nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPp[nAsymmetryBins+1];
  TH1D *helperHistogram;
  double shapeIntegralPbPb[nCentralityBins][nAsymmetryBins+1];
  double shapeIntegralPp[nAsymmetryBins+1];
  
  // Read the pT summed jet shape histograms
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      sumHistogramPbPb[iCentrality][iAsymmetry] = (TH1D*) pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, lowestTrackPtBin)->Clone(Form("sumPbPb%d",iCentrality));
    }
    sumHistogramPp[iAsymmetry] = (TH1D*) ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, 0)->Clone("sumPp");
  }
  
  // Sum the pT:s
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iTrackPt = lowestTrackPtBin+1; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        helperHistogram = (TH1D*)pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,iCentrality,iTrackPt)->Clone();
        sumHistogramPbPb[iCentrality][iAsymmetry]->Add(helperHistogram);
      }
      helperHistogram = (TH1D*)ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,0,iTrackPt)->Clone();
      sumHistogramPp[iAsymmetry]->Add(helperHistogram);
    } // track pT
  }
  
  // Normalize the jet shape histograms
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      shapeIntegralPbPb[iCentrality][iAsymmetry] = sumHistogramPbPb[iCentrality][iAsymmetry]->Integral(1,sumHistogramPbPb[iCentrality][iAsymmetry]->FindBin(0.99),"width");
      sumHistogramPbPb[iCentrality][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iCentrality][iAsymmetry]);
    }
    
    shapeIntegralPp[iAsymmetry] = sumHistogramPp[iAsymmetry]->Integral(1,sumHistogramPp[iAsymmetry]->FindBin(0.99),"width");
    sumHistogramPp[iAsymmetry]->Scale(1.0/shapeIntegralPp[iAsymmetry]);
    
    // Read the uncertainties for the pT summed jet shapes and scale them for jet shapes
    uncertaintyHistogramPp[iAsymmetry] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry);
    uncertaintyHistogramPp[iAsymmetry]->Scale(1.0/shapeIntegralPp[iAsymmetry]);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      uncertaintyHistogramPbPb[iCentrality][iAsymmetry] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, nTrackPtBins, iAsymmetry);
      uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iCentrality][iAsymmetry]);
    }
  }
  
  TH1D* ratioHistogram[nCentralityBins][nAsymmetryBins+1];
  TH1D* asymmetryRatioHistogram[nCentralityBins+1];
  
  // Create uncertainty histograms for the ratio using standard jet shape binning
  TH1D* ratioUncertainty[nCentralityBins][nAsymmetryBins+1];
  TH1D* asymmetryRatioUncertainty[nCentralityBins+1];
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      ratioUncertainty[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("ratioUncertainty%d%d", iCentrality, iAsymmetry));
    }
  }
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    asymmetryRatioUncertainty[iCentrality] = (TH1D*) sumHistogramPbPb[0][0]->Clone(Form("asymmetryRatioUncertainty%d", iCentrality));
  }
  
  JDrawer *drawer = new JDrawer();
  drawer->SetCanvasSize(500,500);
  
  TLegend *legend;
  double ppValue, pbpbValue, ratioValue;
  double ppUncertainty, pbpbUncertainty, ratioUncertaintyValue;
  
  TLine *oneLine = new TLine(0,1,1,1);
  oneLine->SetLineStyle(2);
  
  if(drawJetShapePptoPbPbRatio){
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        ratioHistogram[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("ratio%d%d", iCentrality, iAsymmetry));
        ratioHistogram[iCentrality][iAsymmetry]->Divide(sumHistogramPp[iAsymmetry]);
        
        // Calculate the systemtic uncertainty for the ratio
        for(int iBin = 1; iBin < ratioUncertainty[iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
          ppValue = sumHistogramPp[iAsymmetry]->GetBinContent(iBin);
          pbpbValue = sumHistogramPbPb[iCentrality][iAsymmetry]->GetBinContent(iBin);
          ppUncertainty = uncertaintyHistogramPp[iAsymmetry]->GetBinContent(iBin);
          pbpbUncertainty = uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->GetBinContent(iBin);
          ratioValue = ratioHistogram[iCentrality][iAsymmetry]->GetBinContent(iBin);
          ratioUncertaintyValue = TMath::Sqrt(TMath::Power(pbpbUncertainty/ppValue,2)+TMath::Power(pbpbValue*ppUncertainty/TMath::Power(ppValue,2),2));
          ratioUncertainty[iCentrality][iAsymmetry]->SetBinContent(iBin,ratioValue);
          ratioUncertainty[iCentrality][iAsymmetry]->SetBinError(iBin,ratioUncertaintyValue);
        }
        
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetRangeUser(0,1);
        ratioHistogram[iCentrality][iAsymmetry]->GetYaxis()->SetRangeUser(0,4);
        ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->SetRangeUser(0,1);
        ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetRangeUser(0,4);
        ratioUncertainty[iCentrality][iAsymmetry]->SetFillColorAlpha(29,0.25);
        
        drawer->DrawHistogram(ratioUncertainty[iCentrality][iAsymmetry],"#Deltar","#rho(#Deltar)_{PbPb} / #rho(#Deltar)_{pp}", " ","E2");
        oneLine->Draw();
        ratioHistogram[iCentrality][iAsymmetry]->SetLineColor(kBlack);
        ratioHistogram[iCentrality][iAsymmetry]->Draw("same");
        
        legend = new TLegend(0.25,0.73,0.45,0.88);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->SetHeader(Form("C: %.0f-%.0f%s",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
        legend->AddEntry(ratioHistogram[iCentrality][iAsymmetry],Form("p_{T} > %.1f",trackPtBinBorders[lowestTrackPtBin]),"l");
        legend->Draw();
        
        if(saveFigures){
          if(iAsymmetry == nAsymmetryBins){
            gPad->GetCanvas()->SaveAs(Form("figures/jetShapeRatio_%s_C=%.0f-%.0f_pT%d.pdf", pbpbHistograms->GetJetTrackHistogramName(iJetTrack), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1],lowestTrackPtBin));
          } else {
            saveName = Form("figures/jetShapeRatio_%s_A=%.1f-%.1f_C=%.0f-%.0f", pbpbHistograms->GetJetTrackHistogramName(iJetTrack), asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1], centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]);
            saveName.ReplaceAll(".","v");
            gPad->GetCanvas()->SaveAs(Form("%s.%s",saveName.Data(),figureFormat));
          }
        } // Saving figures
        
      } // Centrality loop
    } // Asymmetry loop
  } // Ratio between pp and PbPb in different asymmetry bins
  
  if(drawJetShapeSymmetricAsymmetricRatio){
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      
      if(iCentrality == nCentralityBins){
        asymmetryRatioHistogram[iCentrality] = (TH1D*) sumHistogramPp[2]->Clone(Form("ratioAsymmetry%d", iCentrality));
        asymmetryRatioHistogram[iCentrality]->Divide(sumHistogramPp[0]);
      } else {
        asymmetryRatioHistogram[iCentrality] = (TH1D*) sumHistogramPbPb[iCentrality][2]->Clone(Form("ratioAsymmetry%d", iCentrality));
        asymmetryRatioHistogram[iCentrality]->Divide(sumHistogramPbPb[iCentrality][0]);
      }
      
      // Calculate the systemtic uncertainty for the ratio
      for(int iBin = 1; iBin < asymmetryRatioUncertainty[iCentrality]->GetNbinsX(); iBin++){
        if(iCentrality == nCentralityBins){
          ppValue = sumHistogramPp[0]->GetBinContent(iBin);
          pbpbValue = sumHistogramPp[2]->GetBinContent(iBin);
          ppUncertainty = uncertaintyHistogramPp[0]->GetBinContent(iBin);
          pbpbUncertainty = uncertaintyHistogramPp[2]->GetBinContent(iBin);
        } else {
          ppValue = sumHistogramPbPb[iCentrality][0]->GetBinContent(iBin);
          pbpbValue = sumHistogramPbPb[iCentrality][2]->GetBinContent(iBin);
          ppUncertainty = uncertaintyHistogramPbPb[iCentrality][0]->GetBinContent(iBin);
          pbpbUncertainty = uncertaintyHistogramPbPb[iCentrality][2]->GetBinContent(iBin);
        }
        ratioValue = asymmetryRatioHistogram[iCentrality]->GetBinContent(iBin);
        ratioUncertaintyValue = TMath::Sqrt(TMath::Power(pbpbUncertainty/ppValue,2)+TMath::Power(pbpbValue*ppUncertainty/TMath::Power(ppValue,2),2));
        asymmetryRatioUncertainty[iCentrality]->SetBinContent(iBin,ratioValue);
        asymmetryRatioUncertainty[iCentrality]->SetBinError(iBin,ratioUncertaintyValue);
      }
      
      asymmetryRatioHistogram[iCentrality]->GetXaxis()->SetRangeUser(0,1);
      asymmetryRatioHistogram[iCentrality]->GetYaxis()->SetRangeUser(0,4);
      asymmetryRatioUncertainty[iCentrality]->GetXaxis()->SetRangeUser(0,1);
      asymmetryRatioUncertainty[iCentrality]->GetYaxis()->SetRangeUser(0,4);
      asymmetryRatioUncertainty[iCentrality]->SetFillColorAlpha(29,0.25);
      
      drawer->DrawHistogram(asymmetryRatioUncertainty[iCentrality],"#Deltar","#rho(#Deltar)_{PbPb} / #rho(#Deltar)_{pp}", " ","E2");
      oneLine->Draw();
      asymmetryRatioHistogram[iCentrality]->SetLineColor(kBlack);
      asymmetryRatioHistogram[iCentrality]->Draw("same");
      
      legend = new TLegend(0.25,0.73,0.45,0.88);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      if(iCentrality == nCentralityBins){
        legend->SetHeader("pp  Track-leading jet");
      } else {
        legend->SetHeader(Form("C: %.0f-%.0f  Track-leading jet",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
      }
      legend->AddEntry(asymmetryRatioHistogram[iCentrality],"0.8 < x_{j} < 1.0 / 0.0 < x_{j} < 0.6","l");
      legend->Draw();
      
      if(saveFigures){
        if(iCentrality == nCentralityBins){
          gPad->GetCanvas()->SaveAs(Form("figures/jetShapeAsymmetryRatio_%s_pp.pdf", pbpbHistograms->GetJetTrackHistogramName(iJetTrack)));
        } else {
          saveName = Form("figures/jetShapeAsymmetryRatio_%s_C=%.0f-%.0f", pbpbHistograms->GetJetTrackHistogramName(iJetTrack), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]);
          saveName.ReplaceAll(".","v");
          gPad->GetCanvas()->SaveAs(Form("%s.%s",saveName.Data(),figureFormat));
        }
      } // Saving figures
      
    } // Centrality loop
  } // Draw asymmetry ratio
  
}
