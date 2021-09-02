#include "JDrawer.h"
#include <tuple>

/*
 * Fit a Gauss function to a histogram and extract parameters from that
 *
 *  TH1* histogram = Histogram from which drawing range in searched
 *
 *  return: Gauss mean, Gauss sigma, Error for Gauss mean, Error for Gauss sigma
 */
std::tuple<double,double,double,double> fitGauss(TH1* histogram, double lowFitRange, double highFitRange, TString title = "", TString bin = "", TString saveName = ""){
  histogram->Fit("gaus","Q","",lowFitRange,highFitRange);
  TF1* gaussFit = histogram->GetFunction("gaus");
  
  double gaussMean = 0;
  double gaussSigma = 0;
  double gaussMeanError = 0;
  double gaussSigmaError = 0;
  
  if(gaussFit){
    gaussMean = gaussFit->GetParameter(1);
    gaussSigma = gaussFit->GetParameter(2);
    gaussMeanError = gaussFit->GetParError(1);
    gaussSigmaError = gaussFit->GetParError(2);
  }
  
  // If title is given, print the fit
  if(!title.EqualTo("")){
    JDrawer *temporaryDrawer = new JDrawer();
    temporaryDrawer->DrawHistogram(histogram,"Peak strip p_{T} [GeV]","Reflected strip p_{T}", " ");
    TLegend *legend = new TLegend(0.62,0.7,0.85,0.85);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->SetHeader(title);
    legend->AddEntry(histogram,bin,"l");
    legend->Draw();
    
    if(!saveName.EqualTo("")){
      gPad->GetCanvas()->SaveAs(saveName);
    }
    
  }
  
  return std::make_tuple(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError);
}


/*
 * Macro for plotting manual JEC debugging histograms
 */
void getAndPlotJecDebug(){

  // Open the data file
  TFile *inputFile = TFile::Open("data/PbPbMC2018_RecoGen_akPfCsJet_onlyJets_4pCentShift_subeNon0_manualJecHalfStrip_2021-08-26_mostStats.root");
  // PbPbMC2018_RecoGen_akPfCsJet_onlyJets_4pCentShift_subeNon0_manualJecHalfStrip_2021-08-26_mostStats.root
  // PbPbMC2018_RecoGen_akCaloJet_onlyRegular_4pCentShift_subeNon0_manualJecDijetScale_2021-08-25.root
  // PbPbMC2018_RecoGen_akPfCsJet_onlyJets_4pCentShift_subeNon0_sAndLJecDebug_lowPtDijet_2021-08-24.root
  // PbPbMC2018_RecoGen_akCaloJet_onlyJets_4pCentShift_subeNon0_sAndLJecDebug_lowPtDijet_2021-08-24.root
  // PbPbMC2018_RecoGen_akPfCsJet_onlyRegular_4pCentShift_subeNon0_etaReflectRandomSeed_2021-08-23_all.root
  // PbPbMC2018_RecoGen_akPfCsJet_onlyRegular_4pCentShift_subeNon0_etaReflectGoodScale_2021-08-20.root
  // PbPbMC2018_RecoGen_akPfCsJet_onlyRegular_4pCentShift_subeNon0_etaReflectNewScale_2021-08-19.root
  // PbPbMC2018_RecoGen_akPfCsJet_onlyRegular_4pCentShift_subeNon0_etaStripCheckNeutral_2021-08-19.root
  // PbPbMC2018_RecoGen_akPfCsJet_onlyRegular_4pCentShift_subeNon0_jecScalingDebug_2021-08-17.root
  // PbPbMC2018_RecoGen_akPfCsJet_onlyJets_4pCentShift_subeNon0_jecDebugHistograms_2021-08-16.root
  // PbPbMC2018_RecoGen_akCaloJet_onlyJets_4pCentShift_matchedPtDifference_etaStrip_2021-08-12.root
  // PbPbMC2018_RecoGen_akCaloJet_onlyJets_4pCentShift_matchedPtDifference_reflectCone_2021-08-12.root
  // PbPbMC2018_RecoGen_akCaloJet_onlyJets_4pCentShift_matchedPtDifference_reflectStrip_2021-08-12.root
  // PbPbMC2018_RecoGen_akPfCsJet_onlyJets_4pCentShift_matchedPtDifference_etaStrip_2021-08-12.root
  // PbPbMC2018_RecoGen_akPfCsJet_onlyJets_4pCentShift_matchedPtDifference_reflectCone_2021-08-12.root
  // PbPbMC2018_RecoGen_akPfCsJet_onlyJets_4pCentShift_matchedPtDifference_reflectStrip_2021-08-12.root
  
  // Figure saving flag
  const bool drawAllJets = false;
  const bool drawLeadingJets = false;
  const bool drawSubleadingJets = false;
  const bool drawJetsAbove100GeV = false;
  const bool drawFlippedJets = false;
  const bool saveFigures = false;
  const char* saveComment = "_scaledReflection";
  const bool drawStripIllustration = true;
  
  // Gauss fit flags
  const bool fitEtaStrip = false;
  const bool fitJetCone = false;
  const bool fitManualCorrection = false;
  
  const bool fitLeading = false;
  const bool fitSubleading = false;
  
  // Possibility to save histograms into an output file
  TString outputFileName = "";
  
  // Read the histograms from the file
  const int nCentralityBins = 4;
  const int nTypes = 3;
  TString typeString[3] = {"leading","subleading","inclusive"};
  TH2D *energyInEtaStrip[nTypes][nCentralityBins];
  TH2D *energyInCone[nTypes][nCentralityBins];
  TH2D *manualJetCorrection[nTypes][nCentralityBins];
  TH2D *manualJetCorrectionRatio[nTypes][nCentralityBins];
  TH2D *energyInEtaStripAbove100GeV[nCentralityBins];
  TH2D *energyInConeAbove100GeV[nCentralityBins];
  TH2D *manualJetCorrectionAbove100GeV[nCentralityBins];
  TH2D *manualJetCorrectionRatioAbove100GeV[nCentralityBins];
  TH2D *energyInEtaStripFlipped[nCentralityBins];
  TH2D *energyInConeFlipped[nCentralityBins];
  TH2D *manualJetCorrectionFlipped[nCentralityBins];
  TH2D *manualJetCorrectionRatioFlipped[nCentralityBins];
  
  // Histograms for projections
  TH1D *currentProjection;
  
  char newName[200];
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iType = 0; iType < nTypes; iType++){
      energyInEtaStrip[iType][iCentrality] = (TH2D*) inputFile->Get(Form("energyInEtaStrip_%s_C%d", typeString[iType].Data(), iCentrality));
      energyInCone[iType][iCentrality] = (TH2D*) inputFile->Get(Form("energyInCone_%s_C%d", typeString[iType].Data(), iCentrality));
      manualJetCorrection[iType][iCentrality] = (TH2D*) inputFile->Get(Form("manualJetCorrection_%s_C%d", typeString[iType].Data(), iCentrality));
      manualJetCorrectionRatio[iType][iCentrality] = (TH2D*) inputFile->Get(Form("manualJetCorrectionRatio_%s_C%d", typeString[iType].Data(), iCentrality));
    }
    energyInEtaStripAbove100GeV[iCentrality] = (TH2D*) inputFile->Get(Form("energyInEtaStripAbove100GeV_C%d",iCentrality));
    energyInConeAbove100GeV[iCentrality] = (TH2D*) inputFile->Get(Form("energyInConeAbove100GeV_C%d",iCentrality));
    manualJetCorrectionAbove100GeV[iCentrality] = (TH2D*) inputFile->Get(Form("manualJetCorrectionAbove100GeV_C%d",iCentrality));
    manualJetCorrectionRatioAbove100GeV[iCentrality] = (TH2D*) inputFile->Get(Form("manualJetCorrectionRatioAbove100GeV_C%d",iCentrality));
    energyInEtaStripFlipped[iCentrality] = (TH2D*) inputFile->Get(Form("energyInEtaStripFlipped_C%d",iCentrality));
    energyInConeFlipped[iCentrality] = (TH2D*) inputFile->Get(Form("energyInConeFlipped_C%d",iCentrality));
    manualJetCorrectionFlipped[iCentrality] = (TH2D*) inputFile->Get(Form("manualJetCorrectionFlipped_C%d",iCentrality));
    manualJetCorrectionRatioFlipped[iCentrality] = (TH2D*) inputFile->Get(Form("manualJetCorrectionRatioFlipped_C%d",iCentrality));
  }
  
//  currentProjection = energyInConeFlipped[0]->ProjectionY("lul",30,30);
//  cout << currentProjection->GetEntries() << endl;
//  currentProjection->Fit("gaus","","",10,80);
//  currentProjection->Draw();
//
//  return;
  
  // Draw the histograms to a canvas
  JDrawer *drawer = new JDrawer();
  
  // Change the margins better suited for 2D-drawing
  drawer->SetLeftMargin(0.15);
  drawer->SetRightMargin(0.13);
  drawer->SetBottomMargin(0.15);
  drawer->SetTopMargin(0.08);
  drawer->SetCanvasSize(700,650);
  
  TLine *line120 = new TLine(0,0,180,180);
  TLine *line160 = new TLine(0,0,220,220);
  TLine *line80 = new TLine(-100,-100,100,100);
  TLine *lineOne = new TLine(-1,-1,1,1);
  TLine *lineHalf = new TLine(-0.7,-0.7,0.7,0.7);
  
  TString centralityString[] = {"Cent: 0-10", "Cent: 10-30", "Cent: 30-50", "Cent: 50-90"};
  
  // Draw the histograms
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    if(drawAllJets){
      drawer->DrawHistogram(energyInEtaStrip[2][iCentrality], "Yield in reflected #eta strip [GeV]", "Yield in #eta strip [GeV]", Form("Inclusive, %s", centralityString[iCentrality].Data()), "colz");
      line120->Draw();
      
      if(saveFigures){
         gPad->GetCanvas()->SaveAs(Form("figures/energyInEtaStripInclusive%s_C%d.pdf", saveComment, iCentrality));
      }
      
      drawer->DrawHistogram(energyInCone[2][iCentrality], "Yield in reflected jet cone [GeV]", "Yield in jet cone [GeV]", Form("Inclusive, %s", centralityString[iCentrality].Data()), "colz");
      line160->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/energyInConeInclusive%s_C%d.pdf", saveComment, iCentrality));
      }
      
      drawer->DrawHistogram(manualJetCorrection[2][iCentrality], "Correction from reflected #eta strip [GeV]", "Correction from #eta strip [GeV]", Form("Inclusive, %s",centralityString[iCentrality].Data()), "colz");
      line80->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/manualJetCorrectionInclusive%s_C%d.pdf", saveComment, iCentrality));
      }
      
      drawer->DrawHistogram(manualJetCorrectionRatio[2][iCentrality], "Correction from reflected #eta strip [GeV]", "Correction from #eta strip [GeV]", Form("Inclusive, %s",centralityString[iCentrality].Data()), "colz");
      lineOne->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/manualJetCorrectionRatioInclusive%s_C%d.pdf", saveComment, iCentrality));
      }
      
    } // Drawing all jets
    
    if(drawLeadingJets){
      drawer->DrawHistogram(energyInEtaStrip[0][iCentrality], "Yield in reflected #eta strip [GeV]", "Yield in #eta strip [GeV]", Form("Leading jets, %s", centralityString[iCentrality].Data()), "colz");
      line120->Draw();
      
      if(saveFigures){
         gPad->GetCanvas()->SaveAs(Form("figures/energyInEtaStripLeading%s_C%d.pdf", saveComment, iCentrality));
      }
      
      drawer->DrawHistogram(energyInCone[0][iCentrality], "Yield in reflected jet cone [GeV]", "Yield in jet cone [GeV]", Form("Leading jets, %s", centralityString[iCentrality].Data()), "colz");
      line160->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/energyInConeLeading%s_C%d.pdf", saveComment, iCentrality));
      }
      
      drawer->DrawHistogram(manualJetCorrection[0][iCentrality], "Correction from reflected #eta strip [GeV]", "Correction from #eta strip [GeV]", Form("Leading jets, %s",centralityString[iCentrality].Data()), "colz");
      line80->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/manualJetCorrectionLeading%s_C%d.pdf", saveComment, iCentrality));
      }
      
      drawer->DrawHistogram(manualJetCorrectionRatio[0][iCentrality], "Correction from reflected #eta strip [GeV]", "Correction from #eta strip [GeV]", Form("Leading jets, %s",centralityString[iCentrality].Data()), "colz");
      lineOne->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/manualJetCorrectionRatioLeading%s_C%d.pdf", saveComment, iCentrality));
      }
      
    } // Drawing leading jets
    
    if(drawSubleadingJets){
      drawer->DrawHistogram(energyInEtaStrip[1][iCentrality], "Yield in reflected #eta strip [GeV]", "Yield in #eta strip [GeV]", Form("Subleading jets, %s", centralityString[iCentrality].Data()), "colz");
      line120->Draw();
      
      if(saveFigures){
         gPad->GetCanvas()->SaveAs(Form("figures/energyInEtaStripSubleading%s_C%d.pdf", saveComment, iCentrality));
      }
      
      drawer->DrawHistogram(energyInCone[1][iCentrality], "Yield in reflected jet cone [GeV]", "Yield in jet cone [GeV]", Form("Subleading jets, %s", centralityString[iCentrality].Data()), "colz");
      line160->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/energyInConeSubleading%s_C%d.pdf", saveComment, iCentrality));
      }
      
      drawer->DrawHistogram(manualJetCorrection[1][iCentrality], "Correction from reflected #eta strip [GeV]", "Correction from #eta strip [GeV]", Form("Subleading jets, %s",centralityString[iCentrality].Data()), "colz");
      line80->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/manualJetCorrectionSubleading%s_C%d.pdf", saveComment, iCentrality));
      }
      
      drawer->DrawHistogram(manualJetCorrectionRatio[1][iCentrality], "Correction from reflected #eta strip [GeV]", "Correction from #eta strip [GeV]", Form("Subleading jets, %s",centralityString[iCentrality].Data()), "colz");
      lineOne->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/manualJetCorrectionRatioSubleading%s_C%d.pdf", saveComment, iCentrality));
      }
      
    } // Drawing subleading jets
    
    if(drawJetsAbove100GeV){
      drawer->DrawHistogram(energyInEtaStripAbove100GeV[iCentrality], "Yield in reflected #eta strip [GeV]", "Yield in #eta strip [GeV]", Form("Jet p_{T} > 100 GeV, %s",centralityString[iCentrality].Data()), "colz");
      line120->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/energyInEtaStripAbove100GeV%s_C%d.pdf", saveComment, iCentrality));
      }
      
      drawer->DrawHistogram(energyInConeAbove100GeV[iCentrality], "Yield in reflected jet cone [GeV]", "Yield in jet cone [GeV]", Form("Jet p_{T} > 100 GeV, %s",centralityString[iCentrality].Data()), "colz");
      line160->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/energyInConeAbove100GeV%s_C%d.pdf", saveComment, iCentrality));
      }
      
      drawer->DrawHistogram(manualJetCorrectionAbove100GeV[iCentrality], "Correction from reflected #eta strip [GeV]", "Correction from #eta strip [GeV]", Form("Jet p_{T} > 100 GeV, %s",centralityString[iCentrality].Data()), "colz");
      line80->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/manualJetCorrectionAbove100GeV%s_C%d.pdf", saveComment, iCentrality));
      }
      
      drawer->DrawHistogram(manualJetCorrectionRatioAbove100GeV[iCentrality], "Correction from reflected #eta strip [GeV]", "Correction from #eta strip [GeV]", Form("Jet p_{T} > 100 GeV, %s",centralityString[iCentrality].Data()), "colz");
      lineHalf->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/manualJetCorrectionRatioAbove100GeV%s_C%d.pdf", saveComment, iCentrality));
      }
      
    } // Drawing jets above 100 GeV
    
    if(drawFlippedJets){
      drawer->DrawHistogram(energyInEtaStripFlipped[iCentrality], "Yield in reflected #eta strip [GeV]", "Yield in #eta strip [GeV]", Form("Jet p_{T} flips around 120 GeV, %s",centralityString[iCentrality].Data()), "colz");
      line120->Draw();
      
      if(saveFigures){
         gPad->GetCanvas()->SaveAs(Form("figures/energyInEtaStripFlipped%s_C%d.pdf", saveComment, iCentrality));
      }
      
      drawer->DrawHistogram(energyInConeFlipped[iCentrality], "Yield in reflected jet cone [GeV]", "Yield in jet cone [GeV]", Form("Jet p_{T} flips around 120 GeV, %s",centralityString[iCentrality].Data()), "colz");
      line160->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/energyInConeFlipped%s_C%d.pdf", saveComment, iCentrality));
      }
      
      drawer->DrawHistogram(manualJetCorrectionFlipped[iCentrality], "Correction from reflected #eta strip [GeV]", "Correction from #eta strip [GeV]", Form("Jet p_{T} flips around 120 GeV, %s",centralityString[iCentrality].Data()), "colz");
      line80->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/manualJetCorrectionFlipped%s_C%d.pdf", saveComment, iCentrality));
      }
      
      drawer->DrawHistogram(manualJetCorrectionRatioFlipped[iCentrality], "Correction from reflected #eta strip [GeV]", "Correction from #eta strip [GeV]", Form("Jet p_{T} flips around 120 GeV, %s",centralityString[iCentrality].Data()), "colz");
      lineHalf->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/manualJetCorrectionRatioFlipped%s_C%d.pdf", saveComment, iCentrality));
      }
      
    } // Drawing jets that flip over or under 120 GeV due to manual correction
    
  } // Centrality loop
  
  if(drawStripIllustration){
    drawer->SetRightMargin(0.04);
    drawer->SetLeftMargin(0.14);
    
    TH2D* emptyHistogram = new TH2D("empty", "empty", 1, -TMath::Pi()/2, 3*TMath::Pi()/2, 1, -2, 2);
    drawer->DrawHistogram(emptyHistogram, "#varphi", "#eta", "Excess yield determination", "colz");
    
    TLine *stripUp = new TLine(-TMath::Pi()/2, -0.4, 3*TMath::Pi()/2, -0.4);
    stripUp->SetLineColor(kRed);
    stripUp->SetLineWidth(2);
    stripUp->Draw();
    
    TLine *stripDown = new TLine(-TMath::Pi()/2, -0.8, 3*TMath::Pi()/2, -0.8);
    stripDown->SetLineColor(kRed);
    stripDown->SetLineWidth(2);
    stripDown->Draw();
    
    TLine *axisX1 = new TLine(-0.05, -0.65, 0.05, -0.55);
    axisX1->SetLineColor(kGreen+4);
    axisX1->SetLineWidth(2);
    axisX1->Draw();
    
    TLine *axisX2 = new TLine(-0.05, -0.55, 0.05, -0.65);
    axisX2->SetLineColor(kGreen+4);
    axisX2->SetLineWidth(2);
    axisX2->Draw();
    
    TEllipse *jetCone = new TEllipse(0, -0.6, 0.4);
    jetCone->SetLineColor(kRed);
    jetCone->SetLineWidth(2);
    jetCone->SetFillStyle(0);
    jetCone->Draw();
    
    TLine *reflectStripUp = new TLine(-TMath::Pi()/2, 0.4, TMath::Pi()/2, 0.4); // 3*TMath::Pi()/2
    reflectStripUp->SetLineColor(kBlue);
    reflectStripUp->SetLineWidth(2);
    reflectStripUp->Draw();
    
    TLine *reflectStripDown = new TLine(-TMath::Pi()/2, 0.8, TMath::Pi()/2, 0.8); // 3*TMath::Pi()/2
    reflectStripDown->SetLineColor(kBlue);
    reflectStripDown->SetLineWidth(2);
    reflectStripDown->Draw();
    
    TEllipse *reflectJetCone = new TEllipse(0, 0.6, 0.4);
    reflectJetCone->SetLineColor(kBlue);
    reflectJetCone->SetLineWidth(2);
    reflectJetCone->SetFillStyle(0);
    reflectJetCone->Draw();
    
    TLegend *upLegend = new TLegend(0.5,0.75,0.8,0.85);
    upLegend->SetFillStyle(0);upLegend->SetBorderSize(0);upLegend->SetTextSize(0.05);upLegend->SetTextFont(62);
    upLegend->AddEntry(reflectStripUp, "#eta reflect method", "l");
    upLegend->Draw();
    
    TLegend *downLegend = new TLegend(0.5,0.47,0.8,0.57);
    downLegend->SetFillStyle(0);downLegend->SetBorderSize(0);downLegend->SetTextSize(0.05);downLegend->SetTextFont(62);
    downLegend->AddEntry(stripUp, "#eta strip method", "l");
    downLegend->Draw();
    
    gPad->GetCanvas()->SaveAs("figures/manualJecEtaStripIllustration.pdf");
  } // Drawing illustrating stripes
  
  
  
  drawer->Reset();
  
  int lowBin, highBin, binCenter, lowFitRange, highFitRange;
  double gaussMean, gaussMeanError, gaussSigma, gaussSigmaError;
  int counter;
  TH2D *culpritHistogram;
  TString typeLegend[3] = {"Leading jets", "Subleading jets", "Jet p_{T} > 100 GeV"};
  int typeIndex = 0;
  
  // Fit Gaussians in slices to the eta strip histograms
  if(fitEtaStrip){
    TH1D *gaussianFitEtaStrip[nCentralityBins];
    double lowBinValueEtaStrip[] = {47,22,4,1};      // Neutral to charged
    double highBinValueEtaStrip[] = {80,44,18,16};   // Neutral to charged
    //double lowBinValueEtaStrip[] = {80,35,15,3};        // Neutral to scaled
    //double highBinValueEtaStrip[] = {130,75,30,10};     // Neutral to scaled
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      if(fitLeading){
        culpritHistogram = energyInEtaStrip[0][iCentrality];
        typeIndex = 0;
      } else if(fitSubleading){
        culpritHistogram = energyInEtaStrip[1][iCentrality];
        typeIndex = 1;
      } else {
        culpritHistogram = energyInEtaStripAbove100GeV[iCentrality];
        typeIndex = 2;
      }
      
      lowBin = culpritHistogram->GetXaxis()->FindBin(lowBinValueEtaStrip[iCentrality]+0.1);
      highBin = culpritHistogram->GetXaxis()->FindBin(highBinValueEtaStrip[iCentrality]-0.1);
      gaussianFitEtaStrip[iCentrality] = new TH1D(Form("etaStripGaussFit%d", iCentrality), Form("etaStripGaussFit%d", iCentrality), highBin-lowBin+1, lowBinValueEtaStrip[iCentrality], highBinValueEtaStrip[iCentrality]);
      gaussianFitEtaStrip[iCentrality]->Sumw2();
      counter = 0;
      
      for(int iBin = lowBin; iBin <= highBin; iBin++){
        binCenter = culpritHistogram->GetXaxis()->GetBinCenter(iBin);
        lowFitRange = binCenter - 20;
        highFitRange = binCenter + 80;
        std::tie(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError) = fitGauss(culpritHistogram->ProjectionY(Form("projection%d",iBin), iBin, iBin), lowFitRange, highFitRange/*, Form("llul%d", iBin), Form("llul%d",iBin)*/);
        //gaussianFitEtaStrip[iCentrality]->SetBinContent(counter,gaussMean/binCenter);
        //gaussianFitEtaStrip[iCentrality]->SetBinError(counter++,gaussMeanError/binCenter);
        gaussianFitEtaStrip[iCentrality]->SetBinContent(counter,gaussMean);
        gaussianFitEtaStrip[iCentrality]->SetBinError(counter++,gaussMeanError);
      }
      
      gaussianFitEtaStrip[iCentrality]->Fit("pol1","","",lowBinValueEtaStrip[iCentrality],highBinValueEtaStrip[iCentrality]);
      
    }
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      drawer->DrawHistogram(gaussianFitEtaStrip[iCentrality], "Yield in reflected #eta strip [GeV]", "#eta strip / reflection", Form("%s, %s", typeLegend[typeIndex].Data(), centralityString[iCentrality].Data()));
    }
  } // Fitting eta strip
  
  // Fit Gaussians in slices to the eta strip histograms
  if(fitJetCone){
    TH1D *gaussianFitJetCone[nCentralityBins];
    double lowBinValueJetCone[] = {45,20,4,3};
    double highBinValueJetCone[] = {75,40,16,22};
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      if(fitLeading){
        culpritHistogram = energyInCone[0][iCentrality];
        typeIndex = 0;
      } else if(fitSubleading){
        culpritHistogram = energyInCone[1][iCentrality];
        typeIndex = 1;
      } else {
        culpritHistogram = energyInConeAbove100GeV[iCentrality];
        typeIndex = 2;
      }
      
      lowBin = culpritHistogram->GetXaxis()->FindBin(lowBinValueJetCone[iCentrality]+0.1);
      highBin = culpritHistogram->GetXaxis()->FindBin(highBinValueJetCone[iCentrality]-0.1);
      gaussianFitJetCone[iCentrality] = new TH1D(Form("jetConeGaussFit%d", iCentrality), Form("jetConeGaussFit%d", iCentrality), highBin-lowBin+1, lowBinValueJetCone[iCentrality], highBinValueJetCone[iCentrality]);
      gaussianFitJetCone[iCentrality]->Sumw2();
      counter = 0;
      
      for(int iBin = lowBin; iBin <= highBin; iBin++){
        binCenter = culpritHistogram->GetXaxis()->GetBinCenter(iBin);
        lowFitRange = binCenter - binCenter*0.7;
        highFitRange = binCenter + 80;
        
        std::tie(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError) = fitGauss(culpritHistogram->ProjectionY(Form("projection%d",iBin), iBin, iBin), lowFitRange, highFitRange, Form("llul%d", iBin), Form("llul%d",iBin));
        //gaussianFitJetCone[iCentrality]->SetBinContent(counter,gaussMean/binCenter);
        //gaussianFitJetCone[iCentrality]->SetBinError(counter++,gaussMeanError/binCenter);
        gaussianFitJetCone[iCentrality]->SetBinContent(counter,gaussMean);
        gaussianFitJetCone[iCentrality]->SetBinError(counter++,gaussMeanError);
      }
      
      gaussianFitJetCone[iCentrality]->Fit("pol1","","",lowBinValueJetCone[iCentrality],highBinValueJetCone[iCentrality]);
      
    }
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      drawer->DrawHistogram(gaussianFitJetCone[iCentrality], "Yield in reflected jet cone [GeV]", "Peak position in jet cone [GeV]", Form("%s, %s", typeLegend[typeIndex].Data(), centralityString[iCentrality].Data()));
    }
  } // Fitting eta strip
  
  // Fit Gaussians in slices to the final correction histograms
  if(fitManualCorrection){
    TH1D *gaussianFitManualCorrection[nCentralityBins];
    double lowBinValueManualCorrection[] = {-40,-30,-20,-15};
    double highBinValueManualCorrection[] = {40,30,20,15};
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      culpritHistogram = manualJetCorrectionAbove100GeV[iCentrality];
      lowBin = culpritHistogram->GetXaxis()->FindBin(lowBinValueManualCorrection[iCentrality]+0.001);
      highBin =culpritHistogram->GetXaxis()->FindBin(highBinValueManualCorrection[iCentrality]-0.001);
      gaussianFitManualCorrection[iCentrality] = new TH1D(Form("manualCorrectionGaussFit%d", iCentrality), Form("manualCorrectionGaussFit%d", iCentrality), highBin-lowBin+1, lowBinValueManualCorrection[iCentrality], highBinValueManualCorrection[iCentrality]);
      gaussianFitManualCorrection[iCentrality]->Sumw2();
      counter = 0;
      
      for(int iBin = lowBin; iBin <= highBin; iBin++){
        binCenter = culpritHistogram->GetXaxis()->GetBinCenter(iBin);
        lowFitRange = binCenter - 45;
        highFitRange = binCenter + 45;
        std::tie(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError) = fitGauss(culpritHistogram->ProjectionY(Form("projection%d",iBin), iBin, iBin), lowFitRange, highFitRange/*, Form("llul%d", iBin), Form("llul%d",iBin)*/);
        gaussianFitManualCorrection[iCentrality]->SetBinContent(counter,gaussMean);
        gaussianFitManualCorrection[iCentrality]->SetBinError(counter++,gaussMeanError);
      }
      
      gaussianFitManualCorrection[iCentrality]->Fit("pol6", "", "", lowBinValueManualCorrection[iCentrality], highBinValueManualCorrection[iCentrality]);
      
    }
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      drawer->DrawHistogram(gaussianFitManualCorrection[iCentrality], "Manual jet energy correction [GeV]", "#eta strip / reflection", Form("Jet p_{T} > 120 GeV, %s", centralityString[iCentrality].Data()));
    }
  } // Fitting eta strip
  
  // Write output files
  if(outputFileName.EndsWith(".root")){
   
    TFile *outputFile = new TFile(outputFileName,"UPDATE");
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      //energyInEtaStrip[iCentrality]->Write("",TObject::kOverwrite);
      //energyInCone[iCentrality]->Write("",TObject::kOverwrite);
      //manualJetCorrection[iCentrality]->Write("",TObject::kOverwrite);
      //manualJetCorrectionRatio[iCentrality]->Write("",TObject::kOverwrite);
      //energyInEtaStripAbove100GeV[iCentrality]->Write("",TObject::kOverwrite);
      //energyInConeAbove100GeV[iCentrality]->Write("",TObject::kOverwrite);
      manualJetCorrectionAbove100GeV[iCentrality]->Write("",TObject::kOverwrite);
      //manualJetCorrectionRatioAbove100GeV[iCentrality]->Write("",TObject::kOverwrite);
      //energyInEtaStripFlipped[iCentrality]->Write("",TObject::kOverwrite);
      //energyInConeFlipped[iCentrality]->Write("",TObject::kOverwrite);
      //manualJetCorrectionFlipped[iCentrality]->Write("",TObject::kOverwrite);
      //manualJetCorrectionRatioFlipped[iCentrality]->Write("",TObject::kOverwrite);
    }
    
    outputFile->Close();
    
  }
}
