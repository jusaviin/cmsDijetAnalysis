#include "JDrawer.h"

/*
 * Macro for plotting manual JEC debugging histograms
 */
void getAndPlotJecDebug(){

  // Open the data file
  TFile *inputFile = TFile::Open("data/PbPbMC2018_RecoGen_akPfCsJet_onlyJets_4pCentShift_matchedPtDifference_reflectCone_2021-08-12.root");
  // PbPbMC2018_RecoGen_akCaloJet_onlyJets_4pCentShift_matchedPtDifference_etaStrip_2021-08-12.root
  // PbPbMC2018_RecoGen_akCaloJet_onlyJets_4pCentShift_matchedPtDifference_reflectCone_2021-08-12.root
  // PbPbMC2018_RecoGen_akCaloJet_onlyJets_4pCentShift_matchedPtDifference_reflectStrip_2021-08-12.root
  // PbPbMC2018_RecoGen_akPfCsJet_onlyJets_4pCentShift_matchedPtDifference_etaStrip_2021-08-12.root
  // PbPbMC2018_RecoGen_akPfCsJet_onlyJets_4pCentShift_matchedPtDifference_reflectCone_2021-08-12.root
  // PbPbMC2018_RecoGen_akPfCsJet_onlyJets_4pCentShift_matchedPtDifference_reflectStrip_2021-08-12.root
  TString methodString = "#eta reflect"; // "#eta strip" "#eta reflect"
  TString methodSaveString = "etaReflect"; // "etaStrip" "etaReflect"
  TString jetTypeString = "PfCs jets";   // "PfCs Jets" "Calo jets"
  TString jetTypeSaveString = "pfCsJet";   // "pfCsJet" "calojet"
  
  // Figure saving flag
  const bool saveFigures = true;
  const bool drawStripIllustration = false;
  
  // Read the histograms from the file
  const int nCentralityBins = 4;
  TH2D *jecDebugHistogram[nCentralityBins];
  
  char newName[200];
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    sprintf(newName,"matchedPtExcessYield_C%d", iCentrality);
    
    jecDebugHistogram[iCentrality] = (TH2D*) inputFile->Get(newName);
  }
  
  // Draw the histograms to a canvas
  JDrawer *drawer = new JDrawer();
  
  // Change the margins better suited for 2D-drawing
  drawer->SetLeftMargin(0.15);
  drawer->SetRightMargin(0.13);
  drawer->SetBottomMargin(0.15);
  drawer->SetTopMargin(0.08);
  
  TString centralityString[] = {"Cent: 0-10", "Cent: 10-30", "Cent: 30-50", "Cent: 50-90"};
  
  // Draw the histograms
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    drawer->DrawHistogram(jecDebugHistogram[iCentrality], "Above average yield [GeV]", "(Reco p_{T} - Gen p_{T}) / Gen p_{T}", Form("%s, %s, %s", jetTypeString.Data(), methodString.Data(), centralityString[iCentrality].Data()), "colz");
    
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/jecDebug2D_%s_%s_C%d.pdf", jetTypeSaveString.Data(), methodSaveString.Data(), iCentrality));
    }
  }
  
  if(drawStripIllustration){
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
    
    TLine *reflectStripUp = new TLine(-TMath::Pi()/2, 0.4, 3*TMath::Pi()/2, 0.4);
    reflectStripUp->SetLineColor(kBlue);
    reflectStripUp->SetLineWidth(2);
    reflectStripUp->Draw();
    
    TLine *reflectStripDown = new TLine(-TMath::Pi()/2, 0.8, 3*TMath::Pi()/2, 0.8);
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
  }
  
}
