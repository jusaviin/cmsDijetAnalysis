#include "JDrawer.h"

/*
 * Macro for deriving weighting functions for MC to match vz and centrality in data
 */
void deriveMcWeightFunctions(){

  // Open files for MC and data
  TFile *dataFile = TFile::Open("data/dijetPbPb2018_akPu4CaloJets_onlyRegular_20eveMix_fixedJEC_eschemeAxis_noCorrections_processed_2021-02-16.root");
  //TFile *dataFile = TFile::Open("data/ppData2017_highForest_pfJets_onlyJets_L2rel_wtaAxis_processed_2019-08-05.root");
  TFile *mcFile = TFile::Open("data/dijetPbPb2018_akCsPfJet_minBias_onlyJets_2020-05-01.root");
  //TFile *mcFile = TFile::Open("data/dijet_ppMC_GenGen_Pythia8_pfJets_wtaAxis_onlyJets_processed_2019-08-06.root");
  // dijetPbPb2018_akCsPfJet_minBias_onlyJets_2020-05-01.root
  // PbPbMC_RecoReco_pfCsJets_onlyEventInfo_2019-07-22.root

  // Read the histograms for vz and centrality
  TH1D *hVz[2];
  TH1D *hCentrality[2];

  hVz[0] = (TH1D*) dataFile->Get("vertexZ");
  hCentrality[0] = (TH1D*) dataFile->Get("centrality");
  hVz[1] = (TH1D*) mcFile->Get("vertexZ");
  hCentrality[1] = (TH1D*) mcFile->Get("centrality");
  
  // Normalize all histograms to one
  for(int i = 0; i < 2; i++){
    hVz[i]->Scale(1/hVz[i]->Integral());
    hCentrality[i]->Scale(1/hCentrality[i]->Integral());
  }
  
  // Draw the histogram before weighting
  JDrawer *drawer = new JDrawer();
  
  drawer->DrawHistogram(hVz[0],"v_{z} (cm)","counts", " ");
  hVz[1]->SetLineColor(kRed);
  hVz[1]->Draw("same");

  drawer->DrawHistogram(hCentrality[0],"centrality","counts", " ");
  hCentrality[1]->SetLineColor(kRed);
  hCentrality[1]->Draw("same");
  
  // Do not draw over centrality
  drawer->CreateCanvas();
  
  // Calculate the ratio of the vz histograms and fit it with pol6
  TH1D *hVzRatio = (TH1D*) hVz[0]->Clone();
  hVzRatio->Divide(hVz[1]);
  hVzRatio->Fit("pol6","","",-15,15);
  TF1* vzFit = hVzRatio->GetFunction("pol6");
  
  // Calculate the ratio of the centrality histograms and fit it with a complicated function
  //TF1* fcent1= new TF1("fcent1","[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4+[7]*exp([5]+[6]*x)",0,90);
  //fcent1->SetParameters(4.40810, -7.75301e-02, 4.91953e-04, -1.34961e-06, 1.44407e-09, -160, 1, 3.68078e-15);
  TH1D *hCentralityRatio = (TH1D*) hCentrality[0]->Clone();
  hCentralityRatio->Divide(hCentrality[1]);
  hCentralityRatio->Fit("pol6","","",0,90);
  TF1 *centralityFit = hCentralityRatio->GetFunction("pol6");
  
  // Draw the fitted ratios
  drawer->DrawHistogram(hVzRatio,"v_{z} (cm)","Ratio fit");
  drawer->DrawHistogram(hCentralityRatio,"centrality","Ratio fit");
  
  // Check with a small simulation that we regain data v_z and centrality
  double randomNumberVz, weigthVz;
  double randomNumberCentrality, weigthCentrality;
  TH1D *vzReconstructed = (TH1D*) hVz[0]->Clone("vzClone");
  vzReconstructed->Reset();
  TH1D *centralityReconstructed = (TH1D*) hCentrality[0]->Clone("centralityClone");
  centralityReconstructed->Reset();
  for(int i = 0; i < 1000000; i++){
    randomNumberVz = hVz[1]->GetRandom();
    weigthVz = vzFit->Eval(randomNumberVz);
    vzReconstructed->Fill(randomNumberVz,weigthVz);
    
    randomNumberCentrality = hCentrality[1]->GetRandom();
    weigthCentrality = centralityFit->Eval(randomNumberCentrality);
    centralityReconstructed->Fill(randomNumberCentrality,weigthCentrality);
  }
  
  vzReconstructed->Scale(1.0/vzReconstructed->Integral());
  drawer->DrawHistogram(hVz[0],"v_{z} (cm)","counts", " ");
  vzReconstructed->SetLineColor(kMagenta);
  vzReconstructed->Draw("same");
  
  centralityReconstructed->Scale(1.0/centralityReconstructed->Integral());
  drawer->DrawHistogram(hCentrality[0],"centrality","counts", " ");
  centralityReconstructed->SetLineColor(kMagenta);
  centralityReconstructed->Draw("same");
}
