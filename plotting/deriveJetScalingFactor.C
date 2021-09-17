#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

/*
 * Combination of third order polynomial and exponential function
 * Functional form: f(x) = a + bx + cx^2 + dx^3 + h/s * e^(-0.5*((x-g)/s)^2)
 */
double jetPolyExp(double *x, double *par){
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]+(par[4]/par[5])*TMath::Exp(-0.5*TMath::Power((x[0]-par[6])/par[5],2));
}

/*
 * Macro for deriving weighting functions for MC to match vz and centrality in data
 */
void deriveJetScalingFactor(){

  // Open data and MC files for pp
  TFile *dataFile = TFile::Open("data/ppData2017_highForest_pfJets_fixedJEC_20EveMixed_wtaAxis_allCorrections_processed_2020-11-04.root");
  TFile *genFile = TFile::Open("data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_onlyJets_triggered_noJetWeight_JECv6_processed_2021-05-11.root");
  // ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_onlyJets_triggered_noJetWeight_JECv6_processed_2021-05-11.root
  // ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_allCorrections_tunedSeagull_processed_2019-10-22.root

  // Read the histograms for inclusive jet pT
  TH1D *hJetPt[2];

  DijetHistogramManager *dataManager = new DijetHistogramManager(dataFile);
  dataManager->SetLoadAllJets(true,true,true,false);
  dataManager->LoadProcessedHistograms();
  
  DijetHistogramManager *genManager = new DijetHistogramManager(genFile);
  genManager->SetLoadAllJets(true,true,true,false);
  genManager->LoadProcessedHistograms();
  
  hJetPt[0] = dataManager->GetHistogramJetPt(DijetHistogramManager::kLeadingJet, 0); // kAnyJet kLeadingJet kSubleadingJet
  hJetPt[1] = genManager->GetHistogramJetPt(DijetHistogramManager::kLeadingJet, 0);  // kAnyJet kLeadingJet kSubleadingJet
  
  // Normalize the histograms to the total yield
  for(int i = 0; i < 2; i++){
      hJetPt[i]->Scale(1.0 / hJetPt[i]->Integral( hJetPt[i]->FindBin(120), hJetPt[i]->GetNbinsX(), "width" ));
  }
  
  // Calculate a ratio between 2015 and 2018 data and fit it with high order polynomial
  TH1D *hRatio = (TH1D*) hJetPt[0]->Clone("ratioHistogram");
  
  //TF1 *ratioFit = new TF1("ratioFit",jetPolyExp,50,500,7);
  //ratioFit->SetParameters(0.723161,0.00236126,3.90984e-06,3.10631e-09,50,10,90);
  
  TF1 *ratioFit = new TF1("ratioFit","pol3",50,500);
  ratioFit->SetParameters(0.723161,0.00236126,-3.90984e-06,3.10631e-09);
  
  hRatio->Divide(hJetPt[1]);
  hRatio->Fit(ratioFit,"","",110,500);
  
  
  // Draw the ratios
  JDrawer *drawer = new JDrawer();
  drawer->DrawHistogram(hRatio,"Inclusive jet p_{T}","Data/MC", " ");
  
  /*
  //drawer->DrawHistogram(hCentrality[0],"centrality","counts", " ");
 // hCentrality[1]->SetLineColor(kRed);
 // hCentrality[1]->Draw("same");
  
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
  //drawer->DrawHistogram(hCentralityRatio,"centrality","Ratio fit");
  
  // Check with a small simulation that we regain data v_z and centrality
  double randomNumberVz, weigthVz;
  //double randomNumberCentrality, weigthCentrality;
  TH1D *vzReconstructed = (TH1D*) hVz[0]->Clone("vzClone");
  vzReconstructed->Reset();
//  TH1D *centralityReconstructed = (TH1D*) hCentrality[0]->Clone("centralityClone");
//  centralityReconstructed->Reset();
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
  centralityReconstructed->Draw("same");*/
}
