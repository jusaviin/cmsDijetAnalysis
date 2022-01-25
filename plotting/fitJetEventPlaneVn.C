#include "DijetMethods.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

void fitJetEventPlaneVn(){

  // Open the data files
  const int nFiles = 3;
  TFile *inputFile[nFiles];
  inputFile[0] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_4pCentShift_manualEventPlane_jetEta1v3_2022-01-21.root");
  if(nFiles > 1) inputFile[1] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_pfCsJets_4pCentShift_manualEventPlane_jetEta1v3_2022-01-21.root");
  if(nFiles > 2) inputFile[2] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_4pCentShift_tunedQweights_manualEventPlane_jetEta1v3_2022-01-21.root");
  if(nFiles > 3) inputFile[3] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_genJets_manualEventPlane_dijetEvents_2021-05-05.root");
  
  TString jetTypeString[4] = {"Calo jet","PFCS jet","Calo tuned Q","Gen jet"};
  //TString jetTypeString[4] = {"Whole #eta","#eta strip","reflected #eta strip","Gen jet"};
  
  double centralityBinBorders[] = {0,10,30,50,90};
  
  bool matchYields = true;  // For comparison purposes, match the average yields of different jet collections
  int referenceYield = 1;   // Choose which jet collection to use as basis for yield matching
  
  bool drawEventPlaneDifference = false;  // Draw difference between manual event plane calculation to the one read from forest
  
  bool drawAllInSamePlot = true;  // True: Draw all three jet collection to the same plot. False: Use separate plots
  bool hideFit = false;
  
  bool printVs = true; // True = Print the jet vn values to the console. False = Do not do that
  
  bool saveFigures = false;
  TString saveComment = "_manualJECpfCs";
  
  // Read the histograms from the data files
  const int nCentralityBins = 4;
  const int nQvectorBins = 4;
  TH1D *jetEventPlaneMidRapidity[nFiles][nQvectorBins][nCentralityBins];
  TF1 *fitFunctionMidRapidity[nFiles][nQvectorBins][nCentralityBins];
  TH1D *jetEventPlaneDifference[nQvectorBins][nCentralityBins];
  double averageYield[nFiles][nQvectorBins][nCentralityBins];
  double scaleFactor[nFiles][nQvectorBins][nCentralityBins];
  
  for(int iJetType = 0; iJetType < nFiles; iJetType++){
    for(int iQvector = 0; iQvector < 1; iQvector++){  // nQvectorBins
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        
        jetEventPlaneMidRapidity[iJetType][iQvector][iCentrality] = (TH1D*) inputFile[iJetType]->Get(Form("jetEventPlaneDeltaPhiForwardRap_Q%dC%d", iQvector, iCentrality));
        
        if(jetEventPlaneMidRapidity[iJetType][iQvector][iCentrality] == NULL){
          jetEventPlaneMidRapidity[iJetType][iQvector][iCentrality] = (TH1D*) inputFile[iJetType]->Get(Form("jetEventPlaneDeltaPhiManual_Q%dC%d", iQvector, iCentrality));
        }
        
      } // Centrality loop
    } // Q-vector loop
  } // Jet type loop
  
  // Difference between manual and HiForest event planes
  if(drawEventPlaneDifference){
    for(int iQvector = 0; iQvector < 1; iQvector++){  // nQvectorBins
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        
        jetEventPlaneDifference[iQvector][iCentrality] = (TH1D*) inputFile[0]->Get(Form("jetEventPlaneDeltaPhiDiff_Q%dC%d", iQvector, iCentrality));
        
        
      } // Centrality loop
    } // Q-vector loop
  }
  
  // Scale the histograms such that yields between different jet collections match
  if(matchYields){
    
    // Find the average yield from each histogram
    for(int iJetType = 0; iJetType < nFiles; iJetType++){
      for(int iQvector = 0; iQvector < 1; iQvector++){ // nQvectorBins
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          
          jetEventPlaneMidRapidity[iJetType][iQvector][iCentrality]->Fit("pol0");
          averageYield[iJetType][iQvector][iCentrality] = jetEventPlaneMidRapidity[iJetType][iQvector][iCentrality]->GetFunction("pol0")->GetParameter(0);
          
        } // Centrality loop
      } // Q-vector loop
    } // Jet type loop
    
    // Scale the yields based on the basis
    
    for(int iJetType = 0; iJetType < nFiles; iJetType++){
      for(int iQvector = 0; iQvector < 1; iQvector++){ // nQvectorBins
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          
          scaleFactor[iJetType][iQvector][iCentrality] =  1 / averageYield[iJetType][iQvector][iCentrality];
          jetEventPlaneMidRapidity[iJetType][iQvector][iCentrality]->Scale(scaleFactor[iJetType][iQvector][iCentrality]);
          
        }
      }
    }
    
  }
  
  
  // Do a fourier fit up to v4 to all histograms
  DijetMethods *fitter = new DijetMethods();
  for(int iJetType = 0; iJetType < nFiles; iJetType++){
    for(int iQvector = 0; iQvector < 1; iQvector++){ // nQvectorBins
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        
        fitter->FourierFit(jetEventPlaneMidRapidity[iJetType][iQvector][iCentrality], 4);
        fitFunctionMidRapidity[iJetType][iQvector][iCentrality] = jetEventPlaneMidRapidity[iJetType][iQvector][iCentrality]->GetFunction("fourier");
        
      } // Centrality loop
    } // Q-vector loop
  } // Jet type loop
  
  // Print values of vn components to console
  if(printVs){
    
    for(int iJetType = 0; iJetType < nFiles; iJetType++){
      
      for(int iFlow = 2; iFlow < 5; iFlow++){
        
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          
          cout << Form("Jet-event plane v%d. ", iFlow) << jetTypeString[iJetType].Data() << Form(". Cent %.1f-%.1f: ", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]) << fitFunctionMidRapidity[iJetType][0][iCentrality]->GetParameter(iFlow) << endl;
          
        } // Centrality loop
      } // Flow component loop
      
    } // Jet type loop
    
  } // Printing values of vn components
    
  JDrawer *drawer = new JDrawer();
  
  TString centralityString;
  TString compactCentralityString;
  TLegend *legend;
  int colors[] = {kBlack, kRed, kBlue, kGreen+3};
  double maxYscale, minYscale;
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    centralityString = Form("Cent: %.0f-%.0f%%",centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]);
    compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]);
    
    if(drawAllInSamePlot){
      legend = new TLegend(0.2,0.68,0.4,1);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*)0, Form("%s, manual jet energy correction",centralityString.Data()), "");
    }
    
    maxYscale = jetEventPlaneMidRapidity[referenceYield][0][iCentrality]->GetMaximum();
    maxYscale = maxYscale + 0.2*maxYscale; // 0.1
    minYscale = jetEventPlaneMidRapidity[referenceYield][0][iCentrality]->GetMinimum();
    minYscale = minYscale - 0.02*minYscale;
    
    for(int iJetType = 0; iJetType < nFiles; iJetType++){
      
      if(hideFit) jetEventPlaneMidRapidity[iJetType][0][iCentrality]->RecursiveRemove(fitFunctionMidRapidity[iJetType][0][iCentrality]);
      
      jetEventPlaneMidRapidity[iJetType][0][iCentrality]->Rebin(2);
      jetEventPlaneMidRapidity[iJetType][0][iCentrality]->Scale(0.5);
      jetEventPlaneMidRapidity[iJetType][0][iCentrality]->GetYaxis()->SetRangeUser(minYscale, maxYscale);
      if(drawAllInSamePlot) {
        jetEventPlaneMidRapidity[iJetType][0][iCentrality]->SetLineColor(colors[iJetType]);
        fitFunctionMidRapidity[iJetType][0][iCentrality]->SetLineColor(colors[iJetType]);
        
      }
      
      if(!drawAllInSamePlot || iJetType == 0){
      drawer->DrawHistogram(jetEventPlaneMidRapidity[iJetType][0][iCentrality], "#Delta#varphi", "A.U.", " ");
      } else {
        jetEventPlaneMidRapidity[iJetType][0][iCentrality]->Draw("same");
      }
      
      if(!drawAllInSamePlot){
        legend = new TLegend(0.2,0.7,0.4,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*)0, Form("%s, Leading jet > 180 GeV",centralityString.Data()), "");
        legend->AddEntry((TObject*)0, jetTypeString[iJetType], "");
        legend->Draw();
        
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/jetEventPlaneDeltaPhi%s_%s%s.pdf", saveComment.Data(), jetTypeString[iJetType].Data(), compactCentralityString.Data()));
        }
        
      } else {
        legend->AddEntry(jetEventPlaneMidRapidity[iJetType][0][iCentrality], Form("%s, v_{2} = %.3f", jetTypeString[iJetType].Data(), fitFunctionMidRapidity[iJetType][0][iCentrality]->GetParameter(2)) ,"l");
      }
      
    } // Jet type loop
    if(drawAllInSamePlot){
      legend->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/jetEventPlaneDeltaPhi%s_jetTypeComparison%s.pdf", saveComment.Data(), compactCentralityString.Data()));
      }
    }
    
    if(drawEventPlaneDifference){
      
      legend = new TLegend(0.5,0.7,0.8,0.85);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*)0, centralityString, "");
      
      jetEventPlaneDifference[0][iCentrality]->Fit("gaus","","",-0.5,0.5);
      
      legend->AddEntry((TObject*)0, Form("#sigma = %.3f", jetEventPlaneDifference[0][iCentrality]->GetFunction("gaus")->GetParameter(2)), "");
      
      jetEventPlaneDifference[0][iCentrality]->SetLineColor(kBlack);
      jetEventPlaneDifference[0][iCentrality]->Rebin(2);
      jetEventPlaneDifference[0][iCentrality]->Scale(0.5);
      
      drawer->DrawHistogram(jetEventPlaneDifference[0][iCentrality], "#Delta#varphi", "A.U.", "Difference between manual and HiForest event plane");
      
      legend->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/jetEventPlaneDeltaPhiDifference%s%s.pdf", saveComment.Data(), compactCentralityString.Data()));
      }
      
    }
  } // Centrality loop
}
