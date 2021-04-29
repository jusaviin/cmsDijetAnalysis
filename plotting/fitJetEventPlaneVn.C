#include "DijetMethods.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

void fitJetEventPlaneVn(){

  // Open the data files
  const int nFiles = 4;
  TFile *inputFile[nFiles];
  inputFile[0] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPb2018_caloJets_180GeVLeadingJet_dijetEvents_2021-04-28.root");
  if(nFiles > 1) inputFile[1] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbData2018_akCaloJet_eschemeAxis_2021-04-12.root");
  if(nFiles > 2) inputFile[2] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPb2018_flowJets_180GeVLeadingJet_dijetEvents_2021-04-28.root");
  if(nFiles > 3) inputFile[3] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbData2018_akFlowJet_wtaAxis_2021-04-12.root");
  
  TString jetTypeString[4] = {"Calo jets > 180 GeV","Calo jets > 120 GeV","PFCS flow jet > 180 GeV","PFCS flow jet > 120 GeV"};
  
  double centralityBinBorders[] = {0,10,30,50,90};
  
  bool matchYields = true;  // For comparison purposes, match the average yields of different jet collections
  int referenceYield = 2;   // Choose which jet collection to use as basis for yield matching
  
  bool drawAllInSamePlot = true;  // True: Draw all three jet collection to the same plot. False: Use separate plots
  bool hideFit = true;
  
  bool saveFigures = false;
  TString saveComment = "_genJets";
  
  // Read the histograms from the data files
  const int nCentralityBins = 4;
  const int nQvectorBins = 4;
  TH1D *jetEventPlaneMidRapidity[nFiles][nQvectorBins][nCentralityBins];
  TF1 *fitFunctionMidRapidity[nFiles][nQvectorBins][nCentralityBins];
  double averageYield[nFiles][nQvectorBins][nCentralityBins];
  double scaleFactor[nFiles][nQvectorBins][nCentralityBins];
  
  for(int iJetType = 0; iJetType < nFiles; iJetType++){
    for(int iQvector = 0; iQvector < 1; iQvector++){  // nQvectorBins
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        
        jetEventPlaneMidRapidity[iJetType][iQvector][iCentrality] = (TH1D*) inputFile[iJetType]->Get(Form("jetEventPlaneDeltaPhiForwardRap_Q%dC%d", iQvector, iCentrality));
        
      } // Centrality loop
    } // Q-vector loop
  } // Jet type loop
  
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
          
          scaleFactor[iJetType][iQvector][iCentrality] =  averageYield[referenceYield][iQvector][iCentrality] / averageYield[iJetType][iQvector][iCentrality];
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
  
//  cout << "Jet-event plane v2. " << jetTypeString[0].Data() << ". No Q-selection. Cent 0-10: " << fitFunctionMidRapidity[0][0][0]->GetParameter(2) << endl;
//  cout << "Jet-event plane v2. " << jetTypeString[0].Data() << ". No Q-selection. Cent 10-30: " << fitFunctionMidRapidity[0][0][1]->GetParameter(2) << endl;
//  cout << "Jet-event plane v2. " << jetTypeString[0].Data() << ". No Q-selection. Cent 30-50: " << fitFunctionMidRapidity[0][0][2]->GetParameter(2) << endl;
//  cout << "Jet-event plane v2. " << jetTypeString[1].Data() << ". No Q-selection. Cent 0-10: " << fitFunctionMidRapidity[1][0][0]->GetParameter(2) << endl;
//  cout << "Jet-event plane v2. " << jetTypeString[1].Data() << ". No Q-selection. Cent 10-30: " << fitFunctionMidRapidity[1][0][1]->GetParameter(2) << endl;
//  cout << "Jet-event plane v2. " << jetTypeString[1].Data() << ". No Q-selection. Cent 30-50: " << fitFunctionMidRapidity[1][0][2]->GetParameter(2) << endl;
//  cout << "Jet-event plane v2. " << jetTypeString[2].Data() << ". No Q-selection. Cent 0-10: " << fitFunctionMidRapidity[2][0][0]->GetParameter(2) << endl;
//  cout << "Jet-event plane v2. " << jetTypeString[2].Data() << ". No Q-selection. Cent 10-30: " << fitFunctionMidRapidity[2][0][1]->GetParameter(2) << endl;
//  cout << "Jet-event plane v2. " << jetTypeString[2].Data() << ". No Q-selection. Cent 30-50: " << fitFunctionMidRapidity[2][0][2]->GetParameter(2) << endl;
//
//  cout << "Jet-event plane v2. " << jetTypeString[0].Data() << ". Integral. Cent 0-10: " << jetEventPlaneMidRapidity[0][0][0]->Integral() << endl;
//  cout << "Jet-event plane v2. " << jetTypeString[0].Data() << ". Integral. Cent 10-30: " << jetEventPlaneMidRapidity[0][0][1]->Integral() << endl;
//  cout << "Jet-event plane v2. " << jetTypeString[0].Data() << ". Integral. Cent 30-50: " << jetEventPlaneMidRapidity[0][0][2]->Integral() << endl;
//  cout << "Jet-event plane v2. " << jetTypeString[1].Data() << ". Integral. Cent 0-10: " << jetEventPlaneMidRapidity[1][0][0]->Integral() << endl;
//  cout << "Jet-event plane v2. " << jetTypeString[1].Data() << ". Integral. Cent 10-30: " << jetEventPlaneMidRapidity[1][0][1]->Integral() << endl;
//  cout << "Jet-event plane v2. " << jetTypeString[1].Data() << ". Integral. Cent 30-50: " << jetEventPlaneMidRapidity[1][0][2]->Integral() << endl;
//  cout << "Jet-event plane v2. " << jetTypeString[2].Data() << ". Integral. Cent 0-10: " << jetEventPlaneMidRapidity[2][0][0]->Integral() << endl;
//  cout << "Jet-event plane v2. " << jetTypeString[2].Data() << ". Integral. Cent 10-30: " << jetEventPlaneMidRapidity[2][0][1]->Integral() << endl;
//  cout << "Jet-event plane v2. " << jetTypeString[2].Data() << ". Integral. Cent 30-50: " << jetEventPlaneMidRapidity[2][0][2]->Integral() << endl;
    
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
      legend = new TLegend(0.2,0.73,0.4,1);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*)0, centralityString, "");
    }
    
    maxYscale = jetEventPlaneMidRapidity[referenceYield][0][iCentrality]->GetMaximum();
    maxYscale = maxYscale + 0.1*maxYscale;
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
        legend->AddEntry((TObject*)0, centralityString, "");
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
  } // Centrality loop
}
