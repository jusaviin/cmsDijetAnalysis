#include "DijetMethods.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

void fitJetEventPlaneVnMultiplicity(){

  // Open the data files
  const int nFiles = 1;
  TFile *inputFile[nFiles];
  inputFile[0] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_subeNon0_multiplicityWeight_allEventPlanes_jetEta1v3_2022-01-30.root");
  if(nFiles > 1) inputFile[1] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multiplicityBins_manualEventPlane_jetEta1v3_2022-01-25.root");
  if(nFiles > 2) inputFile[2] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multiplicityBins_manualEventPlaneV3_jetEta1v3_2022-01-26.root");
  if(nFiles > 3) inputFile[3] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multiplicityBins_manualEventPlaneV4_jetEta1v3_2022-01-26.root");
  
  TString jetTypeString[4] = {"Calo jet","PFCS jet","Calo tuned Q","Gen jet"};
  //TString jetTypeString[4] = {"Whole #eta","#eta strip","reflected #eta strip","Gen jet"};
  
  int eventPlaneOrder = 2; // The vn component that is plotted
  
  double centralityBinBorders[] = {0,10,30,50,90};
  
  bool matchYields = false;  // For comparison purposes, match the average yields of different jet collections
  int referenceYield = 0;   // Choose which jet collection to use as basis for yield matching
  
  bool drawEventPlaneDifference = false;  // Draw difference between manual event plane calculation to the one read from forest
  
  bool drawAllInSamePlot = false;  // True: Draw all three jet collection to the same plot. False: Use separate plots
  bool hideFit = false;
  
  bool printVs = false; // True = Print the jet vn values to the console. False = Do not do that
  
  bool saveFigures = false;
  TString saveComment = "_eventPlane4";
  
  // Read the histograms from the data files
  const int nMultiplicityBins = 20;
  const int firstPlottedBin = 1;
  const int lastPlottedBin = 16;
  double multiplicityBinBorders[] = {0,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3400,3600,3800,4000,5000};
  TH1D *jetEventPlaneMidRapidity[nFiles][nMultiplicityBins];
  TF1 *fitFunctionMidRapidity[nFiles][nMultiplicityBins];
  TH1D *jetEventPlaneDifference[nMultiplicityBins];
  double averageYield[nFiles][nMultiplicityBins];
  double scaleFactor[nFiles][nMultiplicityBins];
  
  double xValues[16] = {500,700,900,1100,1300,1500,1700,1900,2100,2300,2500,2700,2900,3100,3300,3500};
  double xErrors[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double yValues[16];
  double yErrors[16];
  
  for(int iJetType = 0; iJetType < nFiles; iJetType++){
    for(int iMultiplicity = 0; iMultiplicity < nMultiplicityBins; iMultiplicity++){  // nQvectorBins
      
      jetEventPlaneMidRapidity[iJetType][iMultiplicity] = (TH1D*) inputFile[iJetType]->Get(Form("jetEventPlaneDeltaPhiManual_M%d", iMultiplicity));
      
      // If the histogram is not found, try different name
      if(jetEventPlaneMidRapidity[iJetType][iMultiplicity] == NULL){
        jetEventPlaneMidRapidity[iJetType][iMultiplicity] = (TH1D*) inputFile[iJetType]->Get(Form("jetEventPlaneOrder%d_M%d", eventPlaneOrder, iMultiplicity));
      }
    } // Multiplicity loop
  } // Jet type loop
  
  // Difference between manual and HiForest event planes
  if(drawEventPlaneDifference){
    for(int iMultiplicity = 0; iMultiplicity < nMultiplicityBins; iMultiplicity++){  // nQvectorBins
      
      jetEventPlaneDifference[iMultiplicity] = (TH1D*) inputFile[0]->Get(Form("jetEventPlaneDeltaPhiDiff_M%d", iMultiplicity));
      
      // If the histogram is not found, try different name
      if(jetEventPlaneDifference[iMultiplicity] == NULL){
        jetEventPlaneDifference[iMultiplicity] = (TH1D*) inputFile[0]->Get(Form("jetEventPlaneDifferenceOrder%d_M%d", eventPlaneOrder, iMultiplicity));
      }
      
    } // Multiplicity loop
  }
  
  // Scale the histograms such that yields between different jet collections match
  if(matchYields){
    
    // Find the average yield from each histogram
    for(int iJetType = 0; iJetType < nFiles; iJetType++){
      for(int iMultiplicity = 0; iMultiplicity < nMultiplicityBins; iMultiplicity++){
          
          jetEventPlaneMidRapidity[iJetType][iMultiplicity]->Fit("pol0");
          averageYield[iJetType][iMultiplicity] = jetEventPlaneMidRapidity[iJetType][iMultiplicity]->GetFunction("pol0")->GetParameter(0);
          
      } // Multiplicity loop
    } // Jet type loop
    
    // Scale the yields based on the basis
    
    for(int iJetType = 0; iJetType < nFiles; iJetType++){
      for(int iMultiplicity = 0; iMultiplicity < nMultiplicityBins; iMultiplicity++){
          
          scaleFactor[iJetType][iMultiplicity] =  1 / averageYield[iJetType][iMultiplicity];
          jetEventPlaneMidRapidity[iJetType][iMultiplicity]->Scale(scaleFactor[iJetType][iMultiplicity]);
          
      }
    }
    
  }
  
  
  // Do a fourier fit up to v4 to all histograms
  DijetMethods *fitter = new DijetMethods();
  for(int iJetType = 0; iJetType < nFiles; iJetType++){
    for(int iMultiplicity = 0; iMultiplicity < nMultiplicityBins; iMultiplicity++){
        
        fitter->FourierFit(jetEventPlaneMidRapidity[iJetType][iMultiplicity], 4);
        fitFunctionMidRapidity[iJetType][iMultiplicity] = jetEventPlaneMidRapidity[iJetType][iMultiplicity]->GetFunction("fourier");
        
    } // Multiplicity loop
  } // Jet type loop
  
  // Print values of vn components to console
  if(printVs){
    
    for(int iJetType = 0; iJetType < nFiles; iJetType++){
      
      for(int iFlow = 2; iFlow < 5; iFlow++){
        
        for(int iMultiplicity = 0; iMultiplicity < nMultiplicityBins; iMultiplicity++){
          
          cout << Form("Jet-event plane v%d. ", iFlow) << jetTypeString[iJetType].Data() << Form(". Multiplicity %.0f-%.0f: ", multiplicityBinBorders[iMultiplicity], multiplicityBinBorders[iMultiplicity+1]) << fitFunctionMidRapidity[iJetType][iMultiplicity]->GetParameter(iFlow) << endl;
          
        } // Centrality loop
      } // Flow component loop
      
    } // Jet type loop
    
  } // Printing values of vn components
    
  
  // Collect the values for the y-axis arrays to illustrate the multiplicity dependence
  for(int iMultiplicity = firstPlottedBin; iMultiplicity <= lastPlottedBin; iMultiplicity++){
    yValues[iMultiplicity-firstPlottedBin] = fitFunctionMidRapidity[0][iMultiplicity]->GetParameter(eventPlaneOrder);
    yErrors[iMultiplicity-firstPlottedBin] = fitFunctionMidRapidity[0][iMultiplicity]->GetParError(eventPlaneOrder);
  }
  
  
  // Draw a graph illustrating the multiplicity dependence of the jet v2
  TGraphErrors *multiplicityGraph = new TGraphErrors(lastPlottedBin-firstPlottedBin+1, xValues, yValues, xErrors, yErrors);
  multiplicityGraph->SetMarkerStyle(kFullSquare);
  multiplicityGraph->SetMarkerColor(kBlue);
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceGraph();
  
  double highYrange[] = {0.2,0.2,0.14,0.06,0.03,0.02,0.02};
  
  drawer->DrawGraph(multiplicityGraph, 400, 3600, 0.0, highYrange[eventPlaneOrder], "Multiplicity", Form("Jet v_{%d}", eventPlaneOrder), " ", "p");
  
  if(saveFigures){
    gPad->GetCanvas()->SaveAs(Form("figures/jetV%dvsMultiplicity%s.pdf", eventPlaneOrder, saveComment.Data()));
  }
  
  // Draw the jet-event plane fits
  drawer->Reset();

  TString multiplicityString;
  TString compactMultiplicityString;
  TLegend *legend;
  int colors[] = {kBlack, kRed, kBlue, kGreen+3};
  double maxYscale, minYscale;

  for(int iMultiplicity = firstPlottedBin; iMultiplicity <= lastPlottedBin; iMultiplicity++){

    multiplicityString = Form("Multiplicity: %.0f-%.0f%%", multiplicityBinBorders[iMultiplicity], multiplicityBinBorders[iMultiplicity+1]);
    compactMultiplicityString = Form("_M=%.0f-%.0f", multiplicityBinBorders[iMultiplicity], multiplicityBinBorders[iMultiplicity+1]);

    if(drawAllInSamePlot){
      legend = new TLegend(0.2,0.68,0.4,1);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*)0, Form("%s, manual jet energy correction",multiplicityString.Data()), "");
    }

    maxYscale = jetEventPlaneMidRapidity[referenceYield][iMultiplicity]->GetMaximum();
    maxYscale = maxYscale + 0.2*maxYscale; // 0.1
    minYscale = jetEventPlaneMidRapidity[referenceYield][iMultiplicity]->GetMinimum();
    minYscale = minYscale - 0.02*minYscale;

    for(int iJetType = 0; iJetType < nFiles; iJetType++){

      if(hideFit) jetEventPlaneMidRapidity[iJetType][iMultiplicity]->RecursiveRemove(fitFunctionMidRapidity[iJetType][iMultiplicity]);

      jetEventPlaneMidRapidity[iJetType][iMultiplicity]->Rebin(2);
      jetEventPlaneMidRapidity[iJetType][iMultiplicity]->Scale(0.5);
      jetEventPlaneMidRapidity[iJetType][iMultiplicity]->GetYaxis()->SetRangeUser(minYscale, maxYscale);
      if(drawAllInSamePlot) {
        jetEventPlaneMidRapidity[iJetType][iMultiplicity]->SetLineColor(colors[iJetType]);
        fitFunctionMidRapidity[iJetType][iMultiplicity]->SetLineColor(colors[iJetType]);

      }

      if(!drawAllInSamePlot || iJetType == 0){
      drawer->DrawHistogram(jetEventPlaneMidRapidity[iJetType][iMultiplicity], "#Delta#varphi", "A.U.", " ");
      } else {
        jetEventPlaneMidRapidity[iJetType][iMultiplicity]->Draw("same");
      }

      if(!drawAllInSamePlot){
        legend = new TLegend(0.2,0.7,0.4,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*)0, Form("%s, #Psi_{%d}",multiplicityString.Data(), eventPlaneOrder), "");
        legend->AddEntry((TObject*)0, jetTypeString[iJetType], "");
        legend->Draw();

        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/jetEventPlaneDeltaPhi%s%s.pdf", saveComment.Data(), compactMultiplicityString.Data()));
        }

      } else {
        legend->AddEntry(jetEventPlaneMidRapidity[iJetType][iMultiplicity], Form("%s, v_{2} = %.3f", jetTypeString[iJetType].Data(), fitFunctionMidRapidity[iJetType][iMultiplicity]->GetParameter(2)) ,"l");
      }

    } // Jet type loop
    if(drawAllInSamePlot){
      legend->Draw();

      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/jetEventPlaneDeltaPhi%s_jetTypeComparison%s.pdf", saveComment.Data(), compactMultiplicityString.Data()));
      }
    }

    if(drawEventPlaneDifference){

      legend = new TLegend(0.5,0.7,0.8,0.85);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*)0, multiplicityString, "");

      jetEventPlaneDifference[iMultiplicity]->Fit("gaus","","",-0.5,0.5);

      legend->AddEntry((TObject*)0, Form("#sigma = %.3f", jetEventPlaneDifference[iMultiplicity]->GetFunction("gaus")->GetParameter(2)), "");

      jetEventPlaneDifference[iMultiplicity]->SetLineColor(kBlack);
      jetEventPlaneDifference[iMultiplicity]->Rebin(2);
      jetEventPlaneDifference[iMultiplicity]->Scale(0.5);

      drawer->DrawHistogram(jetEventPlaneDifference[iMultiplicity], "#Delta#varphi", "A.U.", "Difference between manual and HiForest event plane");

      legend->Draw();

      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/jetEventPlaneDeltaPhiDifference%s%s.pdf", saveComment.Data(), compactMultiplicityString.Data()));
      }

    }
  } // Multiplicity loop
}
