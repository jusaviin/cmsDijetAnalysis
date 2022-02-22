#include "DijetMethods.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

void fitJetEventPlaneVnMultiplicity(){

  // Open the data files
  const int nFiles = 2;
  TFile *inputFile[nFiles];
  inputFile[0] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multWeight_hadronEventPlane_jetEta1v3_2022-02-21.root");
  if(nFiles > 1) inputFile[1] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multWeight_allEventPlanes_jetEta1v3_2022-02-17.root");
  if(nFiles > 2) inputFile[2] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multWeight_forestEventPlanes_jetEta1v3_2022-02-18.root");
  if(nFiles > 3) inputFile[3] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multiplicityBins_manualEventPlaneV4_jetEta1v3_2022-01-26.root");
  // jetEventPlaneDeltaPhi_PbPb2018_caloJets_forestEventPlanes_jetEta1v3_2022-02-18.root
  // jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multWeight_allEventPlanes_jetEta1v3_2022-02-17.root
  // jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multWeight_forestEventPlanes_jetEta1v3_2022-02-18.root
  // jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multWeight_hadronEventPlane_jetEta1v3_2022-02-21.root
  // jetEventPlaneDeltaPhi_PbPb2018_pfCsJets_forestEventPlanes_jetEta1v6_2022-02-18.root
  // jetEventPlaneDeltaPhi_PbPbMC2018_pfCsJets_multWeight_allEventPlanes_jetEta1v6_2022-02-17.root
  // jetEventPlaneDeltaPhi_PbPbMC2018_pfCsJets_multWeight_forestEventPlanes_jetEta1v6_2022-02-18.root
  // jetEventPlaneDeltaPhi_PbPbMC2018_pfCsJets_multWeight_hadronEventPlane_jetEta1v6_2022-02-21.root
  
  TString jetTypeString[4] = {"Calo jets","PFCS jets","Calo tuned Q","Gen jet"};
  //TString jetTypeString[4] = {"Whole #eta","#eta strip","reflected #eta strip","Gen jet"};
  
  int eventPlaneOrder = 2; // The vn component that is plotted
  
  double centralityBinBorders[] = {0,10,30,50,90};
  
  bool matchYields = false;  // For comparison purposes, match the average yields of different jet collections
  int referenceYield = 0;   // Choose which jet collection to use as basis for yield matching
  
  bool drawEventPlaneDifference = false;  // Draw difference between manual event plane calculation to the one read from forest
  
  bool drawFits = false; // Draw the distributions and fits
  bool drawAllInSamePlot = false;  // True: Draw all three jet collection to the same plot. False: Use separate plots
  bool hideFit = false;
  
  bool printVs = false;    // True = Print the jet vn values to the console. False = Do not do that
  bool printSlides = true; // True = Print slides summarizing multiplicities and vns
  bool printSlidesUsingDifference = true;  // True = Use difference for delta. False = Use ratio for delta.
  bool smoothJetV2 = true;  // Smooth the v2 values for the jet with second order polynomial
  
  bool saveFigures = false;
  TString saveComment = "_jetV2Smoothing";
  TString figureFormat = "pdf";
  
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
  double multiplicityDeltas[nMultiplicityBins] = {0};
  double hadronV2deltas[nMultiplicityBins] = {0};
  double jetV2deltas[nMultiplicityBins] = {0};
  
  double xValues[16] = {501.018, 700.478, 901.016, 1101.64, 1301.44, 1501.33, 1701.07, 1900.85, 2101.17, 2301.97, 2500.19, 2700.93, 2900.64, 3093.29, 3278.72, 3465.86};
  double xErrors[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double yValues[nFiles][16];
  double yErrors[nFiles][16];
  double yValuesRatio[nFiles-1][16];
  double yErrorsRatio[nFiles-1][16];
  
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
  double errorTerm1, errorTerm2;
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iMultiplicity = firstPlottedBin; iMultiplicity <= lastPlottedBin; iMultiplicity++){
      yValues[iFile][iMultiplicity-firstPlottedBin] = fitFunctionMidRapidity[iFile][iMultiplicity]->GetParameter(eventPlaneOrder);
      yErrors[iFile][iMultiplicity-firstPlottedBin] = fitFunctionMidRapidity[iFile][iMultiplicity]->GetParError(eventPlaneOrder);
      
      // Print the numbers to the console
      cout << "Mult: " << multiplicityBinBorders[iMultiplicity] << "-" << multiplicityBinBorders[iMultiplicity+1] << "  v2: " << yValues[iFile][iMultiplicity-firstPlottedBin] << endl;
      
      if(iFile > 0){
        yValuesRatio[iFile-1][iMultiplicity-firstPlottedBin] = yValues[iFile][iMultiplicity-firstPlottedBin] / yValues[0][iMultiplicity-firstPlottedBin];
        errorTerm1 = yErrors[iFile][iMultiplicity-firstPlottedBin] / yValues[0][iMultiplicity-firstPlottedBin];
        errorTerm2 = yErrors[0][iMultiplicity-firstPlottedBin] * yValues[iFile][iMultiplicity-firstPlottedBin] / (yValues[0][iMultiplicity-firstPlottedBin] * yValues[0][iMultiplicity-firstPlottedBin]);
        yErrorsRatio[iFile-1][iMultiplicity-firstPlottedBin] = TMath::Sqrt(errorTerm1*errorTerm1 + errorTerm2*errorTerm2);
      }
      
    }
  }

  int markerColors[] = {kBlue, kRed, kGreen+3, kMagenta};
  int markerStyles[] = {kFullSquare, kFullCircle, kFullDiamond, kFullCross};

  // Draw a graph illustrating the multiplicity dependence of the jet v2
  TGraphErrors *multiplicityGraph[nFiles];
  TGraphErrors *multiplicityRatio[nFiles-1];
  for(int iFile = 0; iFile < nFiles; iFile++){
    multiplicityGraph[iFile] = new TGraphErrors(lastPlottedBin-firstPlottedBin+1, xValues, yValues[iFile], xErrors, yErrors[iFile]);
    multiplicityGraph[iFile]->SetMarkerStyle(markerStyles[iFile]);
    multiplicityGraph[iFile]->SetMarkerColor(markerColors[iFile]);
    
    if(iFile > 0){
      multiplicityRatio[iFile-1] = new TGraphErrors(lastPlottedBin-firstPlottedBin+1, xValues, yValuesRatio[iFile-1], xErrors, yErrorsRatio[iFile-1]);
      multiplicityRatio[iFile-1]->SetMarkerStyle(markerStyles[iFile]);
      multiplicityRatio[iFile-1]->SetMarkerColor(markerColors[iFile]);
    }
  }
  
  TLegend *legend;
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceGraph();
  drawer->SetTitleOffsetY(1.6);
  drawer->SetTitleOffsetX(1.15);
  
  double highYrange[] = {0.2,0.2,0.25,0.06,0.03,0.02,0.02};
  
  // Draw the multiplicity graphs
  drawer->DrawGraph(multiplicityGraph[0], 400, 3600, 0.0, highYrange[eventPlaneOrder], "Multiplicity", Form("Jet v_{%d}", eventPlaneOrder), " ", "p");
  for(int iFile = 1; iFile < nFiles; iFile++){
    multiplicityGraph[iFile]->Draw("same,p");
  }
  
  // Fit the graphs with pol2 function
  if(smoothJetV2){
    for(int iFile = 0; iFile < nFiles; iFile++){
      multiplicityGraph[iFile]->Fit("pol2");
      multiplicityGraph[iFile]->GetFunction("pol2")->SetLineColor(markerColors[iFile]);
    }
  }
  
  // Add a legend to the plot
  legend = new TLegend(0.3,0.2,0.5,0.4);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  //legend->SetHeader("PFCS jets");
  for(int iFile = 0; iFile < nFiles; iFile++){
    legend->AddEntry(multiplicityGraph[iFile], jetTypeString[iFile], "p");
  }
  legend->Draw();
  
  if(saveFigures){
    gPad->GetCanvas()->SaveAs(Form("figures/jetV%dvsMultiplicity%s.%s", eventPlaneOrder, saveComment.Data(), figureFormat.Data()));
  }
  
  // If there is a ratio file, draw it
  if(nFiles > 1){
    drawer->DrawGraph(multiplicityRatio[0], 400, 3600, 0.95, 1.05, "Multiplicity", Form("Jet v_{%d} ratio", eventPlaneOrder), " ", "p");
    for(int iFile = 2; iFile < nFiles; iFile++){
      multiplicityRatio[iFile-1]->Draw("same,p");
    }
  }
  
  // Add a legend to the ratio plot
  legend = new TLegend(0.23,0.2,0.43,0.4);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  for(int iFile = 1; iFile < nFiles; iFile++){
    legend->AddEntry(multiplicityRatio[iFile-1], Form("%s / %s", jetTypeString[iFile].Data(), jetTypeString[0].Data()), "p");
  }
  legend->Draw();
  
  if(saveFigures){
    gPad->GetCanvas()->SaveAs(Form("figures/jetV%dvsMultiplicityRatio%s.%s", eventPlaneOrder, saveComment.Data(), figureFormat.Data()));
  }
  
  
  // =============================================== //
  // Print slides summarizing multiplicities and v2s //
  // =============================================== //
  
  if(printSlides){
    
    // Print a slide with uncertainties for each source and each centrality for each V
    char namer[100];
    double tableMultiplicity, tableHadronV2, tableJetV2, previousMultiplicity, previousHadronV2, previousJetV2;
    double deltaMultiplicity, deltaHadronV2, deltaJetV2;
    double hadronV2scale, multiplicityScale;
    
    sprintf(namer,"\\frametitle{Multiplicity vs. $v_{%d}$}", eventPlaneOrder);
    
    cout << endl;
    cout << "\\begin{frame}" << endl;
    cout << namer << endl;
    cout << "\\small{" << endl;
    cout << "\\begin{center}" << endl;
    cout << "  \\begin{tabular}{cccccc}" << endl;
    cout << "    \\toprule" << endl;
    cout << "    Multiplicity & Hadron $v_{2}$ & Jet $v_{2}$ & $\\Delta(\\mathrm{Mult})$ & $\\Delta(v_{2}^{\\mathrm{hadron}})$ & $\\Delta(v_{2}^{\\mathrm{jet}})$ \\\\" << endl;
    cout << "    \\midrule" << endl;
    
    // Set the correct precision for printing floating point numbers
    cout << fixed << setprecision(4);
    
    for(int iMultiplicity = firstPlottedBin; iMultiplicity <= lastPlottedBin; iMultiplicity++){
      
      // Read the values from this line
      tableMultiplicity = xValues[iMultiplicity - firstPlottedBin];
      tableHadronV2 = yValues[0][iMultiplicity - firstPlottedBin];
      tableJetV2 = yValues[1][iMultiplicity - firstPlottedBin];
      if(smoothJetV2) tableJetV2 = multiplicityGraph[1]->GetFunction("pol2")->Eval(tableMultiplicity);
      
      // Read the values from previous line
      if(iMultiplicity > firstPlottedBin){
        previousMultiplicity = xValues[iMultiplicity - firstPlottedBin - 1];
        previousHadronV2 = yValues[0][iMultiplicity - firstPlottedBin - 1];
        previousJetV2 = yValues[1][iMultiplicity - firstPlottedBin - 1];
        if(smoothJetV2) previousJetV2 = multiplicityGraph[1]->GetFunction("pol2")->Eval(previousMultiplicity);
      } else {
        previousMultiplicity = tableMultiplicity;
        previousHadronV2 = tableHadronV2;
        previousJetV2 = tableJetV2;
      }
      
      // Calculate the delta to previous line, difference or ratio
      if(printSlidesUsingDifference){
        deltaMultiplicity = tableMultiplicity - previousMultiplicity;
        deltaHadronV2 = tableHadronV2 - previousHadronV2;
        deltaJetV2 = tableJetV2 - previousJetV2;
      } else {
        deltaMultiplicity = 1.0 - previousMultiplicity / tableMultiplicity;
        deltaHadronV2 = 1.0 - previousHadronV2 / tableHadronV2;
        deltaJetV2 = 1.0 - previousJetV2 / tableJetV2;
      }
      
      // Collect the deltas to an array for the next step
      multiplicityDeltas[iMultiplicity] = deltaMultiplicity;
      hadronV2deltas[iMultiplicity] = deltaHadronV2;
      jetV2deltas[iMultiplicity] = deltaJetV2;
      
      // Print this line to the console
      cout << "    $" << tableMultiplicity << "$";
      cout << " & $" << tableHadronV2 << "$";
      cout << " & $" << tableJetV2 << "$";
      cout << " & $" << deltaMultiplicity << "$";
      cout << " & $" << deltaHadronV2 << "$";
      cout << " & $" << deltaJetV2 << "$";
      cout << " \\\\" << endl;
      
    } // Multiplicity loop
    
    cout << "    \\bottomrule" << endl;
    cout << "  \\end{tabular}" << endl;
    cout << "\\end{center}" << endl;
    cout << "}" << endl;
    cout << "\\end{frame}" << endl;
    cout << endl;
    
    // After the first slide is printed, calculate scaling factors from different bins assuming linear dependency on hadron v2 and multiplicity
    // X Delta(hadron v2) + Y Delta(Multiplicity) = Delta(jet v2)
    // Use two adjacent lines in the table above to calculate this
    
    cout << endl;
    sprintf(namer,"\\frametitle{Multiplicity and hadron $v_{%d}$ scaling}", eventPlaneOrder);
    
    cout << endl;
    cout << "\\begin{frame}" << endl;
    cout << namer << endl;
    cout << "\\small{" << endl;
    cout << "\\begin{center}" << endl;
    cout << "  \\begin{tabular}{cccccc}" << endl;
    cout << "    \\toprule" << endl;
    cout << "    Multiplicity & $\\Delta(\\mathrm{Mult})$ & $\\Delta(v_{2}^{\\mathrm{hadron}})$ & $\\Delta(v_{2}^{\\mathrm{jet}})$ & Multi scale & $v_{2}^{\\mathrm{hadron}}$ scale  \\\\" << endl;
    cout << "    \\midrule" << endl;
    
    for(int iMultiplicity = firstPlottedBin; iMultiplicity <= lastPlottedBin; iMultiplicity++){
      
      // Read the values from this line
      tableMultiplicity = xValues[iMultiplicity - firstPlottedBin];
      
      // Calculate the scale from the deltas
      if(iMultiplicity > firstPlottedBin && iMultiplicity < lastPlottedBin){
        hadronV2scale = (jetV2deltas[iMultiplicity+1] - multiplicityDeltas[iMultiplicity+1] * jetV2deltas[iMultiplicity] / multiplicityDeltas[iMultiplicity]) / (hadronV2deltas[iMultiplicity+1] - multiplicityDeltas[iMultiplicity+1] * hadronV2deltas[iMultiplicity] / multiplicityDeltas[iMultiplicity]);
        multiplicityScale = (jetV2deltas[iMultiplicity] - hadronV2deltas[iMultiplicity] * hadronV2scale) / multiplicityDeltas[iMultiplicity];
      } else {
        hadronV2scale = 0;
        multiplicityScale = 0;
      }
      
      
      // Print this line to the console
      cout << "    $" << tableMultiplicity << "$";
      cout << " & $" << multiplicityDeltas[iMultiplicity] << "$";
      cout << " & $" << hadronV2deltas[iMultiplicity] << "$";
      cout << " & $" << jetV2deltas[iMultiplicity] << "$";
      cout << fixed << setprecision(6);
      cout << " & $" << multiplicityScale << "$";
      cout << fixed << setprecision(4);
      cout << " & $" << hadronV2scale << "$";
      cout << " \\\\" << endl;
      
    } // Multiplicity loop
    
    cout << "    \\bottomrule" << endl;
    cout << "  \\end{tabular}" << endl;
    cout << "\\end{center}" << endl;
    cout << "}" << endl;
    cout << "\\end{frame}" << endl;
    cout << endl;
    
  } // Printing slides about multiplicity vs. v2
  
  // Draw the jet-event plane fits
  drawer->Reset();
  
  
  TString multiplicityString;
  TString compactMultiplicityString;
  int colors[] = {kBlack, kRed, kBlue, kGreen+3};
  double maxYscale, minYscale;

  if(drawFits){
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
            gPad->GetCanvas()->SaveAs(Form("figures/jetEventPlaneDeltaPhi%s%s.%s", saveComment.Data(), compactMultiplicityString.Data(), figureFormat.Data()));
          }
          
        } else {
          legend->AddEntry(jetEventPlaneMidRapidity[iJetType][iMultiplicity], Form("%s, v_{2} = %.3f", jetTypeString[iJetType].Data(), fitFunctionMidRapidity[iJetType][iMultiplicity]->GetParameter(2)) ,"l");
        }
        
      } // Jet type loop
      if(drawAllInSamePlot){
        legend->Draw();
        
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/jetEventPlaneDeltaPhi%s_jetTypeComparison%s.%s", saveComment.Data(), compactMultiplicityString.Data(), figureFormat.Data()));
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
          gPad->GetCanvas()->SaveAs(Form("figures/jetEventPlaneDeltaPhiDifference%s%s.%s", saveComment.Data(), compactMultiplicityString.Data(), figureFormat.Data()));
        }
        
      }
    } // Multiplicity loop
  } // Drawing fits
}
