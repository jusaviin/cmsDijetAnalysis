#include "DijetMethods.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

void hadronV2ScalingFactorExtractor(){

  // Open the data files
  const int nFiles = 3;
  const int nCorrelationTypes = 2;  // 0 = Jet-event plane correlation. 1 = hadron-event plane correlation
  TFile *inputFile[nCorrelationTypes][nFiles];
  inputFile[0][0] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_updateMultWeight_jetEta1v3_2022-02-23.root");
  if(nFiles > 1) inputFile[0][1] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_qVectorAbove2_updateMultWeight_jetEta1v3_2022-02-23.root");
  if(nFiles > 2) inputFile[0][2] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_qVectorAbove2p5_updateMultWeight_jetEta1v3_2022-02-23.root");
  if(nFiles > 3) inputFile[0][3] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multiplicityBins_manualEventPlaneV4_jetEta1v3_2022-01-26.root");
  
  // jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_updateMultWeight_jetEta1v3_2022-02-23.root
  // jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_qVectorAbove2_updateMultWeight_jetEta1v3_2022-02-23.root
  // jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_qVectorAbove2p5_updateMultWeight_jetEta1v3_2022-02-23.root
  
  // jetEventPlaneDeltaPhi_PbPbMC2018_pfCsJets_updateMultWeight_jetEta1v6_2022-02-23.root
  // jetEventPlaneDeltaPhi_PbPbMC2018_pfCsJets_qVectorAbove2_updateMultWeight_jetEta1v6_2022-02-23.root
  // jetEventPlaneDeltaPhi_PbPbMC2018_pfCsJets_qVectorAbove2p5_updateMultWeight_jetEta1v6_2022-02-23.root
  
  inputFile[1][0] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multWeight_hadronEventPlane_jetEta1v3_2022-02-21.root");
  if(nFiles > 1) inputFile[1][1] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multWeight_hadronEventPlane_qVectorAbove2_jetEta1v3_2022-02-22.root");
  if(nFiles > 2) inputFile[1][2] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multWeight_hadronEventPlane_qVectorAbove2p5_jetEta1v3_2022-02-22.root");
  if(nFiles > 3) inputFile[1][3] = TFile::Open("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multiplicityBins_manualEventPlaneV4_jetEta1v3_2022-01-26.root");
  
  // jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multWeight_hadronEventPlane_jetEta1v3_2022-02-21.root
  // jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multWeight_hadronEventPlane_qVectorAbove2_jetEta1v3_2022-02-22.root
  // jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multWeight_hadronEventPlane_qVectorAbove2p5_jetEta1v3_2022-02-22.root
  
  // jetEventPlaneDeltaPhi_PbPbMC2018_pfCsJets_multWeight_hadronEventPlane_jetEta1v6_2022-02-21.root
  // jetEventPlaneDeltaPhi_PbPbMC2018_pfCsJets_multWeight_hadronEventPlane_qVectorAbove2_jetEta1v6_2022-02-22.root
  // jetEventPlaneDeltaPhi_PbPbMC2018_pfCsJets_multWeight_hadronEventPlane_qVectorAbove2p5_jetEta1v6_2022-02-22.root
  
  TString legendHeader = "Calo jets";
  TString jetTypeString[4] = {"All dijets","Q > 2","Q > 2.5","Gen jet"};
  //TString jetTypeString[4] = {"Whole #eta","#eta strip","reflected #eta strip","Gen jet"};
  
  int eventPlaneOrder = 2; // The vn component that is plotted
  
  double centralityBinBorders[] = {0,10,30,50,90};
  
  bool matchYields = false;  // For comparison purposes, match the average yields of different jet collections
  int referenceYield = 0;   // Choose which jet collection to use as basis for yield matching
  
  bool drawEventPlaneDifference = false;  // Draw difference between manual event plane calculation to the one read from forest
  
  bool drawFits = false; // Draw the distributions and fits
  bool drawAllInSamePlot = false;  // True: Draw all three jet collection to the same plot. False: Use separate plots
  bool hideFit = false;
  
  bool smoothJetV2 = false;  // Smooth the v2 values for the jet with second order polynomial
  bool fitDoubleRatio = true;  // Fit a pol0 to the double ratio
  
  bool saveFigures = true;
  TString saveComment = "_qVectorComparisonCaloJet";
  TString figureFormat = "pdf";
  
  // Read the histograms from the data files
  const int nMultiplicityBins = 20;
  const int firstPlottedBin = 1;
  const int lastPlottedBin = 16;
  double multiplicityBinBorders[] = {0,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3400,3600,3800,4000,5000};
  TH1D *jetEventPlaneCorrelation[nCorrelationTypes][nFiles][nMultiplicityBins];
  TF1 *jetEventPlaneFitFunction[nCorrelationTypes][nFiles][nMultiplicityBins];
  double averageYield[nCorrelationTypes][nFiles][nMultiplicityBins];
  double scaleFactor[nCorrelationTypes][nFiles][nMultiplicityBins];
  
  double xValues[16] = {501.6680, 700.6463, 901.0786, 1101.7435, 1301.0029, 1501.1398, 1701.5657, 1900.8402, 2101.0117, 2301.4069, 2500.8577, 2700.5630, 2899.4207, 3090.2760, 3277.1113, 3463.4320};
  double xErrors[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double yValues[nCorrelationTypes][nFiles][16];
  double yErrors[nCorrelationTypes][nFiles][16];
  double yValuesRatio[nCorrelationTypes][nFiles-1][16];
  double yErrorsRatio[nCorrelationTypes][nFiles-1][16];
  double yValuesDoubleRatio[nFiles-1][16];
  double yErrorsDoubleRatio[nFiles-1][16];
  
  for(int iCorrelationType = 0; iCorrelationType < nCorrelationTypes; iCorrelationType++){
    for(int iFile = 0; iFile < nFiles; iFile++){
      for(int iMultiplicity = 0; iMultiplicity < nMultiplicityBins; iMultiplicity++){  // nQvectorBins
        
        jetEventPlaneCorrelation[iCorrelationType][iFile][iMultiplicity] = (TH1D*) inputFile[iCorrelationType][iFile]->Get(Form("jetEventPlaneOrder%d_M%d", eventPlaneOrder, iMultiplicity));
        
        // If the histogram is not found, try different name
        if(jetEventPlaneCorrelation[iCorrelationType][iFile][iMultiplicity] == NULL){
          jetEventPlaneCorrelation[iCorrelationType][iFile][iMultiplicity] = (TH1D*) inputFile[iCorrelationType][iFile]->Get(Form("jetEventPlaneDeltaPhiManual_M%d", iMultiplicity));
        }
      } // Multiplicity loop
    } // File loop
  } // Correlation type loop
  
  // Scale the histograms such that yields between different jet collections match
  if(matchYields){
    
    // Find the average yield from each histogram
    for(int iCorrelationType = 0; iCorrelationType < nCorrelationTypes; iCorrelationType++){
      for(int iFile = 0; iFile < nFiles; iFile++){
        for(int iMultiplicity = 0; iMultiplicity < nMultiplicityBins; iMultiplicity++){
          
          jetEventPlaneCorrelation[iCorrelationType][iFile][iMultiplicity]->Fit("pol0");
          averageYield[iCorrelationType][iFile][iMultiplicity] = jetEventPlaneCorrelation[iCorrelationType][iFile][iMultiplicity]->GetFunction("pol0")->GetParameter(0);
          
        } // Multiplicity loop
      } // File loop
      
      
      // Scale the yields based on the basis
      
      for(int iFile = 0; iFile < nFiles; iFile++){
        for(int iMultiplicity = 0; iMultiplicity < nMultiplicityBins; iMultiplicity++){
          
          scaleFactor[iCorrelationType][iFile][iMultiplicity] =  1 / averageYield[iCorrelationType][iFile][iMultiplicity];
          jetEventPlaneCorrelation[iCorrelationType][iFile][iMultiplicity]->Scale(scaleFactor[iCorrelationType][iFile][iMultiplicity]);
          
        }
      }
    } // Correlation tpye loop
    
  } // Matching yields
  
  
  // Do a fourier fit up to v4 to all histograms
  DijetMethods *fitter = new DijetMethods();
  for(int iCorrelationType = 0; iCorrelationType < nCorrelationTypes; iCorrelationType++){
    for(int iFile = 0; iFile < nFiles; iFile++){
      for(int iMultiplicity = 0; iMultiplicity < nMultiplicityBins; iMultiplicity++){
        
        fitter->FourierFit(jetEventPlaneCorrelation[iCorrelationType][iFile][iMultiplicity], 4);
        jetEventPlaneFitFunction[iCorrelationType][iFile][iMultiplicity] = jetEventPlaneCorrelation[iCorrelationType][iFile][iMultiplicity]->GetFunction("fourier");
        
      } // Multiplicity loop
    } // File loop
  } // Correlation type loop
  
  // Collect the values for the y-axis arrays to illustrate the multiplicity dependence
  double errorTerm1, errorTerm2;
  for(int iCorrelationType = 0; iCorrelationType < nCorrelationTypes; iCorrelationType++){
    for(int iFile = 0; iFile < nFiles; iFile++){
      for(int iMultiplicity = firstPlottedBin; iMultiplicity <= lastPlottedBin; iMultiplicity++){
        yValues[iCorrelationType][iFile][iMultiplicity-firstPlottedBin] = jetEventPlaneFitFunction[iCorrelationType][iFile][iMultiplicity]->GetParameter(eventPlaneOrder);
        yErrors[iCorrelationType][iFile][iMultiplicity-firstPlottedBin] = jetEventPlaneFitFunction[iCorrelationType][iFile][iMultiplicity]->GetParError(eventPlaneOrder);
        
        if(iFile > 0){
          yValuesRatio[iCorrelationType][iFile-1][iMultiplicity-firstPlottedBin] = yValues[iCorrelationType][iFile][iMultiplicity-firstPlottedBin] / yValues[iCorrelationType][0][iMultiplicity-firstPlottedBin];
          errorTerm1 = yErrors[iCorrelationType][iFile][iMultiplicity-firstPlottedBin] / yValues[iCorrelationType][0][iMultiplicity-firstPlottedBin];
          errorTerm2 = yErrors[iCorrelationType][0][iMultiplicity-firstPlottedBin] * yValues[iCorrelationType][iFile][iMultiplicity-firstPlottedBin] / (yValues[iCorrelationType][0][iMultiplicity-firstPlottedBin] * yValues[iCorrelationType][0][iMultiplicity-firstPlottedBin]);
          yErrorsRatio[iCorrelationType][iFile-1][iMultiplicity-firstPlottedBin] = TMath::Sqrt(errorTerm1*errorTerm1 + errorTerm2*errorTerm2);
        }
        
      } // Multiplicity loop
    } // File loop
  } // Correlation type loop
  
  // Calculate the double ratio to study the hadron v2 scaling factors
  for(int iFile = 0; iFile < nFiles-1; iFile++){
    for(int iMultiplicity = firstPlottedBin; iMultiplicity <= lastPlottedBin; iMultiplicity++){
      yValuesDoubleRatio[iFile][iMultiplicity-firstPlottedBin] = yValuesRatio[1][iFile][iMultiplicity-firstPlottedBin] / yValuesRatio[0][iFile][iMultiplicity-firstPlottedBin];
      errorTerm1 = yErrorsRatio[1][iFile][iMultiplicity-firstPlottedBin] / yValuesRatio[0][iFile][iMultiplicity-firstPlottedBin];
      errorTerm2 = yErrorsRatio[0][iFile][iMultiplicity-firstPlottedBin] * yValuesRatio[1][iFile][iMultiplicity-firstPlottedBin] / (yValuesRatio[0][iFile][iMultiplicity-firstPlottedBin] * yValuesRatio[0][iFile][iMultiplicity-firstPlottedBin]);
      yErrorsDoubleRatio[iFile][iMultiplicity-firstPlottedBin] = TMath::Sqrt(errorTerm1*errorTerm1 + errorTerm2*errorTerm2);
    } // Multiplicity loop
  } // File loop

  int markerColors[] = {kBlue, kRed, kGreen+3, kMagenta};
  int markerStyles[] = {kFullSquare, kFullCircle, kFullDiamond, kFullCross};

  // Create the graphs for jet and hadron v2 as a function of multiplicity, together with ratios and double ratios
  TGraphErrors *multiplicityGraph[nCorrelationTypes][nFiles];
  TGraphErrors *multiplicityRatio[nCorrelationTypes][nFiles-1];
  TGraphErrors *multiplicityDoubleRatio[nFiles-1];
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iCorrelationType = 0; iCorrelationType < nCorrelationTypes; iCorrelationType++){
      multiplicityGraph[iCorrelationType][iFile] = new TGraphErrors(lastPlottedBin-firstPlottedBin+1, xValues, yValues[iCorrelationType][iFile], xErrors, yErrors[iCorrelationType][iFile]);
      multiplicityGraph[iCorrelationType][iFile]->SetMarkerStyle(markerStyles[iFile]);
      multiplicityGraph[iCorrelationType][iFile]->SetMarkerColor(markerColors[iFile]);
      
      if(iFile > 0){
        multiplicityRatio[iCorrelationType][iFile-1] = new TGraphErrors(lastPlottedBin-firstPlottedBin+1, xValues, yValuesRatio[iCorrelationType][iFile-1], xErrors, yErrorsRatio[iCorrelationType][iFile-1]);
        multiplicityRatio[iCorrelationType][iFile-1]->SetMarkerStyle(markerStyles[iFile]);
        multiplicityRatio[iCorrelationType][iFile-1]->SetMarkerColor(markerColors[iFile]);
      }
    }
    
    if(iFile > 0){
      multiplicityDoubleRatio[iFile-1] = new TGraphErrors(lastPlottedBin-firstPlottedBin+1, xValues, yValuesDoubleRatio[iFile-1], xErrors, yErrorsDoubleRatio[iFile-1]);
      multiplicityDoubleRatio[iFile-1]->SetMarkerStyle(markerStyles[iFile]);
      multiplicityDoubleRatio[iFile-1]->SetMarkerColor(markerColors[iFile]);
    }
  }
  
  TLegend *legend;
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceGraph();
  drawer->SetTitleOffsetY(1.6);
  drawer->SetTitleOffsetX(1.15);
  
  double highYrange[] = {0.2,0.2,0.25,0.06,0.03,0.02,0.02};
  TString yAxisName[] = {"Jet","Hadron"};
  TString fileNameAppend[] = {"jet","hadron"};
  
  // Draw the multiplicity graphs
  for(int iCorrelationType = 0; iCorrelationType < nCorrelationTypes; iCorrelationType++){
    drawer->DrawGraph(multiplicityGraph[iCorrelationType][0], 400, 3600, 0.0, highYrange[eventPlaneOrder], "Multiplicity", Form("%s v_{%d}",yAxisName[iCorrelationType].Data(), eventPlaneOrder), " ", "p");
    for(int iFile = 1; iFile < nFiles; iFile++){
      multiplicityGraph[iCorrelationType][iFile]->Draw("same,p");
    } // File loop
    
    // Add a legend to the plot
    legend = new TLegend(0.3,0.65,0.5,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->SetHeader(legendHeader);
    for(int iFile = 0; iFile < nFiles; iFile++){
      legend->AddEntry(multiplicityGraph[iCorrelationType][iFile], jetTypeString[iFile], "p");
    }
    legend->Draw();
    
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/%sV%dvsMultiplicity%s.%s", fileNameAppend[iCorrelationType].Data(), eventPlaneOrder, saveComment.Data(), figureFormat.Data()));
    }
    
  } // Correlation type loop
  
  // Fit the graphs with pol2 function
  if(smoothJetV2){
    for(int iCorrelationType = 0; iCorrelationType < nCorrelationTypes; iCorrelationType++){
      for(int iFile = 0; iFile < nFiles; iFile++){
        multiplicityGraph[iCorrelationType][iFile]->Fit("pol2");
        multiplicityGraph[iCorrelationType][iFile]->GetFunction("pol2")->SetLineColor(markerColors[iFile]);
      }
      
      // Calculate the ratio from smoothed values
      for(int iFile = 1; iFile < nFiles; iFile++){
        for(int iPoint = 0; iPoint < multiplicityRatio[iCorrelationType][iFile-1]->GetN(); iPoint++){
          multiplicityRatio[iCorrelationType][iFile-1]->SetPointY(iPoint, multiplicityGraph[iCorrelationType][iFile]->GetFunction("pol2")->Eval(xValues[iPoint]) / multiplicityGraph[iCorrelationType][iFile-1]->GetFunction("pol2")->Eval(xValues[iPoint]));
          multiplicityRatio[iCorrelationType][iFile-1]->SetPointError(iPoint,0,0.0001);
        }
      }
    }
  }
  
  // If there is a ratio file, draw it
  if(nFiles > 1){
    for(int iCorrelationType = 0; iCorrelationType < nCorrelationTypes; iCorrelationType++){
      drawer->DrawGraph(multiplicityRatio[iCorrelationType][0], 400, 3600, 0.95, 1.3, "Multiplicity", Form("%s v_{%d} ratio", yAxisName[iCorrelationType].Data(), eventPlaneOrder), " ", "p");
      
      for(int iFile = 2; iFile < nFiles; iFile++){
        multiplicityRatio[iCorrelationType][iFile-1]->Draw("same,p");
      }
      
      // Add a legend to the ratio plot
      legend = new TLegend(0.3,0.7,0.5,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->SetHeader(legendHeader);
      for(int iFile = 1; iFile < nFiles; iFile++){
        legend->AddEntry(multiplicityRatio[iCorrelationType][iFile-1], Form("%s / %s", jetTypeString[iFile].Data(), jetTypeString[0].Data()), "p");
      }
      legend->Draw();
      
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/%sV%dvsMultiplicityRatio%s.%s", fileNameAppend[iCorrelationType].Data(), eventPlaneOrder, saveComment.Data(), figureFormat.Data()));
      }
      
    }
  } // Drawing ratios
  
  // Draw the double ratios
  if(nFiles > 1){
    TLine *oneLine = new TLine(400,1,3600,1);
    oneLine->SetLineStyle(2);
    
    drawer->DrawGraph(multiplicityDoubleRatio[0], 400, 3600, 0.8, 1.2, "Multiplicity", Form("%s v_{%d} ratio / %s v_{%d} ratio", yAxisName[1].Data(), eventPlaneOrder, yAxisName[0].Data(), eventPlaneOrder), " ", "p");
    
    for(int iFile = 2; iFile < nFiles; iFile++){
      multiplicityDoubleRatio[iFile-1]->Draw("same,p");
    }
    
    // Add a legend to the ratio plot
    legend = new TLegend(0.3,0.7,0.5,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->SetHeader(legendHeader);
    for(int iFile = 1; iFile < nFiles; iFile++){
      legend->AddEntry(multiplicityDoubleRatio[iFile-1], Form("%s / %s", jetTypeString[iFile].Data(), jetTypeString[0].Data()), "p");
      if(fitDoubleRatio){
        multiplicityDoubleRatio[iFile-1]->Fit("pol0");
        multiplicityDoubleRatio[iFile-1]->GetFunction("pol0")->SetLineColor(markerColors[iFile]);
      }
    }
    legend->Draw();
    oneLine->Draw();
    
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/jetHadronV%dMultiplicityDoubleRatio%s.%s", eventPlaneOrder, saveComment.Data(), figureFormat.Data()));
    }
    
  } // Drawing ratios
  
  // Draw the jet-event plane fits
  drawer->Reset();
  
  
  TString multiplicityString;
  TString compactMultiplicityString;
  int colors[] = {kBlack, kRed, kBlue, kGreen+3};
  double maxYscale, minYscale;
  
  if(drawFits){
    for(int iCorrelationType = 0; iCorrelationType < nCorrelationTypes; iCorrelationType++){
      for(int iMultiplicity = firstPlottedBin; iMultiplicity <= lastPlottedBin; iMultiplicity++){
        
        multiplicityString = Form("Multiplicity: %.0f-%.0f%%", multiplicityBinBorders[iMultiplicity], multiplicityBinBorders[iMultiplicity+1]);
        compactMultiplicityString = Form("_M=%.0f-%.0f", multiplicityBinBorders[iMultiplicity], multiplicityBinBorders[iMultiplicity+1]);
        
        if(drawAllInSamePlot){
          legend = new TLegend(0.2,0.68,0.4,1);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->AddEntry((TObject*)0, Form("%s, manual jet energy correction",multiplicityString.Data()), "");
        }
        
        maxYscale = jetEventPlaneCorrelation[iCorrelationType][referenceYield][iMultiplicity]->GetMaximum();
        maxYscale = maxYscale + 0.2*maxYscale; // 0.1
        minYscale = jetEventPlaneCorrelation[iCorrelationType][referenceYield][iMultiplicity]->GetMinimum();
        minYscale = minYscale - 0.02*minYscale;
        
        for(int iFile = 0; iFile < nFiles; iFile++){
          
          if(hideFit) jetEventPlaneCorrelation[iCorrelationType][iFile][iMultiplicity]->RecursiveRemove(jetEventPlaneFitFunction[iCorrelationType][iFile][iMultiplicity]);
          
          jetEventPlaneCorrelation[iCorrelationType][iFile][iMultiplicity]->Rebin(2);
          jetEventPlaneCorrelation[iCorrelationType][iFile][iMultiplicity]->Scale(0.5);
          jetEventPlaneCorrelation[iCorrelationType][iFile][iMultiplicity]->GetYaxis()->SetRangeUser(minYscale, maxYscale);
          if(drawAllInSamePlot) {
            jetEventPlaneCorrelation[iCorrelationType][iFile][iMultiplicity]->SetLineColor(colors[iFile]);
            jetEventPlaneFitFunction[iCorrelationType][iFile][iMultiplicity]->SetLineColor(colors[iFile]);
            
          }
          
          if(!drawAllInSamePlot || iFile == 0){
            drawer->DrawHistogram(jetEventPlaneCorrelation[iCorrelationType][iFile][iMultiplicity], "#Delta#varphi", "A.U.", " ");
          } else {
            jetEventPlaneCorrelation[iCorrelationType][iFile][iMultiplicity]->Draw("same");
          }
          
          if(!drawAllInSamePlot){
            legend = new TLegend(0.2,0.7,0.4,0.9);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, Form("%s, #Psi_{%d}",multiplicityString.Data(), eventPlaneOrder), "");
            legend->AddEntry((TObject*)0, jetTypeString[iFile], "");
            legend->Draw();
            
            if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/%sEventPlaneDeltaPhi%s%s.%s", fileNameAppend[iCorrelationType].Data(), saveComment.Data(), compactMultiplicityString.Data(), figureFormat.Data()));
            }
            
          } else {
            legend->AddEntry(jetEventPlaneCorrelation[iCorrelationType][iFile][iMultiplicity], Form("%s, v_{2} = %.3f", jetTypeString[iFile].Data(), jetEventPlaneFitFunction[iCorrelationType][iFile][iMultiplicity]->GetParameter(2)) ,"l");
          }
          
        } // Jet type loop
        if(drawAllInSamePlot){
          legend->Draw();
          
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/%sEventPlaneDeltaPhi%s_jetTypeComparison%s.%s", fileNameAppend[iCorrelationType].Data(), saveComment.Data(), compactMultiplicityString.Data(), figureFormat.Data()));
          }
        }
      } // Multiplicity loop
    } // Correlation type loop
  } // Drawing fits
}
