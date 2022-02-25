#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

// Experimental fit function
double expFit(double *x, double *par){
  return par[0]+par[1]*TMath::Exp(x[0]*par[2]);
}

// Experimental fit function
double expWeight(double *x, double *par){
  return 0.0309424+0.110670*TMath::Exp(x[0]*(-0.00138397));
  
  // (0.0309424+0.110670*TMath::Exp(3000*(-0.00138397))) / (0.0685758+0.374178*TMath::Exp(x[0]*(-0.00189837)))
  // (0.0309424+0.110670*TMath::Exp(3000*(-0.00138397))) / (0.0309424+0.110670*TMath::Exp(x[0]*(-0.00138397)))
}

/*
 * Three part fit function to smooth the data
 */
double threePartFit(double *x, double *par){
  
  // Gauss function, if multiplicity > 3000
  //if(x[0] > 3000) return 503/(TMath::Sqrt(2*TMath::Pi())*215)*TMath::Exp(-0.5*TMath::Power((3000-x[0])/215,2));
  if(x[0] > 3000) return TMath::Exp(-0.5*TMath::Power((3000-x[0])/215,2));
  
  // Second order polynomial, if multiplicity < 300
  if(x[0] < 300) return (-10.5233 + 1.47193 * x[0] -0.0024 * x[0] * x[0]) / 503.0;
  
  return (0.10665541 * x[0] + 183.00338) / 503.0;
}

/*
 * Total weight function to match multiplicity in MC to that in data
 */
double totalMultiplicityWeight(double *x, double *par){
  
  double weight = 0;
  double base = 0;
  
  if(x[0] < 600){
    weight = ((0.0309424+0.110670*TMath::Exp(3000*(-0.00138397))) / (0.0685758+0.374178*TMath::Exp(x[0]*(-0.00550382)))) * 1.05;
  } else {
    weight = (0.0309424+0.110670*TMath::Exp(3000*(-0.00138397))) / (0.0309424+0.110670*TMath::Exp(x[0]*(-0.00138397)));
  }
  
  // Gauss function, if multiplicity > 3000
  if(x[0] > 3000) {
    base = TMath::Exp(-0.5*TMath::Power((3000-x[0])/215,2));
  } else if(x[0] < 300) {
    // Second order polynomial, if multiplicity < 300
    base = (-10.5233 + 1.47193 * x[0] -0.0024 * x[0] * x[0]) / 503.0;
  } else {
    base = (0.10665541 * x[0] + 183.00338) / 503.0;
  }
  
  return base * weight;
  
}

/*
 * Macro for studying the different components of the long range fits
 */
void multiplicityPlotter(){
  
  // File from which the long range distributions are plotted
  TString directoryName = "data/";
  const int nFiles = 3;
  TString fileName[] = {"PbPbMC2018_RecoGen_akPfCsJet_onlyJets_multWeight_eventPlanes_jetEta1v6_processed_2022-02-23.root", "PbPbMC2018_RecoGen_akPfCsJet_onlyJets_multWeight_qVectorAbove2_eventPlanes_jetEta1v6_processed_2022-02-23.root", "PbPbMC2018_RecoGen_akPfCsJet_onlyJets_multWeight_qVectorAbove2p5_eventPlanes_jetEta1v6_processed_2022-02-23.root"};
  // dijetPbPb2018_akCaloJet_onlyJets_jet100TriggerEta1v3_withTrackEff_processed_2022-01-19.root
  // PbPbMC2018_RecoGen_akCaloJet_onlyJets_noCentShift_jetEta1v3_reprocessed_2022-01-19.root
  // PbPbMC2018_RecoReco_akCaloJet_onlyJets_noCentShift_withTrackEff_jetEta1v3_processed_2022-01-20.root
  // dijetPbPb2018_akCaloJet_onlyJets_jet100TriggerEta1v3_aprilTrackEff_processed_2022-01-20.root
  // PbPbMC2018_RecoGen_akCaloJet_onlyJets_noCentShift_noCentWeight_jetEta1v3_processed_2022-01-21.root
  // PbPbMC2018_RecoGen_akCaloJet_onlyJets_noCentShift_multiplicityWeight_jetEta1v3_processed_2022-01-24.root
  
  // PbPbMC2018_RecoGen_akCaloJet_onlyJets_multWeight_eventPlanes_jetEta1v3_processed_2022-02-23.root
  // PbPbMC2018_RecoGen_akCaloJet_onlyJets_multWeight_qVectorAbove2_eventPlanes_jetEta1v3_processed_2022-02-23.root
  // PbPbMC2018_RecoGen_akCaloJet_onlyJets_multWeight_qVectorAbove2p5_eventPlanes_jetEta1v3_processed_2022-02-23.root
  
  // PbPbMC2018_RecoGen_akPfCsJet_onlyJets_multWeight_eventPlanes_jetEta1v6_processed_2022-02-23.root
  // PbPbMC2018_RecoGen_akPfCsJet_onlyJets_multWeight_qVectorAbove2_eventPlanes_jetEta1v6_processed_2022-02-23.root
  // PbPbMC2018_RecoGen_akPfCsJet_onlyJets_multWeight_qVectorAbove2p5_eventPlanes_jetEta1v6_processed_2022-02-23.root
  
  TString fileLegend[] = {"PFCS jets", "Q > 2", "Q > 2.5"};
  
  TString systemString = "Data";  // System to be put on a label
  
  // Open the files and create a histogram manager to reaf the histograms
  TFile *multiplicityFile[nFiles];
  DijetHistogramManager *manager[nFiles];
  for(int iFile = 0; iFile < nFiles; iFile++){
    multiplicityFile[iFile] = TFile::Open(directoryName+fileName[iFile]);
    manager[iFile] = new DijetHistogramManager(multiplicityFile[iFile]);
    manager[iFile]->SetLoadEventInformation(true);
    manager[iFile]->LoadProcessedHistograms();
  }
  
  // Binning information
  const int nCentralityBins = 4;
  const int nTrackPtBins = 4;
  const int nAsymmetryBins = 3;
  const int nFlowComponents = 4;
  const int firstDrawnCentralityBin = nCentralityBins;
  const int lastDrawnCentralityBin = nCentralityBins;
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  int maxCentralityBin = manager[0]->GetNCentralityBins();
  
  const bool drawMultiplicity = false;             // Draw multiplicity in each centrality bin
  const bool drawMultiplicityWeighted = false;     // Draw efficiency weighted multiplicity in each centrality bin
  const bool drawMultiplicityMap = false;          // Draw multiplicity vs. centrality map
  const bool drawMultiplicityMapWeighted = false;   // Draw efficiency weighted multiplicity vs. centrality map
  
  const bool drawMultiplicityDijet = true;             // Draw multiplicity in each centrality bin in dijet events
  const bool drawMultiplicityWeightedDijet = false;     // Draw efficiency weighted multiplicity in each centrality bin in dijet events
  const bool drawMultiplicityMapDijet = false;          // Draw multiplicity vs. centrality map in dijet events
  const bool drawMultiplicityMapWeightedDijet = false;   // Draw efficiency weighted multiplicity vs. centrality map in dijet events
  
  const bool drawAllToSameFigure = false;
  
  const bool printMeanSlide = true;
  
  const bool saveFigures = false;
  TString saveComment = "_qComparisonPfCs";
  TString figureFormat = "pdf";
  
  // Define histograms
  TH1D *multiplicity[nFiles][4][nCentralityBins+1];  // Second bin: multiplicity, weighted, dijet, dijet weighted
  TH2D *multiplicityMap[nFiles][4];                  // Second bin: multiplicity, weighted, dijet, dijet weighted
  
  // Read the histograms from the histogram manager
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      multiplicity[iFile][0][iCentrality] = manager[iFile]->GetHistogramMultiplicity(iCentrality);
      multiplicity[iFile][1][iCentrality] = manager[iFile]->GetHistogramMultiplicityWeighted(iCentrality);
      multiplicity[iFile][2][iCentrality] = manager[iFile]->GetHistogramMultiplicityDijet(iCentrality);
      multiplicity[iFile][3][iCentrality] = manager[iFile]->GetHistogramMultiplicityDijetWeighted(iCentrality);
    } // centrality loop
    
    // Load multiplicity histograms without centrality binning
    multiplicity[iFile][0][nCentralityBins] = manager[iFile]->GetHistogramMultiplicity(maxCentralityBin);
    multiplicity[iFile][1][nCentralityBins] = manager[iFile]->GetHistogramMultiplicityWeighted(maxCentralityBin);
    multiplicity[iFile][2][nCentralityBins] = manager[iFile]->GetHistogramMultiplicityDijet(maxCentralityBin);
    multiplicity[iFile][3][nCentralityBins] = manager[iFile]->GetHistogramMultiplicityDijetWeighted(maxCentralityBin);
    
    // Load multiplicity map histograms
    multiplicityMap[iFile][0] = manager[iFile]->GetHistogramMultiplicityMap();
    multiplicityMap[iFile][1] = manager[iFile]->GetHistogramWeightedMultiplicityMap();
    multiplicityMap[iFile][2] = manager[iFile]->GetHistogramMultiplicityMapDijet();
    multiplicityMap[iFile][3] = manager[iFile]->GetHistogramWeightedMultiplicityMapDijet();
  } // File loop
  
  TString saveString[] = {"", "Weighted", "Dijet", "WeightedDijet"};
  TString legendString[] = {"", "", "PFCS dijet events", "PFCS dijet events"};
  TString centralityString[] = {"C = 0-10", "C = 10-30", "C = 30-50", ""};
  TString saveCentralityString[] = {"C=0-10", "C=10-30", "C=30-50", ""};
  bool drawMultiplicityArray[] = {drawMultiplicity, drawMultiplicityWeighted, drawMultiplicityDijet, drawMultiplicityWeightedDijet};
  bool drawMultiplicityMapArray[] = {drawMultiplicityMap, drawMultiplicityMapWeighted, drawMultiplicityMapDijet, drawMultiplicityMapWeightedDijet};
  
  // Configure the histogram drawing class
  JDrawer *drawer = new JDrawer();
  
  double maxYscale, minYscale, yDifference;
  TLegend *legend;
  TString saveName;
  
  const int nMultiplicityBins = 20;
  const int firstPlottedBin = 1;
  const int lastPlottedBin = 16;
  double multiplicityBinBorders[] = {0,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3400,3600,3800,4000,5000};
  double deltaQone, deltaQtwo;
  int lowBin, highBin;
  int fileColors[] = {kBlue, kRed, kGreen+3, kMagenta};
  
  // Define array for mean values within bins
  double meanArray[nFiles][4][nMultiplicityBins] = {{{0}}};
  
  // Draw the distributions
  for(int iType = 0; iType < 4; iType++){
    
    if(!drawMultiplicityArray[iType]) continue;
    
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      
      // Draw the histogram to the canvas
      multiplicity[0][iType][iCentrality]->SetLineColor(fileColors[0]);
      drawer->DrawHistogram(multiplicity[0][iType][iCentrality], "Track multiplicity", "counts", " ");
      
      for(int iFile = 1; iFile < nFiles; iFile++){
        multiplicity[iFile][iType][iCentrality]->SetLineColor(fileColors[iFile]);
        multiplicity[iFile][iType][iCentrality]->Draw("same");
      }
      
      // Draw a legend for the histogram
      legend = new TLegend(0.17,0.7,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      
      if(iCentrality < nCentralityBins) legend->AddEntry((TObject*) 0, Form("C = %.0f-%.0f", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "");
      if(legendString[iType] != "") legend->AddEntry((TObject*) 0, legendString[iType], "");
      for(int iFile = 1; iFile < nFiles; iFile++){
        legend->AddEntry(multiplicity[iFile][iType][iCentrality], fileLegend[iFile], "l");
      }
      legend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        if(iCentrality == nCentralityBins){
          saveName = Form("trackMultiplicity%s%s", saveString[iType].Data(), saveComment.Data());
        } else {
          saveName = Form("trackMultiplicity%s%sC=%.0f-%.0f", saveString[iType].Data(), saveComment.Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]);
        }
        gPad->GetCanvas()->SaveAs(Form("figures/%s.%s", saveName.Data(), figureFormat.Data()));
      } // Saving figures
      
      if(iCentrality == nCentralityBins){
        
        for(int iFile = 0; iFile < nFiles; iFile++){
          for(int iMultiplicity = 0; iMultiplicity < nMultiplicityBins; iMultiplicity++){
            lowBin = multiplicity[iFile][iType][iCentrality]->GetXaxis()->FindBin(multiplicityBinBorders[iMultiplicity]+0.1);
            highBin = multiplicity[iFile][iType][iCentrality]->GetXaxis()->FindBin(multiplicityBinBorders[iMultiplicity+1]-0.1);
            multiplicity[iFile][iType][iCentrality]->GetXaxis()->SetRange(lowBin,highBin);
            meanArray[iFile][iType][iMultiplicity] = multiplicity[iFile][iType][iCentrality]->GetMean();
          }
          multiplicity[iFile][iType][iCentrality]->GetXaxis()->SetRange(0,0);
        }
        
        if(printMeanSlide){
                    
          cout << endl;
          cout << "\\begin{frame}" << endl;
          cout << "\\frametitle{Mean multiplicity values}" << endl;
          cout << "\\small{" << endl;
          cout << "\\begin{center}" << endl;
          cout << "  \\begin{tabular}{cccccc}" << endl;
          cout << "    \\toprule" << endl;
          cout << "    Bin & $\\langle M \\rangle$ & $Q > 2$ & $Q > 2.5$ & $\\Delta(Q > 2)$ & $\\Delta(Q > 2.5)$ \\\\" << endl;
          cout << "    \\midrule" << endl;
          
          // Set the correct precision for printing floating point numbers
          
          for(int iMultiplicity = firstPlottedBin; iMultiplicity <= lastPlottedBin; iMultiplicity++){
            
            // Calculate the deltas from default
            deltaQone = meanArray[1][iType][iMultiplicity] / meanArray[0][iType][iMultiplicity] - 1;
            deltaQtwo = meanArray[2][iType][iMultiplicity] / meanArray[0][iType][iMultiplicity] - 1;
            
            // Print this line to the console
            cout << fixed << setprecision(0);
            cout << "    $" << multiplicityBinBorders[iMultiplicity] << "-" << multiplicityBinBorders[iMultiplicity+1] << "$";
            cout << fixed << setprecision(4);
            cout << " & $" << meanArray[0][iType][iMultiplicity] << "$";
            cout << " & $" << meanArray[1][iType][iMultiplicity] << "$";
            cout << " & $" << meanArray[2][iType][iMultiplicity] << "$";
            cout << fixed << setprecision(4);
            cout << " & $" << deltaQone*100 << "\\%" << "$";
            cout << " & $" << deltaQtwo*100 << "\\%" << "$";
            cout << " \\\\" << endl;
            
          } // Multiplicity loop
          
          cout << "    \\bottomrule" << endl;
          cout << "  \\end{tabular}" << endl;
          cout << "\\end{center}" << endl;
          cout << "}" << endl;
          cout << "\\end{frame}" << endl;
          cout << endl;
        } // Printing slide for mean values
        
      } // Stuff for the integrated centrality bin
    } // Centrality loop
  } // Drawing track multiplicity
  
  // Draw all the centrality bins to the same figure
  if(drawAllToSameFigure){
    int centralityColors[] = {kBlue, kRed, kGreen+3, kCyan, kBlack};
    for(int iType = 0; iType < 4; iType++){

      if(!drawMultiplicityArray[iType]) continue;

      for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){

        multiplicity[0][iType][iCentrality]->SetLineColor(centralityColors[iCentrality]);

        // Draw the histogram to the canvas
        if(iCentrality == 0){
          drawer->DrawHistogram(multiplicity[0][iType][iCentrality], "Track multiplicity", "counts", " ");
          // Draw a legend for the histogram
          legend = new TLegend(0.17,0.6,0.37,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          if(legendString[iType] != "") legend->AddEntry((TObject*) 0, legendString[iType], "");
        } else {
          multiplicity[0][iType][iCentrality]->Draw("same");
        }

        if(iCentrality < nCentralityBins){
          legend->AddEntry(multiplicity[0][iType][iCentrality], Form("C = %.0f-%.0f", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "l");
        } else {
          legend->AddEntry(multiplicity[0][iType][iCentrality], "Total", "l");
        }
        
      } // Centrality loop

      legend->Draw();

      // Save the figures to file
      if(saveFigures){
        saveName = Form("trackMultiplicity%sAllCentralities%s", saveString[iType].Data(), saveComment.Data());
        gPad->GetCanvas()->SaveAs(Form("figures/%s.%s", saveName.Data(), figureFormat.Data()));
      } // Saving figures
      
    } // Drawing track multiplicity
  }
  
  // Change the left margin better suited for 2D-drawing
  drawer->SetLeftMargin(0.12);
  drawer->SetRightMargin(0.14);
  drawer->SetBottomMargin(0.14);
  drawer->SetTopMargin(0.06);
  drawer->SetTitleOffsetY(0.8);
  
  for(int iType = 0; iType < 4; iType++){
    
    if(!drawMultiplicityMapArray[iType]) continue;
    
    // Possibility to zoom to axes
    multiplicityMap[0][iType]->GetYaxis()->SetRangeUser(0,20);
    multiplicityMap[0][iType]->GetXaxis()->SetRangeUser(0,4000);
    
    // Draw the histogram to the canvas
    drawer->DrawHistogram(multiplicityMap[0][iType], "Track multiplicity", "Centrality", " ", "colz");
    
    // Draw a legend for the histogram
    legend = new TLegend(0.67,0.7,0.87,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    
    legend->AddEntry((TObject*) 0, systemString.Data(), "");
    if(legendString[iType] != "") legend->AddEntry((TObject*) 0, legendString[iType], "");
    legend->Draw();
    
    // Save the figures to file
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/multiplicityMap%s%s.%s", saveString[iType].Data(), saveComment.Data(), figureFormat.Data()));
    } // Saving figures
    
  } // Multiplicity vs. centrality map
  
  // Try to fit an exponential function
  
//  TF1 *myExpFit = new TF1("expFit",expFit,200,3000,3);
//  myExpFit->SetParameters(1,1,0);
//
//  drawer->CreateCanvas();
//  multiplicity[3][nCentralityBins]->Fit(myExpFit,"","",200,600);
//
//  drawer->CreateCanvas();
//  TF1* myFit = new TF1("myFit",threePartFit,0,4000,0);
//  myFit->Draw();
//
  drawer->CreateCanvas();
  TF1* myWeight = new TF1("myWeight",totalMultiplicityWeight,0,4000,0);
  myWeight->Draw();
  
  TF1* myExpWeightF = new TF1("myExpWeight",expWeight,0,4000,0);
  cout << "600 match: " << myExpWeightF->Eval(600) << endl;
 
  // Test the weight function
  TH1D *testMultiplicity = (TH1D*) multiplicity[0][2][nCentralityBins]->Clone("myClone");
  
  double currentContent;
  double binCenter;
  double weight;
  
  for(int iBin = 1; iBin <= testMultiplicity->GetNbinsX(); iBin++){
    currentContent = testMultiplicity->GetBinContent(iBin);
    binCenter = testMultiplicity->GetXaxis()->GetBinCenter(iBin);
    weight = myWeight->Eval(binCenter);
    testMultiplicity->SetBinContent(iBin, currentContent*weight);
  }
  
  // Draw the histogram to the canvas
  drawer->DrawHistogram(testMultiplicity, "Weighted multiplicity", "counts", " ");
  
}
