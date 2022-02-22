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
  TString fileName = "PbPbMC2018_RecoGen_akCaloJet_onlyJets_noCentShift_multiplicityWeight_jetEta1v3_processed_2022-01-24.root";
  // dijetPbPb2018_akCaloJet_onlyJets_jet100TriggerEta1v3_withTrackEff_processed_2022-01-19.root
  // PbPbMC2018_RecoGen_akCaloJet_onlyJets_noCentShift_jetEta1v3_reprocessed_2022-01-19.root
  // PbPbMC2018_RecoReco_akCaloJet_onlyJets_noCentShift_withTrackEff_jetEta1v3_processed_2022-01-20.root
  // dijetPbPb2018_akCaloJet_onlyJets_jet100TriggerEta1v3_aprilTrackEff_processed_2022-01-20.root
  // PbPbMC2018_RecoGen_akCaloJet_onlyJets_noCentShift_noCentWeight_jetEta1v3_processed_2022-01-21.root
  // PbPbMC2018_RecoGen_akCaloJet_onlyJets_noCentShift_multiplicityWeight_jetEta1v3_processed_2022-01-24.root
  
  TString systemString = "Data";  // System to be put on a label
  
  // Open the file
  TFile *multiplicityFile = TFile::Open(directoryName+fileName);
  
  // Create a histogram manager for the histograms in the file
  DijetHistogramManager *manager = new DijetHistogramManager(multiplicityFile);
  manager->SetLoadEventInformation(true);
  manager->LoadProcessedHistograms();
  
  // Binning information
  const int nCentralityBins = 4;
  const int nTrackPtBins = 4;
  const int nAsymmetryBins = 3;
  const int nFlowComponents = 4;
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  int maxCentralityBin = manager->GetNCentralityBins();
  
  const bool drawMultiplicity = false;             // Draw multiplicity in each centrality bin
  const bool drawMultiplicityWeighted = false;     // Draw efficiency weighted multiplicity in each centrality bin
  const bool drawMultiplicityMap = false;          // Draw multiplicity vs. centrality map
  const bool drawMultiplicityMapWeighted = false;   // Draw efficiency weighted multiplicity vs. centrality map
  
  const bool drawMultiplicityDijet = true;             // Draw multiplicity in each centrality bin in dijet events
  const bool drawMultiplicityWeightedDijet = false;     // Draw efficiency weighted multiplicity in each centrality bin in dijet events
  const bool drawMultiplicityMapDijet = false;          // Draw multiplicity vs. centrality map in dijet events
  const bool drawMultiplicityMapWeightedDijet = false;   // Draw efficiency weighted multiplicity vs. centrality map in dijet events
  
  const bool drawAllToSameFigure = true;
  
  const bool saveFigures = false;
  TString saveComment = "_data";
  TString figureFormat = "png";
  
  // Define histograms
  TH1D *multiplicity[4][nCentralityBins+1];  // First bin: multiplicity, weighted, dijet, dijet weighted
  TH2D *multiplicityMap[4];                // First bin: multiplicity, weighted, dijet, dijet weighted
  
  // Read the histograms from the histogram manager
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    multiplicity[0][iCentrality] = manager->GetHistogramMultiplicity(iCentrality);
    multiplicity[1][iCentrality] = manager->GetHistogramMultiplicityWeighted(iCentrality);
    multiplicity[2][iCentrality] = manager->GetHistogramMultiplicityDijet(iCentrality);
    multiplicity[3][iCentrality] = manager->GetHistogramMultiplicityDijetWeighted(iCentrality);
  } // centrality loop
  
  // Load multiplicity histograms without centrality binning
  multiplicity[0][nCentralityBins] = manager->GetHistogramMultiplicity(maxCentralityBin);
  multiplicity[1][nCentralityBins] = manager->GetHistogramMultiplicityWeighted(maxCentralityBin);
  multiplicity[2][nCentralityBins] = manager->GetHistogramMultiplicityDijet(maxCentralityBin);
  multiplicity[3][nCentralityBins] = manager->GetHistogramMultiplicityDijetWeighted(maxCentralityBin);
  
  // Load multiplicity map histograms
  multiplicityMap[0] = manager->GetHistogramMultiplicityMap();
  multiplicityMap[1] = manager->GetHistogramWeightedMultiplicityMap();
  multiplicityMap[2] = manager->GetHistogramMultiplicityMapDijet();
  multiplicityMap[3] = manager->GetHistogramWeightedMultiplicityMapDijet();
  
  TString saveString[] = {"", "Weighted", "Dijet", "WeightedDijet"};
  TString legendString[] = {"", "", "Dijet events", "Dijet events"};
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
  double multiplicityBinBorders[] = {0,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3400,3600,3800,4000,5000};
  int lowBin, highBin;
  
  // Draw the distributions
  for(int iType = 0; iType < 4; iType++){
    
    if(!drawMultiplicityArray[iType]) continue;
    
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      
      // Draw the histogram to the canvas
      drawer->DrawHistogram(multiplicity[iType][iCentrality], "Track multiplicity", "counts", " ");
      
      // Draw a legend for the histogram
      legend = new TLegend(0.17,0.7,0.37,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      
      if(iCentrality < nCentralityBins) legend->AddEntry((TObject*) 0, Form("C = %.0f-%.0f", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "");
      if(legendString[iType] != "") legend->AddEntry((TObject*) 0, legendString[iType], "");
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
        
        for(int iMultiplicity = 0; iMultiplicity < nMultiplicityBins; iMultiplicity++){
          lowBin = multiplicity[iType][iCentrality]->GetXaxis()->FindBin(multiplicityBinBorders[iMultiplicity]+0.1);
          highBin = multiplicity[iType][iCentrality]->GetXaxis()->FindBin(multiplicityBinBorders[iMultiplicity+1]-0.1);
          multiplicity[iType][iCentrality]->GetXaxis()->SetRange(lowBin,highBin);
          cout << "Multiplicity bin " << multiplicityBinBorders[iMultiplicity] << "-" << multiplicityBinBorders[iMultiplicity+1] << "  Mean: " << multiplicity[iType][iCentrality]->GetMean(1) << endl;
        }
        multiplicity[iType][iCentrality]->GetXaxis()->SetRange(0,0);
      }
    } // Centrality loop
  } // Drawing track multiplicity
  
  // Draw all the centrality bins to the same figure
  if(drawAllToSameFigure){
    int centralityColors[] = {kBlue, kRed, kGreen+3, kCyan, kBlack};
    for(int iType = 0; iType < 4; iType++){

      if(!drawMultiplicityArray[iType]) continue;

      for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){

        multiplicity[iType][iCentrality]->SetLineColor(centralityColors[iCentrality]);

        // Draw the histogram to the canvas
        if(iCentrality == 0){
          drawer->DrawHistogram(multiplicity[iType][iCentrality], "Track multiplicity", "counts", " ");
          // Draw a legend for the histogram
          legend = new TLegend(0.17,0.6,0.37,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          if(legendString[iType] != "") legend->AddEntry((TObject*) 0, legendString[iType], "");
        } else {
          multiplicity[iType][iCentrality]->Draw("same");
        }

        if(iCentrality < nCentralityBins){
          legend->AddEntry(multiplicity[iType][iCentrality], Form("C = %.0f-%.0f", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "l");
        } else {
          legend->AddEntry(multiplicity[iType][iCentrality], "Total", "l");
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
    multiplicityMap[iType]->GetYaxis()->SetRangeUser(0,20);
    multiplicityMap[iType]->GetXaxis()->SetRangeUser(0,4000);
    
    // Draw the histogram to the canvas
    drawer->DrawHistogram(multiplicityMap[iType], "Track multiplicity", "Centrality", " ", "colz");
    
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
  TH1D *testMultiplicity = (TH1D*) multiplicity[2][nCentralityBins]->Clone("myClone");
  
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
