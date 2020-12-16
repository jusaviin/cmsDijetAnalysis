#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "stackHistXiao.h"
#include "xCanvas.h"
#include "JDrawer.h"

void plotJetShapeXiao(const int nDatasets, DijetHistogramManager *ppHistograms[5], DijetHistogramManager *pbpbHistograms[5], JffCorrector *ppUncertaintyProvider, JffCorrector *pbpbUncertaintyProvider, int iJetTrack){

  const int nTrackPtBins = pbpbHistograms[0]->GetNTrackPtBins();
  const int nCentralityBins = pbpbHistograms[0]->GetNCentralityBins();
  const int nAsymmetryBins = pbpbHistograms[0]->GetNAsymmetryBins();
  const int nRelevantUncertainties = 5;
  
  const bool saveHistogramsForHepData = true; // Save the plotted histograms to a file for HepData submission
  
  const char* jetShapeTitle[] = {"Leading jet shape ratio between x_{j} bins","Subleading jet shape ratio between x_{j} bins","Jet shape ratio between x_{j} bins"};
  const char* jetType[] = {"Leading", "Subleading", "Inclusive"};
  const char* jetShapeSaveName[] = {"trackLeadingJet","trackSubleadingJet","trackInclusiveJet"};
  const char* xjString[] = {"0.0 < x_{j} < 0.6","0.6 < x_{j} < 0.8","0.8 < x_{j} < 1.0","x_{j} inclusive"};
  const char* asymmetrySaveName[] = {"_A=0v0-0v6","_A=0v6-0v8","_A=0v8-1v0",""};
  
  const char* systemLegendString[] = {"WTA", "E-scheme", "Third jet", "", ""};
  
  const bool addPreliminaryTag = false; // Option to add a preliminary tag to the figures
  
  // Number of events in each xj bin. Provided as a table so these do not have to be in the data file
  double xjJetCount[5][4] = {
  // xj = 0.0     0.6      0.8      1.0      all
           { 229806,  190605,  155115,  575526},  // Centrality 0-10 %
           { 234582,  228419,  203306,  666307},  // Centrality 10-30 %
           {  76433,   91979,   91974,  260386},  // Centrality 30-50 %
           {  20115,   27350,   30536,   78001},  // Centrality 50-90 %
           {2690056, 3913971, 4391792, 10995819}  // pp
  };
    
  // Temporary: Get the ratio between pp and PbPb jet shape
  TH1D *sumHistogramPbPb[5][nCentralityBins][nAsymmetryBins+1];
  TH1D *sumHistogramPp[5][nAsymmetryBins+1];
  TH1D *jetShapeArray[5][nCentralityBins+1][nTrackPtBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPbPb[5][nCentralityBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPp[5][nAsymmetryBins+1];
  TH1D *uncertaintySummerPp[5][nAsymmetryBins+1][nRelevantUncertainties];
  TH1D *uncertaintySummerPbPb[5][nCentralityBins][nAsymmetryBins+1][nRelevantUncertainties];
  TH1D *sumUncertainty[5][nCentralityBins+1][nAsymmetryBins+1];
  TH1D *helperHistogram;
  double shapeIntegralPbPb[5][nCentralityBins][nAsymmetryBins+1];
  double shapeIntegralPp[5][nAsymmetryBins+1];
  double uncertaintySum;
  double correlationFactor;
  
  // Read the pT summed jet shape histograms
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        sumHistogramPbPb[iDataset][iCentrality][iAsymmetry] = (TH1D*) pbpbHistograms[iDataset]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, 0)->Clone(Form("sumPbPb%d",iCentrality));
        jetShapeArray[iDataset][iCentrality][0][iAsymmetry] = pbpbHistograms[iDataset]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, 0);
      }
      sumHistogramPp[iDataset][iAsymmetry] = (TH1D*) ppHistograms[iDataset]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, 0)->Clone("sumPp");
      jetShapeArray[iDataset][nCentralityBins][0][iAsymmetry] = ppHistograms[iDataset]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, 0);
    } // Asymmetry loop
  } // Dataset loop
  
  // Sum the pT:s
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iTrackPt = 1; iTrackPt < nTrackPtBins; iTrackPt++){
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          helperHistogram = (TH1D*) pbpbHistograms[iDataset]->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,iCentrality,iTrackPt)->Clone();
          sumHistogramPbPb[iDataset][iCentrality][iAsymmetry]->Add(helperHistogram);
          jetShapeArray[iDataset][iCentrality][iTrackPt][iAsymmetry] = pbpbHistograms[iDataset]->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,iCentrality,iTrackPt);
        }
        helperHistogram = (TH1D*) ppHistograms[iDataset]->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,0,iTrackPt)->Clone();
        sumHistogramPp[iDataset][iAsymmetry]->Add(helperHistogram);
        jetShapeArray[iDataset][nCentralityBins][iTrackPt][iAsymmetry] = ppHistograms[iDataset]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, iTrackPt);
      } // track pT
    } // Asymmetry loop
  } // Dataset loop
  
  // Normalize the jet shape histograms
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        shapeIntegralPbPb[iDataset][iCentrality][iAsymmetry] = sumHistogramPbPb[iDataset][iCentrality][iAsymmetry]->Integral(1,sumHistogramPbPb[iDataset][iCentrality][iAsymmetry]->FindBin(0.99),"width");
        sumHistogramPbPb[iDataset][iCentrality][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iDataset][iCentrality][iAsymmetry]);
        
        // To draw the systematic uncertainties, clone the sum histograms to uncertainty histograms
        sumUncertainty[iDataset][iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iDataset][iCentrality][iAsymmetry]->Clone(Form("sumUncertaintyPbPb%d%d%d", iDataset, iCentrality, iAsymmetry));
        
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          jetShapeArray[iDataset][iCentrality][iTrackPt][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iDataset][iCentrality][iAsymmetry]);
        }
        
      }
      
      shapeIntegralPp[iDataset][iAsymmetry] = sumHistogramPp[iDataset][iAsymmetry]->Integral(1,sumHistogramPp[iDataset][iAsymmetry]->FindBin(0.99),"width");
      sumHistogramPp[iDataset][iAsymmetry]->Scale(1.0/shapeIntegralPp[iDataset][iAsymmetry]);
      
      // To draw the systematic uncertainties, clone the sum histograms to uncertainty histograms
      sumUncertainty[iDataset][nCentralityBins][iAsymmetry] = (TH1D*) sumHistogramPp[iDataset][iAsymmetry]->Clone(Form("sumUncertaintyPp%d%d", iDataset, iAsymmetry));
      
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        jetShapeArray[iDataset][nCentralityBins][iTrackPt][iAsymmetry]->Scale(1.0/shapeIntegralPp[iDataset][iAsymmetry]);
      }
      
      // Read the uncertainties for the pT summed jet shapes but skip a few sources not relevant here
      uncertaintySummerPp[iDataset][iAsymmetry][0] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry, JffCorrector::kBackgroundFluctuation);
      uncertaintySummerPp[iDataset][iAsymmetry][1] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry, JffCorrector::kFragmentationBias);
      uncertaintySummerPp[iDataset][iAsymmetry][2] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry, JffCorrector::kPairAcceptance);
      uncertaintySummerPp[iDataset][iAsymmetry][3] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry, JffCorrector::kBackgroundSubtraction);
      uncertaintySummerPp[iDataset][iAsymmetry][4] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry, JffCorrector::kJetResolution);
      
      uncertaintyHistogramPp[iDataset][iAsymmetry] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry);
      
      
      // Add the relevant uncertainties in quadrature for pp
      for(int iBin = 1; iBin <= uncertaintyHistogramPp[iDataset][iAsymmetry]->GetNbinsX(); iBin++){
        uncertaintySum = 0;
        for(int iUncertainty = 0; iUncertainty < nRelevantUncertainties; iUncertainty++){
          
          // Old files do not have jet resolution uncertainty. Do not crash the code for that
          if(uncertaintySummerPp[iDataset][iAsymmetry][iUncertainty] == NULL) continue;
          
          uncertaintySum += TMath::Power(uncertaintySummerPp[iDataset][iAsymmetry][iUncertainty]->GetBinContent(iBin),2);
        }
        uncertaintyHistogramPp[iDataset][iAsymmetry]->SetBinContent(iBin,TMath::Sqrt(uncertaintySum));
      }
      
      // Scale the histograms from jet shapes
      uncertaintyHistogramPp[iDataset][iAsymmetry]->Scale(1.0/shapeIntegralPp[iDataset][iAsymmetry]);
      
      // Set the bin errors in the sumUncertainty histograms to match the uncertainties read from the file
      for(int iBin = 1; iBin <= sumUncertainty[iDataset][nAsymmetryBins][iAsymmetry]->GetNbinsX(); iBin++){
        sumUncertainty[iDataset][nCentralityBins][iAsymmetry]->SetBinError(iBin, uncertaintyHistogramPp[iDataset][iAsymmetry]->GetBinContent(iBin));
      }
      
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        
        // Read the uncertainties for the pT summed jet shapes but skip a few sources not relevant here
        uncertaintySummerPbPb[iDataset][iCentrality][iAsymmetry][0] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, nTrackPtBins, iAsymmetry, JffCorrector::kBackgroundFluctuation);
        uncertaintySummerPbPb[iDataset][iCentrality][iAsymmetry][1] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, nTrackPtBins, iAsymmetry, JffCorrector::kFragmentationBias);
        uncertaintySummerPbPb[iDataset][iCentrality][iAsymmetry][2] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, nTrackPtBins, iAsymmetry, JffCorrector::kPairAcceptance);
        uncertaintySummerPbPb[iDataset][iCentrality][iAsymmetry][3] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, nTrackPtBins, iAsymmetry, JffCorrector::kBackgroundSubtraction);
        uncertaintySummerPbPb[iDataset][iCentrality][iAsymmetry][4] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, nTrackPtBins, iAsymmetry, JffCorrector::kJetResolution);
        
        uncertaintyHistogramPbPb[iDataset][iCentrality][iAsymmetry] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, nTrackPtBins, iAsymmetry);
        
        // Add the relevant uncertainties in quadrature for PbPb
        for(int iBin = 1; iBin <= uncertaintyHistogramPbPb[iDataset][iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
          uncertaintySum = 0;
          for(int iUncertainty = 0; iUncertainty < nRelevantUncertainties; iUncertainty++){
            
            // Old files do not have jet resolution uncertainty. Do not crash the code for that
            if(uncertaintySummerPbPb[iDataset][iCentrality][iAsymmetry][iUncertainty] == NULL) continue;
            
            uncertaintySum += TMath::Power(uncertaintySummerPbPb[iDataset][iCentrality][iAsymmetry][iUncertainty]->GetBinContent(iBin),2);
          }
          uncertaintyHistogramPbPb[iDataset][iCentrality][iAsymmetry]->SetBinContent(iBin,TMath::Sqrt(uncertaintySum));
        }
        
        uncertaintyHistogramPbPb[iDataset][iCentrality][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iDataset][iCentrality][iAsymmetry]);
        
        // Set the bin errors in the sumUncertainty histograms to match the uncertainties read from the file
        for(int iBin = 1; iBin <= sumUncertainty[iDataset][iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
          sumUncertainty[iDataset][iCentrality][iAsymmetry]->SetBinError(iBin, uncertaintyHistogramPbPb[iDataset][iCentrality][iAsymmetry]->GetBinContent(iBin));
        }
        
      }
    } // Asymmetry loop
  } // Dataset loop

  // Create the asymmetric to nominal and symmetric to nominal ratios
  TH1D* asymmetryRatioHistogram[5][nCentralityBins+1][nAsymmetryBins];
  TH1D* asymmetryRatioUncertainty[5][nCentralityBins+1][nAsymmetryBins];
  double ppValue, pbpbValue, ratioValue;
  double ppUncertainty, pbpbUncertainty, ratioUncertaintyValue;
  double correlationSpillover, correlationJff, correlationResolution;
  double uncertaintyXjBin, uncertaintyXjIntegrated;
  
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
        
        if(iCentrality == nCentralityBins){
          asymmetryRatioHistogram[iDataset][iCentrality][iAsymmetry] = (TH1D*) sumHistogramPp[iDataset][iAsymmetry]->Clone(Form("ratioAsymmetry%d%d%d", iDataset, iCentrality, iAsymmetry));
          asymmetryRatioHistogram[iDataset][iCentrality][iAsymmetry]->Divide(sumHistogramPp[iDataset][nAsymmetryBins]);
        } else {
          asymmetryRatioHistogram[iDataset][iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iDataset][iCentrality][iAsymmetry]->Clone(Form("ratioAsymmetry%d%d%d", iDataset, iCentrality, iAsymmetry));
          asymmetryRatioHistogram[iDataset][iCentrality][iAsymmetry]->Divide(sumHistogramPbPb[iDataset][iCentrality][nAsymmetryBins]);
        }
        
        // The statistical uncertainties are partially correlated, since xj bins are part of the xj integrated sample
        // To take this into account, we need to calculate the statistical uncertainty by hand
        for(int iBin = 1; iBin < asymmetryRatioHistogram[iDataset][iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
          if(iCentrality == nCentralityBins){
            ppValue = sumHistogramPp[iDataset][nAsymmetryBins]->GetBinContent(iBin);
            pbpbValue = sumHistogramPp[iDataset][iAsymmetry]->GetBinContent(iBin);
            ppUncertainty = sumHistogramPp[iDataset][nAsymmetryBins]->GetBinError(iBin);
            pbpbUncertainty = sumHistogramPp[iDataset][iAsymmetry]->GetBinError(iBin);
          } else {
            ppValue = sumHistogramPbPb[iDataset][iCentrality][nAsymmetryBins]->GetBinContent(iBin);
            pbpbValue = sumHistogramPbPb[iDataset][iCentrality][iAsymmetry]->GetBinContent(iBin);
            ppUncertainty = sumHistogramPbPb[iDataset][iCentrality][nAsymmetryBins]->GetBinError(iBin);
            pbpbUncertainty = sumHistogramPbPb[iDataset][iCentrality][iAsymmetry]->GetBinError(iBin);
          }
          
          // Correlation term
          ratioValue = 2 * (1.0/ppValue) * (-1.0 * pbpbValue/TMath::Power(ppValue,2)) * (xjJetCount[iCentrality][iAsymmetry] / xjJetCount[iCentrality][nAsymmetryBins]) * ppUncertainty * pbpbUncertainty;
          
          ratioUncertaintyValue = TMath::Sqrt(TMath::Power(pbpbUncertainty/ppValue,2)+TMath::Power(pbpbValue*ppUncertainty/TMath::Power(ppValue,2),2) + ratioValue);
                              
          asymmetryRatioHistogram[iDataset][iCentrality][iAsymmetry]->SetBinError(iBin,ratioUncertaintyValue);
          
        } // Calculating statistical uncertainty
        
        asymmetryRatioUncertainty[iDataset][iCentrality][iAsymmetry] = (TH1D*) asymmetryRatioHistogram[iDataset][iCentrality][iAsymmetry]->Clone(Form("asymmetryRatioUncertainty%d%d%d", iDataset, iAsymmetry, iCentrality));
        
        // Calculate the systematic uncertainty for the ratio
        for(int iBin = 1; iBin < asymmetryRatioUncertainty[iDataset][iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
          if(iCentrality == nCentralityBins){
            ppValue = sumHistogramPp[iDataset][nAsymmetryBins]->GetBinContent(iBin);
            pbpbValue = sumHistogramPp[iDataset][iAsymmetry]->GetBinContent(iBin);
            ppUncertainty = uncertaintyHistogramPp[iDataset][nAsymmetryBins]->GetBinContent(iBin);
            pbpbUncertainty = uncertaintyHistogramPp[iDataset][iAsymmetry]->GetBinContent(iBin);
            
            // No correlation terms for pp
            correlationSpillover = 0;
            correlationJff = 0;
            correlationResolution = 0;
            
          } else {
            ppValue = sumHistogramPbPb[iDataset][iCentrality][nAsymmetryBins]->GetBinContent(iBin);
            pbpbValue = sumHistogramPbPb[iDataset][iCentrality][iAsymmetry]->GetBinContent(iBin);
            ppUncertainty = uncertaintyHistogramPbPb[iDataset][iCentrality][nAsymmetryBins]->GetBinContent(iBin);
            pbpbUncertainty = uncertaintyHistogramPbPb[iDataset][iCentrality][iAsymmetry]->GetBinContent(iBin);
            
            // JFF and spillover uncertainties share partially the datasets from which they are derived, so we need correlation terms for them
            uncertaintyXjBin = uncertaintySummerPbPb[iDataset][iCentrality][iAsymmetry][0]->GetBinContent(iBin) / shapeIntegralPbPb[iDataset][iCentrality][iAsymmetry];
            uncertaintyXjIntegrated = uncertaintySummerPbPb[iDataset][iCentrality][nAsymmetryBins][0]->GetBinContent(iBin) / shapeIntegralPbPb[iDataset][iCentrality][nAsymmetryBins];
            correlationSpillover = 2 * (1.0/ppValue) * (-1.0 * pbpbValue/TMath::Power(ppValue,2)) * (xjJetCount[iCentrality][iAsymmetry] / xjJetCount[iCentrality][nAsymmetryBins]) * uncertaintyXjBin * uncertaintyXjIntegrated;
            
            uncertaintyXjBin = uncertaintySummerPbPb[iDataset][iCentrality][iAsymmetry][1]->GetBinContent(iBin) / shapeIntegralPbPb[iDataset][iCentrality][iAsymmetry];
            uncertaintyXjIntegrated = uncertaintySummerPbPb[iDataset][iCentrality][nAsymmetryBins][1]->GetBinContent(iBin) / shapeIntegralPbPb[iDataset][iCentrality][nAsymmetryBins];
            correlationJff = 2 * (1.0/ppValue) * (-1.0 * pbpbValue/TMath::Power(ppValue,2)) * (xjJetCount[iCentrality][iAsymmetry] / xjJetCount[iCentrality][nAsymmetryBins]) * uncertaintyXjBin * uncertaintyXjIntegrated;
            
            uncertaintyXjBin = uncertaintySummerPbPb[iDataset][iCentrality][iAsymmetry][4]->GetBinContent(iBin) / shapeIntegralPbPb[iDataset][iCentrality][iAsymmetry];
            uncertaintyXjIntegrated = uncertaintySummerPbPb[iDataset][iCentrality][nAsymmetryBins][4]->GetBinContent(iBin) / shapeIntegralPbPb[iDataset][iCentrality][nAsymmetryBins];
            correlationResolution = 2 * (1.0/ppValue) * (-1.0 * pbpbValue/TMath::Power(ppValue,2)) * (xjJetCount[iCentrality][iAsymmetry] / xjJetCount[iCentrality][nAsymmetryBins]) * uncertaintyXjBin * uncertaintyXjIntegrated;
            
          }
          
          ratioValue = asymmetryRatioHistogram[iDataset][iCentrality][iAsymmetry]->GetBinContent(iBin);
          ratioUncertaintyValue = TMath::Sqrt(TMath::Power(pbpbUncertainty/ppValue,2) + TMath::Power(pbpbValue*ppUncertainty/TMath::Power(ppValue,2),2) + correlationSpillover + correlationJff + correlationResolution);
          asymmetryRatioUncertainty[iDataset][iCentrality][iAsymmetry]->SetBinContent(iBin,ratioValue);
          asymmetryRatioUncertainty[iDataset][iCentrality][iAsymmetry]->SetBinError(iBin,ratioUncertaintyValue);
        }
        
      } // Centrality loop
    } // Asymmetry loop
  } // Dataset loop
  
  TString cent_lab[4] = {"0-10%", "10-30%", "30-50%", "50-90%"};
  int markerColors[5] = {kBlack, kBlue, kRed, kGreen+3, kMagenta};
  int fillColors[5] = {kGray+3, kBlue-6, kRed-6, kGreen-6, kMagenta-6};
  
  // Set a good drawing style for the uncertainty histograms
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    for(int iCentrality = 0; iCentrality < nCentralityBins+1; iCentrality++){
      for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
        asymmetryRatioUncertainty[iDataset][iCentrality][iAsymmetry]->SetFillStyle(1001);
        asymmetryRatioUncertainty[iDataset][iCentrality][iAsymmetry]->SetFillColorAlpha(fillColors[iDataset],0.4);
        asymmetryRatioUncertainty[iDataset][iCentrality][iAsymmetry]->SetMarkerStyle(20);
        asymmetryRatioUncertainty[iDataset][iCentrality][iAsymmetry]->SetMarkerSize(1.6);
        asymmetryRatioUncertainty[iDataset][iCentrality][iAsymmetry]->SetMarkerColor(markerColors[iDataset]);
        asymmetryRatioUncertainty[iDataset][iCentrality][iAsymmetry]->SetLineColor(markerColors[iDataset]);
      } // Asymmetry loop
    } // Centrality loop
  } // Dataset loop
  
  // Draw all the distributions to big canvases
  auxi_canvas *bigCanvas;
  TBox *box;
  TLatex *mainTitle;
  TLine *line;
  
  // Draw a big canvas and put all the plots in it
  bigCanvas = new auxi_canvas(Form("bigCanvas%d",iJetTrack), "", 2500, 1300);
  bigCanvas->SetMargin(0.07, 0.01, 0.10, 0.15);
  bigCanvas->divide(2,5);
  
  mainTitle = new TLatex();
  double axisZoom[] = {1.6,2.4,2.2};
  double axisLow[] = {0.4,0.2,0.5};
  
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    bigCanvas->CD(5-iCentrality);
    
    asymmetryRatioHistogram[0][iCentrality][0]->SetTitle("");
    asymmetryRatioHistogram[0][iCentrality][0]->SetAxisRange(axisLow[iJetTrack/3], axisZoom[iJetTrack/3], "Y");
    asymmetryRatioHistogram[0][iCentrality][0]->SetAxisRange(0, 0.99, "X");
    if(iCentrality < 4)  {
      asymmetryRatioHistogram[0][iCentrality][0]->GetXaxis()->SetTitleOffset(0.7);
      asymmetryRatioHistogram[0][iCentrality][0]->GetXaxis()->SetTitleSize(0.11);
      asymmetryRatioHistogram[0][iCentrality][0]->GetXaxis()->SetNdivisions(505);
      asymmetryRatioHistogram[0][iCentrality][0]->GetXaxis()->SetLabelSize(0.09);
      asymmetryRatioHistogram[0][iCentrality][0]->GetYaxis()->SetNdivisions(505);
    }
    if(iCentrality == 4){
      asymmetryRatioHistogram[0][iCentrality][0]->GetXaxis()->SetTitleOffset(0.92);
      asymmetryRatioHistogram[0][iCentrality][0]->GetXaxis()->SetTitleSize(0.086);
      asymmetryRatioHistogram[0][iCentrality][0]->GetXaxis()->SetNdivisions(505);
      asymmetryRatioHistogram[0][iCentrality][0]->GetXaxis()->SetLabelOffset(0.016);
      asymmetryRatioHistogram[0][iCentrality][0]->GetXaxis()->SetLabelSize(0.064);
      asymmetryRatioHistogram[0][iCentrality][0]->GetYaxis()->SetNdivisions(505);
      asymmetryRatioHistogram[0][iCentrality][0]->GetYaxis()->SetLabelOffset(0.02);
      asymmetryRatioHistogram[0][iCentrality][0]->GetYaxis()->SetLabelSize(0.09);
      asymmetryRatioHistogram[0][iCentrality][0]->GetYaxis()->SetTitleOffset(0.9);
      asymmetryRatioHistogram[0][iCentrality][0]->GetYaxis()->SetTitleSize(0.12);
      asymmetryRatioHistogram[0][iCentrality][0]->GetYaxis()->SetTitle("#rho(#Deltar)_{x_{j} < 0.6 }/#rho(#Deltar)_{All}");
      mainTitle->SetTextSize(0.073);
    }
    
    asymmetryRatioHistogram[0][iCentrality][0]->GetYaxis()->CenterTitle();
    asymmetryRatioHistogram[0][iCentrality][0]->GetXaxis()->CenterTitle();
    asymmetryRatioHistogram[0][iCentrality][0]->SetStats(0);
    asymmetryRatioHistogram[0][iCentrality][0]->SetLineColor(markerColors[0]);
    asymmetryRatioHistogram[0][iCentrality][0]->Draw();
    
    line = new TLine();
    line->SetLineStyle(2);
    line->DrawLine(0, 1, 1, 1);
    
    asymmetryRatioUncertainty[0][iCentrality][0]->Draw("same e2");
    
    for(int iDataset = 1; iDataset < nDatasets; iDataset++){
      asymmetryRatioHistogram[iDataset][iCentrality][0]->SetLineColor(markerColors[iDataset]);
      asymmetryRatioHistogram[iDataset][iCentrality][0]->Draw("same");
      asymmetryRatioUncertainty[iDataset][iCentrality][0]->Draw("same e2");
    }
    
    if(iCentrality == 4){
      mainTitle->SetTextFont(22);
      mainTitle->SetTextSize(0.095);
      mainTitle->DrawLatexNDC(0.35, 0.9, "pp reference");
      //mainTitle->DrawLatexNDC(.35, 0.9, "pp reference");
    }
    else {
      mainTitle->SetTextFont(22);
      mainTitle->SetTextSize(0.1);
      mainTitle->DrawLatexNDC(0.17, 0.9, "PbPb");
      mainTitle->DrawLatexNDC(0.17, 0.82, cent_lab[iCentrality]);
    }
    
    bigCanvas->CD(10-iCentrality);
    
    asymmetryRatioHistogram[0][iCentrality][2]->SetTitle("");
    asymmetryRatioHistogram[0][iCentrality][2]->GetXaxis()->SetTitleOffset(1.1);
    asymmetryRatioHistogram[0][iCentrality][2]->GetXaxis()->SetTitle("#Deltar");
    asymmetryRatioHistogram[0][iCentrality][2]->SetAxisRange(axisLow[iJetTrack/3], axisZoom[iJetTrack/3], "Y");
    asymmetryRatioHistogram[0][iCentrality][2]->SetAxisRange(0, 0.99, "X");
    if( iCentrality<4 )  {
      asymmetryRatioHistogram[0][iCentrality][2]->GetXaxis()->SetTitleOffset(0.65);
      asymmetryRatioHistogram[0][iCentrality][2]->GetXaxis()->SetTitleSize(0.13);
      asymmetryRatioHistogram[0][iCentrality][2]->GetXaxis()->SetNdivisions(505);
      asymmetryRatioHistogram[0][iCentrality][2]->GetXaxis()->SetLabelOffset(0.00001);
      asymmetryRatioHistogram[0][iCentrality][2]->GetXaxis()->SetLabelSize(0.093);
      asymmetryRatioHistogram[0][iCentrality][2]->GetYaxis()->SetNdivisions(505);
    }
    if(iCentrality==4 ){
      asymmetryRatioHistogram[0][iCentrality][2]->GetXaxis()->SetTitleOffset(0.86);
      asymmetryRatioHistogram[0][iCentrality][2]->GetXaxis()->SetTitleSize(0.098);
      asymmetryRatioHistogram[0][iCentrality][2]->GetXaxis()->SetNdivisions(505);
      asymmetryRatioHistogram[0][iCentrality][2]->GetXaxis()->SetLabelOffset(0.016);
      asymmetryRatioHistogram[0][iCentrality][2]->GetXaxis()->SetLabelSize(0.071);
      asymmetryRatioHistogram[0][iCentrality][2]->GetYaxis()->SetNdivisions(505);
      asymmetryRatioHistogram[0][iCentrality][2]->GetYaxis()->SetLabelOffset(0.02);
      asymmetryRatioHistogram[0][iCentrality][2]->GetYaxis()->SetLabelSize(0.07);
      asymmetryRatioHistogram[0][iCentrality][2]->GetYaxis()->SetTitleOffset(1.1);
      asymmetryRatioHistogram[0][iCentrality][2]->GetYaxis()->SetTitleSize(0.1);
      asymmetryRatioHistogram[0][iCentrality][2]->GetYaxis()->SetTitle("#rho(#Deltar)_{x_{j} > 0.8 }/#rho(#Deltar)_{All}");
      mainTitle->SetTextSize(0.073);
    }

    
    asymmetryRatioHistogram[0][iCentrality][2]->GetYaxis()->CenterTitle();
    asymmetryRatioHistogram[0][iCentrality][2]->GetXaxis()->CenterTitle();
    asymmetryRatioHistogram[0][iCentrality][2]->SetStats(0);
    asymmetryRatioHistogram[0][iCentrality][2]->SetLineColor(markerColors[0]);
    asymmetryRatioHistogram[0][iCentrality][2]->Draw("same");
    
    line = new TLine();
    line->SetLineStyle(2);
    line->DrawLine(0, 1, 1, 1);
    
    asymmetryRatioUncertainty[0][iCentrality][2]->Draw("same e2");
    
    for(int iDataset = 1; iDataset < nDatasets; iDataset++){
      asymmetryRatioHistogram[iDataset][iCentrality][2]->SetLineColor(markerColors[iDataset]);
      asymmetryRatioHistogram[iDataset][iCentrality][2]->Draw("same");
      asymmetryRatioUncertainty[iDataset][iCentrality][2]->Draw("same e2");
    }
    
    if(nDatasets > 1 && iCentrality == 4){
      TLegend *systemLegend = new TLegend(0.31,0.65,0.81,0.95);
      systemLegend->SetTextSize(0.06);
      systemLegend->SetLineColor(kWhite);
      systemLegend->SetFillColor(kWhite);
      for(int iDataset = 0; iDataset < nDatasets; iDataset++){
        systemLegend->AddEntry(asymmetryRatioUncertainty[iDataset][iCentrality][0], systemLegendString[iDataset] ,"lp");
      }
      systemLegend->Draw();
    }
    
  } // Centrality loop
  
  TLegend* lt1 = new TLegend(0.01, 0.1, 1, 0.5);
  TLegend* lt2 = new TLegend(0.0, 0.1, 1, 0.5);
  TLegend* lt3 = new TLegend(0.0 ,0.22, 1, 0.5);
  TLegend* lt5 = new TLegend(0.01, 0.5, 1, 0.68);
  TLegend* lt4 = new TLegend(0.25, 0.85, 0.85, 0.95);
  lt1->SetTextSize(0.07);
  lt1->SetLineColor(kWhite);
  lt1->SetFillColor(kWhite);
  lt2->SetTextSize(0.07);
  lt2->SetLineColor(kWhite);
  lt2->SetFillColor(kWhite);
  lt3->SetTextSize(0.07);
  lt3->SetLineColor(kWhite);
  lt3->SetFillColor(kWhite);
  lt4->SetTextSize(0.06);
  lt4->SetLineColor(kWhite);
  lt4->SetFillColor(kWhite);
  lt5->SetTextSize(0.07);
  lt5->SetLineColor(kWhite);
  lt5->SetFillColor(kWhite);
  
  
  /*lt1->AddEntry(jetShapeStack[4][iAsymmetry]->hist_trunk.at(0), "0.7 < p_{T}^{trk}< 1 GeV","f");
  lt1->AddEntry(jetShapeStack[4][iAsymmetry]->hist_trunk.at(1), "1 < p_{T}^{trk}< 2 GeV","f");
  lt1->AddEntry(jetShapeStack[4][iAsymmetry]->hist_trunk.at(2), "2 < p_{T}^{trk}< 3 GeV","f");
  
  lt2->AddEntry(jetShapeStack[4][iAsymmetry]->hist_trunk.at(3), "3 < p_{T}^{trk}< 4 GeV","f");
  lt2->AddEntry(jetShapeStack[4][iAsymmetry]->hist_trunk.at(4), "4 < p_{T}^{trk}< 8 GeV","f");
  lt2->AddEntry(jetShapeStack[4][iAsymmetry]->hist_trunk.at(5), "8 < p_{T}^{trk}< 12 GeV","f");
  
  
  lt3->AddEntry(jetShapeStack[4][iAsymmetry]->hist_trunk.at(6), "12 < p_{T}^{trk}< 300 GeV","f");
  lt3->AddEntry(sumUncertainty[4][iAsymmetry], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
  //lt3->AddEntry(jetShapeStack[4]->hist_trunk.at(8), "20 < p_{T}^{trk}< 300 GeV","f");
  
  lt4->AddEntry(ratioUncertainty[0][iAsymmetry], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
  bigCanvas->CD(9);
  lt4->Draw();*/
  
  //bigCanvas->CD(1);
  //TLegend* lt6 = new TLegend(0.3,0.7,0.98,0.85);
  //lt6->SetTextSize(0.07);
  //lt6->SetLineColor(kWhite);
  //lt6->SetFillColor(kWhite);
  //lt6->AddEntry(js_dr_err_all[0], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
  //lt6->Draw();
  
  //bigCanvas->CD(2);
  //  lt5->AddEntry(js_dr_err_all[0], "0.7 < p_{T}^{track}< 20 GeV","lpfe");
  //  lt5->Draw();
  //line->SetLineStyle(1);
  //line->DrawLineNDC(0, 0, 0, 1);
  //lt1->Draw();
  //bigCanvas->CD(3);
  //lt2->Draw();
  //bigCanvas->CD(4);
  //lt3->Draw();
  
  double cmsPositionX = 0.02;
  double cmsPositionY = 0.9;
  double jetShapeTitlePosition[] = {0.13,0.11,0.13};
  double xjPosition[] = {0.76,0.77,0.76};
  
  if(addPreliminaryTag){
    cmsPositionX = 0.05;
    cmsPositionY = 0.93;
    jetShapeTitlePosition[0] = 0.225;
    jetShapeTitlePosition[1] = 0.21;
    jetShapeTitlePosition[2] = 0.225;
  }
  
  bigCanvas->cd(0);
  
  
  mainTitle->SetTextFont(62);
  mainTitle->SetTextSize(0.065);
  mainTitle->DrawLatexNDC(cmsPositionX, cmsPositionY, "CMS");
  
  if(addPreliminaryTag){
    mainTitle->SetTextFont(42);
    mainTitle->SetTextSize(0.055);
    mainTitle->DrawLatexNDC(0.022, 0.88, "Preliminary");
    
    mainTitle->SetTextFont(42);
    mainTitle->SetTextSize(0.045);
    mainTitle->DrawLatexNDC(jetShapeTitlePosition[iJetTrack/3], 0.935, Form("%s jet shape ratio", jetType[iJetTrack/3]));
    
    mainTitle->SetTextFont(42);
    mainTitle->SetTextSize(0.045);
    mainTitle->DrawLatexNDC(0.26, 0.885, "between x_{j} bins");
    
  } else {
    
    mainTitle->SetTextFont(42);
    mainTitle->SetTextSize(0.045);
    mainTitle->DrawLatexNDC(jetShapeTitlePosition[iJetTrack/3], 0.905, jetShapeTitle[iJetTrack/3]);
    
  }
  
  //mainTitle->SetTextFont(42);
  //mainTitle->SetTextSize(0.035);
  //mainTitle->DrawLatexNDC(xjPosition[iJetTrack/3], 0.94, xjString[iAsymmetry]);
  
  mainTitle->SetTextSize(0.037);
  mainTitle->DrawLatexNDC(0.565, 0.94, "pp 320 pb^{-1} (5.02 TeV)  PbPb 1.7 nb^{-1} (5.02 TeV)");
  mainTitle->SetTextSize(0.034);
  mainTitle->DrawLatexNDC(0.52, 0.885, "anti-k_{T} R = 0.4, |#eta_{jet}| < 1.6, p_{T,1} > 120 GeV, p_{T,2} > 50 GeV, #Delta#varphi_{1,2} > #frac{5#pi}{6}");
  //  lb->drawText("(p_{T}> 120 GeV, |#eta_{jet}|<1.6)", 0.2, 0.25, 4);
  
  box = new TBox();
  box->SetFillColor(kWhite);
  bigCanvas->cd(0);
  box->DrawBox(0.24,0.047, 0.275, 0.09);
  mainTitle->SetTextSize(0.034);
  mainTitle->DrawLatex(0.252, 0.065, "0");
  box->DrawBox(0.42,0.047, 0.45, 0.09);
  box->DrawBox(0.61,0.047, 0.64, 0.09);
  box->DrawBox(0.79,0.047, 0.82, 0.09);
  box->DrawBox(0.97,0.047, 0.99, 0.09);
  //box->DrawBox(0.04,0.46, 0.065, 0.5);
  
  mainTitle->DrawLatex(0.436, 0.065, "0");
  mainTitle->DrawLatex(0.62, 0.065, "0");
  mainTitle->DrawLatex(0.804, 0.065, "0");
  mainTitle->DrawLatex(0.983, 0.065, "1");
  
  //bigCanvas->SaveAs("js_dr_normal_new.eps");
  //bigCanvas->SaveAs(Form("figures/finalJetShapeAsymmetryThirdJetComparison_%s.pdf",jetShapeSaveName[iJetTrack/3]));
  bigCanvas->SaveAs(Form("figures/finalJetShapeAsymmetry_%s_finalStyleUpdates2.pdf",jetShapeSaveName[iJetTrack/3]));
  //bigCanvas->SaveAs("js_dr_normal_v3.eps");
  //bigCanvas->SaveAs("js_dr_normal_v3.pdf");
  
  // Save the histograms to a file for HepData submission
  if(saveHistogramsForHepData){
    TString outputFileName = "hepdata/hepdata_asymmetryRatio_update_hin-19-013.root";
    TFile *outputFile = TFile::Open(outputFileName,"UPDATE");
    TString centralityString[] = {"0-10", "10-30", "30-50", "50-90", "pp"};
    
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      asymmetryRatioHistogram[0][iCentrality][0]->Write(Form("asymmetryRatioImbalanced_%s_%s", jetShapeSaveName[iJetTrack/3], centralityString[iCentrality].Data()), TObject::kOverwrite);
      asymmetryRatioUncertainty[0][iCentrality][0]->GetXaxis()->SetRangeUser(0,1);
      asymmetryRatioUncertainty[0][iCentrality][0]->Write(Form("asymmetryRatioErrorImbalanced_%s_%s", jetShapeSaveName[iJetTrack/3], centralityString[iCentrality].Data()), TObject::kOverwrite);
      asymmetryRatioHistogram[0][iCentrality][2]->Write(Form("asymmetryRatioBalanced_%s_%s", jetShapeSaveName[iJetTrack/3], centralityString[iCentrality].Data()), TObject::kOverwrite);
      asymmetryRatioUncertainty[0][iCentrality][2]->GetXaxis()->SetRangeUser(0,1);
      asymmetryRatioUncertainty[0][iCentrality][2]->Write(Form("asymmetryRatioErrorBalanced_%s_%s", jetShapeSaveName[iJetTrack/3], centralityString[iCentrality].Data()), TObject::kOverwrite);
    }
    
    outputFile->Close();
  }
  
}

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void finalJetShapeRatioPlotter(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Open data files for pp and PbPb data
  TFile *ppFile[5] = { TFile::Open("data/ppData2017_highForest_pfJets_fixedJEC_20EveMixed_wtaAxis_xjBins_allCorrections_processed_2020-11-04.root"), NULL, NULL, NULL, NULL};
  // data/ppData2017_highForest_pfJets_fixedJEC_20EveMixed_wtaAxis_xjBins_allCorrections_processed_2020-11-04.root <-- Fixed pp JEC
  // data/ppData2017_highForest_pfJets_20EveMixed_xjBins_wtaAxis_allCorrections_processed_2020-02-04.root <-- Bug in pp JEC
  // data/ppData2017_highForest_pfJets_20EveMixed_xjBins_finalTrackCorr_JECv4_eschemeAxis_seagullAndJff_processed_2019-10-02.root
  // data/ppData2017_highForest_pfJets_20EventsMixed_finalTrackCorr_xjBins_JECv4_wtaAxis_tunedSeagull_allCorrections_processed_2019-10-17.root
  // data/ppData2017_highForest_pfJets_20EventsMixed_xjBins_finalTrackCorr_JECv4_wtaAxis_allCorrections_processed_2019-09-28.root
  // data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_processed_2019-08-13.root
  // data/dijet_pp_highForest_pfJets_noUncOrInc_allCorrections_wtaAxis_processed_2019-07-13.root
  TFile *pbpbFile[5] = { TFile::Open("data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_allCorrections_finalTuning_onlyJetShapes_processed_2020-02-17.root"), NULL, NULL, NULL, NULL};
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_allCorrections_finalTuning_onlyJetShapes_processed_2020-02-17.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_finalTuning_wtaAxis_processed_2020-02-04.root
  //  data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_finalTuning_wtaAxis_processed_2020-01-29.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrections_seagullTuningProcess_processed_2020-01-15.root
  // dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_JECv6_finalTrack_onlySeagullAndSpillover_correctedCentralityCorrection_eschemeAxis_onlyJetShape_processed_2019-12-05.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_onlySeagullAndSpillover_onlyJetShapes_processed_2019-11-21.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_JECv6_finalTrack_eschemeAxis_onlySeagull_processed_2019-11-21.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_xjBins_JECv6_finalTrack_eschemeAxis_noTrackDeltaRCorrection_firstTry_processed_2019-11-14.root
  // data/PbPbMC2018_RecoGen_akFlowPuCs4PFJet_noUncOrInc_xjBins_5pShiftedCent_5eveMix_jet100Trigger_allCorrections_processed_2019-10-21.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsWithCentShift_trackDeltaRonlyLowPt_processed_2019-10-16.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_allCorrections_noStatErrorFromCorrections_lowPtResidualTrack_processed_2019-10-07_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_allCorrections_processed_2019-10-01_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_calo80Trigger_xjBins_wtaAxis_allCorrections_JECv5b_processed_2019-09-05.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_modifiedSeagull_noErrorJff_averagePeakMixing_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_testNoJffCorr_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_allCorrections_modifiedSeagull_wtaAxis_JECv4_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root
  
  int nFilesPerDataset = 1;
  for(int iFile = 1; iFile < 5; iFile++){
    if(ppFile[iFile] == NULL) break;
    nFilesPerDataset++;
  }
  
  TFile *ppUncertaintyFile = TFile::Open("uncertainties/systematicUncertaintyForPp_20eveMix_xjBins_jffAndJetResolutionUpdate_2020-06-03.root");
  // uncertainties/systematicUncertaintyForPp_20eveMix_xjBins_jffAndJetResolutionUpdate_2020-06-03.root  <-- Final final file
  // uncertainties/systematicUncertaintyForPp_20eveMix_xjBins_fixJES_2020-02-03.root
  // uncertainties/systematicUncertaintyForPp_20percentSpillJff_2019-09-30.root
  // uncertainties/systematicUncertaintyForPythia8RecoGen_mcMode_2019-10-05.root
  TFile *pbpbUncertaintyFile = TFile::Open("uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_manualFluctuationReduction_addTriggerBias_jffUpdate_2020-06-03.root");
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_manualFluctuationReduction_addTriggerBias_jffUpdate_2020-06-03.root <-- Final final file
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_addNewSources_smoothedPairBackground_2020-05-18.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_includeTrackDeltaR_2020-01-27.root
  // uncertainties/systematicUncertaintyForPbPbAsymmetryRatio_2019-10-14.root
  
  // Create histogram managers for pp and PbPb
  DijetHistogramManager *ppHistograms[5];
  DijetHistogramManager *pbpbHistograms[5];
  for(int iFile = 0; iFile < 5; iFile++){
    if(iFile < nFilesPerDataset){
      ppHistograms[iFile] = new DijetHistogramManager(ppFile[iFile]);
      pbpbHistograms[iFile] = new DijetHistogramManager(pbpbFile[iFile]);
    } else {
      ppHistograms[iFile] = NULL;
      pbpbHistograms[iFile] = NULL;
    }
  }
  JffCorrector *ppUncertaintyProvider = new JffCorrector();
  JffCorrector *pbpbUncertaintyProvider = new JffCorrector();
  
  // Choose which figure sets to draw
  bool drawJetShapePptoPbPbRatio = false;
  bool drawJetShapeSymmetricAsymmetricRatio = true;
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = true;
  const char* figureFormat = "pdf";
  TString saveName;
  
  // Get the number of asymmetry bins
  const int nAsymmetryBins = pbpbHistograms[0]->GetNAsymmetryBins();
  const int nTrackPtBins = pbpbHistograms[0]->GetNTrackPtBins();
  const int nCentralityBins = pbpbHistograms[0]->GetNCentralityBins();
  
  double centralityBinBorders[] = {0,10,30,50,90};
  double asymmetryBinBorders[] = {0,0.6,0.8,1};
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};
 
  // Logarithmic scales for figures
  bool logPt = true;          // pT distributions
  bool logCorrelation = true; // track-jet deltaPhi-deltaEta distributions
  bool logJetShape = true;    // Jet shapes
  
  //int iAsymmetry = 0;
  int iJetTrack = DijetHistogramManager::kPtWeightedTrackSubleadingJet;  // DijetHistogramManager::kPtWeightedTrackSubleadingJet
  int lowestTrackPtBin = 0;
  
  TString asymmetryString[nAsymmetryBins+1];
  TString asymmetrySaveString[nAsymmetryBins+1];
  for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
    asymmetryString[iAsymmetry] = Form("   %.1f < x_{j} < %.1f",asymmetryBinBorders[iAsymmetry],asymmetryBinBorders[iAsymmetry+1]);
    asymmetrySaveString[iAsymmetry] = Form("_A=%.1f-%.1f",asymmetryBinBorders[iAsymmetry],asymmetryBinBorders[iAsymmetry+1]);
  }
  asymmetryString[nAsymmetryBins] = "";
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // Load only jet shape histograms from these
  for(int iFile = 0; iFile < nFilesPerDataset; iFile++){
    ppHistograms[iFile]->SetLoadTrackLeadingJetCorrelationsPtWeighted(true);
    ppHistograms[iFile]->SetAsymmetryBinRange(0,nAsymmetryBins);
    ppHistograms[iFile]->LoadProcessedHistograms();
  
    pbpbHistograms[iFile]->SetLoadTrackLeadingJetCorrelationsPtWeighted(true);
    pbpbHistograms[iFile]->SetAsymmetryBinRange(0,nAsymmetryBins);
    pbpbHistograms[iFile]->LoadProcessedHistograms();
  }
  
  // Load the systematic uncertainties from the files
  ppUncertaintyProvider->ReadSystematicFile(ppUncertaintyFile);
  pbpbUncertaintyProvider->ReadSystematicFile(pbpbUncertaintyFile);
  
  // Plot the figures using Xiao's plotting macro
  plotJetShapeXiao(nFilesPerDataset, ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kPtWeightedTrackLeadingJet);
  plotJetShapeXiao(nFilesPerDataset, ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kPtWeightedTrackSubleadingJet);
  
  /*
  // Temporary: Get the ratio between pp and PbPb jet shape
  TH1D *sumHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *sumHistogramPp[nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPp[nAsymmetryBins+1];
  TH1D *helperHistogram;
  double shapeIntegralPbPb[nCentralityBins][nAsymmetryBins+1];
  double shapeIntegralPp[nAsymmetryBins+1];
  
  // Read the pT summed jet shape histograms
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      sumHistogramPbPb[iCentrality][iAsymmetry] = (TH1D*) pbpbHistograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, lowestTrackPtBin)->Clone(Form("sumPbPb%d",iCentrality));
    }
    sumHistogramPp[iAsymmetry] = (TH1D*) ppHistograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, lowestTrackPtBin)->Clone("sumPp");
  }
  
  // Sum the pT:s
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iTrackPt = lowestTrackPtBin+1; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        helperHistogram = (TH1D*)pbpbHistograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,iCentrality,iTrackPt)->Clone();
        sumHistogramPbPb[iCentrality][iAsymmetry]->Add(helperHistogram);
      }
      helperHistogram = (TH1D*)ppHistograms[0]->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,0,iTrackPt)->Clone();
      sumHistogramPp[iAsymmetry]->Add(helperHistogram);
    } // track pT
  }
  
  // Normalize the jet shape histograms
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      shapeIntegralPbPb[iCentrality][iAsymmetry] = sumHistogramPbPb[iCentrality][iAsymmetry]->Integral(1,sumHistogramPbPb[iCentrality][iAsymmetry]->FindBin(0.99),"width");
      sumHistogramPbPb[iCentrality][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iCentrality][iAsymmetry]);
    }
    
    shapeIntegralPp[iAsymmetry] = sumHistogramPp[iAsymmetry]->Integral(1,sumHistogramPp[iAsymmetry]->FindBin(0.99),"width");
    sumHistogramPp[iAsymmetry]->Scale(1.0/shapeIntegralPp[iAsymmetry]);
    
    // Read the uncertainties for the pT summed jet shapes and scale them for jet shapes
    uncertaintyHistogramPp[iAsymmetry] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry);
    uncertaintyHistogramPp[iAsymmetry]->Scale(1.0/shapeIntegralPp[iAsymmetry]);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      uncertaintyHistogramPbPb[iCentrality][iAsymmetry] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, nTrackPtBins, iAsymmetry);
      uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iCentrality][iAsymmetry]);
    }
  }
  
  TH1D* ratioHistogram[nCentralityBins][nAsymmetryBins+1];
  TH1D* asymmetryRatioHistogram[nCentralityBins+1][nAsymmetryBins];
  
  // Create uncertainty histograms for the ratio using standard jet shape binning
  TH1D* ratioUncertainty[nCentralityBins][nAsymmetryBins+1];
  TH1D* asymmetryRatioUncertainty[nCentralityBins+1][nAsymmetryBins];
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      ratioUncertainty[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("ratioUncertainty%d%d", iCentrality, iAsymmetry));
    }
  }
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
      asymmetryRatioUncertainty[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[0][0]->Clone(Form("asymmetryRatioUncertainty%d%d", iCentrality, iAsymmetry));
    }
  }
  
  JDrawer *drawer = new JDrawer();
  drawer->SetCanvasSize(500,500);
  
  TLegend *legend;
  double ppValue, pbpbValue, ratioValue;
  double ppUncertainty, pbpbUncertainty, ratioUncertaintyValue;
  
  TLine *oneLine = new TLine(0,1,1,1);
  oneLine->SetLineStyle(2);
  
  if(drawJetShapePptoPbPbRatio){
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        ratioHistogram[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("ratio%d%d", iCentrality, iAsymmetry));
        ratioHistogram[iCentrality][iAsymmetry]->Divide(sumHistogramPp[iAsymmetry]);
        
        // Calculate the systemtic uncertainty for the ratio
        for(int iBin = 1; iBin < ratioUncertainty[iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
          ppValue = sumHistogramPp[iAsymmetry]->GetBinContent(iBin);
          pbpbValue = sumHistogramPbPb[iCentrality][iAsymmetry]->GetBinContent(iBin);
          ppUncertainty = uncertaintyHistogramPp[iAsymmetry]->GetBinContent(iBin);
          pbpbUncertainty = uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->GetBinContent(iBin);
          ratioValue = ratioHistogram[iCentrality][iAsymmetry]->GetBinContent(iBin);
          ratioUncertaintyValue = TMath::Sqrt(TMath::Power(pbpbUncertainty/ppValue,2)+TMath::Power(pbpbValue*ppUncertainty/TMath::Power(ppValue,2),2));
          ratioUncertainty[iCentrality][iAsymmetry]->SetBinContent(iBin,ratioValue);
          ratioUncertainty[iCentrality][iAsymmetry]->SetBinError(iBin,ratioUncertaintyValue);
        }
        
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetRangeUser(0,1);
        ratioHistogram[iCentrality][iAsymmetry]->GetYaxis()->SetRangeUser(0,3);
        ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->SetRangeUser(0,1);
        ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetRangeUser(0,3);
        ratioUncertainty[iCentrality][iAsymmetry]->SetFillColorAlpha(29,0.25);
        
        drawer->DrawHistogram(ratioUncertainty[iCentrality][iAsymmetry],"#Deltar","#rho(#Deltar)_{PbPb} / #rho(#Deltar)_{pp}", " ","E2");
        oneLine->Draw();
        ratioHistogram[iCentrality][iAsymmetry]->SetLineColor(kBlack);
        ratioHistogram[iCentrality][iAsymmetry]->Draw("same");
        
        legend = new TLegend(0.25,0.73,0.45,0.88);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->SetHeader(Form("C: %.0f-%.0f%s",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
        legend->AddEntry(ratioHistogram[iCentrality][iAsymmetry],Form("p_{T} > %.1f",trackPtBinBorders[lowestTrackPtBin]),"l");
        legend->Draw();
        
        if(saveFigures){
          if(iAsymmetry == nAsymmetryBins){
            gPad->GetCanvas()->SaveAs(Form("figures/jetShapeRatio_%s_C=%.0f-%.0f_pT%d.pdf", pbpbHistograms[0]->GetJetTrackHistogramName(iJetTrack), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1],lowestTrackPtBin));
          } else {
            saveName = Form("figures/jetShapeRatio_%s_A=%.1f-%.1f_C=%.0f-%.0f", pbpbHistograms[0]->GetJetTrackHistogramName(iJetTrack), asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1], centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]);
            saveName.ReplaceAll(".","v");
            gPad->GetCanvas()->SaveAs(Form("%s.%s",saveName.Data(),figureFormat));
          }
        } // Saving figures
        
      } // Centrality loop
    } // Asymmetry loop
  } // Ratio between pp and PbPb in different asymmetry bins
  
  if(drawJetShapeSymmetricAsymmetricRatio){
    for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
        
        if(iCentrality == nCentralityBins){
          asymmetryRatioHistogram[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPp[iAsymmetry]->Clone(Form("ratioAsymmetry%d%d", iCentrality, iAsymmetry));
          asymmetryRatioHistogram[iCentrality][iAsymmetry]->Divide(sumHistogramPp[nAsymmetryBins]);
        } else {
          asymmetryRatioHistogram[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("ratioAsymmetry%d%d", iCentrality, iAsymmetry));
          asymmetryRatioHistogram[iCentrality][iAsymmetry]->Divide(sumHistogramPbPb[iCentrality][nAsymmetryBins]);
        }
        
        // Calculate the systemtic uncertainty for the ratio
        for(int iBin = 1; iBin < asymmetryRatioUncertainty[iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
          if(iCentrality == nCentralityBins){
            ppValue = sumHistogramPp[nAsymmetryBins]->GetBinContent(iBin);
            pbpbValue = sumHistogramPp[iAsymmetry]->GetBinContent(iBin);
            ppUncertainty = uncertaintyHistogramPp[nAsymmetryBins]->GetBinContent(iBin);
            pbpbUncertainty = uncertaintyHistogramPp[iAsymmetry]->GetBinContent(iBin);
          } else {
            ppValue = sumHistogramPbPb[iCentrality][nAsymmetryBins]->GetBinContent(iBin);
            pbpbValue = sumHistogramPbPb[iCentrality][iAsymmetry]->GetBinContent(iBin);
            ppUncertainty = uncertaintyHistogramPbPb[iCentrality][nAsymmetryBins]->GetBinContent(iBin);
            pbpbUncertainty = uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->GetBinContent(iBin);
          }
          ratioValue = asymmetryRatioHistogram[iCentrality][iAsymmetry]->GetBinContent(iBin);
          ratioUncertaintyValue = TMath::Sqrt(TMath::Power(pbpbUncertainty/ppValue,2)+TMath::Power(pbpbValue*ppUncertainty/TMath::Power(ppValue,2),2));
          asymmetryRatioUncertainty[iCentrality][iAsymmetry]->SetBinContent(iBin,ratioValue);
          asymmetryRatioUncertainty[iCentrality][iAsymmetry]->SetBinError(iBin,ratioUncertaintyValue);
        }
        
        asymmetryRatioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetRangeUser(0,1);
        asymmetryRatioHistogram[iCentrality][iAsymmetry]->GetYaxis()->SetRangeUser(0,3);
        asymmetryRatioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->SetRangeUser(0,1);
        asymmetryRatioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetRangeUser(0,3);
        asymmetryRatioUncertainty[iCentrality][iAsymmetry]->SetFillColorAlpha(29,0.25);
        
        drawer->DrawHistogram(asymmetryRatioUncertainty[iCentrality][iAsymmetry],"#Deltar",Form("#rho(#Deltar) %.1f < x_{j} < %.1f / #rho(#Deltar) all x_{j}", asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1]), " ","E2");
        oneLine->Draw();
        asymmetryRatioHistogram[iCentrality][iAsymmetry]->SetLineColor(kBlack);
        asymmetryRatioHistogram[iCentrality][iAsymmetry]->Draw("same");
        
        legend = new TLegend(0.25,0.73,0.45,0.88);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        if(iCentrality == nCentralityBins){
          legend->SetHeader("pp  Track-subleading jet");
        } else {
          legend->SetHeader(Form("C: %.0f-%.0f  Track-subleading jet",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
        }
        legend->AddEntry(asymmetryRatioHistogram[iCentrality][iAsymmetry],Form("%.1f < x_{j} < %.1f / all x_{j}", asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1]),"l");
        legend->Draw();
        
        if(saveFigures){
          if(iCentrality == nCentralityBins){
            saveName = Form("figures/jetShapeAsymmetryRatio_%s_pp%s", pbpbHistograms[0]->GetJetTrackHistogramName(iJetTrack), asymmetrySaveString[iAsymmetry].Data());
            saveName.ReplaceAll(".","v");
            gPad->GetCanvas()->SaveAs(Form("%s.%s",saveName.Data(),figureFormat));
          } else {
            saveName = Form("figures/jetShapeAsymmetryRatio_%s_C=%.0f-%.0f%s", pbpbHistograms[0]->GetJetTrackHistogramName(iJetTrack), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetrySaveString[iAsymmetry].Data());
            saveName.ReplaceAll(".","v");
            gPad->GetCanvas()->SaveAs(Form("%s.%s",saveName.Data(),figureFormat));
          }
        } // Saving figures
        
      } // Centrality loop
    } // Asymmetry loop
  } // Draw asymmetry ratio
  */
}
