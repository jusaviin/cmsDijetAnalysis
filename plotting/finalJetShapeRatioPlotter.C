#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "stackHistXiao.h"
#include "xCanvas.h"
#include "JDrawer.h"

void plotJetShapeXiao(DijetHistogramManager *ppHistograms, DijetHistogramManager *pbpbHistograms, JffCorrector *ppUncertaintyProvider, JffCorrector *pbpbUncertaintyProvider, int iJetTrack){

  const int nTrackPtBins = pbpbHistograms->GetNTrackPtBins();
  const int nCentralityBins = pbpbHistograms->GetNCentralityBins();
  const int nAsymmetryBins = pbpbHistograms->GetNAsymmetryBins();
  const int nRelevantUncertainties = 4;
  
  const char* jetShapeTitle[] = {"Leading jet shape ratio between x_{j} bins","Subleading jet shape ratio between x_{j} bins","Jet shape ratio between x_{j} bins"};
  const char* jetShapeSaveName[] = {"trackLeadingJet","trackSubleadingJet","trackInclusiveJet"};
  const char* xjString[] = {"0.0 < x_{j} < 0.6","0.6 < x_{j} < 0.8","0.8 < x_{j} < 1.0","x_{j} inclusive"};
  const char* asymmetrySaveName[] = {"_A=0v0-0v6","_A=0v6-0v8","_A=0v8-1v0",""};
    
  // Temporary: Get the ratio between pp and PbPb jet shape
  TH1D *sumHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *sumHistogramPp[nAsymmetryBins+1];
  TH1D *jetShapeArray[nCentralityBins+1][nTrackPtBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPp[nAsymmetryBins+1];
  TH1D *uncertaintySummerPp[nAsymmetryBins+1][nRelevantUncertainties];
  TH1D *uncertaintySummerPbPb[nCentralityBins][nAsymmetryBins+1][nRelevantUncertainties];
  TH1D *sumUncertainty[nCentralityBins+1][nAsymmetryBins+1];
  TH1D *helperHistogram;
  double shapeIntegralPbPb[nCentralityBins][nAsymmetryBins+1];
  double shapeIntegralPp[nAsymmetryBins+1];
  double uncertaintySum;
  
  // Read the pT summed jet shape histograms
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      sumHistogramPbPb[iCentrality][iAsymmetry] = (TH1D*) pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, 0)->Clone(Form("sumPbPb%d",iCentrality));
      jetShapeArray[iCentrality][0][iAsymmetry] = pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, 0);
    }
    sumHistogramPp[iAsymmetry] = (TH1D*) ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, 0)->Clone("sumPp");
    jetShapeArray[nCentralityBins][0][iAsymmetry] = ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, 0);
  }
  
  // Sum the pT:s
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iTrackPt = 1; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        helperHistogram = (TH1D*)pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,iCentrality,iTrackPt)->Clone();
        sumHistogramPbPb[iCentrality][iAsymmetry]->Add(helperHistogram);
        jetShapeArray[iCentrality][iTrackPt][iAsymmetry] = pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,iCentrality,iTrackPt);
      }
      helperHistogram = (TH1D*)ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,0,iTrackPt)->Clone();
      sumHistogramPp[iAsymmetry]->Add(helperHistogram);
      jetShapeArray[nCentralityBins][iTrackPt][iAsymmetry] = ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, iTrackPt);
    } // track pT
  }
  
  // Normalize the jet shape histograms
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      shapeIntegralPbPb[iCentrality][iAsymmetry] = sumHistogramPbPb[iCentrality][iAsymmetry]->Integral(1,sumHistogramPbPb[iCentrality][iAsymmetry]->FindBin(0.99),"width");
      sumHistogramPbPb[iCentrality][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iCentrality][iAsymmetry]);
      
      // To draw the systematic uncertainties, clone the sum histograms to uncertainty histograms
      sumUncertainty[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("sumUncertaintyPbPb%d%d", iCentrality, iAsymmetry));
      
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        jetShapeArray[iCentrality][iTrackPt][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iCentrality][iAsymmetry]);
      }
      
    }
    
    shapeIntegralPp[iAsymmetry] = sumHistogramPp[iAsymmetry]->Integral(1,sumHistogramPp[iAsymmetry]->FindBin(0.99),"width");
    sumHistogramPp[iAsymmetry]->Scale(1.0/shapeIntegralPp[iAsymmetry]);
    
    // To draw the systematic uncertainties, clone the sum histograms to uncertainty histograms
    sumUncertainty[nCentralityBins][iAsymmetry] = (TH1D*) sumHistogramPp[iAsymmetry]->Clone(Form("sumUncertaintyPp%d", iAsymmetry));
    
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      jetShapeArray[nCentralityBins][iTrackPt][iAsymmetry]->Scale(1.0/shapeIntegralPp[iAsymmetry]);
    }
     
    // Read the uncertainties for the pT summed jet shapes but skip a few sources not relevant here
    uncertaintySummerPp[iAsymmetry][0] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry, JffCorrector::kBackgroundFluctuation);
    uncertaintySummerPp[iAsymmetry][1] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry, JffCorrector::kFragmentationBias);
    uncertaintySummerPp[iAsymmetry][2] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry, JffCorrector::kPairAcceptance);
    uncertaintySummerPp[iAsymmetry][3] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry, JffCorrector::kBackgroundSubtraction);
    
    uncertaintyHistogramPp[iAsymmetry] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry);
    
    
    // Add the relevant uncertainties in quadrature for pp
    for(int iBin = 1; iBin <= uncertaintyHistogramPp[iAsymmetry]->GetNbinsX(); iBin++){
      uncertaintySum = 0;
      for(int iUncertainty = 0; iUncertainty < nRelevantUncertainties; iUncertainty++){
        uncertaintySum += TMath::Power(uncertaintySummerPp[iAsymmetry][iUncertainty]->GetBinContent(iBin),2);
      }
      uncertaintyHistogramPp[iAsymmetry]->SetBinContent(iBin,TMath::Sqrt(uncertaintySum));
    }
    
    // Scale the histograms from jet shapes
    uncertaintyHistogramPp[iAsymmetry]->Scale(1.0/shapeIntegralPp[iAsymmetry]);
    
    // Set the bin errors in the sumUncertainty histograms to match the uncertainties read from the file
    for(int iBin = 1; iBin <= sumUncertainty[nAsymmetryBins][iAsymmetry]->GetNbinsX(); iBin++){
      sumUncertainty[nCentralityBins][iAsymmetry]->SetBinError(iBin, uncertaintyHistogramPp[iAsymmetry]->GetBinContent(iBin));
    }
    
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      // Read the uncertainties for the pT summed jet shapes but skip a few sources not relevant here
      uncertaintySummerPbPb[iCentrality][iAsymmetry][0] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, nTrackPtBins, iAsymmetry, JffCorrector::kBackgroundFluctuation);
      uncertaintySummerPbPb[iCentrality][iAsymmetry][1] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, nTrackPtBins, iAsymmetry, JffCorrector::kFragmentationBias);
      uncertaintySummerPbPb[iCentrality][iAsymmetry][2] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, nTrackPtBins, iAsymmetry, JffCorrector::kPairAcceptance);
      uncertaintySummerPbPb[iCentrality][iAsymmetry][3] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, nTrackPtBins, iAsymmetry, JffCorrector::kBackgroundSubtraction);
      
      uncertaintyHistogramPbPb[iCentrality][iAsymmetry] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, nTrackPtBins, iAsymmetry);
      
      // Add the relevant uncertainties in quadrature for PbPb
      for(int iBin = 1; iBin <= uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
        uncertaintySum = 0;
        for(int iUncertainty = 0; iUncertainty < nRelevantUncertainties; iUncertainty++){
          uncertaintySum += TMath::Power(uncertaintySummerPbPb[iCentrality][iAsymmetry][iUncertainty]->GetBinContent(iBin),2);
        }
        uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->SetBinContent(iBin,TMath::Sqrt(uncertaintySum));
      }
      
      uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iCentrality][iAsymmetry]);
      
      // Set the bin errors in the sumUncertainty histograms to match the uncertainties read from the file
      for(int iBin = 1; iBin <= sumUncertainty[iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
        sumUncertainty[iCentrality][iAsymmetry]->SetBinError(iBin, uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->GetBinContent(iBin));
      }
      
    }
  }

  // Create the asymmetric to nominal and symmetric to nominal ratios
  TH1D* asymmetryRatioHistogram[nCentralityBins+1][nAsymmetryBins];
  TH1D* asymmetryRatioUncertainty[nCentralityBins+1][nAsymmetryBins];
  double ppValue, pbpbValue, ratioValue;
  double ppUncertainty, pbpbUncertainty, ratioUncertaintyValue;
  
  for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      
      if(iCentrality == nCentralityBins){
        asymmetryRatioHistogram[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPp[iAsymmetry]->Clone(Form("ratioAsymmetry%d%d", iCentrality, iAsymmetry));
        asymmetryRatioHistogram[iCentrality][iAsymmetry]->Divide(sumHistogramPp[nAsymmetryBins]);
      } else {
        asymmetryRatioHistogram[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("ratioAsymmetry%d%d", iCentrality, iAsymmetry));
        asymmetryRatioHistogram[iCentrality][iAsymmetry]->Divide(sumHistogramPbPb[iCentrality][nAsymmetryBins]);
      }
      
      asymmetryRatioUncertainty[iCentrality][iAsymmetry] = (TH1D*) asymmetryRatioHistogram[iCentrality][iAsymmetry]->Clone(Form("asymmetryRatioUncertainty%d%d",iAsymmetry,iCentrality));
      
      // Calculate the systematic uncertainty for the ratio
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
      
    } // Centrality loop
  } // Asymmetry loop
  
  TString cent_lab[4] = {"0-10%", "10-30%", "30-50%", "50-90%"};
  
  // Create the stacked jet shape histograms and set a good drawing style for the uncertainty histgorams
  for(int iCentrality = 0; iCentrality < nCentralityBins+1; iCentrality++){
    for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
      asymmetryRatioUncertainty[iCentrality][iAsymmetry]->SetFillStyle(1001);
      asymmetryRatioUncertainty[iCentrality][iAsymmetry]->SetFillColorAlpha(kGray+3,0.4);
      asymmetryRatioUncertainty[iCentrality][iAsymmetry]->SetMarkerStyle(20);
      asymmetryRatioUncertainty[iCentrality][iAsymmetry]->SetMarkerSize(1.6);
      asymmetryRatioUncertainty[iCentrality][iAsymmetry]->SetMarkerColor(kBlack);
      asymmetryRatioUncertainty[iCentrality][iAsymmetry]->SetLineColor(kBlack);
    } // Asymmetry loop
  } // Centrality loop
  
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
  double axisZoom[] = {2.2,3.2,2.2};
  
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    bigCanvas->CD(5-iCentrality);
    
    asymmetryRatioHistogram[iCentrality][0]->SetTitle("");
    asymmetryRatioHistogram[iCentrality][0]->SetAxisRange(0, axisZoom[iJetTrack/3], "Y");
    asymmetryRatioHistogram[iCentrality][0]->SetAxisRange(0, 0.99, "X");
    if(iCentrality < 4)  {
      asymmetryRatioHistogram[iCentrality][0]->GetXaxis()->SetTitleOffset(0.7);
      asymmetryRatioHistogram[iCentrality][0]->GetXaxis()->SetTitleSize(0.11);
      asymmetryRatioHistogram[iCentrality][0]->GetXaxis()->SetNdivisions(505);
      asymmetryRatioHistogram[iCentrality][0]->GetXaxis()->SetLabelSize(0.09);
      asymmetryRatioHistogram[iCentrality][0]->GetYaxis()->SetNdivisions(505);
    }
    if(iCentrality == 4){
      asymmetryRatioHistogram[iCentrality][0]->GetXaxis()->SetTitleOffset(0.92);
      asymmetryRatioHistogram[iCentrality][0]->GetXaxis()->SetTitleSize(0.086);
      asymmetryRatioHistogram[iCentrality][0]->GetXaxis()->SetNdivisions(505);
      asymmetryRatioHistogram[iCentrality][0]->GetXaxis()->SetLabelOffset(0.016);
      asymmetryRatioHistogram[iCentrality][0]->GetXaxis()->SetLabelSize(0.064);
      asymmetryRatioHistogram[iCentrality][0]->GetYaxis()->SetNdivisions(505);
      asymmetryRatioHistogram[iCentrality][0]->GetYaxis()->SetLabelOffset(0.02);
      asymmetryRatioHistogram[iCentrality][0]->GetYaxis()->SetLabelSize(0.09);
      asymmetryRatioHistogram[iCentrality][0]->GetYaxis()->SetTitleOffset(1.0);
      asymmetryRatioHistogram[iCentrality][0]->GetYaxis()->SetTitleSize(0.1);
      asymmetryRatioHistogram[iCentrality][0]->GetYaxis()->SetTitle("#rho(#Deltar)_{Asymm}/#rho(#Deltar)_{All}");
      mainTitle->SetTextSize(0.073);
    }
    
    asymmetryRatioHistogram[iCentrality][0]->GetYaxis()->CenterTitle();
    asymmetryRatioHistogram[iCentrality][0]->GetXaxis()->CenterTitle();
    asymmetryRatioHistogram[iCentrality][0]->SetStats(0);
    asymmetryRatioHistogram[iCentrality][0]->Draw();
    
    line = new TLine();
    line->SetLineStyle(2);
    line->DrawLine(0, 1, 1, 1);
    
    asymmetryRatioUncertainty[iCentrality][0]->Draw("same e2");
    if(iCentrality == 4){
      mainTitle->SetTextFont(22);
      mainTitle->SetTextSize(.085);
      mainTitle->DrawLatexNDC(.35, 0.9, "pp reference");
    }
    else {
      mainTitle->SetTextFont(22);
      mainTitle->SetTextSize(.09);
      mainTitle->DrawLatexNDC(.19, 0.9, "PbPb");
      mainTitle->DrawLatexNDC(.19, .82, cent_lab[iCentrality]);
    }
    bigCanvas->CD(10-iCentrality);
    asymmetryRatioHistogram[iCentrality][2]->SetTitle("");
    asymmetryRatioHistogram[iCentrality][2]->GetXaxis()->SetTitleOffset(1.1);
    asymmetryRatioHistogram[iCentrality][2]->GetXaxis()->SetTitle("#Deltar");
    asymmetryRatioHistogram[iCentrality][2]->SetAxisRange(0, axisZoom[iJetTrack/3], "Y");
    asymmetryRatioHistogram[iCentrality][2]->SetAxisRange(0, 0.99, "X");
    if( iCentrality<4 )  {
      asymmetryRatioHistogram[iCentrality][2]->GetXaxis()->SetTitleOffset(0.65);
      asymmetryRatioHistogram[iCentrality][2]->GetXaxis()->SetTitleSize(0.115);
      asymmetryRatioHistogram[iCentrality][2]->GetXaxis()->SetNdivisions(505);
      asymmetryRatioHistogram[iCentrality][2]->GetXaxis()->SetLabelOffset(0.00001);
      asymmetryRatioHistogram[iCentrality][2]->GetXaxis()->SetLabelSize(0.092);
      asymmetryRatioHistogram[iCentrality][2]->GetYaxis()->SetNdivisions(505);
    }
    if(iCentrality==4 ){
      asymmetryRatioHistogram[iCentrality][2]->GetXaxis()->SetTitleOffset(0.92);
      asymmetryRatioHistogram[iCentrality][2]->GetXaxis()->SetTitleSize(0.086);
      asymmetryRatioHistogram[iCentrality][2]->GetXaxis()->SetNdivisions(505);
      asymmetryRatioHistogram[iCentrality][2]->GetXaxis()->SetLabelOffset(0.016);
      asymmetryRatioHistogram[iCentrality][2]->GetXaxis()->SetLabelSize(0.071);
      asymmetryRatioHistogram[iCentrality][2]->GetYaxis()->SetNdivisions(505);
      asymmetryRatioHistogram[iCentrality][2]->GetYaxis()->SetLabelOffset(0.02);
      asymmetryRatioHistogram[iCentrality][2]->GetYaxis()->SetLabelSize(0.07);
      asymmetryRatioHistogram[iCentrality][2]->GetYaxis()->SetTitleOffset(1.2);
      asymmetryRatioHistogram[iCentrality][2]->GetYaxis()->SetTitleSize(0.08);
      asymmetryRatioHistogram[iCentrality][2]->GetYaxis()->SetTitle("#rho(#Deltar)_{Symm}/#rho(#Deltar)_{All}");
      mainTitle->SetTextSize(0.073);
    }

    
    asymmetryRatioHistogram[iCentrality][2]->GetYaxis()->CenterTitle();
    asymmetryRatioHistogram[iCentrality][2]->GetXaxis()->CenterTitle();
    asymmetryRatioHistogram[iCentrality][2]->SetStats(0);
    asymmetryRatioHistogram[iCentrality][2]->Draw("same");
    
    line = new TLine();
    line->SetLineStyle(2);
    line->DrawLine(0, 1, 1, 1);
    
    asymmetryRatioUncertainty[iCentrality][2]->Draw("same e2");
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
  
  double cmsPosition[] = {0.42,0.4,0.42};
  double jetShapeTitlePosition[] = {0.18,0.16,0.18};
  double xjPosition[] = {0.76,0.77,0.76};
  
  bigCanvas->cd(0);
  mainTitle->SetTextFont(62);
  mainTitle->SetTextSize(0.065);
  mainTitle->DrawLatexNDC(0.07, 0.9, "CMS");
  
  mainTitle->SetTextFont(42);
  mainTitle->SetTextSize(0.045);
  mainTitle->DrawLatexNDC(jetShapeTitlePosition[iJetTrack/3], 0.905, jetShapeTitle[iJetTrack/3]);
  
  //mainTitle->SetTextFont(42);
  //mainTitle->SetTextSize(0.035);
  //mainTitle->DrawLatexNDC(xjPosition[iJetTrack/3], 0.94, xjString[iAsymmetry]);
  
  mainTitle->SetTextSize(0.03);
  mainTitle->DrawLatexNDC(0.61, 0.94, "pp 320 pb^{-1} (5.02 TeV)  PbPb 1.7 nb^{-1} (5.02 TeV)");
  mainTitle->SetTextSize(0.028);
  mainTitle->DrawLatexNDC(0.57, 0.89, "anti-k_{T} R = 0.4, |#eta_{jet}| < 1.6, p_{T,1} > 120 GeV, p_{T,2} > 50 GeV, #Delta#phi_{1,2} > #frac{5#pi}{6}");
  //  lb->drawText("(p_{T}> 120 GeV, |#eta_{jet}|<1.6)", 0.2, 0.25, 4);
  
  box = new TBox();
  box->SetFillColor(kWhite);
  bigCanvas->cd(0);
  box->DrawBox(0.24,0.047, 0.28, 0.09);
  mainTitle->SetTextSize(0.034);
  mainTitle->DrawLatex(0.252, 0.065, "0");
  box->DrawBox(0.42,0.047, 0.45, 0.09);
  box->DrawBox(0.61,0.047, 0.64, 0.09);
  box->DrawBox(0.79,0.047, 0.82, 0.09);
  box->DrawBox(0.97,0.047, 0.99, 0.09);
  box->DrawBox(0.04,0.46, 0.065, 0.5);
  
  mainTitle->DrawLatex(0.436, 0.065, "0");
  mainTitle->DrawLatex(0.62, 0.065, "0");
  mainTitle->DrawLatex(0.804, 0.065, "0");
  mainTitle->DrawLatex(0.983, 0.065, "1");
  
  //bigCanvas->SaveAs("js_dr_normal_new.eps");
  bigCanvas->SaveAs(Form("figures/finalJetShapeAsymmetryRatio_%s.pdf",jetShapeSaveName[iJetTrack/3]));
  //bigCanvas->SaveAs("js_dr_normal_v3.eps");
  //bigCanvas->SaveAs("js_dr_normal_v3.pdf");
  
}

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void finalJetShapeRatioPlotter(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Open data files for pp and PbPb data
  TFile *ppFile = TFile::Open("data/ppData2017_highForest_pfJets_20EventsMixed_xjBins_finalTrackCorr_JECv4_wtaAxis_allCorrections_processed_2019-09-28.root");
  // data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_processed_2019-08-13.root
  // data/dijet_pp_highForest_pfJets_noUncOrInc_allCorrections_wtaAxis_processed_2019-07-13.root
  TFile *pbpbFile = TFile::Open("data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_allCorrections_processed_2019-10-01_fiveJobsMissing.root");
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_calo80Trigger_xjBins_wtaAxis_allCorrections_JECv5b_processed_2019-09-05.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_modifiedSeagull_noErrorJff_averagePeakMixing_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_testNoJffCorr_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_allCorrections_modifiedSeagull_wtaAxis_JECv4_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root
  
  TFile *ppUncertaintyFile = TFile::Open("systematicUncertaintyForPp_20percentSpillJff_2019-09-30.root");
  TFile *pbpbUncertaintyFile = TFile::Open("systematicUncertaintyForPp_15percentSpill20Jff_2019-10-01.root");
  
  // Create histogram managers for pp and PbPb
  DijetHistogramManager *ppHistograms = new DijetHistogramManager(ppFile);
  DijetHistogramManager *pbpbHistograms = new DijetHistogramManager(pbpbFile);
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
  const int nAsymmetryBins = pbpbHistograms->GetNAsymmetryBins();
  const int nTrackPtBins = pbpbHistograms->GetNTrackPtBins();
  const int nCentralityBins = pbpbHistograms->GetNCentralityBins();
  
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
  ppHistograms->SetLoadTrackLeadingJetCorrelationsPtWeighted(true);
  ppHistograms->SetAsymmetryBinRange(0,nAsymmetryBins);
  ppHistograms->LoadProcessedHistograms();
  
  pbpbHistograms->SetLoadTrackLeadingJetCorrelationsPtWeighted(true);
  pbpbHistograms->SetAsymmetryBinRange(0,nAsymmetryBins);
  pbpbHistograms->LoadProcessedHistograms();
  
  // Load the systematic uncertainties from the files
  ppUncertaintyProvider->ReadSystematicFile(ppUncertaintyFile);
  pbpbUncertaintyProvider->ReadSystematicFile(pbpbUncertaintyFile);
  
  // Plot the figures using Xiao's plotting macro
  plotJetShapeXiao(ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kPtWeightedTrackLeadingJet);
  plotJetShapeXiao(ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kPtWeightedTrackSubleadingJet);
  return;
  
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
      sumHistogramPbPb[iCentrality][iAsymmetry] = (TH1D*) pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, lowestTrackPtBin)->Clone(Form("sumPbPb%d",iCentrality));
    }
    sumHistogramPp[iAsymmetry] = (TH1D*) ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, lowestTrackPtBin)->Clone("sumPp");
  }
  
  // Sum the pT:s
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iTrackPt = lowestTrackPtBin+1; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        helperHistogram = (TH1D*)pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,iCentrality,iTrackPt)->Clone();
        sumHistogramPbPb[iCentrality][iAsymmetry]->Add(helperHistogram);
      }
      helperHistogram = (TH1D*)ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,0,iTrackPt)->Clone();
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
            gPad->GetCanvas()->SaveAs(Form("figures/jetShapeRatio_%s_C=%.0f-%.0f_pT%d.pdf", pbpbHistograms->GetJetTrackHistogramName(iJetTrack), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1],lowestTrackPtBin));
          } else {
            saveName = Form("figures/jetShapeRatio_%s_A=%.1f-%.1f_C=%.0f-%.0f", pbpbHistograms->GetJetTrackHistogramName(iJetTrack), asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1], centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]);
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
            saveName = Form("figures/jetShapeAsymmetryRatio_%s_pp%s", pbpbHistograms->GetJetTrackHistogramName(iJetTrack), asymmetrySaveString[iAsymmetry].Data());
            saveName.ReplaceAll(".","v");
            gPad->GetCanvas()->SaveAs(Form("%s.%s",saveName.Data(),figureFormat));
          } else {
            saveName = Form("figures/jetShapeAsymmetryRatio_%s_C=%.0f-%.0f%s", pbpbHistograms->GetJetTrackHistogramName(iJetTrack), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetrySaveString[iAsymmetry].Data());
            saveName.ReplaceAll(".","v");
            gPad->GetCanvas()->SaveAs(Form("%s.%s",saveName.Data(),figureFormat));
          }
        } // Saving figures
        
      } // Centrality loop
    } // Asymmetry loop
  } // Draw asymmetry ratio
  
}
