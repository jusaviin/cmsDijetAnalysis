#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "stackHistXiao.h"
#include "xCanvas.h"
#include "JDrawer.h"

void plotJetShapeBigAsymmetry(DijetHistogramManager *ppHistograms, DijetHistogramManager *pbpbHistograms, JffCorrector *ppUncertaintyProvider, JffCorrector *pbpbUncertaintyProvider, int iJetTrack){

  const int nTrackPtBins = pbpbHistograms->GetNTrackPtBins();
  const int nCentralityBins = pbpbHistograms->GetNCentralityBins();
  const int nAsymmetryBins = pbpbHistograms->GetNAsymmetryBins();
  
  const char* jetShapeTitle[] = {"Leading jet shape","Subleading jet shape","Jet shape"};
  const char* jetShapeSaveName[] = {"trackLeadingJet","trackSubleadingJet","trackInclusiveJet"};
  const char* xjString[] = {"0.0 < x_{j} < 0.6","0.6 < x_{j} < 0.8","0.8 < x_{j} < 1.0","x_{j} inclusive"};
  const char* asymmetrySaveName[] = {"_A=0v0-0v6","_A=0v6-0v8","_A=0v8-1v0",""};
  const char* saveString = {"JetShape"};
  
  bool drawUncertainties = true;
  bool monteCarloLabels = false;
  bool normalizeJetShape = true;           // True: draw rho. False: draw P.
  
  // Change the titles if the jet shape is not normalized to one
  if(!normalizeJetShape){
    jetShapeTitle[0] = "Leading jet radial momentum distribution";
    jetShapeTitle[1] = "Subleading jet radial momentum distribution";
    jetShapeTitle[2] = "Jet radial momentum distribution";
    saveString = "JetRadialMomentum";
  }
  
  // Temporary: Get the ratio between pp and PbPb jet shape
  TH1D *sumHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *sumHistogramPp[nAsymmetryBins+1];
  TH1D *jetShapeArray[nCentralityBins+1][nTrackPtBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPp[nAsymmetryBins+1];
  TH1D *sumUncertainty[nCentralityBins+1][nAsymmetryBins+1];
  TH1D *helperHistogram;
  double shapeIntegralPbPb[nCentralityBins][nAsymmetryBins+1];
  double shapeIntegralPp[nAsymmetryBins+1];
  
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
      if(normalizeJetShape) sumHistogramPbPb[iCentrality][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iCentrality][iAsymmetry]);
      
      // To draw the systematic uncertainties, clone the sum histograms to uncertainty histograms
      sumUncertainty[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("sumUncertaintyPbPb%d%d", iCentrality, iAsymmetry));
      
      if(normalizeJetShape){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          jetShapeArray[iCentrality][iTrackPt][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iCentrality][iAsymmetry]);
        }
      }
      
    }
    
    shapeIntegralPp[iAsymmetry] = sumHistogramPp[iAsymmetry]->Integral(1,sumHistogramPp[iAsymmetry]->FindBin(0.99),"width");
    if(normalizeJetShape) sumHistogramPp[iAsymmetry]->Scale(1.0/shapeIntegralPp[iAsymmetry]);
    
    // To draw the systematic uncertainties, clone the sum histograms to uncertainty histograms
    sumUncertainty[nCentralityBins][iAsymmetry] = (TH1D*) sumHistogramPp[iAsymmetry]->Clone(Form("sumUncertaintyPp%d", iAsymmetry));
    
    if(normalizeJetShape){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        jetShapeArray[nCentralityBins][iTrackPt][iAsymmetry]->Scale(1.0/shapeIntegralPp[iAsymmetry]);
      }
    }
    
    // Read the uncertainties for the pT summed jet shapes and scale them for jet shapes
    uncertaintyHistogramPp[iAsymmetry] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry);
    if(normalizeJetShape) uncertaintyHistogramPp[iAsymmetry]->Scale(1.0/shapeIntegralPp[iAsymmetry]);
    
    // Set the bin errors in the sumUncertainty histograms to match the uncertainties read from the file
    for(int iBin = 1; iBin <= sumUncertainty[nAsymmetryBins][iAsymmetry]->GetNbinsX(); iBin++){
      sumUncertainty[nCentralityBins][iAsymmetry]->SetBinError(iBin, uncertaintyHistogramPp[iAsymmetry]->GetBinContent(iBin));
    }
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      uncertaintyHistogramPbPb[iCentrality][iAsymmetry] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, nTrackPtBins, iAsymmetry);
      if(normalizeJetShape) uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iCentrality][iAsymmetry]);
      
      // Set the bin errors in the sumUncertainty histograms to match the uncertainties read from the file
      for(int iBin = 1; iBin <= sumUncertainty[iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
        
        sumUncertainty[iCentrality][iAsymmetry]->SetBinError(iBin, uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->GetBinContent(iBin));
      }
      
    }
  }
  
  TH1D* ratioHistogram[nCentralityBins][nAsymmetryBins+1];
  
  // Create uncertainty histograms for the ratio using standard jet shape binning
  TH1D* ratioUncertainty[nCentralityBins][nAsymmetryBins+1];
  double ppValue, pbpbValue, ratioValue;
  double ppUncertainty, pbpbUncertainty, ratioUncertaintyValue;
  
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      //ratioUncertainty[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("ratioUncertainty%d%d", iCentrality, iAsymmetry));
      ratioHistogram[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("ratio%d%d", iCentrality, iAsymmetry));
      ratioHistogram[iCentrality][iAsymmetry]->Divide(sumHistogramPp[iAsymmetry]);
      
      ratioUncertainty[iCentrality][iAsymmetry] = (TH1D*) sumUncertainty[iCentrality][iAsymmetry]->Clone(Form("uncertaintyOfRatio%d%d", iCentrality, iAsymmetry));
      ratioUncertainty[iCentrality][iAsymmetry]->Divide(sumUncertainty[nCentralityBins][iAsymmetry]);
      
    }
  }
  
  TString cent_lab[4] = {"0-10%", "10-30%", "30-50%", "50-90%"};
  stackHist *jetShapeStack[nCentralityBins+1][nAsymmetryBins+1];
  
  // Create the stacked jet shape histograms and set a good drawing style for the uncertainty histgorams
  for(int iCentrality = 0; iCentrality < nCentralityBins+1; iCentrality++){
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      
      sumUncertainty[iCentrality][iAsymmetry]->SetFillStyle(1001);
      sumUncertainty[iCentrality][iAsymmetry]->SetFillColorAlpha(kGray+3,0.4);
      sumUncertainty[iCentrality][iAsymmetry]->SetMarkerStyle(20);
      sumUncertainty[iCentrality][iAsymmetry]->SetMarkerSize(1.6);
      sumUncertainty[iCentrality][iAsymmetry]->SetMarkerColor(kBlack);
      sumUncertainty[iCentrality][iAsymmetry]->SetLineColor(kBlack);
      
      jetShapeStack[iCentrality][iAsymmetry] = new stackHist(Form("st_%d",iCentrality));
      jetShapeStack[iCentrality][iAsymmetry]->setRange(0.0, 0.99, "x");
      if(normalizeJetShape){
        jetShapeStack[iCentrality][iAsymmetry]->setRange(0.005, 30, "y");
      } else {
        jetShapeStack[iCentrality][iAsymmetry]->setRange(0.5, 3000, "y");
      }
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins ; iTrackPt++){
        jetShapeStack[iCentrality][iAsymmetry]->addHist((TH1*) jetShapeArray[iCentrality][iTrackPt][iAsymmetry]);
      }
      //js_dr_err_all[iCentrality]->Scale(1.0/fac);
      //jetShapeSum[iCentrality]->SetMarkerColor(0);
    } // Asymmetry loop
  } // Centrality loop
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins+1; iAsymmetry++){
      ratioUncertainty[iCentrality][iAsymmetry]->SetFillStyle(1001);
      ratioUncertainty[iCentrality][iAsymmetry]->SetFillColorAlpha(kGray+3, 0.4);
      if(monteCarloLabels){
        ratioUncertainty[iCentrality][iAsymmetry]->SetMarkerColor(kBlack);
      } else {
        ratioUncertainty[iCentrality][iAsymmetry]->SetMarkerColor(0);
      }
      ratioUncertainty[iCentrality][iAsymmetry]->SetMarkerStyle(21);
      ratioUncertainty[iCentrality][iAsymmetry]->SetLineColor(kBlack);
    }
  }
  
  // ==============================
  // **  Draw the distributions  **
  // ==============================
  
  // Draw all the distributions to big canvases
  auxi_canvas *bigCanvas = new auxi_canvas(Form("theReallyBigCanvas%d",iJetTrack), "", 2500, 4000);
  bigCanvas->SetMargin(0.06, 0.01, 0.04, 0.16); // Margin order: Left, Right, Bottom, Top
  bigCanvas->divide(6,5);
  
  TBox *box = new TBox();
  TLatex *mainTitle = new TLatex();
  
  TLine *line[nAsymmetryBins];
  const char *pbpbLabel = "PbPb";
  const char *ppLabel = "pp";
  
  if(monteCarloLabels){
    pbpbLabel = "Pythia+Hydjet";
    ppLabel = "Pythia8";
  }
  
  for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
    
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      
      //cout << "In the loop iAsymmetry: " << iAsymmetry << " iCentrality: " << iCentrality << endl;
      
      bigCanvas->CD(5+iAsymmetry*10-iCentrality);
      gPad->SetLogy();
      jetShapeStack[iCentrality][iAsymmetry]->drawStack();
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetNdivisions(505);
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetLabelSize(0.09);
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitleOffset(1.2);
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitleSize(0.1);
      if(normalizeJetShape){
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitle("#rho(#Deltar)");
      } else {
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitle("#Rho(#Deltar)");
      }
      jetShapeStack[iCentrality][iAsymmetry]->hst->Draw();

      if(drawUncertainties) sumUncertainty[iCentrality][iAsymmetry]->Draw("same e2");
      if(iCentrality == nCentralityBins ){
        mainTitle->SetTextFont(22);
        mainTitle->SetTextSize(0.1);
        mainTitle->DrawLatexNDC(0.35, 0.88, ppLabel);
      }
      else {
        mainTitle->SetTextFont(22);
        mainTitle->SetTextSize(0.11);
        mainTitle->DrawLatexNDC(0.19, 0.9, pbpbLabel);
        mainTitle->DrawLatexNDC(0.19, 0.80, cent_lab[iCentrality]);
      }
      
      if(iCentrality == nCentralityBins){
        bigCanvas->CD(10+iAsymmetry*10-iCentrality);
        bigCanvas->pad[9+iAsymmetry*10-iCentrality]->SetLeftMargin(0.99999);
        
        helperHistogram = (TH1D*) ratioHistogram[iCentrality-1][iAsymmetry]->Clone(Form("myHelper%d%d%d",iJetTrack,iCentrality,iAsymmetry));
        helperHistogram->GetXaxis()->SetTitleOffset(2);
        helperHistogram->GetXaxis()->SetLabelOffset(2);
        helperHistogram->GetXaxis()->SetLabelSize(0.064);
        helperHistogram->GetYaxis()->SetNdivisions(505);
        helperHistogram->GetYaxis()->SetLabelOffset(0.02);
        helperHistogram->GetYaxis()->SetLabelSize(0.08);
        helperHistogram->GetYaxis()->SetTitleOffset(0.9);
        helperHistogram->GetYaxis()->SetTitleSize(0.1);
        
        if(iAsymmetry == 2){
          helperHistogram->GetYaxis()->SetTitleOffset(0.9);
          helperHistogram->GetYaxis()->SetTitleSize(0.09);
        }
        
        if(normalizeJetShape){
          helperHistogram->GetYaxis()->SetTitle("#rho(#Deltar)_{PbPb}/#rho(#Deltar)_{pp}");
        } else {
          helperHistogram->GetYaxis()->SetTitle("#Rho(#Deltar)_{PbPb}/#Rho(#Deltar)_{pp}");
        }
        
        helperHistogram->GetYaxis()->CenterTitle();
        helperHistogram->GetXaxis()->CenterTitle();
        helperHistogram->SetStats(0);
        helperHistogram->Draw("same");
        
      }
      
      // The ratios are only drawn to four last columns
      if(iCentrality < nCentralityBins){
        
        bigCanvas->CD(10+iAsymmetry*10-iCentrality);
        //ratio[iCentrality]->GetYaxis()->SetNdivisions(505);
        ratioHistogram[iCentrality][iAsymmetry]->SetTitle("");
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetTitleOffset(0.5);
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetTitle("#Deltar");
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetTitleSize(0.15);
        if(monteCarloLabels){
          ratioHistogram[iCentrality][iAsymmetry]->SetAxisRange(0.0, 2.2, "Y");
        } else if(normalizeJetShape){
          ratioHistogram[iCentrality][iAsymmetry]->SetAxisRange(0, 3.2, "Y");
        } else {
          ratioHistogram[iCentrality][iAsymmetry]->SetAxisRange(0, 4.2, "Y");
        }
        ratioHistogram[iCentrality][iAsymmetry]->SetAxisRange(0, 0.99, "X");
        
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetNdivisions(505);
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetLabelOffset(-0.015);
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetLabelSize(0.1);
        ratioHistogram[iCentrality][iAsymmetry]->GetYaxis()->SetNdivisions(505);
        
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->CenterTitle();
        ratioHistogram[iCentrality][iAsymmetry]->SetStats(0);
        ratioHistogram[iCentrality][iAsymmetry]->Draw("same");
        
        line[iAsymmetry] = new TLine();
        line[iAsymmetry]->SetLineStyle(2);
        line[iAsymmetry]->DrawLine(0, 1, 1, 1);
        
        if(drawUncertainties) ratioUncertainty[iCentrality][iAsymmetry]->Draw("same e2");
        
      }
      
    } // Centrality loop
    
  } // Asymmetry loop for drawing distributions to the really big canvas
  
  TLegend* ptLegend1 = new TLegend(0.04, 0.85, 0.27, 0.90);
  TLegend* ptLegend2 = new TLegend(0.28, 0.85, 0.51, 0.90);
  TLegend* ptLegend3 = new TLegend(0.53 ,0.85, 0.76, 0.90);
  TLegend* ptLegend4 = new TLegend(0.77, 0.85, 1.00, 0.90);
  ptLegend1->SetTextSize(0.02);
  ptLegend1->SetLineColor(kWhite);
  ptLegend1->SetFillColor(kWhite);
  ptLegend2->SetTextSize(0.02);
  ptLegend2->SetLineColor(kWhite);
  ptLegend2->SetFillColor(kWhite);
  ptLegend3->SetTextSize(0.02);
  ptLegend3->SetLineColor(kWhite);
  ptLegend3->SetFillColor(kWhite);
  ptLegend4->SetTextSize(0.02);
  ptLegend4->SetLineColor(kWhite);
  ptLegend4->SetFillColor(kWhite);
  
  TLegend* lt4 = new TLegend(0.25, 0.85, 0.85, 0.95);
  lt4->SetTextSize(0.06);
  lt4->SetLineColor(kWhite);
  lt4->SetFillColor(kWhite);
  
  
  
  ptLegend1->AddEntry(jetShapeStack[4][0]->hist_trunk.at(0), "0.7 < p_{T}^{trk}< 1 GeV","f");
  ptLegend1->AddEntry(jetShapeStack[4][0]->hist_trunk.at(1), "1 < p_{T}^{trk}< 2 GeV","f");
  
  ptLegend2->AddEntry(jetShapeStack[4][0]->hist_trunk.at(2), "2 < p_{T}^{trk}< 3 GeV","f");
  ptLegend2->AddEntry(jetShapeStack[4][0]->hist_trunk.at(3), "3 < p_{T}^{trk}< 4 GeV","f");
  
  ptLegend3->AddEntry(jetShapeStack[4][0]->hist_trunk.at(4), "4 < p_{T}^{trk}< 8 GeV","f");
  ptLegend3->AddEntry(jetShapeStack[4][0]->hist_trunk.at(5), "8 < p_{T}^{trk}< 12 GeV","f");
  
  ptLegend4->AddEntry(jetShapeStack[4][0]->hist_trunk.at(6), "12 < p_{T}^{trk}< 300 GeV","f");
  ptLegend4->AddEntry(sumUncertainty[4][0], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
  
  lt4->AddEntry(ratioUncertainty[0][0], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
  //bigCanvas->CD(9);
  //lt4->Draw();
  
  //bigCanvas->CD(1);
  //TLegend* lt6 = new TLegend(0.3,0.7,0.98,0.85);
  //lt6->SetTextSize(0.07);
  //lt6->SetLineColor(kWhite);
  //lt6->SetFillColor(kWhite);
  //lt6->AddEntry(js_dr_err_all[0], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
  //lt6->Draw();
  
  //bigCanvas->CD(2);
  //line[0]->SetLineStyle(1);
  //line[0]->DrawLineNDC(0, 0, 0, 1);

  double cmsPosition[] = {0.22,0.2,0.22};
  double cmsPositionY = 0.97;
  double jetShapeTitlePosition[] = {0.56,0.54,0.56};
  double xjPosition[] = {0.76,0.77,0.76};
  double systemPosition = 0.26;
  double selectionPosition = 0.205;
  double cmsSize = 0.035;
  
  // If we draw radial momentum distribution instead of jet shape, need to change and reposition the labels
  if(!normalizeJetShape){
    cmsPosition[iJetTrack/3] = 0.33;
    jetShapeTitlePosition[iJetTrack/3] -= 0.17;
    xjPosition[iJetTrack/3] += 0.08;
    cmsPositionY = 0.88;
    systemPosition += 0.04;
    selectionPosition += 0.04;
    cmsSize = 0.04;
  }
  
  // Draw the labels for different xj bins
  bigCanvas->CD(6);
  mainTitle->SetTextFont(42);
  mainTitle->SetTextSize(0.13);
  mainTitle->DrawLatexNDC(0.08,0.44,"0.0 < x_{j} < 0.6");
  
  bigCanvas->CD(16);
  mainTitle->DrawLatexNDC(0.08,0.44,"0.6 < x_{j} < 0.8");
  
  bigCanvas->CD(26);
  mainTitle->SetTextSize(0.12);
  mainTitle->DrawLatexNDC(0.09,0.55,"0.8 < x_{j} < 1.0");
  
  bigCanvas->cd(0);
  
  ptLegend1->Draw();
  ptLegend2->Draw();
  ptLegend3->Draw();
  ptLegend4->Draw();
  
  mainTitle->SetTextFont(62);
  mainTitle->SetTextSize(cmsSize);
  mainTitle->DrawLatexNDC(cmsPosition[iJetTrack/3], cmsPositionY, "CMS Preliminary");
  
  mainTitle->SetTextFont(42);
  mainTitle->SetTextSize(0.035);
  mainTitle->DrawLatexNDC(jetShapeTitlePosition[iJetTrack/3], 0.97, jetShapeTitle[iJetTrack/3]);
  
  //mainTitle->SetTextFont(42);
  //mainTitle->SetTextSize(0.035);
  //mainTitle->DrawLatexNDC(xjPosition[iJetTrack/3], 0.94, xjString[iAsymmetry]);
  
  mainTitle->SetTextSize(0.03);
  mainTitle->DrawLatexNDC(systemPosition, 0.94, "5.02 TeV   pp 320 pb^{-1}   PbPb 1.7 nb^{-1}");
  mainTitle->SetTextSize(0.022);
  mainTitle->DrawLatexNDC(selectionPosition, 0.915, "anti-k_{T} R = 0.4, |#eta_{jet}| < 1.6, p_{T,1} > 120 GeV, p_{T,2} > 50 GeV, #Delta#phi_{1,2} > #frac{5#pi}{6}");
  //  lb->drawText("(p_{T}> 120 GeV, |#eta_{jet}|<1.6)", 0.2, 0.25, 4);
  
  box = new TBox();
  box->SetFillColor(kWhite);
  bigCanvas->cd(0);
  box->DrawBox(0.24,0.017, 0.253, 0.039);
  mainTitle->SetTextSize(0.02);
  mainTitle->DrawLatex(0.242, 0.0285, "0");
  box->DrawBox(0.42,0.017, 0.437, 0.039);
  box->DrawBox(0.605,0.017, 0.63, 0.039);
  box->DrawBox(0.78,0.017, 0.81, 0.039);
  box->DrawBox(0.98,0.017, 0.99, 0.039);
  
  box->DrawBox(0.23,0.3067, 0.243, 0.315);
  box->DrawBox(0.23,0.5735, 0.243, 0.58);
  
  mainTitle->DrawLatex(0.428, 0.0285, "0");
  mainTitle->DrawLatex(0.614, 0.0285, "0");
  mainTitle->DrawLatex(0.801, 0.0285, "0");
  mainTitle->DrawLatex(0.982, 0.0285, "1");
  
  //bigCanvas->SaveAs("js_dr_normal_new.eps");
  bigCanvas->SaveAs(Form("figures/final%s_%s_DesDng.png", saveString, jetShapeSaveName[iJetTrack/3]));
  //bigCanvas->SaveAs("js_dr_normal_v3.eps");
  //bigCanvas->SaveAs("js_dr_normal_v3.pdf");
  
}

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void finalBigAsymmetryPlotter(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Open data files for pp and PbPb data
  TFile *ppFile = TFile::Open("data/ppData2017_highForest_pfJets_20EventsMixed_finalTrackCorr_xjBins_JECv4_wtaAxis_tunedSeagull_allCorrections_processed_2019-10-17.root");
  // data/ppData2017_highForest_pfJets_20EventsMixed_finalTrackCorr_xjBins_JECv4_wtaAxis_tunedSeagull_allCorrections_processed_2019-10-17.root"
  // data/ppData2017_highForest_pfJets_20EveMixed_finalTrackCorr_JECv4_eschemeAxis_tunedSeagull_allCorrections_processed_2019-10-14.root
  // data/ppData2017_highForest_pfJets_20EveMixed_finalTrackCorr_JECv4_eschemeAxis_allCorrections_processed_2019-10-02.root
  // data/ppData2017_highForest_pfJets_20EventsMixed_xjBins_finalTrackCorr_JECv4_wtaAxis_allCorrections_processed_2019-09-28.root
  // data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_processed_2019-08-13.root
  // data/dijet_pp_highForest_pfJets_noUncOrInc_allCorrections_wtaAxis_processed_2019-07-13.root
  // data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_allCorrections_processed_2019-09-28.root
  // data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_testSeagull_allCorrections_processed_2019-09-28.root
  // data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_allCorrections_processed_2019-09-28.root
  TFile *pbpbFile = TFile::Open("data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsWithShifterCentrality_trackDeltaRonlyInLowPt_processed_2019-10-17.root");
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsWithShifterCentrality_trackDeltaRonlyInLowPt_processed_2019-10-17.root
  // data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_5eveMix_smootheMixing_allCorrections_processed_2019-10-16.root
  // data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_allCorrections_processed_2019-10-07.root
  // data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncorr_5eveMix_widerMixingPeakFinding_allCorrections_processed_2019-10-07.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_JECv6_finalTrack_eschemeAxis_allCorrections_allPtTrackDeltaR_processed_2019-10-16.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsWithCentShift_trackDeltaRonlyLowPt_processed_2019-10-16.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_JECv6_finalTrack_eschemeAxis_allCorrections_trackDeltaRFromWta_processed_2019-10-02.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trig_JECv6_finalTrack_eschemeAxis_allCorrections_shortRangeJff_processed_2019-10-16.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncOrInc_improvisedMixingFromSubeNon0AtLowerPt_wtaAxis_JECv6_noSymmetrySpilloverLooseCutUntil8_processed_2019-09-26.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0AtLowPt_wtaAxis_JECv6_newSpilloverLooseCutUntil8_processed_2019-09-26.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0AtLowPt_wtaAxis_JECv6_spilloverVeryLooseCutUntil8_processed_2019-09-26.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_wtaAxis_selectiveSeagull_modifiedSpillover_JECv6_processed_2019-09-26.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_allCorrections_noStatErrorFromCorrections_lowPtResidualTrack_processed_2019-10-07_fiveJobsMissing.root
  // data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_allCorrections_processed_2019-10-07.root
  
  TFile *ppUncertaintyFile = TFile::Open("uncertainties/systematicUncertaintyForPp_20percentSpillJff_2019-09-30.root");
  // uncertainties/systematicUncertaintyForPp_20percentSpillJff_2019-09-30.root
  // uncertainties/systematicUncertaintyForPythia8RecoGen_mcMode_2019-10-05.root
  TFile *pbpbUncertaintyFile = TFile::Open("uncertainties/systematicUncertaintyForPbPb_15percentSpill5Jff_2019-10-14.root");
  // uncertainties/systematicUncertaintyForPythiaHydjetRecoGen_mcMode_2019-10-06.root
  // uncertainties/systematicUncertaintyForPythiaHydjetRecoGen_mcMode_2019-10-05.root
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_oldJES_15percentSpill10Jff_2019-10-17.root
  // uncertainties/systematicUncertaintyForPbPb_15percentSpill5Jff_2019-10-14.root
  // uncertainties/systematicUncertaintyForPp_15percentSpill20Jff_2019-10-01.root
  
  // Create histogram managers for pp and PbPb
  DijetHistogramManager *ppHistograms = new DijetHistogramManager(ppFile);
  DijetHistogramManager *pbpbHistograms = new DijetHistogramManager(pbpbFile);
  JffCorrector *ppUncertaintyProvider = new JffCorrector();
  JffCorrector *pbpbUncertaintyProvider = new JffCorrector();
  
  // Choose which figure sets to draw
  bool drawJetShapePptoPbPbRatio = true;
  bool drawJetShapeSymmetricAsymmetricRatio = false;
  
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
  int jetTrackIndex[2] = {DijetHistogramManager::kPtWeightedTrackLeadingJet, DijetHistogramManager::kPtWeightedTrackSubleadingJet};
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
  plotJetShapeBigAsymmetry(ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kPtWeightedTrackLeadingJet);
  plotJetShapeBigAsymmetry(ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kPtWeightedTrackSubleadingJet);
}
