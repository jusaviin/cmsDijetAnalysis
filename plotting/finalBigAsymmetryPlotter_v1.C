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
    for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
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
    if(monteCarloLabels){
    }
  }
  
  // ==============================
  // **  Draw the distributions  **
  // ==============================
  
  // Draw all the distributions to big canvases
  auxi_canvas *bigCanvas = new auxi_canvas(Form("theReallyBigCanvas%d",iJetTrack), "", 2500, 2500);
  bigCanvas->SetMargin(0.06, 0.01, 0.08, 0.2); // Margin order: Left, Right, Bottom, Top
  bigCanvas->divide(4,5);
  
  
  TBox *box = new TBox();
  TLatex *mainTitle = new TLatex();
  
  TLine *line[nAsymmetryBins];
  const char *pbpbLabel = "PbPb";
  const char *ppLabel = "pp";
  
  if(monteCarloLabels){
    pbpbLabel = "Pythia+Hydjet";
    ppLabel = "Pythia8";
  }
  
  TString xjbin[] = {"0.0 < x_{j} < 0.6", "0.6 < x_{j} < 0.8", "0.8 < x_{j} < 1.0", "all dijets"};
  
  for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins+1; iAsymmetry++){
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      
      //cout << "In the loop iAsymmetry: " << iAsymmetry << " iCentrality: " << iCentrality << endl;
      
      if( iAsymmetry == 3) {
        bigCanvas->CD(5-iCentrality);
      }else
        bigCanvas->CD(10+iAsymmetry*5-iCentrality);
      gPad->SetLogy();
      jetShapeStack[iCentrality][iAsymmetry]->drawStack();
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetNdivisions(505);
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetLabelSize(0.09);
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitleOffset(1.2);
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitleSize(0.1);
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitleOffset(0.5);
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitle("#Deltar");
      if(iCentrality == nCentralityBins && iAsymmetry == 2){
        //most left-bottom conner
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitleSize(0.11);
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitleOffset(0.83);
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetLabelSize(0.08);
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetLabelOffset(0.03);
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetLabelSize(0.07);
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetLabelOffset(0.01);
      }else {
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitleOffset(0.8);
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitleSize(0.15);
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitleSize(0.14);
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitleOffset(0.71);
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetLabelSize(0.11);
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetLabelOffset(0.011);
      }
      if(normalizeJetShape){
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitle("#rho(#Deltar)");
      } else {
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitle("#Rho(#Deltar)");
      }
      jetShapeStack[iCentrality][iAsymmetry]->hst->Draw();
      
      if(drawUncertainties) sumUncertainty[iCentrality][iAsymmetry]->Draw("same e2");
      if(iCentrality == nCentralityBins){
        mainTitle->SetTextFont(22);
        if( iAsymmetry == 2)
          mainTitle->SetTextSize(0.08);
        else
          mainTitle->SetTextSize(0.11);
        mainTitle->DrawLatexNDC(0.35, 0.88, ppLabel);
        mainTitle->DrawLatexNDC(0.55, 0.86, xjbin[iAsymmetry]);
      }
      else {
        mainTitle->SetTextFont(22);
        mainTitle->SetTextSize(0.11);
        mainTitle->DrawLatexNDC(0.19, 0.9, pbpbLabel);
        mainTitle->DrawLatexNDC(0.19, 0.80, cent_lab[iCentrality]);
      }
      
      
    } // Centrality loop
    
  } // Asymmetry loop for drawing distributions to the really big canvas
  
  TLegend* ptLegend1 = new TLegend(0.04, 0.82, 0.27, 0.89);
  TLegend* ptLegend2 = new TLegend(0.28, 0.82, 0.51, 0.89);
  TLegend* ptLegend3 = new TLegend(0.53 ,0.82, 0.76, 0.89);
  TLegend* ptLegend4 = new TLegend(0.77, 0.82, 1.00, 0.89);
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
  /*
   bigCanvas->CD(6);
   mainTitle->SetTextFont(42);
   mainTitle->SetTextSize(0.13);
   mainTitle->DrawLatexNDC(0.08,0.44,"0.0 < x_{j} < 0.6");
   
   bigCanvas->CD(16);
   mainTitle->DrawLatexNDC(0.08,0.44,"0.6 < x_{j} < 0.8");
   
   bigCanvas->CD(26);
   mainTitle->SetTextSize(0.12);
   mainTitle->DrawLatexNDC(0.09,0.55,"0.8 < x_{j} < 1.0");
   
   */
  
  bigCanvas->cd(0);
  
  ptLegend1->Draw();
  ptLegend2->Draw();
  ptLegend3->Draw();
  ptLegend4->Draw();
  
  mainTitle->SetTextFont(62);
  mainTitle->SetTextSize(cmsSize);
  mainTitle->DrawLatexNDC(cmsPosition[iJetTrack/3]+0.14, cmsPositionY, "CMS");
  
  mainTitle->SetTextFont(42);
  mainTitle->SetTextSize(0.035);
  mainTitle->DrawLatexNDC(jetShapeTitlePosition[iJetTrack/3]-0.08, 0.97, jetShapeTitle[iJetTrack/3]);
  
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
  mainTitle->SetTextSize(0.02);
  box->DrawBox(0.24,0.017, 0.253, 0.068);
  box->DrawBox(0.42,0.017, 0.45, 0.068);
  box->DrawBox(0.605,0.017,0.63, 0.068);
  box->DrawBox(0.78,0.017, 0.81, 0.068);
  box->DrawBox(0.98,0.017, 0.99, 0.068);
  
  box->DrawBox(0.23,0.3067, 0.243, 0.315);
  box->DrawBox(0.23,0.5735, 0.243, 0.58);
  
  //zeroMark
  mainTitle->DrawLatex(0.242, 0.0545, "0");
  mainTitle->DrawLatex(0.428, 0.0545, "0");
  mainTitle->DrawLatex(0.614, 0.0545, "0");
  mainTitle->DrawLatex(0.801, 0.0545, "0");
  mainTitle->DrawLatex(0.982, 0.0545, "1");
  
  //=======================
  // canvas2
  auto *ratioCanvas = new TCanvas(Form("ratioCanvas_%d",iJetTrack), "", 1600, 1400);
  ratioCanvas->SetMargin(0.2, 0.05, 0.25, 0.25); // Margin order: Left, Right, Bottom, Top
  ratioCanvas->Divide(4,3,0 ,0);
  TH1D* auxi_hist[4];
  for(int i=0 ; i< 4; i++){
    auxi_hist[i] = (TH1D*) ratioUncertainty[i][3]->Clone(Form("auxi_hist_%d", i));
    auxi_hist[i]->SetFillColorAlpha(kWhite, 0 );
    auxi_hist[i]->SetMarkerColor(kWhite);
    auxi_hist[i]->SetMarkerStyle(20);
    auxi_hist[i]->SetMarkerSize(1.1);
  }
  for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      ratioCanvas->cd(4+iAsymmetry*4-iCentrality);
      if(iCentrality==3) {
        ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->SetTitleOffset(0.7);
        ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->SetTitleSize(0.125);
        ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->SetLabelOffset(0.01);
        ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->SetLabelSize(0.08);
      }else {
        ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->SetTitleOffset(0.58);
        ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->SetTitleSize(0.15);
        ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->SetLabelOffset(-0.003);
        ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->SetLabelSize(0.1);
      }
      if(iAsymmetry==2) {
        if(iCentrality == 3)
          ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetTitleOffset(1);
        ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetTitleSize(0.08);
        ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetLabelOffset(0.008);
        ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetLabelSize(0.07);
      }else {
        ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetTitleSize(0.1);
        ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetLabelOffset(0.01);
        ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetLabelSize(0.08);
      }
      
      ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetTitle("#rho(#Deltar)_{PbPb}/#rho(#Deltar)_{pp}");
      ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->CenterTitle();
      ratioUncertainty[iCentrality][iAsymmetry]->SetStats(0);
      ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetNdivisions(505);
      ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->CenterTitle();
      ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->SetTitle("#Deltar");
      ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->CenterTitle();
      ratioUncertainty[iCentrality][iAsymmetry]->SetAxisRange(0.01, 3.7, "Y");
      ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->SetNdivisions(505);
      ratioUncertainty[iCentrality][iAsymmetry]->SetAxisRange(0, 0.99, "X");
      
      
      ratioUncertainty[iCentrality][iAsymmetry]->SetMarkerStyle(20);
      ratioUncertainty[iCentrality][iAsymmetry]->SetMarkerSize(1.8);
      ratioHistogram[iCentrality][iAsymmetry]->SetMarkerStyle(20);
      ratioHistogram[iCentrality][iAsymmetry]->SetMarkerSize(1.8);
      ratioHistogram[iCentrality][iAsymmetry]->SetMarkerSize(1.8);
      if(iAsymmetry==0){
        ratioUncertainty[iCentrality][iAsymmetry]->SetFillColorAlpha(kRed+2, 0.4);
        ratioUncertainty[iCentrality][iAsymmetry]->SetLineColor(kRed+2);
        ratioHistogram[iCentrality][iAsymmetry]->SetLineColor(kRed+2);
        ratioUncertainty[iCentrality][iAsymmetry]->SetMarkerColor(kRed+2);
        ratioHistogram[iCentrality][iAsymmetry]->SetMarkerColor(kRed+2);
      }else if(iAsymmetry==1){
        Color_t color = kViolet-5;
        ratioUncertainty[iCentrality][iAsymmetry]->SetFillColorAlpha(color, 0.4);
        ratioUncertainty[iCentrality][iAsymmetry]->SetLineColor(color);
        ratioUncertainty[iCentrality][iAsymmetry]->SetMarkerColor(color);
        ratioHistogram[iCentrality][iAsymmetry]->SetLineColor(color);
        ratioHistogram[iCentrality][iAsymmetry]->SetMarkerColor(color-1);
      }else if(iAsymmetry==2){
        ratioUncertainty[iCentrality][iAsymmetry]->SetFillColorAlpha(kAzure+2, 0.4);
        ratioUncertainty[iCentrality][iAsymmetry]->SetLineColor(kAzure+2);
        ratioUncertainty[iCentrality][iAsymmetry]->SetMarkerColor(kAzure+2);
        ratioHistogram[iCentrality][iAsymmetry]->SetMarkerColor(kAzure+2);
        ratioHistogram[iCentrality][iAsymmetry]->SetLineColor(kAzure+2);
      }
      ratioHistogram[iCentrality][nAsymmetryBins]->SetLineColor(kBlack);
      ratioHistogram[iCentrality][nAsymmetryBins]->SetMarkerSize(1.3);
      ratioHistogram[iCentrality][nAsymmetryBins]->SetMarkerColor(kBlack);
      ratioHistogram[iCentrality][nAsymmetryBins]->SetMarkerStyle(24);
      ratioHistogram[iCentrality][nAsymmetryBins]->SetFillColorAlpha(kWhite, 0);
      ratioUncertainty[iCentrality][nAsymmetryBins]->SetLineColor(kBlack);
      ratioUncertainty[iCentrality][nAsymmetryBins]->SetFillColorAlpha(kGray+2, 0.3);
      ratioUncertainty[iCentrality][nAsymmetryBins]->SetMarkerColor(kWhite);
      ratioUncertainty[iCentrality][nAsymmetryBins]->SetMarkerStyle(20);
      ratioUncertainty[iCentrality][nAsymmetryBins]->SetMarkerSize(1.1);
      
      if(iAsymmetry== 0 && iCentrality == 0 ){
      }
      ratioUncertainty[iCentrality][iAsymmetry]->Draw("same e2");
      ratioUncertainty[iCentrality][nAsymmetryBins]->Draw("same e2");
      ratioHistogram[iCentrality][iAsymmetry]->Draw("same");
      ratioHistogram[iCentrality][nAsymmetryBins]->Draw("same");
      auxi_hist[iCentrality]->Draw("same e2");
      
      line[iAsymmetry] = new TLine();
      line[iAsymmetry]->SetLineStyle(2);
      line[iAsymmetry]->DrawLine(0, 1, 1, 1);
      
    }
  }
  //canvas2caption
  ratioCanvas->cd(0);
  mainTitle->SetTextFont(62);
  mainTitle->SetTextSize(0.04);
  mainTitle->DrawLatexNDC(0.06, 0.95, "CMS");
  
  mainTitle->SetTextFont(42);
  mainTitle->SetTextSize(0.035);
  mainTitle->DrawLatexNDC(.15, 0.95, Form("%s ratio",jetShapeTitle[iJetTrack/3]));
  mainTitle->DrawLatexNDC(.52, 0.95, "5.02 TeV   pp 320 pb^{-1}   PbPb 1.7 nb^{-1}");
  
  auto ltc2 = new TLegend(0.27, 0.45, 0.72, 0.8);
  ltc2->SetLineColor(0);
  ltc2->AddEntry(ratioUncertainty[0][3], xjbin[3], "lpf");
  for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
    ltc2->AddEntry(ratioUncertainty[0][iAsymmetry], xjbin[iAsymmetry], "lpf");
    ratioCanvas->cd(4+iAsymmetry*4-3);
  }
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    ratioCanvas->cd(4-iCentrality);
    mainTitle->SetTextFont(22);
    mainTitle->SetTextSize(0.09);
    if(iCentrality == 3){
      mainTitle->DrawLatexNDC(0.3, 0.88, "Cent: "+cent_lab[iCentrality]);
    }else {
      mainTitle->DrawLatexNDC(0.1, 0.88, "Cent: "+cent_lab[iCentrality]);
    }
  }
  
  ratioCanvas->cd(1);
  ltc2->Draw();
  
  //bigCanvas->SaveAs("js_dr_normal_new.eps");
  bigCanvas->SaveAs(Form("figures/final%s_%s_DesDng.pdf", saveString, jetShapeSaveName[iJetTrack/3]));
  ratioCanvas->SaveAs(Form("figures/final%s_%s_ratioDesDng.pdf", saveString, jetShapeSaveName[iJetTrack/3]));
  //bigCanvas->SaveAs("js_dr_normal_v3.eps");
  //bigCanvas->SaveAs("js_dr_normal_v3.pdf");
  
}

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void finalBigAsymmetryPlotter_v1(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Open data files for pp and PbPb data
  TFile *ppFile = TFile::Open("data/ppData2017_highForest_pfJets_20EventsMixed_finalTrackCorr_xjBins_JECv4_wtaAxis_tunedSeagull_allCorrections_processed_2019-10-17.root");
  TFile *pbpbFile = TFile::Open("data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsWithShifterCentrality_trackDeltaRonlyInLowPt_processed_2019-10-17.root");
  
  TFile *ppUncertaintyFile = TFile::Open("uncertainties/systematicUncertaintyForPp_20percentSpillJff_2019-09-30.root");
  TFile *pbpbUncertaintyFile = TFile::Open("uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_includeTrackDeltaR_2020-01-27.root");
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
