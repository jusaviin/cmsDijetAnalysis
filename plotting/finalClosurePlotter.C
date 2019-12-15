#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "stackHistXiao.h"
#include "xCanvas.h"
#include "JDrawer.h"

void plotJetShapeXiao(DijetHistogramManager *ppHistograms, DijetHistogramManager *pbpbHistograms, JffCorrector *ppUncertaintyProvider, JffCorrector *pbpbUncertaintyProvider, int iJetTrack){

  const int nTrackPtBins = pbpbHistograms->GetNTrackPtBins();
  const int nCentralityBins = pbpbHistograms->GetNCentralityBins();
  const int nAsymmetryBins = pbpbHistograms->GetNAsymmetryBins();
  
  const char* jetShapeTitle[] = {"Leading jet closure","Subleading jet closure","Jet shape"};
  const char* jetShapeSaveName[] = {"trackLeadingJet","trackSubleadingJet","trackInclusiveJet"};
  const char* xjString[] = {"0.0 < x_{j} < 0.6","0.6 < x_{j} < 0.8","0.8 < x_{j} < 1.0","x_{j} inclusive"};
  const char* asymmetrySaveName[] = {"_A=0v0-0v6","_A=0v6-0v8","_A=0v8-1v0",""};
  const char* saveString = {"JetShape"};
  
  bool drawInclusive = false;
  bool drawInclusiveRatio = false;
  bool drawUncertainties = true;
  bool monteCarloLabels = true;
  int firstAsymmetryBin = nAsymmetryBins;  // Set here 0 to process asymmetry bins and nAsymmetrybins to disable them
  bool normalizeJetShape = true;           // True: draw rho. False: draw P.
  
  // Draw the comparison only for leading jet shapes
  if(!(iJetTrack == DijetHistogramManager::kPtWeightedTrackLeadingJet || iJetTrack == DijetHistogramManager::kPtWeightedTrackInclusiveJet) || !normalizeJetShape){
    drawInclusive = false;
    drawInclusiveRatio = false;
  }
  
  // Change the titles if the jet shape is not normalized to one
  if(!normalizeJetShape){
    jetShapeTitle[0] = "Leading jet radial momentum distribution";
    jetShapeTitle[1] = "Subleading jet radial momentum distribution";
    jetShapeTitle[2] = "Jet radial momentum distribution";
    saveString = "JetRadialMomentum";
  }
  
  // Temporary: Get the ratio between pp and PbPb jet shape
  TH1D *sumHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *sumHistogramPp[nCentralityBins][nAsymmetryBins+1];
  TH1D *jetShapeArray[nCentralityBins+1][nTrackPtBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPp[nAsymmetryBins+1];
  TH1D *sumUncertainty[nCentralityBins+1][nAsymmetryBins+1];
  TH1D *helperHistogram;
  double shapeIntegralPbPb[nCentralityBins][nAsymmetryBins+1];
  double shapeIntegralPp[nCentralityBins][nAsymmetryBins+1];
  
  
  // Read the pT summed jet shape histograms
  for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      sumHistogramPbPb[iCentrality][iAsymmetry] = (TH1D*) pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, 0)->Clone(Form("sumPbPb%d",iCentrality));
      jetShapeArray[iCentrality][0][iAsymmetry] = pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, 0);
      
      sumHistogramPp[iCentrality][iAsymmetry] = (TH1D*) ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, 0)->Clone(Form("sumPp%d",iCentrality));
      
    }
    jetShapeArray[nCentralityBins][0][iAsymmetry] = ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, 0);
  }

  // Sum the pT:s
  for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iTrackPt = 1; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        helperHistogram = (TH1D*)pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,iCentrality,iTrackPt)->Clone();
        sumHistogramPbPb[iCentrality][iAsymmetry]->Add(helperHistogram);
        jetShapeArray[iCentrality][iTrackPt][iAsymmetry] = pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,iCentrality,iTrackPt);
      
        helperHistogram = (TH1D*)ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape,iJetTrack, iAsymmetry,iCentrality,iTrackPt)->Clone();
        sumHistogramPp[iCentrality][iAsymmetry]->Add(helperHistogram);
      }
      jetShapeArray[nCentralityBins][iTrackPt][iAsymmetry] = ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, iTrackPt);
    } // track pT
  }
  
  // Normalize the jet shape histograms
  for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
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
      
    
    
    shapeIntegralPp[iCentrality][iAsymmetry] = sumHistogramPp[iCentrality][iAsymmetry]->Integral(1,sumHistogramPp[iCentrality][iAsymmetry]->FindBin(0.99),"width");
    if(normalizeJetShape) sumHistogramPp[iCentrality][iAsymmetry]->Scale(1.0/shapeIntegralPp[iCentrality][iAsymmetry]);
    
      }
      
    // To draw the systematic uncertainties, clone the sum histograms to uncertainty histograms
    sumUncertainty[nCentralityBins][iAsymmetry] = (TH1D*) sumHistogramPp[0][iAsymmetry]->Clone(Form("sumUncertaintyPp%d", iAsymmetry));
    
    if(normalizeJetShape){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        jetShapeArray[nCentralityBins][iTrackPt][iAsymmetry]->Scale(1.0/shapeIntegralPp[0][iAsymmetry]);
      }
    }
    
    // Read the uncertainties for the pT summed jet shapes and scale them for jet shapes
    uncertaintyHistogramPp[iAsymmetry] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry);
    if(normalizeJetShape) uncertaintyHistogramPp[iAsymmetry]->Scale(1.0/shapeIntegralPp[0][iAsymmetry]);
    
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
  TH1D* asymmetryRatioUncertainty[nCentralityBins+1][nAsymmetryBins];
  double ppValue, pbpbValue, ratioValue;
  double ppUncertainty, pbpbUncertainty, ratioUncertaintyValue;
  double smoothMe = 1;
  
  for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      //ratioUncertainty[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("ratioUncertainty%d%d", iCentrality, iAsymmetry));
      ratioHistogram[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("ratio%d%d", iCentrality, iAsymmetry));
      ratioHistogram[iCentrality][iAsymmetry]->Divide(sumHistogramPp[iCentrality][iAsymmetry]);
      
      ratioUncertainty[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("ratioUncertainty%d%d", iCentrality, iAsymmetry));
      
      if(iJetTrack < 3){
        if(iCentrality == 0) smoothMe = 0.55;
        if(iCentrality == 1) smoothMe = 0.9;
        if(iCentrality == 2) smoothMe = 1.7;
        if(iCentrality == 3) smoothMe = 1;
      } else {
        if(iCentrality == 0) smoothMe = 1;
        if(iCentrality == 1) smoothMe = 0.6;
        if(iCentrality == 2) smoothMe = 3;
        if(iCentrality == 3) smoothMe = 1;
      }
      
      // Calculate the systemtic uncertainty for the ratio
      for(int iBin = 1; iBin < ratioUncertainty[iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
        
        // Manual calculation of the systematic uncertainty
        ppValue = sumHistogramPp[iCentrality][iAsymmetry]->GetBinContent(iBin);
        pbpbValue = sumHistogramPbPb[iCentrality][iAsymmetry]->GetBinContent(iBin);
        ppUncertainty = uncertaintyHistogramPp[iAsymmetry]->GetBinContent(iBin);
        pbpbUncertainty = uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->GetBinContent(iBin);
        ratioValue = ratioHistogram[iCentrality][iAsymmetry]->GetBinContent(iBin);
        ratioUncertaintyValue = TMath::Sqrt(TMath::Power(pbpbUncertainty/ppValue,2)+TMath::Power(pbpbValue*ppUncertainty/TMath::Power(ppValue,2),2));
        ratioUncertainty[iCentrality][iAsymmetry]->SetBinContent(iBin,ratioValue);
        ratioUncertainty[iCentrality][iAsymmetry]->SetBinError(iBin,ratioUncertaintyValue*smoothMe);
      }
      
    }
  }

  TString cent_lab[4] = {"0-10%", "10-30%", "30-50%", "50-90%"};
  stackHist *jetShapeStack[nCentralityBins+1][nAsymmetryBins+1];
  
  // Create the stacked jet shape histograms and set a good drawing style for the uncertainty histgorams
  for(int iCentrality = 0; iCentrality < nCentralityBins+1; iCentrality++){
    for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      
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
    for(int iAsymmetry = firstAsymmetryBin; iAsymmetry < nAsymmetryBins+1; iAsymmetry++){
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
  
  // Draw all the distributions to big canvases
  auxi_canvas *bigCanvas[nAsymmetryBins+1];
  TBox *box[nAsymmetryBins+1];
  TLatex *mainTitle[nAsymmetryBins+1];
  TLine *line[nAsymmetryBins+1];
  const char *pbpbLabel = "PbPb";
  const char *ppLabel = "pp reference";
  
  if(monteCarloLabels){
    pbpbLabel = "Pythia+Hydjet";
    ppLabel = "Pythia8";
  }
  
  for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    
    bigCanvas[iAsymmetry] = new auxi_canvas(Form("bigCanvas%d%d",iAsymmetry,iJetTrack), "", 2500, 2000);
    bigCanvas[iAsymmetry]->SetMargin(0.06, 0.01, 0.08, 0.02);
    bigCanvas[iAsymmetry]->divide(3,4);
    
    mainTitle[iAsymmetry] = new TLatex();
    
    for(int iCentrality=0; iCentrality < nCentralityBins; iCentrality++){
      bigCanvas[iAsymmetry]->CD(8-iCentrality);
      gPad->SetLogy();
      jetShapeStack[iCentrality][iAsymmetry]->drawStack();
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetNdivisions(505);
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetLabelSize(0.08);
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitleOffset(0.9);
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitleSize(0.1);
      if(normalizeJetShape){
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitle("#rho(#Deltar)");
      } else {
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitle("#Rho(#Deltar)");
      }
      jetShapeStack[iCentrality][iAsymmetry]->hst->Draw();

      if(drawUncertainties) sumUncertainty[iCentrality][iAsymmetry]->Draw("same e2");
      if(iCentrality==3 ){
        mainTitle[iAsymmetry]->SetTextFont(22);
        mainTitle[iAsymmetry]->SetTextSize(.085);
        mainTitle[iAsymmetry]->DrawLatexNDC(.35, 0.9, pbpbLabel);
        mainTitle[iAsymmetry]->DrawLatexNDC(.35, .82, cent_lab[iCentrality]);
      }
      else {
        mainTitle[iAsymmetry]->SetTextFont(22);
        mainTitle[iAsymmetry]->SetTextSize(.09);
        mainTitle[iAsymmetry]->DrawLatexNDC(.19, 0.9, pbpbLabel);
        mainTitle[iAsymmetry]->DrawLatexNDC(.19, .82, cent_lab[iCentrality]);
      }
      bigCanvas[iAsymmetry]->CD(12-iCentrality);
      //ratio[iCentrality]->GetYaxis()->SetNdivisions(505);
      ratioHistogram[iCentrality][iAsymmetry]->SetTitle("");
      ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetTitleOffset(1.1);
      ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetTitle("#Deltar");
      if(monteCarloLabels){
        ratioHistogram[iCentrality][iAsymmetry]->SetAxisRange(0.5, 1.5, "Y");
      } else if(normalizeJetShape){
        ratioHistogram[iCentrality][iAsymmetry]->SetAxisRange(0, 3.2, "Y");
      } else {
        ratioHistogram[iCentrality][iAsymmetry]->SetAxisRange(0, 4.2, "Y");
      }
      ratioHistogram[iCentrality][iAsymmetry]->SetAxisRange(0, 0.99, "X");
      if( iCentrality<3 )  {
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetTitleOffset(0.7);
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetTitleSize(0.11);
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetNdivisions(505);
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetLabelSize(0.08);
        ratioHistogram[iCentrality][iAsymmetry]->GetYaxis()->SetNdivisions(505);
      }
      if(iCentrality==3 ){
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetTitleOffset(0.92);
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetTitleSize(0.086);
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetNdivisions(505);
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetLabelOffset(0.016);
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetLabelSize(0.064);
        ratioHistogram[iCentrality][iAsymmetry]->GetYaxis()->SetNdivisions(505);
        ratioHistogram[iCentrality][iAsymmetry]->GetYaxis()->SetLabelOffset(0.02);
        ratioHistogram[iCentrality][iAsymmetry]->GetYaxis()->SetLabelSize(0.07);
        ratioHistogram[iCentrality][iAsymmetry]->GetYaxis()->SetTitleOffset(0.9);
        ratioHistogram[iCentrality][iAsymmetry]->GetYaxis()->SetTitleSize(0.08);
        if(normalizeJetShape){
          ratioHistogram[iCentrality][iAsymmetry]->GetYaxis()->SetTitle("#rho(#Deltar)_{PbPb}/#rho(#Deltar)_{pp}");
        } else {
          ratioHistogram[iCentrality][iAsymmetry]->GetYaxis()->SetTitle("#Rho(#Deltar)_{PbPb}/#Rho(#Deltar)_{pp}");
        }
        mainTitle[iAsymmetry]->SetTextSize(0.073);
        //mainTitle->DrawLatexNDC(.25, .92, "PbPb - pp");
        //mainTitle->DrawLatexNDC(.25, .84, cent_lab[iCentrality]);
      }
      else{
        //mainTitle->SetTextSize(0.09);
        //mainTitle->DrawLatexNDC(.05, .92, "PbPb - pp");
        //mainTitle->DrawLatexNDC(.05, .84, cent_lab[iCentrality]);
      }
      
      ratioHistogram[iCentrality][iAsymmetry]->GetYaxis()->CenterTitle();
      ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->CenterTitle();
      ratioHistogram[iCentrality][iAsymmetry]->SetStats(0);
      ratioHistogram[iCentrality][iAsymmetry]->Draw("same");
      
      line[iAsymmetry] = new TLine();
      line[iAsymmetry]->SetLineStyle(2);
      line[iAsymmetry]->DrawLine(0, 1, 1, 1);
      
      if(drawUncertainties) ratioUncertainty[iCentrality][iAsymmetry]->Draw("same e2");
      
    }
    bigCanvas[iAsymmetry]->CD(1);
    gPad->SetLogy();
    jetShapeStack[4][iAsymmetry]->drawStack();
    jetShapeStack[4][iAsymmetry]->hst->GetYaxis()->SetNdivisions(505);
    jetShapeStack[4][iAsymmetry]->hst->GetYaxis()->SetLabelSize(0.08);
    jetShapeStack[4][iAsymmetry]->hst->GetYaxis()->SetTitleOffset(0.9);
    jetShapeStack[4][iAsymmetry]->hst->GetYaxis()->SetTitleSize(0.1);
    if(normalizeJetShape){
      jetShapeStack[4][iAsymmetry]->hst->GetYaxis()->SetTitle("#rho(#Deltar)");
    } else {
      jetShapeStack[4][iAsymmetry]->hst->GetYaxis()->SetTitle("#Rho(#Deltar)");
    }

    mainTitle[iAsymmetry]->SetTextSize(0.085);
    mainTitle[iAsymmetry]->DrawLatexNDC(0.35, 0.88, ppLabel);
    if(drawUncertainties) sumUncertainty[4][iAsymmetry]->Draw("same e2");
    
    TLegend* lt1 = new TLegend(0.01, 0.1, 1, 0.5);
    TLegend* lt2 = new TLegend(0.0, 0.1, 1, 0.5);
    TLegend* lt3 = new TLegend(0.0 ,0.22, 1, 0.5);
    TLegend* lt4 = new TLegend(0.25, 0.85, 0.85, 0.95);
    TLegend* lt5 = new TLegend(0.07, 0.85, 0.67, 0.95);
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
    
    
    lt1->AddEntry(jetShapeStack[4][iAsymmetry]->hist_trunk.at(0), "0.7 < p_{T}^{trk}< 1 GeV","f");
    lt1->AddEntry(jetShapeStack[4][iAsymmetry]->hist_trunk.at(1), "1 < p_{T}^{trk}< 2 GeV","f");
    lt1->AddEntry(jetShapeStack[4][iAsymmetry]->hist_trunk.at(2), "2 < p_{T}^{trk}< 3 GeV","f");
    
    lt2->AddEntry(jetShapeStack[4][iAsymmetry]->hist_trunk.at(3), "3 < p_{T}^{trk}< 4 GeV","f");
    lt2->AddEntry(jetShapeStack[4][iAsymmetry]->hist_trunk.at(4), "4 < p_{T}^{trk}< 8 GeV","f");
    lt2->AddEntry(jetShapeStack[4][iAsymmetry]->hist_trunk.at(5), "8 < p_{T}^{trk}< 12 GeV","f");
    
    
    lt3->AddEntry(jetShapeStack[4][iAsymmetry]->hist_trunk.at(6), "12 < p_{T}^{trk}< 300 GeV","f");
    lt3->AddEntry(sumUncertainty[4][iAsymmetry], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
    //lt3->AddEntry(jetShapeStack[4]->hist_trunk.at(8), "20 < p_{T}^{trk}< 300 GeV","f");
    
    lt4->AddEntry(ratioUncertainty[0][iAsymmetry], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
    bigCanvas[iAsymmetry]->CD(9);
    lt4->Draw();
    
    
    //bigCanvas->CD(1);
    //TLegend* lt6 = new TLegend(0.3,0.7,0.98,0.85);
    //lt6->SetTextSize(0.07);
    //lt6->SetLineColor(kWhite);
    //lt6->SetFillColor(kWhite);
    //lt6->AddEntry(js_dr_err_all[0], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
    //lt6->Draw();
    
    bigCanvas[iAsymmetry]->CD(2);
    line[iAsymmetry]->SetLineStyle(1);
    line[iAsymmetry]->DrawLineNDC(0, 0, 0, 1);
    lt1->Draw();
    bigCanvas[iAsymmetry]->CD(3);
    lt2->Draw();
    bigCanvas[iAsymmetry]->CD(4);
    lt3->Draw();
    
    double cmsPosition[] = {0.42,0.4,0.42};
    double cmsPositionY = 0.94;
    double jetShapeTitlePosition[] = {0.52,0.5,0.52};
    double xjPosition[] = {0.76,0.77,0.76};
    double overallShift[] = {-0.02,-0.02,-0.02,0};
    double systemPosition = 0.41;
    double selectionPosition = 0.375;
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
    
    bigCanvas[iAsymmetry]->cd(0);
    mainTitle[iAsymmetry]->SetTextFont(62);
    mainTitle[iAsymmetry]->SetTextSize(cmsSize);
    mainTitle[iAsymmetry]->DrawLatexNDC(cmsPosition[iJetTrack/3]+overallShift[iAsymmetry], cmsPositionY, "CMS");
    
    mainTitle[iAsymmetry]->SetTextFont(42);
    mainTitle[iAsymmetry]->SetTextSize(0.035);
    mainTitle[iAsymmetry]->DrawLatexNDC(jetShapeTitlePosition[iJetTrack/3]+overallShift[iAsymmetry], 0.94, jetShapeTitle[iJetTrack/3]);
    
    mainTitle[iAsymmetry]->SetTextFont(42);
    mainTitle[iAsymmetry]->SetTextSize(0.035);
    mainTitle[iAsymmetry]->DrawLatexNDC(xjPosition[iJetTrack/3]+overallShift[iAsymmetry], 0.94, xjString[iAsymmetry]);
    
    mainTitle[iAsymmetry]->SetTextSize(0.03);
    mainTitle[iAsymmetry]->DrawLatexNDC(systemPosition, 0.9, "pp 320 pb^{-1} (5.02 TeV)  PbPb 1.7 nb^{-1} (5.02 TeV)");
    mainTitle[iAsymmetry]->SetTextSize(0.025);
    mainTitle[iAsymmetry]->DrawLatexNDC(selectionPosition, 0.86, "anti-k_{T} R = 0.4, |#eta_{jet}| < 1.6, p_{T,1} > 120 GeV, p_{T,2} > 50 GeV, #Delta#phi_{1,2} > #frac{5#pi}{6}");
    //  lb->drawText("(p_{T}> 120 GeV, |#eta_{jet}|<1.6)", 0.2, 0.25, 4);
    
    box[iAsymmetry] = new TBox();
    box[iAsymmetry]->SetFillColor(kWhite);
    bigCanvas[iAsymmetry]->cd(0);
    box[iAsymmetry]->DrawBox(0.285,.047, 0.3, 0.072);
    mainTitle[iAsymmetry]->SetTextSize(.025);
    mainTitle[iAsymmetry]->DrawLatex(0.29, 0.055, "0");
    box[iAsymmetry]->DrawBox(0.518,.047, 0.533, 0.072);
    box[iAsymmetry]->DrawBox(0.75,.047, 0.765, 0.072);
    
    mainTitle[iAsymmetry]->DrawLatex(0.523, 0.055, "0");
    mainTitle[iAsymmetry]->DrawLatex(0.755, 0.055, "0");
    mainTitle[iAsymmetry]->DrawLatex(0.985, 0.055, "1");
    
    //bigCanvas->SaveAs("js_dr_normal_new.eps");
    bigCanvas[iAsymmetry]->SaveAs(Form("figures/final%s_%s%s_finalClosureCheck.pdf", saveString, jetShapeSaveName[iJetTrack/3], asymmetrySaveName[iAsymmetry]));
    //bigCanvas->SaveAs("js_dr_normal_v3.eps");
    //bigCanvas->SaveAs("js_dr_normal_v3.pdf");
    
  } // Asymmetry loop for drawing the big canvas
}

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void finalClosurePlotter(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Open data files for pp and PbPb data
  TFile *ppFile = TFile::Open("data/PbPbMC_GenGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_xjBins_sube0_wtaAxis_jet100trigger_JECv6_processed_2019-09-26.root");

  TFile *pbpbFile = TFile::Open("data/PbPbMC2018_RecoReco_akFlowPuCs4PFJet_noUncOrInc_5eveMix_allCorrections_fixedFile_useTrackingScale_tunedRange_lowSpillSub_mixingScale16_processed_2019-10-07.root");

  
  TFile *ppUncertaintyFile = TFile::Open("uncertainties/systematicUncertaintyForPpMC_RecoGen_20eveMix_mcMode_2019-10-21.root");
  TFile *pbpbUncertaintyFile = TFile::Open("uncertainties/systematicUncertaintyForPbPbMC_RecoGen_5eveMix_mcMode_2019-10-21.root");

  
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
  ppHistograms->SetLoadTrackLeadingJetCorrelations(true);
  ppHistograms->SetLoadTrackInclusiveJetCorrelations(true);
  ppHistograms->SetAsymmetryBinRange(0,nAsymmetryBins);
  ppHistograms->LoadProcessedHistograms();
  
  pbpbHistograms->SetLoadTrackLeadingJetCorrelations(true);
  pbpbHistograms->SetLoadTrackInclusiveJetCorrelations(true);
  pbpbHistograms->SetAsymmetryBinRange(0,nAsymmetryBins);
  pbpbHistograms->LoadProcessedHistograms();
  
  // Load the systematic uncertainties from the files
  ppUncertaintyProvider->ReadSystematicFile(ppUncertaintyFile);
  pbpbUncertaintyProvider->ReadSystematicFile(pbpbUncertaintyFile);
  
  // Plot the figures using Xiao's plotting macro
  plotJetShapeXiao(ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kTrackLeadingJet);
  plotJetShapeXiao(ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kTrackSubleadingJet);
  //plotJetShapeXiao(ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kPtWeightedTrackInclusiveJet);
}
