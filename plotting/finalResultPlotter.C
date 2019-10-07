#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "stackHistXiao.h"
#include "xCanvas.h"
#include "JDrawer.h"

void plotJetShapeXiao(DijetHistogramManager *ppHistograms, DijetHistogramManager *pbpbHistograms, JffCorrector *ppUncertaintyProvider, JffCorrector *pbpbUncertaintyProvider, int iJetTrack){

  const int nTrackPtBins = pbpbHistograms->GetNTrackPtBins();
  const int nCentralityBins = pbpbHistograms->GetNCentralityBins();
  const int nAsymmetryBins = pbpbHistograms->GetNAsymmetryBins();
  
  const char* jetShapeTitle[] = {"Leading jet shape","Subleading jet shape","Jet shape"};
  const char* jetShapeSaveName[] = {"trackLeadingJet","trackSubleadingJet","trackInclusiveJet"};
  const char* xjString[] = {"0.0 < x_{j} < 0.6","0.6 < x_{j} < 0.8","0.8 < x_{j} < 1.0","x_{j} inclusive"};
  const char* asymmetrySaveName[] = {"_A=0v0-0v6","_A=0v6-0v8","_A=0v8-1v0",""};
  
  bool drawInclusive = false;
  bool drawInclusiveRatio = false;
  bool drawUncertainties = true;
  bool monteCarloLabels = false;
  int firstAsymmetryBin = 0;  // Set here 0 to process asymmetry bins and nAsymmetrybins to disable them
    
  // Draw the comparison only for leading jets
  if(iJetTrack != DijetHistogramManager::kPtWeightedTrackLeadingJet){
    drawInclusive = false;
    drawInclusiveRatio = false;
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
  
  // Open a file to include results from inclusive jet analysis HIN-16-020
  TFile *inclusiveResultFile = TFile::Open("data/publishedResults/JS5TeV_HIN_16_020.root");
  TH1D *sumHistogramPbPbInclusive[nCentralityBins];
  TH1D *sumHistogramPpInclusive;
  TH1D *inclusiveRatio[nCentralityBins];
  
  // Read the pT summed jet shape histograms
  for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      sumHistogramPbPb[iCentrality][iAsymmetry] = (TH1D*) pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, 0)->Clone(Form("sumPbPb%d",iCentrality));
      jetShapeArray[iCentrality][0][iAsymmetry] = pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, 0);
    }
    sumHistogramPp[iAsymmetry] = (TH1D*) ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, 0)->Clone("sumPp");
    jetShapeArray[nCentralityBins][0][iAsymmetry] = ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, 0);
  }
  
  // Sum the pT:s
  for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
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
  
  // Read the histograms from the inclusive jet shape analysis
  sumHistogramPpInclusive = (TH1D*) inclusiveResultFile->Get("JS_pp_0")->Clone(Form("inclusiveJetShape%d", iJetTrack));
  
  for(int iTrackPt = 1; iTrackPt < 9; iTrackPt++){
    helperHistogram = (TH1D*) inclusiveResultFile->Get(Form("JS_pp_%d",iTrackPt));
    sumHistogramPpInclusive->Add(helperHistogram);
  }
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    sumHistogramPbPbInclusive[iCentrality] = (TH1D*) inclusiveResultFile->Get(Form("JS_pb_0_%d",iCentrality))->Clone(Form("inclusiveJetShape%d%d", iJetTrack, iCentrality));
    
    for(int iTrackPt = 1; iTrackPt < 9; iTrackPt++){
      helperHistogram = (TH1D*) inclusiveResultFile->Get(Form("JS_pb_%d_%d",iTrackPt,iCentrality));
      sumHistogramPbPbInclusive[iCentrality]->Add(helperHistogram);
    }
  }
  
  // Normalize the histograms from the inclusive jet shape analysis
  sumHistogramPpInclusive->Scale(1.0/sumHistogramPpInclusive->Integral(1, sumHistogramPpInclusive->FindBin(0.99) ,"width"));
  sumHistogramPpInclusive->SetLineColor(kMagenta);
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    sumHistogramPbPbInclusive[iCentrality]->Scale(1.0/sumHistogramPbPbInclusive[iCentrality]->Integral(1, sumHistogramPbPbInclusive[iCentrality]->FindBin(0.99) ,"width"));
    sumHistogramPbPbInclusive[iCentrality]->SetLineColor(kMagenta);
    
    // Also calculate the ratio of the inclusive histograms
    inclusiveRatio[iCentrality] = (TH1D*) sumHistogramPbPbInclusive[iCentrality]->Clone(Form("inclusiveRatio%d",iCentrality));
    inclusiveRatio[iCentrality]->Divide(sumHistogramPpInclusive);
  }
  
  // Normalize the jet shape histograms
  for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
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
    
    // Read the uncertainties for the pT summed jet shapes and scale them for jet shapes
    uncertaintyHistogramPp[iAsymmetry] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, nTrackPtBins, iAsymmetry);
    uncertaintyHistogramPp[iAsymmetry]->Scale(1.0/shapeIntegralPp[iAsymmetry]);
    
    // Set the bin errors in the sumUncertainty histograms to match the uncertainties read from the file
    for(int iBin = 1; iBin <= sumUncertainty[nAsymmetryBins][iAsymmetry]->GetNbinsX(); iBin++){
      sumUncertainty[nCentralityBins][iAsymmetry]->SetBinError(iBin, uncertaintyHistogramPp[iAsymmetry]->GetBinContent(iBin));
    }
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      uncertaintyHistogramPbPb[iCentrality][iAsymmetry] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, nTrackPtBins, iAsymmetry);
      uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iCentrality][iAsymmetry]);
      
      // Set the bin errors in the sumUncertainty histograms to match the uncertainties read from the file
      for(int iBin = 1; iBin <= sumUncertainty[iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
        
        sumUncertainty[iCentrality][iAsymmetry]->SetBinError(iBin, uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->GetBinContent(iBin));
      }
      
    }
  }
  
  TH1D* ratioHistogram[nCentralityBins][nAsymmetryBins+1];
  TH1D* asymmetryRatioHistogram[nCentralityBins+1][nAsymmetryBins];
  
  // Create uncertainty histograms for the ratio using standard jet shape binning
  TH1D* ratioUncertainty[nCentralityBins][nAsymmetryBins+1];
  TH1D* asymmetryRatioUncertainty[nCentralityBins+1][nAsymmetryBins];
  double ppValue, pbpbValue, ratioValue;
  double ppUncertainty, pbpbUncertainty, ratioUncertaintyValue;
  
  for(int iAsymmetry = firstAsymmetryBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      ratioUncertainty[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("ratioUncertainty%d%d", iCentrality, iAsymmetry));
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
      
    }
  }
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    for(int iAsymmetry = firstAsymmetryBin; iAsymmetry < nAsymmetryBins; iAsymmetry++){
      asymmetryRatioUncertainty[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[0][0]->Clone(Form("asymmetryRatioUncertainty%d%d", iCentrality, iAsymmetry));
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
      jetShapeStack[iCentrality][iAsymmetry]->setRange(0.005, 30, "y");
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
      ratioUncertainty[iCentrality][iAsymmetry]->SetMarkerColor(0);
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
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitle("#rho(#Deltar)");
      jetShapeStack[iCentrality][iAsymmetry]->hst->Draw();
      if(iAsymmetry == nAsymmetryBins && iJetTrack == 2 && drawInclusive){
        sumHistogramPbPbInclusive[iCentrality]->Draw("same");
      }
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
      ratioHistogram[iCentrality][iAsymmetry]->SetAxisRange(0., 3.2, "Y");
      ratioHistogram[iCentrality][iAsymmetry]->SetAxisRange(0, .99, "X");
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
        ratioHistogram[iCentrality][iAsymmetry]->GetYaxis()->SetTitle("#rho(#Deltar)_{PbPb}/#rho(#Deltar)_{pp}");
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
      
      if(drawInclusiveRatio && iAsymmetry == nAsymmetryBins){
        inclusiveRatio[iCentrality]->Draw("same");
      }
    }
    bigCanvas[iAsymmetry]->CD(1);
    gPad->SetLogy();
    jetShapeStack[4][iAsymmetry]->drawStack();
    jetShapeStack[4][iAsymmetry]->hst->GetYaxis()->SetNdivisions(505);
    jetShapeStack[4][iAsymmetry]->hst->GetYaxis()->SetLabelSize(0.08);
    jetShapeStack[4][iAsymmetry]->hst->GetYaxis()->SetTitleOffset(0.9);
    jetShapeStack[4][iAsymmetry]->hst->GetYaxis()->SetTitleSize(0.1);
    jetShapeStack[4][iAsymmetry]->hst->GetYaxis()->SetTitle("#rho(#Deltar)");
    if(iAsymmetry == nAsymmetryBins && iJetTrack == 2 && drawInclusive){
      sumHistogramPpInclusive->Draw("same");
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
    
    if(drawInclusiveRatio && iAsymmetry == nAsymmetryBins){
      lt5->AddEntry(inclusiveRatio[0], "JHEP 05(2018)006","lpfe");
      bigCanvas[iAsymmetry]->CD(10);
      lt5->Draw();
    }
    
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
    double jetShapeTitlePosition[] = {0.52,0.5,0.52};
    double xjPosition[] = {0.76,0.77,0.76};
    double overallShift[] = {-0.02,-0.02,-0.02,0};
    
    bigCanvas[iAsymmetry]->cd(0);
    mainTitle[iAsymmetry]->SetTextFont(62);
    mainTitle[iAsymmetry]->SetTextSize(0.035);
    mainTitle[iAsymmetry]->DrawLatexNDC(cmsPosition[iJetTrack/3]+overallShift[iAsymmetry], 0.94, "CMS");
    
    mainTitle[iAsymmetry]->SetTextFont(42);
    mainTitle[iAsymmetry]->SetTextSize(0.035);
    mainTitle[iAsymmetry]->DrawLatexNDC(jetShapeTitlePosition[iJetTrack/3]+overallShift[iAsymmetry], 0.94, jetShapeTitle[iJetTrack/3]);
    
    mainTitle[iAsymmetry]->SetTextFont(42);
    mainTitle[iAsymmetry]->SetTextSize(0.035);
    mainTitle[iAsymmetry]->DrawLatexNDC(xjPosition[iJetTrack/3]+overallShift[iAsymmetry], 0.94, xjString[iAsymmetry]);
    
    mainTitle[iAsymmetry]->SetTextSize(0.03);
    mainTitle[iAsymmetry]->DrawLatexNDC(0.41, 0.9, "pp 320 pb^{-1} (5.02 TeV)  PbPb 1.7 nb^{-1} (5.02 TeV)");
    mainTitle[iAsymmetry]->SetTextSize(0.025);
    mainTitle[iAsymmetry]->DrawLatexNDC(0.375, 0.86, "anti-k_{T} R = 0.4, |#eta_{jet}| < 1.6, p_{T,1} > 120 GeV, p_{T,2} > 50 GeV, #Delta#phi_{1,2} > #frac{5#pi}{6}");
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
    bigCanvas[iAsymmetry]->SaveAs(Form("figures/finalJetShape_%s%s_newErrors.pdf",jetShapeSaveName[iJetTrack/3],asymmetrySaveName[iAsymmetry]));
    //bigCanvas->SaveAs("js_dr_normal_v3.eps");
    //bigCanvas->SaveAs("js_dr_normal_v3.pdf");
    
  } // Asymmetry loop for drawing the big canvas
}

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void finalResultPlotter(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Open data files for pp and PbPb data
  TFile *ppFile = TFile::Open("data/ppData2017_highForest_pfJets_20EventsMixed_xjBins_finalTrackCorr_JECv4_wtaAxis_allCorrections_processed_2019-09-28.root");
  // data/ppData2017_highForest_pfJets_20EventsMixed_xjBins_finalTrackCorr_JECv4_wtaAxis_allCorrections_processed_2019-09-28.root
  // data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_processed_2019-08-13.root
  // data/dijet_pp_highForest_pfJets_noUncOrInc_allCorrections_wtaAxis_processed_2019-07-13.root
  // data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_noUncorr_20EventsMixed_JECv4_allCorrections_processed_2019-09-28.root
  TFile *pbpbFile = TFile::Open("data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_allCorrections_noStatErrorFromCorrections_lowPtResidualTrack_processed_2019-10-07_fiveJobsMissing.root");
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_forSubleading_allCorrections_wtaAxis_JECv6_processed_2019-09-26.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixingFromSubeNon0_wtaAxis_allCorrections_JECv6_processed_2019-09-26.root
  // data/PbPbMC_RecoGen_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_JECv6_allCorrections_processed_2019-10-04.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_5eveMix_calo100Trigger_JECv6_finalTrack_allCorrections_wtaAxis_processed_2019-09-26.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_allCorrections_lowPtResidualTrack_processed_2019-10-01_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_calo80Trigger_xjBins_wtaAxis_allCorrections_JECv5b_processed_2019-09-05.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_modifiedSeagull_noErrorJff_averagePeakMixing_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_testNoJffCorr_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_allCorrections_modifiedSeagull_wtaAxis_JECv4_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_processed_2019-07-05.root
  
  TFile *ppUncertaintyFile = TFile::Open("systematicUncertaintyForPp_20percentSpillJff_2019-09-30.root");
  // systematicUncertaintyForPp_20percentSpillJff_2019-09-30.root
  TFile *pbpbUncertaintyFile = TFile::Open("systematicUncertaintyForPp_15percentSpill20Jff_2019-10-01.root");
  // systematicUncertaintyForPythiaHydjetRecoGen_mcMode_2019-10-06.root
  // systematicUncertaintyForPythiaHydjetRecoGen_mcMode_2019-10-05.root
  // systematicUncertaintyForPp_15percentSpill20Jff_2019-10-01.root
  
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
  plotJetShapeXiao(ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kPtWeightedTrackLeadingJet);
  plotJetShapeXiao(ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kPtWeightedTrackSubleadingJet);
  return;
  
  // Temporary: Get the ratio between pp and PbPb jet shape
  TH1D *sumHistogramPbPb[2][nCentralityBins][nAsymmetryBins+1];
  TH1D *sumHistogramPp[2][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPbPb[2][nCentralityBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPp[2][nAsymmetryBins+1];
  TH1D *helperHistogram;
  double shapeIntegralPbPb[2][nCentralityBins][nAsymmetryBins+1];
  double shapeIntegralPp[2][nAsymmetryBins+1];
  
  // Read the pT summed jet shape histograms
  for(int iJetTrack = 0; iJetTrack < 2; iJetTrack++){
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        sumHistogramPbPb[iJetTrack][iCentrality][iAsymmetry] = (TH1D*) pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, jetTrackIndex[iJetTrack], iAsymmetry, iCentrality, lowestTrackPtBin)->Clone(Form("sumPbPb%d",iCentrality));
      }
      sumHistogramPp[iJetTrack][iAsymmetry] = (TH1D*) ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, jetTrackIndex[iJetTrack], iAsymmetry, 0, lowestTrackPtBin)->Clone("sumPp");
    }
  }
  
  // Sum the pT:s
  for(int iJetTrack = 0; iJetTrack < 2; iJetTrack++){
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iTrackPt = lowestTrackPtBin+1; iTrackPt < nTrackPtBins; iTrackPt++){
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          helperHistogram = (TH1D*)pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, jetTrackIndex[iJetTrack], iAsymmetry, iCentrality, iTrackPt)->Clone();
          sumHistogramPbPb[iJetTrack][iCentrality][iAsymmetry]->Add(helperHistogram);
        }
        helperHistogram = (TH1D*)ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, jetTrackIndex[iJetTrack], iAsymmetry, 0, iTrackPt)->Clone();
        sumHistogramPp[iJetTrack][iAsymmetry]->Add(helperHistogram);
      } // track pT
    } // Asymmetry loop
  } // Jet track type loop
  
  // Normalize the jet shape histograms
  for(int iJetTrack = 0; iJetTrack < 2; iJetTrack++){
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        shapeIntegralPbPb[iJetTrack][iCentrality][iAsymmetry] = sumHistogramPbPb[iJetTrack][iCentrality][iAsymmetry]->Integral(1,sumHistogramPbPb[iJetTrack][iCentrality][iAsymmetry]->FindBin(0.99),"width");
        sumHistogramPbPb[iJetTrack][iCentrality][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iJetTrack][iCentrality][iAsymmetry]);
      }
      
      shapeIntegralPp[iJetTrack][iAsymmetry] = sumHistogramPp[iJetTrack][iAsymmetry]->Integral(1,sumHistogramPp[iJetTrack][iAsymmetry]->FindBin(0.99),"width");
      sumHistogramPp[iJetTrack][iAsymmetry]->Scale(1.0/shapeIntegralPp[iJetTrack][iAsymmetry]);
      
      // Read the uncertainties for the pT summed jet shapes and scale them for jet shapes
      uncertaintyHistogramPp[iJetTrack][iAsymmetry] = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(jetTrackIndex[iJetTrack], 0, nTrackPtBins, iAsymmetry);
      uncertaintyHistogramPp[iJetTrack][iAsymmetry]->Scale(1.0/shapeIntegralPp[iJetTrack][iAsymmetry]);
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        uncertaintyHistogramPbPb[iJetTrack][iCentrality][iAsymmetry] = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(jetTrackIndex[iJetTrack], iCentrality, nTrackPtBins, iAsymmetry);
        uncertaintyHistogramPbPb[iJetTrack][iCentrality][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iJetTrack][iCentrality][iAsymmetry]);
      }
    } // Asymmetry loop
  } // Jet track type loop
  
  TH1D* ratioHistogram[2][nCentralityBins][nAsymmetryBins+1];
  TH1D* asymmetryRatioHistogram[2][nCentralityBins+1][nAsymmetryBins];
  
  // Create uncertainty histograms for the ratio using standard jet shape binning
  TH1D* ratioUncertainty[2][nCentralityBins][nAsymmetryBins+1];
  TH1D* asymmetryRatioUncertainty[2][nCentralityBins+1][nAsymmetryBins];
  for(int iJetTrack = 0; iJetTrack < 2; iJetTrack++){
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        ratioUncertainty[iJetTrack][iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iJetTrack][iCentrality][iAsymmetry]->Clone(Form("ratioUncertainty%d%d%d", iJetTrack, iCentrality, iAsymmetry));
      } // Centrality loop
    } // Asymmetry loop
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
        asymmetryRatioUncertainty[iJetTrack][iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iJetTrack][0][0]->Clone(Form("asymmetryRatioUncertainty%d%d%d", iJetTrack, iCentrality, iAsymmetry));
      } // Asymmetry loop
    } // Centrality loop
  } // Jet track type loop
  
  // Open a file to include results from inclusive jet analysis HIN-16-020
  TFile *inclusiveResultFile = TFile::Open("data/publishedResults/JS5TeV_HIN_16_020.root");
  TH1D *sumHistogramPbPbInclusive[nCentralityBins];
  TH1D *sumHistogramPpInclusive;
  TH1D *inclusiveRatio[nCentralityBins];
  
  // Read the histograms from the inclusive jet shape analysis
  sumHistogramPpInclusive = (TH1D*) inclusiveResultFile->Get("JS_pp_0")->Clone("inclusiveJetShapePp");
  
  for(int iTrackPt = 1; iTrackPt < 9; iTrackPt++){
    helperHistogram = (TH1D*) inclusiveResultFile->Get(Form("JS_pp_%d",iTrackPt));
    sumHistogramPpInclusive->Add(helperHistogram);
  }
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    sumHistogramPbPbInclusive[iCentrality] = (TH1D*) inclusiveResultFile->Get(Form("JS_pb_0_%d",iCentrality))->Clone(Form("inclusiveJetShapePbPb%d", iCentrality));
    
    for(int iTrackPt = 1; iTrackPt < 9; iTrackPt++){
      helperHistogram = (TH1D*) inclusiveResultFile->Get(Form("JS_pb_%d_%d",iTrackPt,iCentrality));
      sumHistogramPbPbInclusive[iCentrality]->Add(helperHistogram);
    }
  }
  
  // Normalize the histograms from the inclusive jet shape analysis
  sumHistogramPpInclusive->Scale(1.0/sumHistogramPpInclusive->Integral(1, sumHistogramPpInclusive->FindBin(0.99) ,"width"));
  sumHistogramPpInclusive->SetLineColor(kMagenta);
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    sumHistogramPbPbInclusive[iCentrality]->Scale(1.0/sumHistogramPbPbInclusive[iCentrality]->Integral(1, sumHistogramPbPbInclusive[iCentrality]->FindBin(0.99) ,"width"));
    sumHistogramPbPbInclusive[iCentrality]->SetLineColor(kMagenta);
    
    // Also calculate the ratio of the inclusive histograms
    inclusiveRatio[iCentrality] = (TH1D*) sumHistogramPbPbInclusive[iCentrality]->Clone(Form("inclusiveRatio%d",iCentrality));
    inclusiveRatio[iCentrality]->Divide(sumHistogramPpInclusive);
  }
  
  JDrawer *drawer = new JDrawer();
  drawer->SetCanvasSize(500,500);
  
  TLegend *legend;
  double ppValue, pbpbValue, ratioValue;
  double ppUncertainty, pbpbUncertainty, ratioUncertaintyValue;
  
  TLine *oneLine = new TLine(0,1,1,1);
  oneLine->SetLineStyle(2);
  
  if(drawJetShapePptoPbPbRatio){
    
    // Prepare the uncertainty histograms for the ratio
    for(int iJetTrack = 0; iJetTrack < 2; iJetTrack++){
      for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          ratioHistogram[iJetTrack][iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iJetTrack][iCentrality][iAsymmetry]->Clone(Form("ratio%d%d%d", iJetTrack, iCentrality, iAsymmetry));
          ratioHistogram[iJetTrack][iCentrality][iAsymmetry]->Divide(sumHistogramPp[iJetTrack][iAsymmetry]);
          
          // Calculate the systemtic uncertainty for the ratio
          for(int iBin = 1; iBin < ratioUncertainty[iJetTrack][iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
            ppValue = sumHistogramPp[iJetTrack][iAsymmetry]->GetBinContent(iBin);
            pbpbValue = sumHistogramPbPb[iJetTrack][iCentrality][iAsymmetry]->GetBinContent(iBin);
            ppUncertainty = uncertaintyHistogramPp[iJetTrack][iAsymmetry]->GetBinContent(iBin);
            pbpbUncertainty = uncertaintyHistogramPbPb[iJetTrack][iCentrality][iAsymmetry]->GetBinContent(iBin);
            ratioValue = ratioHistogram[iJetTrack][iCentrality][iAsymmetry]->GetBinContent(iBin);
            ratioUncertaintyValue = TMath::Sqrt(TMath::Power(pbpbUncertainty/ppValue,2)+TMath::Power(pbpbValue*ppUncertainty/TMath::Power(ppValue,2),2));
            ratioUncertainty[iJetTrack][iCentrality][iAsymmetry]->SetBinContent(iBin,ratioValue);
            ratioUncertainty[iJetTrack][iCentrality][iAsymmetry]->SetBinError(iBin,ratioUncertaintyValue);
          } // Bin loop
        } // Centrality loop
      } // Asymmetry loop
    } // Jet track type loop
    
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        
        ratioHistogram[0][iCentrality][iAsymmetry]->GetXaxis()->SetRangeUser(0,1);
        ratioHistogram[0][iCentrality][iAsymmetry]->GetYaxis()->SetRangeUser(0,3);
        ratioHistogram[0][iCentrality][iAsymmetry]->SetLineColor(kBlack);
        ratioUncertainty[0][iCentrality][iAsymmetry]->GetXaxis()->SetRangeUser(0,1);
        ratioUncertainty[0][iCentrality][iAsymmetry]->GetYaxis()->SetRangeUser(0,3);
        ratioUncertainty[0][iCentrality][iAsymmetry]->SetFillColorAlpha(29,0.25);
        
        ratioHistogram[1][iCentrality][iAsymmetry]->SetLineColor(kRed);
        ratioUncertainty[1][iCentrality][iAsymmetry]->SetFillColorAlpha(46,0.25);
        
        drawer->DrawHistogram(ratioUncertainty[0][iCentrality][iAsymmetry],"#Deltar","#rho(#Deltar)_{PbPb} / #rho(#Deltar)_{pp}", " ","E2");
        ratioUncertainty[1][iCentrality][iAsymmetry]->Draw("same,E2");
        oneLine->Draw();
        
        ratioHistogram[0][iCentrality][iAsymmetry]->Draw("same");
        ratioHistogram[1][iCentrality][iAsymmetry]->Draw("same");
        
        if(iAsymmetry == nAsymmetryBins){
          inclusiveRatio[iCentrality]->SetLineColor(kBlue);
          inclusiveRatio[iCentrality]->Draw("same");
        }
        
        legend = new TLegend(0.21,0.63,0.41,0.88);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->SetHeader(Form("C: %.0f-%.0f%s",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1], asymmetryString[iAsymmetry].Data()));
        legend->AddEntry(ratioHistogram[0][iCentrality][iAsymmetry],"Leading jets","l");
        legend->AddEntry(ratioHistogram[1][iCentrality][iAsymmetry],"Subleading jets","l");
        if(iAsymmetry == nAsymmetryBins) legend->AddEntry(inclusiveRatio[iCentrality], "JHEP 05(2018)006");
        legend->Draw();
        
        if(saveFigures){
          if(iAsymmetry == nAsymmetryBins){
            gPad->GetCanvas()->SaveAs(Form("figures/jetShapeRatioComparison_%s_C=%.0f-%.0f_pT%d.pdf", pbpbHistograms->GetJetTrackHistogramName(jetTrackIndex[0]), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1],lowestTrackPtBin));
          } else {
            saveName = Form("figures/jetShapeRatioComparison_%s_A=%.1f-%.1f_C=%.0f-%.0f", pbpbHistograms->GetJetTrackHistogramName(jetTrackIndex[0]), asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1], centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]);
            saveName.ReplaceAll(".","v");
            gPad->GetCanvas()->SaveAs(Form("%s.%s",saveName.Data(),figureFormat));
          }
        } // Saving figures
        
      } // Centrality loop
    } // Asymmetry loop
  } // Ratio between pp and PbPb in different asymmetry bins
  
  if(drawJetShapeSymmetricAsymmetricRatio){
    
    // Prepare the uncertainties for asymmetry ratios
    for(int iJetTrack = 0; iJetTrack < 2; iJetTrack++){
      for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
        for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
          
          if(iCentrality == nCentralityBins){
            asymmetryRatioHistogram[iJetTrack][iCentrality][iAsymmetry] = (TH1D*) sumHistogramPp[iJetTrack][iAsymmetry]->Clone(Form("ratioAsymmetry%d%d%d", iJetTrack, iCentrality, iAsymmetry));
            asymmetryRatioHistogram[iJetTrack][iCentrality][iAsymmetry]->Divide(sumHistogramPp[iJetTrack][nAsymmetryBins]);
          } else {
            asymmetryRatioHistogram[iJetTrack][iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iJetTrack][iCentrality][iAsymmetry]->Clone(Form("ratioAsymmetry%d%d%d", iJetTrack, iCentrality, iAsymmetry));
            asymmetryRatioHistogram[iJetTrack][iCentrality][iAsymmetry]->Divide(sumHistogramPbPb[iJetTrack][iCentrality][nAsymmetryBins]);
          }
          
          // Calculate the systemtic uncertainty for the ratio
          for(int iBin = 1; iBin < asymmetryRatioUncertainty[iJetTrack][iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
            if(iCentrality == nCentralityBins){
              ppValue = sumHistogramPp[iJetTrack][nAsymmetryBins]->GetBinContent(iBin);
              pbpbValue = sumHistogramPp[iJetTrack][iAsymmetry]->GetBinContent(iBin);
              ppUncertainty = uncertaintyHistogramPp[iJetTrack][nAsymmetryBins]->GetBinContent(iBin);
              pbpbUncertainty = uncertaintyHistogramPp[iJetTrack][iAsymmetry]->GetBinContent(iBin);
            } else {
              ppValue = sumHistogramPbPb[iJetTrack][iCentrality][nAsymmetryBins]->GetBinContent(iBin);
              pbpbValue = sumHistogramPbPb[iJetTrack][iCentrality][iAsymmetry]->GetBinContent(iBin);
              ppUncertainty = uncertaintyHistogramPbPb[iJetTrack][iCentrality][nAsymmetryBins]->GetBinContent(iBin);
              pbpbUncertainty = uncertaintyHistogramPbPb[iJetTrack][iCentrality][iAsymmetry]->GetBinContent(iBin);
            }
            ratioValue = asymmetryRatioHistogram[iJetTrack][iCentrality][iAsymmetry]->GetBinContent(iBin);
            ratioUncertaintyValue = TMath::Sqrt(TMath::Power(pbpbUncertainty/ppValue,2)+TMath::Power(pbpbValue*ppUncertainty/TMath::Power(ppValue,2),2));
            asymmetryRatioUncertainty[iJetTrack][iCentrality][iAsymmetry]->SetBinContent(iBin,ratioValue);
            asymmetryRatioUncertainty[iJetTrack][iCentrality][iAsymmetry]->SetBinError(iBin,ratioUncertaintyValue);
          }
        } // Centrality loop
      } // Asymmetry loop
    } // Jet track loop
    
    // Draw the asymmetry ratios
    for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
        
        asymmetryRatioHistogram[0][iCentrality][iAsymmetry]->GetXaxis()->SetRangeUser(0,1);
        asymmetryRatioHistogram[0][iCentrality][iAsymmetry]->GetYaxis()->SetRangeUser(0,3);
        asymmetryRatioUncertainty[0][iCentrality][iAsymmetry]->GetXaxis()->SetRangeUser(0,1);
        asymmetryRatioUncertainty[0][iCentrality][iAsymmetry]->GetYaxis()->SetRangeUser(0,3);
        asymmetryRatioUncertainty[0][iCentrality][iAsymmetry]->SetFillColorAlpha(29,0.25);
        
        drawer->DrawHistogram(asymmetryRatioUncertainty[0][iCentrality][iAsymmetry],"#Deltar",Form("#rho(#Deltar) %.1f < x_{j} < %.1f / #rho(#Deltar) all x_{j}", asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1]), " ","E2");
        oneLine->Draw();
        asymmetryRatioHistogram[0][iCentrality][iAsymmetry]->SetLineColor(kBlack);
        asymmetryRatioHistogram[0][iCentrality][iAsymmetry]->Draw("same");
        
        legend = new TLegend(0.25,0.73,0.45,0.88);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        if(iCentrality == nCentralityBins){
          legend->SetHeader("pp  Track-subleading jet");
        } else {
          legend->SetHeader(Form("C: %.0f-%.0f  Track-subleading jet",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
        }
        legend->AddEntry(asymmetryRatioHistogram[0][iCentrality][iAsymmetry],Form("%.1f < x_{j} < %.1f / all x_{j}", asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1]),"l");
        legend->Draw();
        
        if(saveFigures){
          if(iCentrality == nCentralityBins){
            saveName = Form("figures/jetShapeAsymmetryRatio_%s_pp%s", pbpbHistograms->GetJetTrackHistogramName(jetTrackIndex[0]), asymmetrySaveString[iAsymmetry].Data());
            saveName.ReplaceAll(".","v");
            gPad->GetCanvas()->SaveAs(Form("%s.%s",saveName.Data(),figureFormat));
          } else {
            saveName = Form("figures/jetShapeAsymmetryRatio_%s_C=%.0f-%.0f%s", pbpbHistograms->GetJetTrackHistogramName(jetTrackIndex[0]), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], asymmetrySaveString[iAsymmetry].Data());
            saveName.ReplaceAll(".","v");
            gPad->GetCanvas()->SaveAs(Form("%s.%s",saveName.Data(),figureFormat));
          }
        } // Saving figures
        
      } // Centrality loop
    } // Asymmetry loop
  } // Draw asymmetry ratio
  
}
