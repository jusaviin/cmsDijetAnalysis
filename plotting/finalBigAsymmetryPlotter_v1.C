#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "stackHistXiao.h"
#include "xCanvas.h"
#include "JDrawer.h"

void plotJetShapeBigAsymmetry(DijetHistogramManager *ppHistograms, DijetHistogramManager *pbpbHistograms, JffCorrector *ppUncertaintyProvider, JffCorrector *pbpbUncertaintyProvider, int iJetTrack){
  
  //  Note, changing the first track pT bin properly requires some manual work. Not all is done automagically.
  const int iFirstTrackPtBin = 0;  // Can change the lower bound of the drawn track pT bins. Default is 0.
  bool printIntegrals = false;     // Print the integrals of all the distributions to slides
  
  const int nTrackPtBins = pbpbHistograms->GetNTrackPtBins();
  const int nCentralityBins = pbpbHistograms->GetNCentralityBins();
  const int nAsymmetryBins = pbpbHistograms->GetNAsymmetryBins();
  
  const char* jetShapeTitle[] = {"Leading jet shape","Subleading jet shape","Jet shape"};
  const char* jetShapeSaveName[] = {"trackLeadingJet","trackSubleadingJet","trackInclusiveJet"};
  const char* xjString[] = {"0.0 < x_{j} < 0.6","0.6 < x_{j} < 0.8","0.8 < x_{j} < 1.0","x_{j} inclusive"};
  const char* asymmetrySaveName[] = {"_A=0v0-0v6","_A=0v6-0v8","_A=0v8-1v0",""};
  const char* saveString = {"JetShape"};
  
  const bool drawUncertainties = true;
  const bool monteCarloLabels = false;
  const bool normalizeJetShape = true;         // True: draw rho. False: draw P.
  const bool drawExtraRatio = false;           // Draw illustration of third jet effects to the ratio
  const bool saveHistogramsForHepData = false; // Save the plotted histograms to a file for HepData submission
  
  const bool addPreliminaryTag = false;
  
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
  TH1D *jetShapeForIntegral[nCentralityBins+1][nTrackPtBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPp[nAsymmetryBins+1];
  TH1D *sumUncertainty[nCentralityBins+1][nAsymmetryBins+1];
  TH1D *helperHistogram;
  double shapeIntegralPbPb[nCentralityBins][nAsymmetryBins+1];
  double shapeIntegralPp[nAsymmetryBins+1];
  double integralAfterNormalization[nCentralityBins+1][nTrackPtBins][nAsymmetryBins+1];
  double integralErrorAfterNormalization[nCentralityBins+1][nTrackPtBins][nAsymmetryBins+1];
  
  // Read the pT summed jet shape histograms
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      sumHistogramPbPb[iCentrality][iAsymmetry] = (TH1D*) pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iFirstTrackPtBin)->Clone(Form("sumPbPb%d",iCentrality));
      jetShapeArray[iCentrality][iFirstTrackPtBin][iAsymmetry] = pbpbHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, iCentrality, iFirstTrackPtBin);
    }
    sumHistogramPp[iAsymmetry] = (TH1D*) ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, iFirstTrackPtBin)->Clone("sumPp");
    jetShapeArray[nCentralityBins][iFirstTrackPtBin][iAsymmetry] = ppHistograms->GetHistogramJetShape(DijetHistogramManager::kJetShape, iJetTrack, iAsymmetry, 0, iFirstTrackPtBin);
  }
  
  // Sum the pT:s
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iTrackPt = iFirstTrackPtBin+1; iTrackPt < nTrackPtBins; iTrackPt++){
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
        for(int iTrackPt = iFirstTrackPtBin; iTrackPt < nTrackPtBins; iTrackPt++){
          jetShapeArray[iCentrality][iTrackPt][iAsymmetry]->Scale(1.0/shapeIntegralPbPb[iCentrality][iAsymmetry]);
        }
      }
      
    }
    
    shapeIntegralPp[iAsymmetry] = sumHistogramPp[iAsymmetry]->Integral(1,sumHistogramPp[iAsymmetry]->FindBin(0.99),"width");
    if(normalizeJetShape) sumHistogramPp[iAsymmetry]->Scale(1.0/shapeIntegralPp[iAsymmetry]);
    
    // To draw the systematic uncertainties, clone the sum histograms to uncertainty histograms
    sumUncertainty[nCentralityBins][iAsymmetry] = (TH1D*) sumHistogramPp[iAsymmetry]->Clone(Form("sumUncertaintyPp%d", iAsymmetry));
    
    if(normalizeJetShape){
      for(int iTrackPt = iFirstTrackPtBin; iTrackPt < nTrackPtBins; iTrackPt++){
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
  
  // Calculate the jet shape integrals and errors taking into account the systematic uncertainties
  TH1D *systematicErrorHistogram;
  TH1D *resultPlusSystematicError;
  TH1D *resultMinusSystematicError;
  int oneBin;
  double binContent;
  double statisticalError, systematicError, totalError;
  double systematicPlusIntegral, systematicMinusIntegral;
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        // Clone the jet shape histogram including statistical uncertainties
        jetShapeForIntegral[iCentrality][iTrackPt][iAsymmetry] = (TH1D*)jetShapeArray[iCentrality][iTrackPt][iAsymmetry]->Clone(Form("jetShapeForIntegral%d%d%d", iAsymmetry, iCentrality, iTrackPt));
        
        resultPlusSystematicError = (TH1D*)jetShapeArray[iCentrality][iTrackPt][iAsymmetry]->Clone(Form("resultPlus%d%d%d", iAsymmetry, iCentrality, iTrackPt));
        
        resultMinusSystematicError = (TH1D*)jetShapeArray[iCentrality][iTrackPt][iAsymmetry]->Clone(Form("resultMinus%d%d%d", iAsymmetry, iCentrality, iTrackPt));
        
        // Find the systematic uncertainties for this bin
        if(iCentrality == nCentralityBins){
          systematicErrorHistogram = ppUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, 0, iTrackPt, iAsymmetry);
          if(normalizeJetShape) systematicErrorHistogram->Scale(1.0/shapeIntegralPp[iAsymmetry]);
        } else {
          systematicErrorHistogram = pbpbUncertaintyProvider->GetJetShapeSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry);
          if(normalizeJetShape) systematicErrorHistogram->Scale(1.0/shapeIntegralPbPb[iCentrality][iAsymmetry]);
        }
        
        // Add systematic errors to the main histogram to evaluate the systematic error of the integral
        oneBin = jetShapeArray[iCentrality][iTrackPt][iAsymmetry]->FindBin(0.99);
        for(int iBin = 1; iBin <= oneBin; iBin++){
          
          // Find the bin content and systematic error
          binContent = jetShapeForIntegral[iCentrality][iTrackPt][iAsymmetry]->GetBinContent(iBin);
          systematicError = systematicErrorHistogram->GetBinError(iBin);
          
          // Add and subtract the systematic error to/from the distribution
          resultPlusSystematicError->SetBinContent(iBin, binContent+systematicError);
          resultMinusSystematicError->SetBinContent(iBin, binContent-systematicError);
          
        } // Loop over bins of the histogram
        
        // Calculate the integral with statistical errors errors
        integralAfterNormalization[iCentrality][iTrackPt][iAsymmetry] = jetShapeArray[iCentrality][iTrackPt][iAsymmetry]->IntegralAndError(1, oneBin, statisticalError, "width");
        
        // Calculate the integrals for distributions shifted by the systematic uncertainty
        systematicPlusIntegral = resultPlusSystematicError->Integral(1, oneBin, "width");
        systematicMinusIntegral = resultMinusSystematicError->Integral(1, oneBin, "width");
        
        // Systematic uncertainty of the integral is the bigger deviation from the nominal
        systematicError = TMath::Max(systematicPlusIntegral-integralAfterNormalization[iCentrality][iTrackPt][iAsymmetry], integralAfterNormalization[iCentrality][iTrackPt][iAsymmetry]-systematicMinusIntegral);
        
        // Combine the statistical and systematic uncertainty for the total error
        totalError = TMath::Sqrt(systematicError*systematicError+statisticalError*statisticalError);
        
        // Set the total error as the error of the integral
        integralErrorAfterNormalization[iCentrality][iTrackPt][iAsymmetry] = totalError;
        
      } // Track pT loop
    } // Centrality loop
  } // Asymmetry loop
    
  TH1D* ratioHistogram[nCentralityBins][nAsymmetryBins+1];
  TH1D* extraRatio[nCentralityBins];   // Extra ratio histograms to illustrate third jet effects
  
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
      
      jetShapeStack[iCentrality][iAsymmetry] = new stackHist(Form("st_%d_%d",iCentrality, iAsymmetry));
      jetShapeStack[iCentrality][iAsymmetry]->setRange(0.0, 0.99, "x");
      if(normalizeJetShape){
        jetShapeStack[iCentrality][iAsymmetry]->setRange(0.005, 30, "y");
      } else {
        jetShapeStack[iCentrality][iAsymmetry]->setRange(0.5, 3000, "y");
      }
      for(int iTrackPt = iFirstTrackPtBin; iTrackPt < nTrackPtBins ; iTrackPt++){
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
  
  // Create the extra ratio histograms
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    extraRatio[iCentrality] = (TH1D*) sumHistogramPbPb[iCentrality][0]->Clone(Form("extraRatio%d", iCentrality)); // For PbPb use 0.0 < xj < 0.6 bin
    extraRatio[iCentrality]->Divide(sumHistogramPp[3]); // For pp, use the xj integrated bin
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
  
  TString xjbin[] = {"0.0 < x_{j} < 0.6", "0.6 < x_{j} < 0.8", "0.8 < x_{j} < 1.0", "All dijets"};
  
  for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins+1; iAsymmetry++){
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      if( iAsymmetry == 3) {
        bigCanvas->CD(5-iCentrality);
      }else
        bigCanvas->CD(10+iAsymmetry*5-iCentrality);
      gPad->SetLogy();
      jetShapeStack[iCentrality][iAsymmetry]->drawStack();
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetNdivisions(505);
      jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetLabelSize(0.11);
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
        jetShapeStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetLabelSize(0.08);
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
  
  // Original
  ptLegend1->AddEntry(jetShapeStack[4][0]->hist_trunk.at(0), "0.7 < p_{T}^{ch}< 1 GeV","f");
  ptLegend1->AddEntry(jetShapeStack[4][0]->hist_trunk.at(1), "1 < p_{T}^{ch}< 2 GeV","f");
  
  ptLegend2->AddEntry(jetShapeStack[4][0]->hist_trunk.at(2), "2 < p_{T}^{ch}< 3 GeV","f");
  ptLegend2->AddEntry(jetShapeStack[4][0]->hist_trunk.at(3), "3 < p_{T}^{ch}< 4 GeV","f");
  
  ptLegend3->AddEntry(jetShapeStack[4][0]->hist_trunk.at(4), "4 < p_{T}^{ch}< 8 GeV","f");
  ptLegend3->AddEntry(jetShapeStack[4][0]->hist_trunk.at(5), "8 < p_{T}^{ch}< 12 GeV","f");
  
  ptLegend4->AddEntry(jetShapeStack[4][0]->hist_trunk.at(6), "12 < p_{T}^{ch}< 300 GeV","f");
  ptLegend4->AddEntry(sumUncertainty[4][0], "0.7 < p_{T}^{ch}< 300 GeV","lpfe");
  
  lt4->AddEntry(ratioUncertainty[0][0], "0.7 < p_{T}^{ch}< 300 GeV","lpfe");
  
  // No low pT
  /*ptLegend1->AddEntry(jetShapeStack[4][0]->hist_trunk.at(0), "1 < p_{T}^{trk}< 2 GeV","f");
  ptLegend1->AddEntry(jetShapeStack[4][0]->hist_trunk.at(1), "2 < p_{T}^{trk}< 3 GeV","f");
  
  ptLegend2->AddEntry(jetShapeStack[4][0]->hist_trunk.at(2), "3 < p_{T}^{trk}< 4 GeV","f");
  ptLegend2->AddEntry(jetShapeStack[4][0]->hist_trunk.at(3), "4 < p_{T}^{trk}< 8 GeV","f");
  
  ptLegend3->AddEntry(jetShapeStack[4][0]->hist_trunk.at(4), "8 < p_{T}^{trk}< 12 GeV","f");
  ptLegend3->AddEntry(jetShapeStack[4][0]->hist_trunk.at(5), "12 < p_{T}^{trk}< 300 GeV","f");
  
  ptLegend4->AddEntry(sumUncertainty[4][0], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
  
  //lt4->AddEntry(ratioUncertainty[0][0], "1 < p_{T}^{trk}< 300 GeV","lpfe"); */
  
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
  
  double cmsPositionX = 0.28;
  double cmsPositionY = 0.97;
  double jetShapeTitlePosition[] = {0.4,0.4,0.4};
  double xjPosition[] = {0.76,0.77,0.76};
  double systemPosition = 0.26;
  double selectionPosition = 0.205;
  double cmsSize = 0.035;
  
  if(addPreliminaryTag){
    cmsPositionX = 0.065;
    cmsPositionY = 0.96;
    jetShapeTitlePosition[0] = 0.36;
    jetShapeTitlePosition[1] = 0.34;
    jetShapeTitlePosition[2] = 0.36;
  }
  
  // If we draw radial momentum distribution instead of jet shape, need to change and reposition the labels
  if(!normalizeJetShape){
    cmsPositionX = 0.33;
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
  mainTitle->SetTextSize(0.04);
  mainTitle->DrawLatexNDC(cmsPositionX, cmsPositionY, "CMS");
  
  if(addPreliminaryTag){
    mainTitle->SetTextFont(42);
    mainTitle->SetTextSize(0.035);
    mainTitle->DrawLatexNDC(0.03, 0.93, "Preliminary");
  }
  
  mainTitle->SetTextFont(42);
  mainTitle->SetTextSize(0.035);
  mainTitle->DrawLatexNDC(jetShapeTitlePosition[iJetTrack/3], 0.97, jetShapeTitle[iJetTrack/3]);
  
  //mainTitle->SetTextFont(42);
  //mainTitle->SetTextSize(0.035);
  //mainTitle->DrawLatexNDC(xjPosition[iJetTrack/3], 0.94, xjString[iAsymmetry]);
  
  mainTitle->SetTextSize(0.03);
  if(monteCarloLabels){
    mainTitle->DrawLatexNDC(systemPosition, 0.94, "5.02 TeV   Pythia8   Pythia+Hydjet");
  } else {
    mainTitle->DrawLatexNDC(systemPosition, 0.94, "5.02 TeV   pp 320 pb^{-1}   PbPb 1.7 nb^{-1}");
  }

  mainTitle->SetTextSize(0.022);
  mainTitle->DrawLatexNDC(selectionPosition, 0.915, "anti-k_{T} R = 0.4, |#eta_{jet}| < 1.6, p_{T,1} > 120 GeV, p_{T,2} > 50 GeV, #Delta#varphi_{1,2} > #frac{5#pi}{6}");
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
  
  //box->DrawBox(0.23,0.3067, 0.243, 0.315);
  //box->DrawBox(0.23,0.5735, 0.243, 0.58);
  
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
        ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetLabelOffset(0.02);
        ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetLabelSize(0.08);
      }else {
        ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetTitleSize(0.1);
        ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetLabelOffset(0.015);
        ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetLabelSize(0.1);
      }
      
      ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetTitle("#rho(#Deltar)_{PbPb}/#rho(#Deltar)_{pp}");
      ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->CenterTitle();
      ratioUncertainty[iCentrality][iAsymmetry]->SetStats(0);
      ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->SetNdivisions(505);
      ratioUncertainty[iCentrality][iAsymmetry]->GetYaxis()->CenterTitle();
      ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->SetTitle("#Deltar");
      ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->CenterTitle();
      ratioUncertainty[iCentrality][iAsymmetry]->SetAxisRange(0.01, 3.1 - (iJetTrack/3)*0.5, "Y");
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
      
      ratioUncertainty[iCentrality][iAsymmetry]->SetTitle("");
      ratioUncertainty[iCentrality][iAsymmetry]->Draw("same e2");
      ratioUncertainty[iCentrality][nAsymmetryBins]->Draw("same e2");
      ratioHistogram[iCentrality][iAsymmetry]->Draw("same");
      ratioHistogram[iCentrality][nAsymmetryBins]->Draw("same");
      auxi_hist[iCentrality]->Draw("same e2");
      if(iAsymmetry == 0 && iCentrality != 3 && drawExtraRatio){
        extraRatio[iCentrality]->SetLineColor(kRed+2);
        extraRatio[iCentrality]->SetLineStyle(2);
        extraRatio[iCentrality]->SetLineWidth(2);
        extraRatio[iCentrality]->Draw("l,hist,same");
      }
      
      line[iAsymmetry] = new TLine();
      line[iAsymmetry]->SetLineStyle(2);
      line[iAsymmetry]->DrawLine(0, 1, 1, 1);
      
    }
  }
  //canvas2caption
  ratioCanvas->cd(0);
  
  double cmsPositionXratio = 0.05;
  double cmsPositionYratio = 0.95;
  double jetShapeRatioTitlePosition[] = {0.14, 0.14, 0.14};
  
  if(addPreliminaryTag){
    cmsPositionXratio = 0.04;
    cmsPositionYratio = 0.97;
    jetShapeRatioTitlePosition[0] = 0.18;
    jetShapeRatioTitlePosition[1] = 0.158;
    jetShapeRatioTitlePosition[2] = 0.18;
  }
  
  mainTitle->SetTextFont(62);
  mainTitle->SetTextSize(0.04);
  mainTitle->DrawLatexNDC(cmsPositionXratio, cmsPositionYratio, "CMS");
  
  if(addPreliminaryTag){
    mainTitle->SetTextFont(42);
    mainTitle->SetTextSize(0.035);
    mainTitle->DrawLatexNDC(0.01, 0.94, "Preliminary");
  }
  
  mainTitle->SetTextFont(42);
  mainTitle->SetTextSize(0.035);
  mainTitle->DrawLatexNDC(jetShapeRatioTitlePosition[iJetTrack/3], 0.95, Form("%s ratios",jetShapeTitle[iJetTrack/3]));
  
  mainTitle->SetTextSize(0.024);
  if(monteCarloLabels){
    mainTitle->DrawLatexNDC(0.58, 0.975, "5.02 TeV   Pythia8   Pythia+Hydjet");
  } else {
    mainTitle->DrawLatexNDC(0.58, 0.975, "5.02 TeV   pp 320 pb^{-1}   PbPb 1.7 nb^{-1}");
  }
  
  mainTitle->SetTextSize(0.02);
  mainTitle->DrawLatexNDC(0.505, 0.945, "anti-k_{T} R = 0.4, |#eta_{jet}| < 1.6, p_{T,1} > 120 GeV, p_{T,2} > 50 GeV, #Delta#varphi_{1,2} > #frac{5#pi}{6}");
  
  mainTitle->SetTextFont(42);
  mainTitle->SetTextSize(0.03);
  box->SetFillColor(kWhite);
  float boxwidth = 0.02, boxhight = 0.03;
  float x= 0.28, y=0.058, offsetx = 0.005, offsety = 0.005;
  box->DrawBox(x, y, x+boxwidth, y+boxhight); mainTitle->DrawLatex(x+offsetx, y+offsety, "0");
  //		box->SetFillColor(kAzure+6);
  x= 0.51, y=0.058; box->DrawBox(x, y, x+boxwidth, y+boxhight); mainTitle->DrawLatex(x+offsetx, y+offsety, "0");
  x= 0.746,y=0.058; box->DrawBox(x, y, x+boxwidth, y+boxhight); mainTitle->DrawLatex(x+offsetx, y+offsety, "0");
  x= 0.973, y=0.058; box->DrawBox(x, y, x+boxwidth, y+boxhight); mainTitle->DrawLatex(x+offsetx, y+offsety, "1");
  
  TLegend *ltc2 = new TLegend(0.27, 0.6, 0.72, 0.8);
  ltc2->SetLineColor(0);
  ltc2->SetTextSize(0.085);
  ltc2->AddEntry(ratioUncertainty[0][3], xjbin[3], "lpf");
  ltc2->AddEntry(ratioUncertainty[0][0], xjbin[0], "lpf");

  TLegend *ltc3 = new TLegend(0.27, 0.75, 0.72, 0.95);
  ltc3->SetLineColor(0);
  ltc3->SetTextSize(0.085);
  ltc3->AddEntry(ratioUncertainty[0][3], xjbin[3], "lpf");
  ltc3->AddEntry(ratioUncertainty[0][1], xjbin[1], "lpf");
  
  TLegend *ltc4 = new TLegend(0.27, 0.78, 0.72, 0.95);
  ltc4->SetLineColor(0);
  ltc4->SetTextSize(0.07);
  ltc4->AddEntry(ratioUncertainty[0][3], xjbin[3], "lpf");
  ltc4->AddEntry(ratioUncertainty[0][2], xjbin[2], "lpf");
  
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
  ratioCanvas->cd(5);
  ltc3->Draw();
  ratioCanvas->cd(9);
  ltc4->Draw();
  
  // Legend for extra ratio
  if(drawExtraRatio){
    ratioCanvas->cd(2);
    TLegend *extraLegend = new TLegend(0.05, 0.68, 0.72, 0.84);
    extraLegend->SetFillStyle(0);extraLegend->SetBorderSize(0);extraLegend->SetTextSize(0.07);extraLegend->SetTextFont(62);
    extraLegend->AddEntry(extraRatio[0], "PbPb(x_{j} < 0.6) / pp (all x_{j})","l");
    extraLegend->Draw();
  } // Drawing extra ratio for illustration
  
  bigCanvas->SaveAs(Form("figures/final%s_%s_preliminaryTag.pdf", saveString, jetShapeSaveName[iJetTrack/3]));
  ratioCanvas->SaveAs(Form("figures/final%s_%s_ratioPreliminaryTag.pdf", saveString, jetShapeSaveName[iJetTrack/3]));
  
  // Save the histograms to a file for HepData submission
  if(saveHistogramsForHepData){
    
    TString outputFileName = "hepdata/hepdata_jetShapes_hin-19-013.root";
    if(!normalizeJetShape) outputFileName = "hepdata/hepdata_jetRadialMomentum_hin-19-013.root";
    TString shapeName = "Shape";
    if(!normalizeJetShape) shapeName = "RadialMomentum";
    TFile *outputFile = TFile::Open(outputFileName,"UPDATE");
    TString centralityString[] = {"0-10", "10-30", "30-50", "50-90", "pp"};
    TString asymmetryString[] = {"0<xj<06","06<xj<08","08<xj<1","allXj"};
    
    int firstSaveBin = 0;
    if(!normalizeJetShape) firstSaveBin = nAsymmetryBins;
    
    for(int iAsymmetry = firstSaveBin; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        
        sumHistogramPbPb[iCentrality][iAsymmetry]->GetXaxis()->SetRangeUser(0,1);
        sumHistogramPbPb[iCentrality][iAsymmetry]->Write(Form("jet%s_%s_%s_%s", shapeName.Data(), jetShapeSaveName[iJetTrack/3], centralityString[iCentrality].Data(), asymmetryString[iAsymmetry].Data()), TObject::kOverwrite);;
        sumUncertainty[iCentrality][iAsymmetry]->GetXaxis()->SetRangeUser(0,1);
        sumUncertainty[iCentrality][iAsymmetry]->Write(Form("jet%sError_%s_%s_%s", shapeName.Data(), jetShapeSaveName[iJetTrack/3], centralityString[iCentrality].Data(), asymmetryString[iAsymmetry].Data()), TObject::kOverwrite);
        ratioHistogram[iCentrality][iAsymmetry]->GetXaxis()->SetRangeUser(0,1);
        ratioHistogram[iCentrality][iAsymmetry]->Write(Form("jet%sRatio_%s_%s_%s", shapeName.Data(), jetShapeSaveName[iJetTrack/3], centralityString[iCentrality].Data(), asymmetryString[iAsymmetry].Data()), TObject::kOverwrite);
        ratioUncertainty[iCentrality][iAsymmetry]->GetXaxis()->SetRangeUser(0,1);
        ratioUncertainty[iCentrality][iAsymmetry]->Write(Form("jet%sRatioError_%s_%s_%s", shapeName.Data(), jetShapeSaveName[iJetTrack/3], centralityString[iCentrality].Data(), asymmetryString[iAsymmetry].Data()), TObject::kOverwrite);
      } // Centrality loop
      
      // No centrality binning for pp
      
      sumHistogramPp[iAsymmetry]->GetXaxis()->SetRangeUser(0,1);
      sumHistogramPp[iAsymmetry]->Write(Form("jet%s_%s_pp_%s", shapeName.Data(), jetShapeSaveName[iJetTrack/3], asymmetryString[iAsymmetry].Data()), TObject::kOverwrite);
      sumUncertainty[nCentralityBins][iAsymmetry]->GetXaxis()->SetRangeUser(0,1);
      sumUncertainty[nCentralityBins][iAsymmetry]->Write(Form("jet%sError_%s_pp_%s", shapeName.Data(), jetShapeSaveName[iJetTrack/3], asymmetryString[iAsymmetry].Data()), TObject::kOverwrite);
      
    } // Asymmetry loop
    
    outputFile->Close();
    
  } // Save HepData
  
  // Print the integrals for different bins to the console
  if(printIntegrals){
    char namer[200];
    double ptBinBorders[] = {0.7, 1, 2, 3, 4, 8, 12, 300};
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      
      sprintf(namer,"\\frametitle{%s integral, $%s$}",jetShapeTitle[iJetTrack/3], xjbin[iAsymmetry].Data());
      
      cout << endl;
      cout << "\\begin{frame}" << endl;
      cout << namer << endl;
      cout << "\\begin{center}" << endl;
      cout << "  \\begin{tabular}{cccccc}" << endl;
      cout << "    \\toprule" << endl;
      cout << "     p_{\\mathrm{T}} (GeV) & C: 0-10 & C: 10-30 & C: 30-50 & C: 50-90 & pp \\\\" << endl;
      cout << "    \\midrule" << endl;
      
      // Set the correct precision for printing floating point numbers
      cout << fixed << setprecision(3);
      
      for(int iTrackPt = iFirstTrackPtBin; iTrackPt < nTrackPtBins; iTrackPt++){
        
        sprintf(namer, "$%.1f < p_{\\mathrm{T}} < %.1f$", ptBinBorders[iTrackPt], ptBinBorders[iTrackPt+1]);
        cout << namer;
        for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
          cout << " & $" << integralAfterNormalization[iCentrality][iTrackPt][iAsymmetry] << " \\pm " <<  integralErrorAfterNormalization[iCentrality][iTrackPt][iAsymmetry] << "$";
        }
        cout << " \\\\" << endl;
      }
      
      cout << "    \\bottomrule" << endl;
      cout << "  \\end{tabular}" << endl;
      cout << "\\end{center}" << endl;
      cout << "\\end{frame}" << endl;
      cout << endl;
      
    }
  }
  
}

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void finalBigAsymmetryPlotter_v1(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Open data files for pp and PbPb data
  TFile *ppFile = TFile::Open("data/ppData2017_highForest_pfJets_20EveMixed_xjBins_wtaAxis_allCorrections_processed_2020-02-04.root");
  // data/ppData2017_highForest_pfJets_20EveMixed_xjBins_wtaAxis_allCorrections_processed_2020-02-04.root <--- Nominal result file
  TFile *pbpbFile = TFile::Open("data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_allCorrections_finalTuning_onlyJetShapes_processed_2020-02-17.root");
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_allCorrections_finalTuning_onlyJetShapes_processed_2020-02-17.root <-- File with manual fluctuation reduction
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_onlyJetShapa_manualTuning_wtaAxis_processed_2020-02-04.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_tuningForSeagull_wtaAxis_processed_2020-02-04.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrectionsExceptSpillover_wtaAxis_processed_2020-02-06.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_oldSpillover_finalTuning_wtaAxis_processed_2020-02-04.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_scaledSpillover_finalTuning_wtaAxis_processed_2020-02-04.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_finalTuning_wtaAxis_processed_2020-01-29.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrections_seagullTuningProcess_processed_2020-01-15.root
  // data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_allCorrectionsWithShifterCentrality_trackDeltaRonlyInLowPt_processed_2019-10-17.root
  
  TFile *ppUncertaintyFile = TFile::Open("uncertainties/systematicUncertaintyForPp_20eveMix_xjBins_jffAndJetResolutionUpdate_2020-06-03.root");
  // uncertainties/systematicUncertaintyForPp_20eveMix_xjBins_jffAndJetResolutionUpdate_2020-06-03.root  <-- Final final file
  // uncertainties/systematicUncertaintyForPp_20eveMix_xjBins_fixJES_2020-02-03.root
  TFile *pbpbUncertaintyFile = TFile::Open("uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_manualFluctuationReduction_addTriggerBias_jffUpdate_2020-06-03.root");
  // uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_manualFluctuationReduction_addTriggerBias_jffUpdate_2020-06-03.root <-- Final final file
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
