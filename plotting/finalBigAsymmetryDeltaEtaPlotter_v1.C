#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "stackHist.h"
#include "xCanvas.h"
#include "JDrawer.h"
#include "DijetMethods.h"

void plotDeltaEtaBigAsymmetry(DijetHistogramManager *ppHistograms, DijetHistogramManager *pbpbHistograms, JffCorrector *ppUncertaintyProvider, JffCorrector *pbpbUncertaintyProvider, int iJetTrack){
  
  // Determine the number of bins from the histogram manager
  const int nTrackPtBins = pbpbHistograms->GetNTrackPtBins();
  const int nCentralityBins = pbpbHistograms->GetNCentralityBins();
  const int nAsymmetryBins = pbpbHistograms->GetNAsymmetryBins();
  
  // Define labels and plotting layout
  const char* deltaEtaTitle[] = {"Particle yields associated with leading jets","Particle yields associated with subleading jets","Particle yield wrt. inclusive jet"};
  //const char* deltaEtaTitle[] = {"Particle yield w.r.t. the leading jet","Particle yield w.r.t. the subleading jet","Particle yield wrt. inclusive jet"};
  const char* deltaEtaSaveName[] = {"trackLeadingJet","trackSubleadingJet","trackInclusiveJet"};
  const char* xjString[] = {"0.0 < x_{j} < 0.6","0.6 < x_{j} < 0.8","0.8 < x_{j} < 1.0","x_{j} inclusive"};
  const char* asymmetrySaveName[] = {"_A=0v0-0v6","_A=0v6-0v8","_A=0v8-1v0",""};
  double deltaEtaZoom[] = {37,37,30};
  double subtractionZoom[] = {23,23,15};
  
  const bool drawInclusive = false;
  
  // Open a file to include results from inclusive jet analysis HIN-16-020
  //  TFile *inclusiveResultFile = TFile::Open("data/publishedResults/officialHist_py_deta_16_020.root");
  
  // Methods for doing stuff
  DijetMethods *projector = new DijetMethods();
  double projectionRegion = 1;  // Region in deltaPhi which is used to project the deltaEta peak
  const int nDeltaEtaBinsRebin = 21;
  const double deltaEtaBinBordersRebin[nDeltaEtaBinsRebin+1] = {-4,-3,-2,-1.5,-1,-0.8,-0.6,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.6,0.8,1,1.5,2,3,4};
  
  // Histograms needed to calculate the stacked deltaEta distributions
  TH1D *sumHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *sumHistogramPp[nAsymmetryBins+1];
  TH1D *deltaEtaArray[nCentralityBins+1][nTrackPtBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPbPb[nCentralityBins][nAsymmetryBins+1];
  TH1D *uncertaintyHistogramPp[nAsymmetryBins+1];
  TH1D *sumUncertainty[nCentralityBins+1][nAsymmetryBins+1];
  TH2D *helperHistogram;
  TH1D *addedHistogram;
  TH1D *rebinnedHistogram;
  double shapeIntegralPbPb[nCentralityBins][nAsymmetryBins+1];
  double shapeIntegralPp[nAsymmetryBins+1];
  double binContent;
  
  // Histgrams for comparison with inclusive results
  TH1D *sumHistogramInclusive[nCentralityBins+1];
  
  // Read two dimensional deltaEta-deltaPhi histograms and project the deltaEta yields out from them
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      helperHistogram = pbpbHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, 0);
      
      deltaEtaArray[iCentrality][0][iAsymmetry] = (TH1D*) projector->ProjectAnalysisYieldDeltaEta(helperHistogram, pbpbHistograms->GetTrackPtBinBorder(0), pbpbHistograms->GetTrackPtBinBorder(1), true)->Clone(Form("deltaEtaArray%d%d%d", iCentrality, iAsymmetry, 0));
      
      sumHistogramPbPb[iCentrality][iAsymmetry] = (TH1D*) deltaEtaArray[iCentrality][0][iAsymmetry]->Clone(Form("sumPbPb%d%d",iCentrality,iAsymmetry));
      
      /*
       addedHistogram = projector->ProjectRegionDeltaEta(helperHistogram, -projectionRegion, projectionRegion, Form("DeltaEtaStack%d%d%d%d",iJetTrack,iAsymmetry,iCentrality,0));
       
       // Since we want to plot yield, we do not want to normalize over the number of bins projected over
       // but simply look at the yield in certain region
       addedHistogram->Scale(projector->GetNBinsProjectedOver());
       
       // The different pT bins can have different size, so we need to normalize the yield with the pT bin width
       addedHistogram->Scale(1/(pbpbHistograms->GetTrackPtBinBorder(1) - pbpbHistograms->GetTrackPtBinBorder(0)));
       
       // Rebin the histogram to match the binning in the inclusive jet shape paper
       sumHistogramPbPb[iCentrality][iAsymmetry] = (TH1D*) projector->RebinAsymmetric(addedHistogram, nDeltaEtaBinsRebin, deltaEtaBinBordersRebin)->Clone(Form("sumPbPb%d%d",iCentrality,iAsymmetry));
       
       deltaEtaArray[iCentrality][0][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("deltaEtaArray%d%d%d", iCentrality, iAsymmetry, 0));*/
    } // Centrality loop
    
    helperHistogram = ppHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, 0, 0);
    
    deltaEtaArray[nCentralityBins][0][iAsymmetry] = (TH1D*) projector->ProjectAnalysisYieldDeltaEta(helperHistogram, ppHistograms->GetTrackPtBinBorder(0), ppHistograms->GetTrackPtBinBorder(1), true)->Clone(Form("deltaEtaArray%d%d%d", nCentralityBins, iAsymmetry, 0));
    
    sumHistogramPp[iAsymmetry] = (TH1D*) deltaEtaArray[nCentralityBins][0][iAsymmetry]->Clone(Form("sumPp%d",iAsymmetry));
    
    /*addedHistogram = projector->ProjectRegionDeltaEta(helperHistogram, -projectionRegion, projectionRegion, Form("DeltaEtaStack%d%d%d%d",iJetTrack,iAsymmetry,nCentralityBins,0));
     
     // Since we want to plot yield, we do not want to normalize over the number of bins projected over
     // but simply look at the yield in certain region
     addedHistogram->Scale(projector->GetNBinsProjectedOver());
     
     // The different pT bins can have different size, so we need to normalize the yield with the pT bin width
     addedHistogram->Scale(1/(ppHistograms->GetTrackPtBinBorder(1) - ppHistograms->GetTrackPtBinBorder(0)));
     
     // Rebin the histogram to match the binning in the inclusive jet shape paper
     sumHistogramPp[iAsymmetry] = (TH1D*) projector->RebinAsymmetric(addedHistogram, nDeltaEtaBinsRebin, deltaEtaBinBordersRebin)->Clone(Form("sumPp%d",iAsymmetry));
     
     deltaEtaArray[nCentralityBins][0][iAsymmetry] = (TH1D*) sumHistogramPp[iAsymmetry]->Clone(Form("deltaEtaArray%d%d%d", nCentralityBins, iAsymmetry, 0));*/
    
  } // Asymmetry loop
  
  // Sum the pT:s
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iTrackPt = 1; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        
        helperHistogram = pbpbHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, iTrackPt);
        
        deltaEtaArray[iCentrality][iTrackPt][iAsymmetry] = (TH1D*) projector->ProjectAnalysisYieldDeltaEta(helperHistogram, pbpbHistograms->GetTrackPtBinBorder(iTrackPt), pbpbHistograms->GetTrackPtBinBorder(iTrackPt+1), true)->Clone(Form("deltaEtaArray%d%d%d", iCentrality, iAsymmetry, iTrackPt));
        
        sumHistogramPbPb[iCentrality][iAsymmetry]->Add(deltaEtaArray[iCentrality][iTrackPt][iAsymmetry]);
        
      } // Centrality loop
      
      helperHistogram = ppHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, 0, iTrackPt);
      
      deltaEtaArray[nCentralityBins][iTrackPt][iAsymmetry] = (TH1D*) projector->ProjectAnalysisYieldDeltaEta(helperHistogram, ppHistograms->GetTrackPtBinBorder(iTrackPt), ppHistograms->GetTrackPtBinBorder(iTrackPt+1), true)->Clone(Form("deltaEtaArray%d%d%d", nCentralityBins, iAsymmetry, iTrackPt));
      
      sumHistogramPp[iAsymmetry]->Add(deltaEtaArray[nCentralityBins][iTrackPt][iAsymmetry]);
      
    } // Track pT loop
  } // Asymmetry loop
  
  // Read the inclusive histograms to make a comparison to the published results
  int inclusiveCentralityBins[] = {0,10,30,50,100};
  int inclusiveTrackPtBins[] = {0,1,2,3,4,8,12};
  
  // Find the uncertainties for the deltaEta histograms
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    
    // First, read the uncertainties in the lowest pT bin to the uncertainty histograms
    uncertaintyHistogramPp[iAsymmetry] = ppUncertaintyProvider->GetDeltaEtaSystematicUncertainty(iJetTrack, 0, 0, iAsymmetry);
    uncertaintyHistogramPp[iAsymmetry]->Multiply(uncertaintyHistogramPp[iAsymmetry]);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      uncertaintyHistogramPbPb[iCentrality][iAsymmetry] = pbpbUncertaintyProvider->GetDeltaEtaSystematicUncertainty(iJetTrack, iCentrality, 0, iAsymmetry);
      uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->Multiply(uncertaintyHistogramPbPb[iCentrality][iAsymmetry]);
    }
    
    // For the pT summed uncertainties, add the uncertainties in each pT bin in quadrature
    for(int iTrackPt = 1; iTrackPt < 6; iTrackPt++){
      addedHistogram = ppUncertaintyProvider->GetDeltaEtaSystematicUncertainty(iJetTrack, 0, iTrackPt, iAsymmetry);
      addedHistogram->Multiply(addedHistogram);
      uncertaintyHistogramPp[iAsymmetry]->Add(addedHistogram);
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        addedHistogram = pbpbUncertaintyProvider->GetDeltaEtaSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry);
        addedHistogram->Multiply(addedHistogram);
        uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->Add(addedHistogram);
      }
    }
    
    // In the end, take the square root bin by bin
    for(int iBin = 1; iBin <= uncertaintyHistogramPp[iAsymmetry]->GetNbinsX(); iBin++){
      binContent = uncertaintyHistogramPp[iAsymmetry]->GetBinContent(iBin);
      uncertaintyHistogramPp[iAsymmetry]->SetBinContent(iBin, TMath::Sqrt(binContent));
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        binContent = uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->GetBinContent(iBin);
        uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->SetBinContent(iBin, TMath::Sqrt(binContent));
      }
    }
    
    // Create the uncertainties for summed histograms by first cloning the summed histgorams
    sumUncertainty[nCentralityBins][iAsymmetry] = (TH1D*) sumHistogramPp[iAsymmetry]->Clone(Form("summedUncertaintyPpA%d",iAsymmetry));
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      sumUncertainty[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("summedUncertaintyPbPbC%dA%d", iCentrality, iAsymmetry));
    }
    
    // Set the errors for the summed uncertainty histograms from the extracted and summed error from the file
    for(int iBin = 1; iBin <= sumUncertainty[nCentralityBins][iAsymmetry]->GetNbinsX(); iBin++){
      sumUncertainty[nCentralityBins][iAsymmetry]->SetBinError(iBin,uncertaintyHistogramPp[iAsymmetry]->GetBinContent(iBin));
    }
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iBin = 1; iBin <= sumUncertainty[iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
        sumUncertainty[iCentrality][iAsymmetry]->SetBinError(iBin,uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->GetBinContent(iBin));
      }
    }
  }
  
  TH1D* subtractionHistogram[nCentralityBins][nTrackPtBins][nAsymmetryBins+1];
  
  // Create uncertainty histograms for the ratio using standard jet shape binning
  TH1D* subtractionUncertainty[nCentralityBins][nAsymmetryBins+1];
  double ppUncertainty, pbpbUncertainty, subtractionUncertaintyValue;
  
  // Subtract pp from PbPb in each bin and find the uncertainty for the total subtraction
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      // Subtract pp from PbPb in each track pT bin
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        subtractionHistogram[iCentrality][iTrackPt][iAsymmetry] = (TH1D*) deltaEtaArray[iCentrality][iTrackPt][iAsymmetry]->Clone(Form("subtraction%d%d%d", iCentrality, iTrackPt, iAsymmetry));
        subtractionHistogram[iCentrality][iTrackPt][iAsymmetry]->Add(deltaEtaArray[nCentralityBins][iTrackPt][iAsymmetry],-1);
      } // Track pt loop
      
      subtractionUncertainty[iCentrality][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("subtractionUncertainty%d%d", iCentrality, iAsymmetry));
      subtractionUncertainty[iCentrality][iAsymmetry]->Add(sumHistogramPp[iAsymmetry],-1);
      
      // Calculate the systematic uncertainty for the ratio
      for(int iBin = 1; iBin < subtractionUncertainty[iCentrality][iAsymmetry]->GetNbinsX(); iBin++){
        ppUncertainty = uncertaintyHistogramPp[iAsymmetry]->GetBinContent(iBin);
        pbpbUncertainty = uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->GetBinContent(iBin);
        subtractionUncertaintyValue = TMath::Sqrt(TMath::Power(pbpbUncertainty,2)+TMath::Power(ppUncertainty,2));
        subtractionUncertainty[iCentrality][iAsymmetry]->SetBinError(iBin,subtractionUncertaintyValue);
      } // Bin loop
    } // Centrality loop
  } // Asymmetry loop
  
  TString cent_lab[4] = {"0-10%", "10-30%", "30-50%", "50-90%"};
  stackHist *deltaEtaStack[nCentralityBins+1][nAsymmetryBins+1];
  stackHist *subtractedStack[nCentralityBins][nAsymmetryBins+1];
  
  // Stack the delta eta histograms
  for(int iCentrality = 0; iCentrality < nCentralityBins+1; iCentrality++){
    for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
      sumUncertainty[iCentrality][iAsymmetry]->SetFillStyle(1001);
      sumUncertainty[iCentrality][iAsymmetry]->SetFillColorAlpha(kGray+3, 0.4);
      sumUncertainty[iCentrality][iAsymmetry]->SetMarkerStyle(20);
      //	sumUncertainty[iCentrality][iAsymmetry]->SetMarkerSize(.0);
      sumUncertainty[iCentrality][iAsymmetry]->SetMarkerSize(1.4);
      sumUncertainty[iCentrality][iAsymmetry]->SetMarkerColor(0);
      sumUncertainty[iCentrality][iAsymmetry]->SetLineColor(kBlack);
      
      deltaEtaStack[iCentrality][iAsymmetry] = new stackHist(Form("st_%d",iCentrality));
      deltaEtaStack[iCentrality][iAsymmetry]->setRange(-1.5, 1.5, "x");
      deltaEtaStack[iCentrality][iAsymmetry]->setRange(-4, deltaEtaZoom[iJetTrack/3], "y");
      for(int iTrackPt = nTrackPtBins-2; iTrackPt >= 0; iTrackPt--){
        deltaEtaStack[iCentrality][iAsymmetry]->addHist((TH1*) deltaEtaArray[iCentrality][iTrackPt][iAsymmetry]);
      }
      //jetShapeSum[iCentrality]->SetMarkerColor(0);
    } // Asymmetry loop
  } // Centrality loop
  
  // Stack the subtracted (PbPb - pp) deltaEta histograms
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins+1; iAsymmetry++){
      
      subtractedStack[iCentrality][iAsymmetry] = new stackHist(Form("st_%d",iCentrality));
      subtractedStack[iCentrality][iAsymmetry]->setRange(-1.5, 1.5, "x");
      subtractedStack[iCentrality][iAsymmetry]->setRange(-1, subtractionZoom[iJetTrack/3], "y");
      
      for(int iTrackPt = nTrackPtBins-2; iTrackPt >= 0; iTrackPt--){
        subtractedStack[iCentrality][iAsymmetry]->addHist((TH1*) subtractionHistogram[iCentrality][iTrackPt][iAsymmetry]);
      }
      
      subtractionUncertainty[iCentrality][iAsymmetry]->SetFillStyle(1001);
      subtractionUncertainty[iCentrality][iAsymmetry]->SetFillColorAlpha(kGray+3, 0.4);
      subtractionUncertainty[iCentrality][iAsymmetry]->SetMarkerColor(0);
      subtractionUncertainty[iCentrality][iAsymmetry]->SetMarkerStyle(21);
      subtractionUncertainty[iCentrality][iAsymmetry]->SetLineColor(kBlack);
    }
  }
  
  // ==============================
  // **  Draw the distributions  **
  // ==============================
  //canvas2
  //
  TH1D* auxi_hist[5][4];
  for(int j=0 ;j<4; j++){
    for(int i=0 ;i<5; i++){
      //						auxi_hist[i][j] = (TH1D*) sumUncertainty[i][j]->Clone(Form("auxi_hist_%d_%d", i,j));
      // maybe later for state error
      if(i<4){
        auxi_hist[i][j] = (TH1D*) sumHistogramPbPb[i][j]->Clone(Form("auxi_hist_%d_%d", i,j));
      }else {
        auxi_hist[i][j] = (TH1D*) sumHistogramPp[j]->Clone(Form("auxi_hist_%d_%d", i,j));
      }
    }
  }
  
  TH1D *axisDrawer;
  
  // Draw all the distributions to big canvases
  auxi_canvas *bigCanvas = new auxi_canvas(Form("theReallyBigCanvas%d",iJetTrack), "", 2500, 2500);
  bigCanvas->SetMargin(0.08, 0.01, 0.08, 0.2); // Margin order: Left, Right, Bottom, Top
  bigCanvas->divide(4,5);
  
  TBox *box = new TBox();
  TLatex *mainTitle = new TLatex();
  TString xjbin[] = {"0.0 < x_{j} < 0.6", "0.6 < x_{j} < 0.8", "0.8 < x_{j} < 1.0", "all dijets"};
  
  TLine *line[nAsymmetryBins];
  const char *pbpbLabel = "PbPb";
  const char *ppLabel = "pp";
  
  for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins+1; iAsymmetry++){
    
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      
      cout << "In the loop iAsymmetry: " << iAsymmetry << " iCentrality: " << iCentrality << endl;
      
      if(iAsymmetry == 3 ) bigCanvas->CD(5-iCentrality);
      else bigCanvas->CD(10+iAsymmetry*5-iCentrality);
      
      deltaEtaStack[iCentrality][iAsymmetry]->drawStack("","hist",true);
      deltaEtaStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetNdivisions(505);
      deltaEtaStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetLabelSize(0.11);
      deltaEtaStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitleOffset(1.2);
      deltaEtaStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitleSize(0.1);
      deltaEtaStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitle("Y = #frac{1}{N_{dijet}} #frac{dN}{d#Delta#eta}");
      deltaEtaStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitle("#Delta#eta");
      if(iCentrality == nCentralityBins && iAsymmetry == 2){
        //most left-bottom conner
        deltaEtaStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitleSize(0.1);
        deltaEtaStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitleOffset(1);
        deltaEtaStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetLabelSize(0.08);
        deltaEtaStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetLabelOffset(0.03);
        deltaEtaStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetLabelSize(0.08);
        deltaEtaStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetLabelOffset(0.01);
      }else {
        deltaEtaStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitleOffset(0.8);
        deltaEtaStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitleSize(0.15);
        deltaEtaStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitleSize(0.14);
        deltaEtaStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitleOffset(0.71);
        deltaEtaStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetLabelSize(0.11);
        deltaEtaStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetLabelOffset(0.011);
      }
      auxi_hist[iCentrality][iAsymmetry]->SetLineColor(kBlack);
      auxi_hist[iCentrality][iAsymmetry]->SetFillColorAlpha(kWhite, 0);
      auxi_hist[iCentrality][iAsymmetry]->SetMarkerStyle(24);
      auxi_hist[iCentrality][iAsymmetry]->SetMarkerColor(kBlack);
      //auxi_hist[iCentrality][iAsymmetry]->SetMarkerSize(.0);
      auxi_hist[iCentrality][iAsymmetry]->SetMarkerSize(1.7);
      
      deltaEtaStack[iCentrality][iAsymmetry]->hst->Draw();
      auxi_hist[iCentrality][iAsymmetry]->Draw("same");
      
      sumUncertainty[iCentrality][iAsymmetry]->Draw("same e2");
      
      if(iCentrality == nCentralityBins){
        mainTitle->SetTextFont(22);
        if( iAsymmetry == 2)
          mainTitle->SetTextSize(0.08);
        else
          mainTitle->SetTextSize(0.11);
        mainTitle->DrawLatex(-1.15,deltaEtaZoom[0]-5, ppLabel);
        mainTitle->DrawLatex(-0.4,deltaEtaZoom[0]-5, xjbin[iAsymmetry]);
      }
      else {
        mainTitle->SetTextFont(22);
        mainTitle->SetTextSize(0.11);
        mainTitle->DrawLatex(-1.15, deltaEtaZoom[0]-3.15, pbpbLabel);
        mainTitle->DrawLatex(-1.15, deltaEtaZoom[0]-7, cent_lab[iCentrality]);
      }
    } // Centrality loop
  } // Asymmetry loop for drawing distributions to the really big canvas
  cout << "Loopity loop done" << endl;
  TLegend* ptLegend1 = new TLegend(0.04, 0.82, 0.27, 0.89);
  TLegend* ptLegend2 = new TLegend(0.28, 0.82, 0.51, 0.89);
  TLegend* ptLegend3 = new TLegend(0.53 ,0.82, 0.76, 0.89);
  TLegend* ptLegend4 = new TLegend(0.77, 0.855, 1.00,0.89);
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
  
  ptLegend1->AddEntry(deltaEtaStack[4][0]->hist_trunk.at(5), "0.7 < p_{T}^{trk}< 1 GeV","f");
  ptLegend1->AddEntry(deltaEtaStack[4][0]->hist_trunk.at(4), "1 < p_{T}^{trk}< 2 GeV","f");
  
  ptLegend2->AddEntry(deltaEtaStack[4][0]->hist_trunk.at(3), "2 < p_{T}^{trk}< 3 GeV","f");
  ptLegend2->AddEntry(deltaEtaStack[4][0]->hist_trunk.at(2), "3 < p_{T}^{trk}< 4 GeV","f");
  
  ptLegend3->AddEntry(deltaEtaStack[4][0]->hist_trunk.at(1), "4 < p_{T}^{trk}< 8 GeV","f");
  ptLegend3->AddEntry(deltaEtaStack[4][0]->hist_trunk.at(0), "8 < p_{T}^{trk}< 12 GeV","f");
  
  //ptLegend4->AddEntry(deltaEtaStack[4][0]->hist_trunk.at(6), "12 < p_{T}^{trk}< 300 GeV","f");
  ptLegend4->AddEntry(sumUncertainty[4][0], "0.7 < p_{T}^{trk}< 12 GeV","lpfe");
  
  
  double cmsPosition[] = {0.12,0.1,0.12};
  double cmsPositionY = 0.97;
  double jetShapeTitlePosition[] = {0.46,0.44,0.46};
  double xjPosition[] = {0.76,0.77,0.76};
  double systemPosition = 0.26;
  double selectionPosition = 0.205;
  double cmsSize = 0.035;
  
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
  mainTitle->DrawLatexNDC(.13, cmsPositionY, "CMS");
  
  mainTitle->SetTextFont(42);
  mainTitle->SetTextSize(0.035);
  mainTitle->DrawLatexNDC(0.24, 0.97, deltaEtaTitle[iJetTrack/3]);
  
  //mainTitle->SetTextFont(42);
  //mainTitle->SetTextSize(0.035);
  //mainTitle->DrawLatexNDC(xjPosition[iJetTrack/3], 0.94, xjString[iAsymmetry]);
  
  mainTitle->SetTextSize(0.03);
  mainTitle->DrawLatexNDC(systemPosition, 0.94, "5.02 TeV   pp 320 pb^{-1}   PbPb 1.7 nb^{-1}");
  mainTitle->SetTextSize(0.022);
  mainTitle->DrawLatexNDC(selectionPosition, 0.915, "anti-k_{T} R = 0.4, |#eta_{jet}| < 1.6, p_{T,1} > 120 GeV, p_{T,2} > 50 GeV, #Delta#phi_{1,2} > #frac{5#pi}{6}");
  //  lb->drawText("(p_{T}> 120 GeV, |#eta_{jet}|<1.6)", 0.2, 0.25, 4);
  //
  
  /*box = new TBox();
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
   mainTitle->DrawLatex(0.982, 0.0285, "1");*/
  
  //bigCanvas->SaveAs("js_dr_normal_new.eps");
  bigCanvas->SaveAs(Form("figures/final%s_bigCanvas_closeToFinal.pdf", deltaEtaSaveName[iJetTrack/3]));
  //bigCanvas->SaveAs("deltaEta_normal_v3.png");
  //bigCanvas->SaveAs("js_dr_normal_v3.pdf");
  
}

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void finalBigAsymmetryDeltaEtaPlotter_v1(){
  TString dir = "../data";
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Open data files for pp and PbPb data
  TFile *ppFile = TFile::Open("data/ppData2017_highForest_pfJets_20EveMixed_xjBins_wtaAxis_allCorrections_processed_2020-02-04.root");
  // data/ppData2017_highForest_pfJets_20EveMixed_xjBins_wtaAxis_allCorrections_processed_2020-01-31.root
  // data/ppData2017_highForest_pfJets_20EventsMixed_xjBins_finalTrackCorr_JECv4_wtaAxis_allCorrections_processed_2019-09-28.root
  // data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_processed_2019-08-13.root
  // data/dijet_pp_highForest_pfJets_noUncOrInc_allCorrections_wtaAxis_processed_2019-07-13.root
  TFile *pbpbFile = TFile::Open("data/dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_allCorrections_tuningForSeagull_wtaAxis_processed_2020-02-04.root");
  
  TFile *ppUncertaintyFile = TFile::Open("uncertainties/systematicUncertaintyForPp_20eveMix_xjBins_fixJES_2020-02-03.root");
  TFile *pbpbUncertaintyFile = TFile::Open("uncertainties/systematicUncertaintyForPbPb_25eveMix_xjBins_includeTrackDeltaR_2020-01-27.root");
  // uncertainties/systematicUncertaintyForPbPb_15percentSpill5Jff_2019-10-14.root
  // uncertainties/systematicUncertaintyForPbPb_15percentSpill20Jff_2019-10-01.root
  
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
  int iJetTrack = DijetHistogramManager::kTrackLeadingJet;  // DijetHistogramManager::kTrackSubleadingJet
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
  ppHistograms->SetLoad2DHistograms(true);
  ppHistograms->SetAsymmetryBinRange(0,nAsymmetryBins);
  ppHistograms->LoadProcessedHistograms();
  
  pbpbHistograms->SetLoadTrackLeadingJetCorrelations(true);
  pbpbHistograms->SetLoad2DHistograms(true);
  pbpbHistograms->SetAsymmetryBinRange(0,nAsymmetryBins);
  pbpbHistograms->LoadProcessedHistograms();
  
  // Load the systematic uncertainties from the files
  ppUncertaintyProvider->ReadSystematicFile(ppUncertaintyFile);
  pbpbUncertaintyProvider->ReadSystematicFile(pbpbUncertaintyFile);
  
  // Plot the figures using Xiao's plotting macro
  plotDeltaEtaBigAsymmetry(ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kTrackLeadingJet);
  plotDeltaEtaBigAsymmetry(ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kTrackSubleadingJet);
  
}
