#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JffCorrector.h"
#include "stackHist.h"
#include "xCanvas.h"
#include "JDrawer.h"
#include "DijetMethods.h"

void plotDeltaEtaXiao(DijetHistogramManager *ppHistograms, DijetHistogramManager *pbpbHistograms, JffCorrector *ppUncertaintyProvider, JffCorrector *pbpbUncertaintyProvider, int iJetTrack){

  // Determine the number of bins from the histogram manager
  const int nTrackPtBins = pbpbHistograms->GetNTrackPtBins();
  const int nCentralityBins = pbpbHistograms->GetNCentralityBins();
  const int nAsymmetryBins = pbpbHistograms->GetNAsymmetryBins();
  
  // Define labels and plotting layout
  const char* deltaEtaTitle[] = {"Particle yield vs. leading #Delta#eta","Particle yield vs. subleading #Delta#eta","Particle yield vs. #Delta#eta"};
  const char* deltaEtaSaveName[] = {"trackLeadingJet","trackSubleadingJet","trackInclusiveJet"};
  const char* xjString[] = {"0.0 < x_{j} < 0.6","0.6 < x_{j} < 0.8","0.8 < x_{j} < 1.0","x_{j} inclusive"};
  const char* asymmetrySaveName[] = {"_A=0v0-0v6","_A=0v6-0v8","_A=0v8-1v0",""};
  double deltaEtaZoom[] = {30,40,30};
  double subtractionZoom[] = {15,23,15};
  
  const bool drawInclusive = false;
  
  // Open a file to include results from inclusive jet analysis HIN-16-020
  TFile *inclusiveResultFile = TFile::Open("data/publishedResults/officialHist_py_deta_16_020.root");
  
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
  
  // Histgrams for comparison with inclusive results
  TH1D *sumHistogramInclusive[nCentralityBins+1];
  
  // Read two dimensional deltaEta-deltaPhi histograms and project the deltaEta yields out from them
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      helperHistogram = pbpbHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, 0);
      
      addedHistogram = projector->ProjectRegionDeltaEta(helperHistogram, -projectionRegion, projectionRegion, Form("DeltaEtaStack%d%d%d%d",iJetTrack,iAsymmetry,iCentrality,0));
      
      // Since we want to plot yield, we do not want to normalize over the number of bins projected over
      // but simply look at the yield in certain region
      addedHistogram->Scale(projector->GetNBinsProjectedOver());
      
      // The different pT bins can have different size, so we need to normalize the yield with the pT bin width
      addedHistogram->Scale(1/(pbpbHistograms->GetTrackPtBinBorder(1) - pbpbHistograms->GetTrackPtBinBorder(0)));
      
      // Rebin the histogram to match the binning in the inclusive jet shape paper
      sumHistogramPbPb[iCentrality][iAsymmetry] = (TH1D*) projector->RebinAsymmetric(addedHistogram, nDeltaEtaBinsRebin, deltaEtaBinBordersRebin)->Clone(Form("sumPbPb%d%d",iCentrality,iAsymmetry));

      deltaEtaArray[iCentrality][0][iAsymmetry] = (TH1D*) sumHistogramPbPb[iCentrality][iAsymmetry]->Clone(Form("deltaEtaArray%d%d%d", iCentrality, iAsymmetry, 0));
    }
    
    helperHistogram = ppHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, 0, 0);
    
    addedHistogram = projector->ProjectRegionDeltaEta(helperHistogram, -projectionRegion, projectionRegion, Form("DeltaEtaStack%d%d%d%d",iJetTrack,iAsymmetry,nCentralityBins,0));
    
    // Since we want to plot yield, we do not want to normalize over the number of bins projected over
    // but simply look at the yield in certain region
    addedHistogram->Scale(projector->GetNBinsProjectedOver());
    
    // The different pT bins can have different size, so we need to normalize the yield with the pT bin width
    addedHistogram->Scale(1/(ppHistograms->GetTrackPtBinBorder(1) - ppHistograms->GetTrackPtBinBorder(0)));
    
    // Rebin the histogram to match the binning in the inclusive jet shape paper
    sumHistogramPp[iAsymmetry] = (TH1D*) projector->RebinAsymmetric(addedHistogram, nDeltaEtaBinsRebin, deltaEtaBinBordersRebin)->Clone(Form("sumPp%d",iAsymmetry));
    
    deltaEtaArray[nCentralityBins][0][iAsymmetry] = (TH1D*) sumHistogramPp[iAsymmetry]->Clone(Form("deltaEtaArray%d%d%d", nCentralityBins, iAsymmetry, 0));
    
  }
  
  // Sum the pT:s
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    for(int iTrackPt = 1; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        
        helperHistogram = pbpbHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, iCentrality, iTrackPt);
        
        addedHistogram = projector->ProjectRegionDeltaEta(helperHistogram, -projectionRegion, projectionRegion, Form("DeltaEtaStack%d%d%d%d",iJetTrack,iAsymmetry,iCentrality,iTrackPt));
        
        // Since we want to plot yield, we do not want to normalize over the number of bins projected over
        // but simply look at the yield in certain region
        addedHistogram->Scale(projector->GetNBinsProjectedOver());
        
        // The different pT bins can have different size, so we need to normalize the yield with the pT bin width
        addedHistogram->Scale(1/(pbpbHistograms->GetTrackPtBinBorder(iTrackPt+1) - pbpbHistograms->GetTrackPtBinBorder(iTrackPt)));
        
        // Rebin the histogram to match the binning in the inclusive jet shape paper
        rebinnedHistogram = projector->RebinAsymmetric(addedHistogram,nDeltaEtaBinsRebin,deltaEtaBinBordersRebin);
        
        sumHistogramPbPb[iCentrality][iAsymmetry]->Add(rebinnedHistogram);
        
        deltaEtaArray[iCentrality][iTrackPt][iAsymmetry] = (TH1D*) rebinnedHistogram->Clone(Form("deltaEtaArray%d%d%d", iCentrality, iAsymmetry, iTrackPt));
      }
      
      helperHistogram = ppHistograms->GetHistogramJetTrackDeltaEtaDeltaPhi(iJetTrack, DijetHistogramManager::kBackgroundSubtracted, iAsymmetry, 0, iTrackPt);
      
      addedHistogram = projector->ProjectRegionDeltaEta(helperHistogram, -projectionRegion, projectionRegion, Form("DeltaEtaStack%d%d%d%d",iJetTrack,iAsymmetry,nCentralityBins,iTrackPt));
      
      // Since we want to plot yield, we do not want to normalize over the number of bins projected over
      // but simply look at the yield in certain region
      addedHistogram->Scale(projector->GetNBinsProjectedOver());
      
      // The different pT bins can have different size, so we need to normalize the yield with the pT bin width
      addedHistogram->Scale(1/(ppHistograms->GetTrackPtBinBorder(iTrackPt+1) - ppHistograms->GetTrackPtBinBorder(iTrackPt)));
      
      // Rebin the histogram to match the binning in the inclusive jet shape paper
      rebinnedHistogram = projector->RebinAsymmetric(addedHistogram,nDeltaEtaBinsRebin,deltaEtaBinBordersRebin);
      
      sumHistogramPp[iAsymmetry]->Add(rebinnedHistogram);
      
      deltaEtaArray[nCentralityBins][iTrackPt][iAsymmetry] = (TH1D*) rebinnedHistogram->Clone(Form("deltaEtaArray%d%d%d", nCentralityBins, iAsymmetry, iTrackPt));
    } // track pT
  }
  
  // Read the inclusive histograms to make a comparison to the published results
  int inclusiveCentralityBins[] = {0,10,30,50,100};
  int inclusiveTrackPtBins[] = {0,1,2,3,4,8,12};
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    sumHistogramInclusive[iCentrality] = (TH1D*) inclusiveResultFile->Get(Form("Proj_dEta_PbPb_Cent%d_Cent%d_TrkPt07_TrkPt1_Rebin_rebin", inclusiveCentralityBins[iCentrality], inclusiveCentralityBins[iCentrality+1]))->Clone(Form("inclusiveResults%d", iCentrality));
    
    for(int iTrackPt = 1; iTrackPt < 6; iTrackPt++){
      
      addedHistogram = (TH1D*) inclusiveResultFile->Get(Form("Proj_dEta_PbPb_Cent%d_Cent%d_TrkPt%d_TrkPt%d_Rebin_rebin", inclusiveCentralityBins[iCentrality], inclusiveCentralityBins[iCentrality+1], inclusiveTrackPtBins[iTrackPt], inclusiveTrackPtBins[iTrackPt+1]));
      
      sumHistogramInclusive[iCentrality]->Add(addedHistogram);
      
    } // Track pT loop for published inclusive results
    
  } // Centrality loop for published inclusive results
  
  sumHistogramInclusive[nCentralityBins] = (TH1D*) inclusiveResultFile->Get("Proj_dEta_pp_TrkPt07_TrkPt1_Rebin_rebin")->Clone("inclusiveResultsPp");
  
  for(int iTrackPt = 1; iTrackPt < 6; iTrackPt++){
    
    addedHistogram = (TH1D*) inclusiveResultFile->Get(Form("Proj_dEta_pp_TrkPt%d_TrkPt%d_Rebin_rebin", inclusiveTrackPtBins[iTrackPt], inclusiveTrackPtBins[iTrackPt+1]));
    
    sumHistogramInclusive[nCentralityBins]->Add(addedHistogram);
    
  } // Track pT loop for published inclusive results
  cout << "Giong to read uncertainties" << endl;
  // Find the uncertainties for the deltaEta histograms
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){

    // First, read the uncertainties in the lowest pT bin to the uncertainty histograms
    uncertaintyHistogramPp[iAsymmetry] = ppUncertaintyProvider->GetDeltaEtaSystematicUncertainty(iJetTrack, 0, 0, iAsymmetry);
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      uncertaintyHistogramPbPb[iCentrality][iAsymmetry] = pbpbUncertaintyProvider->GetDeltaEtaSystematicUncertainty(iJetTrack, iCentrality, 0, iAsymmetry);
    }

    // For the pT summed uncertainties, add the uncertainties in each pT bin
    for(int iTrackPt = 1; iTrackPt < 6; iTrackPt++){
      addedHistogram = ppUncertaintyProvider->GetDeltaEtaSystematicUncertainty(iJetTrack, 0, iTrackPt, iAsymmetry);
      uncertaintyHistogramPp[iAsymmetry]->Add(addedHistogram);
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        addedHistogram = pbpbUncertaintyProvider->GetDeltaEtaSystematicUncertainty(iJetTrack, iCentrality, iTrackPt, iAsymmetry);
        uncertaintyHistogramPbPb[iCentrality][iAsymmetry]->Add(addedHistogram);
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
      sumUncertainty[iCentrality][iAsymmetry]->SetMarkerSize(1.6);
      sumUncertainty[iCentrality][iAsymmetry]->SetMarkerColor(0);
      sumUncertainty[iCentrality][iAsymmetry]->SetLineColor(kBlack);
      
      deltaEtaStack[iCentrality][iAsymmetry] = new stackHist(Form("st_%d",iCentrality));
      deltaEtaStack[iCentrality][iAsymmetry]->setRange(-1.5, 1.5, "x");
      deltaEtaStack[iCentrality][iAsymmetry]->setRange(-1, deltaEtaZoom[iJetTrack/3], "y");
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
  
  // Draw all the distributions to big canvases
  auxi_canvas *bigCanvas[nAsymmetryBins+1];
  TBox *box[nAsymmetryBins+1];
  TLatex *mainTitle[nAsymmetryBins+1];
  TLine *line[nAsymmetryBins+1];
  
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
    
    bigCanvas[iAsymmetry] = new auxi_canvas(Form("bigCanvas%d%d",iAsymmetry,iJetTrack), "", 2500, 2000);
    bigCanvas[iAsymmetry]->SetMargin(0.075, 0.01, 0.08, 0.02);
    bigCanvas[iAsymmetry]->divide(3,4);
    
    mainTitle[iAsymmetry] = new TLatex();
    
    for(int iCentrality=0; iCentrality < nCentralityBins; iCentrality++){
      bigCanvas[iAsymmetry]->CD(8-iCentrality);
      deltaEtaStack[iCentrality][iAsymmetry]->drawStack("","hist",true);
      deltaEtaStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetNdivisions(505);
      deltaEtaStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetLabelSize(0.08);
      deltaEtaStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitleOffset(1.0);
      deltaEtaStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitleSize(0.1);
      deltaEtaStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitle("Y = #frac{1}{N_{dijet}} #frac{dN}{d#Delta#eta}");
      deltaEtaStack[iCentrality][iAsymmetry]->hst->Draw();
      if(iAsymmetry == nAsymmetryBins && iJetTrack == 0 && drawInclusive){
        sumHistogramInclusive[iCentrality]->Draw("same");
      }
      sumUncertainty[iCentrality][iAsymmetry]->Draw("same e2");
      if(iCentrality==3 ){
        mainTitle[iAsymmetry]->SetTextFont(22);
        mainTitle[iAsymmetry]->SetTextSize(.085);
        mainTitle[iAsymmetry]->DrawLatexNDC(.35, 0.9, "PbPb");
        mainTitle[iAsymmetry]->DrawLatexNDC(.35, .82, cent_lab[iCentrality]);
      }
      else {
        mainTitle[iAsymmetry]->SetTextFont(22);
        mainTitle[iAsymmetry]->SetTextSize(.09);
        mainTitle[iAsymmetry]->DrawLatexNDC(.19, 0.9, "PbPb");
        mainTitle[iAsymmetry]->DrawLatexNDC(.19, .82, cent_lab[iCentrality]);
      }
      bigCanvas[iAsymmetry]->CD(12-iCentrality);
      //ratio[iCentrality]->GetYaxis()->SetNdivisions(505);
      //subtractionHistogram[iCentrality][iAsymmetry]->SetTitle("");
      subtractedStack[iCentrality][iAsymmetry]->drawStack("","hist",true);
      subtractedStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitleOffset(1.1);
      subtractedStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitle("#Delta#eta");
      //subtractionHistogram[iCentrality][iAsymmetry]->SetAxisRange(0., 10, "Y");
      //subtractionHistogram[iCentrality][iAsymmetry]->SetAxisRange(0, .99, "X");
      if( iCentrality<3 )  {
        subtractedStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitleOffset(0.7);
        subtractedStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitleSize(0.11);
        subtractedStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetNdivisions(505);
        subtractedStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetLabelSize(0.08);
        subtractedStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetNdivisions(505);
      }
      if(iCentrality==3 ){
        subtractedStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitleOffset(0.92);
        subtractedStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetTitleSize(0.086);
        subtractedStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetNdivisions(505);
        subtractedStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetLabelOffset(0.016);
        subtractedStack[iCentrality][iAsymmetry]->hst->GetXaxis()->SetLabelSize(0.064);
        subtractedStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetNdivisions(505);
        subtractedStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetLabelOffset(0.02);
        subtractedStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetLabelSize(0.07);
        subtractedStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitleOffset(1.0);
        subtractedStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitleSize(0.08);
        subtractedStack[iCentrality][iAsymmetry]->hst->GetYaxis()->SetTitle("Y_{PbPb} - Y_{pp}");
        mainTitle[iAsymmetry]->SetTextSize(0.073);
        mainTitle[iAsymmetry]->DrawLatexNDC(.25, .92, "PbPb - pp");
        mainTitle[iAsymmetry]->DrawLatexNDC(.25, .84, cent_lab[iCentrality]);
      }
      else{
        mainTitle[iAsymmetry]->SetTextSize(0.09);
        mainTitle[iAsymmetry]->DrawLatexNDC(.05, .92, "PbPb - pp");
        mainTitle[iAsymmetry]->DrawLatexNDC(.05, .84, cent_lab[iCentrality]);
      }
      
      subtractedStack[iCentrality][iAsymmetry]->hst->Draw();
      
      //subtractionHistogram[iCentrality][iAsymmetry]->GetYaxis()->CenterTitle();
      //subtractionHistogram[iCentrality][iAsymmetry]->GetXaxis()->CenterTitle();
      //subtractionHistogram[iCentrality][iAsymmetry]->SetStats(0);
      //subtractionHistogram[iCentrality][iAsymmetry]->Draw("same");
      
      //line[iAsymmetry] = new TLine();
      //line[iAsymmetry]->SetLineStyle(2);
      //line[iAsymmetry]->DrawLine(0, 1, 1, 1);
      
      subtractionUncertainty[iCentrality][iAsymmetry]->Draw("same e2");
    }
    bigCanvas[iAsymmetry]->CD(1);
    deltaEtaStack[4][iAsymmetry]->drawStack("","hist",true);
    deltaEtaStack[4][iAsymmetry]->hst->GetYaxis()->SetNdivisions(505);
    deltaEtaStack[4][iAsymmetry]->hst->GetYaxis()->SetLabelSize(0.08);
    deltaEtaStack[4][iAsymmetry]->hst->GetYaxis()->SetTitleOffset(1.0);
    deltaEtaStack[4][iAsymmetry]->hst->GetYaxis()->SetTitleSize(0.1);
    deltaEtaStack[4][iAsymmetry]->hst->GetYaxis()->SetTitle("Y = #frac{1}{N_{dijet}} #frac{dN}{d#Delta#eta}");
    mainTitle[iAsymmetry]->SetTextSize(0.085);
    mainTitle[iAsymmetry]->DrawLatexNDC(0.35, 0.88, "pp reference");
    if(iAsymmetry == nAsymmetryBins && iJetTrack == 0 && drawInclusive){
      sumHistogramInclusive[nCentralityBins]->Draw("same");
    }
    sumUncertainty[4][iAsymmetry]->Draw("same e2");
    
    TLegend* lt1 = new TLegend(0.01, 0.1, 1, 0.5);
    TLegend* lt2 = new TLegend(0.0, 0.1, 1, 0.5);
    TLegend* lt3 = new TLegend(0.0 ,0.35, 1, 0.5);
    TLegend* lt5 = new TLegend(0.01, 0.5, 1, 0.68);
    //TLegend* lt4 = new TLegend(0.25, 0.85, 0.85, 0.95);
    lt1->SetTextSize(0.07);
    lt1->SetLineColor(kWhite);
    lt1->SetFillColor(kWhite);
    lt2->SetTextSize(0.07);
    lt2->SetLineColor(kWhite);
    lt2->SetFillColor(kWhite);
    lt3->SetTextSize(0.07);
    lt3->SetLineColor(kWhite);
    lt3->SetFillColor(kWhite);
    //lt4->SetTextSize(0.06);
    //lt4->SetLineColor(kWhite);
    //lt4->SetFillColor(kWhite);
    lt5->SetTextSize(0.07);
    lt5->SetLineColor(kWhite);
    lt5->SetFillColor(kWhite);
    
    lt1->AddEntry(deltaEtaStack[4][iAsymmetry]->hist_trunk.at(5), "0.7 < p_{T}^{trk}< 1 GeV","f");
    lt1->AddEntry(deltaEtaStack[4][iAsymmetry]->hist_trunk.at(4), "1 < p_{T}^{trk}< 2 GeV","f");
    lt1->AddEntry(deltaEtaStack[4][iAsymmetry]->hist_trunk.at(3), "2 < p_{T}^{trk}< 3 GeV","f");
    
    lt2->AddEntry(deltaEtaStack[4][iAsymmetry]->hist_trunk.at(2), "3 < p_{T}^{trk}< 4 GeV","f");
    lt2->AddEntry(deltaEtaStack[4][iAsymmetry]->hist_trunk.at(1), "4 < p_{T}^{trk}< 8 GeV","f");
    lt2->AddEntry(deltaEtaStack[4][iAsymmetry]->hist_trunk.at(0), "8 < p_{T}^{trk}< 12 GeV","f");
    
    
    lt3->AddEntry(sumUncertainty[4][iAsymmetry], "0.7 < p_{T}^{trk}< 12 GeV","lpfe");
    //lt3->AddEntry(deltaEtaStack[4]->hist_trunk.at(7), "16 < p_{T}^{trk}< 20 GeV","f");
    //lt3->AddEntry(deltaEtaStack[4]->hist_trunk.at(8), "20 < p_{T}^{trk}< 300 GeV","f");
    
    //lt4->AddEntry(ratioUncertainty[0][iAsymmetry], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
    //bigCanvas[iAsymmetry]->CD(9);
    //lt4->Draw();
    
    //bigCanvas->CD(1);
    //TLegend* lt6 = new TLegend(0.3,0.7,0.98,0.85);
    //lt6->SetTextSize(0.07);
    //lt6->SetLineColor(kWhite);
    //lt6->SetFillColor(kWhite);
    //lt6->AddEntry(js_dr_err_all[0], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
    //lt6->Draw();
    
    
    
    bigCanvas[iAsymmetry]->CD(2);
    //  lt5->AddEntry(js_dr_err_all[0], "0.7 < p_{T}^{track}< 20 GeV","lpfe");
    //  lt5->Draw();
    //line[iAsymmetry]->SetLineStyle(1);
    //line[iAsymmetry]->DrawLineNDC(0, 0, 0, 1);
    lt1->Draw();
    bigCanvas[iAsymmetry]->CD(3);
    lt2->Draw();
    bigCanvas[iAsymmetry]->CD(4);
    lt3->Draw();
    
    double cmsPosition[] = {0.37,0.35,0.37};
    double deltaEtaTitlePosition[] = {0.46,0.44,0.46};
    double xjPosition[] = {0.81,0.82,0.81};
    double overallShift[] = {-0.02,-0.02,-0.02,0};
    
    bigCanvas[iAsymmetry]->cd(0);
    mainTitle[iAsymmetry]->SetTextFont(62);
    mainTitle[iAsymmetry]->SetTextSize(0.035);
    mainTitle[iAsymmetry]->DrawLatexNDC(cmsPosition[iJetTrack/3]+overallShift[iAsymmetry], 0.94, "CMS");
    
    mainTitle[iAsymmetry]->SetTextFont(42);
    mainTitle[iAsymmetry]->SetTextSize(0.035);
    mainTitle[iAsymmetry]->DrawLatexNDC(deltaEtaTitlePosition[iJetTrack/3]+overallShift[iAsymmetry], 0.94, deltaEtaTitle[iJetTrack/3]);
    
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
    box[iAsymmetry]->DrawBox(0.518,.047, 0.533, 0.072);
    box[iAsymmetry]->DrawBox(0.75,.047, 0.765, 0.072);
    
    //bigCanvas->SaveAs("js_dr_normal_new.eps");
    bigCanvas[iAsymmetry]->SaveAs(Form("figures/finalDeltaEta_%s%s.pdf",deltaEtaSaveName[iJetTrack/3],asymmetrySaveName[iAsymmetry]));
    //bigCanvas->SaveAs("js_dr_normal_v3.eps");
    //bigCanvas->SaveAs("js_dr_normal_v3.pdf");
    
  } // Asymmetry loop for drawing the big canvas
}

/*
 * Macro for configuring the DijetDrawer and defining which histograms are drawn
 */
void finalDeltaEtaPlotter(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Open data files for pp and PbPb data
  TFile *ppFile = TFile::Open("data/ppData2017_highForest_pfJets_20EventsMixed_xjBins_finalTrackCorr_JECv4_wtaAxis_allCorrections_processed_2019-09-28.root");
  // data/ppData2017_highForest_pfJets_20eventsMixed_xjBins_JECv2_averagePeakMixing_wtaAxis_allCorrections_processed_2019-08-13.root
  // data/dijet_pp_highForest_pfJets_noUncOrInc_allCorrections_wtaAxis_processed_2019-07-13.root
  TFile *pbpbFile = TFile::Open("data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_allCorrections_lowPtResidualTrack_processed_2019-10-01_fiveJobsMissing.root");
  // data/dijetPbPb2018_highForest_akPu4CaloJets_jet80trigger_5eveMix_xjBins_wtaAxis_allCorrections_JECv5b_processed_2019-09-09.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_allCorrections_modifiedSeagull_wtaAxis_JECv4_processed_2019-08-13_fiveJobsMissing.root
  // data/dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_wtaAxis_JECv4_modifiedSeagull_noErrorJff_averagePeakMixing_processed_2019-08-13_fiveJobsMissing.root
  // dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_xjBins_allCorrections_modifiedSeagull_wtaAxis_JECv4_processed_2019-08-13_fiveJobsMissing.root
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
  plotDeltaEtaXiao(ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kTrackLeadingJet);
  plotDeltaEtaXiao(ppHistograms, pbpbHistograms, ppUncertaintyProvider, pbpbUncertaintyProvider, DijetHistogramManager::kTrackSubleadingJet);
  
}
