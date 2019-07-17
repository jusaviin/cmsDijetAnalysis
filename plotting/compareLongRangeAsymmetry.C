#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "DijetMethods.h"
#include "JDrawer.h"

/*
 * Print a slide with fit chi2/NDF information to the console
 */
void printChi2Slide(double (*chi2PerNdf)[5][7], int nAsymmetryBins){
  
  // Create a histogram manager to facilitate binning info exchange
  const int nTrackPtBins = 6;
  const int nCentralityBins = 4;
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12};  // Bin borders for track pT
  double asymmetryBinBorders[] = {0,0.6,0.8,1};     // Bin borders for xj
  
  char namer[100];
  
  for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
  
    if(iAsymmetry == nAsymmetryBins){
      sprintf(namer,"\\frametitle{$\\chi^{2}$ / NDF for long range}");
    } else {
      sprintf(namer,"\\frametitle{$\\chi^{2}$ / NDF for long range, $%.1f < x_{j} < %.1f$}", asymmetryBinBorders[iAsymmetry], asymmetryBinBorders[iAsymmetry+1]);
    }
    
    cout << endl;
    cout << "\\begin{frame}" << endl;
    cout << namer << endl;
    cout << "\\begin{center}" << endl;
    cout << "  \\begin{tabular}{cccccc}" << endl;
    cout << "    \\toprule" << endl;
    cout << "    $p_{\\mathrm{T}} (GeV)$ & C: 0-10 & C: 10-30 & C: 30-50 & C: 50-100 & pp \\\\" << endl;
    cout << "    \\midrule" << endl;
    
    // Set the correct precision for printing floating point numbers
    cout << fixed << setprecision(3);
    
    // Print one line for each track pT bin
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      sprintf(namer,"    %.1f-%.1f",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
      cout << namer;
      for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
        cout << " & $" << chi2PerNdf[iAsymmetry][iCentrality][iTrackPt] << "$";
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

/*
 * Macro for producing the JFF correction from PYTHIA simulation results
 */
void compareLongRangeAsymmetry(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString pbpbFileName = "data/dijetPbPb_pfCsJets_xjBins_wtaAxis_noUncOrInc_improvisedMixing_allCorrections_adjustedBackground_processed_2019-07-05.root";
  TString ppFileName = "data/dijet_pp_highForest_pfJets_noUncOrInc_allCorrections_wtaAxis_processed_2019-07-13.root";
  
  const bool printChi2 = false;
  const bool saveFigures = false;
  TString saveComment = "";
  
  const int fourierV = 0;  // Select which vn component to draw. 0 = All, 1...4 = v1...v4
  TString vString = "";
  TString vHeader = "v_{1} - v_{4}";
  if(fourierV > 0) {
    vString = Form("_v%d",fourierV);
    vHeader = Form("v_{%d}",fourierV);
  }
  
  // Open the input files
  TFile *pbpbDataFile = TFile::Open(pbpbFileName);
  TFile *ppDataFile = TFile::Open(ppFileName);
  
  // Create readers for the histograms and define binning
  DijetHistogramManager *pbpbReader = new DijetHistogramManager(pbpbDataFile);
  DijetHistogramManager *ppReader = new DijetHistogramManager(ppDataFile);
  
  const int nCentralityBins = pbpbReader->GetNCentralityBins();
  const int nTrackPtBins = pbpbReader->GetNTrackPtBins();
  const int nAsymmetryBins = pbpbReader->GetNAsymmetryBins();
  double centralityBinBorders[] = {0,10,30,50,100};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  
  // Define a fitter to test fitting the background distribution with different amount of vn:s
  DijetMethods *refitter = new DijetMethods();
  const int nRefit = 4; // Number of vn:s included in the refit
  bool refitBackground = true; // Refit the background distribution
  
  if(refitBackground && fourierV == 0) vHeader = Form("v_{1} - v_{%d}",nRefit);
  
  // Define arrays for the histograms and supporting numbers
  TH1D *background[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  TF1 *backgroundFit[nAsymmetryBins+1][nCentralityBins+1][nTrackPtBins];
  
  // To pass the array as an argument, the numbers will need to be given by hand.
  if(nAsymmetryBins != 3 || nCentralityBins != 4 || nTrackPtBins != 7){
    cout << "Wrong binning for chi2/ndf histogram!" << endl;
    cout << "nAsymmetry = " << nAsymmetryBins << " expected = 3" << endl;
    cout << "nCentralityBins = " << nCentralityBins << " expected = 4" << endl;
    cout << "nTrackPtBins = " << nTrackPtBins << " expected = 7" << endl;
    cout << "Fix the binning and the code will run again!" << endl;
    return;
  }
  double chi2PerNdf[4][5][7];
  
  // Read the histograms from the input files
  pbpbReader->SetLoadTrackLeadingJetCorrelations(true);
  pbpbReader->SetAsymmetryBinRange(0,nAsymmetryBins);
  pbpbReader->LoadProcessedHistograms();
  
  ppReader->SetLoadTrackLeadingJetCorrelations(true);
  ppReader->SetAsymmetryBinRange(0,nAsymmetryBins);
  ppReader->LoadProcessedHistograms();
  
  double chi2;  // Chi2 for fit quality estimation
  int ndf;      // Number of degrees of freedom for fit quality estimation
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iAsymmetry = 0; iAsymmetry <= nAsymmetryBins; iAsymmetry++){
        if(iCentrality == nCentralityBins){
          background[iAsymmetry][nCentralityBins][iTrackPt] = ppReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, iAsymmetry, 0, iTrackPt, DijetHistogramManager::kWholeEta);
        } else {
          background[iAsymmetry][iCentrality][iTrackPt] = pbpbReader->GetHistogramJetTrackDeltaPhi(DijetHistogramManager::kTrackLeadingJet, DijetHistogramManager::kBackground, iAsymmetry, iCentrality, iTrackPt, DijetHistogramManager::kWholeEta);
        }
        backgroundFit[iAsymmetry][iCentrality][iTrackPt] = background[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        
        // Option to refit the background distribution
        if(refitBackground){
          background[iAsymmetry][iCentrality][iTrackPt]->RecursiveRemove(backgroundFit[iAsymmetry][iCentrality][iTrackPt]);
          background[iAsymmetry][iCentrality][iTrackPt]->Rebin(4);
          background[iAsymmetry][iCentrality][iTrackPt]->Scale(1/4.0);
          refitter->FourierFit(background[iAsymmetry][iCentrality][iTrackPt], nRefit);
          backgroundFit[iAsymmetry][iCentrality][iTrackPt] = background[iAsymmetry][iCentrality][iTrackPt]->GetFunction("fourier");
        }
        
        chi2 = backgroundFit[iAsymmetry][iCentrality][iTrackPt]->GetChisquare();
        ndf = backgroundFit[iAsymmetry][iCentrality][iTrackPt]->GetNDF();
        chi2PerNdf[iAsymmetry][iCentrality][iTrackPt] = chi2/ndf;
      } // asymmetry loop
    } // track pT loop
  } // centrality loop
  
  if(printChi2) printChi2Slide(chi2PerNdf,nAsymmetryBins);
  
  // Draw all asymmery bins deltaR curves to a single plot to see if the spillover is similar regardless of asymmetry bin
  JDrawer *drawer = new JDrawer();
  
  // Make subeNon0 to spillover comparison in all bins
  TLegend *legend;
  TLegend *vLegend;
  int colors[] = {kBlue,kRed,kGreen+4,kBlack};
  double vx1 = 0.2; double vx2 = 0.25;  // x position of vn information legend
  if(fourierV == 0){
    vx1 = 0.183; vx2 = 0.233;
  }
  
  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      
      // Setup a legend for the figure
      legend = new TLegend(0.35,0.75,0.75,0.99);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      if(iCentrality == nCentralityBins){
        legend->SetHeader(Form("%.1f < p_{T} < %.1f, pp",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]));
      } else {
        legend->SetHeader(Form("%.1f < p_{T} < %.1f, C: %.0f-%.0f %%",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1], centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]));
      }

      // Draw the first histogram to set scale
      drawer->DrawHistogram(background[0][iCentrality][iTrackPt],"#Delta#phi","#frac{dN}{d#Delta#phi}"," ");
      
      // Draw all the asymmetry bins to the same figure
      for(int iAsymmetry = 0; iAsymmetry < nAsymmetryBins; iAsymmetry++){
        background[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(colors[iAsymmetry]);
        background[iAsymmetry][iCentrality][iTrackPt]->Draw("same");
        backgroundFit[iAsymmetry][iCentrality][iTrackPt]->SetLineColor(colors[iAsymmetry]);
        
        // Select the wanted vn component
        if(fourierV > 0){
          for(int iFitComponent = 1; iFitComponent < backgroundFit[iAsymmetry][iCentrality][iTrackPt]->GetNpar(); iFitComponent++){
            if(iFitComponent != fourierV) backgroundFit[iAsymmetry][iCentrality][iTrackPt]->SetParameter(iFitComponent,0);
          }
        }
        
        // Add the histogram to legend
        legend->AddEntry(background[iAsymmetry][iCentrality][iTrackPt],Form("%.1f < x_{j} < %.1f",xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]),"l");
      } // asymmetry loop
      
      // Draw the legend
      legend->Draw();
      
      // Make a legend for the fourier component
      vLegend = new TLegend(vx1,0.83,vx2,0.86);
      vLegend->SetFillStyle(0);vLegend->SetBorderSize(0);vLegend->SetTextSize(0.06);vLegend->SetTextFont(62);
      vLegend->SetHeader(vHeader);
      vLegend->Draw();
      
      // Save the figures to file
      if(saveFigures){
        if(iCentrality == nCentralityBins){
          gPad->GetCanvas()->SaveAs(Form("figures/longRangeAsymmetryComparison%s%s_pp_pT=%.0f-%.0f.pdf", saveComment.Data(), vString.Data(), trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
        } else {
          gPad->GetCanvas()->SaveAs(Form("figures/longRangeAsymmetryComparison%s%s_C=%.0f-%.0f_pT=%.0f-%.0f.pdf", saveComment.Data(), vString.Data(), centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1], trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]));
        }
      }
    } // Track pt Loop
  } // Centrality loop
  
}

