#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "../src/DijetHistograms.h"
#include "DijetCard.h"
#include "JDrawer.h"
#include <tuple>

/*
 * Fit a Gauss function to a histogram and extract parameters from that
 *
 *  TH1* histogram = Histogram from which drawing range in searched
 *
 *  return: Gauss mean, Gauss sigma, Error for Gauss mean, Error for Gauss sigma
 */
std::tuple<double,double,double,double> fitGauss(TH1* histogram, TString title = "", TString bin = "", TString saveName = ""){
  histogram->Fit("gaus","","",0.5,1.5);
  TF1 * gaussFit = histogram->GetFunction("gaus");
  
  double gaussMean = 0;
  double gaussSigma = 0;
  double gaussMeanError = 0;
  double gaussSigmaError = 0;
  
  if(gaussFit){
    gaussMean = gaussFit->GetParameter(1);
    gaussSigma = gaussFit->GetParameter(2);
    gaussMeanError = gaussFit->GetParError(1);
    gaussSigmaError = gaussFit->GetParError(2);
  }
  
  // If title is given, print the fit
  if(!title.EqualTo("")){
    JDrawer *temporaryDrawer = new JDrawer();
    temporaryDrawer->DrawHistogram(histogram,"Reco p_{T} / Gen p_{T}","Counts", " ");
    TLegend *legend = new TLegend(0.62,0.7,0.85,0.85);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->SetHeader(title);
    legend->AddEntry(histogram,bin,"l");
    legend->Draw();
    
    if(!saveName.EqualTo("")){
      gPad->GetCanvas()->SaveAs(saveName);
    }
    
  }
  
  return std::make_tuple(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError);
}

/*
 * Draw a closure histogram with quark/gluon discrimination to the plot
 *
 *  TH1D *histogram[DijetHistograms::knClosureTypes+1] = Array of histograms for each closure particle type
 *  const char* xTitle = Title given to the x-axis
 *  const char* yTitle = Title given to the y-axis
 *  bool ppData = Are we using pp od PbPb data in closure plotting
 *  int iClosureType = Leading/subleading/inclusive closures
 *  int iCentrality = Index of the centrality bin
 *  int iAsymmetry = Index of the xj bin
 *  int legendNzoom = Define the y-axis zoom and the legend position in the plot
 *  const char* saveComment = Comment given to the save name file
 *  bool saveFigures = Choose whether to save the figures or not
 */
void drawClosureHistogram(TH1D *histogram[DijetHistograms::knClosureTypes+1], const char* xTitle, const char* yTitle, bool ppData, int iClosureType, int iCentrality, int iAsymmetry, int legendNzoom, const char* saveComment, bool saveFigures){
  
  // Create a new drawer and define bin borders and drawing style
  JDrawer *drawer = new JDrawer();
  drawer->SetCanvasSize(700,700);
  double centralityBinBorders[5] = {0,10,30,50,90};
  double xjBinBorders[4] = {0.0, 0.6, 0.8, 1.0};
  char centralityString[100];
  char centralitySaveName[100];
  char namer[200];
  int lineColors[2] = {kBlue,kRed};
  const char* particleNames[2] = {"Quark","Gluon"};
  const char* closureTypeName[3] = {"Leading","Subleading","Inclusive"};
  const char* asymmetrySaveName[4] = {"_A0","_A1","_A2",""};
  
  // Zooming and legend position options
  double yZoomLow = 0.9;
  double yZoomHigh = 1.1;
  double legendX1 = 0.39;
  double legendX2 = 0.81;
  double legendY1 = 0.63;
  double legendY2 = 0.9;
  
  if(legendNzoom == 1){
    yZoomLow = 0;
    yZoomHigh = 0.34; // 0.24
    //legendX1 = 0.19;
    //legendX2 = 0.61;
    //legendY1 = 0.18;
    //legendY2 = 0.45;
    legendX1 = 0.46;
    legendX2 = 0.88;
    legendY1 = 0.65;
    legendY2 = 0.92;
  }
  
  // Set a good style for the inclusive histogram
  histogram[DijetHistograms::knClosureParticleTypes]->SetLineColor(kBlack);
  histogram[DijetHistograms::knClosureParticleTypes]->GetYaxis()->SetRangeUser(yZoomLow,yZoomHigh);
  drawer->DrawHistogram(histogram[DijetHistograms::knClosureParticleTypes],xTitle,yTitle," ");
  
  // Create a legend to the plot
  TLegend *legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
  
  if(ppData){
    sprintf(centralityString,", pp");
    sprintf(centralitySaveName,"");
  } else {
    sprintf(centralityString,", Cent:%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    sprintf(centralitySaveName,"_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
  }
  
  sprintf(namer,"%s jet%s",closureTypeName[iClosureType],centralityString);
  legend->SetHeader(namer);
  if(iAsymmetry < 3) legend->AddEntry((TObject*)0, Form("%.1f < x_{j} < %.1f", xjBinBorders[iAsymmetry], xjBinBorders[iAsymmetry+1]), "");
  legend->AddEntry(histogram[DijetHistograms::knClosureParticleTypes],"Inclusive","l");
  
  for(int iClosureParticle = 0; iClosureParticle < DijetHistograms::knClosureParticleTypes; iClosureParticle++){
    histogram[iClosureParticle]->SetLineColor(lineColors[iClosureParticle]);
    histogram[iClosureParticle]->Draw("same");
    legend->AddEntry(histogram[iClosureParticle],particleNames[iClosureParticle],"l");
  } // Closure particle loop (quark/gluon)
  
  legend->Draw();
  
  // Save the figures if selected to do so
  if(saveFigures){
    sprintf(namer,"figures/jet%s%sJet%s%s.pdf",saveComment,closureTypeName[iClosureType],centralitySaveName,asymmetrySaveName[iAsymmetry]);
    gPad->GetCanvas()->SaveAs(namer);
  }
  
}

/*
 * Macro for producing the JFF correction from PYTHIA simulation results
 */
void constructJetPtClosures(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString closureFileName = "data/PbPbMC2018_GenGen_akFlowJet_onlyJets_5pCentShift_jetClosuresWithXj_processed_2020-04-16.root";  // File from which the RecoGen histograms are read for the correction
  // data/PbPbMC2018_GenGen_akFlowJet_onlyJets_5pCentShift_smearCheck_jetPtClosure_processed_2020-05-14.root
  // data/PbPbMC2018_RecoReco_akFlowJet_onlyJets_5pCentShift_jetClosuresWithXj_processed_2020-04-16.root
  // data/PbPbMC2018_GenGen_akFlowJet_onlyJets_5pCentShift_jetClosuresWithXj_processed_2020-04-16.root
  // data/PbPbMC2018_GenGen_akFlowJet_onlyJets_5pCentShift_jet100trigger_jetClosuresWithXj_processed_2020-04-21.root
  // data/PbPbMC_GenGen_akPu4CaloJet_jetClosures_noMixing_matchedJets_eschemeAxis_oldJEC_processed_2019-09-20.root
  // data/PbPbMC_GenGen_akFlowPuCs4PfJets_jetsNtracks_JECv5b_jetClosures_processed_2019-09-03.root
  // data/PbPbMC_GenGen_skims_pfJets_noCorrelations_jetPtClosure_processed_2019-02-07.root
  // "data/PbPbMC_GenGen_pfJets_noCorrelations_jetPtClosure_processed_2019-06-25.root"
  // "data/PbPbMC_GenGen_pfCsJets_noCorrelations_jetPtClosure_eta1v3_2019-06-26_processed.root"
  // data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_onlyJets_jetPtClosure_JECv4_processed_2020-05-20.root
  // data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_onlyJets_20pSmear_jetPtClosure_JECv4_processed_2020-05-20.root
  
  bool drawLeadingClosure = true;       // Produce the closure plots for leading jet pT
  bool drawSubleadingClosure = true;    // Produce the closure plots for subleading jet pT
  bool drawInclusiveClosure = false;     // Produce the closure plots for inclusive jet pT
  
  bool drawPtClosure = true;
  bool drawEtaClosure = false;
  
  bool closureSelector[DijetHistograms::knClosureTypes] = {drawLeadingClosure,drawSubleadingClosure,drawInclusiveClosure};
  
  bool fitResolution = false;  // Fit the jet pT resolution histograms
  
  bool saveFigures = true ;  // Save the figures to file
  
  // ==================================================================
  // =================== Configuration ready ==========================
  // ==================================================================
  
  // Open the input files
  TFile *closureFile = TFile::Open(closureFileName);
  
  // Check if we are using pp or PbPb data
  DijetCard *card = new DijetCard(closureFile);
  TString collisionSystem = card->GetDataType();
  bool ppData = (collisionSystem.Contains("pp") || collisionSystem.Contains("localTest"));
  
  // Create histogram managers to provide the histograms for the correction
  DijetHistogramManager *closureHistograms = new DijetHistogramManager(closureFile);
  closureHistograms->SetLoadJetPtClosureHistograms(true);
  if(ppData) closureHistograms->SetCentralityBinRange(0,0);  // Disable centrality binning for pp data
  closureHistograms->SetAsymmetryBinRange(0,3);
  closureHistograms->LoadProcessedHistograms();
  
  // Find the correct number of centrality and track pT bins
  const int nCentralityBins = ppData ? 1 : closureHistograms->GetNCentralityBins();
  const int nTrackPtBins = closureHistograms->GetNTrackPtBins();
  const int nAsymmetryBins = closureHistograms->GetNAsymmetryBins();
  double centralityBinBorders[5] = {0,10,30,50,90};  // Bin borders for centrality
  
  // Choose which xj bins to draw
  int firstDrawnAsymmetryBin = 0;
  int lastDrawAsymmetryBin = 3;
  
  // Initialize reco/gen ratio and closure histograms
  TH1D *hRecoGenRatio[DijetHistograms::knClosureTypes][DijetHistogramManager::knGenJetPtBins][nCentralityBins][nAsymmetryBins+1][DijetHistograms::knClosureParticleTypes+1];
  TH1D *hRecoGenRatioEta[DijetHistograms::knClosureTypes][DijetHistogramManager::knJetEtaBins][nCentralityBins][nAsymmetryBins+1][DijetHistograms::knClosureParticleTypes+1];
  TH1D *hJetPtClosure[DijetHistograms::knClosureTypes][nCentralityBins][nAsymmetryBins+1][DijetHistograms::knClosureParticleTypes+1];
  TH1D *hJetPtClosureSigma[DijetHistograms::knClosureTypes][nCentralityBins][nAsymmetryBins+1][DijetHistograms::knClosureParticleTypes+1];
  TH1D *hJetPtClosureEta[DijetHistograms::knClosureTypes][nCentralityBins][nAsymmetryBins+1][DijetHistograms::knClosureParticleTypes+1];
  TH1D *hJetPtClosureSigmaEta[DijetHistograms::knClosureTypes][nCentralityBins][nAsymmetryBins+1][DijetHistograms::knClosureParticleTypes+1];
  char histogramNamer[100];
  
  for(int iClosureType = 0; iClosureType < DijetHistograms::knClosureTypes; iClosureType++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawAsymmetryBin; iAsymmetry++){
        for(int iClosureParticle = 0; iClosureParticle < DijetHistograms::knClosureParticleTypes+1; iClosureParticle++){
          for(int iGenJetPt = 0; iGenJetPt < DijetHistogramManager::knGenJetPtBins; iGenJetPt++){
            hRecoGenRatio[iClosureType][iGenJetPt][iCentrality][iAsymmetry][iClosureParticle] = NULL;
          } // Generator level jet pT loop
          for(int iJetEta = 0; iJetEta < DijetHistogramManager::knJetEtaBins; iJetEta++){
            hRecoGenRatioEta[iClosureType][iJetEta][iCentrality][iAsymmetry][iClosureParticle] = NULL;
          }
          sprintf(histogramNamer,"jetPtClosure_Corr%d_Cent%d_xj%d_Part%d", iClosureType, iCentrality, iAsymmetry, iClosureParticle);
          hJetPtClosure[iClosureType][iCentrality][iAsymmetry][iClosureParticle] = new TH1D(histogramNamer,histogramNamer,45,50,500);
          sprintf(histogramNamer,"jetPtClosureSigma_Corr%d_Cent%d_xj%d_Part%d",iClosureType, iCentrality, iAsymmetry, iClosureParticle);
          hJetPtClosureSigma[iClosureType][iCentrality][iAsymmetry][iClosureParticle] = new TH1D(histogramNamer,histogramNamer,45,50,500);
          sprintf(histogramNamer,"jetPtClosureEta_Corr%d_Cent%d_xj%d_Part%d",iClosureType, iCentrality, iAsymmetry, iClosureParticle);
          hJetPtClosureEta[iClosureType][iCentrality][iAsymmetry][iClosureParticle] = new TH1D(histogramNamer,histogramNamer,50,-2.5,2.5);
          sprintf(histogramNamer,"jetPtClosureSigmaEta_Corr%d_Cent%d_xj%d_Part%d",iClosureType, iCentrality, iAsymmetry, iClosureParticle);
          hJetPtClosureSigmaEta[iClosureType][iCentrality][iAsymmetry][iClosureParticle] = new TH1D(histogramNamer,histogramNamer,50,-2.5,2.5);
        } // Closure particle loop (quark/gluon/no selection)
      } // Asymmetry loop
    } // Centrality loop
  } // Closure type loop (leading/subleading/inclusive)
  
  JDrawer *drawer = new JDrawer();
  TF1* gaussFit;
  double gaussMean = 0;
  double gaussSigma = 0;
  double gaussMeanError = 0;
  double gaussSigmaError = 0;
  char namer[100];
  int minGenPt = 0;  // Skip bins below 120 GeV for leading jet
  
  // Read the reco/gen histograms from the file and fit them to construct the closure plots
  for(int iClosureType = 0; iClosureType < DijetHistograms::knClosureTypes; iClosureType++){
    if(!closureSelector[iClosureType]) continue;
    minGenPt = (iClosureType == 0) ? 7 : 0; // Skip the first jet pT bins for leading jet closure
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawAsymmetryBin; iAsymmetry++){
        for(int iClosureParticle = 0; iClosureParticle < DijetHistograms::knClosureParticleTypes+1; iClosureParticle++){
          if(drawPtClosure){
            for(int iGenJetPt = minGenPt; iGenJetPt < DijetHistogramManager::knGenJetPtBins; iGenJetPt++){
              // cout << "Fitting for pT: " << iGenJetPt << " closureParticle: " << iClosureParticle << " centrality: " << iCentrality << " closureType:" << iClosureType << endl;
              // Read the reco/gen histogram from the file
              
              hRecoGenRatio[iClosureType][iGenJetPt][iCentrality][iAsymmetry][iClosureParticle] = closureHistograms->GetHistogramJetPtClosure(iClosureType, iGenJetPt, DijetHistogramManager::knJetEtaBins, iCentrality, iAsymmetry, iClosureParticle);
              
              // Fit a gauss to the histogram
              //std::tie(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError) = fitGauss(hRecoGenRatio[iClosureType][iGenJetPt][iCentrality][iAsymmetry][iClosureParticle],"Leading jet 0-10 %",Form("%.0f < Gen pT < %.0f", hJetPtClosure[iClosureType][iCentrality][iAsymmetry][iClosureParticle]->GetBinLowEdge(iGenJetPt+1), hJetPtClosure[iClosureType][iCentrality][iAsymmetry][iClosureParticle]->GetBinLowEdge(iGenJetPt+2)), Form("figures/jetLeadingJetPtFit_GenPt%d.pdf",iGenJetPt));
              std::tie(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError) = fitGauss(hRecoGenRatio[iClosureType][iGenJetPt][iCentrality][iAsymmetry][iClosureParticle]);
              
              // Fill the histogram with the fit parameters
              hJetPtClosure[iClosureType][iCentrality][iAsymmetry][iClosureParticle]->SetBinContent(iGenJetPt+1,gaussMean);
              hJetPtClosure[iClosureType][iCentrality][iAsymmetry][iClosureParticle]->SetBinError(iGenJetPt+1,gaussMeanError);
              hJetPtClosureSigma[iClosureType][iCentrality][iAsymmetry][iClosureParticle]->SetBinContent(iGenJetPt+1,gaussSigma);
              hJetPtClosureSigma[iClosureType][iCentrality][iAsymmetry][iClosureParticle]->SetBinError(iGenJetPt+1,gaussSigmaError);
              
            } // Generator level jet pT loop
          } // pT closure if
          
          if(drawEtaClosure){
            // For eta, bins from 9 to nBins-9 cover the ragion -1.6 < eta < 1.6
            for(int iJetEta = 9; iJetEta < DijetHistogramManager::knJetEtaBins-9; iJetEta++){
              
              // Read the reco/gen histogram from the file
              hRecoGenRatioEta[iClosureType][iJetEta][iCentrality][iAsymmetry][iClosureParticle] = closureHistograms->GetHistogramJetPtClosure(iClosureType, DijetHistogramManager::knGenJetPtBins, iJetEta, iCentrality, iAsymmetry, iClosureParticle);
              
              // Fit a gauss to the histogram
              std::tie(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError) = fitGauss(hRecoGenRatioEta[iClosureType][iJetEta][iCentrality][iAsymmetry][iClosureParticle]);
              
              // Fill the histogram with the fit parameters
              hJetPtClosureEta[iClosureType][iCentrality][iAsymmetry][iClosureParticle]->SetBinContent(iJetEta+1,gaussMean);
              hJetPtClosureEta[iClosureType][iCentrality][iAsymmetry][iClosureParticle]->SetBinError(iJetEta+1,gaussMeanError);
              hJetPtClosureSigmaEta[iClosureType][iCentrality][iAsymmetry][iClosureParticle]->SetBinContent(iJetEta+1,gaussSigma);
              hJetPtClosureSigmaEta[iClosureType][iCentrality][iAsymmetry][iClosureParticle]->SetBinError(iJetEta+1,gaussSigmaError);
              
            } // Jet eta loop
          } // eta closure if
          
        } // Closure particle loop (quark/gluon/no selection)
      } // Asymmetry loop
    } // Centrality loop
  } // Closure type loop (leading/subleading/inclusive)
  
  double minFitPt = 120;
  double maxFitPt = 500;
  
  // Fit the resolution plots with a polynomial function
  if(fitResolution){
    for(int iClosureType = 0; iClosureType < DijetHistograms::knClosureTypes; iClosureType++){
      if(!closureSelector[iClosureType]) continue;
      
      minFitPt = 50;
      if(iClosureType == DijetHistograms::kLeadingClosure) minFitPt = 120;
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawAsymmetryBin; iAsymmetry++){
          
          cout << "Fitting iClosureType: " << iClosureType << " iCentrality: " << iCentrality << endl;
          hJetPtClosureSigma[iClosureType][iCentrality][iAsymmetry][DijetHistograms::knClosureParticleTypes]->Fit("pol4","","",minFitPt,maxFitPt);
          
        } // Asymmetry loop
      } // Centrality loop
    } // Closure type loop (leading/subleading/inclusive)
  } // Fitting the resolution
  
  // Draw the closure plots
  int lineColors[2] = {kBlue,kRed};
  const char* particleNames[2] = {"Quark","Gluon"};
  const char* closureTypeName[3] = {"Leading","Subleading","Inclusive"};
  double adders[3] = {0,0.05,0};
  TLegend *legend;
  char centralityString[100];
  char centralitySaveName[100];
  
  // Create the output file
  TFile *outputFile = new TFile("closureLul.root","UPDATE");
  
  for(int iClosureType = 0; iClosureType < DijetHistograms::knClosureTypes; iClosureType++){
    if(!closureSelector[iClosureType]) continue;
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iAsymmetry = firstDrawnAsymmetryBin; iAsymmetry <= lastDrawAsymmetryBin; iAsymmetry++){
        if(drawPtClosure){
          drawClosureHistogram(hJetPtClosure[iClosureType][iCentrality][iAsymmetry], "Gen p_{T} (GeV)", "#mu(reco p_{T} / gen p_{T})", ppData, iClosureType, iCentrality, iAsymmetry, 0, "PtClosureGenDijet", saveFigures);
          
          drawClosureHistogram(hJetPtClosureSigma[iClosureType][iCentrality][iAsymmetry], "Gen p_{T} (GeV)", "#sigma(reco p_{T} / gen p_{T})", ppData, iClosureType, iCentrality, iAsymmetry, 1, "PtResolution", saveFigures);
          
          hJetPtClosureSigma[iClosureType][iCentrality][iAsymmetry][2]->Write();
        }
        
        if(drawEtaClosure){
          drawClosureHistogram(hJetPtClosureEta[iClosureType][iCentrality][iAsymmetry], "#eta", "#mu(reco p_{T} / gen p_{T})", ppData, iClosureType, iCentrality, iAsymmetry, 0, "EtaClosure", saveFigures);
          
          drawClosureHistogram(hJetPtClosureSigmaEta[iClosureType][iCentrality][iAsymmetry], "Gen p_{T} (GeV)", "#sigma(reco p_{T} / gen p_{T})", ppData, iClosureType, iCentrality, iAsymmetry, 1, "EtaResolution", saveFigures);
          
        }
        
      } // Asymmetry loop
    } // Centrality loop
  } // Closure type loop (leading/subleading/inclusive)
  
  outputFile->Close();
  
}
