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
std::tuple<double,double,double,double> fitGauss(TH1* histogram){
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
 *  int legendNzoom = Define the y-axis zoom and the legend position in the plot
 *  const char* saveComment = Comment given to the save name file
 *  bool saveFigures = Choose whether to save the figures or not
 */
void drawClosureHistogram(TH1D *histogram[DijetHistograms::knClosureTypes+1], const char* xTitle, const char* yTitle, bool ppData, int iClosureType, int iCentrality, int legendNzoom, const char* saveComment, bool saveFigures){
  
  // Create a new drawer and define bin borders and drawing style
  JDrawer *drawer = new JDrawer();
  drawer->SetCanvasSize(700,700);
  double centralityBinBorders[5] = {0,10,30,50,100};
  char centralityString[100];
  char centralitySaveName[100];
  char namer[200];
  int lineColors[2] = {kBlue,kRed};
  const char* particleNames[2] = {"Quark","Gluon"};
  const char* closureTypeName[3] = {"Leading","Subleading","Inclusive"};
  
  // Zooming and legend position options
  double yZoomLow = 0.9;
  double yZoomHigh = 1.1;
  double legendX1 = 0.39;
  double legendX2 = 0.81;
  double legendY1 = 0.63;
  double legendY2 = 0.9;
  
  if(legendNzoom == 1){
    yZoomLow = 0;
    yZoomHigh = 0.2;
    legendY1 = 0.18;
    legendY2 = 0.45;
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
  legend->AddEntry(histogram[DijetHistograms::knClosureParticleTypes],"Inclusive","l");
  
  for(int iClosureParticle = 0; iClosureParticle < DijetHistograms::knClosureParticleTypes; iClosureParticle++){
    histogram[iClosureParticle]->SetLineColor(lineColors[iClosureParticle]);
    histogram[iClosureParticle]->Draw("same");
    legend->AddEntry(histogram[iClosureParticle],particleNames[iClosureParticle],"l");
  } // Closure particle loop (quark/gluon)
  
  legend->Draw();
  
  // Save the figures if selected to do so
  if(saveFigures){
    sprintf(namer,"figures/jetPtClosure%s%sJet%s.pdf",saveComment,closureTypeName[iClosureType],centralitySaveName);
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
  
  TString closureFileName = "data/PbPbMC_GenGen_skims_pfJets_noCorrelations_jetPtClosure_processed_2019-02-07.root";  // File from which the RecoGen histograms are read for the correction
  // "data/PbPbMC_GenGen_pfJets_noCorrelations_jetPtClosure_processed_2019-06-25.root"
  // "data/PbPbMC_GenGen_pfCsJets_noCorrelations_jetPtClosure_eta1v3_2019-06-26_processed.root"
  
  bool drawLeadingClosure = true;       // Produce the closure plots for leading jet pT
  bool drawSubleadingClosure = false;    // Produce the closure plots for subleading jet pT
  bool drawInclusiveClosure = false;     // Produce the closure plots for inclusive jet pT
  
  bool drawPtClosure = true;
  bool drawEtaClosure = false;
  
  bool closureSelector[DijetHistograms::knClosureTypes] = {drawLeadingClosure,drawSubleadingClosure,drawInclusiveClosure};
  
  bool saveFigures = false;  // Save the figures to file
  
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
  closureHistograms->LoadProcessedHistograms();
  
  // Find the correct number of centrality and track pT bins
  const int nCentralityBins = ppData ? 1 : closureHistograms->GetNCentralityBins();
  const int nTrackPtBins = closureHistograms->GetNTrackPtBins();
  double centralityBinBorders[5] = {0,10,30,50,100};  // Bin borders for centrality
  
  // Initialize reco/gen ratio and closure histograms
  TH1D *hRecoGenRatio[DijetHistograms::knClosureTypes][DijetHistogramManager::knGenJetPtBins][nCentralityBins][DijetHistograms::knClosureParticleTypes+1];
  TH1D *hRecoGenRatioEta[DijetHistograms::knClosureTypes][DijetHistogramManager::knJetEtaBins][nCentralityBins][DijetHistograms::knClosureParticleTypes+1];
  TH1D *hJetPtClosure[DijetHistograms::knClosureTypes][nCentralityBins][DijetHistograms::knClosureParticleTypes+1];
  TH1D *hJetPtClosureSigma[DijetHistograms::knClosureTypes][nCentralityBins][DijetHistograms::knClosureParticleTypes+1];
  TH1D *hJetPtClosureEta[DijetHistograms::knClosureTypes][nCentralityBins][DijetHistograms::knClosureParticleTypes+1];
  TH1D *hJetPtClosureSigmaEta[DijetHistograms::knClosureTypes][nCentralityBins][DijetHistograms::knClosureParticleTypes+1];
  char histogramNamer[100];
  
  for(int iClosureType = 0; iClosureType < DijetHistograms::knClosureTypes; iClosureType++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iClosureParticle = 0; iClosureParticle < DijetHistograms::knClosureParticleTypes+1; iClosureParticle++){
        for(int iGenJetPt = 0; iGenJetPt < DijetHistogramManager::knGenJetPtBins; iGenJetPt++){
          hRecoGenRatio[iClosureType][iGenJetPt][iCentrality][iClosureParticle] = NULL;
        } // Generator level jet pT loop
        for(int iJetEta = 0; iJetEta < DijetHistogramManager::knJetEtaBins; iJetEta++){
          hRecoGenRatioEta[iClosureType][iJetEta][iCentrality][iClosureParticle] = NULL;
        }
        sprintf(histogramNamer,"jetPtClosure_Corr%d_Cent%d_Part%d",iClosureType,iCentrality,iClosureParticle);
        hJetPtClosure[iClosureType][iCentrality][iClosureParticle] = new TH1D(histogramNamer,histogramNamer,45,50,500);
        sprintf(histogramNamer,"jetPtClosureSigma_Corr%d_Cent%d_Part%d",iClosureType,iCentrality,iClosureParticle);
        hJetPtClosureSigma[iClosureType][iCentrality][iClosureParticle] = new TH1D(histogramNamer,histogramNamer,45,50,500);
        sprintf(histogramNamer,"jetPtClosureEta_Corr%d_Cent%d_Part%d",iClosureType,iCentrality,iClosureParticle);
        hJetPtClosureEta[iClosureType][iCentrality][iClosureParticle] = new TH1D(histogramNamer,histogramNamer,50,-2.5,2.5);
        sprintf(histogramNamer,"jetPtClosureSigmaEta_Corr%d_Cent%d_Part%d",iClosureType,iCentrality,iClosureParticle);
        hJetPtClosureSigmaEta[iClosureType][iCentrality][iClosureParticle] = new TH1D(histogramNamer,histogramNamer,50,-2.5,2.5);
      } // Closure particle loop (quark/gluon/no selection)
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
      for(int iClosureParticle = 0; iClosureParticle < DijetHistograms::knClosureParticleTypes+1; iClosureParticle++){
        if(drawPtClosure){
          for(int iGenJetPt = minGenPt; iGenJetPt < DijetHistogramManager::knGenJetPtBins; iGenJetPt++){
            cout << "Fitting for pT: " << iGenJetPt << endl;
            // Read the reco/gen histogram from the file
            hRecoGenRatio[iClosureType][iGenJetPt][iCentrality][iClosureParticle] = closureHistograms->GetHistogramJetPtClosure(iClosureType, iGenJetPt, DijetHistogramManager::knJetEtaBins, iCentrality, iClosureParticle);
            
            // Fit a gauss to the histogram
            std::tie(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError) = fitGauss(hRecoGenRatio[iClosureType][iGenJetPt][iCentrality][iClosureParticle]);
            
            // Fill the histogram with the fit parameters
            hJetPtClosure[iClosureType][iCentrality][iClosureParticle]->SetBinContent(iGenJetPt+1,gaussMean);
            hJetPtClosure[iClosureType][iCentrality][iClosureParticle]->SetBinError(iGenJetPt+1,gaussMeanError);
            hJetPtClosureSigma[iClosureType][iCentrality][iClosureParticle]->SetBinContent(iGenJetPt+1,gaussSigma);
            hJetPtClosureSigma[iClosureType][iCentrality][iClosureParticle]->SetBinError(iGenJetPt+1,gaussSigmaError);
            
          } // Generator level jet pT loop
        } // pT closure if
        
        if(drawEtaClosure){
          for(int iJetEta = 0; iJetEta < DijetHistogramManager::knJetEtaBins; iJetEta++){
            
            // Read the reco/gen histogram from the file
            hRecoGenRatioEta[iClosureType][iJetEta][iCentrality][iClosureParticle] = closureHistograms->GetHistogramJetPtClosure(iClosureType, DijetHistogramManager::knGenJetPtBins, iJetEta, iCentrality, iClosureParticle);
            
            // Fit a gauss to the histogram
            std::tie(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError) = fitGauss(hRecoGenRatioEta[iClosureType][iJetEta][iCentrality][iClosureParticle]);
            
            // Fill the histogram with the fit parameters
            hJetPtClosureEta[iClosureType][iCentrality][iClosureParticle]->SetBinContent(iJetEta+1,gaussMean);
            hJetPtClosureEta[iClosureType][iCentrality][iClosureParticle]->SetBinError(iJetEta+1,gaussMeanError);
            hJetPtClosureSigmaEta[iClosureType][iCentrality][iClosureParticle]->SetBinContent(iJetEta+1,gaussSigma);
            hJetPtClosureSigmaEta[iClosureType][iCentrality][iClosureParticle]->SetBinError(iJetEta+1,gaussSigmaError);
            
          } // Jet eta loop
        } // eta closure if
        
      } // Closure particle loop (quark/gluon/no selection)
    } // Centrality loop
  } // Closure type loop (leading/subleading/inclusive)
 
  // Draw the closure plots
  int lineColors[2] = {kBlue,kRed};
  const char* particleNames[2] = {"Quark","Gluon"};
  const char* closureTypeName[3] = {"Leading","Subleading","Inclusive"};
  double adders[3] = {0,0.05,0};
  TLegend *legend;
  char centralityString[100];
  char centralitySaveName[100];
  for(int iClosureType = 0; iClosureType < DijetHistograms::knClosureTypes; iClosureType++){
    if(!closureSelector[iClosureType]) continue;
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      if(drawPtClosure){
        drawClosureHistogram(hJetPtClosure[iClosureType][iCentrality], "Gen p_{T} (GeV)", "#mu(reco p_{T} / gen p_{T})", ppData, iClosureType, iCentrality, 0, "", saveFigures);
        
        drawClosureHistogram(hJetPtClosureSigma[iClosureType][iCentrality], "Gen p_{T} (GeV)", "#sigma(reco p_{T} / gen p_{T})", ppData, iClosureType, iCentrality, 1, "", saveFigures);
      }
      
      if(drawEtaClosure){
        drawClosureHistogram(hJetPtClosureEta[iClosureType][iCentrality], "#eta", "#mu(reco p_{T} / gen p_{T})", ppData, iClosureType, iCentrality, 0, "", saveFigures);
        
        drawClosureHistogram(hJetPtClosureSigmaEta[iClosureType][iCentrality], "Gen p_{T} (GeV)", "#sigma(reco p_{T} / gen p_{T})", ppData, iClosureType, iCentrality, 1, "", saveFigures);
        
      }

    } // Centrality loop
  } // Closure type loop (leading/subleading/inclusive)
  
}
