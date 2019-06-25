#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "../src/DijetHistograms.h"
#include "DijetCard.h"
#include "JDrawer.h"

/*
 * Macro for producing the JFF correction from PYTHIA simulation results
 */
void constructJetPtClosures(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString closureFileName = "data/PbPbMC_RecoReco_skims_pfJets_noUncorr_noCorrelations_jetPtClosure_noJetLimit_processed_2019-02-02.root";  // File from which the RecoGen histograms are read for the correction
  // "data/PbPbMC_RecoReco_skims_pfJets_noUncorr_noCorrelations_jetPtClosure_noJetLimit_processed_2019-02-02.root"
  // "data/dijet_ppMC_RecoReco_mergedSkims_Pythia6_pfJets_noCorrelations_jetPtClosure_processed_2019-02-06.root"
  // "data/PbPbMC_GenGen_skims_pfJets_noCorrelations_jetPtClosure_processed_2019-02-07.root"
  
  bool drawLeadingClosure = false;       // Produce the closure plots for leading jet pT
  bool drawSubleadingClosure = true;    // Produce the closure plots for subleading jet pT
  bool drawInclusiveClosure = true;     // Produce the closure plots for inclusive jet pT
  
  // TODO: If leading jet closure, min GenPtBin should be 7
  
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
  TH1D *hJetPtClosure[DijetHistograms::knClosureTypes][nCentralityBins][DijetHistograms::knClosureParticleTypes+1];
  TH1D *hJetPtClosureSigma[DijetHistograms::knClosureTypes][nCentralityBins][DijetHistograms::knClosureParticleTypes+1];
  char histogramNamer[100];
  
  for(int iClosureType = 0; iClosureType < DijetHistograms::knClosureTypes; iClosureType++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iClosureParticle = 0; iClosureParticle < DijetHistograms::knClosureParticleTypes+1; iClosureParticle++){
        for(int iGenJetPt = 0; iGenJetPt < DijetHistogramManager::knGenJetPtBins; iGenJetPt++){
          hRecoGenRatio[iClosureType][iGenJetPt][iCentrality][iClosureParticle] = NULL;
        } // Generator level jet pT loop
        sprintf(histogramNamer,"jetPtClosure_Corr%d_Cent%d_Part%d",iClosureType,iCentrality,iClosureParticle);
        hJetPtClosure[iClosureType][iCentrality][iClosureParticle] = new TH1D(histogramNamer,histogramNamer,45,50,500);
        sprintf(histogramNamer,"jetPtClosureSigma_Corr%d_Cent%d_Part%d",iClosureType,iCentrality,iClosureParticle);
        hJetPtClosureSigma[iClosureType][iCentrality][iClosureParticle] = new TH1D(histogramNamer,histogramNamer,45,50,500);
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
  
  // Read the reco/gen histograms from the file and fit them to construct the closure plots
  for(int iClosureType = 0; iClosureType < DijetHistograms::knClosureTypes; iClosureType++){
    if(!closureSelector[iClosureType]) continue;
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iClosureParticle = 0; iClosureParticle < DijetHistograms::knClosureParticleTypes+1; iClosureParticle++){
        for(int iGenJetPt = 0; iGenJetPt < DijetHistogramManager::knGenJetPtBins; iGenJetPt++){
          
          // Read the reco/gen histogram from the file
          hRecoGenRatio[iClosureType][iGenJetPt][iCentrality][iClosureParticle] = closureHistograms->GetHistogramJetPtClosure(iClosureType, iGenJetPt, DijetHistogramManager::knJetEtaBins, iCentrality, iClosureParticle);
          
          // Fit the histogram with a gauss function
          hRecoGenRatio[iClosureType][iGenJetPt][iCentrality][iClosureParticle]->Fit("gaus","","",0.5,1.5);
          
          // Draw the function with the fit for debunking purposes
          if(iCentrality == 0 && iClosureParticle == DijetHistograms::knClosureParticleTypes){
            //sprintf(namer,"Reco/Gen GenPt%d",iGenJetPt);
            //drawer->DrawHistogram(hRecoGenRatio[iClosureType][iGenJetPt][iCentrality][iClosureParticle],namer,"N");
          }
          
          // Read the fit parameters from the function
          gaussFit = hRecoGenRatio[iClosureType][iGenJetPt][iCentrality][iClosureParticle]->GetFunction("gaus");
          gaussMean = gaussFit->GetParameter(1);
          gaussSigma = gaussFit->GetParameter(2);
          gaussMeanError = gaussFit->GetParError(1);
          gaussSigmaError = gaussFit->GetParError(2);
          
          // Fill the histogram with the fit parameters
          hJetPtClosure[iClosureType][iCentrality][iClosureParticle]->SetBinContent(iGenJetPt+1,gaussMean);
          hJetPtClosure[iClosureType][iCentrality][iClosureParticle]->SetBinError(iGenJetPt+1,gaussMeanError);
          hJetPtClosureSigma[iClosureType][iCentrality][iClosureParticle]->SetBinContent(iGenJetPt+1,gaussSigma);
          hJetPtClosureSigma[iClosureType][iCentrality][iClosureParticle]->SetBinError(iGenJetPt+1,gaussSigmaError);
          
        } // Generator level jet pT loop
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
      
      hJetPtClosure[iClosureType][iCentrality][DijetHistograms::knClosureParticleTypes]->SetLineColor(kBlack);
      hJetPtClosure[iClosureType][iCentrality][DijetHistograms::knClosureParticleTypes]->GetYaxis()->SetRangeUser(0.9,1.1);
      drawer->DrawHistogram(hJetPtClosure[iClosureType][iCentrality][DijetHistograms::knClosureParticleTypes],"Gen p_{T} (GeV)","#mu(reco p_{T} / gen p_{T})"," ");
      
      legend = new TLegend(0.54-adders[iClosureType],0.63,0.91,0.9);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      
      if(ppData){
        sprintf(centralityString,", pp");
        sprintf(centralitySaveName,"");
      } else {
        sprintf(centralityString,", Cent:%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
        sprintf(centralitySaveName,"_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
      }
      
      sprintf(namer,"%s jet%s",closureTypeName[iClosureType],centralityString);
      legend->SetHeader(namer);
      legend->AddEntry(hJetPtClosure[iClosureType][iCentrality][DijetHistograms::knClosureParticleTypes],"Inclusive","l");
      
      for(int iClosureParticle = 0; iClosureParticle < DijetHistograms::knClosureParticleTypes; iClosureParticle++){
        hJetPtClosure[iClosureType][iCentrality][iClosureParticle]->SetLineColor(lineColors[iClosureParticle]);
        hJetPtClosure[iClosureType][iCentrality][iClosureParticle]->Draw("same");
        legend->AddEntry(hJetPtClosure[iClosureType][iCentrality][iClosureParticle],particleNames[iClosureParticle],"l");
      } // Closure particle loop (quark/gluon)
      
      legend->Draw();
      
      // Save the figures if selected to do so
      if(saveFigures){
        sprintf(namer,"figures/jetPtClosure%sJet%s.pdf",closureTypeName[iClosureType],centralitySaveName);
        gPad->GetCanvas()->SaveAs(namer);
      }
      
    } // Centrality loop
  } // Closure type loop (leading/subleading/inclusive)
  
}
