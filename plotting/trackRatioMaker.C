#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)

/*
 * Macro for finding ratio of uncorrected minimum bias tracks from 2018 and 2015 to be used in a temporary correction while waiting for final tracking corrections to arrive
 */
void trackRatioMaker(){
  
  // Load the minimum bias files
  TFile *minBiasFile2018 = TFile::Open("data/dijetPbPb2018_akPu4CaloJets_minBiasTrackScan_processed_2019-09-20.root");
  TFile *minBiasFile2015 = TFile::Open("data/dijetPbPb2015_minBiasTrackScan_processed_2019-09-21.root");

  // Create histogram managers to read the histograms from the files
  DijetHistogramManager *histogramProvider2018 = new DijetHistogramManager(minBiasFile2018);
  DijetHistogramManager *histogramProvider2015 = new DijetHistogramManager(minBiasFile2015);
  
  // Load the uncorrected track histograms to the histogram managers
  histogramProvider2018->SetLoadEventInformation(true);
  histogramProvider2018->SetLoadInclusiveTracksUncorrected(true);
  histogramProvider2018->SetLoad2DHistograms(true);
  histogramProvider2018->LoadProcessedHistograms();

  histogramProvider2015->SetLoadEventInformation(true);
  histogramProvider2015->SetLoadInclusiveTracksUncorrected(true);
  histogramProvider2015->SetLoad2DHistograms(true);
  histogramProvider2015->LoadProcessedHistograms();
  
  // Read the two dimensional eta-phi maps in all the analysis centrality and track pT bins
  const int nCentralityBins = histogramProvider2018->GetNCentralityBins();
  const int nTrackPtBins = histogramProvider2018->GetNTrackPtBins();
  TH2D *minBiasTracks2018[nCentralityBins][nTrackPtBins];
  TH2D *minBiasTracks2015[nCentralityBins][nTrackPtBins];
  TH2D *minBiasTrackRatio[nCentralityBins][nTrackPtBins];
  TH1D *eventCount2018;
  TH1D *eventCount2015;
  int nEvents2018;
  int nEvents2015;
  
  eventCount2018 = histogramProvider2018->GetHistogramCentralityWeighted();
  eventCount2015 = histogramProvider2015->GetHistogramCentralityWeighted();
  
  double binBorders2018[] = {0,10,30,50,90};
  double binBorders2015[] = {0,10,30,50,90};
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    nEvents2018 = eventCount2018->Integral(eventCount2018->FindBin(binBorders2018[iCentrality]+0.01), eventCount2018->FindBin(binBorders2018[iCentrality+1]-0.01));
    nEvents2015 = eventCount2015->Integral(eventCount2018->FindBin(binBorders2015[iCentrality]+0.01), eventCount2015->FindBin(binBorders2015[iCentrality+1]-0.01));
    
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      
      minBiasTracks2018[iCentrality][iTrackPt] = histogramProvider2018->GetHistogramTrackEtaPhi(DijetHistogramManager::kUncorrectedInclusiveTrack, DijetHistogramManager::kSameEvent, iCentrality, iTrackPt);
      
      minBiasTracks2015[iCentrality][iTrackPt] = histogramProvider2015->GetHistogramTrackEtaPhi(DijetHistogramManager::kUncorrectedInclusiveTrack, DijetHistogramManager::kSameEvent, iCentrality, iTrackPt);
      
      // Normalize both distribution by the weighted number of events in a centrality bin
      minBiasTracks2018[iCentrality][iTrackPt]->Scale(1.0/nEvents2018);
      minBiasTracks2015[iCentrality][iTrackPt]->Scale(1.0/nEvents2015);
      
      // Take the ratio of the distributions
      minBiasTrackRatio[iCentrality][iTrackPt] = (TH2D*) minBiasTracks2015[iCentrality][iTrackPt]->Clone(Form("minimumBiasTrackRatio_C%dT%d", iCentrality, iTrackPt));
      
      minBiasTrackRatio[iCentrality][iTrackPt]->Divide(minBiasTracks2018[iCentrality][iTrackPt]);
      
    } // Track pT loop
  } // Centrality loop
  
  // For some bad statistics bins, reset everything to one
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = nTrackPtBins-2; iTrackPt < nTrackPtBins; iTrackPt++){
      if((iTrackPt == nTrackPtBins-2) && iCentrality != nCentralityBins-1) continue;
      
      for(int iPhi = 1; iPhi <= minBiasTrackRatio[iCentrality][iTrackPt]->GetNbinsX(); iPhi++){
        for(int iEta = 1; iEta <= minBiasTrackRatio[iCentrality][iTrackPt]->GetNbinsY(); iEta++){
          minBiasTrackRatio[iCentrality][iTrackPt]->SetBinContent(iPhi,iEta,1);
        }
      }
    }
  }
  
  // If there are zeros left in the ratio histograms, set them to one
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iPhi = 1; iPhi <= minBiasTrackRatio[iCentrality][iTrackPt]->GetNbinsX(); iPhi++){
        for(int iEta = 1; iEta <= minBiasTrackRatio[iCentrality][iTrackPt]->GetNbinsY(); iEta++){
          
          if(minBiasTrackRatio[iCentrality][iTrackPt]->GetBinContent(iPhi,iEta) == 0){
            minBiasTrackRatio[iCentrality][iTrackPt]->SetBinContent(iPhi,iEta,1);
          }
        }
      }
    }
  }
  
  // Write the ratio histogram to a file
  TFile *outputFile = new TFile("minimumBiasTrackRatio.root","RECREATE");
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      minBiasTrackRatio[iCentrality][iTrackPt]->Write();
    }
  }
  
  // Close the output file
  outputFile->Close();
}
