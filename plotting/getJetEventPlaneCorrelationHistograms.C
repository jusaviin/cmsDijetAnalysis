void getJetEventPlaneCorrelationHistograms(){

  // Open the data file
  TFile *inputFile = TFile::Open("data/PbPbMC2018_RecoGen_akCaloJet_onlyJets_multEventPlane_multWeight_jetEta1v3_2022-01-25.root");
  
  // Read the histogram with the given name from the file
  THnSparseD *histogramArray[2];
  histogramArray[0] = (THnSparseD*) inputFile->Get("jetEventPlaneForwardRap");
  histogramArray[1] = (THnSparseD*) inputFile->Get("jetEventPlaneMidRap");
  
  // If cannot find histogram, inform that it could not be found and return null
  if(histogramArray[0] == nullptr || histogramArray[1] == nullptr){
    cout << "Could not find histograms. Will not compute." << endl;
    return;
  }
  
  // Apply the restrictions in the set of axes
  const int nCentralityBins = 4;
  const int nQvectorBins = 9;
  TH1D *jetEventPlane[2][nQvectorBins][nCentralityBins];
  
  const double lowQvectorBin[nQvectorBins] = {1,5,6,7};
  const double highQvectorBin = 9;
  
  const char* eventString[2] = {"jetEventPlaneDeltaPhiManual","jetEventPlaneDeltaPhiDiff"};
  
  char newName[200];
  // Regular centrality based binning
//  for(int iType = 0; iType < 2; iType++){
//    for(int iQvector = 0; iQvector < nQvectorBins; iQvector++){  // nQvectorBins
//      histogramArray[iType]->GetAxis(1)->SetRange(iQvector+1,iQvector+1);
//      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
//        histogramArray[iType]->GetAxis(2)->SetRange(iCentrality+2,iCentrality+2);
//        sprintf(newName,"%s_Q%dC%d",eventString[iType], iQvector, iCentrality);
//
//        jetEventPlane[iType][iQvector][iCentrality] = (TH1D*) histogramArray[iType]->Projection(0);
//        jetEventPlane[iType][iQvector][iCentrality]->SetName(newName);
//
//      }
//    }
//  }
  
  // Binning in multiplicity bins
    for(int iType = 0; iType < 2; iType++){
      for(int iQvector = 0; iQvector < nQvectorBins; iQvector++){  // nQvectorBins
        histogramArray[iType]->GetAxis(1)->SetRange(iQvector+1,iQvector+1);
          sprintf(newName,"%s_M%d",eventString[iType], iQvector);
  
          jetEventPlane[iType][iQvector][0] = (TH1D*) histogramArray[iType]->Projection(0);
          jetEventPlane[iType][iQvector][0]->SetName(newName);
      }
    }
  
  // Save the histogram to a file
  TFile *outputFile = new TFile("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_caloJets_multiplicityBins_manualEventPlane_jetEta1v3_2022-01-25.root","UPDATE");
  for(int iType = 0; iType < 2; iType++){
    for(int iQvector = 0; iQvector < nQvectorBins; iQvector++){ // nQvectorBins
      for(int iCentrality = 0; iCentrality < 1; iCentrality++){
        jetEventPlane[iType][iQvector][iCentrality]->Write("",TObject::kOverwrite);
      }
    }
  }
  
  outputFile->Close();
}
