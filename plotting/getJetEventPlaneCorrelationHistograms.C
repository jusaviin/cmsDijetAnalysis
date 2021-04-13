void getJetEventPlaneCorrelationHistograms(){

  // Open the data file
  TFile *inputFile = TFile::Open("data/dijetPbPb2018_akPfCsJet_onlyJets_eventPlane_noMixing_2021-04-12.root");
  
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
  const int nQvectorBins = 4;
  TH1D *jetEventPlane[2][nQvectorBins][nCentralityBins];
  
  const double lowQvectorBin[nQvectorBins] = {1,5,6,7};
  const double highQvectorBin = 9;
  
  const char* eventString[2] = {"jetEventPlaneDeltaPhiForwardRap","jetEventPlaneDeltaPhiMidRap"};
  
  char newName[200];
  for(int iType = 0; iType < 2; iType++){
    for(int iQvector = 0; iQvector < 1; iQvector++){  // nQvectorBins
      histogramArray[iType]->GetAxis(1)->SetRange(lowQvectorBin[iQvector],highQvectorBin);
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        histogramArray[iType]->GetAxis(2)->SetRange(iCentrality+2,iCentrality+2);
        sprintf(newName,"%s_Q%dC%d",eventString[iType], iQvector, iCentrality);
        
        jetEventPlane[iType][iQvector][iCentrality] = (TH1D*) histogramArray[iType]->Projection(0);
        jetEventPlane[iType][iQvector][iCentrality]->SetName(newName);
        
      }
    }
  }
  
  // Save the histogram to a file
  TFile *outputFile = new TFile("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbData2018_akPfCsJet_wtaAxis_2021-04-12.root","UPDATE");
  for(int iType = 0; iType < 2; iType++){
    for(int iQvector = 0; iQvector < 1; iQvector++){ // nQvectorBins
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        jetEventPlane[iType][iQvector][iCentrality]->Write("",TObject::kOverwrite);
      }
    }
  }
  
  outputFile->Close();
}
