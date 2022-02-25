void getJetEventPlaneCorrelationHistograms(){

  // Open the data file
  TFile *inputFile = TFile::Open("data/PbPbMC2018_RecoGen_akPfCsJet_onlyRegular_multWeight_subeNon0_fakeJetV2p8_jetEta1v6_2022-02-21.root");
  
  // Configuration
  const int nEventPlaneOrder = 3;
  
  // Read the histogram with the given name from the file
  THnSparseD *jetEventPlaneArray[nEventPlaneOrder];
  THnSparseD *jetEventPlaneDifferenceArray[nEventPlaneOrder];
  for(int iOrder = 0; iOrder < nEventPlaneOrder; iOrder++){
    jetEventPlaneArray[iOrder] = (THnSparseD*) inputFile->Get(Form("jetEventPlaneOrder%d", iOrder+2));
    jetEventPlaneDifferenceArray[iOrder] = (THnSparseD*) inputFile->Get(Form("jetEventPlaneDifferenceOrder%d", iOrder+2));
  }
  
  // If cannot find histogram, inform that it could not be found and return null
  for(int iOrder = 0; iOrder < nEventPlaneOrder; iOrder++){
    if(jetEventPlaneArray[iOrder] == nullptr || jetEventPlaneDifferenceArray[iOrder] == nullptr){
      cout << "Could not find histograms of order " << iOrder << ". Will not compute." << endl;
      return;
    }
  }
  
  // Apply the restrictions in the set of axes
  const int nCentralityBins = 4;
  const int nMultiplicityBins = 20;
  TH1D *jetEventPlaneCentrality[nEventPlaneOrder][nCentralityBins];
  TH1D *jetEventPlaneMultiplicity[nEventPlaneOrder][nMultiplicityBins];
  TH1D *jetEventPlaneDifferenceCentrality[nEventPlaneOrder][nCentralityBins];
  TH1D *jetEventPlaneDifferenceMultiplicity[nEventPlaneOrder][nMultiplicityBins];
  
  char newName[200];
  
  // Regular centrality based binning
  for(int iOrder = 0; iOrder < nEventPlaneOrder; iOrder++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      jetEventPlaneArray[iOrder]->GetAxis(2)->SetRange(iCentrality+2,iCentrality+2);
      sprintf(newName,"jetEventPlaneOrder%d_C%d", iOrder+2, iCentrality);
      
      jetEventPlaneCentrality[iOrder][iCentrality] = (TH1D*) jetEventPlaneArray[iOrder]->Projection(0);
      jetEventPlaneCentrality[iOrder][iCentrality]->SetName(newName);
      
      jetEventPlaneDifferenceArray[iOrder]->GetAxis(2)->SetRange(iCentrality+2,iCentrality+2);
      sprintf(newName,"jetEventPlaneDifferenceOrder%d_C%d", iOrder+2, iCentrality);
      
      jetEventPlaneDifferenceCentrality[iOrder][iCentrality] = (TH1D*) jetEventPlaneDifferenceArray[iOrder]->Projection(0);
      jetEventPlaneDifferenceCentrality[iOrder][iCentrality]->SetName(newName);
      
    }
  }
  
  // Reset the range in histogram arrays
  for(int iOrder = 0; iOrder < nEventPlaneOrder; iOrder++){
    jetEventPlaneArray[iOrder]->GetAxis(2)->SetRange(0,0);
    jetEventPlaneDifferenceArray[iOrder]->GetAxis(2)->SetRange(0,0);
  }
  
  
  // Binning in multiplicity bins
  for(int iOrder = 0; iOrder < nEventPlaneOrder; iOrder++){
    for(int iMultiplicity = 0; iMultiplicity < nMultiplicityBins; iMultiplicity++){
      jetEventPlaneArray[iOrder]->GetAxis(1)->SetRange(iMultiplicity+1,iMultiplicity+1);
      sprintf(newName,"jetEventPlaneOrder%d_M%d", iOrder+2, iMultiplicity);
      
      jetEventPlaneMultiplicity[iOrder][iMultiplicity] = (TH1D*) jetEventPlaneArray[iOrder]->Projection(0);
      jetEventPlaneMultiplicity[iOrder][iMultiplicity]->SetName(newName);
      
      jetEventPlaneDifferenceArray[iOrder]->GetAxis(1)->SetRange(iMultiplicity+1,iMultiplicity+1);
      sprintf(newName,"jetEventPlaneDifferenceOrder%d_M%d", iOrder+2, iMultiplicity);
      
      jetEventPlaneDifferenceMultiplicity[iOrder][iMultiplicity] = (TH1D*) jetEventPlaneDifferenceArray[iOrder]->Projection(0);
      jetEventPlaneDifferenceMultiplicity[iOrder][iMultiplicity]->SetName(newName);
    }
  }
  
  // Save the histogram to a file
  TFile *outputFile = new TFile("eventPlaneCorrelation/jetEventPlaneDeltaPhi_PbPbMC2018_pfCsJets_multWeight_fakeJetV2p8_jetEta1v6_2022-02-24.root","UPDATE");
  for(int iOrder = 0; iOrder < nEventPlaneOrder; iOrder++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      jetEventPlaneCentrality[iOrder][iCentrality]->Write("",TObject::kOverwrite);
      jetEventPlaneDifferenceCentrality[iOrder][iCentrality]->Write("",TObject::kOverwrite);
    }
    
    for(int iMultiplicity = 0; iMultiplicity < nMultiplicityBins; iMultiplicity++){
      jetEventPlaneMultiplicity[iOrder][iMultiplicity]->Write("",TObject::kOverwrite);
      jetEventPlaneDifferenceMultiplicity[iOrder][iMultiplicity]->Write("",TObject::kOverwrite);
    }
  }
  
  outputFile->Close();
}
