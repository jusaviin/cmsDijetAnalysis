/*
 * Toy simulation to see a shape resulting from a hole in the acceptance for event mixing
 */
void toySmearing(){

  // Open the data file containing the generator level jet spectra from generator level dijets
  TFile *inputFile = TFile::Open("../data/PbPbMC_GenGen_skims_pfJets_noCorrelations_matchedJetsWithFlavor_processed_2019-02-04.root");
  
  // Configuration for the simulation
  const int nSmearedJets = 1000;
  const double smearingSigma = 0.2;  // Gaussian sigma in percentage of jet pT
  const int nCentralityBins = 4;
  
  const char* outputFileName = "toySmearedGenJetPtTeest.root";
  
  // Create a random number generator
  TRandom3 *rng = new TRandom3();
  rng->SetSeed(0);
  gRandom->SetSeed(0);
  
  // Create histograms to store the generated distributions
  const double minJetPt = 0;    // Minimum pT for the jets
  const double maxJetPt = 500;  // Maximum pT for the jets
  const int nBinsJetPt = 100;   // Number of jet pT bins
  
  TH1D *smearedJetPt[nCentralityBins];
  TH1D *genJetPt[nCentralityBins];
  char namer[100];
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    // Create histograms for smeared jet pT
    sprintf(namer,"smearedJetPt_C%d",iCentrality);
    smearedJetPt[iCentrality] = new TH1D(namer,namer,nBinsJetPt,minJetPt,maxJetPt); smearedJetPt[iCentrality]->Sumw2();
    
    // Read histograms for gen jet pT
    sprintf(namer,"subleadingJet/subleadingJetPt_C%d",iCentrality);
    genJetPt[iCentrality] = (TH1D*) inputFile->Get(namer);
    
  }
  
  // Start the event loop
  double jetPt, jetPtAfterSmearing;
  
  for(int iJets = 0; iJets < nSmearedJets; iJets++){
    
    if(iJets % 10000 == 0){
      cout << "Doing jet: " << iJets << " / " << nSmearedJets << endl;
    }
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
      // Get a random jet pT from the generator level jet pT distribution
      jetPt = genJetPt[iCentrality]->GetRandom();
    
      // Do a gaussian smearing for the jet pT
      jetPtAfterSmearing = rng->Gaus(jetPt,jetPt*smearingSigma);
      
      // Fill the newly obtained jet pT to the smeared distribution
      smearedJetPt[iCentrality]->Fill(jetPtAfterSmearing);
      
    } // Centrality loop
    
  } // Jet loop
  
  // Write the histograms to the output file
  TFile *outfile = new TFile(outputFileName, "RECREATE");
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    smearedJetPt[iCentrality]->Write();
  }
  outfile->Close();
  
}
