/*
 * Toy simulation to see a shape resulting from a hole in the acceptance for event mixing
 */
void distributionMixing(){

  // Open the data file and read the jet phi from same events
  TFile *inputFile = TFile::Open("../data/PbPbMC_RecoGen_skims_pfJets_noUncorr_5eveImprovedMix_subeNon0_fixedCentality_processed_2019-02-15.root");
  
  // Configuration for the simulation
  const int nEvents = 1000000;
  const int nTracksInEvent = 1000;
  
  const char* outputFileName = "mixingFromSubeNon0Distribution1M.root";
  
  // Create a random number generator
  gRandom->SetSeed(0);
  
  // Create histograms to store the generated distributions
  const double minDeltaPhiJetTrack = -TMath::Pi()/2.0;    // Minimum deltaPhi for jet-track correlations
  const double maxDeltaPhiJetTrack = 3.0*TMath::Pi()/2.0; // Maximum deltaPhi for jet-track correlations
  const int nDeltaPhiBinsJetTrack = 200;                  // Number of deltaPhi bins for jet-track correlations (match the common number in UIC group)
  const double minDeltaEtaJetTrack = -5.0;   // Minimum deltaEta for jet-track correlations
  const double maxDeltaEtaJetTrack = 5.0;    // Maximum deltaEta for jet-track correlations
  const int nDeltaEtaBinsJetTrack = 500;     // Number of deltaEta bins for jet-track correlations (match the common number in UIC group)
  
  TH2D *jetTrackDeltaPhiDeltaEta = new TH2D("jetDeltaPhiDeltaEta","jetDeltaPhiDeltaEta",nDeltaPhiBinsJetTrack,minDeltaPhiJetTrack,maxDeltaPhiJetTrack,nDeltaEtaBinsJetTrack,minDeltaEtaJetTrack,maxDeltaEtaJetTrack); jetTrackDeltaPhiDeltaEta->Sumw2();
  TH2D *jetTrackDeltaPhiDeltaEtaSame = new TH2D("jetDeltaPhiDeltaEtaSame","jetDeltaPhiDeltaEtaSame",nDeltaPhiBinsJetTrack,minDeltaPhiJetTrack,maxDeltaPhiJetTrack,nDeltaEtaBinsJetTrack,minDeltaEtaJetTrack,maxDeltaEtaJetTrack); jetTrackDeltaPhiDeltaEtaSame->Sumw2();

  TH2D *jetPhiEta = (TH2D*) inputFile->Get("leadingJet/leadingJetEtaPhi_C0");
  TH2D *trackPhiEtaSame = (TH2D*) inputFile->Get("track/trackEtaPhi_SameEvent_C0T1");
  TH2D *trackPhiEta = (TH2D*) inputFile->Get("track/trackEtaPhi_MixedEvent_C0T1");
  
  // Start the event loop
  double jetEta, jetPhi;
  double trackEta, trackPhi, trackEtaSame, trackPhiSame;
  double deltaEta, deltaPhi, deltaEtaSame, deltaPhiSame;
  for(int iEvents = 0; iEvents < nEvents; iEvents++){
    
    if(iEvents % 10000 == 0){
      cout << "Doing event: " << iEvents << " / " << nEvents << endl;
    }
    
    // For jets without hole, just get the values from a uniform distribution
    jetPhiEta->GetRandom2(jetPhi,jetEta);
    
    // Loop over tracks in an event
    for(int iTrack = 0; iTrack < nTracksInEvent; iTrack++){
      
      trackPhiEta->GetRandom2(trackPhi,trackEta);
      trackPhiEtaSame->GetRandom2(trackPhiSame,trackEtaSame);
      
      
      // Calculate the jet-track correlation
      deltaEta = jetEta - trackEta;
      deltaPhi = jetPhi - trackPhi;
      deltaEtaSame = jetEta - trackEtaSame;
      deltaPhiSame = jetPhi - trackPhiSame;
      
      // Transform deltaPhis to interval [-pi/2,3pi/2]
      while(deltaPhi > (1.5*TMath::Pi())){deltaPhi += -2*TMath::Pi();}
      while(deltaPhiSame > (1.5*TMath::Pi())){deltaPhiSame += -2*TMath::Pi();}
      while(deltaPhi < (-0.5*TMath::Pi())){deltaPhi += 2*TMath::Pi();}
      while(deltaPhiSame < (-0.5*TMath::Pi())){deltaPhiSame += 2*TMath::Pi();}
      
      // Fill the track and correlation histograms
      jetTrackDeltaPhiDeltaEta->Fill(deltaPhi,deltaEta);
      jetTrackDeltaPhiDeltaEtaSame->Fill(deltaPhiSame,deltaEtaSame);
      
    } // Track loop
    
  } // Event loop
  
  // Write the histograms to the output file
  TFile *outfile = new TFile(outputFileName, "RECREATE");
  jetTrackDeltaPhiDeltaEta->Write();
  jetTrackDeltaPhiDeltaEtaSame->Write();
  outfile->Close();
  
}
