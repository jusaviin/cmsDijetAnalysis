/*
 * Check if given coordinates fall to the hole in acceptance
 */
bool isInHole(double eta, double phi, double minHoleEta, double maxHoleEta, double minHolePhi, double maxHolePhi){
  
  if(eta > minHoleEta && eta < maxHoleEta && phi > minHolePhi && phi < maxHolePhi) return true;
  return false;
  
}

/*
 * Toy simulation to see a shape resulting from a hole in the acceptance for event mixing
 */
void acceptanceHoleMixing(){

  // Configuration for the simulation
  const double maxJetEta = 1.6;
  const double maxTrackEta = 2.4;
  const double minHoleEta = -1.8;
  const double maxHoleEta = -0.9;
  const double minHolePhi = 2.3;
  const double maxHolePhi = 2.7;
  const int nEvents = 10000000;
  const int nTracksInEvent = 1000;
  
  const char* outputFileName = "toyHoleInAcceptance.root";
  
  // Create a random number generator
  TRandom3 *rng = new TRandom3();
  rng->SetSeed(0);
  
  // Create histograms to store the generated distributions
  const double minDeltaPhiJetTrack = -TMath::Pi()/2.0;    // Minimum deltaPhi for jet-track correlations
  const double maxDeltaPhiJetTrack = 3.0*TMath::Pi()/2.0; // Maximum deltaPhi for jet-track correlations
  const int nDeltaPhiBinsJetTrack = 200;                  // Number of deltaPhi bins for jet-track correlations (match the common number in UIC group)
  const double minDeltaEtaJetTrack = -5.0;   // Minimum deltaEta for jet-track correlations
  const double maxDeltaEtaJetTrack = 5.0;    // Maximum deltaEta for jet-track correlations
  const int nDeltaEtaBinsJetTrack = 500;     // Number of deltaEta bins for jet-track correlations (match the common number in UIC group)
  const double minPhi = -TMath::Pi();  // Minimum phi
  const double maxPhi = TMath::Pi();   // Maximum phi
  const int nPhiBins = 64;             // Number of phi bins
  const double minEta = -2.5;          // Minimum eta (current eta cut for tracks = 2.4)
  const double maxEta = 2.5;           // Maximum eta (current eta cut for tracks = 2.4)
  const int nEtaBins = 50;             // Number of eta bins
  
  TH2D *jetTrackDeltaPhiDeltaEta = new TH2D("jetDeltaPhiDeltaEta","jetDeltaPhiDeltaEta",nDeltaPhiBinsJetTrack,minDeltaPhiJetTrack,maxDeltaPhiJetTrack,nDeltaEtaBinsJetTrack,minDeltaEtaJetTrack,maxDeltaEtaJetTrack); jetTrackDeltaPhiDeltaEta->Sumw2();
  TH2D *jetTrackDeltaPhiDeltaEtaNoJetHole =  new TH2D("jetDeltaPhiDeltaEtaNoJetHole","jetDeltaPhiDeltaEtaNoJetHole",nDeltaPhiBinsJetTrack,minDeltaPhiJetTrack,maxDeltaPhiJetTrack,nDeltaEtaBinsJetTrack,minDeltaEtaJetTrack,maxDeltaEtaJetTrack); jetTrackDeltaPhiDeltaEtaNoJetHole->Sumw2();
  TH2D *jetTrackDeltaPhiDeltaEtaNoHole = new TH2D("jetDeltaPhiDeltaEtaNoHole","jetDeltaPhiDeltaEtaNoHole",nDeltaPhiBinsJetTrack,minDeltaPhiJetTrack,maxDeltaPhiJetTrack,nDeltaEtaBinsJetTrack,minDeltaEtaJetTrack,maxDeltaEtaJetTrack); jetTrackDeltaPhiDeltaEtaNoHole->Sumw2();
  TH2D *jetPhiEta = new TH2D("jetPhiEta","jetPhiEta",nPhiBins,minPhi,maxPhi,nEtaBins,minEta,maxEta); jetPhiEta->Sumw2();
  TH2D *jetPhiEtaNoHole = new TH2D("jetPhiEtaNoHole","jetPhiEtaNoHole",nPhiBins,minPhi,maxPhi,nEtaBins,minEta,maxEta); jetPhiEtaNoHole->Sumw2();
  TH2D *trackPhiEta = new TH2D("trackPhiEta","trackPhiEta",nPhiBins,minPhi,maxPhi,nEtaBins,minEta,maxEta); trackPhiEta->Sumw2();
  TH2D *trackPhiEtaNoHole = new TH2D("trackPhiEtaNoHole","trackPhiEtaNoHole",nPhiBins,minPhi,maxPhi,nEtaBins,minEta,maxEta); trackPhiEtaNoHole->Sumw2();
  
  // Start the event loop
  double jetEta, jetPhi, jetEtaNoHole, jetPhiNoHole;
  double trackEta, trackPhi, trackEtaNoHole, trackPhiNoHole;
  double deltaEta, deltaPhi, deltaEtaNoJetHole, deltaPhiNoJetHole, deltaEtaNoHole, deltaPhiNoHole;
  for(int iEvents = 0; iEvents < nEvents; iEvents++){
    
    if(iEvents % 10000 == 0){
      cout << "Doing event: " << iEvents << " / " << nEvents << endl;
    }
    
    // For jets without hole, just get the values from a uniform distribution
    jetEtaNoHole = rng->Uniform(-maxJetEta,maxJetEta);
    jetPhiNoHole = rng->Uniform(-TMath::Pi(),TMath::Pi());
    jetEta = jetEtaNoHole;
    jetPhi = jetPhiNoHole;
    
    // With the hole, do not accept events inside the hole
    while (isInHole(jetEta,jetPhi,minHoleEta,maxHoleEta,minHolePhi,maxHolePhi)) {
      jetEta = rng->Uniform(-maxJetEta,maxJetEta);
      jetPhi = rng->Uniform(-TMath::Pi(),TMath::Pi());
    }
    
    // Fill the jet histograms
    jetPhiEta->Fill(jetPhi,jetEta);
    jetPhiEtaNoHole->Fill(jetPhiNoHole,jetEtaNoHole);
    
    // Loop over tracks in an event
    for(int iTrack = 0; iTrack < nTracksInEvent; iTrack++){
      
      trackEtaNoHole = rng->Uniform(-maxTrackEta,maxTrackEta);
      trackPhiNoHole = rng->Uniform(-TMath::Pi(),TMath::Pi());
      trackEta = trackEtaNoHole;
      trackPhi = trackPhiNoHole;
      
      // Discard all the tracks that fall inside the hole
      while (isInHole(trackEta,trackPhi,minHoleEta,maxHoleEta,minHolePhi,maxHolePhi)){
        trackEta = rng->Uniform(-maxTrackEta,maxTrackEta);
        trackPhi = rng->Uniform(-TMath::Pi(),TMath::Pi());
      }
      
      // Calculate the jet-track correlation
      deltaEta = jetEta - trackEta;
      deltaPhi = jetPhi - trackPhi;
      deltaEtaNoJetHole = jetEtaNoHole - trackEta;
      deltaPhiNoJetHole = jetPhiNoHole - trackPhi;
      deltaEtaNoHole = jetEtaNoHole - trackEtaNoHole;
      deltaPhiNoHole = jetPhiNoHole - trackPhiNoHole;
      
      // Transform deltaPhis to interval [-pi/2,3pi/2]
      while(deltaPhi > (1.5*TMath::Pi())){deltaPhi += -2*TMath::Pi();}
      while(deltaPhiNoJetHole > (1.5*TMath::Pi())){deltaPhiNoJetHole += -2*TMath::Pi();}
      while(deltaPhiNoHole > (1.5*TMath::Pi())){deltaPhiNoHole += -2*TMath::Pi();}
      while(deltaPhi < (-0.5*TMath::Pi())){deltaPhi += 2*TMath::Pi();}
      while(deltaPhiNoJetHole < (-0.5*TMath::Pi())){deltaPhiNoJetHole += 2*TMath::Pi();}
      while(deltaPhiNoHole < (-0.5*TMath::Pi())){deltaPhiNoHole += 2*TMath::Pi();}
      
      // Fill the track and correlation histograms
      trackPhiEta->Fill(trackPhi,trackEta);
      trackPhiEtaNoHole->Fill(trackPhiNoHole,trackEtaNoHole);
      jetTrackDeltaPhiDeltaEta->Fill(deltaPhi,deltaEta);
      jetTrackDeltaPhiDeltaEtaNoJetHole->Fill(deltaPhiNoJetHole,deltaEtaNoJetHole);
      jetTrackDeltaPhiDeltaEtaNoHole->Fill(deltaPhiNoHole,deltaEtaNoHole);
      
    } // Track loop
    
  } // Event loop
  
  // Write the histograms to the output file
  TFile *outfile = new TFile(outputFileName, "RECREATE");
  jetPhiEta->Write();
  jetPhiEtaNoHole->Write();
  trackPhiEta->Write();
  trackPhiEtaNoHole->Write();
  jetTrackDeltaPhiDeltaEta->Write();
  jetTrackDeltaPhiDeltaEtaNoJetHole->Write();
  jetTrackDeltaPhiDeltaEtaNoHole->Write();
  outfile->Close();
  
}
