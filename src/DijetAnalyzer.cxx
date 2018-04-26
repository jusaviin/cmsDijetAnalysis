// Class for the main analysis algorithms for the leading-subleading jet analysis

// Root includes
#include <TFile.h>
#include <TMath.h>

// Own includes
#include "DijetAnalyzer.h"
#include "ForestReader.h"

using namespace std;

/*
 * Default constructor
 */
DijetAnalyzer::DijetAnalyzer() :
  fFileNames(0),
  fCard(0),
  fHistograms(0),
  fTrackCorrection()
{
  // Default constructor
  fHistograms = new DijetHistograms();
  fHistograms->CreateHistograms();
}

/*
 * Custom constructor
 */
DijetAnalyzer::DijetAnalyzer(std::vector<TString> fileNameVector, ConfigurationCard *newCard) :
  fFileNames(fileNameVector),
  fCard(newCard),
  fHistograms(0)
{
  // Custom constructor
  fHistograms = new DijetHistograms();
  fHistograms->CreateHistograms();
  
  // Find the correct folder for track correction tables based on data type
  Int_t dataType = fCard->Get("DataType");
  if(dataType == ForestReader::kPp || dataType == ForestReader::kPpMC || dataType == ForestReader::kLocalTest){
    fTrackCorrection = new TrkCorr("trackCorrectionTables/TrkCorr_July22_Iterative_pp_eta2p4/");
  } else if (dataType == ForestReader::kPbPb || dataType == ForestReader::kPbPbMC){
    fTrackCorrection = new TrkCorr("trackCorrectionTables/TrkCorr_Jun7_Iterative_PbPb_etaLT2p4/");
  } else {
    fTrackCorrection = new TrkCorr(""); // Bad data type, no corrections initialized
  }
}

/*
 * Copy constructor
 */
DijetAnalyzer::DijetAnalyzer(const DijetAnalyzer& in) :
  fFileNames(in.fFileNames),
  fCard(in.fCard),
  fHistograms(in.fHistograms),
  fTrackCorrection(in.fTrackCorrection)
{
  // Copy constructor
}

/*
 * Assingment operator
 */
DijetAnalyzer& DijetAnalyzer::operator=(const DijetAnalyzer& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fFileNames = in.fFileNames;
  fCard = in.fCard;
  fHistograms = in.fHistograms;
  fTrackCorrection = in.fTrackCorrection;
  
  return *this;
}

/*
 * Destructor
 */
DijetAnalyzer::~DijetAnalyzer(){
  // destructor
  delete fHistograms;
  delete fTrackCorrection;
}

/*
 * Main analysis loop
 */
void DijetAnalyzer::RunAnalysis(){
  
  //****************************************
  //        Event selection cuts
  //****************************************
  
  const Double_t vzCut = fCard->Get("ZVertexCut");  // Event cut vor the z-posiiton of the primary vertex
  
  //****************************************
  //        Dijet selection cuts
  //****************************************
  
  const Double_t jetEtaCut = fCard->Get("JetEtaCut");           // Eta cut around midrapidity
  const Double_t jetSearchEtaCut = fCard->Get("SearchEtaCut");  // Eta cut when searching for a dijet
  const Double_t jetMaximumPtCut = fCard->Get("MaxPtCut");      // Maximum pT accepted for leading jet (and tracks)
  const Double_t leadingJetMinPtCut = fCard->Get("MinPtCut");   // Minimum pT cut for leading jet
  const Double_t subleadingJetMinPtCut = fCard->Get("SubleadingPtCut"); // Minimum pT cut for subleading jet
  const Double_t deltaPhiCut = fCard->Get("DeltaPhiCut");       // DeltaPhi cut for the dijet system
  const Double_t minimumMaxTrackPtFraction = fCard->Get("MinMaxTrackPtFraction");  // Cut for jets consisting only from soft particles
  const Double_t maximumMaxTrackPtFraction = fCard->Get("MaxMaxTrackPtFraction");  // Cut for jets consisting only from one high pT particle
  
  //****************************************
  //        Traxk selection cuts
  //****************************************
  
  const Double_t trackEtaCut = fCard->Get("TrackEtaCut");     // Eta cut around midrapidity
  const Double_t trackMinPtCut = fCard->Get("MinTrackPtCut"); // Minimum pT cut
  const Double_t maxTrackPtRelativeError = fCard->Get("MaxTrackPtRelativeError");   // Maximum relative error for pT
  const Double_t maxTrackDistanceToVertex = fCard->Get("VertexMaxDistance");        // Maximum distance to primary vetrex
  const Double_t calorimeterSignalLimitPt = fCard->Get("CalorimeterSignalLimitPt"); // Require signal in calorimeters for track above this pT
  const Double_t highPtEtFraction = fCard->Get("HighPtEtFraction"); // For high pT tracks, minimum required Et as a fraction of track pT
  const Double_t chi2QualityCut = fCard->Get("Chi2QualityCut");     // Quality cut for track reconstruction
  const Double_t minimumTrackHits = fCard->Get("MinimumTrackHits"); // Quality cut for track hits
  
  //****************************************
  //            All cuts set!
  //****************************************

  // Event variables
  Int_t nEvents = 0;            // Number of events
  Bool_t dijetFound = false;    // Is there a dijet in the event?
  Bool_t twoJetsFound = false;  // Are there two jets in the event?
  Double_t vz = 0;               // Vertex z-position
  Double_t centrality = 0;       // Event centrality
  Int_t hiBin = 0;              // CMS hiBin (centrality * 2)
  
  // Variables for jets
  Double_t dijetAsymmetry = -99;   // Dijet asymmetry
  Double_t leadingPt = 0;          // Leading jet pT
  Double_t subleadingPt = 0;       // Subleading jet pT
  Int_t secondHighestIndex = -1;  // Index of the subleading jet in the event
  Int_t highestIndex = -1;        // Index of the leading jet in the event
  Double_t jetPt = 0;              // pT of the i:th jet in the event
  Double_t dphi = 0;               // deltaPhi for the considered jets
  
  // Variables for tracks
  Double_t trackPt;       // Track pT
  Double_t trackEta;      // Track eta
  Double_t trackPhi;      // Track phi
  Double_t trackR;        // Distance of a track to the current jet
  Double_t trackRMin;     // Minimum distance to a jet
  Double_t trackEt;       // Track transverse energy
  Double_t trackEfficiencyCorrection;  // Efficiency correction for the track
  Double_t deltaPhiTrackLeadingJet;    // DeltaPhi between track and leading jet
  Double_t deltaEtaTrackLeadingJet;    // DeltaEta between track and leading jet
  Double_t deltaPhiTrackSubleadingJet; // DeltaPhi between track and subleading jet
  Double_t deltaEtaTrackSubleadingJet; // DeltaEta between track and subleading jet
  
  // File name helper variables
  TString currentFile;
  
  // Fillers for THnSparses
  Double_t fillerAnyJet[4];
  Double_t fillerDijet[9];
  Double_t fillerTrack[9];
  Double_t fillerTrackPtWeight[6];
  
  // Amount of debugging messages
  Int_t debugLevel = fCard->Get("DebugLevel");
  
  // Loop over files
  for(Int_t iFile = 0; iFile < (Int_t) fFileNames.size(); iFile++) {
    
    // Find the filename
    currentFile = fFileNames.at(iFile);
    
    // Open the file and check that everything goes fine
    TFile *inputFile = new TFile(currentFile); // Create file with new so we can delete it at the end of file loop

    // Check that the file exists
    if(!inputFile){
      cout << "Warning! Could not open the file: " << currentFile.Data() << endl;
      delete inputFile;
      continue;
    }
    
    // Check that the file is not zombie
    if(inputFile->IsZombie()){
      cout << "Warning! The following file is a zombie: " << currentFile.Data() << endl;
      delete inputFile;
      continue;
    }
    
    // Debug message, if wanted
    if(debugLevel > 0) cout << "Reading from file: " << currentFile.Data() << endl;
    
    // If file is good, read the forest from the file
    ForestReader *treeReader = new ForestReader(fCard->Get("DataType")); // Creating this with new and deleting at the end of file loop makes memory leak slower, but does not prevent it completely
    treeReader->ReadForestFromFile(inputFile);  // There seems to be a memory leak here...
    nEvents = treeReader->GetNEvents();
    
    // Event loop
    for(Int_t iEvent = 0; iEvent < nEvents; iEvent++){
      
      // Read the event to memory
      treeReader->GetEvent(iEvent);
      
      // Get vz and centrality information from all events
      vz = treeReader->GetVz();
      centrality = treeReader->GetCentrality();
      hiBin = treeReader->GetHiBin();
      fHistograms->fhVertexZ->Fill(vz);                   // z vertex distribution from all events
      fHistograms->fhEvents->Fill(DijetHistograms::kAll); // All the events looped over
      fHistograms->fhCentrality->Fill(centrality);        // Centrality filled from all events
      
      //  ===== Apply all the event quality cuts =====
      
      // Cut for primary vertex. Only applied for data.
      if(treeReader->GetPrimaryVertexFilterBit() == 0) continue;
      fHistograms->fhEvents->Fill(DijetHistograms::kPrimaryVertex);
      
      // Cut for HB/HE noise. Only applied for data.
      if(treeReader->GetHBHENoiseFilterBit() == 0) continue;
      fHistograms->fhEvents->Fill(DijetHistograms::kHBHENoise);
      
      // Cut for collision event selection. Only applied for PbPb data.
      if(treeReader->GetCollisionEventSelectionFilterBit() == 0) continue;
      fHistograms->fhEvents->Fill(DijetHistograms::kCollisionEventSelection);
      
      // Cut for beam scraping. Only applied for pp data.
      if(treeReader->GetBeamScrapingFilterBit() == 0) continue;
      fHistograms->fhEvents->Fill(DijetHistograms::kBeamScraping);
      
      // Cut for calirimeter jet quality. Only applied for data.
      if(treeReader->GetCaloJetFilterBit() == 0) continue;
      fHistograms->fhEvents->Fill(DijetHistograms::kCaloJet);
      
      // Cut for vertex z-position
      if(TMath::Abs(vz) > vzCut) continue;
      fHistograms->fhEvents->Fill(DijetHistograms::kVzCut);
      
      // ===== Event quality cuts applied =====
      
      // Reset the variables used in dijet finding
      twoJetsFound = false;
      dijetFound = false;
      highestIndex = -1;
      secondHighestIndex = -1;
      leadingPt = 0;
      subleadingPt = 0;
      
      // Search for leading jet and fill histograms for all jets within the eta vut
      for(Int_t jetIndex = 0; jetIndex < treeReader->GetNJets(); jetIndex++) {
        jetPt = treeReader->GetJetPt(jetIndex);
        if(TMath::Abs(treeReader->GetJetEta(jetIndex)) >= jetSearchEtaCut) continue; // Cut for search eta range
        if(minimumMaxTrackPtFraction >= treeReader->GetJetMaxTrackPt(jetIndex)/treeReader->GetJetRawPt(jetIndex)) continue; // Cut for jets with only very low pT particles
        if(maximumMaxTrackPtFraction <= treeReader->GetJetMaxTrackPt(jetIndex)/treeReader->GetJetRawPt(jetIndex)) continue; // Cut for jets where all the pT is taken by one track
        
        // Fill the histogram for all jets within eta range
        if(TMath::Abs(treeReader->GetJetEta(jetIndex)) < jetEtaCut){
          
          // Fill the axes in correct order
          fillerAnyJet[0] = jetPt;                             // Axis 0 = any jet pT
          fillerAnyJet[1] = treeReader->GetJetPhi(jetIndex);   // Axis 1 = any jet phi
          fillerAnyJet[2] = treeReader->GetJetEta(jetIndex);   // Axis 2 = any jet eta
          fillerAnyJet[3] = centrality;                        // Axis 3 = centrality
          fHistograms->fhAnyJet->Fill(fillerAnyJet);           // Fill the data point to histogram
          
        }
        
        if(jetPt <= leadingJetMinPtCut) continue; // Minimum leading jet pT cut
        if(jetPt > leadingPt){
          leadingPt = jetPt;
          highestIndex = jetIndex;
        }
      } // End of search for leading jet loop
      
      // Search for subleading jet
      for(Int_t jetIndex = 0 ; jetIndex < treeReader->GetNJets(); jetIndex++){
        jetPt = treeReader->GetJetPt(jetIndex);
        if(jetIndex == highestIndex) continue; // Do not consider leading particle
        if(TMath::Abs(treeReader->GetJetEta(jetIndex)) >= jetSearchEtaCut) continue; // Cut for search eta range
        if(minimumMaxTrackPtFraction >= treeReader->GetJetMaxTrackPt(jetIndex)/treeReader->GetJetRawPt(jetIndex)) continue; // Cut for jets with only very low pT particles
        if(maximumMaxTrackPtFraction <= treeReader->GetJetMaxTrackPt(jetIndex)/treeReader->GetJetRawPt(jetIndex)) continue; // Cut for jets where all the pT is taken by one track
        if(jetPt <= subleadingJetMinPtCut) continue; // Minimum subleading jet pT cut
        if(jetPt > subleadingPt){
          subleadingPt = jetPt;
          secondHighestIndex = jetIndex;
        }
      }  //End of subleading jet search
      
      // Check if at least two jets were found
      if(highestIndex > -1 && secondHighestIndex > -1) twoJetsFound = true;
      
      // Only apply the dijet cuts for events with at least two jets
      if(twoJetsFound){
        
        dijetFound = true;
        dphi =  treeReader->GetJetPhi(highestIndex) - treeReader->GetJetPhi(secondHighestIndex);
        if(dphi < 0) dphi = -dphi;
        if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
        
        
        if((treeReader->GetJetPt(highestIndex) >= jetMaximumPtCut) ||              // Maximum leading jet pT cut
           (TMath::Abs(treeReader->GetJetEta(highestIndex)) >= jetEtaCut) ||       // Leading jet eta cut
           (TMath::Abs(treeReader->GetJetEta(secondHighestIndex)) >= jetEtaCut)||  // Subleading jet eta cut
           (TMath::Abs(dphi) <= deltaPhiCut)){                                     // DeltaPhi cut
          dijetFound = false;
        }
      } // End of dijet cuts
      
      // If a dijet is found, fill some information to fHistograms
      if(dijetFound){
        fHistograms->fhEvents->Fill(DijetHistograms::kDijet);
        
        // Calculate the asymmetry
        dijetAsymmetry = (treeReader->GetJetPt(highestIndex) - treeReader->GetJetPt(secondHighestIndex))/(treeReader->GetJetPt(highestIndex) + treeReader->GetJetPt(secondHighestIndex));
        
        // Fill the dijet histogram axes in correct order
        fillerDijet[0] = treeReader->GetJetPt(highestIndex);        // Axis 0: Leading jet pT
        fillerDijet[1] = treeReader->GetJetPhi(highestIndex);       // Axis 1: Leading jet phi
        fillerDijet[2] = treeReader->GetJetEta(highestIndex);       // Axis 2: Leading jet eta
        fillerDijet[3] = treeReader->GetJetPt(secondHighestIndex);  // Axis 3: Subleading jet pT
        fillerDijet[4] = treeReader->GetJetPhi(secondHighestIndex); // Axis 4: Subleading jet phi
        fillerDijet[5] = treeReader->GetJetEta(secondHighestIndex); // Axis 5: Subleading jet eta
        fillerDijet[6] = TMath::Abs(dphi);                          // Axis 6: deltaPhi
        fillerDijet[7] = dijetAsymmetry;                            // Axis 7: Asymmetry
        fillerDijet[8] = centrality;                                // Axis 8: Centrality
        fHistograms->fhDijet->Fill(fillerDijet);                    // Fill the data point to dijet histogram
        
        // Correlate jets with tracks in dijet events
        Int_t nTracks = treeReader->GetNTracks();
        for(Int_t iTrack = 0; iTrack <nTracks; iTrack++){
          trackPt = treeReader->GetTrackPt(iTrack);
          trackPhi = treeReader->GetTrackPhi(iTrack);
          trackEta = treeReader->GetTrackEta(iTrack);
          trackEt = (treeReader->GetTrackEnergyEcal(iTrack)+treeReader->GetTrackEnergyHcal(iTrack))/TMath::CosH(trackEta);
          
          //  ==== Apply cuts for tracks and collect information on how much track are cut in each step ====
          fHistograms->fhTrackCuts->Fill(DijetHistograms::kAllTracks);
          
          // Cut for track pT
          if(trackPt <= trackMinPtCut) continue;                     // Minimum pT cut
          if(trackPt >= jetMaximumPtCut) continue;                   // Maximum pT cut (same as for leading jets)
          fHistograms->fhTrackCuts->Fill(DijetHistograms::kPtCuts);
          
          // Cut for track eta
          if(TMath::Abs(trackEta) >= trackEtaCut) continue;          // Eta cut
          fHistograms->fhTrackCuts->Fill(DijetHistograms::kEtaCut);
          
          // Cut for high purity
          if(!treeReader->GetTrackHighPurity(iTrack)) continue;     // High purity cut
          fHistograms->fhTrackCuts->Fill(DijetHistograms::kHighPurity);
          
          // Cut for realtive error for track pT
          if(treeReader->GetTrackPtError(iTrack)/trackPt >= maxTrackPtRelativeError) continue; // Cut for track pT relative error
          fHistograms->fhTrackCuts->Fill(DijetHistograms::kPtError);
          
          // Cut for track distance from primary vertex
          if(treeReader->GetTrackVertexDistanceZ(iTrack)/treeReader->GetTrackVertexDistanceZError(iTrack) >= maxTrackDistanceToVertex) continue; // Mysterious cut about track proximity to vertex in z-direction
          if(treeReader->GetTrackVertexDistanceXY(iTrack)/treeReader->GetTrackVertexDistanceXYError(iTrack) >= maxTrackDistanceToVertex) continue; // Mysterious cut about track proximity to vertex in xy-direction
          fHistograms->fhTrackCuts->Fill(DijetHistograms::kVertexDistance);
          
          // Cut for energy deposition in calorimeters for high pT tracks
          if(!(trackPt < calorimeterSignalLimitPt || (trackEt >= highPtEtFraction*trackPt))) continue;  // For high pT tracks, require signal also in calorimeters
          fHistograms->fhTrackCuts->Fill(DijetHistograms::kCaloSignal);
          
          // Cuts for track reconstruction quality
          if(treeReader->GetTrackChi2(iTrack)/(1.0*treeReader->GetNTrackDegreesOfFreedom(iTrack))/(1.0*treeReader->GetNHitsTrackerLayer(iTrack)) >= chi2QualityCut) continue; // Track reconstruction quality cut
          if(treeReader->GetNHitsTrack(iTrack) < minimumTrackHits) continue; // Cut for minimum number of hits per track
          fHistograms->fhTrackCuts->Fill(DijetHistograms::kReconstructionQuality);
          
          //     ==== Track cuts done ====
          
          // Calculate minimum distance of a track to closest jet. This is needed for track efficiency correction
          trackRMin = 666;   // Initialize the minimum distance to a jet to some very large value
          for(Int_t iJet = 0; iJet < treeReader->GetNJets(); iJet++){
            
            // Require the same jet quality cuts as when searching for dijets
            if(TMath::Abs(treeReader->GetJetEta(iJet)) >= jetSearchEtaCut) continue; // Require jet eta to be in the search range for dijets
            if(treeReader->GetJetPt(iJet) <= subleadingJetMinPtCut) continue; // Only consider jets that pass the pT cut for subleading jets
            if(minimumMaxTrackPtFraction >= treeReader->GetJetMaxTrackPt(iJet)/treeReader->GetJetRawPt(iJet)) continue; // Cut for jets with only very low pT particles
            if(maximumMaxTrackPtFraction <= treeReader->GetJetMaxTrackPt(iJet)/treeReader->GetJetRawPt(iJet)) continue; // Cut for jets where all the pT is taken by one track
            // if(TMath::Abs(chargedSum[k]/rawpt[k]) < 0.01) continue; // Jet quality cut from readme file. TODO: Check if should be applied
            
            // Note: The ACos(Cos(jetPhi-trackPhi)) structure transforms deltaPhi to interval [0,Pi]
            trackR = TMath::Power(treeReader->GetJetEta(iJet)-trackEta,2)+TMath::Power(TMath::ACos(TMath::Cos(treeReader->GetJetPhi(iJet)-trackPhi)),2);
            if(trackRMin*trackRMin>trackR) trackRMin=TMath::Power(trackR,0.5);
            
          } // Loop for calculating Rmin
          
          // Get the track efficiency corrections
          trackEfficiencyCorrection = fTrackCorrection->getTrkCorr(trackPt, trackEta, trackPhi, hiBin, trackRMin);
          
          // Calculate deltaEta and deltaPhi between track and leading and subleading jets
          deltaEtaTrackLeadingJet = treeReader->GetJetEta(highestIndex) - trackEta;
          deltaPhiTrackLeadingJet = treeReader->GetJetPhi(highestIndex) - trackPhi;
          deltaEtaTrackSubleadingJet = treeReader->GetJetEta(secondHighestIndex) - trackEta;
          deltaPhiTrackSubleadingJet = treeReader->GetJetPhi(secondHighestIndex) - trackPhi;
          
          // Tranform deltaPhis to interval [-pi/2,3pi/2]
          while(deltaPhiTrackLeadingJet > (1.5*TMath::Pi())){deltaPhiTrackLeadingJet += -2*TMath::Pi();}
          while(deltaPhiTrackSubleadingJet > (1.5*TMath::Pi())){deltaPhiTrackSubleadingJet += -2*TMath::Pi();}
          while(deltaPhiTrackLeadingJet < (-0.5*TMath::Pi())){deltaPhiTrackLeadingJet += 2*TMath::Pi();}
          while(deltaPhiTrackSubleadingJet < (-0.5*TMath::Pi())){deltaPhiTrackSubleadingJet += 2*TMath::Pi();}
          
          // Fill track histograms
          
          // ===== Nominal tracks =====
          fillerTrack[0] = trackPt;                    // Axis 0: Track pT
          fillerTrack[1] = trackPhi;                   // Axis 1: Track phi
          fillerTrack[2] = trackEta;                   // Axis 2: Track eta
          fillerTrack[3] = deltaPhiTrackLeadingJet;    // Axis 3: DeltaPhi between track and leading jet
          fillerTrack[4] = deltaEtaTrackLeadingJet;    // Axis 4: DeltaEta between track and leading jet
          fillerTrack[5] = deltaPhiTrackSubleadingJet; // Axis 5: DeltaPhi between track and subleading jet
          fillerTrack[6] = deltaEtaTrackSubleadingJet; // Axis 6: DeltaEta between track and subleading jet
          fillerTrack[7] = dijetAsymmetry;             // Axis 7: Dijet asymmetry
          fillerTrack[8] = centrality;                 // Axis 8: Centrality

          // ===== Track for pT weighting =====
          fillerTrackPtWeight[0] = deltaPhiTrackLeadingJet;    // Axis 0: DeltaPhi between track and leading jet
          fillerTrackPtWeight[1] = deltaEtaTrackLeadingJet;    // Axis 1: DeltaEta between track and leading jet
          fillerTrackPtWeight[2] = deltaPhiTrackSubleadingJet; // Axis 2: DeltaPhi between track and subleading jet
          fillerTrackPtWeight[3] = deltaEtaTrackSubleadingJet; // Axis 3: DeltaEta between track and subleading jet
          fillerTrackPtWeight[4] = dijetAsymmetry;             // Axis 4: Dijet asymmetry
          fillerTrackPtWeight[5] = centrality;                 // Axis 5: Centrality
          
          // Fill the histogram for nominal tracks
          fHistograms->fhTrack->Fill(fillerTrack,trackEfficiencyCorrection);
          
          // Fill the histogram for undorrected tracks
          fHistograms->fhTrackUncorrected->Fill(fillerTrack);
          
          // Fill the histogram for pT weighted tracks
          fHistograms->fhTrackPtWeighted->Fill(fillerTrackPtWeight,trackEfficiencyCorrection*trackPt);
          
        } // Loop over tracks
        
      } // Dijet in event
      
    } // Event loop
    
    // Close the input file after the event has been read
    inputFile->Close();
    delete treeReader;
    delete inputFile;
    
  } // File loop
  
}

/*
 * Getter for dijet histograms
 */
DijetHistograms* DijetAnalyzer::GetHistograms() const{
  return fHistograms;
}
