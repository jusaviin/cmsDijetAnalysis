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
  int dataType = fCard->Get("DataType");
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
  
  TFile *inputFile;
  ForestReader *treeReader = new ForestReader(fCard->Get("DataType"));
  
  //****************************************
  //        Event selection cuts
  //****************************************
  
  const double vzCut = fCard->Get("ZVertexCut");  // Event cut vor the z-posiiton of the primary vertex
  
  //****************************************
  //        Dijet selection cuts
  //****************************************
  
  const double jetEtaCut = fCard->Get("JetEtaCut");           // Eta cut around midrapidity
  const double jetSearchEtaCut = fCard->Get("SearchEtaCut");  // Eta cut when searching for a dijet
  const double jetMaximumPtCut = fCard->Get("MaxPtCut");      // Maximum pT accepted for leading jet (and tracks)
  const double leadingJetMinPtCut = fCard->Get("MinPtCut");   // Minimum pT cut for leading jet
  const double subleadingJetMinPtCut = fCard->Get("SubleadingPtCut"); // Minimum pT cut for subleading jet
  const double deltaPhiCut = fCard->Get("DeltaPhiCut");       // DeltaPhi cut for the dijet system
  const double minimumMaxTrackPtFraction = fCard->Get("MinMaxTrackPtFraction");  // Cut for jets consisting only from soft particles
  const double maximumMaxTrackPtFraction = fCard->Get("MaxMaxTrackPtFraction");  // Cut for jets consisting only from one high pT particle
  
  //****************************************
  //        Traxk selection cuts
  //****************************************
  
  const double trackEtaCut = fCard->Get("TrackEtaCut");     // Eta cut around midrapidity
  const double trackMinPtCut = fCard->Get("MinTrackPtCut"); // Minimum pT cut
  const double maxTrackPtRelativeError = fCard->Get("MaxTrackPtRelativeError");   // Maximum relative error for pT
  const double maxTrackDistanceToVertex = fCard->Get("VertexMaxDistance");        // Maximum distance to primary vetrex
  const double calorimeterSignalLimitPt = fCard->Get("CalorimeterSignalLimitPt"); // Require signal in calorimeters for track above this pT
  const double highPtEtFraction = fCard->Get("HighPtEtFraction"); // For high pT tracks, minimum required Et as a fraction of track pT
  const double chi2QualityCut = fCard->Get("Chi2QualityCut");     // Quality cut for track reconstruction
  const double minimumTrackHits = fCard->Get("MinimumTrackHits"); // Quality cut for track hits
  
  //****************************************
  //            All cuts set!
  //****************************************

  // Event variables
  int nEvents = 0;            // Number of events
  bool dijetFound = false;    // Is there a dijet in the event?
  bool twoJetsFound = false;  // Are there two jets in the event?
  double vz = 0;              // Vertex z-position
  double centrality = 0;      // Event centrality
  int hiBin = 0;              // CMS hiBin (centrality * 2)
  
  // Variables for jets
  double Aj = -99;              // Jet asymmetry
  double leadingPt = 0;         // Leading jet pT
  double subleadingPt = 0;      // Subleading jet pT
  int secondHighestIndex = -1;  // Index of the subleading jet in the event
  int highestIndex = -1;        // Index of the leading jet in the event
  double jetPt = 0;             // pT of the i:th jet in the event
  double dphi = 0;              // deltaPhi for the considered jets
  
  // Variables for tracks
  double trackPt;       // Track pT
  double trackEta;      // Track eta
  double trackPhi;      // Track phi
  double trackR;        // Distance of a track to the current jet
  double trackRMin;     // Minimum distance to a jet
  double trackEt;       // Track transverse energy
  double trackEfficiencyCorrection; // Efficiency correction for the track
  
  // File name helper variables
  TString currentFile;
  
  // Fillers for THnSparses
  double fillerAnyJet[4];
  double fillerDijet[9];
  
  // Amount of debugging messages
  int debugLevel = fCard->Get("DebugLevel");
  
  // Loop over files
  for(int iFile = 0; iFile < (int) fFileNames.size(); iFile++) {
    
    // Find the filename
    currentFile = fFileNames.at(iFile);
    
    // Open the file and check that everything goes fine
    inputFile = TFile::Open(currentFile);

    // Check that the file exists
    if(!inputFile){
      cout << "Warning! Could not open the file: " << currentFile.Data() << endl;
      continue;
    }
    
    // Check that the file is not zombie
    if(inputFile->IsZombie()){
      cout << "Warning! The following file is a zombie: " << currentFile.Data() << endl;
      continue;
    }
    
    // Debug message, if wanted
    if(debugLevel > 0) cout << "Reading from file: " << currentFile.Data() << endl;
    
    // If file is good, read the forest from the file
    treeReader->ReadForestFromFile(inputFile);
    nEvents = treeReader->GetNEvents();
    
    // Event loop
    for(int iEvent = 0; iEvent < nEvents; iEvent++){
      
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
      for(int jetIndex = 0; jetIndex < treeReader->GetNJets(); jetIndex++) {
        jetPt = treeReader->GetJetPt(jetIndex);
        if(TMath::Abs(treeReader->GetJetEta(jetIndex)) >= jetSearchEtaCut) continue; // Cut for search eta range
        if(minimumMaxTrackPtFraction > treeReader->GetJetMaxTrackPt(jetIndex)/treeReader->GetJetRawPt(jetIndex)) continue; // Cut for jets with only very low pT particles
        if(maximumMaxTrackPtFraction < treeReader->GetJetMaxTrackPt(jetIndex)/treeReader->GetJetRawPt(jetIndex)) continue; // Cut for jets where all the pT is taken by one track
        
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
      for(int jetIndex = 0 ; jetIndex < treeReader->GetNJets(); jetIndex++){
        jetPt = treeReader->GetJetPt(jetIndex);
        if(jetIndex == highestIndex) continue; // Do not consider leading particle
        if(TMath::Abs(treeReader->GetJetEta(jetIndex)) >= jetSearchEtaCut) continue; // Cut for search eta range
        if(minimumMaxTrackPtFraction > treeReader->GetJetMaxTrackPt(jetIndex)/treeReader->GetJetRawPt(jetIndex)) continue; // Cut for jets with only very low pT particles
        if(maximumMaxTrackPtFraction < treeReader->GetJetMaxTrackPt(jetIndex)/treeReader->GetJetRawPt(jetIndex)) continue; // Cut for jets where all the pT is taken by one track
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
        Aj = (treeReader->GetJetPt(highestIndex) - treeReader->GetJetPt(secondHighestIndex))/(treeReader->GetJetPt(highestIndex) + treeReader->GetJetPt(secondHighestIndex));
        
        // Fill the dijet histogram axes in correct order
        fillerDijet[0] = treeReader->GetJetPt(highestIndex);        // Axis 0: Leading jet pT
        fillerDijet[1] = treeReader->GetJetPhi(highestIndex);       // Axis 1: Leading jet phi
        fillerDijet[2] = treeReader->GetJetEta(highestIndex);       // Axis 2: Leading jet eta
        fillerDijet[3] = treeReader->GetJetPt(secondHighestIndex);  // Axis 3: Subleading jet pT
        fillerDijet[4] = treeReader->GetJetPhi(secondHighestIndex); // Axis 4: Subleading jet phi
        fillerDijet[5] = treeReader->GetJetEta(secondHighestIndex); // Axis 5: Subleading jet eta
        fillerDijet[6] = TMath::Abs(dphi);                          // Axis 6: deltaPhi
        fillerDijet[7] = Aj;                                        // Axis 7: Asymmetry
        fillerDijet[8] = centrality;                                // Axis 8: Centrality
        fHistograms->fhDijet->Fill(fillerDijet);                    // Fill the data point to dijet histogram
        
        // Correlate jets with tracks in dijet events
        int nTracks = treeReader->GetNTracks();
        for(int iTrack = 0; iTrack <nTracks; iTrack++){
          trackEta = treeReader->GetTrackEta(iTrack);
          trackPt = treeReader->GetTrackPt(iTrack);
          trackPhi = treeReader->GetTrackPhi(iTrack);
          trackEt = (treeReader->GetTrackEnergyEcal(iTrack)+treeReader->GetTrackEnergyHcal(iTrack))/TMath::CosH(trackEta);
          
          //  ==== Apply cuts for tracks ====
          if(trackPt < trackMinPtCut) continue;                     // Minimum pT cut
          if(trackPt > jetMaximumPtCut) continue;                   // Maximum pT cut (same as for leading jets)
          if(TMath::Abs(trackEta) > trackEtaCut) continue;          // Eta cut
          if(!treeReader->GetTrackHighPurity(iTrack)) continue;     // High purity cut
          if(treeReader->GetTrackPtError(iTrack)/trackPt > maxTrackPtRelativeError) continue; // Cut for track pT relative error
          if(treeReader->GetTrackVertexDistanceZ(iTrack)/treeReader->GetTrackVertexDistanceZError(iTrack) > maxTrackDistanceToVertex) continue; // Mysterious cut about track proximity to vertex in z-direction
          if(treeReader->GetTrackVertexDistanceXY(iTrack)/treeReader->GetTrackVertexDistanceXYError(iTrack) > maxTrackDistanceToVertex) continue; // Mysterious cut about track proximity to vertex in xy-direction
          if(!(trackPt < calorimeterSignalLimitPt || (trackEt > highPtEtFraction*trackPt))) continue;  // For high pT tracks, require signal also in calorimeters
          if(treeReader->GetTrackChi2(iTrack)/(1.0*treeReader->GetNTrackDegreesOfFreedom(iTrack))/(1.0*treeReader->GetNHitsTrackerLayer(iTrack)) > chi2QualityCut) continue; // Track reconstruction quality cut
          if(treeReader->GetNHitsTrack(iTrack) < minimumTrackHits) continue; // Cut for minimum number of hits per track
          //     ==== Track cuts done ====
          
          // Calculate minimum distance of a track to closest jet. This is needed for track efficiency correction
          trackRMin = 666;   // Initialize the minimum distance to a jet to some very large value
          for(int iJet = 0; iJet < treeReader->GetNJets(); iJet++){
            
            // Require the same jet quality cuts as when searching for dijets
            if(TMath::Abs(treeReader->GetJetEta(iJet)) > jetSearchEtaCut) continue; // Require jet eta to be in the search range for dijets
            if(treeReader->GetJetPt(iJet) < subleadingJetMinPtCut) continue; // Only consider jets that pass the pT cut for subleading jets
            if(minimumMaxTrackPtFraction > treeReader->GetJetMaxTrackPt(iJet)/treeReader->GetJetRawPt(iJet)) continue; // Cut for jets with only very low pT particles
            if(maximumMaxTrackPtFraction < treeReader->GetJetMaxTrackPt(iJet)/treeReader->GetJetRawPt(iJet)) continue; // Cut for jets where all the pT is taken by one track
            // if(TMath::Abs(chargedSum[k]/rawpt[k]) < 0.01) continue; // Jet quality cut from readme file. TODO: Check if should be applied
            
            // Note: The ACos(Cos(jetPhi-trackPhi)) structure transforms deltaPhi to interval [0,Pi]
            trackR = TMath::Power(treeReader->GetJetEta(iJet)-trackEta,2)+TMath::Power(TMath::ACos(TMath::Cos(treeReader->GetJetPhi(iJet)-trackPhi)),2);
            if(trackRMin*trackRMin>trackR) trackRMin=TMath::Power(trackR,0.5);
          } // Loop for calculating Rmin
          
          trackEfficiencyCorrection = fTrackCorrection->getTrkCorr(trackPt, trackEta, trackPhi, hiBin, trackRMin);
          
        } // Loop over tracks
        
      } // Dijet in event
      
    } // Event loop
    
    // Close the input file after the event has been read
    inputFile->Close();
    
  } // File loop
  
  delete treeReader;  // Delete the created tree reader after the analysis is done
}

/*
 * Getter for dijet histograms
 */
DijetHistograms* DijetAnalyzer::GetHistograms() const{
  return fHistograms;
}
