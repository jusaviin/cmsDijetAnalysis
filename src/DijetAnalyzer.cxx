// Class for the main analysis algorithms for the leading-subleading jet analysis

// Root includes
#include <TFile.h>
#include <TMath.h>
#include <TRandom3.h>

// Own includes
#include "DijetAnalyzer.h"

using namespace std;

/*
 * Default constructor
 */
DijetAnalyzer::DijetAnalyzer() :
  fFileNames(0),
  fCard(0),
  fHistograms(0),
  fTrackCorrection(),
  fVzCut(0),
  fJetEtaCut(0),
  fJetSearchEtaCut(0),
  fJetMaximumPtCut(0),
  fLeadingJetMinPtCut(0),
  fSubleadingJetMinPtCut(0),
  fDeltaPhiCut(0),
  fMinimumMaxTrackPtFraction(0),
  fMaximumMaxTrackPtFraction(0),
  fTrackEtaCut(0),
  fTrackMinPtCut(0),
  fMaxTrackPtRelativeError(0),
  fMaxTrackDistanceToVertex(0),
  fCalorimeterSignalLimitPt(0),
  fHighPtEtFraction(0),
  fChi2QualityCut(0),
  fMinimumTrackHits(0)
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
  fHistograms(0),
  fVzCut(0),
  fJetEtaCut(0),
  fJetSearchEtaCut(0),
  fJetMaximumPtCut(0),
  fLeadingJetMinPtCut(0),
  fSubleadingJetMinPtCut(0),
  fDeltaPhiCut(0),
  fMinimumMaxTrackPtFraction(0),
  fMaximumMaxTrackPtFraction(0),
  fTrackEtaCut(0),
  fTrackMinPtCut(0),
  fMaxTrackPtRelativeError(0),
  fMaxTrackDistanceToVertex(0),
  fCalorimeterSignalLimitPt(0),
  fHighPtEtFraction(0),
  fChi2QualityCut(0),
  fMinimumTrackHits(0)
{
  // Custom constructor
  fHistograms = new DijetHistograms(fCard);
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
  fTrackCorrection(in.fTrackCorrection),
  fVzCut(in.fVzCut),
  fJetEtaCut(in.fJetEtaCut),
  fJetSearchEtaCut(in.fJetSearchEtaCut),
  fJetMaximumPtCut(in.fJetMaximumPtCut),
  fLeadingJetMinPtCut(in.fLeadingJetMinPtCut),
  fSubleadingJetMinPtCut(in.fSubleadingJetMinPtCut),
  fDeltaPhiCut(in.fDeltaPhiCut),
  fMinimumMaxTrackPtFraction(in.fMinimumMaxTrackPtFraction),
  fMaximumMaxTrackPtFraction(in.fMaximumMaxTrackPtFraction),
  fTrackEtaCut(in.fTrackEtaCut),
  fTrackMinPtCut(in.fTrackMinPtCut),
  fMaxTrackPtRelativeError(in.fMaxTrackPtRelativeError),
  fMaxTrackDistanceToVertex(in.fMaxTrackDistanceToVertex),
  fCalorimeterSignalLimitPt(in.fCalorimeterSignalLimitPt),
  fHighPtEtFraction(in.fHighPtEtFraction),
  fChi2QualityCut(in.fChi2QualityCut),
  fMinimumTrackHits(in.fMinimumTrackHits)
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
  fVzCut = in.fVzCut;
  fJetEtaCut = in.fJetEtaCut;
  fJetSearchEtaCut = in.fJetSearchEtaCut;
  fJetMaximumPtCut = in.fJetMaximumPtCut;
  fLeadingJetMinPtCut = in.fLeadingJetMinPtCut;
  fSubleadingJetMinPtCut = in.fSubleadingJetMinPtCut;
  fDeltaPhiCut = in.fDeltaPhiCut;
  fMinimumMaxTrackPtFraction = in.fMinimumMaxTrackPtFraction;
  fMaximumMaxTrackPtFraction = in.fMaximumMaxTrackPtFraction;
  fTrackEtaCut = in.fTrackEtaCut;
  fTrackMinPtCut = in.fTrackMinPtCut;
  fMaxTrackPtRelativeError = in.fMaxTrackPtRelativeError;
  fMaxTrackDistanceToVertex = in.fMaxTrackDistanceToVertex;
  fCalorimeterSignalLimitPt = in.fCalorimeterSignalLimitPt;
  fHighPtEtFraction = in.fHighPtEtFraction;
  fChi2QualityCut = in.fChi2QualityCut;
  fMinimumTrackHits = in.fMinimumTrackHits;
  
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
  
  fVzCut = fCard->Get("ZVertexCut");  // Event cut vor the z-posiiton of the primary vertex
  
  //****************************************
  //        Dijet selection cuts
  //****************************************
  
  fJetEtaCut = fCard->Get("JetEtaCut");           // Eta cut around midrapidity
  fJetSearchEtaCut = fCard->Get("SearchEtaCut");  // Eta cut when searching for a dijet
  fJetMaximumPtCut = fCard->Get("MaxPtCut");      // Maximum pT accepted for leading jet (and tracks)
  fLeadingJetMinPtCut = fCard->Get("MinPtCut");   // Minimum pT cut for leading jet
  fSubleadingJetMinPtCut = fCard->Get("SubleadingPtCut"); // Minimum pT cut for subleading jet
  fDeltaPhiCut = fCard->Get("DeltaPhiCut");       // DeltaPhi cut for the dijet system
  fMinimumMaxTrackPtFraction = fCard->Get("MinMaxTrackPtFraction");  // Cut for jets consisting only from soft particles
  fMaximumMaxTrackPtFraction = fCard->Get("MaxMaxTrackPtFraction");  // Cut for jets consisting only from one high pT particle
  
  //****************************************
  //        Track selection cuts
  //****************************************
  
  fTrackEtaCut = fCard->Get("TrackEtaCut");     // Eta cut around midrapidity
  fTrackMinPtCut = fCard->Get("MinTrackPtCut"); // Minimum pT cut
  fMaxTrackPtRelativeError = fCard->Get("MaxTrackPtRelativeError");   // Maximum relative error for pT
  fMaxTrackDistanceToVertex = fCard->Get("VertexMaxDistance");        // Maximum distance to primary vetrex
  fCalorimeterSignalLimitPt = fCard->Get("CalorimeterSignalLimitPt"); // Require signal in calorimeters for track above this pT
  fHighPtEtFraction = fCard->Get("HighPtEtFraction"); // For high pT tracks, minimum required Et as a fraction of track pT
  fChi2QualityCut = fCard->Get("Chi2QualityCut");     // Quality cut for track reconstruction
  fMinimumTrackHits = fCard->Get("MinimumTrackHits"); // Quality cut for track hits
  
  //****************************************
  //            All cuts set!
  //****************************************

  // Input file and forest reader for analysis
  TFile *inputFile;
  TFile *mixedEventFile;
  HighForestReader *treeReader = new HighForestReader(fCard->Get("DataType"));
  ForestReader *mixedEventReader;
  
  // For PbPb, the Forest in mixing file is in different format as for other datasets
  if(fCard->Get("DataType") == ForestReader::kPbPb){
    mixedEventReader = new SkimForestReader(ForestReader::kPbPb);
  } else {
    mixedEventReader = new HighForestReader(fCard->Get("DataType"));
  }
  
  // Event variables
  Int_t nEvents = 0;            // Number of events
  Bool_t dijetFound = false;    // Is there a dijet in the event?
  Bool_t twoJetsFound = false;  // Are there two jets in the event?
  Double_t vz = 0;              // Vertex z-position
  Double_t centrality = 0;      // Event centrality
  Int_t hiBin = 0;              // CMS hiBin (centrality * 2)
  
  // Event mixing information
  Bool_t mixEvents = (fCard->Get("DoEventMixing") == 1);           // Do or do not do event mixing
  Int_t nMixedEventsPerDijet = fCard->Get("NMixedEventsPerDijet"); // Number of events mixed with each dijet event
  TRandom3 *mixedEventRandomizer = new TRandom3();                 // Randomizer for starting point in the mixed event file
  Int_t firstMixingEvent;                                          // Event index from which we start the mixing
  Double_t mixingVzTolerance = fCard->Get("VzTolerance");          // Maximum vz distance between mixed event and dijet event
  Int_t mixedEventIndex;                                           // Index of current event in mixing loop
  Int_t eventsMixed;                                               // Number of events mixed
  Int_t nEventsInMixingFile;                                       // Number of events in the mixing file
  Bool_t allEventsWentThrough;                                     // Have we checked all the events in the mixing file
  std::vector<Double_t> mixedEventVz;                              // Vector for vz in mixing events
  std::vector<Int_t> mixedEventHiBin;                              // Vector for hiBin in mixing events
  
  // Initialize the mixed event randomizer
  mixedEventRandomizer->SetSeed(0);
  
  // Variables for jets
  Double_t dijetAsymmetry = -99;   // Dijet asymmetry
  Double_t leadingJetPt = 0;       // Leading jet pT
  Double_t leadingJetPhi = 0;      // Leading jet phi
  Double_t leadingJetEta = 0;      // Leading jet eta
  Double_t subleadingJetPt = 0;    // Subleading jet pT
  Double_t subleadingJetPhi = 0;   // Subleading jet phi
  Double_t subleadingJetEta = 0;   // Subleading jet eta
  Int_t secondHighestIndex = -1;   // Index of the subleading jet in the event
  Int_t highestIndex = -1;         // Index of the leading jet in the event
  Double_t jetPt = 0;              // pT of the i:th jet in the event
  Double_t dphi = 0;               // deltaPhi for the considered jets
  Double_t leadingJetInfo[3];      // Array for leading jet pT, phi and eta
  Double_t subleadingJetInfo[3];   // Array for subleading jet pT, phi and eta
  
  // File name helper variables
  TString currentFile;
  TString currentMixedEventFile;
  
  // Fillers for THnSparses
  Double_t fillerJet[4];
  Double_t fillerDijet[5];
  
  
  // Amount of debugging messages
  Int_t debugLevel = fCard->Get("DebugLevel");
  
  // Loop over files
  for(Int_t iFile = 0; iFile < (Int_t) fFileNames.size(); iFile++) {
    
    // Find the filename
    currentFile = fFileNames.at(iFile);
    
    // PbPb data has different data file for mixing, other data sets use the regular data files for mixing
    if(fCard->Get("DataType") == ForestReader::kPbPb){
      currentMixedEventFile = "root://cmsxrootd.fnal.gov///store/user/kjung/PbPb_5TeV_MinBiasSkim/Data2015_finalTrkCut_1Mevts.root";
    } else {
      currentMixedEventFile = fFileNames.at(iFile);
    }
    
    // Open the file and check that everything goes fine
    inputFile = TFile::Open(currentFile);
    mixedEventFile = TFile::Open(currentMixedEventFile);
    
    // Check that the file exists
    if(!inputFile){
      cout << "Warning! Could not open the file: " << currentFile.Data() << endl;
      continue;
    }
    
    // Check that the mixing file exists
    if(!mixedEventFile && mixEvents){
      cout << "Warning! Could not open the mixing file: " << currentMixedEventFile.Data() << endl;
      continue;
    }
    
    // Check that the file is not zombie
    if(inputFile->IsZombie()){
      cout << "Warning! The following file is a zombie: " << currentFile.Data() << endl;
      continue;
    }
    
    // Check that the file is not zombie
    if(mixedEventFile->IsZombie() && mixEvents){
      cout << "Warning! The following mixing file is a zombie: " << currentMixedEventFile.Data() << endl;
      continue;
    }
    
    // Debug message, if wanted
    if(debugLevel > 0) cout << "Reading from file: " << currentFile.Data() << endl;
    
    // If file is good, read the forest from the file
    treeReader->ReadForestFromFile(inputFile);  // There seems to be a memory leak here...
    nEvents = treeReader->GetNEvents();
    
    // Read also the forest for event mixing
    if(mixEvents){
      mixedEventReader->ReadForestFromFile(mixedEventFile);  // Read the mixed event forest
      nEventsInMixingFile = mixedEventReader->GetNEvents(); // Read the number of events in the mixing file
      firstMixingEvent = nEventsInMixingFile*mixedEventRandomizer->Rndm();  // Start mixing from random spot in file
      if(firstMixingEvent == nEventsInMixingFile) firstMixingEvent--;  // Move the index to allowed range
      
      // Read vz and hiBin from each event in event mixing file to memory.
      // This way we avoid loading different mixed events in a loop several times
      mixedEventVz.clear();     // Clear the vectors for any possible
      mixedEventHiBin.clear();  // contents they might have
      for(Int_t iMixedEvent = 0; iMixedEvent < nEventsInMixingFile; iMixedEvent++){
        mixedEventReader->GetEvent(iMixedEvent);
        mixedEventVz.push_back(mixedEventReader->GetVz());
        mixedEventHiBin.push_back(mixedEventReader->GetHiBin());
      }
    }

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
      
      // Cut for energy deposition in at least 3 hadronic forward towers. Only applied for PbPb data.
      if(treeReader->GetHfCoincidenceFilterBit() == 0) continue;
      fHistograms->fhEvents->Fill(DijetHistograms::kHfCoincidence);
      
      // Cut for cluster compatibility. Only applied for PbPb data.
      if(treeReader->GetClusterCompatibilityFilterBit() == 0) continue;
      fHistograms->fhEvents->Fill(DijetHistograms::kClusterCompatibility);
      
      // Cut for calirimeter jet quality. Only applied for data.
      if(treeReader->GetCaloJetFilterBit() == 0) continue;
      fHistograms->fhEvents->Fill(DijetHistograms::kCaloJet);
      
      // Cut for vertex z-position
      if(TMath::Abs(vz) > fVzCut) continue;
      fHistograms->fhEvents->Fill(DijetHistograms::kVzCut);
      
      // ===== Event quality cuts applied =====
      
      // Reset the variables used in dijet finding
      twoJetsFound = false;
      dijetFound = false;
      highestIndex = -1;
      secondHighestIndex = -1;
      leadingJetPt = 0;
      subleadingJetPt = 0;
      
      // Search for leading jet and fill histograms for all jets within the eta vut
      for(Int_t jetIndex = 0; jetIndex < treeReader->GetNJets(); jetIndex++) {
        jetPt = treeReader->GetJetPt(jetIndex);
        if(TMath::Abs(treeReader->GetJetEta(jetIndex)) >= fJetSearchEtaCut) continue; // Cut for search eta range
        if(fMinimumMaxTrackPtFraction >= treeReader->GetJetMaxTrackPt(jetIndex)/treeReader->GetJetRawPt(jetIndex)) continue; // Cut for jets with only very low pT particles
        if(fMaximumMaxTrackPtFraction <= treeReader->GetJetMaxTrackPt(jetIndex)/treeReader->GetJetRawPt(jetIndex)) continue; // Cut for jets where all the pT is taken by one track
        
        // Fill the histogram for all jets within eta range
        if(TMath::Abs(treeReader->GetJetEta(jetIndex)) < fJetEtaCut){
          
          // Fill the axes in correct order
          fillerJet[0] = jetPt;                             // Axis 0 = any jet pT
          fillerJet[1] = treeReader->GetJetPhi(jetIndex);   // Axis 1 = any jet phi
          fillerJet[2] = treeReader->GetJetEta(jetIndex);   // Axis 2 = any jet eta
          fillerJet[3] = centrality;                        // Axis 3 = centrality
          fHistograms->fhAnyJet->Fill(fillerJet);           // Fill the data point to histogram
          
        }
        
        if(jetPt <= fLeadingJetMinPtCut) continue; // Minimum leading jet pT cut
        if(jetPt > leadingJetPt){
          leadingJetPt = jetPt;
          highestIndex = jetIndex;
        }
      } // End of search for leading jet loop
      
      // Search for subleading jet
      for(Int_t jetIndex = 0 ; jetIndex < treeReader->GetNJets(); jetIndex++){
        jetPt = treeReader->GetJetPt(jetIndex);
        if(jetIndex == highestIndex) continue; // Do not consider leading particle
        if(TMath::Abs(treeReader->GetJetEta(jetIndex)) >= fJetSearchEtaCut) continue; // Cut for search eta range
        if(fMinimumMaxTrackPtFraction >= treeReader->GetJetMaxTrackPt(jetIndex)/treeReader->GetJetRawPt(jetIndex)) continue; // Cut for jets with only very low pT particles
        if(fMaximumMaxTrackPtFraction <= treeReader->GetJetMaxTrackPt(jetIndex)/treeReader->GetJetRawPt(jetIndex)) continue; // Cut for jets where all the pT is taken by one track
        if(jetPt <= fSubleadingJetMinPtCut) continue; // Minimum subleading jet pT cut
        if(jetPt > subleadingJetPt){
          subleadingJetPt = jetPt;
          secondHighestIndex = jetIndex;
        }
      }  //End of subleading jet search
      
      // Check if at least two jets were found
      if(highestIndex > -1 && secondHighestIndex > -1) twoJetsFound = true;
      
      // Only apply the dijet cuts for events with at least two jets
      if(twoJetsFound){
        
        // Read the eta and phi values for leading and subleading jets
        leadingJetPhi = treeReader->GetJetPhi(highestIndex);
        leadingJetEta = treeReader->GetJetEta(highestIndex);
        subleadingJetPhi = treeReader->GetJetPhi(secondHighestIndex);
        subleadingJetEta = treeReader->GetJetEta(secondHighestIndex);
        
        dijetFound = true;
        dphi =  leadingJetPhi - subleadingJetPhi;
        if(dphi < 0) dphi = -dphi;
        if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
        
        
        if((leadingJetPt >= fJetMaximumPtCut) ||              // Maximum leading jet pT cut
           (TMath::Abs(leadingJetEta) >= fJetEtaCut) ||       // Leading jet eta cut
           (TMath::Abs(subleadingJetEta) >= fJetEtaCut)||  // Subleading jet eta cut
           (TMath::Abs(dphi) <= fDeltaPhiCut)){                                     // DeltaPhi cut
          dijetFound = false;
        }
      } // End of dijet cuts
      
      // If a dijet is found, fill some information to fHistograms
      if(dijetFound){
        fHistograms->fhEvents->Fill(DijetHistograms::kDijet);
        fHistograms->fhCentralityDijet->Fill(centrality);
        
        // Calculate the asymmetry
        dijetAsymmetry = (leadingJetPt - subleadingJetPt)/(leadingJetPt + subleadingJetPt);
        
        // Fill the leading jet histogram
        fillerDijet[0] = leadingJetPt;                   // Axis 0: Leading jet pT
        fillerDijet[1] = leadingJetPhi;                  // Axis 1: Leading jet phi
        fillerDijet[2] = leadingJetEta;                  // Axis 2: Leading jet eta
        fillerDijet[3] = dijetAsymmetry;                 // Axis 3: Asymmetry
        fillerDijet[4] = centrality;                     // Axis 4: Centrality
        fHistograms->fhLeadingJet->Fill(fillerDijet);    // Fill the data point to leading jet histogram
        
        // Fill the subleading jet histogram
        fillerDijet[0] = subleadingJetPt;                // Axis 0: Subleading jet pT
        fillerDijet[1] = subleadingJetPhi;               // Axis 1: Subleading jet phi
        fillerDijet[2] = subleadingJetEta;               // Axis 2: Subleading jet eta
        fillerDijet[3] = dijetAsymmetry;                 // Axis 3: Asymmetry
        fillerDijet[4] = centrality;                     // Axis 4: Centrality
        fHistograms->fhSubleadingJet->Fill(fillerDijet); // Fill the data point to subleading jet histogram

        // Fill the dijet histogram
        fillerDijet[0] = leadingJetPt;                   // Axis 0: Leading jet pT
        fillerDijet[1] = subleadingJetPt;                // Axis 1: Subleading jet pT
        fillerDijet[2] = TMath::Abs(dphi);               // Axis 2: deltaPhi
        fillerDijet[3] = dijetAsymmetry;                 // Axis 3: Asymmetry
        fillerDijet[4] = centrality;                     // Axis 4: Centrality
        fHistograms->fhDijet->Fill(fillerDijet);         // Fill the data point to dijet histogram
        
        // Fill the arrays with leading and subleading jet information for correlation with tracks
        leadingJetInfo[0] = leadingJetPt;
        leadingJetInfo[1] = leadingJetPhi;
        leadingJetInfo[2] = leadingJetEta;
        subleadingJetInfo[0] = subleadingJetPt;
        subleadingJetInfo[1] = subleadingJetPhi;
        subleadingJetInfo[2] = subleadingJetEta;
        
        // Correlate jets with tracks in dijet events
        CorrelateTracksAndJets(treeReader,leadingJetInfo,subleadingJetInfo,DijetHistograms::kSameEvent);
        
        // Do event mixing
        if(mixEvents){
          
          // Start mixing from the first event index
          mixedEventIndex = firstMixingEvent;
          eventsMixed = 0;
          allEventsWentThrough = false;
          
          // Continue mixing until we have reached required number of event or no event candidates remain in the file
          while (eventsMixed < nMixedEventsPerDijet && !allEventsWentThrough) {
            
            // Increment the counter for event index to be mixed with the current event
            mixedEventIndex++;
            
            // If we are out of bounds from the event in data file, reset the counter
            if(mixedEventIndex == nEventsInMixingFile) {
              mixedEventIndex = -1;
              continue;
            }
            
            // If we come back to the first event index, we have gone through all the events without finding 20 similar events form the file
            if(mixedEventIndex == firstMixingEvent) allEventsWentThrough = true;
            
            // Do not mix with the same event
            if((mixedEventIndex == iEvent) && (currentFile == currentMixedEventFile)) continue;
            
            // Match vz and hiBin between the dijet event and mixed event
            if(TMath::Abs(mixedEventVz.at(mixedEventIndex) - vz) > mixingVzTolerance) continue;
            if(mixedEventHiBin.at(mixedEventIndex) != hiBin) continue;
            
            // If match vz and hiBin, then load the event from the mixed event tree
            mixedEventReader->GetEvent(mixedEventIndex);
            
            // Do the correlations with the dijet from current event and track from mixing event
            CorrelateTracksAndJets(mixedEventReader,leadingJetInfo,subleadingJetInfo,DijetHistograms::kMixedEvent);
            eventsMixed++;
            
          } // While loop for finding events to mix
          
          // Print out a message if we could not find 20 events to mix with the current event
          if(allEventsWentThrough && (debugLevel > 0)){
            cout << "Could only find " << eventsMixed << " events to mix with event " << iEvent << " in file " << currentFile.Data() << endl;
          }
          
          // For the next event, start mixing the events from where we were left with in the previous event
          firstMixingEvent = mixedEventIndex;
          
        } // Event mixing
      } // Dijet in event
      
    } // Event loop
    
    // Close the input files after the event has been read
    inputFile->Close();
    mixedEventFile->Close();
    
  } // File loop
  
  // When we are done, delete ForestReaders
  delete treeReader;
  delete mixedEventReader;
}

/*
 * Method for all jet-track correlations
 *
 *  ForestReader *jetReader = Forest from which the track information is read
 *  Double_t leadingJetInfo[3] = Array containing leading jet pT, phi and eta
 *  Double_t subleadingJetInfo[3] = Array containing subleading jet pT, phi and eta
 *  Int_t correlationType = DijetHistograms::kSameEvent for same event correlations, DijetHistograms::kMixedEvent for mixed event correlations
 */
void DijetAnalyzer::CorrelateTracksAndJets(ForestReader *treeReader, Double_t leadingJetInfo[3], Double_t subleadingJetInfo[3], Int_t correlationType){
  
  // Define a filler for THnSparses
  Double_t fillerJetTrack[6];
  Double_t fillerTrack[5];
  
  // Event information
  Int_t hiBin = treeReader->GetHiBin();
  Int_t centrality = treeReader->GetCentrality();
  
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
  
  // Read the leading jet and subleading jet information from arrays
  Double_t leadingJetPt = leadingJetInfo[0];
  Double_t leadingJetPhi = leadingJetInfo[1];
  Double_t leadingJetEta = leadingJetInfo[2];
  Double_t subleadingJetPt = subleadingJetInfo[0];
  Double_t subleadingJetPhi = subleadingJetInfo[1];
  Double_t subleadingJetEta = subleadingJetInfo[2];
  
  // Calculate the dijet asymmetry
  Double_t dijetAsymmetry = (leadingJetPt - subleadingJetPt)/(leadingJetPt + subleadingJetPt);
  
  // Loop over all track in the event
  Int_t nTracks = treeReader->GetNTracks();
  for(Int_t iTrack = 0; iTrack <nTracks; iTrack++){
    trackPt = treeReader->GetTrackPt(iTrack);
    trackPhi = treeReader->GetTrackPhi(iTrack);
    trackEta = treeReader->GetTrackEta(iTrack);
    trackEt = (treeReader->GetTrackEnergyEcal(iTrack)+treeReader->GetTrackEnergyHcal(iTrack))/TMath::CosH(trackEta);
    
    //  ==== Apply cuts for tracks and collect information on how much track are cut in each step ====
    
    // Only fill the track cut histograms for same event data
    if(correlationType == DijetHistograms::kSameEvent) fHistograms->fhTrackCuts->Fill(DijetHistograms::kAllTracks);
    
    // Cut for track pT
    if(trackPt <= fTrackMinPtCut) continue;                     // Minimum pT cut
    if(trackPt >= fJetMaximumPtCut) continue;                   // Maximum pT cut (same as for leading jets)
    if(correlationType == DijetHistograms::kSameEvent) fHistograms->fhTrackCuts->Fill(DijetHistograms::kPtCuts);
    
    // Cut for track eta
    if(TMath::Abs(trackEta) >= fTrackEtaCut) continue;          // Eta cut
    if(correlationType == DijetHistograms::kSameEvent) fHistograms->fhTrackCuts->Fill(DijetHistograms::kEtaCut);
    
    // Cut for high purity
    if(!treeReader->GetTrackHighPurity(iTrack)) continue;     // High purity cut
    if(correlationType == DijetHistograms::kSameEvent) fHistograms->fhTrackCuts->Fill(DijetHistograms::kHighPurity);
    
    // Cut for relative error for track pT
    if(treeReader->GetTrackPtError(iTrack)/trackPt >= fMaxTrackPtRelativeError) continue; // Cut for track pT relative error
    if(correlationType == DijetHistograms::kSameEvent) fHistograms->fhTrackCuts->Fill(DijetHistograms::kPtError);
    
    // Cut for track distance from primary vertex
    if(treeReader->GetTrackVertexDistanceZ(iTrack)/treeReader->GetTrackVertexDistanceZError(iTrack) >= fMaxTrackDistanceToVertex) continue; // Mysterious cut about track proximity to vertex in z-direction
    if(treeReader->GetTrackVertexDistanceXY(iTrack)/treeReader->GetTrackVertexDistanceXYError(iTrack) >= fMaxTrackDistanceToVertex) continue; // Mysterious cut about track proximity to vertex in xy-direction
    if(correlationType == DijetHistograms::kSameEvent) fHistograms->fhTrackCuts->Fill(DijetHistograms::kVertexDistance);
    
    // Cut for energy deposition in calorimeters for high pT tracks
    if(!(trackPt < fCalorimeterSignalLimitPt || (trackEt >= fHighPtEtFraction*trackPt))) continue;  // For high pT tracks, require signal also in calorimeters
    if(correlationType == DijetHistograms::kSameEvent) fHistograms->fhTrackCuts->Fill(DijetHistograms::kCaloSignal);
    
    // Cuts for track reconstruction quality
    if(treeReader->GetTrackChi2(iTrack)/(1.0*treeReader->GetNTrackDegreesOfFreedom(iTrack))/(1.0*treeReader->GetNHitsTrackerLayer(iTrack)) >= fChi2QualityCut) continue; // Track reconstruction quality cut
    if(treeReader->GetNHitsTrack(iTrack) < fMinimumTrackHits) continue; // Cut for minimum number of hits per track
    if(correlationType == DijetHistograms::kSameEvent) fHistograms->fhTrackCuts->Fill(DijetHistograms::kReconstructionQuality);
    
    //     ==== Track cuts done ====
    
    // Calculate minimum distance of a track to closest jet. This is needed for track efficiency correction
    trackRMin = 666;   // Initialize the minimum distance to a jet to some very large value
    for(Int_t iJet = 0; iJet < treeReader->GetNJets(); iJet++){
      
      // Require the same jet quality cuts as when searching for dijets
      if(TMath::Abs(treeReader->GetJetEta(iJet)) >= fJetSearchEtaCut) continue; // Require jet eta to be in the search range for dijets
      if(treeReader->GetJetPt(iJet) <= fSubleadingJetMinPtCut) continue; // Only consider jets that pass the pT cut for subleading jets
      if(fMinimumMaxTrackPtFraction >= treeReader->GetJetMaxTrackPt(iJet)/treeReader->GetJetRawPt(iJet)) continue; // Cut for jets with only very low pT particles
      if(fMaximumMaxTrackPtFraction <= treeReader->GetJetMaxTrackPt(iJet)/treeReader->GetJetRawPt(iJet)) continue; // Cut for jets where all the pT is taken by one track
      // if(TMath::Abs(chargedSum[k]/rawpt[k]) < 0.01) continue; // Jet quality cut from readme file. TODO: Check if should be applied
      
      // Note: The ACos(Cos(jetPhi-trackPhi)) structure transforms deltaPhi to interval [0,Pi]
      trackR = TMath::Power(treeReader->GetJetEta(iJet)-trackEta,2)+TMath::Power(TMath::ACos(TMath::Cos(treeReader->GetJetPhi(iJet)-trackPhi)),2);
      if(trackRMin*trackRMin>trackR) trackRMin=TMath::Power(trackR,0.5);
      
    } // Loop for calculating Rmin
    
    // Get the track efficiency corrections
    trackEfficiencyCorrection = fTrackCorrection->getTrkCorr(trackPt, trackEta, trackPhi, hiBin, trackRMin);
    
    // Calculate deltaEta and deltaPhi between track and leading and subleading jets
    deltaEtaTrackLeadingJet = leadingJetEta - trackEta;
    deltaPhiTrackLeadingJet = leadingJetPhi - trackPhi;
    deltaEtaTrackSubleadingJet = subleadingJetEta - trackEta;
    deltaPhiTrackSubleadingJet = subleadingJetPhi - trackPhi;
    
    // Tranform deltaPhis to interval [-pi/2,3pi/2]
    while(deltaPhiTrackLeadingJet > (1.5*TMath::Pi())){deltaPhiTrackLeadingJet += -2*TMath::Pi();}
    while(deltaPhiTrackSubleadingJet > (1.5*TMath::Pi())){deltaPhiTrackSubleadingJet += -2*TMath::Pi();}
    while(deltaPhiTrackLeadingJet < (-0.5*TMath::Pi())){deltaPhiTrackLeadingJet += 2*TMath::Pi();}
    while(deltaPhiTrackSubleadingJet < (-0.5*TMath::Pi())){deltaPhiTrackSubleadingJet += 2*TMath::Pi();}
    
    // Fill track histograms
    fillerTrack[0] = trackPt;                    // Axis 0: Track pT
    fillerTrack[1] = trackPhi;                   // Axis 1: Track phi
    fillerTrack[2] = trackEta;                   // Axis 2: Track eta
    fillerTrack[3] = centrality;                 // Axis 3: Centrality
    fillerTrack[4] = correlationType;            // Axis 4: Correlation type (same or mixed event)
    fHistograms->fhTrack->Fill(fillerTrack,trackEfficiencyCorrection);  // Fill the track histogram
    fHistograms->fhTrackUncorrected->Fill(fillerTrack);                 // Fill the uncorrected track histogram
    
    // Fill the track-leading jet correlation histograms
    fillerJetTrack[0] = trackPt;                    // Axis 0: Track pT
    fillerJetTrack[1] = deltaPhiTrackLeadingJet;    // Axis 1: DeltaPhi between track and leading jet
    fillerJetTrack[2] = deltaEtaTrackLeadingJet;    // Axis 2: DeltaEta between track and leading jet
    fillerJetTrack[3] = dijetAsymmetry;             // Axis 3: Dijet asymmetry
    fillerJetTrack[4] = centrality;                 // Axis 4: Centrality
    fillerJetTrack[5] = correlationType;            // Axis 5: Correlation type (same or mixed event)
    fHistograms->fhTrackLeadingJet->Fill(fillerJetTrack,trackEfficiencyCorrection); // Fill the track-leading jet correlation histogram
    fHistograms->fhTrackLeadingJetUncorrected->Fill(fillerJetTrack);                // Fill the uncorrected track-leading jet correlation histogram
    fHistograms->fhTrackLeadingJetPtWeighted->Fill(fillerJetTrack,trackEfficiencyCorrection*trackPt); // Fill the pT weighted track-leading jet correlation histogram
    
    // Fill the track-subleading jet correlation histograms
    fillerJetTrack[0] = trackPt;                    // Axis 0: Track pT
    fillerJetTrack[1] = deltaPhiTrackSubleadingJet; // Axis 1: DeltaPhi between track and subleading jet
    fillerJetTrack[2] = deltaEtaTrackSubleadingJet; // Axis 2: DeltaEta between track and subleading jet
    fillerJetTrack[3] = dijetAsymmetry;             // Axis 3: Dijet asymmetry
    fillerJetTrack[4] = centrality;                 // Axis 4: Centrality
    fillerJetTrack[5] = correlationType;            // Axis 5: Correlation type (same or mixed event)
    fHistograms->fhTrackSubleadingJet->Fill(fillerJetTrack,trackEfficiencyCorrection); // Fill the track-subleading jet correlation histogram
    fHistograms->fhTrackSubleadingJetUncorrected->Fill(fillerJetTrack);                // Fill the uncorrected track-subleading jet correlation histogram
    fHistograms->fhTrackSubleadingJetPtWeighted->Fill(fillerJetTrack,trackEfficiencyCorrection*trackPt); // Fill the pT weighted track-subleading jet correlation histogram
    
  } // Loop over tracks
}

/*
 * Getter for dijet histograms
 */
DijetHistograms* DijetAnalyzer::GetHistograms() const{
  return fHistograms;
}
