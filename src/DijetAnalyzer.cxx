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
  fJffCorrection(),
  fVzWeightFunction(0),
  fCentralityWeightFunction(0),
  fDataType(-1),
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
  fMinimumTrackHits(0),
  fSubeventCut(0),
  fMcCorrelationType(0),
  fFilledHistograms(kFillAll)
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
  fMinimumTrackHits(0),
  fSubeventCut(0),
  fMcCorrelationType(0),
  fFilledHistograms(kFillAll)
{
  // Custom constructor
  fHistograms = new DijetHistograms(fCard);
  fHistograms->CreateHistograms();
  
  // Find the correct folder for track correction tables based on data type
  fDataType = fCard->Get("DataType");
  bool ppData = true;
  if(fDataType == ForestReader::kPp || fDataType == ForestReader::kPpMC || fDataType == ForestReader::kLocalTest){
    
    // Track correction
    fTrackCorrection = new TrkCorr("trackCorrectionTables/TrkCorr_July22_Iterative_pp_eta2p4/");
    
    // Common vz weight function used by UIC group for pp MC
    fVzWeightFunction = new TF1("fvz","gaus",-15,15);
    fVzWeightFunction->SetParameter(0,1.10477);
    fVzWeightFunction->SetParameter(1,2.52738);
    fVzWeightFunction->SetParameter(2,1.30296e1);
    
    fCentralityWeightFunction = NULL;
    
  } else if (fDataType == ForestReader::kPbPb || fDataType == ForestReader::kPbPbMC){
    
    // Track correction
    fTrackCorrection = new TrkCorr("trackCorrectionTables/TrkCorr_Jun7_Iterative_PbPb_etaLT2p4/");
    
    // Flag for PbPb data
    ppData = false;
    
    // Common vz weight function used by UIC group for PbPb MC
    fVzWeightFunction = new TF1("fvz","pol6",-15,15);
    fVzWeightFunction->SetParameters(1.18472,-0.132675,0.00857998,-0.000326085,-1.48786e-06,4.68665e-07,-7.32942e-09);
    
    // Comment centrality weight function used by UIC group for PbPb MC
    fCentralityWeightFunction = new TF1("fcent1","[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4+[7]*exp([5]+[6]*x)",0,180);
    fCentralityWeightFunction->SetParameters(4.40810, -7.75301e-02, 4.91953e-04, -1.34961e-06, 1.44407e-09, -160, 1, 3.68078e-15);
    
  } else {
    fTrackCorrection = new TrkCorr(""); // Bad data type, no corrections initialized
    fVzWeightFunction = NULL;
    fCentralityWeightFunction = NULL;
  }
  
  // Initialize the class for JFF correction
  fJffCorrection = new JffCorrection(ppData);
  
  // Read which histograms to fill this run
  fFilledHistograms = fCard->Get("FilledHistograms");
  if(fFilledHistograms < 0) {
    cout << "WARNING: Negative value given for FilledHistograms in ConfigurationCard!" << endl;
    cout << "Filling all histograms" << endl;
    fFilledHistograms = kFillAll;
  }
  if(fFilledHistograms >= knFillModes){
    cout << "WARNING: Integer larger than " << knFillModes-1 << " given for FilledHistograms in ConfigurationCard!" << endl;
    cout << "Please give smaller positive integer. Will fill all histograms now" << endl;
    fFilledHistograms = kFillAll;
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
  fJffCorrection(in.fJffCorrection),
  fVzWeightFunction(in.fVzWeightFunction),
  fCentralityWeightFunction(in.fCentralityWeightFunction),
  fDataType(in.fDataType),
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
  fMinimumTrackHits(in.fMinimumTrackHits),
  fSubeventCut(in.fSubeventCut),
  fMcCorrelationType(in.fMcCorrelationType)
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
  fJffCorrection = in.fJffCorrection;
  fVzWeightFunction = in.fVzWeightFunction;
  fCentralityWeightFunction = in.fCentralityWeightFunction;
  fDataType = in.fDataType;
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
  fSubeventCut = in.fSubeventCut;
  fMcCorrelationType = in.fMcCorrelationType;
  
  return *this;
}

/*
 * Destructor
 */
DijetAnalyzer::~DijetAnalyzer(){
  // destructor
  delete fHistograms;
  delete fTrackCorrection;
  delete fJffCorrection;
  if(fVzWeightFunction) delete fVzWeightFunction;
  if(fCentralityWeightFunction) delete fCentralityWeightFunction;
}

/*
 * Get the number of particle flow candidates within a range of 0.4 of the given eta-phi angle
 */
Int_t DijetAnalyzer::GetNParticleFlowCandidatesInJet(ForestReader *treeReader, Double_t jetPhi, Double_t jetEta){

  // No correction for generator level jets
  if(treeReader->GetNParticleFlowCandidates() < 0) return -1;
  
  // Start counting from zero
  Int_t nParticleFlowCandidatesInThisJet = 0;
  Double_t distanceToThisJet = 0;
  
  // Loop over all particle flow candidates
  for(Int_t iParticleFlowCandidate = 0; iParticleFlowCandidate < treeReader->GetNParticleFlowCandidates(); iParticleFlowCandidate++){
    if(treeReader->GetParticleFlowCandidateId(iParticleFlowCandidate) != 1) continue; // Require ID 1 for candidates
    if(treeReader->GetParticleFlowCandidatePt(iParticleFlowCandidate) < 2) continue; // Require minimum pT of 2 GeV for candidates
    distanceToThisJet = TMath::Power(TMath::Power(jetPhi-treeReader->GetParticleFlowCandidatePhi(iParticleFlowCandidate),2)+TMath::Power(jetEta-treeReader->GetParticleFlowCandidateEta(iParticleFlowCandidate),2),0.5);
    if(distanceToThisJet > 0.4) continue;  // Require the particle to be inside the jet cone of R = 0.4
    nParticleFlowCandidatesInThisJet++;
  }
  
  // Return all we find
  return nParticleFlowCandidatesInThisJet;
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
  
  fSubeventCut = fCard->Get("SubeventCut");     // Required subevent type
  
  //****************************************
  //    Correlation type for Monte Carlo
  //****************************************
  fMcCorrelationType = fCard->Get("McCorrelationType");         // Correlation type for Monte Carlo
  
  //****************************************
  //            All cuts set!
  //****************************************

  // Input files and forest readers for analysis
  TFile *inputFile;
  TFile *copyInputFile; // If we read forest for tracks and jets with different readers, we need to give different file pointers to them
  TFile *mixedEventFile;
  ForestReader *jetReader;
  ForestReader *trackReader;
  ForestReader *mixedEventReader;
  
  // Event variables
  Int_t nEvents = 0;            // Number of events
  Bool_t dijetFound = false;    // Is there a dijet in the event?
  Bool_t twoJetsFound = false;  // Are there two jets in the event?
  Double_t vz = 0;              // Vertex z-position
  Double_t centrality = 0;      // Event centrality
  Int_t hiBin = 0;              // CMS hiBin (centrality * 2)
  Int_t pthat = 0;              // pT hat for MC events
  
  // Combining bools to make the code more readable
  Bool_t fillMainLoopHistograms = (fFilledHistograms == kFillAll || fFilledHistograms == kFillAllButJetTrack); // Select to fill the histograms in the main loop
  Bool_t fillEventInformationHistograms = (fFilledHistograms == kFillAll || fFilledHistograms == kFillJetTrack || fFilledHistograms == kFillOnlyEventInformation); // Need number of jets for normalizing jet-track correlation histograms
  Bool_t useDifferentReaderFotJetsAndTracks = (fMcCorrelationType == kRecoGen || fMcCorrelationType == kGenReco); // Use different forest reader for jets and tracks
  
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
  Bool_t checkForDuplicates;                                       // Check if we have already mixed with the corresponding event
  Bool_t skipEvent;                                                // Skip the current mixing event
  std::vector<Double_t> mixedEventVz;                              // Vector for vz in mixing events
  std::vector<Int_t> mixedEventHiBin;                              // Vector for hiBin in mixing events
  std::vector<Int_t> mixedEventIndices;                            // Indices for events that have already been mixed
  Double_t additionalTolerance;                                    // Increased vz tolerance, if not enough events found from the mixing file
  
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
  Double_t swapJetPt = 0;          // Swapping helper variable
  Double_t swapJetPhi = 0;         // Swapping helper variable
  Double_t swapJetEta = 0;         // Swapping helper variable
  Int_t secondHighestIndex = -1;   // Index of the subleading jet in the event
  Int_t highestIndex = -1;         // Index of the leading jet in the event
  Int_t swapIndex = -1;            // Swapping helper variable
  Double_t jetPt = 0;              // pT of the i:th jet in the event
  Double_t jetPhi = 0;             // phi of the i:th jet in the event
  Double_t jetEta = 0;             // eta of the i:th jet in the event
  Double_t dphi = 0;               // deltaPhi for the considered jets
  Double_t leadingJetInfo[3];      // Array for leading jet pT, phi and eta
  Double_t subleadingJetInfo[3];   // Array for subleading jet pT, phi and eta
  
  // Variables for particle flow candidates
  Int_t nParticleFlowCandidatesInThisJet;  // Number of particle flow candidates in the current jet
  
  // File name helper variables
  TString currentFile;
  TString currentMixedEventFile;
  
  // Fillers for THnSparses
  Double_t fillerJet[4];
  Double_t fillerDijet[5];
  
  // Select the reader for jets
  if(fMcCorrelationType == kGenReco || fMcCorrelationType == kGenGen){
    jetReader = new GeneratorLevelForestReader(fDataType);
  } else {
    jetReader = new HighForestReader(fDataType);
  }
  
  // Select the reader for tracks
  if(fMcCorrelationType == kRecoGen){
    trackReader = new GeneratorLevelForestReader(fDataType);
  } else if (fMcCorrelationType == kGenReco){
    trackReader = new HighForestReader(fDataType);
  } else {
    trackReader = jetReader;
  }
  
  // If mixing events, create ForestReader for that. For PbPb, the Forest in mixing file is in different format as for other datasets
  if(mixEvents){
    if(fDataType == ForestReader::kPbPb){
      mixedEventReader = new SkimForestReader(ForestReader::kPbPb);
    } else if (fMcCorrelationType == kRecoGen || fMcCorrelationType ==kGenGen) { // Mixed event reader for generator tracks
      mixedEventReader = new GeneratorLevelForestReader(fDataType);
    } else {
      mixedEventReader = new HighForestReader(fDataType);
    }
  }
  
  // Amount of debugging messages
  Int_t debugLevel = fCard->Get("DebugLevel");
  
  // Loop over files
  for(Int_t iFile = 0; iFile < (Int_t) fFileNames.size(); iFile++) {
    
    // Find the filename
    currentFile = fFileNames.at(iFile);
    
    // PbPb data has different data file for mixing, other data sets use the regular data files for mixing
    if(fDataType == ForestReader::kPbPb){
      currentMixedEventFile = "root://cmsxrootd.fnal.gov///store/user/kjung/PbPb_5TeV_MinBiasSkim/Data2015_finalTrkCut_1Mevts.root";
    } else {
      currentMixedEventFile = fFileNames.at(iFile);
    }
    
    // Open the file and check that everything goes fine
    inputFile = TFile::Open(currentFile);
    if(useDifferentReaderFotJetsAndTracks) copyInputFile = TFile::Open(currentFile);
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
    jetReader->ReadForestFromFile(inputFile);  // There might be a memory leak in handling the forest...
    if(useDifferentReaderFotJetsAndTracks) trackReader->ReadForestFromFile(copyInputFile); // If we mix reco and gen, the reader for jets and tracks is different
    nEvents = jetReader->GetNEvents();
    
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
      jetReader->GetEvent(iEvent);
      
      // If track reader is not the same as jet reader, read the event to memory in trackReader
      if(useDifferentReaderFotJetsAndTracks) trackReader->GetEvent(iEvent);

      // Get vz and centrality information from all events
      vz = jetReader->GetVz();
      centrality = jetReader->GetCentrality();
      hiBin = jetReader->GetHiBin();
      pthat = jetReader->GetPtHat();
      if(fillEventInformationHistograms){
        fHistograms->fhVertexZ->Fill(vz);                   // z vertex distribution from all events
        fHistograms->fhEvents->Fill(DijetHistograms::kAll); // All the events looped over
        fHistograms->fhCentrality->Fill(centrality);        // Centrality filled from all events
        fHistograms->fhPtHat->Fill(pthat);                  // pT hat histogram
      }
      //  ===== Apply all the event quality cuts =====
      
      // Cut for primary vertex. Only applied for data.
      if(jetReader->GetPrimaryVertexFilterBit() == 0) continue;
      if(fillEventInformationHistograms) fHistograms->fhEvents->Fill(DijetHistograms::kPrimaryVertex);
      
      // Cut for HB/HE noise. Only applied for data.
      if(jetReader->GetHBHENoiseFilterBit() == 0) continue;
      if(fillEventInformationHistograms) fHistograms->fhEvents->Fill(DijetHistograms::kHBHENoise);
      
      // Cut for collision event selection. Only applied for PbPb data.
      if(jetReader->GetCollisionEventSelectionFilterBit() == 0) continue;
      if(fillEventInformationHistograms) fHistograms->fhEvents->Fill(DijetHistograms::kCollisionEventSelection);
      
      // Cut for beam scraping. Only applied for pp data.
      if(jetReader->GetBeamScrapingFilterBit() == 0) continue;
      if(fillEventInformationHistograms) fHistograms->fhEvents->Fill(DijetHistograms::kBeamScraping);
      
      // Cut for energy deposition in at least 3 hadronic forward towers. Only applied for PbPb data.
      if(jetReader->GetHfCoincidenceFilterBit() == 0) continue;
      if(fillEventInformationHistograms) fHistograms->fhEvents->Fill(DijetHistograms::kHfCoincidence);
      
      // Cut for cluster compatibility. Only applied for PbPb data.
      if(jetReader->GetClusterCompatibilityFilterBit() == 0) continue;
      if(fillEventInformationHistograms) fHistograms->fhEvents->Fill(DijetHistograms::kClusterCompatibility);
      
      // Cut for calirimeter jet quality. Only applied for data.
      if(jetReader->GetCaloJetFilterBit() == 0) continue;
      if(fillEventInformationHistograms) fHistograms->fhEvents->Fill(DijetHistograms::kCaloJet);
      
      // Cut for vertex z-position
      if(TMath::Abs(vz) > fVzCut) continue;
      if(fillEventInformationHistograms) fHistograms->fhEvents->Fill(DijetHistograms::kVzCut);
      
      // ===== Event quality cuts applied =====
      
      // Reset the variables used in dijet finding
      twoJetsFound = false;
      dijetFound = false;
      highestIndex = -1;
      secondHighestIndex = -1;
      leadingJetPt = 0;
      subleadingJetPt = 0;
      
      // Search for leading jet and fill histograms for all jets within the eta vut
      for(Int_t jetIndex = 0; jetIndex < jetReader->GetNJets(); jetIndex++) {
        jetPt = jetReader->GetJetPt(jetIndex);
        jetPhi = jetReader->GetJetPhi(jetIndex);
        jetEta = jetReader->GetJetEta(jetIndex);
        if(TMath::Abs(jetEta) >= fJetSearchEtaCut) continue; // Cut for search eta range
        if(fMinimumMaxTrackPtFraction >= jetReader->GetJetMaxTrackPt(jetIndex)/jetReader->GetJetRawPt(jetIndex)) continue; // Cut for jets with only very low pT particles
        if(fMaximumMaxTrackPtFraction <= jetReader->GetJetMaxTrackPt(jetIndex)/jetReader->GetJetRawPt(jetIndex)) continue; // Cut for jets where all the pT is taken by one track
        
        // Only fill the any jet histogram if selected
        if(fillMainLoopHistograms){
          // Fill the histogram for all jets within eta range
          if(TMath::Abs(jetEta) < fJetEtaCut){
          
            // Fill the axes in correct order
            nParticleFlowCandidatesInThisJet = GetNParticleFlowCandidatesInJet(jetReader,jetPhi,jetEta);       // Apply JFF correction for jet pT
            fillerJet[0] = fJffCorrection->GetCorrection(nParticleFlowCandidatesInThisJet,hiBin,jetPt,jetEta);  // Axis 0 = any jet pT
            fillerJet[1] = jetPhi;                  // Axis 1 = any jet phi
            fillerJet[2] = jetEta;                  // Axis 2 = any jet eta
            fillerJet[3] = centrality;              // Axis 3 = centrality
            fHistograms->fhAnyJet->Fill(fillerJet); // Fill the data point to histogram
        
          } // Eta cut
        } // Check if we want to fill any jet histograms
        
        if(jetPt > leadingJetPt){
          leadingJetPt = jetPt;
          highestIndex = jetIndex;
        }
      } // End of search for leading jet loop
      
      // Search for subleading jet
      for(Int_t jetIndex = 0 ; jetIndex < jetReader->GetNJets(); jetIndex++){
        jetPt = jetReader->GetJetPt(jetIndex);
        jetPhi = jetReader->GetJetPhi(jetIndex);
        jetEta = jetReader->GetJetEta(jetIndex);
        if(jetIndex == highestIndex) continue; // Do not consider leading particle
        if(TMath::Abs(jetEta) >= fJetSearchEtaCut) continue; // Cut for search eta range
        if(fMinimumMaxTrackPtFraction >= jetReader->GetJetMaxTrackPt(jetIndex)/jetReader->GetJetRawPt(jetIndex)) continue; // Cut for jets with only very low pT particles
        if(fMaximumMaxTrackPtFraction <= jetReader->GetJetMaxTrackPt(jetIndex)/jetReader->GetJetRawPt(jetIndex)) continue; // Cut for jets where all the pT is taken by one track
        
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
        leadingJetPhi = jetReader->GetJetPhi(highestIndex);
        leadingJetEta = jetReader->GetJetEta(highestIndex);
        subleadingJetPhi = jetReader->GetJetPhi(secondHighestIndex);
        subleadingJetEta = jetReader->GetJetEta(secondHighestIndex);
        
        // Apply the JFF correction for leading and subleading jet pT
        nParticleFlowCandidatesInThisJet = GetNParticleFlowCandidatesInJet(jetReader,leadingJetPhi,leadingJetEta);
        leadingJetPt = fJffCorrection->GetCorrection(nParticleFlowCandidatesInThisJet,hiBin,leadingJetPt,leadingJetEta);
        nParticleFlowCandidatesInThisJet = GetNParticleFlowCandidatesInJet(jetReader,subleadingJetPhi,subleadingJetEta);
        subleadingJetPt = fJffCorrection->GetCorrection(nParticleFlowCandidatesInThisJet,hiBin,subleadingJetPt,subleadingJetEta);
        
        // If after the correction subleading jet becomes leading jet and vice versa, swap the leading/subleading info
        if(subleadingJetPt > leadingJetPt){
          swapJetPt = leadingJetPt;   leadingJetPt = subleadingJetPt;    subleadingJetPt = swapJetPt;
          swapJetPhi = leadingJetPhi; leadingJetPhi = subleadingJetPhi;  subleadingJetPhi = swapJetPhi;
          swapJetEta = leadingJetEta; leadingJetEta = subleadingJetEta;  subleadingJetEta = swapJetEta;
          swapIndex = highestIndex;   highestIndex = secondHighestIndex; secondHighestIndex = swapIndex;
        }
        
        dijetFound = true;
        dphi =  leadingJetPhi - subleadingJetPhi;
        if(dphi < 0) dphi = -dphi;
        if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
        
        // Apply dijet cuts
        if((leadingJetPt >= fJetMaximumPtCut) ||              // Maximum leading jet pT cut
           (leadingJetPt <= fLeadingJetMinPtCut) ||           // Leading jet minimum pT cut
           (subleadingJetPt <= fSubleadingJetMinPtCut) ||     // Subleading jet minimum pT cut
           (TMath::Abs(leadingJetEta) >= fJetEtaCut) ||       // Leading jet eta cut
           (TMath::Abs(subleadingJetEta) >= fJetEtaCut)||     // Subleading jet eta cut
           (TMath::Abs(dphi) <= fDeltaPhiCut)){               // DeltaPhi cut
          dijetFound = false;
        }
      } // End of dijet cuts
      
      // If a dijet is found, fill some information to fHistograms
      if(dijetFound){
        
        // Only fill the dijet histograms if selected
        if(fillEventInformationHistograms){
          fHistograms->fhEvents->Fill(DijetHistograms::kDijet);
          fHistograms->fhCentralityDijet->Fill(centrality);
        }
        if(fillMainLoopHistograms){
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
          
        }
        
        // Fill the arrays with leading and subleading jet information for correlation with tracks
        leadingJetInfo[0] = leadingJetPt;
        leadingJetInfo[1] = leadingJetPhi;
        leadingJetInfo[2] = leadingJetEta;
        subleadingJetInfo[0] = subleadingJetPt;
        subleadingJetInfo[1] = subleadingJetPhi;
        subleadingJetInfo[2] = subleadingJetEta;
        
        // Do not do the jet-track correlation if we are only filling the event information histograms
        if(fFilledHistograms == kFillOnlyEventInformation) continue;
        
        // Correlate jets with tracks in dijet events
        CorrelateTracksAndJets(trackReader,leadingJetInfo,subleadingJetInfo,DijetHistograms::kSameEvent);
        
        // Do event mixing
        if(mixEvents){
          
          // Start mixing from the first event index
          mixedEventIndex = firstMixingEvent;
          eventsMixed = 0;
          allEventsWentThrough = false;
          checkForDuplicates = false;   // Start checking for duplicates in the second round over the file
          skipEvent = false;
          additionalTolerance = 0;
          mixedEventIndices.clear();
          
          // Continue mixing until we have reached required number of evens
          while (eventsMixed < nMixedEventsPerDijet) {
            
            // By default, do not skip an event
            skipEvent = false;
            
            // If we have checked all the events but not found enough event to mix with, increase vz tolerance
            if(allEventsWentThrough){
              if(debugLevel > 0){
                cout << "Could only find " << eventsMixed << " events to mix with event " << iEvent << " in file " << currentFile.Data() << endl;
                cout << "Increasing vz tolerance by 0.25" << endl;
              }
              
              additionalTolerance += 0.25;
              allEventsWentThrough = false;
              checkForDuplicates = true;
            }
            
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
            
            // Do not mix with events used in the previous iteration over the file
            if(checkForDuplicates){
              for(Int_t iMixedEvent : mixedEventIndices) {
                if(iMixedEvent == mixedEventIndex) skipEvent = true;
              }
            }
            if(skipEvent) continue;
            
            // Match vz and hiBin between the dijet event and mixed event
            if(TMath::Abs(mixedEventVz.at(mixedEventIndex) - vz) > (mixingVzTolerance + additionalTolerance)) continue;
            if(mixedEventHiBin.at(mixedEventIndex) != hiBin) continue;
            
            // If match vz and hiBin, then load the event from the mixed event tree and mark that we have mixed this event
            mixedEventReader->GetEvent(mixedEventIndex);
            mixedEventIndices.push_back(mixedEventIndex);
            
            // Do the correlations with the dijet from current event and track from mixing event
            CorrelateTracksAndJets(mixedEventReader,leadingJetInfo,subleadingJetInfo,DijetHistograms::kMixedEvent);
            eventsMixed++;
            
          } // While loop for finding events to mix
          
          // For the next event, start mixing the events from where we were left with in the previous event
          firstMixingEvent = mixedEventIndex;
          
        } // Event mixing
      } // Dijet in event
      
    } // Event loop
        
    // Close the input files after the event has been read
    inputFile->Close();
    if(useDifferentReaderFotJetsAndTracks) copyInputFile->Close();
    mixedEventFile->Close();
    
  } // File loop
  
  // When we are done, delete ForestReaders
  delete jetReader;
  if(fMcCorrelationType == kRecoGen || fMcCorrelationType == kGenReco) delete trackReader;
  if(mixEvents) delete mixedEventReader;
  delete mixedEventRandomizer;
}

/*
 * Method for all jet-track correlations
 *
 *  ForestReader *jetReader = Forest from which the track information is read
 *  Double_t leadingJetInfo[3] = Array containing leading jet pT, phi and eta
 *  Double_t subleadingJetInfo[3] = Array containing subleading jet pT, phi and eta
 *  Int_t correlationType = DijetHistograms::kSameEvent for same event correlations, DijetHistograms::kMixedEvent for mixed event correlations
 */
void DijetAnalyzer::CorrelateTracksAndJets(ForestReader *trackReader, Double_t leadingJetInfo[3], Double_t subleadingJetInfo[3], Int_t correlationType){
  
  // Define a filler for THnSparses
  Double_t fillerJetTrack[6];
  Double_t fillerTrack[5];
  
  // Select which histograms to fill
  Bool_t fillTrackHistograms = (fFilledHistograms == kFillAll || fFilledHistograms == kFillAllButJetTrack);
  Bool_t fillJetTrackCorrelationHistograms = (fFilledHistograms == kFillAll || fFilledHistograms == kFillJetTrack);
  
  // Event information
  Int_t hiBin = trackReader->GetHiBin();
  Double_t centrality = trackReader->GetCentrality();
  
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
  Int_t nTracks = trackReader->GetNTracks();
  for(Int_t iTrack = 0; iTrack <nTracks; iTrack++){
    
    // Only look at charged tracks
    if(trackReader->GetTrackCharge(iTrack) == 0) continue;
    if(!PassSubeventCut(trackReader->GetTrackSubevent(iTrack))) continue;
    
    trackPt = trackReader->GetTrackPt(iTrack);
    trackPhi = trackReader->GetTrackPhi(iTrack);
    trackEta = trackReader->GetTrackEta(iTrack);
    trackEt = (trackReader->GetTrackEnergyEcal(iTrack)+trackReader->GetTrackEnergyHcal(iTrack))/TMath::CosH(trackEta);
    
    //  ==== Apply cuts for tracks and collect information on how much track are cut in each step ====
    
    // Only fill the track cut histograms for same event data
    if(correlationType == DijetHistograms::kSameEvent && fillTrackHistograms) fHistograms->fhTrackCuts->Fill(DijetHistograms::kAllTracks);
    
    // Cut for track pT
    if(trackPt <= fTrackMinPtCut) continue;                     // Minimum pT cut
    if(trackPt >= fJetMaximumPtCut) continue;                   // Maximum pT cut (same as for leading jets)
    if(correlationType == DijetHistograms::kSameEvent && fillTrackHistograms) fHistograms->fhTrackCuts->Fill(DijetHistograms::kPtCuts);
    
    // Cut for track eta
    if(TMath::Abs(trackEta) >= fTrackEtaCut) continue;          // Eta cut
    if(correlationType == DijetHistograms::kSameEvent && fillTrackHistograms) fHistograms->fhTrackCuts->Fill(DijetHistograms::kEtaCut);
    
    // Cut for high purity
    if(!trackReader->GetTrackHighPurity(iTrack)) continue;     // High purity cut
    if(correlationType == DijetHistograms::kSameEvent && fillTrackHistograms) fHistograms->fhTrackCuts->Fill(DijetHistograms::kHighPurity);
    
    // Cut for relative error for track pT
    if(trackReader->GetTrackPtError(iTrack)/trackPt >= fMaxTrackPtRelativeError) continue; // Cut for track pT relative error
    if(correlationType == DijetHistograms::kSameEvent && fillTrackHistograms) fHistograms->fhTrackCuts->Fill(DijetHistograms::kPtError);
    
    // Cut for track distance from primary vertex
    if(trackReader->GetTrackVertexDistanceZ(iTrack)/trackReader->GetTrackVertexDistanceZError(iTrack) >= fMaxTrackDistanceToVertex) continue; // Mysterious cut about track proximity to vertex in z-direction
    if(trackReader->GetTrackVertexDistanceXY(iTrack)/trackReader->GetTrackVertexDistanceXYError(iTrack) >= fMaxTrackDistanceToVertex) continue; // Mysterious cut about track proximity to vertex in xy-direction
    if(correlationType == DijetHistograms::kSameEvent && fillTrackHistograms) fHistograms->fhTrackCuts->Fill(DijetHistograms::kVertexDistance);
    
    // Cut for energy deposition in calorimeters for high pT tracks
    if(!(trackPt < fCalorimeterSignalLimitPt || (trackEt >= fHighPtEtFraction*trackPt))) continue;  // For high pT tracks, require signal also in calorimeters
    if(correlationType == DijetHistograms::kSameEvent && fillTrackHistograms) fHistograms->fhTrackCuts->Fill(DijetHistograms::kCaloSignal);
    
    // Cuts for track reconstruction quality
    if(trackReader->GetTrackChi2(iTrack)/(1.0*trackReader->GetNTrackDegreesOfFreedom(iTrack))/(1.0*trackReader->GetNHitsTrackerLayer(iTrack)) >= fChi2QualityCut) continue; // Track reconstruction quality cut
    if(trackReader->GetNHitsTrack(iTrack) < fMinimumTrackHits) continue; // Cut for minimum number of hits per track
    if(correlationType == DijetHistograms::kSameEvent && fillTrackHistograms) fHistograms->fhTrackCuts->Fill(DijetHistograms::kReconstructionQuality);
    
    //     ==== Track cuts done ====
    
    // Calculate minimum distance of a track to closest jet. This is needed for track efficiency correction
    trackRMin = 666;   // Initialize the minimum distance to a jet to some very large value
    
    // No efficiency correction for generator level tracks
    if(fMcCorrelationType == kRecoGen || fMcCorrelationType == kGenGen){
      trackEfficiencyCorrection = 1;
    } else { // Apply the correction for real data and reconstructed Monte Carlo tracks
      
      // Need distance to nearest jet for the efficiency correction
      for(Int_t iJet = 0; iJet < trackReader->GetNJets(); iJet++){
        
        // Require the same jet quality cuts as when searching for dijets
        if(TMath::Abs(trackReader->GetJetEta(iJet)) >= fJetSearchEtaCut) continue; // Require jet eta to be in the search range for dijets
        if(trackReader->GetJetPt(iJet) <= fSubleadingJetMinPtCut) continue; // Only consider jets that pass the pT cut for subleading jets
        if(fMinimumMaxTrackPtFraction >= trackReader->GetJetMaxTrackPt(iJet)/trackReader->GetJetRawPt(iJet)) continue; // Cut for jets with only very low pT particles
        if(fMaximumMaxTrackPtFraction <= trackReader->GetJetMaxTrackPt(iJet)/trackReader->GetJetRawPt(iJet)) continue; // Cut for jets where all the pT is taken by one track
        // if(TMath::Abs(chargedSum[k]/rawpt[k]) < 0.01) continue; // Jet quality cut from readme file. TODO: Check if should be applied
        
        // Note: The ACos(Cos(jetPhi-trackPhi)) structure transforms deltaPhi to interval [0,Pi]
        trackR = TMath::Power(trackReader->GetJetEta(iJet)-trackEta,2)+TMath::Power(TMath::ACos(TMath::Cos(trackReader->GetJetPhi(iJet)-trackPhi)),2);
        if(trackRMin*trackRMin>trackR) trackRMin=TMath::Power(trackR,0.5);
        
      } // Loop for calculating Rmin
      
      // Get the track efficiency corrections
      trackEfficiencyCorrection = fTrackCorrection->getTrkCorr(trackPt, trackEta, trackPhi, hiBin, trackRMin);
    } // Track efficiency correction if
    
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
    
    // Fill track histograms if selected to do so
    if(fillTrackHistograms){
      fillerTrack[0] = trackPt;                    // Axis 0: Track pT
      fillerTrack[1] = trackPhi;                   // Axis 1: Track phi
      fillerTrack[2] = trackEta;                   // Axis 2: Track eta
      fillerTrack[3] = centrality;                 // Axis 3: Centrality
      fillerTrack[4] = correlationType;            // Axis 4: Correlation type (same or mixed event)
      fHistograms->fhTrack->Fill(fillerTrack,trackEfficiencyCorrection);  // Fill the track histogram
      fHistograms->fhTrackUncorrected->Fill(fillerTrack);                 // Fill the uncorrected track histogram
    }
    
    // Fill all the jet-track correlation histograms if selected
    if(fillJetTrackCorrelationHistograms){
      
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
    }
    
  } // Loop over tracks
}
/*
 * Get the proper vz weighting depending on analyzed system
 *
 *  Arguments:
 *   const Double_t vz = Vertex z position for the event
 *
 *   return: Multiplicative correction factor for vz
 */
Double_t DijetAnalyzer::GetVzWeight(const Double_t vz) const{
  if(fDataType == ForestReader::kPp || fDataType == ForestReader::kPbPb) return 1;  // No correction for real data
  if(fDataType == ForestReader::kPpMC || fDataType == ForestReader::kLocalTest) return 1.0/fVzWeightFunction->Eval(vz); // Weight for pp MC
  if(fDataType == ForestReader::kPbPbMC) return fVzWeightFunction->Eval(vz); // Weight for PbPb MC
  return -1; // Return crazy value for unknown data types, so user will not miss it
}

/*
 * Get the proper centrality weighting depending on analyzed system
 *
 *  Arguments:
 *   const Int_t hiBin = CMS hiBin
 *
 *   return: Multiplicative correction factor for the given CMS hiBin
 */
Double_t DijetAnalyzer::GetCentralityWeight(const Int_t hiBin) const{
  if(fDataType != ForestReader::kPbPbMC) return 1;  // Centrality weighting only for PbPb MC
  return (hiBin < 194) ? fCentralityWeightFunction->Eval(hiBin) : 1;  // No weighting for the most peripheral centrality bins
}

/*
 * Get the proper pT hat weighting for MC
 *
 *  Arguments:
 *   const Int_t ptHat = pT hat value in the event
 *
 *   return: Multiplicative correction factor for the given pT hat
 */
Double_t DijetAnalyzer::GetPtHatWeight(const Int_t ptHat) const{
  if(fDataType == ForestReader::kPp || fDataType == ForestReader::kPbPb || fDataType == ForestReader::kLocalTest) return 1;
  
  // The event are counted for specific bin edges, so those edges will be used for correction
  const Int_t nBins = 9;
  Int_t correctionBinEdges[nBins] = {30, 50, 80, 120, 170, 220, 280, 370, 460};
  
  // Search the correct bin for the given pT hat value
  Int_t currentBin = -1;
  for(Int_t iBin = 0; iBin < nBins; iBin++){
    if(ptHat < correctionBinEdges[iBin]){
      currentBin = iBin;
      break;
    }
  }
  
  // Cross sections for each bin are given in the twiki https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiForest2015
  Double_t crossSections[nBins+1] = {2.403,3.378e-2,3.778e-3,4.423e-4,6.147e-5,1.018e-5,2.477e-6,6.160e-7,1.088e-7,2.537e-8};
  
  // Number of events for different pT hat bins in the forest file list ppMC_Pythia6_forest_5TeV.txt
  Int_t ppMcEvents[nBins] = {0,444104,322347,383263,468748,447937,259209,234447,39275};
  
  // Return the weight for pp
  if(fDataType == ForestReader::kPpMC) {
    if(ppMcEvents[currentBin] == 0){ // This should never happen
      cout << "WARNING! No events found for pT hat bin " << currentBin << endl;
      return 0;
    }
    return (crossSections[currentBin]-crossSections[currentBin+1])/ppMcEvents[currentBin];
  }
  
  // Number of events for different pT hat bins in the forest file list PbPbMC_5TeVPythia6+Hydjet_forests.txt
  Int_t PbPbMcEvents[nBins] = {0,0,0,2571563,2850815,2680567,2891375,781744,129417};
  
  // Return the weight for PbPb
  if(PbPbMcEvents[currentBin] == 0) { // This should never happen
    cout << "WARNING! No events found for pT hat bin " << currentBin << endl;
    return 0;
  }
  return (crossSections[currentBin]-crossSections[currentBin+1])/PbPbMcEvents[currentBin];
}

/*
 * Check is a track passes the required subevent cut
 *
 *  Arguments:
 *   const Int_t subeventIndex = Subevent index for the track in consideration
 *
 *  return: true if subevent cut is passes, false if not
 */
Bool_t DijetAnalyzer::PassSubeventCut(const Int_t subeventIndex) const{
  if(fSubeventCut == kSubeventAny) return true;
  if(fSubeventCut == kSubeventZero && subeventIndex == 0) return true;
  if(fSubeventCut == kSubeventNonZero && subeventIndex > 0) return true;
  return false;
}

/*
 * Getter for dijet histograms
 */
DijetHistograms* DijetAnalyzer::GetHistograms() const{
  return fHistograms;
}
