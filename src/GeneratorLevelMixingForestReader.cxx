// Implementation for GeneratorLevelMixingForestReader

// Own includes
#include "GeneratorLevelMixingForestReader.h"

/*
 * Default constructor
 */
GeneratorLevelMixingForestReader::GeneratorLevelMixingForestReader() :
  ForestReader(),
  fHeavyIonTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fTrackPtArray(0),
  fTrackPhiArray(0),
  fTrackEtaArray(0),
  fTrackChargeArray(0),
  fTrackSubeventArray(0)
{
  // Default constructor
}

/*
 * Custom constructor
 *
 *  Arguments:
 *   Int_t dataType: 0 = pp, 1 = PbPb, 2 = pp MC, 3 = PbPb MC, 4 = Local Test
 *   Int_t readMode: 0 = Regular forests, 1 = Official PYTHIA8 forest
 *   Int_t jetType: 0 = Calo jets, 1 = PF jets
 *   Int_t jetAxis: 0 = Anti-kT axis, 1 = Leading particle flow candidate axis, 2 = WTA axis
 *   Bool_t matchJets: True = Do matching for reco and gen jets. False = Do not require matching
 *   Bool_t doEventPlane: Read the event plane branches from the tree. Branches not included in older trees.
 *   Bool_t minimumBiasMode: There is no HLT tree in minimum bias files, it cannot be loaded for those runs
 */
GeneratorLevelMixingForestReader::GeneratorLevelMixingForestReader(Int_t dataType, Int_t readMode, Int_t jetType, Int_t jetAxis, Bool_t matchJets, Bool_t doEventPlane, Bool_t minimumBiasMode) :
  ForestReader(dataType,readMode,jetType,jetAxis,matchJets,doEventPlane,minimumBiasMode),
  fHeavyIonTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fTrackPtArray(0),
  fTrackPhiArray(0),
  fTrackEtaArray(0),
  fTrackChargeArray(0),
  fTrackSubeventArray(0)
{
  // Custom constructor
  
}

/*
 * Copy constructor
 */
GeneratorLevelMixingForestReader::GeneratorLevelMixingForestReader(const GeneratorLevelMixingForestReader& in) :
  ForestReader(in),
  fHeavyIonTree(in.fHeavyIonTree),
  fSkimTree(in.fSkimTree),
  fTrackTree(in.fTrackTree),
  fTrackPtArray(in.fTrackPtArray),
  fTrackPhiArray(in.fTrackPhiArray),
  fTrackEtaArray(in.fTrackEtaArray),
  fTrackChargeArray(in.fTrackChargeArray),
  fTrackSubeventArray(in.fTrackSubeventArray)
{
  // Copy constructor  

}

/*
 * Assignment operator
 */
GeneratorLevelMixingForestReader& GeneratorLevelMixingForestReader::operator=(const GeneratorLevelMixingForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  ForestReader::operator=(in);
  
  fHeavyIonTree = in.fHeavyIonTree;
  fSkimTree = in.fSkimTree;
  fTrackTree = in.fTrackTree;
  
  fTrackPtArray = in.fTrackPtArray;
  fTrackPhiArray = in.fTrackPhiArray;
  fTrackEtaArray = in.fTrackEtaArray;
  fTrackChargeArray = in.fTrackChargeArray;
  fTrackSubeventArray = in.fTrackSubeventArray;
  
  return *this;
}

/*
 * Destructor
 */
GeneratorLevelMixingForestReader::~GeneratorLevelMixingForestReader(){
  // destructor
}

/*
 * Initialization, meaning that the branches are connected to the tree
 */
void GeneratorLevelMixingForestReader::Initialize(){
  
  // Connect the branches of the heavy ion tree
  fHeavyIonTree->SetBranchStatus("*",0);
  fHeavyIonTree->SetBranchStatus("vz",1);
  fHeavyIonTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  fHeavyIonTree->SetBranchStatus("hiBin",1);
  fHeavyIonTree->SetBranchAddress("hiBin",&fHiBin,&fHiBinBranch);
  fPtHat = 1;       // pT hat information not needed for mixing
  fEventWeight = 1; // Event weight not needed for mixing
    
  // Mixing from minimum bias file -> no trigger selection
  fCaloJetFilterBit = 1;  // No filter for local test of PbPb MC forest
  fCaloJetFilterBitPrescale = 1; // Set the prescaled filter bit to 1. Only relevant for minimum bias PbPb (data skim)
  
  // Connect the branches to the skim tree (different for pp and PbPb data and Monte Carlo)
  fSkimTree->SetBranchStatus("*",0);
  if(fDataType == kPp || fDataType == kPpMC){ // pp data or MC
    fSkimTree->SetBranchStatus("pPAprimaryVertexFilter",1);
    fSkimTree->SetBranchAddress("pPAprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fSkimTree->SetBranchStatus("pBeamScrapingFilter",1);
    fSkimTree->SetBranchAddress("pBeamScrapingFilter",&fBeamScrapingFilterBit,&fBeamScrapingBranch);
    fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
    fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fCollisionEventSelectionFilterBit = 1;  // No collision event selection filter for pp
    fHfCoincidenceFilterBit = 1; // No HF energy coincidence requirement for pp
    fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pp
  } else if (fDataType == kPbPb || fDataType == kPbPbMC){ // PbPb data or MC
    
    fSkimTree->SetBranchStatus("pprimaryVertexFilter",1);
    fSkimTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
    fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    
    // Branches with 2018 syntax
    if(fReadMode > 2000){
      fSkimTree->SetBranchStatus("collisionEventSelectionAODv2",1); // 2018 syntax (or not v2?)
      fSkimTree->SetBranchAddress("collisionEventSelectionAODv2", &fCollisionEventSelectionFilterBit, &fCollisionEventSelectionBranch); // 2018 syntax (or not v2?)
      fSkimTree->SetBranchStatus("phfCoincFilter3Th3",1); // 2018 syntax (of Th4?)
      fSkimTree->SetBranchAddress("phfCoincFilter3Th3", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch); // 2018 syntax
    } else { // Banches with 2015 syntax
      fSkimTree->SetBranchStatus("pcollisionEventSelection",1); // 2015 syntax
      fSkimTree->SetBranchAddress("pcollisionEventSelection", &fCollisionEventSelectionFilterBit, &fCollisionEventSelectionBranch); // 2015 syntax
      fSkimTree->SetBranchStatus("phfCoincFilter3",1); // 2015 syntax
      fSkimTree->SetBranchAddress("phfCoincFilter3", &fHfCoincidenceFilterBit ,&fHfCoincidenceBranch); // 2015 syntax
    }
    
    fSkimTree->SetBranchStatus("pclusterCompatibilityFilter",1);
    fSkimTree->SetBranchAddress("pclusterCompatibilityFilter",&fClusterCompatibilityFilterBit,&fClusterCompatibilityBranch);
    fBeamScrapingFilterBit = 1;  // No beam scraping filter for PbPb
  } else { // Local test
    fPrimaryVertexFilterBit = 1;
    fBeamScrapingFilterBit = 1;
    fCollisionEventSelectionFilterBit = 1;
    fHBHENoiseFilterBit = 1;
    fHfCoincidenceFilterBit = 1;
    fClusterCompatibilityFilterBit = 1;
  }
  
  // Connect the branches to the track tree
  fTrackTree->SetBranchStatus("*",0);
  fTrackTree->SetBranchStatus("pt",1);
  fTrackTree->SetBranchAddress("pt",&fTrackPtArray,&fTrackPtBranch);
  fTrackTree->SetBranchStatus("phi",1);
  fTrackTree->SetBranchAddress("phi",&fTrackPhiArray,&fTrackPhiBranch);
  fTrackTree->SetBranchStatus("eta",1);
  fTrackTree->SetBranchAddress("eta",&fTrackEtaArray,&fTrackEtaBranch);
  fTrackTree->SetBranchStatus("chg",1);
  fTrackTree->SetBranchAddress("chg",&fTrackChargeArray,&fTrackPtErrorBranch);  // Reuse a branch from ForestReader that is not otherwise needed here
  fTrackTree->SetBranchStatus("sube",1);
  fTrackTree->SetBranchAddress("sube",&fTrackSubeventArray,&fTrackChi2Branch);  // Reuse a branch from ForestReader that is not otherwise needed here
  
}

/*
 * Connect a new tree to the reader
 */
void GeneratorLevelMixingForestReader::ReadForestFromFile(TFile* inputFile){
  std::vector<TString> fileList;
  fileList.push_back(inputFile->GetName());
  ReadForestFromFileList(fileList);
}

/*
 * Connect a new tree to the reader
 */
void GeneratorLevelMixingForestReader::ReadForestFromFileList(std::vector<TString> inputFileList){
  
  // Connect a trees from the file to the reader
  fHeavyIonTree = new TChain("hiEvtAnalyzer/HiTree");
  fSkimTree = new TChain("skimanalysis/HltTree");
  fTrackTree = new TChain("HiGenParticleAna/hi");
  
  for(std::vector<TString>::iterator listIterator = inputFileList.begin(); listIterator != inputFileList.end(); listIterator++){
    fHeavyIonTree->Add(*listIterator);
    fSkimTree->Add(*listIterator);
    fTrackTree->Add(*listIterator);
  }
  
  Initialize();
}

/*
 * Burn the current forest.
 */
void GeneratorLevelMixingForestReader::BurnForest(){
  fHeavyIonTree->Delete();
  fSkimTree->Delete();
  fTrackTree->Delete();
}

/*
 * Load an event to memory
 */
void GeneratorLevelMixingForestReader::GetEvent(Int_t nEvent){
  fHeavyIonTree->GetEntry(nEvent);
  fSkimTree->GetEntry(nEvent);
  fTrackTree->GetEntry(nEvent);
  
  // Read the numbers of tracks for this event
  fnTracks = fTrackPtArray->size();
}

// Getter for number of events in the tree
Int_t GeneratorLevelMixingForestReader::GetNEvents() const{
  return fHeavyIonTree->GetEntries();
}

// Getter for jet pT
Float_t GeneratorLevelMixingForestReader::GetJetPt(Int_t iJet) const{
  return 0; // No jet information required for mixing
}

// Getter for jet phi
Float_t GeneratorLevelMixingForestReader::GetJetPhi(Int_t iJet) const{
  return 0; // No jet information required for mixing
}

// Getter for jet eta
Float_t GeneratorLevelMixingForestReader::GetJetEta(Int_t iJet) const{
   return 0; // No jet information required for mixing
}

// Getter for jet raw pT
Float_t GeneratorLevelMixingForestReader::GetJetRawPt(Int_t iJet) const{
   return 0; // No jet information required for mixing
}

// Getter for maximum track pT inside a jet
Float_t GeneratorLevelMixingForestReader::GetJetMaxTrackPt(Int_t iJet) const{
   return 0; // No jet information required for mixing
}

// Getter for track pT
Float_t GeneratorLevelMixingForestReader::GetTrackPt(Int_t iTrack) const{
  return fTrackPtArray->at(iTrack);
}

// Getter for track pT error
Float_t GeneratorLevelMixingForestReader::GetTrackPtError(Int_t iTrack) const{
  return 0; // No error for generator level tracks
}

// Getter for track phi
Float_t GeneratorLevelMixingForestReader::GetTrackPhi(Int_t iTrack) const{
  return fTrackPhiArray->at(iTrack);
}

// Getter for track eta
Float_t GeneratorLevelMixingForestReader::GetTrackEta(Int_t iTrack) const{
  return fTrackEtaArray->at(iTrack);
}

// Getter for high purity of the track
Bool_t GeneratorLevelMixingForestReader::GetTrackHighPurity(Int_t iTrack) const{
  return true; // Generatod tracks always have high purity
}

// Getter for track distance from primary vertex in z-direction
Float_t GeneratorLevelMixingForestReader::GetTrackVertexDistanceZ(Int_t iTrack) const{
  return 0; // The cut is on distance/error. Setting this to 0 and error to 1 always passes the cut
}

// Getter for error of track distance from primary vertex in z-direction
Float_t GeneratorLevelMixingForestReader::GetTrackVertexDistanceZError(Int_t iTrack) const{
  return 1; // The cut is on distance/error. Setting this to 1 and distance to 0 always passes the cut
}

// Getter for track distance from primary vertex in xy-direction
Float_t GeneratorLevelMixingForestReader::GetTrackVertexDistanceXY(Int_t iTrack) const{
  return 0; // The cut is on distance/error. Setting this to 0 and error to 1 always passes the cut
}

// Getter for error of track distance from primary vertex in xy-direction
Float_t GeneratorLevelMixingForestReader::GetTrackVertexDistanceXYError(Int_t iTrack) const{
  return 1; // The cut is on distance/error. Setting this to 1 and distance to 0 always passes the cut
}

// Getter for track chi2 value from reconstruction fit
Float_t GeneratorLevelMixingForestReader::GetTrackChi2(Int_t iTrack) const{
  return 1; // Track cuts disabled for generator level tracks
}

// Getter for number of degrees of freedom in reconstruction fit
Int_t GeneratorLevelMixingForestReader::GetNTrackDegreesOfFreedom(Int_t iTrack) const{
  return 1; // Track cuts disabled for generator level tracks
}

// Getter for number of hits in tracker layers
Int_t GeneratorLevelMixingForestReader::GetNHitsTrackerLayer(Int_t iTrack) const{
  return 1; // Track cuts disabled for generator level tracks
}

// Getter for number of hits for the track
Int_t GeneratorLevelMixingForestReader::GetNHitsTrack(Int_t iTrack) const{
  return 1; // Track cuts disabled for generator level tracks
}

// Getter for track energy in ECal
Float_t GeneratorLevelMixingForestReader::GetTrackEnergyEcal(Int_t iTrack) const{
  return 1; // Track cuts disabled for generator level tracks
}

// Getter for track energy in HCal
Float_t GeneratorLevelMixingForestReader::GetTrackEnergyHcal(Int_t iTrack) const{
  return 1; // Track cuts disabled for generator level tracks
}

// Getter for track charge.
Int_t GeneratorLevelMixingForestReader::GetTrackCharge(Int_t iTrack) const{
  return fTrackChargeArray->at(iTrack);
}

// Getter for track subevent index.
Int_t GeneratorLevelMixingForestReader::GetTrackSubevent(Int_t iTrack) const{
  return fTrackSubeventArray->at(iTrack);
}

// Getter for track MC status. Relevant only for some 2015 datasets.
Int_t GeneratorLevelMixingForestReader::GetTrackMCStatus(Int_t iTrack) const{
  return 1;
}

// Getter for track algorithm
Int_t GeneratorLevelMixingForestReader::GetTrackAlgorithm(Int_t iTrack) const{
  return 1; // Track cuts disabled for generator level tracks
}

// Getter for track MVA
Float_t GeneratorLevelMixingForestReader::GetTrackMVA(Int_t iTrack) const{
  return 1; // Track cuts disabled for generator level tracks
}

// Check if generator level jet has a matching reconstructed jet
Bool_t GeneratorLevelMixingForestReader::HasMatchingJet(Int_t iJet) const{
  return false;  // No matching is done for mixing
}

// Getter for matched generator level jet pT
Float_t GeneratorLevelMixingForestReader::GetMatchedPt(Int_t iJet) const{
  return 0; // No matching is done for mixing
}

// Parton flavor for the parton initiating the jet
Int_t GeneratorLevelMixingForestReader::GetPartonFlavor(Int_t iJet) const{
  return 0; // No matching is done for mixing
}

// Get the eta of the matched generator level jet
Float_t GeneratorLevelMixingForestReader::GetMatchedEta(Int_t iJet) const{
  return 0; // No matching is done for mixing
}

// Get the pT of the matched reconstructed jet
Float_t GeneratorLevelMixingForestReader::GetMatchedPhi(Int_t iJet) const{
  return 0; // No matching is done for mixing
}
