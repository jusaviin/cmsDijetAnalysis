// Implementation for MixingForestReader

// Own includes
#include "MixingForestReader.h"

/*
 * Default constructor
 */
MixingForestReader::MixingForestReader() :
  ForestReader(),
  fHeavyIonTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fnTracksBranch(0),
  fTrackAlgorithmBranch(0),
  fTrackMVABranch(0),
  fTrackPtArray(),
  fTrackPtErrorArray(),
  fTrackPhiArray(),
  fTrackEtaArray(),
  fHighPurityTrackArray(),
  fTrackVertexDistanceZArray(),
  fTrackVertexDistanceZErrorArray(),
  fTrackVertexDistanceXYArray(),
  fTrackVertexDistanceXYErrorArray(),
  fTrackChi2Array(),
  fnTrackDegreesOfFreedomArray(),
  fnHitsTrackerLayerArray(),
  fnHitsTrackArray(),
  fTrackEnergyEcalArray(),
  fTrackEnergyHcalArray(),
  fTrackAlgorithmArray(),
  fTrackMVAArray()
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
 */
MixingForestReader::MixingForestReader(Int_t dataType, Int_t readMode, Int_t jetType, Int_t jetAxis, Bool_t matchJets) :
  ForestReader(dataType,readMode,jetType,jetAxis,matchJets),
  fHeavyIonTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fnTracksBranch(0),
  fTrackAlgorithmBranch(0),
  fTrackMVABranch(0),
  fTrackPtArray(),
  fTrackPtErrorArray(),
  fTrackPhiArray(),
  fTrackEtaArray(),
  fHighPurityTrackArray(),
  fTrackVertexDistanceZArray(),
  fTrackVertexDistanceZErrorArray(),
  fTrackVertexDistanceXYArray(),
  fTrackVertexDistanceXYErrorArray(),
  fTrackChi2Array(),
  fnTrackDegreesOfFreedomArray(),
  fnHitsTrackerLayerArray(),
  fnHitsTrackArray(),
  fTrackEnergyEcalArray(),
  fTrackEnergyHcalArray(),
  fTrackAlgorithmArray(),
  fTrackMVAArray()
{
  // Custom constructor
  
}

/*
 * Copy constructor
 */
MixingForestReader::MixingForestReader(const MixingForestReader& in) :
  ForestReader(in),
  fHeavyIonTree(in.fHeavyIonTree),
  fSkimTree(in.fSkimTree),
  fTrackTree(in.fTrackTree),
  fnTracksBranch(in.fnTracksBranch),
  fTrackAlgorithmBranch(in.fTrackAlgorithmBranch),
  fTrackMVABranch(in.fTrackMVABranch)
{
  // Copy constructor  
  for(Int_t i = 0; i < fnMaxTrack; i++){
    fTrackPtArray[i] = in.fTrackPtArray[i];
    fTrackPtErrorArray[i] = in.fTrackPtErrorArray[i];
    fTrackPhiArray[i] = in.fTrackPhiArray[i];
    fTrackEtaArray[i] = in.fTrackEtaArray[i];
    fHighPurityTrackArray[i] = in.fHighPurityTrackArray[i];
    fTrackVertexDistanceZArray[i] = in.fTrackVertexDistanceZArray[i];
    fTrackVertexDistanceZErrorArray[i] = in.fTrackVertexDistanceZErrorArray[i];
    fTrackVertexDistanceXYArray[i] = in.fTrackVertexDistanceXYArray[i];
    fTrackVertexDistanceXYErrorArray[i] = in.fTrackVertexDistanceXYErrorArray[i];
    fTrackChi2Array[i] = in.fTrackChi2Array[i];
    fnTrackDegreesOfFreedomArray[i] = in.fnTrackDegreesOfFreedomArray[i];
    fnHitsTrackerLayerArray[i] = in.fnHitsTrackerLayerArray[i];
    fnHitsTrackArray[i] = in.fnHitsTrackArray[i];
    fTrackEnergyEcalArray[i] = in.fTrackEnergyEcalArray[i];
    fTrackEnergyHcalArray[i] = in.fTrackEnergyHcalArray[i];
    fTrackAlgorithmArray[i] = in.fTrackAlgorithmArray[i];
    fTrackMVAArray[i] = in.fTrackMVAArray[i];
  }
}

/*
 * Assignment operator
 */
MixingForestReader& MixingForestReader::operator=(const MixingForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  ForestReader::operator=(in);
  
  fHeavyIonTree = in.fHeavyIonTree;
  fSkimTree = in.fSkimTree;
  fTrackTree = in.fTrackTree;
  fnTracksBranch = in.fnTracksBranch;
  fTrackAlgorithmBranch = in.fTrackAlgorithmBranch;
  fTrackMVABranch = in.fTrackMVABranch;
  
  for(Int_t i = 0; i < fnMaxTrack; i++){
    fTrackPtArray[i] = in.fTrackPtArray[i];
    fTrackPtErrorArray[i] = in.fTrackPtErrorArray[i];
    fTrackPhiArray[i] = in.fTrackPhiArray[i];
    fTrackEtaArray[i] = in.fTrackEtaArray[i];
    fHighPurityTrackArray[i] = in.fHighPurityTrackArray[i];
    fTrackVertexDistanceZArray[i] = in.fTrackVertexDistanceZArray[i];
    fTrackVertexDistanceZErrorArray[i] = in.fTrackVertexDistanceZErrorArray[i];
    fTrackVertexDistanceXYArray[i] = in.fTrackVertexDistanceXYArray[i];
    fTrackVertexDistanceXYErrorArray[i] = in.fTrackVertexDistanceXYErrorArray[i];
    fTrackChi2Array[i] = in.fTrackChi2Array[i];
    fnTrackDegreesOfFreedomArray[i] = in.fnTrackDegreesOfFreedomArray[i];
    fnHitsTrackerLayerArray[i] = in.fnHitsTrackerLayerArray[i];
    fnHitsTrackArray[i] = in.fnHitsTrackArray[i];
    fTrackEnergyEcalArray[i] = in.fTrackEnergyEcalArray[i];
    fTrackEnergyHcalArray[i] = in.fTrackEnergyHcalArray[i];
    fTrackAlgorithmArray[i] = in.fTrackAlgorithmArray[i];
    fTrackMVAArray[i] = in.fTrackMVAArray[i];
  }
  
  return *this;
}

/*
 * Destructor
 */
MixingForestReader::~MixingForestReader(){
  // destructor
}

/*
 * Initialization, meaning that the branches are connected to the tree
 */
void MixingForestReader::Initialize(){
  
  // Connect the branches of the heavy ion tree
  fHeavyIonTree->SetBranchStatus("*",0);
  fHeavyIonTree->SetBranchStatus("vz",1);
  fHeavyIonTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  fHeavyIonTree->SetBranchStatus("hiBin",1);
  fHeavyIonTree->SetBranchAddress("hiBin",&fHiBin,&fHiBinBranch);
  fPtHat = 0; // pT hat information not needed for mixing
  fEventWeight = 0; // Event weight not needed for mixing
    
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
  fTrackTree->SetBranchStatus("trkPt",1);
  fTrackTree->SetBranchAddress("trkPt",&fTrackPtArray,&fTrackPtBranch);
  fTrackTree->SetBranchStatus("trkPtError",1);
  fTrackTree->SetBranchAddress("trkPtError",&fTrackPtErrorArray,&fTrackPtErrorBranch);
  fTrackTree->SetBranchStatus("trkPhi",1);
  fTrackTree->SetBranchAddress("trkPhi",&fTrackPhiArray,&fTrackPhiBranch);
  fTrackTree->SetBranchStatus("trkEta",1);
  fTrackTree->SetBranchAddress("trkEta",&fTrackEtaArray,&fTrackEtaBranch);
  fTrackTree->SetBranchStatus("nTrk",1);
  fTrackTree->SetBranchAddress("nTrk",&fnTracks,&fnTracksBranch);
  fTrackTree->SetBranchStatus("highPurity",1);
  fTrackTree->SetBranchAddress("highPurity",&fHighPurityTrackArray,&fHighPurityTrackBranch);
  fTrackTree->SetBranchStatus("trkDz1",1);
  fTrackTree->SetBranchAddress("trkDz1",&fTrackVertexDistanceZArray,&fTrackVertexDistanceZBranch);
  fTrackTree->SetBranchStatus("trkDzError1",1);
  fTrackTree->SetBranchAddress("trkDzError1",&fTrackVertexDistanceZErrorArray,&fTrackVertexDistanceZErrorBranch);
  fTrackTree->SetBranchStatus("trkDxy1",1);
  fTrackTree->SetBranchAddress("trkDxy1",&fTrackVertexDistanceXYArray,&fTrackVertexDistanceXYBranch);
  fTrackTree->SetBranchStatus("trkDxyError1",1);
  fTrackTree->SetBranchAddress("trkDxyError1",&fTrackVertexDistanceXYErrorArray,&fTrackVertexDistanceXYErrorBranch);
  fTrackTree->SetBranchStatus("trkChi2",1);
  fTrackTree->SetBranchAddress("trkChi2",&fTrackChi2Array,&fTrackChi2Branch);
  fTrackTree->SetBranchStatus("trkNdof",1);
  fTrackTree->SetBranchAddress("trkNdof",&fnTrackDegreesOfFreedomArray,&fnTrackDegreesOfFreedomBranch);
  fTrackTree->SetBranchStatus("trkNlayer",1);
  fTrackTree->SetBranchAddress("trkNlayer",&fnHitsTrackerLayerArray,&fnHitsTrackerLayerBranch);
  fTrackTree->SetBranchStatus("trkNHit",1);
  fTrackTree->SetBranchAddress("trkNHit",&fnHitsTrackArray,&fnHitsTrackBranch);
  fTrackTree->SetBranchStatus("pfEcal",1);
  fTrackTree->SetBranchAddress("pfEcal",&fTrackEnergyEcalArray,&fTrackEnergyEcalBranch);
  fTrackTree->SetBranchStatus("pfHcal",1);
  fTrackTree->SetBranchAddress("pfHcal",&fTrackEnergyHcalArray,&fTrackEnergyHcalBranch);
  
  // Additional information needed for 2018 track cuts
  fTrackTree->SetBranchStatus("trkAlgo",1);
  fTrackTree->SetBranchAddress("trkAlgo",&fTrackAlgorithmArray,&fTrackAlgorithmBranch);
  fTrackTree->SetBranchStatus("trkMVA",1);
  fTrackTree->SetBranchAddress("trkMVA",&fTrackMVAArray,&fTrackMVABranch);
  
}

/*
 * Connect a new tree to the reader
 */
void MixingForestReader::ReadForestFromFile(TFile* inputFile){
  std::vector<TString> fileList;
  fileList.push_back(inputFile->GetName());
  ReadForestFromFileList(fileList);
}

/*
 * Connect a new tree to the reader
 */
void MixingForestReader::ReadForestFromFileList(std::vector<TString> inputFileList){
  
  // Connect a trees from the file to the reader
  fHeavyIonTree = new TChain("hiEvtAnalyzer/HiTree");
  fSkimTree = new TChain("skimanalysis/HltTree");
  
  // The track tree and the particle flow candidate tree have different names for differant datasets
  if(fDataType == kPp || fDataType == kPpMC || fDataType == kLocalTest || fReadMode > 2000){
    fTrackTree = new TChain("ppTrack/trackTree");
  } else {
    fTrackTree = new TChain("anaTrack/trackTree"); // 2015 PbPb syntax
  }
  
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
void MixingForestReader::BurnForest(){
  fHeavyIonTree->Delete();
  fSkimTree->Delete();
  fTrackTree->Delete();
}

/*
 * Load an event to memory
 */
void MixingForestReader::GetEvent(Int_t nEvent){
  fHeavyIonTree->GetEntry(nEvent);
  fSkimTree->GetEntry(nEvent);
  fTrackTree->GetEntry(nEvent);
}

// Getter for number of events in the tree
Int_t MixingForestReader::GetNEvents() const{
  return fHeavyIonTree->GetEntries();
}

// Getter for jet pT
Float_t MixingForestReader::GetJetPt(Int_t iJet) const{
  return 0; // No jet information required for mixing
}

// Getter for jet phi
Float_t MixingForestReader::GetJetPhi(Int_t iJet) const{
  return 0; // No jet information required for mixing
}

// Getter for jet eta
Float_t MixingForestReader::GetJetEta(Int_t iJet) const{
   return 0; // No jet information required for mixing
}

// Getter for jet raw pT
Float_t MixingForestReader::GetJetRawPt(Int_t iJet) const{
   return 0; // No jet information required for mixing
}

// Getter for maximum track pT inside a jet
Float_t MixingForestReader::GetJetMaxTrackPt(Int_t iJet) const{
   return 0; // No jet information required for mixing
}

// Getter for track pT
Float_t MixingForestReader::GetTrackPt(Int_t iTrack) const{
  return fTrackPtArray[iTrack];
}

// Getter for track pT error
Float_t MixingForestReader::GetTrackPtError(Int_t iTrack) const{
  return fTrackPtErrorArray[iTrack];
}

// Getter for track phi
Float_t MixingForestReader::GetTrackPhi(Int_t iTrack) const{
  return fTrackPhiArray[iTrack];
}

// Getter for track eta
Float_t MixingForestReader::GetTrackEta(Int_t iTrack) const{
  return fTrackEtaArray[iTrack];
}

// Getter for high purity of the track
Bool_t MixingForestReader::GetTrackHighPurity(Int_t iTrack) const{
  return fHighPurityTrackArray[iTrack];
}

// Getter for track distance from primary vertex in z-direction
Float_t MixingForestReader::GetTrackVertexDistanceZ(Int_t iTrack) const{
  return fTrackVertexDistanceZArray[iTrack];
}

// Getter for error of track distance from primary vertex in z-direction
Float_t MixingForestReader::GetTrackVertexDistanceZError(Int_t iTrack) const{
  return fTrackVertexDistanceZErrorArray[iTrack];
}

// Getter for track distance from primary vertex in xy-direction
Float_t MixingForestReader::GetTrackVertexDistanceXY(Int_t iTrack) const{
  return fTrackVertexDistanceXYArray[iTrack];
}

// Getter for error of track distance from primary vertex in xy-direction
Float_t MixingForestReader::GetTrackVertexDistanceXYError(Int_t iTrack) const{
  return fTrackVertexDistanceXYErrorArray[iTrack];
}

// Getter for track chi2 value from reconstruction fit
Float_t MixingForestReader::GetTrackChi2(Int_t iTrack) const{
  return fTrackChi2Array[iTrack];
}

// Getter for number of degrees of freedom in reconstruction fit
Int_t MixingForestReader::GetNTrackDegreesOfFreedom(Int_t iTrack) const{
  return fnTrackDegreesOfFreedomArray[iTrack];
}

// Getter for number of hits in tracker layers
Int_t MixingForestReader::GetNHitsTrackerLayer(Int_t iTrack) const{
  return fnHitsTrackerLayerArray[iTrack];
}

// Getter for number of hits for the track
Int_t MixingForestReader::GetNHitsTrack(Int_t iTrack) const{
  return fnHitsTrackArray[iTrack];
}

// Getter for track energy in ECal
Float_t MixingForestReader::GetTrackEnergyEcal(Int_t iTrack) const{
  return fTrackEnergyEcalArray[iTrack];
}

// Getter for track energy in HCal
Float_t MixingForestReader::GetTrackEnergyHcal(Int_t iTrack) const{
  return fTrackEnergyHcalArray[iTrack];
}

// Getter for track charge. Charge is only relevant for generator level tracks.
Int_t MixingForestReader::GetTrackCharge(Int_t iTrack) const{
  return 1;
}

// Getter for track subevent index. Relevant only for generator level tracks.
Int_t MixingForestReader::GetTrackSubevent(Int_t iTrack) const{
  return -1;
}

// Getter for track MC status. Relevant only for generator level tracks.
Int_t MixingForestReader::GetTrackMCStatus(Int_t iTrack) const{
  return 1;
}

// Getter for track algorithm
Int_t MixingForestReader::GetTrackAlgorithm(Int_t iTrack) const{
  return fTrackAlgorithmArray[iTrack];
}

// Getter for track MVA
Float_t MixingForestReader::GetTrackMVA(Int_t iTrack) const{
  return fTrackMVAArray[iTrack];
}

// Check if generator level jet has a matching reconstructed jet
Bool_t MixingForestReader::HasMatchingJet(Int_t iJet) const{
  return false;  // No matching is done for mixing
}

// Getter for matched generator level jet pT
Float_t MixingForestReader::GetMatchedPt(Int_t iJet) const{
  return 0; // No matching is done for mixing
}

// Parton flavor for the parton initiating the jet
Int_t MixingForestReader::GetPartonFlavor(Int_t iJet) const{
  return 0; // No matching is done for mixing
}

// Get the eta of the matched generator level jet
Float_t MixingForestReader::GetMatchedEta(Int_t iJet) const{
  return 0; // No matching is done for mixing
}

// Get the pT of the matched reconstructed jet
Float_t MixingForestReader::GetMatchedPhi(Int_t iJet) const{
  return 0; // No matching is done for mixing
}
