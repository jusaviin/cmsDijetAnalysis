// Implementation for SkimForestReader

// Own includes
#include "SkimForestReader.h"

/*
 * Default constructor
 */
SkimForestReader::SkimForestReader() :
  ForestReader(),
  fEventTree(0),
  fJetPtArray(0),
  fJetPhiArray(0),
  fJetEtaArray(0),
  fJetRawPtArray(0),
  fJetMaxTrackPtArray(0),
  fTrackPtArray(0),
  fTrackPtErrorArray(0),
  fTrackPhiArray(0),
  fTrackEtaArray(0),
  fHighPurityTrackArray(0),
  fTrackVertexDistanceZArray(0),
  fTrackVertexDistanceZErrorArray(0),
  fTrackVertexDistanceXYArray(0),
  fTrackVertexDistanceXYErrorArray(0),
  fTrackChi2Array(0),
  fnTrackDegreesOfFreedomArray(0),
  fnHitsTrackerLayerArray(0),
  fnHitsTrackArray(0),
  fTrackEnergyEcalArray(0),
  fTrackEnergyHcalArray(0)
{
  // Default constructor
}

/*
 * Custom constructor
 *
 *  Arguments:
 *   Int_t dataType: 0 = pp, 1 = PbPb, 2 = pp MC, 3 = PbPb MC, 4 = Local Test
 */
SkimForestReader::SkimForestReader(Int_t dataType) :
  ForestReader(dataType),
  fEventTree(0),
  fJetPtArray(0),
  fJetPhiArray(0),
  fJetEtaArray(0),
  fJetRawPtArray(0),
  fJetMaxTrackPtArray(0),
  fTrackPtArray(0),
  fTrackPtErrorArray(0),
  fTrackPhiArray(0),
  fTrackEtaArray(0),
  fHighPurityTrackArray(0),
  fTrackVertexDistanceZArray(0),
  fTrackVertexDistanceZErrorArray(0),
  fTrackVertexDistanceXYArray(0),
  fTrackVertexDistanceXYErrorArray(0),
  fTrackChi2Array(0),
  fnTrackDegreesOfFreedomArray(0),
  fnHitsTrackerLayerArray(0),
  fnHitsTrackArray(0),
  fTrackEnergyEcalArray(0),
  fTrackEnergyHcalArray(0)
{
  // Custom constructor
  
}

/*
 * Copy constructor
 */
SkimForestReader::SkimForestReader(const SkimForestReader& in) :
  ForestReader(in),
  fEventTree(in.fEventTree),
  fJetPtArray(in.fJetPtArray),
  fJetPhiArray(in.fJetPhiArray),
  fJetEtaArray(in.fJetEtaArray),
  fJetRawPtArray(in.fJetRawPtArray),
  fJetMaxTrackPtArray(in.fJetMaxTrackPtArray),
  fTrackPtArray(in.fTrackPtArray),
  fTrackPtErrorArray(in.fTrackPtErrorArray),
  fTrackPhiArray(in.fTrackPhiArray),
  fTrackEtaArray(in.fTrackEtaArray),
  fHighPurityTrackArray(in.fHighPurityTrackArray),
  fTrackVertexDistanceZArray(in.fTrackVertexDistanceZArray),
  fTrackVertexDistanceZErrorArray(in.fTrackVertexDistanceZErrorArray),
  fTrackVertexDistanceXYArray(in.fTrackVertexDistanceXYArray),
  fTrackVertexDistanceXYErrorArray(in.fTrackVertexDistanceXYErrorArray),
  fTrackChi2Array(in.fTrackChi2Array),
  fnTrackDegreesOfFreedomArray(in.fnTrackDegreesOfFreedomArray),
  fnHitsTrackerLayerArray(in.fnHitsTrackerLayerArray),
  fnHitsTrackArray(in.fnHitsTrackArray),
  fTrackEnergyEcalArray(in.fTrackEnergyEcalArray),
  fTrackEnergyHcalArray(in.fTrackEnergyHcalArray)
{
  // Copy constructor
  
}

/*
 * Assignment operator
 */
SkimForestReader& SkimForestReader::operator=(const SkimForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  ForestReader::operator=(in);
  
  fEventTree = in.fEventTree;
  
  fJetPtArray = in.fJetPtArray;
  fJetPhiArray = in.fJetPhiArray;
  fJetEtaArray = in.fJetEtaArray;
  fJetRawPtArray = in.fJetRawPtArray;
  fJetMaxTrackPtArray = in.fJetMaxTrackPtArray;
  fTrackPtArray = in.fTrackPtArray;
  fTrackPtErrorArray = in.fTrackPtErrorArray;
  fTrackPhiArray = in.fTrackPhiArray;
  fTrackEtaArray = in.fTrackEtaArray;
  fHighPurityTrackArray = in.fHighPurityTrackArray;
  fTrackVertexDistanceZArray = in.fTrackVertexDistanceZArray;
  fTrackVertexDistanceZErrorArray = in.fTrackVertexDistanceZErrorArray;
  fTrackVertexDistanceXYArray = in.fTrackVertexDistanceXYArray;
  fTrackVertexDistanceXYErrorArray = in.fTrackVertexDistanceXYErrorArray;
  fTrackChi2Array = in.fTrackChi2Array;
  fnTrackDegreesOfFreedomArray = in.fnTrackDegreesOfFreedomArray;
  fnHitsTrackerLayerArray = in.fnHitsTrackerLayerArray;
  fnHitsTrackArray = in.fnHitsTrackArray;
  fTrackEnergyEcalArray = in.fTrackEnergyEcalArray;
  fTrackEnergyHcalArray = in.fTrackEnergyHcalArray;
  
  return *this;
}

/*
 * Destructor
 */
SkimForestReader::~SkimForestReader(){
  // destructor
}

/*
 * Initialization, meaning that the branches are connected to the tree
 */
void SkimForestReader::Initialize(){
  
  // Connect the branches related to event information
  fEventTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  fEventTree->SetBranchAddress("hiBin",&fHiBin,&fHiBinBranch);
  if(fDataType == kPpMC || fDataType == kPbPbMC){
    fEventTree->SetBranchAddress("pthat",&fPtHat,&fPtHatBranch); // pT hat only for MC
  } else {
    fPtHat = 0; // We do not have pT hat information for real data
  }
  
  // Connect the branches to jet properties
  fEventTree->SetBranchAddress("calo_jtpt",&fJetPtArray,&fJetPtBranch);
  fEventTree->SetBranchAddress("calo_jtphi",&fJetPhiArray,&fJetPhiBranch);
  fEventTree->SetBranchAddress("calo_jteta",&fJetEtaArray,&fJetEtaBranch);
  fEventTree->SetBranchAddress("calo_rawpt",&fJetRawPtArray,&fJetRawPtBranch);
  fEventTree->SetBranchAddress("calo_trackMax",&fJetMaxTrackPtArray,&fJetMaxTrackPtBranch);
  
  // Connect the branches regarding jet trigger requirements
  if(fDataType == kPp || fDataType == kPpMC){ // pp data of MC
    fEventTree->SetBranchAddress("HLT_ak4CaloJet80",&fCaloJetFilterBit,&fCaloJetFilterBranch);
  } else if (fDataType == kPbPb || fDataType == kPbPbMC){ // PbPb data or MC
    fEventTree->SetBranchAddress("HLT_ak4CaloJet100",&fCaloJetFilterBit,&fCaloJetFilterBranch);
  } else { // Local test
    fCaloJetFilterBit = 1;  // No filter for local test
  }
  
  // Connect the branches containing event selection filter bits
  if(fDataType == kPp || fDataType == kPpMC){ // pp data pr MC
    fEventTree->SetBranchAddress("pPAprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fEventTree->SetBranchAddress("pBeamScrapingFilter",&fBeamScrapingFilterBit,&fBeamScrapingBranch);
    fEventTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fCollisionEventSelectionFilterBit = 1;  // No collision event selection filter for pp
    fHfCoincidenceFilterBit = 1; // No HF energy coincidence requirement for pp
    fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pp
  } else if (fDataType == kPbPb || fDataType == kPbPbMC){ // PbPb data or MC
    fEventTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fEventTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fEventTree->SetBranchAddress("pcollisionEventSelection",&fCollisionEventSelectionFilterBit,&fCollisionEventSelectionBranch);
    fEventTree->SetBranchAddress("phfCoincFilter3",&fHfCoincidenceFilterBit,&fHfCoincidenceBranch);
    fEventTree->SetBranchAddress("pclusterCompatibilityFilter",&fClusterCompatibilityFilterBit,&fClusterCompatibilityBranch);
    fBeamScrapingFilterBit = 1;  // No beam scraping filter for PbPb
  } else { // Local test
    fPrimaryVertexFilterBit = 1;
    fBeamScrapingFilterBit = 1;
    fCollisionEventSelectionFilterBit = 1;
    fHBHENoiseFilterBit = 1;
    fHfCoincidenceFilterBit = 1;
    fClusterCompatibilityFilterBit = 1;
  }
  
  // Connect the branches related to tracks
  fEventTree->SetBranchAddress("trkPt",&fTrackPtArray,&fTrackPtBranch);
  fEventTree->SetBranchAddress("trkPtError",&fTrackPtErrorArray,&fTrackPtErrorBranch);
  fEventTree->SetBranchAddress("trkPhi",&fTrackPhiArray,&fTrackPhiBranch);
  fEventTree->SetBranchAddress("trkEta",&fTrackEtaArray,&fTrackEtaBranch);
  fEventTree->SetBranchAddress("highPurity",&fHighPurityTrackArray,&fHighPurityTrackBranch);
  fEventTree->SetBranchAddress("trkDz",&fTrackVertexDistanceZArray,&fTrackVertexDistanceZBranch);
  fEventTree->SetBranchAddress("trkDzError",&fTrackVertexDistanceZErrorArray,&fTrackVertexDistanceZErrorBranch);
  fEventTree->SetBranchAddress("trkDxy",&fTrackVertexDistanceXYArray,&fTrackVertexDistanceXYBranch);
  fEventTree->SetBranchAddress("trkDxyError",&fTrackVertexDistanceXYErrorArray,&fTrackVertexDistanceXYErrorBranch);
  fEventTree->SetBranchAddress("trkChi2",&fTrackChi2Array,&fTrackChi2Branch);
  fEventTree->SetBranchAddress("trkNdof",&fnTrackDegreesOfFreedomArray,&fnTrackDegreesOfFreedomBranch);
  fEventTree->SetBranchAddress("trkNlayer",&fnHitsTrackerLayerArray,&fnHitsTrackerLayerBranch);
  fEventTree->SetBranchAddress("trkNHit",&fnHitsTrackArray,&fnHitsTrackBranch);
  fEventTree->SetBranchAddress("pfEcal",&fTrackEnergyEcalArray,&fTrackEnergyEcalBranch);
  fEventTree->SetBranchAddress("pfHcal",&fTrackEnergyHcalArray,&fTrackEnergyHcalBranch);
  
  // Connect branches related to particle flow candidates
  fEventTree->SetBranchAddress("pfId",&fParticleFlowCandidateIdArray,&fParticleFlowCandidateIdBranch);
  fEventTree->SetBranchAddress("pfPt",&fParticleFlowCandidatePtArray,&fParticleFlowCandidatePtBranch);
  fEventTree->SetBranchAddress("pfPhi",&fParticleFlowCandidatePhiArray,&fParticleFlowCandidatePhiBranch);
  fEventTree->SetBranchAddress("pfEta",&fParticleFlowCandidateEtaArray,&fParticleFlowCandidateEtaBranch);
}

/*
 * Connect a new tree to the reader
 */
void SkimForestReader::ReadForestFromFile(TFile *inputFile){
  
  // Connect a trees from the file to the reader
  fEventTree = (TTree*)inputFile->Get("mixing_tree");
  
  Initialize();
}

/*
 * Burn the current forest.
 */
void SkimForestReader::BurnForest(){
  fEventTree->Delete();
}

/*
 * Load an event to memory
 */
void SkimForestReader::GetEvent(Int_t nEvent){
  
  fEventTree->GetEntry(nEvent);
  
  // Read the numbers of tracks and jets for this event
  fnJets = fJetPtArray->size();
  fnTracks = fTrackPtArray->size();
}

// Getter for jet pT
Float_t SkimForestReader::GetJetPt(Int_t iJet) const{
  return fJetPtArray->at(iJet);
}

// Getter for jet phi
Float_t SkimForestReader::GetJetPhi(Int_t iJet) const{
  return fJetPhiArray->at(iJet);
}

// Getter for jet eta
Float_t SkimForestReader::GetJetEta(Int_t iJet) const{
  return fJetEtaArray->at(iJet);
}

// Getter for jet raw pT
Float_t SkimForestReader::GetJetRawPt(Int_t iJet) const{
  return fJetRawPtArray->at(iJet);
}

// Getter for maximum track pT inside a jet
Float_t SkimForestReader::GetJetMaxTrackPt(Int_t iJet) const{
  return fJetMaxTrackPtArray->at(iJet);
}

// Getter for track pT
Float_t SkimForestReader::GetTrackPt(Int_t iTrack) const{
  return fTrackPtArray->at(iTrack);
}

// Getter for track pT error
Float_t SkimForestReader::GetTrackPtError(Int_t iTrack) const{
  return fTrackPtErrorArray->at(iTrack);
}

// Getter for track phi
Float_t SkimForestReader::GetTrackPhi(Int_t iTrack) const{
  return fTrackPhiArray->at(iTrack);
}

// Getter for track eta
Float_t SkimForestReader::GetTrackEta(Int_t iTrack) const{
  return fTrackEtaArray->at(iTrack);
}

// Getter for high purity of the track
Bool_t SkimForestReader::GetTrackHighPurity(Int_t iTrack) const{
  return fHighPurityTrackArray->at(iTrack);
}

// Getter for track distance from primary vertex in z-direction
Float_t SkimForestReader::GetTrackVertexDistanceZ(Int_t iTrack) const{
  return fTrackVertexDistanceZArray->at(iTrack);
}

// Getter for error of track distance from primary vertex in z-direction
Float_t SkimForestReader::GetTrackVertexDistanceZError(Int_t iTrack) const{
  return fTrackVertexDistanceZErrorArray->at(iTrack);
}

// Getter for track distance from primary vertex in xy-direction
Float_t SkimForestReader::GetTrackVertexDistanceXY(Int_t iTrack) const{
  return fTrackVertexDistanceXYArray->at(iTrack);
}

// Getter for error of track distance from primary vertex in xy-direction
Float_t SkimForestReader::GetTrackVertexDistanceXYError(Int_t iTrack) const{
  return fTrackVertexDistanceXYErrorArray->at(iTrack);
}

// Getter for track chi2 value from reconstruction fit
Float_t SkimForestReader::GetTrackChi2(Int_t iTrack) const{
  return fTrackChi2Array->at(iTrack);
}

// Getter for number of degrees of freedom in reconstruction fit
Int_t SkimForestReader::GetNTrackDegreesOfFreedom(Int_t iTrack) const{
  return fnTrackDegreesOfFreedomArray->at(iTrack);
}

// Getter for number of hits in tracker layers
Int_t SkimForestReader::GetNHitsTrackerLayer(Int_t iTrack) const{
  return fnHitsTrackerLayerArray->at(iTrack);
}

// Getter for number of hits for the track
Int_t SkimForestReader::GetNHitsTrack(Int_t iTrack) const{
  return fnHitsTrackArray->at(iTrack);
}

// Getter for track energy in ECal
Float_t SkimForestReader::GetTrackEnergyEcal(Int_t iTrack) const{
  return fTrackEnergyEcalArray->at(iTrack);
}

// Getter for track energy in HCal
Float_t SkimForestReader::GetTrackEnergyHcal(Int_t iTrack) const{
  return fTrackEnergyHcalArray->at(iTrack);
}

// Getter for track charge. Charge is only relevant for generator level tracks
Int_t SkimForestReader::GetTrackCharge(Int_t iTrack) const{
  return 1;
}

// Getter for track subevent index. Relevant only for generator level tracks.
Int_t SkimForestReader::GetTrackSubevent(Int_t iTrack) const{
  return -1;
}

