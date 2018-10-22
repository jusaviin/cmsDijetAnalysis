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
 *   Int_t readMode: 0 = Regular forests, 1 = Official PYTHIA8 forest
 *   Int_t jetType: 0 = Calo jets, 1 = PF jets
 */
SkimForestReader::SkimForestReader(Int_t dataType, Int_t readMode, Int_t jetType) :
  ForestReader(dataType,readMode,jetType),
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
  
  // Only enable the branches that are actually read
  fEventTree->SetBranchStatus("*",0);
  
  // Connect the branches related to event information
  fEventTree->SetBranchStatus("vz",1);
  fEventTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  if(fDataType == kPp || fDataType == kPpMC){
    fHiBin = -1;  // The skims for pp do not have hiBin
  } else {
    fEventTree->SetBranchStatus("hiBin",1);
    fEventTree->SetBranchAddress("hiBin",&fHiBin,&fHiBinBranch);
  }
  if(fDataType == kPpMC || fDataType == kPbPbMC){
    fEventTree->SetBranchStatus("pthat",1);
    fEventTree->SetBranchAddress("pthat",&fPtHat,&fPtHatBranch); // pT hat only for MC
  } else {
    fPtHat = 0; // We do not have pT hat information for real data
  }
  
  // Connect the branches to jet properties
  const char * jetType[2] = {"calo","pf"};
  char branchName[20];
  sprintf(branchName,"%s_jtpt",jetType[fJetType]);
  fEventTree->SetBranchStatus(branchName,1);
  fEventTree->SetBranchAddress(branchName,&fJetPtArray,&fJetPtBranch);
  sprintf(branchName,"%s_jtphi",jetType[fJetType]);
  fEventTree->SetBranchStatus(branchName,1);
  fEventTree->SetBranchAddress(branchName,&fJetPhiArray,&fJetPhiBranch);
  sprintf(branchName,"%s_jteta",jetType[fJetType]);
  fEventTree->SetBranchStatus(branchName,1);
  fEventTree->SetBranchAddress(branchName,&fJetEtaArray,&fJetEtaBranch);
  sprintf(branchName,"%s_rawpt",jetType[fJetType]);
  fEventTree->SetBranchStatus(branchName,1);
  fEventTree->SetBranchAddress(branchName,&fJetRawPtArray,&fJetRawPtBranch);
  sprintf(branchName,"%s_trackMax",jetType[fJetType]);
  fEventTree->SetBranchStatus(branchName,1);
  fEventTree->SetBranchAddress(branchName,&fJetMaxTrackPtArray,&fJetMaxTrackPtBranch);
  
  // Connect the branches to the HLT tree
  if(fDataType == kPp){ // pp data
    fEventTree->SetBranchStatus("HLT_AK4CaloJet80_Eta5p1_v1",1);
    fEventTree->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1",&fCaloJetFilterBit,&fCaloJetFilterBranch);
  } else if (fDataType == kPbPb){ // PbPb
    fEventTree->SetBranchStatus("HLT_ak4CaloJet100_Prescale",1);
    fEventTree->SetBranchAddress("HLT_ak4CaloJet100_Prescale",&fCaloJetFilterBit,&fCaloJetFilterBranch);
  } else { // Local test or MC
    fCaloJetFilterBit = 1;  // No filter for local test or MC
  }
  
  // Connect the branches containing event selection filter bits
  if(fDataType == kPp){ // pp data pr MC
    fEventTree->SetBranchStatus("pPAprimaryVertexFilter",1);
    fEventTree->SetBranchAddress("pPAprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fEventTree->SetBranchStatus("pBeamScrapingFilter",1);
    fEventTree->SetBranchAddress("pBeamScrapingFilter",&fBeamScrapingFilterBit,&fBeamScrapingBranch);
    fEventTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
    fEventTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fCollisionEventSelectionFilterBit = 1;  // No collision event selection filter for pp
    fHfCoincidenceFilterBit = 1; // No HF energy coincidence requirement for pp
    fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pp
  } else if (fDataType == kPpMC){
    fEventTree->SetBranchStatus("pprimaryVertexFilter",1);
    fEventTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch); // Naming in Dhanush's skim
    //fEventTree->SetBranchStatus("pPAprimaryVertexFilter",1);
    //fEventTree->SetBranchAddress("pPAprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch); // Naming in Kurt's skim
    fBeamScrapingFilterBit = 1; // No beam scraping filter for MC
    fEventTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
    fEventTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fCollisionEventSelectionFilterBit = 1;  // No collision event selection filter for pp MC
    fHfCoincidenceFilterBit = 1; // No HF energy coincidence requirement for pp MC
    fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pp MC
  } else if (fDataType == kPbPb || fDataType == kPbPbMC){ // PbPb data or MC
    fEventTree->SetBranchStatus("pprimaryVertexFilter",1);
    fEventTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fEventTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
    fEventTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fEventTree->SetBranchStatus("pcollisionEventSelection",1);
    fEventTree->SetBranchAddress("pcollisionEventSelection",&fCollisionEventSelectionFilterBit,&fCollisionEventSelectionBranch);
    fEventTree->SetBranchStatus("phfCoincFilter3",1);
    fEventTree->SetBranchAddress("phfCoincFilter3",&fHfCoincidenceFilterBit,&fHfCoincidenceBranch);
    fEventTree->SetBranchStatus("pclusterCompatibilityFilter",1);
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
  fEventTree->SetBranchStatus("trkPt",1);
  fEventTree->SetBranchAddress("trkPt",&fTrackPtArray,&fTrackPtBranch);
  fEventTree->SetBranchStatus("trkPtError",1);
  fEventTree->SetBranchAddress("trkPtError",&fTrackPtErrorArray,&fTrackPtErrorBranch);
  fEventTree->SetBranchStatus("trkPhi",1);
  fEventTree->SetBranchAddress("trkPhi",&fTrackPhiArray,&fTrackPhiBranch);
  fEventTree->SetBranchStatus("trkEta",1);
  fEventTree->SetBranchAddress("trkEta",&fTrackEtaArray,&fTrackEtaBranch);
  fEventTree->SetBranchStatus("highPurity",1);
  fEventTree->SetBranchAddress("highPurity",&fHighPurityTrackArray,&fHighPurityTrackBranch);
  fEventTree->SetBranchStatus("trkDz",1);
  fEventTree->SetBranchAddress("trkDz",&fTrackVertexDistanceZArray,&fTrackVertexDistanceZBranch);
  fEventTree->SetBranchStatus("trkDzError",1);
  fEventTree->SetBranchAddress("trkDzError",&fTrackVertexDistanceZErrorArray,&fTrackVertexDistanceZErrorBranch);
  fEventTree->SetBranchStatus("trkDxy",1);
  fEventTree->SetBranchAddress("trkDxy",&fTrackVertexDistanceXYArray,&fTrackVertexDistanceXYBranch);
  fEventTree->SetBranchStatus("trkDxyError",1);
  fEventTree->SetBranchAddress("trkDxyError",&fTrackVertexDistanceXYErrorArray,&fTrackVertexDistanceXYErrorBranch);
  fEventTree->SetBranchStatus("trkChi2",1);
  fEventTree->SetBranchAddress("trkChi2",&fTrackChi2Array,&fTrackChi2Branch);
  fEventTree->SetBranchStatus("trkNdof",1);
  fEventTree->SetBranchAddress("trkNdof",&fnTrackDegreesOfFreedomArray,&fnTrackDegreesOfFreedomBranch);
  fEventTree->SetBranchStatus("trkNlayer",1);
  fEventTree->SetBranchAddress("trkNlayer",&fnHitsTrackerLayerArray,&fnHitsTrackerLayerBranch);
  fEventTree->SetBranchStatus("trkNHit",1);
  fEventTree->SetBranchAddress("trkNHit",&fnHitsTrackArray,&fnHitsTrackBranch);
  fEventTree->SetBranchStatus("pfEcal",1);
  fEventTree->SetBranchAddress("pfEcal",&fTrackEnergyEcalArray,&fTrackEnergyEcalBranch);
  fEventTree->SetBranchStatus("pfHcal",1);
  fEventTree->SetBranchAddress("pfHcal",&fTrackEnergyHcalArray,&fTrackEnergyHcalBranch);
  
  // Connect branches related to particle flow candidates
  fEventTree->SetBranchStatus("pfId",1);
  fEventTree->SetBranchAddress("pfId",&fParticleFlowCandidateIdVector,&fParticleFlowCandidateIdBranch);
  fEventTree->SetBranchStatus("pfPt",1);
  fEventTree->SetBranchAddress("pfPt",&fParticleFlowCandidatePtVector,&fParticleFlowCandidatePtBranch);
  fEventTree->SetBranchStatus("pfPhi",1);
  fEventTree->SetBranchAddress("pfPhi",&fParticleFlowCandidatePhiVector,&fParticleFlowCandidatePhiBranch);
  fEventTree->SetBranchStatus("pfEta",1);
  fEventTree->SetBranchAddress("pfEta",&fParticleFlowCandidateEtaVector,&fParticleFlowCandidateEtaBranch);
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

// Getter for track MC status. Relevant only for generator level tracks.
Int_t SkimForestReader::GetTrackMCStatus(Int_t iTrack) const{
  return 1;
}
