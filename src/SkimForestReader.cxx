// Implementation for SkimForestReader

// Own includes
#include "SkimForestReader.h"

/*
 * Default constructor
 */
SkimForestReader::SkimForestReader() :
  ForestReader(),
  fEventTree(0),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fJetRawPtArray(),
  fJetMaxTrackPtArray(),
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
  fTrackEnergyHcalArray()
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
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
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
  fTrackEnergyHcalArray()
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
  
  // Connect the branches to jet properties
  fEventTree->SetBranchAddress("calo_jtpt",&fJetPtArray,&fJetPtBranch);
  fEventTree->SetBranchAddress("calo_jtphi",&fJetPhiArray,&fJetPhiBranch);
  fEventTree->SetBranchAddress("calo_jteta",&fJetEtaArray,&fJetEtaBranch);
  fEventTree->SetBranchAddress("calo_rawpt",&fJetRawPtArray,&fJetRawPtBranch);
  fEventTree->SetBranchAddress("calo_trackMax",&fJetMaxTrackPtArray,&fJetMaxTrackPtBranch);
  
  // Connect the branches regarding jet trigger requirements
  if(fDataType == kPp){ // pp data
    fEventTree->SetBranchAddress("HLT_ak4CaloJet80",&fCaloJetFilterBit,&fCaloJetFilterBranch);
  } else if (fDataType == kPbPb){ // PbPb data
    fEventTree->SetBranchAddress("HLT_ak4CaloJet100",&fCaloJetFilterBit,&fCaloJetFilterBranch);
  } else { // Monte Carlo
    fCaloJetFilterBit = 1;  // No filter for Monte Carlo
  }
  
  // Connect the branches containing event selection filter bits
  if(fDataType == kPp){ // pp data
    fEventTree->SetBranchAddress("pPAprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fEventTree->SetBranchAddress("pBeamScrapingFilter",&fBeamScrapingFilterBit,&fBeamScrapingBranch);
    fEventTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fCollisionEventSelectionFilterBit = 1;  // No collision event selection filter for pp
    fHfCoincidenceFilterBit = 1; // No HF energy coincidence requirement for pp
    fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pp
  } else if (fDataType == kPbPb){ // PbPb data
    fEventTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fEventTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fEventTree->SetBranchAddress("pcollisionEventSelection",&fCollisionEventSelectionFilterBit,&fCollisionEventSelectionBranch);
    fEventTree->SetBranchAddress("phfCoincFilter3",&fHfCoincidenceFilterBit,&fHfCoincidenceBranch);
    fEventTree->SetBranchAddress("pclusterCompatibilityFilter",&fClusterCompatibilityFilterBit,&fClusterCompatibilityBranch);
    fBeamScrapingFilterBit = 1;  // No beam scraping filter for PbPb
  } else { // Monte Carlo
    fPrimaryVertexFilterBit = 1;            // No filter for Monte Carlo
    fBeamScrapingFilterBit = 1;             // No filter for Monte Carlo
    fCollisionEventSelectionFilterBit = 1;  // No filter for Monte Carlo
    fHBHENoiseFilterBit = 1;                // No filter for Monte Carlo
    fHfCoincidenceFilterBit = 1;            // No HF energy coincidence requirement for Monte Carlo
    fClusterCompatibilityFilterBit = 1;     // No cluster compatibility requirement for Monte Carlo
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
