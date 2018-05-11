// Implementation for ForestReader

// Own includes
#include "ForestReader.h"

/*
 * Default constructor
 */
ForestReader::ForestReader() :
  fDataType(0),
  fHiVzBranch(0),
  fHiBinBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetEtaBranch(0),
  fJetRawPtBranch(0),
  fJetMaxTrackPtBranch(0),
  fCaloJetFilterBranch(0),
  fPrimaryVertexBranch(0),
  fBeamScrapingBranch(0),
  fCollisionEventSelectionBranch(0),
  fHBHENoiseBranch(0),
  fHfCoincidenceBranch(0),
  fClusterCompatibilityBranch(0),
  fTrackPtBranch(0),
  fTrackPtErrorBranch(0),
  fTrackPhiBranch(0),
  fTrackEtaBranch(0),
  fHighPurityTrackBranch(0),
  fTrackVertexDistanceZBranch(0),
  fTrackVertexDistanceZErrorBranch(0),
  fTrackVertexDistanceXYBranch(0),
  fTrackVertexDistanceXYErrorBranch(0),
  fTrackChi2Branch(0),
  fnTrackDegreesOfFreedomBranch(0),
  fnHitsTrackerLayerBranch(0),
  fnHitsTrackBranch(0),
  fTrackEnergyEcalBranch(0),
  fTrackEnergyHcalBranch(0),
  fParticleFlowCandidateIdBranch(0),
  fParticleFlowCandidatePtBranch(0),
  fParticleFlowCandidatePhiBranch(0),
  fParticleFlowCandidateEtaBranch(0),
  fVertexZ(-100),
  fHiBin(-1),
  fnJets(0),
  fCaloJetFilterBit(0),
  fPrimaryVertexFilterBit(0),
  fBeamScrapingFilterBit(0),
  fCollisionEventSelectionFilterBit(0),
  fHBHENoiseFilterBit(0),
  fHfCoincidenceFilterBit(0),
  fClusterCompatibilityFilterBit(0),
  fnTracks(0)
{
  // Default constructor
}

/*
 * Custom constructor
 *
 *  Arguments:
 *   Int_t dataType: 0 = pp, 1 = PbPb, 2 = pp MC, 3 = PbPb MC, 4 = Local Test
 */
ForestReader::ForestReader(Int_t dataType) :
  fDataType(0),
  fHiVzBranch(0),
  fHiBinBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetEtaBranch(0),
  fJetRawPtBranch(0),
  fJetMaxTrackPtBranch(0),
  fCaloJetFilterBranch(0),
  fPrimaryVertexBranch(0),
  fBeamScrapingBranch(0),
  fCollisionEventSelectionBranch(0),
  fHBHENoiseBranch(0),
  fHfCoincidenceBranch(0),
  fClusterCompatibilityBranch(0),
  fTrackPtBranch(0),
  fTrackPtErrorBranch(0),
  fTrackPhiBranch(0),
  fTrackEtaBranch(0),
  fHighPurityTrackBranch(0),
  fTrackVertexDistanceZBranch(0),
  fTrackVertexDistanceZErrorBranch(0),
  fTrackVertexDistanceXYBranch(0),
  fTrackVertexDistanceXYErrorBranch(0),
  fTrackChi2Branch(0),
  fnTrackDegreesOfFreedomBranch(0),
  fnHitsTrackerLayerBranch(0),
  fnHitsTrackBranch(0),
  fTrackEnergyEcalBranch(0),
  fTrackEnergyHcalBranch(0),
  fParticleFlowCandidateIdBranch(0),
  fParticleFlowCandidatePtBranch(0),
  fParticleFlowCandidatePhiBranch(0),
  fParticleFlowCandidateEtaBranch(0),
  fVertexZ(-100),
  fHiBin(-1),
  fnJets(0),
  fCaloJetFilterBit(0),
  fPrimaryVertexFilterBit(0),
  fBeamScrapingFilterBit(0),
  fCollisionEventSelectionFilterBit(0),
  fHBHENoiseFilterBit(0),
  fHfCoincidenceFilterBit(0),
  fClusterCompatibilityFilterBit(0),
  fnTracks(0)
{
  // Custom constructor
  
  SetDataType(dataType);
  
}

/*
 * Copy constructor
 */
ForestReader::ForestReader(const ForestReader& in) :
  fDataType(in.fDataType),
  fHiVzBranch(in.fHiVzBranch),
  fHiBinBranch(in.fHiBinBranch),
  fJetPtBranch(in.fJetPtBranch),
  fJetPhiBranch(in.fJetPhiBranch),
  fJetEtaBranch(in.fJetEtaBranch),
  fJetRawPtBranch(in.fJetRawPtBranch),
  fJetMaxTrackPtBranch(in.fJetMaxTrackPtBranch),
  fCaloJetFilterBranch(in.fCaloJetFilterBranch),
  fPrimaryVertexBranch(in.fPrimaryVertexBranch),
  fBeamScrapingBranch(in.fBeamScrapingBranch),
  fCollisionEventSelectionBranch(in.fCollisionEventSelectionBranch),
  fHBHENoiseBranch(in.fHBHENoiseBranch),
  fHfCoincidenceBranch(in.fHfCoincidenceBranch),
  fClusterCompatibilityBranch(in.fClusterCompatibilityBranch),
  fTrackPtBranch(in.fTrackPtBranch),
  fTrackPtErrorBranch(in.fTrackPtErrorBranch),
  fTrackPhiBranch(in.fTrackPhiBranch),
  fTrackEtaBranch(in.fTrackEtaBranch),
  fHighPurityTrackBranch(in.fHighPurityTrackBranch),
  fTrackVertexDistanceZBranch(in.fTrackVertexDistanceZBranch),
  fTrackVertexDistanceZErrorBranch(in.fTrackVertexDistanceZErrorBranch),
  fTrackVertexDistanceXYBranch(in.fTrackVertexDistanceXYBranch),
  fTrackVertexDistanceXYErrorBranch(in.fTrackVertexDistanceXYErrorBranch),
  fTrackChi2Branch(in.fTrackChi2Branch),
  fnTrackDegreesOfFreedomBranch(in.fnTrackDegreesOfFreedomBranch),
  fnHitsTrackerLayerBranch(in.fnHitsTrackerLayerBranch),
  fnHitsTrackBranch(in.fnHitsTrackBranch),
  fTrackEnergyEcalBranch(in.fTrackEnergyEcalBranch),
  fTrackEnergyHcalBranch(in.fTrackEnergyHcalBranch),
  fParticleFlowCandidateIdBranch(in.fParticleFlowCandidateIdBranch),
  fParticleFlowCandidatePtBranch(in.fParticleFlowCandidatePtBranch),
  fParticleFlowCandidatePhiBranch(in.fParticleFlowCandidatePhiBranch),
  fParticleFlowCandidateEtaBranch(in.fParticleFlowCandidateEtaBranch),
  fVertexZ(in.fVertexZ),
  fHiBin(in.fHiBin),
  fnJets(in.fnJets),
  fCaloJetFilterBit(in.fCaloJetFilterBit),
  fPrimaryVertexFilterBit(in.fPrimaryVertexFilterBit),
  fBeamScrapingFilterBit(in.fBeamScrapingFilterBit),
  fCollisionEventSelectionFilterBit(in.fCollisionEventSelectionFilterBit),
  fHBHENoiseFilterBit(in.fHBHENoiseFilterBit),
  fHfCoincidenceFilterBit(in.fHfCoincidenceFilterBit),
  fClusterCompatibilityFilterBit(in.fClusterCompatibilityFilterBit),
  fnTracks(in.fnTracks)
{

}

/*
 * Assignment operator
 */
ForestReader& ForestReader::operator=(const ForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  fDataType = in.fDataType;
  fHiVzBranch = in.fHiVzBranch;
  fHiBinBranch = in.fHiBinBranch;
  fJetPtBranch = in.fJetPtBranch;
  fJetPhiBranch = in.fJetPhiBranch;
  fJetEtaBranch = in.fJetEtaBranch;
  fJetRawPtBranch = in.fJetRawPtBranch;
  fJetMaxTrackPtBranch = in.fJetMaxTrackPtBranch;
  fCaloJetFilterBranch = in.fCaloJetFilterBranch;
  fPrimaryVertexBranch = in.fPrimaryVertexBranch;
  fBeamScrapingBranch = in.fBeamScrapingBranch;
  fCollisionEventSelectionBranch = in.fCollisionEventSelectionBranch;
  fHBHENoiseBranch = in.fHBHENoiseBranch;
  fHfCoincidenceBranch = in.fHfCoincidenceBranch;
  fClusterCompatibilityBranch = in.fClusterCompatibilityBranch;
  fTrackPtBranch = in.fTrackPtBranch;
  fTrackPtErrorBranch = in.fTrackPtErrorBranch;
  fTrackPhiBranch = in.fTrackPhiBranch;
  fTrackEtaBranch = in.fTrackEtaBranch;
  fHighPurityTrackBranch = in.fHighPurityTrackBranch;
  fTrackVertexDistanceZBranch = in.fTrackVertexDistanceZBranch;
  fTrackVertexDistanceZErrorBranch = in.fTrackVertexDistanceZErrorBranch;
  fTrackVertexDistanceXYBranch = in.fTrackVertexDistanceXYBranch;
  fTrackVertexDistanceXYErrorBranch = in.fTrackVertexDistanceXYErrorBranch;
  fTrackChi2Branch = in.fTrackChi2Branch;
  fnTrackDegreesOfFreedomBranch = in.fnTrackDegreesOfFreedomBranch;
  fnHitsTrackerLayerBranch = in.fnHitsTrackerLayerBranch;
  fnHitsTrackBranch = in.fnHitsTrackBranch;
  fTrackEnergyEcalBranch = in.fTrackEnergyEcalBranch;
  fTrackEnergyHcalBranch = in.fTrackEnergyHcalBranch;
  fParticleFlowCandidateIdBranch = in.fParticleFlowCandidateIdBranch;
  fParticleFlowCandidatePtBranch = in.fParticleFlowCandidatePtBranch;
  fParticleFlowCandidatePhiBranch = in.fParticleFlowCandidatePhiBranch;
  fParticleFlowCandidateEtaBranch = in.fParticleFlowCandidateEtaBranch;
  fVertexZ = in.fVertexZ;
  fHiBin = in.fHiBin;
  fnJets = in.fnJets;
  fCaloJetFilterBit = in.fCaloJetFilterBit;
  fPrimaryVertexFilterBit = in.fPrimaryVertexFilterBit;
  fBeamScrapingFilterBit = in.fBeamScrapingFilterBit;
  fCollisionEventSelectionFilterBit = in.fCollisionEventSelectionFilterBit;
  fHBHENoiseFilterBit = in.fHBHENoiseFilterBit;
  fHfCoincidenceFilterBit = in.fHfCoincidenceFilterBit;
  fClusterCompatibilityFilterBit = in.fClusterCompatibilityFilterBit;
  fnTracks = in.fnTracks;
  
  return *this;
}

/*
 * Destructor
 */
ForestReader::~ForestReader(){
  // destructor
}

/*
 * Setter for fDataType
 */
void ForestReader::SetDataType(Int_t dataType){
  
  //Sanity check for given data type
  if(dataType < 0 || dataType > knDataTypes-1){
    cout << "ERROR: Data type input " << dataType << " is invalid in ForestReader.cxx!" << endl;
    cout << "Please give integer between 0 and " << knDataTypes-1 << "." << endl;
    cout << "Setting data type to 0 (pp)." << endl;
    fDataType = 0;
  } else {
    
    // If the sanity check passes, set the given data type
    fDataType = dataType;
  }
}

// Getter for number of events in the tree
Int_t ForestReader::GetNEvents() const{
  return fJetPtBranch->GetEntries();
}

// Getter for number of jets in an event
Int_t ForestReader::GetNJets() const{
  return fnJets;
}

// Getter for vertex z position
Float_t ForestReader::GetVz() const{
  return fVertexZ;
}

// Getter for centrality. CMS has integer centrality bins from 0 to 200, thus division by 2.
Float_t ForestReader::GetCentrality() const{
  return fHiBin/2.0;
}

// Getter for hiBin. Return 1 for negative values (for easier handling of tracking efficiency correction)
Int_t ForestReader::GetHiBin() const{
  if(fHiBin < 0) return 1;
  return fHiBin;
}

// Getter for calorimeter jet filter bit. Always 1 for MC (set in the initializer).
Int_t ForestReader::GetCaloJetFilterBit() const{
  return fCaloJetFilterBit;
}

// Getter for primary vertex filter bit. Always 1 for MC (set in the initializer).
Int_t ForestReader::GetPrimaryVertexFilterBit() const{
  return fPrimaryVertexFilterBit;
}

// Getter for beam scraping filter bit. Always 1 for MC and PbPb (set in the initializer).
Int_t ForestReader::GetBeamScrapingFilterBit() const{
  return fBeamScrapingFilterBit;
}

// Getter for HB/HE noise filter bit. Always 1 for MC (set in the initializer).
Int_t ForestReader::GetHBHENoiseFilterBit() const{
  return fHBHENoiseFilterBit;
}

// Getter for collision event selection filter bit. Always 1 for MC and pp (set in the initializer).
Int_t ForestReader::GetCollisionEventSelectionFilterBit() const{
  return fCollisionEventSelectionFilterBit;
}

// Getter for HF energy coincidence filter bit. Always 1 for MC and pp (set in the initializer).
Int_t ForestReader::GetHfCoincidenceFilterBit() const{
  return fHfCoincidenceFilterBit;
}

// Getter for cluster compatibility filter bit. Always 1 for MC and pp (set in the initializer).
Int_t ForestReader::GetClusterCompatibilityFilterBit() const{
  return fClusterCompatibilityFilterBit;
}

// Getter for number of tracks in an event
Int_t ForestReader::GetNTracks() const{
  return fnTracks;
}

// Getter for particle flow candidate ID
Int_t ForestReader::GetParticleFlowCandidateId(Int_t iCandidate) const{
  return fParticleFlowCandidateIdArray->at(iCandidate);
}

// Getter for particle flow candidate pT
Float_t ForestReader::GetParticleFlowCandidatePt(Int_t iCandidate) const{
  return fParticleFlowCandidatePtArray->at(iCandidate);
}

// Getter for particle flow candidate phi
Float_t ForestReader::GetParticleFlowCandidatePhi(Int_t iCandidate) const{
  return fParticleFlowCandidatePhiArray->at(iCandidate);
}

// Getter for particle flow candidate eta
Float_t ForestReader::GetParticleFlowCandidateEta(Int_t iCandidate) const{
  return fParticleFlowCandidateEtaArray->at(iCandidate);
}

// Getter number of particle flow candidates in an event
Int_t ForestReader::GetNParticleFlowCandidates() const{
  return fParticleFlowCandidateIdArray->size();
}
