// Implementation for ForestReader

// Own includes
#include "ForestReader.h"

/*
 * Default constructor
 */
ForestReader::ForestReader() :
  fDataType(0),
  fReadMode(0),
  fJetType(0),
  fJetAxis(0),
  fMatchJets(false),
  fDoEventPlane(false),
  fHiVzBranch(0),
  fHiBinBranch(0),
  fPtHatBranch(0),
  fnEventPlaneBranch(0),
  fEventPlaneAngleBranch(0),
  fEventPlaneQBranch(0),
  fEventPlaneMultiplicityBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetEtaBranch(0),
  fJetRawPtBranch(0),
  fJetMaxTrackPtBranch(0),
  fJetRefPtBranch(0),
  fJetRefFlavorBranch(0),
  fJetMatchedPtBranch(0),
  fJetMatchedEtaBranch(0),
  fJetMatchedPhiBranch(0),
  fEventWeightBranch(0),
  fCaloJetFilterBranch(0),
  fCaloJetFilterPrescaleBranch(0),
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
  fPtHat(0),
  fnEventPlane(0),
  fEventPlaneAngle(),
  fEventPlaneQ(),
  fEventPlaneMultiplicity(),
  fnJets(0),
  fnMatchedJets(0),
  fEventWeight(1),
  fCaloJetFilterBit(0),
  fCaloJetFilterBitPrescale(0),
  fPrimaryVertexFilterBit(0),
  fBeamScrapingFilterBit(0),
  fCollisionEventSelectionFilterBit(0),
  fHBHENoiseFilterBit(0),
  fHfCoincidenceFilterBit(0),
  fClusterCompatibilityFilterBit(0),
  fnTracks(0),
  fnParticleFlowCandidates(0),
  fParticleFlowCandidateIdArray(),
  fParticleFlowCandidatePtArray(),
  fParticleFlowCandidatePhiArray(),
  fParticleFlowCandidateEtaArray(),
  fParticleFlowCandidateIdVector(0),
  fParticleFlowCandidatePtVector(0),
  fParticleFlowCandidatePhiVector(0),
  fParticleFlowCandidateEtaVector(0)

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
 */
ForestReader::ForestReader(Int_t dataType, Int_t readMode, Int_t jetType, Int_t jetAxis, Bool_t matchJets, Bool_t doEventPlane) :
  fDataType(0),
  fReadMode(readMode),
  fJetType(jetType),
  fJetAxis(jetAxis),
  fMatchJets(matchJets),
  fDoEventPlane(doEventPlane),
  fHiVzBranch(0),
  fHiBinBranch(0),
  fPtHatBranch(0),
  fnEventPlaneBranch(0),
  fEventPlaneAngleBranch(0),
  fEventPlaneQBranch(0),
  fEventPlaneMultiplicityBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetEtaBranch(0),
  fJetRawPtBranch(0),
  fJetMaxTrackPtBranch(0),
  fJetRefPtBranch(0),
  fJetRefFlavorBranch(0),
  fJetMatchedPtBranch(0),
  fJetMatchedEtaBranch(0),
  fJetMatchedPhiBranch(0),
  fEventWeightBranch(0),
  fCaloJetFilterBranch(0),
  fCaloJetFilterPrescaleBranch(0),
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
  fPtHat(0),
  fnEventPlane(0),
  fEventPlaneAngle(),
  fEventPlaneQ(),
  fEventPlaneMultiplicity(),
  fnJets(0),
  fnMatchedJets(0),
  fEventWeight(1),
  fCaloJetFilterBit(0),
  fCaloJetFilterBitPrescale(0),
  fPrimaryVertexFilterBit(0),
  fBeamScrapingFilterBit(0),
  fCollisionEventSelectionFilterBit(0),
  fHBHENoiseFilterBit(0),
  fHfCoincidenceFilterBit(0),
  fClusterCompatibilityFilterBit(0),
  fnTracks(0),
  fnParticleFlowCandidates(0),
  fParticleFlowCandidateIdArray(),
  fParticleFlowCandidatePtArray(),
  fParticleFlowCandidatePhiArray(),
  fParticleFlowCandidateEtaArray(),
  fParticleFlowCandidateIdVector(0),
  fParticleFlowCandidatePtVector(0),
  fParticleFlowCandidatePhiVector(0),
  fParticleFlowCandidateEtaVector(0)
{
  // Custom constructor
  
  SetDataType(dataType);
  
}

/*
 * Copy constructor
 */
ForestReader::ForestReader(const ForestReader& in) :
  fDataType(in.fDataType),
  fReadMode(in.fReadMode),
  fJetType(in.fJetType),
  fJetAxis(in.fJetAxis),
  fMatchJets(in.fMatchJets),
  fDoEventPlane(in.fDoEventPlane),
  fHiVzBranch(in.fHiVzBranch),
  fHiBinBranch(in.fHiBinBranch),
  fPtHatBranch(in.fPtHatBranch),
  fnEventPlaneBranch(in.fnEventPlaneBranch),
  fEventPlaneAngleBranch(in.fEventPlaneAngleBranch),
  fEventPlaneQBranch(in.fEventPlaneQBranch),
  fEventPlaneMultiplicityBranch(in.fEventPlaneMultiplicityBranch),
  fJetPtBranch(in.fJetPtBranch),
  fJetPhiBranch(in.fJetPhiBranch),
  fJetEtaBranch(in.fJetEtaBranch),
  fJetRawPtBranch(in.fJetRawPtBranch),
  fJetMaxTrackPtBranch(in.fJetMaxTrackPtBranch),
  fJetRefPtBranch(in.fJetRefPtBranch),
  fJetRefFlavorBranch(in.fJetRefFlavorBranch),
  fJetMatchedPtBranch(in.fJetMatchedPtBranch),
  fJetMatchedEtaBranch(in.fJetMatchedEtaBranch),
  fJetMatchedPhiBranch(in.fJetMatchedPhiBranch),
  fEventWeightBranch(in.fEventWeightBranch),
  fCaloJetFilterBranch(in.fCaloJetFilterBranch),
  fCaloJetFilterPrescaleBranch(in.fCaloJetFilterPrescaleBranch),
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
  fPtHat(in.fPtHat),
  fnEventPlane(in.fnEventPlane),
  fnJets(in.fnJets),
  fnMatchedJets(in.fnMatchedJets),
  fEventWeight(in.fEventWeight),
  fCaloJetFilterBit(in.fCaloJetFilterBit),
  fCaloJetFilterBitPrescale(in.fCaloJetFilterBitPrescale),
  fPrimaryVertexFilterBit(in.fPrimaryVertexFilterBit),
  fBeamScrapingFilterBit(in.fBeamScrapingFilterBit),
  fCollisionEventSelectionFilterBit(in.fCollisionEventSelectionFilterBit),
  fHBHENoiseFilterBit(in.fHBHENoiseFilterBit),
  fHfCoincidenceFilterBit(in.fHfCoincidenceFilterBit),
  fClusterCompatibilityFilterBit(in.fClusterCompatibilityFilterBit),
  fnTracks(in.fnTracks),
  fnParticleFlowCandidates(in.fnParticleFlowCandidates),
  fParticleFlowCandidateIdVector(in.fParticleFlowCandidateIdVector),
  fParticleFlowCandidatePtVector(in.fParticleFlowCandidatePtVector),
  fParticleFlowCandidatePhiVector(in.fParticleFlowCandidatePhiVector),
  fParticleFlowCandidateEtaVector(in.fParticleFlowCandidateEtaVector)
{

  for(Int_t iParticleFlowCandidate = 0; iParticleFlowCandidate < fnMaxParticleFlowCandidates; iParticleFlowCandidate++){
    fParticleFlowCandidateIdArray[iParticleFlowCandidate] = in.fParticleFlowCandidateIdArray[iParticleFlowCandidate];
    fParticleFlowCandidatePtArray[iParticleFlowCandidate] = in.fParticleFlowCandidatePtArray[iParticleFlowCandidate];
    fParticleFlowCandidatePhiArray[iParticleFlowCandidate] = in.fParticleFlowCandidatePhiArray[iParticleFlowCandidate];
    fParticleFlowCandidateEtaArray[iParticleFlowCandidate] = in.fParticleFlowCandidateEtaArray[iParticleFlowCandidate];
  }
  
  for(Int_t iEventPlane = 0; iEventPlane < fMaxEventPlanes; iEventPlane++){
    fEventPlaneAngle[iEventPlane] = in.fEventPlaneAngle[iEventPlane];
    fEventPlaneQ[iEventPlane] = in.fEventPlaneQ[iEventPlane];
    fEventPlaneMultiplicity[iEventPlane] = in.fEventPlaneMultiplicity[iEventPlane];
  }
  
}

/*
 * Assignment operator
 */
ForestReader& ForestReader::operator=(const ForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  fDataType = in.fDataType;
  fReadMode = in.fReadMode;
  fJetType = in.fJetType;
  fJetAxis = in.fJetAxis;
  fMatchJets = in.fMatchJets;
  fDoEventPlane = in.fDoEventPlane;
  fHiVzBranch = in.fHiVzBranch;
  fHiBinBranch = in.fHiBinBranch;
  fPtHatBranch = in.fPtHatBranch;
  fnEventPlaneBranch = in.fnEventPlaneBranch;
  fEventPlaneAngleBranch = in.fEventPlaneAngleBranch;
  fEventPlaneQBranch = in.fEventPlaneQBranch;
  fEventPlaneMultiplicityBranch = in.fEventPlaneMultiplicityBranch;
  fJetPtBranch = in.fJetPtBranch;
  fJetPhiBranch = in.fJetPhiBranch;
  fJetEtaBranch = in.fJetEtaBranch;
  fJetRawPtBranch = in.fJetRawPtBranch;
  fJetMaxTrackPtBranch = in.fJetMaxTrackPtBranch;
  fJetRefPtBranch = in.fJetRefPtBranch;
  fJetRefFlavorBranch = in.fJetRefFlavorBranch;
  fJetMatchedPtBranch = in.fJetMatchedPtBranch;
  fJetMatchedEtaBranch = in.fJetMatchedEtaBranch;
  fJetMatchedPhiBranch = in.fJetMatchedPhiBranch;
  fEventWeightBranch = in.fEventWeightBranch;
  fCaloJetFilterBranch = in.fCaloJetFilterBranch;
  fCaloJetFilterPrescaleBranch = in.fCaloJetFilterPrescaleBranch;
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
  fPtHat = in.fPtHat;
  fnEventPlane = in.fnEventPlane;
  fnJets = in.fnJets;
  fnMatchedJets = in.fnMatchedJets;
  fEventWeight = in.fEventWeight;
  fCaloJetFilterBit = in.fCaloJetFilterBit;
  fCaloJetFilterBitPrescale = in.fCaloJetFilterBitPrescale;
  fPrimaryVertexFilterBit = in.fPrimaryVertexFilterBit;
  fBeamScrapingFilterBit = in.fBeamScrapingFilterBit;
  fCollisionEventSelectionFilterBit = in.fCollisionEventSelectionFilterBit;
  fHBHENoiseFilterBit = in.fHBHENoiseFilterBit;
  fHfCoincidenceFilterBit = in.fHfCoincidenceFilterBit;
  fClusterCompatibilityFilterBit = in.fClusterCompatibilityFilterBit;
  fnTracks = in.fnTracks;
  fnParticleFlowCandidates = in.fnParticleFlowCandidates;
  fParticleFlowCandidateIdVector = in.fParticleFlowCandidateIdVector;
  fParticleFlowCandidatePtVector = in.fParticleFlowCandidatePtVector;
  fParticleFlowCandidatePhiVector = in.fParticleFlowCandidatePhiVector;
  fParticleFlowCandidateEtaVector = in.fParticleFlowCandidateEtaVector;
  
  for(Int_t iParticleFlowCandidate = 0; iParticleFlowCandidate < fnMaxParticleFlowCandidates; iParticleFlowCandidate++){
    fParticleFlowCandidateIdArray[iParticleFlowCandidate] = in.fParticleFlowCandidateIdArray[iParticleFlowCandidate];
    fParticleFlowCandidatePtArray[iParticleFlowCandidate] = in.fParticleFlowCandidatePtArray[iParticleFlowCandidate];
    fParticleFlowCandidatePhiArray[iParticleFlowCandidate] = in.fParticleFlowCandidatePhiArray[iParticleFlowCandidate];
    fParticleFlowCandidateEtaArray[iParticleFlowCandidate] = in.fParticleFlowCandidateEtaArray[iParticleFlowCandidate];
  }
  
  for(Int_t iEventPlane = 0; iEventPlane < fMaxEventPlanes; iEventPlane++){
    fEventPlaneAngle[iEventPlane] = in.fEventPlaneAngle[iEventPlane];
    fEventPlaneQ[iEventPlane] = in.fEventPlaneQ[iEventPlane];
    fEventPlaneMultiplicity[iEventPlane] = in.fEventPlaneMultiplicity[iEventPlane];
  }
  
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

// Getter for pT hat
Float_t ForestReader::GetPtHat() const{
  return fPtHat;
}

// Getter for pT hat
Float_t ForestReader::GetEventWeight() const{
  return fEventWeight;
}

// Getter for the number of event planes
Int_t ForestReader::GenNEventPlane() const{
  return fnEventPlane;
}

// Getter for the event plane angle for the i:th event plane
Float_t ForestReader::GetEventPlaneAngle(Int_t iEventPlane) const{
  return fEventPlaneAngle[iEventPlane];
}

// Getter for the magnitude of the q-vector for the i:th event plane
Float_t ForestReader::GetEventPlaneQ(Int_t iEventPlane) const{
  return fEventPlaneQ[iEventPlane];
}

// Getter for the particle multiplicity in the i:th event plane
Float_t ForestReader::GetEventPlaneMultiplicity(Int_t iEventPlane) const{
  return fEventPlaneMultiplicity[iEventPlane];
}

// Getter for calorimeter jet filter bit. Always 1 for MC (set in the initializer).
Int_t ForestReader::GetCaloJetFilterBit() const{
  return fCaloJetFilterBit;
}

// Getter for prescaled calorimeter jet filter bit. Always 1 for MC (set in the initializer).
Int_t ForestReader::GetPrescaledCaloJetFilterBit() const{
  return fCaloJetFilterBitPrescale;
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
  if(fReadMode > 2000) return 0; // No PF candidates for 2018 data
  if(fReadMode == 1 && fDataType == kPpMC) return fParticleFlowCandidateIdArray[iCandidate];   // For PYTHIA8 forest
  return fParticleFlowCandidateIdVector->at(iCandidate);                                       // Regular
}

// Getter for particle flow candidate pT
Float_t ForestReader::GetParticleFlowCandidatePt(Int_t iCandidate) const{
  if(fReadMode > 2000) return 0; // No PF candidates for 2018 data
  if(fReadMode == 1 && fDataType == kPpMC) return fParticleFlowCandidatePtArray[iCandidate];  // For PYTHIA8 forest
  return fParticleFlowCandidatePtVector->at(iCandidate);                                      // Regular
}

// Getter for particle flow candidate phi
Float_t ForestReader::GetParticleFlowCandidatePhi(Int_t iCandidate) const{
  if(fReadMode > 2000) return 0; // No PF candidates for 2018 data
  if(fReadMode == 1 && fDataType == kPpMC) return fParticleFlowCandidatePhiArray[iCandidate];  // For PYTHIA8 forest
  return fParticleFlowCandidatePhiVector->at(iCandidate);                                      // Regular
}

// Getter for particle flow candidate eta
Float_t ForestReader::GetParticleFlowCandidateEta(Int_t iCandidate) const{
  if(fReadMode > 2000) return 0; // No PF candidates for 2018 data
  if(fReadMode == 1 && fDataType == kPpMC)  return fParticleFlowCandidateEtaArray[iCandidate];  // For PYTHIA8 forest
  return fParticleFlowCandidateEtaVector->at(iCandidate);                                       // Regular
}

// Getter number of particle flow candidates in an event
Int_t ForestReader::GetNParticleFlowCandidates() const{
  if(fReadMode > 2000) return 0; // No PF candidates for 2018 data
  if(fReadMode == 1 && fDataType == kPpMC) return fnParticleFlowCandidates; // For PYTHIA8 forest
  return fParticleFlowCandidateIdVector->size();                            // Regular
}

// Getter for track algorithm
Int_t ForestReader::GetTrackAlgorithm(Int_t iTrack) const{
  return 0;
}

// Getter for the original track algorithm
Int_t ForestReader::GetTrackOriginalAlgorithm(Int_t iTrack) const{
  return 0;
}

// Getter for track MVA
Float_t ForestReader::GetTrackMVA(Int_t iTrack) const{
  return 0;
}
