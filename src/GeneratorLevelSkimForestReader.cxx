// Implementation for GeneratorLevelSkimForestReader

// Own includes
#include "GeneratorLevelSkimForestReader.h"

/*
 * Default constructor
 */
GeneratorLevelSkimForestReader::GeneratorLevelSkimForestReader() :
  SkimForestReader(),
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
 */
GeneratorLevelSkimForestReader::GeneratorLevelSkimForestReader(Int_t dataType) :
  SkimForestReader(dataType),
  fTrackChargeArray(0),
  fTrackSubeventArray(0)
{
  // Custom constructor
  
}

/*
 * Copy constructor
 */
GeneratorLevelSkimForestReader::GeneratorLevelSkimForestReader(const GeneratorLevelSkimForestReader& in) :
  SkimForestReader(in),
  fTrackChargeArray(in.fTrackChargeArray),
  fTrackSubeventArray(in.fTrackSubeventArray)
{
  // Copy constructor
  
}

/*
 * Assignment operator
 */
GeneratorLevelSkimForestReader& GeneratorLevelSkimForestReader::operator=(const GeneratorLevelSkimForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  SkimForestReader::operator=(in);
  
  fTrackChargeArray = in.fTrackChargeArray;
  fTrackSubeventArray = in.fTrackSubeventArray;
  
  return *this;
}

/*
 * Destructor
 */
GeneratorLevelSkimForestReader::~GeneratorLevelSkimForestReader(){
  // destructor
}

/*
 * Initialization, meaning that the branches are connected to the tree
 */
void GeneratorLevelSkimForestReader::Initialize(){
  
  // Connect the branches related to event information
  fEventTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  if(fDataType == kPpMC){
    fHiBin = -1;  // The skims for pp do not have hiBin
  } else {
    fEventTree->SetBranchAddress("hiBin",&fHiBin,&fHiBinBranch);
  }
  fEventTree->SetBranchAddress("pthat",&fPtHat,&fPtHatBranch); // pT hat only for MC
  
  // Connect the branches to jet properties
  fEventTree->SetBranchAddress("genpt",&fJetPtArray,&fJetPtBranch);
  fEventTree->SetBranchAddress("genphi",&fJetPhiArray,&fJetPhiBranch);
  fEventTree->SetBranchAddress("geneta",&fJetEtaArray,&fJetEtaBranch);
  
  // Connect the branches to the HLT tree
  fCaloJetFilterBit = 1;  // No calorimeter filter bit is present in the skims
  
  // Connect the branches containing event selection filter bits
  if (fDataType == kPpMC){
    //fEventTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch); // Naming in Dhanush's skim
    fEventTree->SetBranchAddress("pPAprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch); // Naming in Kurt's skim
    fBeamScrapingFilterBit = 1; // No beam scraping filter for MC
    fEventTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fCollisionEventSelectionFilterBit = 1;  // No collision event selection filter for pp MC
    fHfCoincidenceFilterBit = 1; // No HF energy coincidence requirement for pp MC
    fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pp MC
  } else if (fDataType == kPbPbMC){ // PbPb data or MC
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
  fEventTree->SetBranchAddress("pt",&fTrackPtArray,&fTrackPtBranch);
  fEventTree->SetBranchAddress("phi",&fTrackPhiArray,&fTrackPhiBranch);
  fEventTree->SetBranchAddress("eta",&fTrackEtaArray,&fTrackEtaBranch);
  fEventTree->SetBranchAddress("chg",&fTrackChargeArray,&fTrackPtErrorBranch); // Reuse a branch from ForestReader that is not otherwise needed here
  fEventTree->SetBranchAddress("sube",&fTrackSubeventArray,&fTrackChi2Branch); // Reuse a branch from ForestReader that is not otherwise needed here
  
}

// Getter for jet raw pT (not relevant for generator jets, just return value that passes cuts)
Float_t GeneratorLevelSkimForestReader::GetJetRawPt(Int_t iJet) const{
  return 2; // The cut is on maxTrackPt/rawPt. Giving 1 and 2 here passes the analysis cuts
}

// Getter for maximum track pT inside a jet (not relevant for generator jets, just return value that passes cuts)
Float_t GeneratorLevelSkimForestReader::GetJetMaxTrackPt(Int_t iJet) const{
  return 1; // The cut is on maxTrackPt/rawPt. Giving 1 and 2 here passes the analysis cuts
}


// Getter for track charge.
Int_t GeneratorLevelSkimForestReader::GetTrackCharge(Int_t iTrack) const{
  return fTrackChargeArray->at(iTrack);
}

// Getter for track subevent index.
Int_t GeneratorLevelSkimForestReader::GetTrackSubevent(Int_t iTrack) const{
  return fTrackSubeventArray->at(iTrack);
}

// Getter for track pT error (not relevant for generator tracks)
Float_t GeneratorLevelSkimForestReader::GetTrackPtError(Int_t iTrack) const{
  return 0; // Setting all errors to 0 always passes the track quality cut
}

// Getter for high purity of the track (not relevant for generator tracks)
Bool_t GeneratorLevelSkimForestReader::GetTrackHighPurity(Int_t iTrack) const{
  return true; // All the generator tracks are of high purity
}

// Getter for track distance from primary vertex in z-direction (not relevant for generator tracks)
Float_t GeneratorLevelSkimForestReader::GetTrackVertexDistanceZ(Int_t iTrack) const{
  return 0; // The cut is on distance/error. Setting this to 0 and error to 1 always passes the cut
}

// Getter for error of track distance from primary vertex in z-direction (not relevant for generator tracks)
Float_t GeneratorLevelSkimForestReader::GetTrackVertexDistanceZError(Int_t iTrack) const{
  return 1; // The cut is on distance/error. Setting this to 1 and distance to 0 always passes the cut
}

// Getter for track distance from primary vertex in xy-direction (not relevant for generator tracks)
Float_t GeneratorLevelSkimForestReader::GetTrackVertexDistanceXY(Int_t iTrack) const{
  return 0; // The cut is on distance/error. Setting this to 0 and error to 1 always passes the cut
}

// Getter for error of track distance from primary vertex in xy-direction (not relevant for generator tracks)
Float_t GeneratorLevelSkimForestReader::GetTrackVertexDistanceXYError(Int_t iTrack) const{
  return 1; // The cut is on distance/error. Setting this to 1 and distance to 0 always passes the cut
}

// Getter for track chi2 value from reconstruction fit (not relevant for generator tracks)
Float_t GeneratorLevelSkimForestReader::GetTrackChi2(Int_t iTrack) const{
  return 1; // Note: The cut on chi2 quality should be disabled in the card for generator tracks
}

// Getter for number of degrees of freedom in reconstruction fit (not relevant for generator tracks)
Int_t GeneratorLevelSkimForestReader::GetNTrackDegreesOfFreedom(Int_t iTrack) const{
  return 1; // Note: The cut on chi2 quality should be disabled in the card for generator tracks
}

// Getter for number of hits in tracker layers (not relevant for generator tracks)
Int_t GeneratorLevelSkimForestReader::GetNHitsTrackerLayer(Int_t iTrack) const{
  return 1; // Note: The cut on chi2 quality should be disabled in the card for generator tracks
}

// Getter for number of hits for the track (not relevant for generator tracks)
Int_t GeneratorLevelSkimForestReader::GetNHitsTrack(Int_t iTrack) const{
  return 1; // Note: The cut on NHits should be disabled in the card for generator tracks
}

// Getter for track energy in ECal (not relevant for generator tracks)
Float_t GeneratorLevelSkimForestReader::GetTrackEnergyEcal(Int_t iTrack) const{
  return 1; // Note: The cut on Et should be disabled on the ConfigurationCard for generator tracks
}

// Getter for track energy in HCal (not relevant for generator tracks)
Float_t GeneratorLevelSkimForestReader::GetTrackEnergyHcal(Int_t iTrack) const{
  return 1; // Note: The cut on Et should be disabled on the ConfigurationCard for generator tracks
}

// Getter for particle flow candidate ID (not relevant for generator tracks)
Int_t GeneratorLevelSkimForestReader::GetParticleFlowCandidateId(Int_t iCandidate) const{
  return -1; // Return negative value to show that correction is not needed
}

// Getter for particle flow candidate pT (not relevant for generator tracks)
Float_t GeneratorLevelSkimForestReader::GetParticleFlowCandidatePt(Int_t iCandidate) const{
  return -1; // Return negative value to show that correction is not needed
}

// Getter for particle flow candidate phi (not relevant for generator tracks)
Float_t GeneratorLevelSkimForestReader::GetParticleFlowCandidatePhi(Int_t iCandidate) const{
  return -1; // Return negative value to show that correction is not needed
}

// Getter for particle flow candidate eta (not relevant for generator tracks)
Float_t GeneratorLevelSkimForestReader::GetParticleFlowCandidateEta(Int_t iCandidate) const{
  return -1; // Return negative value to show that correction is not needed
}

// Getter number of particle flow candidates in an event (not relevant for generator tracks)
Int_t GeneratorLevelSkimForestReader::GetNParticleFlowCandidates() const{
  return -1; // Return negative value to show that correction is not needed
}
