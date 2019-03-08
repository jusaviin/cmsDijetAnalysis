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
 *   Int_t readMode: 0 = Regular forests, 1 = Official PYTHIA8 forest
 *   Int_t jetType: 0 = Calo jets, 1 = PF jets
 *   Int_t jetAxis: 0 = Anti-kT axis, 1 = Leading particle flow candidate axis, 2 = WTA axis
 *   Bool_t matchJets: True = Do matching for reco and gen jets. False = Do not require matching
 */
GeneratorLevelSkimForestReader::GeneratorLevelSkimForestReader(Int_t dataType, Int_t readMode, Int_t jetType, Int_t jetAxis, Bool_t matchJets) :
  SkimForestReader(dataType,readMode,jetType,jetAxis,matchJets),
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
  
  // Only enable the branches that are actually read
  fEventTree->SetBranchStatus("*",0);
  
  // Connect the branches related to event information
  fEventTree->SetBranchStatus("vz",1);
  fEventTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  if(fDataType == kPpMC){
    fHiBin = -1;  // The skims for pp do not have hiBin
  } else {
    fEventTree->SetBranchStatus("hiBin",1);
    fEventTree->SetBranchAddress("hiBin",&fHiBin,&fHiBinBranch);
  }
  fEventTree->SetBranchStatus("pthat",1);
  fEventTree->SetBranchAddress("pthat",&fPtHat,&fPtHatBranch); // pT hat only for MC
  
  // Connect the branches to jet properties
  const char *jetAxis[3] = {"","","_wta_"};
  char branchName[30];
  fEventTree->SetBranchStatus("genpt",1);
  fEventTree->SetBranchAddress("genpt",&fJetPtArray,&fJetPtBranch);
  sprintf(branchName,"gen%sphi",jetAxis[fJetAxis]);
  fEventTree->SetBranchStatus(branchName,1);
  fEventTree->SetBranchAddress(branchName,&fJetPhiArray,&fJetPhiBranch);
  sprintf(branchName,"gen%seta",jetAxis[fJetAxis]);
  fEventTree->SetBranchStatus(branchName,1);
  fEventTree->SetBranchAddress(branchName,&fJetEtaArray,&fJetEtaBranch);
  
  if(fMatchJets){
    const char * jetType[2] = {"calo","pf"};
    jetAxis[0] = "jt"; jetAxis[1] = "jt"; jetAxis[2] = "wta_";
    sprintf(branchName,"%s_refpt",jetType[fJetType]);
    fEventTree->SetBranchStatus(branchName,1);
    fEventTree->SetBranchAddress(branchName,&fJetRefPtArray,&fJetRefPtBranch);
    sprintf(branchName,"%s_refparton_flavor",jetType[fJetType]);
    fEventTree->SetBranchStatus(branchName,1);
    fEventTree->SetBranchAddress(branchName,&fJetRefFlavorArray,&fJetRefFlavorBranch);
    sprintf(branchName,"%s_jtpt",jetType[fJetType]);
    fEventTree->SetBranchStatus(branchName,1);
    fEventTree->SetBranchAddress(branchName,&fMatchedJetPtArray,&fJetMatchedPtBranch);
    sprintf(branchName,"%s_%seta",jetType[fJetType],jetAxis[fJetAxis]);
    fEventTree->SetBranchStatus(branchName,1);
    fEventTree->SetBranchAddress(branchName,&fMatchedJetEtaArray,&fJetMatchedEtaBranch);
    sprintf(branchName,"%s_%sphi",jetType[fJetType],jetAxis[fJetAxis]);
    fEventTree->SetBranchStatus(branchName,1);
    fEventTree->SetBranchAddress(branchName,&fMatchedJetPhiArray,&fJetMatchedPhiBranch);
  }
  
  // Connect the branches to the HLT tree
  fCaloJetFilterBit = 1;         // No calorimeter filter bit is present in the skims
  fCaloJetFilterBitPrescale = 1; // Set the prescaled filter bit to 1. Only relevant for minimum bias PbPb (data skim)
  
  // Connect the branches containing event selection filter bits
  if (fDataType == kPpMC){
    if(fJetAxis == 2){
      fEventTree->SetBranchStatus("pPAprimaryVertexFilter",1);
      fEventTree->SetBranchAddress("pPAprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch); // Naming in Xiao's skim
    } else {
      fEventTree->SetBranchStatus("pprimaryVertexFilter",1);
      fEventTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch); // Naming in regular skim
    }
    fBeamScrapingFilterBit = 1; // No beam scraping filter for MC
    fEventTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
    fEventTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fCollisionEventSelectionFilterBit = 1;  // No collision event selection filter for pp MC
    fHfCoincidenceFilterBit = 1; // No HF energy coincidence requirement for pp MC
    fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pp MC
  } else if (fDataType == kPbPbMC){ // PbPb data or MC
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
  fEventTree->SetBranchStatus("pt",1);
  fEventTree->SetBranchAddress("pt",&fTrackPtArray,&fTrackPtBranch);
  fEventTree->SetBranchStatus("phi",1);
  fEventTree->SetBranchAddress("phi",&fTrackPhiArray,&fTrackPhiBranch);
  fEventTree->SetBranchStatus("eta",1);
  fEventTree->SetBranchAddress("eta",&fTrackEtaArray,&fTrackEtaBranch);
  fEventTree->SetBranchStatus("chg",1);
  fEventTree->SetBranchAddress("chg",&fTrackChargeArray,&fTrackPtErrorBranch); // Reuse a branch from ForestReader that is not otherwise needed here
  fEventTree->SetBranchStatus("sube",1);
  fEventTree->SetBranchAddress("sube",&fTrackSubeventArray,&fTrackChi2Branch); // Reuse a branch from ForestReader that is not otherwise needed here
  
  // Need to check track status for Xiao's skims
  if(fJetAxis == 2){
    fEventTree->SetBranchStatus("status",1);
    fEventTree->SetBranchAddress("status",&fTrackStatusArray,&fnHitsTrackerLayerBranch); // Reuse a branch from ForestReader that is not otherwise needed here
  }
  
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

// Getter for track MC status.
Int_t GeneratorLevelSkimForestReader::GetTrackMCStatus(Int_t iTrack) const{
  if(fJetAxis != 2) return 1;  // Need to check this for Xiao's skims, not for Kurt's skims
  return fTrackStatusArray->at(iTrack);
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
  return 1; // Note: The chi2 quality cut is disabled for generator level tracks in the main analysis code
}

// Getter for number of degrees of freedom in reconstruction fit (not relevant for generator tracks)
Int_t GeneratorLevelSkimForestReader::GetNTrackDegreesOfFreedom(Int_t iTrack) const{
  return 1; // Note: The chi2 quality cut is disabled for generator level tracks in the main analysis code
}

// Getter for number of hits in tracker layers (not relevant for generator tracks)
Int_t GeneratorLevelSkimForestReader::GetNHitsTrackerLayer(Int_t iTrack) const{
  return 1; // Note: The chi2 quality cut is disabled for generator level tracks in the main analysis code
}

// Getter for number of hits for the track (not relevant for generator tracks)
Int_t GeneratorLevelSkimForestReader::GetNHitsTrack(Int_t iTrack) const{
  return 1; // Note: The cut on NHits is disabled for generator level tracks in the main analysis code
}

// Getter for track energy in ECal (not relevant for generator tracks)
Float_t GeneratorLevelSkimForestReader::GetTrackEnergyEcal(Int_t iTrack) const{
  return 1; // Note: The cut on Et is disabled for generator level tracks in the main analysis code
}

// Getter for track energy in HCal (not relevant for generator tracks)
Float_t GeneratorLevelSkimForestReader::GetTrackEnergyHcal(Int_t iTrack) const{
  return 1; // Note: The cut on Et is disabled for generator level tracks in the main analysis code
}

// Getter for particle flow candidate ID (not relevant for generator tracks)
Int_t GeneratorLevelSkimForestReader::GetParticleFlowCandidateId(Int_t iCandidate) const{
  return 1; // Return 1 as we use regular generator lavel tracks to determine leading particle flow candidate
}

// Getter for particle flow candidate pT (just regular tracks in generator level)
Float_t GeneratorLevelSkimForestReader::GetParticleFlowCandidatePt(Int_t iCandidate) const{
  return fTrackPtArray->at(iCandidate);  // Use regular generated particle pT for particle flow cnadidates
}

// Getter for particle flow candidate phi (just regular tracks in generator level)
Float_t GeneratorLevelSkimForestReader::GetParticleFlowCandidatePhi(Int_t iCandidate) const{
  return fTrackPhiArray->at(iCandidate); // Use regular track phi for particle flow candidates
}

// Getter for particle flow candidate eta (just regular tracks in generator level)
Float_t GeneratorLevelSkimForestReader::GetParticleFlowCandidateEta(Int_t iCandidate) const{
  return fTrackEtaArray->at(iCandidate); // Use regular track eta for particle flow candidates
}

// Getter number of particle flow candidates in an event (just regular tracks in generator level)
Int_t GeneratorLevelSkimForestReader::GetNParticleFlowCandidates() const{
  return fnTracks; // Use regular tracks in generator level
}

// Check if generator level jet has a matching reconstructed jet
Bool_t GeneratorLevelSkimForestReader::HasMatchingJet(Int_t iJet) const{
  
  // If not matching jets, just tell that everything is fine
  if(!fMatchJets) return true;
  
  // Ref pT array has pT for all the generator level jets that are matched with reconstructed jets
  // If our generator level pT is found from this array, it has a matching reconstructed jet
  Double_t jetPt = GetJetPt(iJet);
  for(Int_t iRef = 0; iRef < fJetRefPtArray->size(); iRef++){
    if(TMath::Abs(jetPt - fJetRefPtArray->at(iRef)) < 0.001) return true;
  }
  
  return false;
}

// Get the matching reconstructed jet index for the given generator level jet
Int_t GeneratorLevelSkimForestReader::GetMatchingIndex(Int_t iJet) const{
  
  // Ref pT array has pT for all the generator level jets that are matched with reconstructed jets
  // If our generator level pT is found from this array, it has a matching reconstructed jet
  Double_t jetPt = GetJetPt(iJet);
  Int_t matchedIndex = -1;
  for(Int_t iRef = 0; iRef < fJetRefPtArray->size(); iRef++){
    if(TMath::Abs(jetPt - fJetRefPtArray->at(iRef)) < 0.001){
      matchedIndex = iRef;
      break;
    }
  }
  
  return matchedIndex;
  
}

// Get the pT of the matched reconstructed jet
Float_t GeneratorLevelSkimForestReader::GetMatchedPt(Int_t iJet) const{
  
  // If not matching jets, just return something because this has no meaning
  if(!fMatchJets) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchedIndex = GetMatchingIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchedIndex == -1) return -999;
  
  // Return the pT of the matching reconstructed jet
  return fMatchedJetPtArray->at(matchedIndex);
  
}

// Get the eta of the matched reconstructed jet
Float_t GeneratorLevelSkimForestReader::GetMatchedEta(Int_t iJet) const{
  
  // If not matching jets, just return something because this has no meaning
  if(!fMatchJets) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchedIndex = GetMatchingIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchedIndex == -1) return -999;
  
  // Return the pT of the matching reconstructed jet
  return fMatchedJetEtaArray->at(matchedIndex);
  
}

// Get the pT of the matched reconstructed jet
Float_t GeneratorLevelSkimForestReader::GetMatchedPhi(Int_t iJet) const{
  
  // If not matching jets, just return something because this has no meaning
  if(!fMatchJets) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchedIndex = GetMatchingIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchedIndex == -1) return -999;
  
  // Return the pT of the matching reconstructed jet
  return fMatchedJetPhiArray->at(matchedIndex);
  
}

// Get the flavor of the matched jet
Int_t GeneratorLevelSkimForestReader::GetPartonFlavor(Int_t iJet) const{
  
  // If not matching jets, just return something because this has no meaning
  if(!fMatchJets) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchedIndex = GetMatchingIndex(iJet);

  // If we did not find macth, something went wrong. Return -999
  if(matchedIndex == -1) return -999;
  
  // Return the matching parton flavor
  return fJetRefFlavorArray->at(matchedIndex);

}
