// Implementation for GeneratorLevelForestReader

// Own includes
#include "GeneratorLevelForestReader.h"

/*
 * Default constructor
 */
GeneratorLevelForestReader::GeneratorLevelForestReader() :
  ForestReader(),
  fHeavyIonTree(0),
  fJetTree(0),
  fTrackTree(0),
  fJetPtArray(0),
  fJetPhiArray(0),
  fJetEtaArray(0),
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
 */
GeneratorLevelForestReader::GeneratorLevelForestReader(Int_t dataType) :
  ForestReader(dataType),
  fHeavyIonTree(0),
  fJetTree(0),
  fTrackTree(0),
  fJetPtArray(0),
  fJetPhiArray(0),
  fJetEtaArray(0),
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
GeneratorLevelForestReader::GeneratorLevelForestReader(const GeneratorLevelForestReader& in) :
  ForestReader(in),
  fHeavyIonTree(in.fHeavyIonTree),
  fJetTree(in.fJetTree),
  fTrackTree(in.fTrackTree),
  fJetPtArray(in.fJetPtArray),
  fJetPhiArray(in.fJetPhiArray),
  fJetEtaArray(in.fJetEtaArray),
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
GeneratorLevelForestReader& GeneratorLevelForestReader::operator=(const GeneratorLevelForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  ForestReader::operator=(in);
  
  fHeavyIonTree = in.fHeavyIonTree;
  fJetTree = in.fJetTree;
  fTrackTree = in.fTrackTree;
  fJetPtArray = in.fJetPtArray;
  fJetPhiArray = in.fJetPhiArray;
  fJetEtaArray = in.fJetEtaArray;
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
GeneratorLevelForestReader::~GeneratorLevelForestReader(){
  // destructor
}

/*
 * Initialization, meaning that the branches are connected to the tree
 */
void GeneratorLevelForestReader::Initialize(){
  
  // Connect the branches of the heavy ion tree
  fHeavyIonTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  fHeavyIonTree->SetBranchAddress("hiBin",&fHiBin,&fHiBinBranch);
  fHeavyIonTree->SetBranchAddress("pthat",&fPtHat,&fPtHatBranch); // pT hat only for MC
  
  // Connect the branches to the jet tree
  fJetTree->SetBranchAddress("jtpt",&fJetPtArray,&fJetPtBranch);
  fJetTree->SetBranchAddress("jtphi",&fJetPhiArray,&fJetPhiBranch);
  fJetTree->SetBranchAddress("jteta",&fJetEtaArray,&fJetEtaBranch);
  
  // Connect the branches to the HLT tree (only for real data)
  fCaloJetFilterBit = 1;  // No filter for Monte Carlo
  
  // Connect the branches to the skim tree (different for pp and PbPb data, no connection for Monte Carlo)
  fPrimaryVertexFilterBit = 1;            // No filter for Monte Carlo
  fBeamScrapingFilterBit = 1;             // No filter for Monte Carlo
  fCollisionEventSelectionFilterBit = 1;  // No filter for Monte Carlo
  fHBHENoiseFilterBit = 1;                // No filter for Monte Carlo
  fHfCoincidenceFilterBit = 1;            // No HF energy coincidence requirement for Monte Carlo
  fClusterCompatibilityFilterBit = 1;     // No cluster compatibility requirement for Monte Carlo
  
  // Connect the branches to the track tree
  fTrackTree->SetBranchAddress("pt",&fTrackPtArray,&fTrackPtBranch);
  fTrackTree->SetBranchAddress("phi",&fTrackPhiArray,&fTrackPhiBranch);
  fTrackTree->SetBranchAddress("eta",&fTrackEtaArray,&fTrackEtaBranch);
  fTrackTree->SetBranchAddress("chg",&fTrackChargeArray,&fTrackPtErrorBranch); // Reuse a branch from ForestReader that is not otherwise needed here
  fTrackTree->SetBranchAddress("sube",&fTrackSubeventArray,&fTrackChi2Branch); // Reuse a branch from ForestReader that is not otherwise needed here
}

/*
 * Connect a new tree to the reader
 */
void GeneratorLevelForestReader::ReadForestFromFile(TFile *inputFile){
  
  // Connect a trees from the file to the reader
  fHeavyIonTree = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
  
  // The jet tree has different name in different datasets
  if(fDataType == kPpMC){
    fJetTree = (TTree*)inputFile->Get("ak4CaloJetAnalyzer/t");
  } else if (fDataType == kPbPbMC){
    fJetTree = (TTree*)inputFile->Get("akPu4CaloJetAnalyzer/t");
  }
  
  // The track tree is different than in other types of forests
  fTrackTree = (TTree*)inputFile->Get("HiGenParticleAna/hi");
  
  // Connect branches to trees
  Initialize();
}

/*
 * Burn the current forest.
 */
void GeneratorLevelForestReader::BurnForest(){
  fHeavyIonTree->Delete();
  fJetTree->Delete();
  fTrackTree->Delete();
}

/*
 * Load an event to memory
 */
void GeneratorLevelForestReader::GetEvent(Int_t nEvent){
  fHeavyIonTree->GetEntry(nEvent);
  fJetTree->GetEntry(nEvent);
  fTrackTree->GetEntry(nEvent);
  
  // Read the numbers of tracks and jets for this event
  fnJets = fJetPtArray->size();
  fnTracks = fTrackPtArray->size();
}

// Getter for jet pT
Float_t GeneratorLevelForestReader::GetJetPt(Int_t iJet) const{
  return fJetPtArray->at(iJet);
}

// Getter for jet phi
Float_t GeneratorLevelForestReader::GetJetPhi(Int_t iJet) const{
  return fJetPhiArray->at(iJet);
}

// Getter for jet eta
Float_t GeneratorLevelForestReader::GetJetEta(Int_t iJet) const{
  return fJetEtaArray->at(iJet);
}

// Getter for jet raw pT (not relevant for generator jets, just return value that passes cuts)
Float_t GeneratorLevelForestReader::GetJetRawPt(Int_t iJet) const{
  return 2; // The cut is on maxTrackPt/rawPt. Giving 1 and 2 here passes the analysis cuts
}

// Getter for maximum track pT inside a jet (not relevant for generator jets, just return value that passes cuts)
Float_t GeneratorLevelForestReader::GetJetMaxTrackPt(Int_t iJet) const{
  return 1; // The cut is on maxTrackPt/rawPt. Giving 1 and 2 here passes the analysis cuts
}

// Getter for track pT
Float_t GeneratorLevelForestReader::GetTrackPt(Int_t iTrack) const{
  return fTrackPtArray->at(iTrack);
}

// Getter for track phi
Float_t GeneratorLevelForestReader::GetTrackPhi(Int_t iTrack) const{
  return fTrackPhiArray->at(iTrack);
}

// Getter for track eta
Float_t GeneratorLevelForestReader::GetTrackEta(Int_t iTrack) const{
  return fTrackEtaArray->at(iTrack);
}

// Getter for track charge
Int_t GeneratorLevelForestReader::GetTrackCharge(Int_t iTrack) const{
  return fTrackChargeArray->at(iTrack);
}

// Getter for track subevent index
Int_t GeneratorLevelForestReader::GetTrackSubevent(Int_t iTrack) const{
  return fTrackSubeventArray->at(iTrack);
}

// Getter for track pT error (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackPtError(Int_t iTrack) const{
  return 0; // Setting all errors to 0 always passes the track quality cut
}

// Getter for high purity of the track (not relevant for generator tracks)
Bool_t GeneratorLevelForestReader::GetTrackHighPurity(Int_t iTrack) const{
  return true; // All the generator tracks are of high purity
}

// Getter for track distance from primary vertex in z-direction (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackVertexDistanceZ(Int_t iTrack) const{
  return 0; // The cut is on distance/error. Setting this to 0 and error to 1 always passes the cut
}

// Getter for error of track distance from primary vertex in z-direction (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackVertexDistanceZError(Int_t iTrack) const{
  return 1; // The cut is on distance/error. Setting this to 1 and distance to 0 always passes the cut
}

// Getter for track distance from primary vertex in xy-direction (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackVertexDistanceXY(Int_t iTrack) const{
  return 0; // The cut is on distance/error. Setting this to 0 and error to 1 always passes the cut
}

// Getter for error of track distance from primary vertex in xy-direction (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackVertexDistanceXYError(Int_t iTrack) const{
  return 1; // The cut is on distance/error. Setting this to 1 and distance to 0 always passes the cut
}

// Getter for track chi2 value from reconstruction fit (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackChi2(Int_t iTrack) const{
  return 1; // Note: The cut on chi2 quality should be disabled in the card for generator tracks
}

// Getter for number of degrees of freedom in reconstruction fit (not relevant for generator tracks)
Int_t GeneratorLevelForestReader::GetNTrackDegreesOfFreedom(Int_t iTrack) const{
  return 1; // Note: The cut on chi2 quality should be disabled in the card for generator tracks
}

// Getter for number of hits in tracker layers (not relevant for generator tracks)
Int_t GeneratorLevelForestReader::GetNHitsTrackerLayer(Int_t iTrack) const{
  return 1; // Note: The cut on chi2 quality should be disabled in the card for generator tracks
}

// Getter for number of hits for the track (not relevant for generator tracks)
Int_t GeneratorLevelForestReader::GetNHitsTrack(Int_t iTrack) const{
  return 1; // Note: The cut on NHits should be disabled in the card for generator tracks
}

// Getter for track energy in ECal (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackEnergyEcal(Int_t iTrack) const{
  return 1; // Note: The cut on Et should be disabled on the ConfigurationCard for generator tracks
}

// Getter for track energy in HCal (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackEnergyHcal(Int_t iTrack) const{
  return 1; // Note: The cut on Et should be disabled on the ConfigurationCard for generator tracks
}

// Getter for particle flow candidate ID (not relevant for generator tracks)
Int_t GeneratorLevelForestReader::GetParticleFlowCandidateId(Int_t iCandidate) const{
  return -1; // Return negative value to show that correction is not needed
}

// Getter for particle flow candidate pT (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetParticleFlowCandidatePt(Int_t iCandidate) const{
  return -1; // Return negative value to show that correction is not needed
}

// Getter for particle flow candidate phi (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetParticleFlowCandidatePhi(Int_t iCandidate) const{
  return -1; // Return negative value to show that correction is not needed
}

// Getter for particle flow candidate eta (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetParticleFlowCandidateEta(Int_t iCandidate) const{
  return -1; // Return negative value to show that correction is not needed
}

// Getter number of particle flow candidates in an event (not relevant for generator tracks)
Int_t GeneratorLevelForestReader::GetNParticleFlowCandidates() const{
  return -1; // Return negative value to show that correction is not needed
}
