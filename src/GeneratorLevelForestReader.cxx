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
  fHltTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fTrackPtArray(0),
  fTrackPhiArray(0),
  fTrackEtaArray(0),
  fTrackChargeArray(0),
  fTrackSubeventArray(0),
  fTrackStatusArray(0)
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
GeneratorLevelForestReader::GeneratorLevelForestReader(Int_t dataType, Int_t readMode, Int_t jetType) :
  ForestReader(dataType,readMode,jetType),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fTrackPtArray(0),
  fTrackPhiArray(0),
  fTrackEtaArray(0),
  fTrackChargeArray(0),
  fTrackSubeventArray(0),
  fTrackStatusArray(0)
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
  fHltTree(in.fHltTree),
  fSkimTree(in.fSkimTree),
  fTrackTree(in.fTrackTree),
  fTrackPtArray(in.fTrackPtArray),
  fTrackPhiArray(in.fTrackPhiArray),
  fTrackEtaArray(in.fTrackEtaArray),
  fTrackChargeArray(in.fTrackChargeArray),
  fTrackSubeventArray(in.fTrackSubeventArray),
  fTrackStatusArray(in.fTrackStatusArray)
{
  // Copy constructor
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
  }
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
  fHltTree = in.fHltTree;
  fSkimTree = in.fSkimTree;
  fTrackTree = in.fTrackTree;
  fTrackPtArray = in.fTrackPtArray;
  fTrackPhiArray = in.fTrackPhiArray;
  fTrackEtaArray = in.fTrackEtaArray;
  fTrackChargeArray = in.fTrackChargeArray;
  fTrackSubeventArray = in.fTrackSubeventArray;
  fTrackStatusArray = in.fTrackStatusArray;
  
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
  }
  
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
  
  // Helper variable for choosing correct branches
  const char *branchName[2] = {"none","none"};
  
  // Connect the branches of the heavy ion tree
  fHeavyIonTree->SetBranchStatus("*",0);
  fHeavyIonTree->SetBranchStatus("vz",1);
  fHeavyIonTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  fHeavyIonTree->SetBranchStatus("hiBin",1);
  fHeavyIonTree->SetBranchAddress("hiBin",&fHiBin,&fHiBinBranch);
  fHeavyIonTree->SetBranchStatus("pthat",1);
  fHeavyIonTree->SetBranchAddress("pthat",&fPtHat,&fPtHatBranch); // pT hat only for MC
  
  // Connect the branches to the jet tree
  fJetTree->SetBranchStatus("*",0);
  fJetTree->SetBranchStatus("genpt",1);
  fJetTree->SetBranchAddress("genpt",&fJetPtArray,&fJetPtBranch);
  fJetTree->SetBranchStatus("genphi",1);
  fJetTree->SetBranchAddress("genphi",&fJetPhiArray,&fJetPhiBranch);
  fJetTree->SetBranchStatus("geneta",1);
  fJetTree->SetBranchAddress("geneta",&fJetEtaArray,&fJetEtaBranch);
  fJetTree->SetBranchStatus("ngen",1);
  fJetTree->SetBranchAddress("ngen",&fnJets,&fJetRawPtBranch); // Reuse a branch from ForestReader that is not otherwise needed here
  
  // Connect the branches to the HLT tree
  fHltTree->SetBranchStatus("*",0);
  if(fDataType == kPp){ // pp data
    branchName[0] = "HLT_AK4CaloJet80_Eta5p1_v1";
    branchName[1] = "HLT_AK4PFJet80_Eta5p1_v1";
    fHltTree->SetBranchStatus(branchName[fJetType],1);
    fHltTree->SetBranchAddress(branchName[fJetType],&fCaloJetFilterBit,&fCaloJetFilterBranch);
  } else if (fDataType == kPpMC){
    branchName[0] = "HLT_AK4CaloJet80_Eta5p1ForPPRef_v1";
    branchName[1] = "HLT_AK4PFJet80_Eta5p1ForPPRef_v1";
    if(fReadMode == 0 || fReadMode == 2){
      fHltTree->SetBranchStatus(branchName[fJetType],1);
      fHltTree->SetBranchAddress(branchName[fJetType],&fCaloJetFilterBit,&fCaloJetFilterBranch);  // For Purdue high forest
    } else {
      fCaloJetFilterBit = 1; // This filter bit does not exist in the official PYTHIA8 dijet forest
    }
  } else if (fDataType == kPbPb){ // PbPb
    fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet100_Eta5p1_v1",1);
    fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1",&fCaloJetFilterBit,&fCaloJetFilterBranch);
  } else if (fDataType == kPbPbMC){
    fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet100_Eta5p1_v2",1);
    fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v2",&fCaloJetFilterBit,&fCaloJetFilterBranch);
  } else { // Local test
    fCaloJetFilterBit = 1;  // No filter for local test
  }
  
  // Connect the branches to the skim tree (different for pp and PbPb Monte Carlo)
  fSkimTree->SetBranchStatus("*",0);
  if(fDataType == kPpMC){ // pp MC
    fSkimTree->SetBranchStatus("pPAprimaryVertexFilter",1);
    fSkimTree->SetBranchAddress("pPAprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fSkimTree->SetBranchStatus("pBeamScrapingFilter",1);
    fSkimTree->SetBranchAddress("pBeamScrapingFilter",&fBeamScrapingFilterBit,&fBeamScrapingBranch);
    fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
    fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fCollisionEventSelectionFilterBit = 1;  // No collision event selection filter for pp
    fHfCoincidenceFilterBit = 1; // No HF energy coincidence requirement for pp
    fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pp
  } else if (fDataType == kPbPbMC){ // PbPb MC
    fSkimTree->SetBranchStatus("pprimaryVertexFilter",1);
    fSkimTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
    fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fSkimTree->SetBranchStatus("pcollisionEventSelection",1);
    fSkimTree->SetBranchAddress("pcollisionEventSelection",&fCollisionEventSelectionFilterBit,&fCollisionEventSelectionBranch);
    fSkimTree->SetBranchStatus("phfCoincFilter3",1);
    fSkimTree->SetBranchAddress("phfCoincFilter3",&fHfCoincidenceFilterBit,&fHfCoincidenceBranch);
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
  if(fDataType != kLocalTest && fReadMode == 0) {
    fTrackTree->SetBranchStatus("sta",1);
    fTrackTree->SetBranchAddress("sta",&fTrackStatusArray,&fTrackEnergyEcalBranch); // Reuse a branch from ForestReader that is not otherwise needed here. Not available for local test or PYTHIA8 forest
  }

}

/*
 * Connect a new tree to the reader
 */
void GeneratorLevelForestReader::ReadForestFromFile(TFile *inputFile){
  
  // Helper variable for finding the correct tree
  const char *treeName[2] = {"none","none"};
  
  // Connect a trees from the file to the reader
  fHeavyIonTree = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
  fHltTree = (TTree*)inputFile->Get("hltanalysis/HltTree");
  fSkimTree = (TTree*)inputFile->Get("skimanalysis/HltTree");
  
  // The jet tree has different name in different datasets
  if(fDataType == kPp || fDataType == kPpMC){
    treeName[0] = "ak4CaloJetAnalyzer/t"; // Tree for calo jets
    treeName[1] = "ak4PFJetAnalyzer/t";   // Tree for PF jets
  } else if (fDataType == kPbPb || fDataType == kPbPbMC){
    treeName[0] = "akPu4CaloJetAnalyzer/t";  // Tree for calo jets
    treeName[1] = "akPu4PFJetAnalyzer/t";    // Tree for PF jets
  } else if (fDataType == kLocalTest){
    treeName[0] = "ak4PFJetAnalyzer/t";  // Only PF jets in local test file
    treeName[1] = "ak4PFJetAnalyzer/t";  // Only PF jets in local test file
  }
  fJetTree = (TTree*)inputFile->Get(treeName[fJetType]);
  
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
  fHltTree->Delete();
  fSkimTree->Delete();
  fTrackTree->Delete();
}

/*
 * Load an event to memory
 */
void GeneratorLevelForestReader::GetEvent(Int_t nEvent){
  fHeavyIonTree->GetEntry(nEvent);
  fJetTree->GetEntry(nEvent);
  fHltTree->GetEntry(nEvent);
  fSkimTree->GetEntry(nEvent);
  fTrackTree->GetEntry(nEvent);
  
  // Read the numbers of tracks and jets for this event
  fnTracks = fTrackPtArray->size();
}

// Getter for jet pT
Float_t GeneratorLevelForestReader::GetJetPt(Int_t iJet) const{
  return fJetPtArray[iJet];
}

// Getter for jet phi
Float_t GeneratorLevelForestReader::GetJetPhi(Int_t iJet) const{
  return fJetPhiArray[iJet];
}

// Getter for jet eta
Float_t GeneratorLevelForestReader::GetJetEta(Int_t iJet) const{
  return fJetEtaArray[iJet];
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

// Getter for track MC status
Int_t GeneratorLevelForestReader::GetTrackMCStatus(Int_t iTrack) const{
  if(fDataType == kLocalTest || fReadMode == 1 || fReadMode == 2) return 1;
  return fTrackStatusArray->at(iTrack);
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
  return 1; // Note: The chi2 quality cut is disabled for generator level tracks in the main analysis code
}

// Getter for number of degrees of freedom in reconstruction fit (not relevant for generator tracks)
Int_t GeneratorLevelForestReader::GetNTrackDegreesOfFreedom(Int_t iTrack) const{
  return 1; // Note: The chi2 quality cut is disabled for generator level tracks in the main analysis code
}

// Getter for number of hits in tracker layers (not relevant for generator tracks)
Int_t GeneratorLevelForestReader::GetNHitsTrackerLayer(Int_t iTrack) const{
  return 1; // Note: The chi2 quality cut is disabled for generator level tracks in the main analysis code
}

// Getter for number of hits for the track (not relevant for generator tracks)
Int_t GeneratorLevelForestReader::GetNHitsTrack(Int_t iTrack) const{
  return 1; // Note: The cut on NHits is disabled for generator level tracks in the main analysis code
}

// Getter for track energy in ECal (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackEnergyEcal(Int_t iTrack) const{
  return 1; // Note: The cut on Et is disabled for generator level tracks in the main analysis code
}

// Getter for track energy in HCal (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackEnergyHcal(Int_t iTrack) const{
  return 1; // Note: The cut on Et is disabled for generator level tracks in the main analysis code
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
