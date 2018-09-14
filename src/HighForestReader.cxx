// Implementation for HighForestReader

// Own includes
#include "HighForestReader.h"

/*
 * Default constructor
 */
HighForestReader::HighForestReader() :
  ForestReader(),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fParticleFlowCandidateTree(0),
  fnJetsBranch(0),
  fnTracksBranch(0),
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
 *   Int_t readMode: 0 = Regular forests, 1 = Official PYTHIA8 forest
 *   Int_t jetType: 0 = Calo jets, 1 = PF jets
 */
HighForestReader::HighForestReader(Int_t dataType, Int_t readMode, Int_t jetType) :
  ForestReader(dataType,readMode,jetType),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fParticleFlowCandidateTree(0),
  fnJetsBranch(0),
  fnTracksBranch(0),
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
HighForestReader::HighForestReader(const HighForestReader& in) :
  ForestReader(in),
  fHeavyIonTree(in.fHeavyIonTree),
  fJetTree(in.fJetTree),
  fHltTree(in.fHltTree),
  fSkimTree(in.fSkimTree),
  fTrackTree(in.fTrackTree),
  fParticleFlowCandidateTree(in.fParticleFlowCandidateTree),
  fnJetsBranch(in.fnJetsBranch),
  fnTracksBranch(in.fnTracksBranch)
{
  // Copy constructor
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetRawPtArray[i] = in.fJetRawPtArray[i];
    fJetMaxTrackPtArray[i] = in.fJetMaxTrackPtArray[i];
  }
  
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
  }
}

/*
 * Assignment operator
 */
HighForestReader& HighForestReader::operator=(const HighForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  ForestReader::operator=(in);
  
  fHeavyIonTree = in.fHeavyIonTree;
  fJetTree = in.fJetTree;
  fHltTree = in.fHltTree;
  fSkimTree = in.fSkimTree;
  fTrackTree = in.fTrackTree;
  fParticleFlowCandidateTree = in.fParticleFlowCandidateTree;
  fnJetsBranch = in.fnJetsBranch;
  fnTracksBranch = in.fnTracksBranch;
  
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetRawPtArray[i] = in.fJetRawPtArray[i];
    fJetMaxTrackPtArray[i] = in.fJetMaxTrackPtArray[i];
  }
  
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
  }
  
  return *this;
}

/*
 * Destructor
 */
HighForestReader::~HighForestReader(){
  // destructor
}

/*
 * Initialization, meaning that the branches are connected to the tree
 */
void HighForestReader::Initialize(){
  
  // Helper variable for choosing correct branches
  const char *branchName[2] = {"none","none"};
  
  // Connect the branches of the heavy ion tree
  fHeavyIonTree->SetBranchStatus("*",0);
  fHeavyIonTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  fHeavyIonTree->SetBranchStatus("vz",1);
  fHeavyIonTree->SetBranchAddress("hiBin",&fHiBin,&fHiBinBranch);
  fHeavyIonTree->SetBranchStatus("hiBin",1);
  if(fDataType == kPpMC || fDataType == kPbPbMC || fDataType == kLocalTest){
    fHeavyIonTree->SetBranchAddress("pthat",&fPtHat,&fPtHatBranch); // pT hat only for MC
    fHeavyIonTree->SetBranchStatus("pthat",1);
  } else {
    fPtHat = 0; // We do not have pT hat information for real data
  }
  
  // Connect the branches to the jet tree
  fJetTree->SetBranchStatus("*",0);
  fJetTree->SetBranchAddress("jtpt",&fJetPtArray,&fJetPtBranch);
  fJetTree->SetBranchStatus("jtpt",1);
  fJetTree->SetBranchAddress("jtphi",&fJetPhiArray,&fJetPhiBranch);
  fJetTree->SetBranchStatus("jtphi",1);
  fJetTree->SetBranchAddress("jteta",&fJetEtaArray,&fJetEtaBranch);
  fJetTree->SetBranchStatus("jteta",1);
  fJetTree->SetBranchAddress("nref",&fnJets,&fnJetsBranch);
  fJetTree->SetBranchStatus("nref",1);
  fJetTree->SetBranchAddress("rawpt",&fJetRawPtArray,&fJetRawPtBranch);
  fJetTree->SetBranchStatus("rawpt",1);
  fJetTree->SetBranchAddress("trackMax",&fJetMaxTrackPtArray,&fJetMaxTrackPtBranch);
  fJetTree->SetBranchStatus("trackMax",1);
  
  // Event selection summary
  //
  //         tree                      branch                         What it is
  //  hltanalysis/HltTree   HLT_HIPuAK4CaloJet100_Eta5p1_v1      Event selection for PbPb
  //  hltanalysis/HltTree      HLT_AK4CaloJet80_Eta5p1_v1         Event selection for pp
  // skimanalysis/HltTree         pprimaryVertexFilter           Event selection for PbPb
  // skimanalysis/HltTree       pcollisionEventSelection         Event selection for PbPb
  // skimanalysis/HltTree    HBHENoiseFilterResultRun2Loose   Event selection for pp and PbPb
  // skimanalysis/HltTree         pPAprimaryVertexFilter          Event selection for pp
  // skimanalysis/HltTree           pBeamScrapingFilter           Event selection for pp
  
  // Connect the branches to the HLT tree
  fHltTree->SetBranchStatus("*",0);
  if(fDataType == kPp){ // pp data
    branchName[0] = "HLT_AK4CaloJet80_Eta5p1_v1";
    branchName[1] = "HLT_AK4PFJet80_Eta5p1_v1";
    fHltTree->SetBranchAddress(branchName[fJetType],&fCaloJetFilterBit,&fCaloJetFilterBranch);
    fHltTree->SetBranchStatus(branchName[fJetType],1);
  } else if (fDataType == kPpMC){
    branchName[0] = "HLT_AK4CaloJet80_Eta5p1ForPPRef_v1";
    branchName[1] = "HLT_AK4PFJet80_Eta5p1ForPPRef_v1";
    if(fReadMode == 0 || fReadMode == 2){
      fHltTree->SetBranchAddress(branchName[fJetType],&fCaloJetFilterBit,&fCaloJetFilterBranch);  // For Purdue high forest
      fHltTree->SetBranchStatus(branchName[fJetType],1);
    } else {
      fCaloJetFilterBit = 1; // This filter bit does not exist in the official PYTHIA8 dijet forest
    }
  } else if (fDataType == kPbPb){ // PbPb
    fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1",&fCaloJetFilterBit,&fCaloJetFilterBranch);
    fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet100_Eta5p1_v1",1);
  } else if (fDataType == kPbPbMC){
    fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v2",&fCaloJetFilterBit,&fCaloJetFilterBranch);
    fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet100_Eta5p1_v2",1);
  } else { // Local test
    fCaloJetFilterBit = 1;  // No filter for local test
  }
  
  // Connect the branches to the skim tree (different for pp and PbPb data and Monte Carlo)
  fSkimTree->SetBranchStatus("*",0);
  if(fDataType == kPp || fDataType == kPpMC){ // pp data or MC
    fSkimTree->SetBranchAddress("pPAprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fSkimTree->SetBranchStatus("pPAprimaryVertexFilter",1);
    fSkimTree->SetBranchAddress("pBeamScrapingFilter",&fBeamScrapingFilterBit,&fBeamScrapingBranch);
    fSkimTree->SetBranchStatus("pBeamScrapingFilter",1);
    fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
    fCollisionEventSelectionFilterBit = 1;  // No collision event selection filter for pp
    fHfCoincidenceFilterBit = 1; // No HF energy coincidence requirement for pp
    fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pp
  } else if (fDataType == kPbPb || fDataType == kPbPbMC){ // PbPb data or MC
    fSkimTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fSkimTree->SetBranchStatus("pprimaryVertexFilter",1);
    fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
    fSkimTree->SetBranchAddress("pcollisionEventSelection",&fCollisionEventSelectionFilterBit,&fCollisionEventSelectionBranch);
    fSkimTree->SetBranchStatus("pcollisionEventSelection",1);
    fSkimTree->SetBranchAddress("phfCoincFilter3",&fHfCoincidenceFilterBit,&fHfCoincidenceBranch);
    fSkimTree->SetBranchStatus("phfCoincFilter3",1);
    fSkimTree->SetBranchAddress("pclusterCompatibilityFilter",&fClusterCompatibilityFilterBit,&fClusterCompatibilityBranch);
    fSkimTree->SetBranchStatus("pclusterCompatibilityFilter",1);
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
  fTrackTree->SetBranchAddress("trkPt",&fTrackPtArray,&fTrackPtBranch);
  fTrackTree->SetBranchStatus("trkPt",1);
  fTrackTree->SetBranchAddress("trkPtError",&fTrackPtErrorArray,&fTrackPtErrorBranch);
  fTrackTree->SetBranchStatus("trkPtError",1);
  fTrackTree->SetBranchAddress("trkPhi",&fTrackPhiArray,&fTrackPhiBranch);
  fTrackTree->SetBranchStatus("trkPhi",1);
  fTrackTree->SetBranchAddress("trkEta",&fTrackEtaArray,&fTrackEtaBranch);
  fTrackTree->SetBranchStatus("trkEta",1);
  fTrackTree->SetBranchAddress("nTrk",&fnTracks,&fnTracksBranch);
  fTrackTree->SetBranchStatus("nTrk",1);
  fTrackTree->SetBranchAddress("highPurity",&fHighPurityTrackArray,&fHighPurityTrackBranch);
  fTrackTree->SetBranchStatus("highPurity",1);
  fTrackTree->SetBranchAddress("trkDz1",&fTrackVertexDistanceZArray,&fTrackVertexDistanceZBranch);
  fTrackTree->SetBranchStatus("trkDz1",1);
  fTrackTree->SetBranchAddress("trkDzError1",&fTrackVertexDistanceZErrorArray,&fTrackVertexDistanceZErrorBranch);
  fTrackTree->SetBranchStatus("trkDzError1",1);
  fTrackTree->SetBranchAddress("trkDxy1",&fTrackVertexDistanceXYArray,&fTrackVertexDistanceXYBranch);
  fTrackTree->SetBranchStatus("trkDxy1",1);
  fTrackTree->SetBranchAddress("trkDxyError1",&fTrackVertexDistanceXYErrorArray,&fTrackVertexDistanceXYErrorBranch);
  fTrackTree->SetBranchStatus("trkDxyError1",1);
  fTrackTree->SetBranchAddress("trkChi2",&fTrackChi2Array,&fTrackChi2Branch);
  fTrackTree->SetBranchStatus("trkChi2",1);
  fTrackTree->SetBranchAddress("trkNdof",&fnTrackDegreesOfFreedomArray,&fnTrackDegreesOfFreedomBranch);
  fTrackTree->SetBranchStatus("trkNdof",1);
  fTrackTree->SetBranchAddress("trkNlayer",&fnHitsTrackerLayerArray,&fnHitsTrackerLayerBranch);
  fTrackTree->SetBranchStatus("trkNlayer",1);
  fTrackTree->SetBranchAddress("trkNHit",&fnHitsTrackArray,&fnHitsTrackBranch);
  fTrackTree->SetBranchStatus("trkNHit",1);
  fTrackTree->SetBranchAddress("pfEcal",&fTrackEnergyEcalArray,&fTrackEnergyEcalBranch);
  fTrackTree->SetBranchStatus("pfEcal",1);
  fTrackTree->SetBranchAddress("pfHcal",&fTrackEnergyHcalArray,&fTrackEnergyHcalBranch);
  fTrackTree->SetBranchStatus("pfHcal",1);
  
  // Connect the branches to the particle flow candidate tree
  fParticleFlowCandidateTree->SetBranchStatus("*",0);
  if(fReadMode == 0 || fReadMode == 2){ // Regular forests have vectors for particle flow candidate tree
    fParticleFlowCandidateTree->SetBranchAddress("pfId",&fParticleFlowCandidateIdVector,&fParticleFlowCandidateIdBranch);
    fParticleFlowCandidateTree->SetBranchStatus("pfId",1);
    fParticleFlowCandidateTree->SetBranchAddress("pfPt",&fParticleFlowCandidatePtVector,&fParticleFlowCandidatePtBranch);
    fParticleFlowCandidateTree->SetBranchStatus("pfPt",1);
    fParticleFlowCandidateTree->SetBranchAddress("pfPhi",&fParticleFlowCandidatePhiVector,&fParticleFlowCandidatePhiBranch);
    fParticleFlowCandidateTree->SetBranchStatus("pfPhi",1);
    fParticleFlowCandidateTree->SetBranchAddress("pfEta",&fParticleFlowCandidateEtaVector,&fParticleFlowCandidateEtaBranch);
    fParticleFlowCandidateTree->SetBranchStatus("pfEta",1);
  } else { // PYTHIA8 forest has arrays instead of vectors for particle flow candidate tree
    fParticleFlowCandidateTree->SetBranchAddress("nPFpart",&fnParticleFlowCandidates,&nfParticleFlowCandidateBranch);
    fParticleFlowCandidateTree->SetBranchStatus("nPFpart",1);
    fParticleFlowCandidateTree->SetBranchAddress("pfId",&fParticleFlowCandidateIdArray,&fParticleFlowCandidateIdBranch);
    fParticleFlowCandidateTree->SetBranchStatus("pfId",1);
    fParticleFlowCandidateTree->SetBranchAddress("pfPt",&fParticleFlowCandidatePtArray,&fParticleFlowCandidatePtBranch);
    fParticleFlowCandidateTree->SetBranchStatus("pfPt",1);
    fParticleFlowCandidateTree->SetBranchAddress("pfPhi",&fParticleFlowCandidatePhiArray,&fParticleFlowCandidatePhiBranch);
    fParticleFlowCandidateTree->SetBranchStatus("pfPhi",1);
    fParticleFlowCandidateTree->SetBranchAddress("pfEta",&fParticleFlowCandidateEtaArray,&fParticleFlowCandidateEtaBranch);
    fParticleFlowCandidateTree->SetBranchStatus("pfEta",1);
  }
  
}

/*
 * Connect a new tree to the reader
 */
void HighForestReader::ReadForestFromFile(TFile *inputFile){
  
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
  
  // The track tree and the particle flow candidate tree have different names for differant datasets
  if(fDataType == kPp || fDataType == kPpMC || fDataType == kLocalTest){
    fTrackTree = (TTree*)inputFile->Get("ppTrack/trackTree");
    fParticleFlowCandidateTree = (TTree*)inputFile->Get("pfcandAnalyzer/pfTree");
  } else if (fDataType == kPbPb || fDataType == kPbPbMC){
    fTrackTree = (TTree*)inputFile->Get("anaTrack/trackTree");
    fParticleFlowCandidateTree = (TTree*)inputFile->Get("pfcandAnalyzerCS/pfTree");
  }
  
  
  Initialize();
}

/*
 * Burn the current forest.
 */
void HighForestReader::BurnForest(){
  fHeavyIonTree->Delete();
  fHltTree->Delete();
  fSkimTree->Delete();
  fJetTree->Delete();
  fTrackTree->Delete();
  fParticleFlowCandidateTree->Delete();
}

/*
 * Load an event to memory
 */
void HighForestReader::GetEvent(Int_t nEvent){
  fHeavyIonTree->GetEntry(nEvent);
  fJetTree->GetEntry(nEvent);
  fHltTree->GetEntry(nEvent);
  fSkimTree->GetEntry(nEvent);
  fTrackTree->GetEntry(nEvent);
  fParticleFlowCandidateTree->GetEntry(nEvent);
}

// Getter for jet pT
Float_t HighForestReader::GetJetPt(Int_t iJet) const{
  return fJetPtArray[iJet];
}

// Getter for jet phi
Float_t HighForestReader::GetJetPhi(Int_t iJet) const{
  return fJetPhiArray[iJet];
}

// Getter for jet eta
Float_t HighForestReader::GetJetEta(Int_t iJet) const{
  return fJetEtaArray[iJet];
}

// Getter for jet raw pT
Float_t HighForestReader::GetJetRawPt(Int_t iJet) const{
  return fJetRawPtArray[iJet];
}

// Getter for maximum track pT inside a jet
Float_t HighForestReader::GetJetMaxTrackPt(Int_t iJet) const{
  return fJetMaxTrackPtArray[iJet];
}

// Getter for track pT
Float_t HighForestReader::GetTrackPt(Int_t iTrack) const{
  return fTrackPtArray[iTrack];
}

// Getter for track pT error
Float_t HighForestReader::GetTrackPtError(Int_t iTrack) const{
  return fTrackPtErrorArray[iTrack];
}

// Getter for track phi
Float_t HighForestReader::GetTrackPhi(Int_t iTrack) const{
  return fTrackPhiArray[iTrack];
}

// Getter for track eta
Float_t HighForestReader::GetTrackEta(Int_t iTrack) const{
  return fTrackEtaArray[iTrack];
}

// Getter for high purity of the track
Bool_t HighForestReader::GetTrackHighPurity(Int_t iTrack) const{
  return fHighPurityTrackArray[iTrack];
}

// Getter for track distance from primary vertex in z-direction
Float_t HighForestReader::GetTrackVertexDistanceZ(Int_t iTrack) const{
  return fTrackVertexDistanceZArray[iTrack];
}

// Getter for error of track distance from primary vertex in z-direction
Float_t HighForestReader::GetTrackVertexDistanceZError(Int_t iTrack) const{
  return fTrackVertexDistanceZErrorArray[iTrack];
}

// Getter for track distance from primary vertex in xy-direction
Float_t HighForestReader::GetTrackVertexDistanceXY(Int_t iTrack) const{
  return fTrackVertexDistanceXYArray[iTrack];
}

// Getter for error of track distance from primary vertex in xy-direction
Float_t HighForestReader::GetTrackVertexDistanceXYError(Int_t iTrack) const{
  return fTrackVertexDistanceXYErrorArray[iTrack];
}

// Getter for track chi2 value from reconstruction fit
Float_t HighForestReader::GetTrackChi2(Int_t iTrack) const{
  return fTrackChi2Array[iTrack];
}

// Getter for number of degrees of freedom in reconstruction fit
Int_t HighForestReader::GetNTrackDegreesOfFreedom(Int_t iTrack) const{
  return fnTrackDegreesOfFreedomArray[iTrack];
}

// Getter for number of hits in tracker layers
Int_t HighForestReader::GetNHitsTrackerLayer(Int_t iTrack) const{
  return fnHitsTrackerLayerArray[iTrack];
}

// Getter for number of hits for the track
Int_t HighForestReader::GetNHitsTrack(Int_t iTrack) const{
  return fnHitsTrackArray[iTrack];
}

// Getter for track energy in ECal
Float_t HighForestReader::GetTrackEnergyEcal(Int_t iTrack) const{
  return fTrackEnergyEcalArray[iTrack];
}

// Getter for track energy in HCal
Float_t HighForestReader::GetTrackEnergyHcal(Int_t iTrack) const{
  return fTrackEnergyHcalArray[iTrack];
}

// Getter for track charge. Charge is only relevant for generator level tracks.
Int_t HighForestReader::GetTrackCharge(Int_t iTrack) const{
  return 1;
}

// Getter for track subevent index. Relevant only for generator level tracks.
Int_t HighForestReader::GetTrackSubevent(Int_t iTrack) const{
  return -1;
}

// Getter for track MC status. Relevant only for generator level tracks.
Int_t HighForestReader::GetTrackMCStatus(Int_t iTrack) const{
  return 1;
}
