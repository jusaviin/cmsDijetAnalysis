// Implementation for ForestReader

// Own includes
#include "ForestReader.h"

/*
 * Default constructor
 */
ForestReader::ForestReader() :
  fDataType(0),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fHiVzBranch(0),
  fHiBinBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetEtaBranch(0),
  fnJetsBranch(0),
  fJetRawPtBranch(0),
  fJetMaxTrackPtBranch(0),
  fCaloJetFilterBranch(0),
  fPrimaryVertexBranch(0),
  fBeamScrapingBranch(0),
  fCollisionEventSelectionBranch(0),
  fHBHENoiseBranch(0),
  fTrackPtBranch(0),
  fTrackPtErrorBranch(0),
  fTrackPhiBranch(0),
  fTrackEtaBranch(0),
  fnTracksBranch(0),
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
  fVertexZ(-100),
  fHiBin(-1),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fnJets(0),
  fJetRawPtArray(),
  fJetMaxTrackPtArray(),
  fCaloJetFilterBit(0),
  fPrimaryVertexFilterBit(0),
  fBeamScrapingFilterBit(0),
  fCollisionEventSelectionFilterBit(0),
  fHBHENoiseFilterBit(0),
  fTrackPtArray(),
  fTrackPtErrorArray(),
  fTrackPhiArray(),
  fTrackEtaArray(),
  fnTracks(0),
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
ForestReader::ForestReader(Int_t dataType) :
  fDataType(0),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fHiVzBranch(0),
  fHiBinBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetEtaBranch(0),
  fnJetsBranch(0),
  fJetRawPtBranch(0),
  fJetMaxTrackPtBranch(0),
  fCaloJetFilterBranch(0),
  fPrimaryVertexBranch(0),
  fBeamScrapingBranch(0),
  fCollisionEventSelectionBranch(0),
  fHBHENoiseBranch(0),
  fTrackPtBranch(0),
  fTrackPtErrorBranch(0),
  fTrackPhiBranch(0),
  fTrackEtaBranch(0),
  fnTracksBranch(0),
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
  fVertexZ(-100),
  fHiBin(-1),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fnJets(0),
  fJetRawPtArray(),
  fJetMaxTrackPtArray(),
  fCaloJetFilterBit(0),
  fPrimaryVertexFilterBit(0),
  fBeamScrapingFilterBit(0),
  fCollisionEventSelectionFilterBit(0),
  fHBHENoiseFilterBit(0),
  fTrackPtArray(),
  fTrackPtErrorArray(),
  fTrackPhiArray(),
  fTrackEtaArray(),
  fnTracks(0),
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
  
  SetDataType(dataType);
  
}

/*
 * Copy constructor
 */
ForestReader::ForestReader(const ForestReader& in) :
  fDataType(in.fDataType),
  fHeavyIonTree(in.fHeavyIonTree),
  fJetTree(in.fJetTree),
  fHltTree(in.fHltTree),
  fSkimTree(in.fSkimTree),
  fTrackTree(in.fTrackTree),
  fHiVzBranch(in.fHiVzBranch),
  fHiBinBranch(in.fHiBinBranch),
  fJetPtBranch(in.fJetPtBranch),
  fJetPhiBranch(in.fJetPhiBranch),
  fJetEtaBranch(in.fJetEtaBranch),
  fnJetsBranch(in.fnJetsBranch),
  fJetRawPtBranch(in.fJetRawPtBranch),
  fJetMaxTrackPtBranch(in.fJetMaxTrackPtBranch),
  fCaloJetFilterBranch(in.fCaloJetFilterBranch),
  fPrimaryVertexBranch(in.fPrimaryVertexBranch),
  fBeamScrapingBranch(in.fBeamScrapingBranch),
  fCollisionEventSelectionBranch(in.fCollisionEventSelectionBranch),
  fHBHENoiseBranch(in.fHBHENoiseBranch),
  fTrackPtBranch(in.fTrackPtBranch),
  fTrackPtErrorBranch(in.fTrackPtErrorBranch),
  fTrackPhiBranch(in.fTrackPhiBranch),
  fTrackEtaBranch(in.fTrackEtaBranch),
  fnTracksBranch(in.fnTracksBranch),
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
  fVertexZ(in.fVertexZ),
  fHiBin(in.fHiBin),
  fnJets(in.fnJets),
  fCaloJetFilterBit(in.fCaloJetFilterBit),
  fPrimaryVertexFilterBit(in.fPrimaryVertexFilterBit),
  fBeamScrapingFilterBit(in.fBeamScrapingFilterBit),
  fCollisionEventSelectionFilterBit(in.fCollisionEventSelectionFilterBit),
  fHBHENoiseFilterBit(in.fHBHENoiseFilterBit),
  fnTracks(in.fnTracks)
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
ForestReader& ForestReader::operator=(const ForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  fDataType = in.fDataType;
  fHeavyIonTree = in.fHeavyIonTree;
  fJetTree = in.fJetTree;
  fHltTree = in.fHltTree;
  fSkimTree = in.fSkimTree;
  fTrackTree = in.fTrackTree;
  fHiVzBranch = in.fHiVzBranch;
  fHiBinBranch = in.fHiBinBranch;
  fJetPtBranch = in.fJetPtBranch;
  fJetPhiBranch = in.fJetPhiBranch;
  fJetEtaBranch = in.fJetEtaBranch;
  fnJetsBranch = in.fnJetsBranch;
  fJetRawPtBranch = in.fJetRawPtBranch;
  fJetMaxTrackPtBranch = in.fJetMaxTrackPtBranch;
  fCaloJetFilterBranch = in.fCaloJetFilterBranch;
  fPrimaryVertexBranch = in.fPrimaryVertexBranch;
  fBeamScrapingBranch = in.fBeamScrapingBranch;
  fCollisionEventSelectionBranch = in.fCollisionEventSelectionBranch;
  fHBHENoiseBranch = in.fHBHENoiseBranch;
  fTrackPtBranch = in.fTrackPtBranch;
  fTrackPtErrorBranch = in.fTrackPtErrorBranch;
  fTrackPhiBranch = in.fTrackPhiBranch;
  fTrackEtaBranch = in.fTrackEtaBranch;
  fnTracksBranch = in.fnTracksBranch;
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
  fVertexZ = in.fVertexZ;
  fHiBin = in.fHiBin;
  fnJets = in.fnJets;
  fCaloJetFilterBit = in.fCaloJetFilterBit;
  fPrimaryVertexFilterBit = in.fPrimaryVertexFilterBit;
  fBeamScrapingFilterBit = in.fBeamScrapingFilterBit;
  fCollisionEventSelectionFilterBit = in.fCollisionEventSelectionFilterBit;
  fHBHENoiseFilterBit = in.fHBHENoiseFilterBit;
  fnTracks = in.fnTracks;
  
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

/*
 * Initialization, meaning that the branches are connected to the tree
 */
void ForestReader::Initialize(){
  
  // Connect the branches of the heavy ion tree
  fHeavyIonTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  fHeavyIonTree->SetBranchAddress("hiBin",&fHiBin,&fHiBinBranch);
  
  // Connect the branches to the jet tree
  fJetTree->SetBranchAddress("jtpt",&fJetPtArray,&fJetPtBranch);
  fJetTree->SetBranchAddress("jtphi",&fJetPhiArray,&fJetPhiBranch);
  fJetTree->SetBranchAddress("jteta",&fJetEtaArray,&fJetEtaBranch);
  fJetTree->SetBranchAddress("nref",&fnJets,&fnJetsBranch);
  fJetTree->SetBranchAddress("rawpt",&fJetRawPtArray,&fJetRawPtBranch);
  fJetTree->SetBranchAddress("trackMax",&fJetMaxTrackPtArray,&fJetMaxTrackPtBranch);
  
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
  
  // Connect the branches to the HLT tree (only for real data)
  if(fDataType == kPp){ // pp data
    fHltTree->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1",&fCaloJetFilterBit,&fCaloJetFilterBranch);
  } else if (fDataType == kPbPb){ // PbPb data
    fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1",&fCaloJetFilterBit,&fCaloJetFilterBranch);
  } else { // Monte Carlo
    fCaloJetFilterBit = 1;  // No filter for Monte Carlo
  }
  
  // TODO: ADD cuts for phfCoincFilter3 and pclusterCompatibilityFilter for HI collisions!!
  
  // Connect the branches to the skim tree (different for pp and PbPb data, no connection for Monte Carlo)
  if(fDataType == kPp){ // pp data
    fSkimTree->SetBranchAddress("pPAprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fSkimTree->SetBranchAddress("pBeamScrapingFilter",&fBeamScrapingFilterBit,&fBeamScrapingBranch);
    fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fCollisionEventSelectionFilterBit = 1;  // No collision event selection filter for pp
  } else if (fDataType == kPbPb){ // PbPb data
    fSkimTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fSkimTree->SetBranchAddress("pcollisionEventSelection",&fCollisionEventSelectionFilterBit,&fCollisionEventSelectionBranch);
    fBeamScrapingFilterBit = 1;  // No beam scraping filter for PbPb
  } else { // Monte Carlo
    fPrimaryVertexFilterBit = 1;            // No filter for Monte Carlo
    fBeamScrapingFilterBit = 1;             // No filter for Monte Carlo
    fCollisionEventSelectionFilterBit = 1;  // No filter for Monte Carlo
    fHBHENoiseFilterBit = 1;                // No filter for Monte Carlo
  }
  
  // Connect the branches to the track tree
  fTrackTree->SetBranchAddress("trkPt",&fTrackPtArray,&fTrackPtBranch);
  fTrackTree->SetBranchAddress("trkPtError",&fTrackPtErrorArray,&fTrackPtErrorBranch);
  fTrackTree->SetBranchAddress("trkPhi",&fTrackPhiArray,&fTrackPhiBranch);
  fTrackTree->SetBranchAddress("trkEta",&fTrackEtaArray,&fTrackEtaBranch);
  fTrackTree->SetBranchAddress("nTrk",&fnTracks,&fnTracksBranch);
  fTrackTree->SetBranchAddress("highPurity",&fHighPurityTrackArray,&fHighPurityTrackBranch);
  fTrackTree->SetBranchAddress("trkDz1",&fTrackVertexDistanceZArray,&fTrackVertexDistanceZBranch);
  fTrackTree->SetBranchAddress("trkDzError1",&fTrackVertexDistanceZErrorArray,&fTrackVertexDistanceZErrorBranch);
  fTrackTree->SetBranchAddress("trkDxy1",&fTrackVertexDistanceXYArray,&fTrackVertexDistanceXYBranch);
  fTrackTree->SetBranchAddress("trkDxyError1",&fTrackVertexDistanceXYErrorArray,&fTrackVertexDistanceXYErrorBranch);
  fTrackTree->SetBranchAddress("trkChi2",&fTrackChi2Array,&fTrackChi2Branch);
  fTrackTree->SetBranchAddress("trkNdof",&fnTrackDegreesOfFreedomArray,&fnTrackDegreesOfFreedomBranch);
  fTrackTree->SetBranchAddress("trkNlayer",&fnHitsTrackerLayerArray,&fnHitsTrackerLayerBranch);
  fTrackTree->SetBranchAddress("trkNHit",&fnHitsTrackArray,&fnHitsTrackBranch);
  fTrackTree->SetBranchAddress("pfEcal",&fTrackEnergyEcalArray,&fTrackEnergyEcalBranch);
  fTrackTree->SetBranchAddress("pfHcal",&fTrackEnergyHcalArray,&fTrackEnergyHcalBranch);
}

/*
 * Connect a new tree to the reader
 */
void ForestReader::ReadForestFromFile(TFile *inputFile){
  
  // Connect a trees from the file to the reader
  fHeavyIonTree = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
  fHltTree = (TTree*)inputFile->Get("hltanalysis/HltTree");
  fSkimTree = (TTree*)inputFile->Get("skimanalysis/HltTree");
  
  // The jet tree has different name in different datasets
  if(fDataType == kPp || fDataType == kPpMC){
    fJetTree = (TTree*)inputFile->Get("ak4CaloJetAnalyzer/t");
  } else if (fDataType == kPbPb || fDataType == kPbPbMC){
    fJetTree = (TTree*)inputFile->Get("akPu4CaloJetAnalyzer/t");
  } else if (fDataType == kLocalTest){
    fJetTree = (TTree*)inputFile->Get("ak4PFJetAnalyzer/t");
  }
  
  // The track tree has different name for differant datasets
  if(fDataType == kPp || fDataType == kPpMC || fDataType == kLocalTest){
    fTrackTree = (TTree*)inputFile->Get("ppTrack/trackTree");
  } else if (fDataType == kPbPb || fDataType == kPbPbMC){
    fTrackTree = (TTree*)inputFile->Get("anaTrack/trackTree");
  }
  
  Initialize();
}

/*
 * Load an event to memory
 */
void ForestReader::GetEvent(Int_t nEvent) const{
  fHeavyIonTree->GetEntry(nEvent);
  fJetTree->GetEntry(nEvent);
  fHltTree->GetEntry(nEvent);
  fSkimTree->GetEntry(nEvent);
  fTrackTree->GetEntry(nEvent);
}

// Getter for number of events in the tree
Int_t ForestReader::GetNEvents() const{
  return fJetPtBranch->GetEntries();
}

// Getter for jet pT
Float_t ForestReader::GetJetPt(Int_t iJet) const{
  return fJetPtArray[iJet];
}

// Getter for jet phi
Float_t ForestReader::GetJetPhi(Int_t iJet) const{
  return fJetPhiArray[iJet];
}

// Getter for jet eta
Float_t ForestReader::GetJetEta(Int_t iJet) const{
  return fJetEtaArray[iJet];
}

// Getter for number of jets in an event
Int_t ForestReader::GetNJets() const{
  return fnJets;
}

// Getter for jet raw pT
Float_t ForestReader::GetJetRawPt(Int_t iJet) const{
  return fJetRawPtArray[iJet];
}

// Getter for maximum track pT inside a jet
Float_t ForestReader::GetJetMaxTrackPt(Int_t iJet) const{
  return fJetMaxTrackPtArray[iJet];
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

// Getter for HB/HE noisr filter bit. Always 1 for MC (set in the initializer).
Int_t ForestReader::GetHBHENoiseFilterBit() const{
  return fHBHENoiseFilterBit;
}

// Getter for HB/HE noisr filter bit. Always 1 for MC and pp (set in the initializer).
Int_t ForestReader::GetCollisionEventSelectionFilterBit() const{
  return fCollisionEventSelectionFilterBit;
}

// Getter for track pT
Float_t ForestReader::GetTrackPt(Int_t iTrack) const{
  return fTrackPtArray[iTrack];
}

// Getter for track pT error
Float_t ForestReader::GetTrackPtError(Int_t iTrack) const{
  return fTrackPtErrorArray[iTrack];
}

// Getter for track phi
Float_t ForestReader::GetTrackPhi(Int_t iTrack) const{
  return fTrackPhiArray[iTrack];
}

// Getter for track eta
Float_t ForestReader::GetTrackEta(Int_t iTrack) const{
  return fTrackEtaArray[iTrack];
}

// Getter for number of tracks in an event
Int_t ForestReader::GetNTracks() const{
  return fnTracks;
}

// Getter for high purity of the track
bool ForestReader::GetTrackHighPurity(Int_t iTrack) const{
  return fHighPurityTrackArray[iTrack];
}

// Getter for track distance from primary vertex in z-direction
Float_t ForestReader::GetTrackVertexDistanceZ(Int_t iTrack) const{
  return fTrackVertexDistanceZArray[iTrack];
}

// Getter for error of track distance from primary vertex in z-direction
Float_t ForestReader::GetTrackVertexDistanceZError(Int_t iTrack) const{
  return fTrackVertexDistanceZErrorArray[iTrack];
}

// Getter for track distance from primary vertex in xy-direction
Float_t ForestReader::GetTrackVertexDistanceXY(Int_t iTrack) const{
  return fTrackVertexDistanceXYArray[iTrack];
}

// Getter for error of track distance from primary vertex in xy-direction
Float_t ForestReader::GetTrackVertexDistanceXYError(Int_t iTrack) const{
  return fTrackVertexDistanceXYErrorArray[iTrack];
}

// Getter for track chi2 value from reconstruction fit
Float_t ForestReader::GetTrackChi2(Int_t iTrack) const{
  return fTrackChi2Array[iTrack];
}

// Getter for number of degrees of freedom in reconstruction fit
UChar_t ForestReader::GetNTrackDegreesOfFreedom(Int_t iTrack) const{
  return fnTrackDegreesOfFreedomArray[iTrack];
}

// Getter for number of hits in tracker layers
UChar_t ForestReader::GetNHitsTrackerLayer(Int_t iTrack) const{
  return fnHitsTrackerLayerArray[iTrack];
}

// Getter for number of hits for the track
UChar_t ForestReader::GetNHitsTrack(Int_t iTrack) const{
  return fnHitsTrackArray[iTrack];
}

// Getter for track energy in ECal
Float_t ForestReader::GetTrackEnergyEcal(Int_t iTrack) const{
  return fTrackEnergyEcalArray[iTrack];
}

// Getter for track energy in HCal
Float_t ForestReader::GetTrackEnergyHcal(Int_t iTrack) const{
  return fTrackEnergyHcalArray[iTrack];
}
