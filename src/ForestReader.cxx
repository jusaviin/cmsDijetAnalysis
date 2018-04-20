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
  fHiVzBranch(0),
  fHiCentralityBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetEtaBranch(0),
  fnJetsBranch(0),
  fCaloJetFilterBranch(0),
  fPrimaryVertexBranch(0),
  fBeamScrapingBranch(0),
  fCollisionEventSelectionBranch(0),
  fHBHENoiseBranch(0),
  fVertexZ(-100),
  fCentrality(-1),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fnJets(0),
  fCaloJetFilterBit(0),
  fPrimaryVertexFilterBit(0),
  fBeamScrapingFilterBit(0),
  fCollisionEventSelectionFilterBit(0),
  fHBHENoiseFilterBit(0)
{
  // Default constructor
}

/*
 * Custom constructor
 *
 *  Arguments:
 *   int dataType: 0 = pp, 1 = PbPb, 2 = pp MC, 3 = PbPb MC, 4 = Local Test
 */
ForestReader::ForestReader(int dataType) :
  fDataType(0),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fHiVzBranch(0),
  fHiCentralityBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetEtaBranch(0),
  fnJetsBranch(0),
  fCaloJetFilterBranch(0),
  fPrimaryVertexBranch(0),
  fBeamScrapingBranch(0),
  fCollisionEventSelectionBranch(0),
  fHBHENoiseBranch(0),
  fVertexZ(-100),
  fCentrality(-1),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fnJets(0),
  fCaloJetFilterBit(0),
  fPrimaryVertexFilterBit(0),
  fBeamScrapingFilterBit(0),
  fCollisionEventSelectionFilterBit(0),
  fHBHENoiseFilterBit(0)
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
  fHiVzBranch(in.fHiVzBranch),
  fHiCentralityBranch(in.fHiCentralityBranch),
  fJetPtBranch(in.fJetPtBranch),
  fJetPhiBranch(in.fJetPhiBranch),
  fJetEtaBranch(in.fJetEtaBranch),
  fnJetsBranch(in.fnJetsBranch),
  fCaloJetFilterBranch(in.fCaloJetFilterBranch),
  fPrimaryVertexBranch(in.fPrimaryVertexBranch),
  fBeamScrapingBranch(in.fBeamScrapingBranch),
  fCollisionEventSelectionBranch(in.fCollisionEventSelectionBranch),
  fHBHENoiseBranch(in.fHBHENoiseBranch),
  fVertexZ(in.fVertexZ),
  fCentrality(in.fCentrality),
  fnJets(in.fnJets),
  fCaloJetFilterBit(in.fCaloJetFilterBit),
  fPrimaryVertexFilterBit(in.fPrimaryVertexFilterBit),
  fBeamScrapingFilterBit(in.fBeamScrapingFilterBit),
  fCollisionEventSelectionFilterBit(in.fCollisionEventSelectionFilterBit),
  fHBHENoiseFilterBit(in.fHBHENoiseFilterBit)
{
  // Copy constructor
  for(int i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
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
  fHiVzBranch = in.fHiVzBranch;
  fHiCentralityBranch = in.fHiCentralityBranch;
  fJetPtBranch = in.fJetPtBranch;
  fJetPhiBranch = in.fJetPhiBranch;
  fJetEtaBranch = in.fJetEtaBranch;
  fnJetsBranch = in.fnJetsBranch;
  fCaloJetFilterBranch = in.fCaloJetFilterBranch;
  fPrimaryVertexBranch = in.fPrimaryVertexBranch;
  fBeamScrapingBranch = in.fBeamScrapingBranch;
  fCollisionEventSelectionBranch = in.fCollisionEventSelectionBranch;
  fHBHENoiseBranch = in.fHBHENoiseBranch;
  fVertexZ = in.fVertexZ;
  fCentrality = in.fCentrality;
  fnJets = in.fnJets;
  fCaloJetFilterBit = in.fCaloJetFilterBit;
  fPrimaryVertexFilterBit = in.fPrimaryVertexFilterBit;
  fBeamScrapingFilterBit = in.fBeamScrapingFilterBit;
  fCollisionEventSelectionFilterBit = in.fCollisionEventSelectionFilterBit;
  fHBHENoiseFilterBit = in.fHBHENoiseFilterBit;
  for(int i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
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
void ForestReader::SetDataType(int dataType){
  
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
  fHeavyIonTree->SetBranchAddress("hiBin",&fCentrality,&fHiCentralityBranch);
  
  // Connect the branches to the jet tree
  fJetTree->SetBranchAddress("jtpt",&fJetPtArray,&fJetPtBranch);
  fJetTree->SetBranchAddress("jtphi",&fJetPhiArray,&fJetPhiBranch);
  fJetTree->SetBranchAddress("jteta",&fJetEtaArray,&fJetEtaBranch);
  fJetTree->SetBranchAddress("nref",&fnJets,&fnJetsBranch);
  
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
  
  Initialize();
}

/*
 * Load an event to memory
 */
void ForestReader::GetEvent(int nEvent) const{
  fHeavyIonTree->GetEntry(nEvent);
  fJetTree->GetEntry(nEvent);
  fHltTree->GetEntry(nEvent);
  fSkimTree->GetEntry(nEvent);
}

// Getter for number of events in the tree
int ForestReader::GetNEvents() const{
  return fJetPtBranch->GetEntries();
}

// Getter for jet pT
float ForestReader::GetJetPt(int iJet) const{
  return fJetPtArray[iJet];
}

// Getter for jet phi
float ForestReader::GetJetPhi(int iJet) const{
  return fJetPhiArray[iJet];
}

// Getter for jet eta
float ForestReader::GetJetEta(int iJet) const{
  return fJetEtaArray[iJet];
}

// Getter for number of jets in an event
int ForestReader::GetNJets() const{
  return fnJets;
}

// Getter for vertex z position
float ForestReader::GetVz() const{
  return fVertexZ;
}

// Getter for centrality. For some reason CMS normalizes centrality from 0 to 200, thus division by 2.
float ForestReader::GetCentrality() const{
  return fCentrality/2.0;
}

// Getter for calorimeter jet filter bit. Always 1 for MC (set in the initializer).
int ForestReader::GetCaloJetFilterBit() const{
  return fCaloJetFilterBit;
}

// Getter for primary vertex filter bit. Always 1 for MC (set in the initializer).
int ForestReader::GetPrimaryVertexFilterBit() const{
  return fPrimaryVertexFilterBit;
}

// Getter for beam scraping filter bit. Always 1 for MC and PbPb (set in the initializer).
int ForestReader::GetBeamScrapingFilterBit() const{
  return fBeamScrapingFilterBit;
}

// Getter for HB/HE noisr filter bit. Always 1 for MC (set in the initializer).
int ForestReader::GetHBHENoiseFilterBit() const{
  return fHBHENoiseFilterBit;
}

// Getter for HB/HE noisr filter bit. Always 1 for MC and pp (set in the initializer).
int ForestReader::GetCollisionEventSelectionFilterBit() const{
  return fCollisionEventSelectionFilterBit;
}
