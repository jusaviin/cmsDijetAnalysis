// Implementation for ForestReader

#include "ForestReader.h"

/*
 * Default constructor
 */
ForestReader::ForestReader() :
  fHeavyIonTree(0),
  fJetTree(0),
  fHiVzBranch(0),
  fHiCentralityBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetEtaBranch(0),
  fnJetsBranch(0),
  fVertexZ(-100),
  fCentrality(-1),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fnJets(0)
{
  // Default constructor
}

/*
 * Copy constructor
 */
ForestReader::ForestReader(const ForestReader& in) :
  fHeavyIonTree(in.fHeavyIonTree),
  fJetTree(in.fJetTree),
  fHiVzBranch(in.fHiVzBranch),
  fHiCentralityBranch(in.fHiCentralityBranch),
  fJetPtBranch(in.fJetPtBranch),
  fJetPhiBranch(in.fJetPhiBranch),
  fJetEtaBranch(in.fJetEtaBranch),
  fnJetsBranch(in.fnJetsBranch),
  fVertexZ(in.fVertexZ),
  fCentrality(in.fCentrality),
  fnJets(in.fnJets)
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
  
  fHeavyIonTree = in.fHeavyIonTree;
  fJetTree = in.fJetTree;
  fHiVzBranch = in.fHiVzBranch;
  fHiCentralityBranch = in.fHiCentralityBranch;
  fJetPtBranch = in.fJetPtBranch;
  fJetPhiBranch = in.fJetPhiBranch;
  fJetEtaBranch = in.fJetEtaBranch;
  fnJetsBranch = in.fnJetsBranch;
  fVertexZ = in.fVertexZ;
  fCentrality = in.fCentrality;
  fnJets = in.fnJets;
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
}

/*
 * Connect a new tree to the reader
 */
void ForestReader::ReadForestFromFile(TFile *inputFile){
  // Connect a trees from the file to the reader
  fHeavyIonTree = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
  fJetTree = (TTree*)inputFile->Get("ak4PFJetAnalyzer/t");
  Initialize();
}

/*
 * Load an event to memory
 */
void ForestReader::GetEvent(int nEvent) const{
  fHeavyIonTree->GetEntry(nEvent);
  fJetTree->GetEntry(nEvent);
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
