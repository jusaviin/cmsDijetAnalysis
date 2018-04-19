// Reader for jet trees from CMS data
//
//===========================================================
// ForestReader.h
//
// Author: Jussi Viinikainen
//===========================================================

#ifndef FORESTREADER_H
#define FORESTREADER_H

// C++ includes
#include <iostream>
#include <assert.h>

// Root includes
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>

using namespace std;

class ForestReader{
  
private:
  static const int fnMaxJet = 60;
  
public:
  
  // Possible data types to be read with the reader class
  enum enumDataTypes{kPp, kPbPb, kPpMC, kPbPbMC, kLocalTest, knDataTypes};
  
  // Constructors and destructors
  ForestReader();                                   // Default constructor
  ForestReader(int dataType);                       // Custom constructor
  ForestReader(const ForestReader& in);             // Copy constructor
  virtual ~ForestReader();                          // Destructor
  ForestReader& operator=(const ForestReader& obj); // Equal sign operator
  
  // Methods
  void GetEvent(int nEvent) const;             // Get the nth event in tree
  int GetNEvents() const;                      // Get the number of events
  void ReadForestFromFile(TFile *inputFile);   // Read the forest from a file
  
  // Getters for leaves in heavy ion tree
  float GetVz() const;              // Getter for vertex z position
  float GetCentrality() const;      // Getter for centrality
  
  // Getters for leaves in jet tree
  float GetJetPt(int iJet) const;   // Getter for jet pT
  float GetJetPhi(int iJet) const;  // Getter for jet phi
  float GetJetEta(int iJet) const;  // Getter for jet eta
  int GetNJets() const;             // Getter for number of jets
  
  // Getters for leaves in HLT tree
  int GetCaloJetFilterBit() const;  // Getter for calorimeter jet filter bit
  
  // Getters for leaves in skim tree
  int GetPrimaryVertexFilterBit() const;           // Getter for primary vertex filter bit
  int GetBeamScrapingFilterBit() const;            // Getter got beam scraping filter bit
  int GetHBHENoiseFilterBit() const;               // Getter for HB/HE noise filter bit
  int GetCollisionEventSelectionFilterBit() const; // Getter for collision event selection filter bit
  
  // Setter for data type
  void SetDataType(int dataType); // Setter for data type
  
private:
  
  // Methods
  void Initialize();  // Connect the branches to the tree
  
  int fDataType;  // Type of data read with the tree. 0 = pp, 1 = PbPb, 2 = ppMC, 3 = PbPbMC, 4 = LocalTest
  
  // Trees in the forest
  TTree *fHeavyIonTree;    // Tree for heavy ion event information
  TTree *fJetTree;         // Tree for jet information
  TTree *fHltTree;         // Tree for high level trigger information
  TTree *fSkimTree;        // Tree for event selection information
  
  // Branches for heavy ion tree
  TBranch *fHiVzBranch;            // Branch for vertex z-position
  TBranch *fHiCentralityBranch;    // Branch for centrality
  
  // Branches for jet tree
  TBranch *fJetPtBranch;    // Branch for jet pT
  TBranch *fJetPhiBranch;   // Branch for jet phi
  TBranch *fJetEtaBranch;   // Branch for jet eta
  TBranch *fnJetsBranch;    // Branch for number of jets in an event
  
  // Branches for HLT tree
  TBranch *fCaloJetFilterBranch;    // Branch for calo jet filter bit
  
  // Branches for Skim Tree
  TBranch *fPrimaryVertexBranch;           // Branch for primary vertex filter bit
  TBranch *fBeamScrapingBranch;            // Branch for beam scraping filter bit
  TBranch *fCollisionEventSelectionBranch; // Branch for collision event selection filter bit
  TBranch *fHBHENoiseBranch;               // Branch for HB/HE noise filter bit
  
  // Leaves for heavy ion tree
  float fVertexZ;     // Vertex z-position
  int fCentrality;    // Centrality bin
  
  // Leaves for jet tree
  float fJetPtArray[fnMaxJet] = {0};     // pT:s of all the jets in an event
  float fJetPhiArray[fnMaxJet] = {0};    // phis of all the jets in an event
  float fJetEtaArray[fnMaxJet] = {0};    // etas of all the jets in an event
  int fnJets;                            // number of jets in an event
  
  // Leaves for the HLT tree
  int fCaloJetFilterBit;    // Filter bit for calorimeter jets
  
  // Leaves for the skim tree
  int fPrimaryVertexFilterBit;           // Filter bit for primary vertex
  int fBeamScrapingFilterBit;            // Filter bit for beam scraping
  int fCollisionEventSelectionFilterBit; // Filter bit for collision event selection
  int fHBHENoiseFilterBit;               // Filter bit for HB/HE noise
  
};

#endif






















