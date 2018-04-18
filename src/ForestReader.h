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
  
  // Constructors and destructors
  ForestReader(); // Default constructor
  ForestReader(const ForestReader& in); // Copy constructor
  virtual ~ForestReader(); // Destructor
  ForestReader& operator=(const ForestReader& obj); // Equal sign operator
  
  // Methods
  void GetEvent(int nEvent) const;             // Get the nth event in tree
  int GetNEvents() const;                      // Get the number of events
  void ReadForestFromFile(TFile *inputFile);   // Read the forest from a file
  
  // Getters for leaves for heavy ion tree
  float GetVz() const;              // Getter for vertex z position
  float GetCentrality() const;      // Getter for centrality
  
  // Getters for leaves for jet tree
  float GetJetPt(int iJet) const;   // Getter for jet pT
  float GetJetPhi(int iJet) const;  // Getter for jet phi
  float GetJetEta(int iJet) const;  // Getter for jet eta
  int GetNJets() const;             // Getter for number of jets
  
private:
  
  // Methods
  void Initialize();  // Connect the branches to the tree
  
  // Trees in the forest
  TTree *fHeavyIonTree;
  TTree *fJetTree;
  
  // Branches for heavy ion tree
  TBranch *fHiVzBranch;
  TBranch *fHiCentralityBranch;
  
  // Branches for jet tree
  TBranch *fJetPtBranch;
  TBranch *fJetPhiBranch;
  TBranch *fJetEtaBranch;
  TBranch *fnJetsBranch;
  
  // Leaves for heavy ion tree
  float fVertexZ;
  int fCentrality;
  
  // Leaves for jet tree
  float fJetPtArray[fnMaxJet] = {0};
  float fJetPhiArray[fnMaxJet] = {0};
  float fJetEtaArray[fnMaxJet] = {0};
  int fnJets;
  
};

#endif






















