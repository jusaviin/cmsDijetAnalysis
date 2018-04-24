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
  static const int fnMaxJet = 60;        // Maximum number of jets in an event
  static const int fnMaxTrack = 11000;   // Maximum number of tracks in an event
  
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
  int GetHiBin() const;             // Getter for CMS hiBin
  
  // Getters for leaves in jet tree
  float GetJetPt(int iJet) const;         // Getter for jet pT
  float GetJetPhi(int iJet) const;        // Getter for jet phi
  float GetJetEta(int iJet) const;        // Getter for jet eta
  int GetNJets() const;                   // Getter for number of jets
  float GetJetRawPt(int iJet) const;      // Getter for jet raw pT
  float GetJetMaxTrackPt(int iJet) const; // Getter for maximum track pT inside a jet
  
  // Getters for leaves in HLT tree
  int GetCaloJetFilterBit() const;  // Getter for calorimeter jet filter bit
  
  // Getters for leaves in skim tree
  int GetPrimaryVertexFilterBit() const;           // Getter for primary vertex filter bit
  int GetBeamScrapingFilterBit() const;            // Getter got beam scraping filter bit
  int GetHBHENoiseFilterBit() const;               // Getter for HB/HE noise filter bit
  int GetCollisionEventSelectionFilterBit() const; // Getter for collision event selection filter bit
  
  // Getters for leaves in the track tree
  float GetTrackPt(int iTrack) const;       // Getter for track pT
  float GetTrackPtError(int iTrack) const;  // Getter for track pT error
  float GetTrackPhi(int iTrack) const;      // Getter for track phi
  float GetTrackEta(int iTrack) const;      // Getter for track eta
  int GetNTracks() const;                   // Getter for number of track
  int GetTrackHighPurity(int iTrack) const; // Getter for the high purity of the track
  
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
  TTree *fTrackTree;       // Tree for tracks  PbPb: anaTrack/trackTree pp: ppTrack/trackTree GenParticles: HiGenParticleAna/hi
  
  // Branches for heavy ion tree
  TBranch *fHiVzBranch;            // Branch for vertex z-position
  TBranch *fHiBinBranch;    // Branch for centrality
  
  // Branches for jet tree
  TBranch *fJetPtBranch;         // Branch for jet pT
  TBranch *fJetPhiBranch;        // Branch for jet phi
  TBranch *fJetEtaBranch;        // Branch for jet eta
  TBranch *fnJetsBranch;         // Branch for number of jets in an event
  TBranch *fJetRawPtBranch;      // Branch for raw jet pT
  TBranch *fJetMaxTrackPtBranch; // Maximum pT for a track inside a jet
  
  // Branches for HLT tree
  TBranch *fCaloJetFilterBranch;    // Branch for calo jet filter bit
  
  // Branches for skim tree
  TBranch *fPrimaryVertexBranch;           // Branch for primary vertex filter bit
  TBranch *fBeamScrapingBranch;            // Branch for beam scraping filter bit
  TBranch *fCollisionEventSelectionBranch; // Branch for collision event selection filter bit
  TBranch *fHBHENoiseBranch;               // Branch for HB/HE noise filter bit
  
  // Branches for track tree
  TBranch *fTrackPtBranch;                    // Branch for track pT
  TBranch *fTrackPtErrorBranch;               // Branch for track pT error
  TBranch *fTrackPhiBranch;                   // Branch for track phi
  TBranch *fTrackEtaBranch;                   // Branch for track eta
  TBranch *fnTracksBranch;                    // Branch for number of tracks
  TBranch *fHighPurityTrackBranch;            // Branch for high purity of the track
  TBranch *fTrackVertexDistanceZBranch;       // Branch for track distance from primary vertex in z-direction
  TBranch *fTrackVertexDistanceZErrorBranch;  // Branch for error for track distance from primary vertex in z-direction
  TBranch *fTrackVertexDistanceXYBranch;      // Branch for track distance from primary vertex in xy-direction
  TBranch *fTrackVertexDistanceXYErrorBranch; // Branch for error for track distance from primary vertex in xy-direction
  TBranch *fTrackChi2Branch;                  // Branch for track chi2 value from reconstruction fit
  TBranch *fnTrackDegreesOfFreedomBranch;     // Branch for number of degrees of freedom in reconstruction fit
  TBranch *fnHitsTrackerLayerBranch;          // Branch for number of hits in tracker layers
  
  // Leaves for heavy ion tree
  float fVertexZ;    // Vertex z-position
  int fHiBin;        // HiBin = Centrality percentile * 2
  
  // Leaves for jet tree
  float fJetPtArray[fnMaxJet] = {0};         // pT:s of all the jets in an event
  float fJetPhiArray[fnMaxJet] = {0};        // phis of all the jets in an event
  float fJetEtaArray[fnMaxJet] = {0};        // etas of all the jets in an event
  int fnJets;                                // number of jets in an event
  float fJetRawPtArray[fnMaxJet] = {0};      // raw jet pT for all the jets in an event
  float fJetMaxTrackPtArray[fnMaxJet] = {0}; // maximum track pT inside a jet for all the jets in an event
  
  // Leaves for the HLT tree
  int fCaloJetFilterBit;    // Filter bit for calorimeter jets
  
  // Leaves for the skim tree
  int fPrimaryVertexFilterBit;           // Filter bit for primary vertex
  int fBeamScrapingFilterBit;            // Filter bit for beam scraping
  int fCollisionEventSelectionFilterBit; // Filter bit for collision event selection
  int fHBHENoiseFilterBit;               // Filter bit for HB/HE noise
  
  // Leaves for the track tree
  float fTrackPtArray[fnMaxTrack] = {0};                    // Array for track pT:s
  float fTrackPtErrorArray[fnMaxTrack] = {0};               // Array for track pT errors
  float fTrackPhiArray[fnMaxTrack] = {0};                   // Array for track phis
  float fTrackEtaArray[fnMaxTrack] = {0};                   // Array for track etas
  int fnTracks;                                             // Number of tracks
  int fHighPurityTrackArray[fnMaxTrack] = {0};              // Array for the high purity of tracks
  float fTrackVertexDistanceZArray[fnMaxTrack] = {0};       // Array for track distance from primary vertex in z-direction
  float fTrackVertexDistanceZErrorArray[fnMaxTrack] = {0};  // Array for error for track distance from primary vertex in z-direction
  float fTrackVertexDistanceXYArray[fnMaxTrack] = {0};      // Array for track distance from primary vertex in xy-direction
  float fTrackVertexDistanceXYErrorArray[fnMaxTrack] = {0}; // Array for error for track distance from primary vertex in xy-direction
  float fTrackChi2Array[fnMaxTrack] = {0};                  // Array for track chi2 value from reconstruction fit
  int fnTrackDegreesOfFreedomArray[fnMaxTrack] = {0};       // Array for number of degrees of freedom in reconstruction fit
  int fnHitsTrackerLayerArray[fnMaxTrack] = {0};            // Array for number of hits in tracker layers
  
};

#endif






















