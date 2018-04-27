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
  static const Int_t fnMaxJet = 60;        // Maximum number of jets in an event
  static const Int_t fnMaxTrack = 15000;   // Maximum number of tracks in an event
  
public:
  
  // Possible data types to be read with the reader class
  enum enumDataTypes{kPp, kPbPb, kPpMC, kPbPbMC, kLocalTest, knDataTypes};
  
  // Constructors and destructors
  ForestReader();                                   // Default constructor
  ForestReader(Int_t dataType);                       // Custom constructor
  ForestReader(const ForestReader& in);             // Copy constructor
  virtual ~ForestReader();                          // Destructor
  ForestReader& operator=(const ForestReader& obj); // Equal sign operator
  
  // Methods
  void GetEvent(Int_t nEvent) const;             // Get the nth event in tree
  Int_t GetNEvents() const;                      // Get the number of events
  void ReadForestFromFile(TFile *inputFile);   // Read the forest from a file
  
  // Getters for leaves in heavy ion tree
  Float_t GetVz() const;              // Getter for vertex z position
  Float_t GetCentrality() const;      // Getter for centrality
  Int_t GetHiBin() const;             // Getter for CMS hiBin
  
  // Getters for leaves in jet tree
  Float_t GetJetPt(Int_t iJet) const;         // Getter for jet pT
  Float_t GetJetPhi(Int_t iJet) const;        // Getter for jet phi
  Float_t GetJetEta(Int_t iJet) const;        // Getter for jet eta
  Int_t GetNJets() const;                     // Getter for number of jets
  Float_t GetJetRawPt(Int_t iJet) const;      // Getter for jet raw pT
  Float_t GetJetMaxTrackPt(Int_t iJet) const; // Getter for maximum track pT inside a jet
  
  // Getters for leaves in HLT tree
  Int_t GetCaloJetFilterBit() const;  // Getter for calorimeter jet filter bit
  
  // Getters for leaves in skim tree
  Int_t GetPrimaryVertexFilterBit() const;           // Getter for primary vertex filter bit
  Int_t GetBeamScrapingFilterBit() const;            // Getter got beam scraping filter bit
  Int_t GetHBHENoiseFilterBit() const;               // Getter for HB/HE noise filter bit
  Int_t GetCollisionEventSelectionFilterBit() const; // Getter for collision event selection filter bit
  
  // Getters for leaves in the track tree
  Float_t GetTrackPt(Int_t iTrack) const;                    // Getter for track pT
  Float_t GetTrackPtError(Int_t iTrack) const;               // Getter for track pT error
  Float_t GetTrackPhi(Int_t iTrack) const;                   // Getter for track phi
  Float_t GetTrackEta(Int_t iTrack) const;                   // Getter for track eta
  Int_t GetNTracks() const;                                  // Getter for number of track
  Bool_t GetTrackHighPurity(Int_t iTrack) const;             // Getter for the high purity of the track
  Float_t GetTrackVertexDistanceZ(Int_t iTrack) const;       // Getter for track distance from primary vertex in z-direction
  Float_t GetTrackVertexDistanceZError(Int_t iTrack) const;  // Getter for error of track distance from primary vertex in z-direction
  Float_t GetTrackVertexDistanceXY(Int_t iTrack) const;      // Getter for track distance from primary vertex in xy-direction
  Float_t GetTrackVertexDistanceXYError(Int_t iTrack) const; // Getter for error of track distance from primary vertex in xy-direction
  Float_t GetTrackChi2(Int_t iTrack) const;                  // Getter for track chi2 value from reconstruction fit
  UChar_t GetNTrackDegreesOfFreedom(Int_t iTrack) const;     // Getter for number of degrees of freedom in reconstruction fit
  UChar_t GetNHitsTrackerLayer(Int_t iTrack) const;          // Getter for number of hits in tracker layers
  UChar_t GetNHitsTrack(Int_t iTrack) const;                 // Getter for number of hits for the track
  Float_t GetTrackEnergyEcal(Int_t iTrack) const;            // Getter for track energy in ECal
  Float_t GetTrackEnergyHcal(Int_t iTrack) const;            // Getter for track energy in HCal
  
  // Setter for data type
  void SetDataType(Int_t dataType); // Setter for data type
  
private:
  
  // Methods
  void Initialize();  // Connect the branches to the tree
  
  Int_t fDataType;  // Type of data read with the tree. 0 = pp, 1 = PbPb, 2 = ppMC, 3 = PbPbMC, 4 = LocalTest
  
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
  TBranch *fnHitsTrackBranch;                 // Branch for number of hits for the track
  TBranch *fTrackEnergyEcalBranch;            // Branch for track energy in ECal
  TBranch *fTrackEnergyHcalBranch;            // Branch for track energy in HCal
  
  // TODO: Change all leaf types to root defined ones to avoid losing precision!!!!
  
  // Leaves for heavy ion tree
  Float_t fVertexZ;    // Vertex z-position
  Int_t fHiBin;        // HiBin = Centrality percentile * 2
  
  // Leaves for jet tree
  Float_t fJetPtArray[fnMaxJet] = {0};         // pT:s of all the jets in an event
  Float_t fJetPhiArray[fnMaxJet] = {0};        // phis of all the jets in an event
  Float_t fJetEtaArray[fnMaxJet] = {0};        // etas of all the jets in an event
  Int_t fnJets;                                // number of jets in an event
  Float_t fJetRawPtArray[fnMaxJet] = {0};      // raw jet pT for all the jets in an event
  Float_t fJetMaxTrackPtArray[fnMaxJet] = {0}; // maximum track pT inside a jet for all the jets in an event
  
  // Leaves for the HLT tree
  Int_t fCaloJetFilterBit;    // Filter bit for calorimeter jets
  
  // Leaves for the skim tree
  Int_t fPrimaryVertexFilterBit;           // Filter bit for primary vertex
  Int_t fBeamScrapingFilterBit;            // Filter bit for beam scraping
  Int_t fCollisionEventSelectionFilterBit; // Filter bit for collision event selection
  Int_t fHBHENoiseFilterBit;               // Filter bit for HB/HE noise
  
  // Leaves for the track tree
  Float_t fTrackPtArray[fnMaxTrack] = {0};                    // Array for track pT:s
  Float_t fTrackPtErrorArray[fnMaxTrack] = {0};               // Array for track pT errors
  Float_t fTrackPhiArray[fnMaxTrack] = {0};                   // Array for track phis
  Float_t fTrackEtaArray[fnMaxTrack] = {0};                   // Array for track etas
  Int_t fnTracks;                                             // Number of tracks
  Bool_t fHighPurityTrackArray[fnMaxTrack] = {0};             // Array for the high purity of tracks
  Float_t fTrackVertexDistanceZArray[fnMaxTrack] = {0};       // Array for track distance from primary vertex in z-direction
  Float_t fTrackVertexDistanceZErrorArray[fnMaxTrack] = {0};  // Array for error for track distance from primary vertex in z-direction
  Float_t fTrackVertexDistanceXYArray[fnMaxTrack] = {0};      // Array for track distance from primary vertex in xy-direction
  Float_t fTrackVertexDistanceXYErrorArray[fnMaxTrack] = {0}; // Array for error for track distance from primary vertex in xy-direction
  Float_t fTrackChi2Array[fnMaxTrack] = {0};                  // Array for track chi2 value from reconstruction fit
  UChar_t fnTrackDegreesOfFreedomArray[fnMaxTrack] = {0};     // Array for number of degrees of freedom in reconstruction fit
  UChar_t fnHitsTrackerLayerArray[fnMaxTrack] = {0};          // Array for number of hits in tracker layers
  UChar_t fnHitsTrackArray[fnMaxTrack] = {0};                 // Array for number of hits for the track
  Float_t fTrackEnergyEcalArray[fnMaxTrack] = {0};            // Array for track energy in ECal
  Float_t fTrackEnergyHcalArray[fnMaxTrack] = {0};            // Array for track energy in HCal
  
};

#endif






















