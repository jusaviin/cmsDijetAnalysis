// Reader for jet trees from CMS data
//
//===========================================================
// HighForestReader.h
//
// Author: Jussi Viinikainen
//===========================================================

#ifndef HIGHFORESTREADER_H
#define HIGHFORESTREADER_H

// Own includes
#include "ForestReader.h"

using namespace std;

class HighForestReader : public ForestReader{
  
private:
  static const Int_t fnMaxJet = 60;        // Maximum number of jets in an event
  static const Int_t fnMaxTrack = 60000;   // Maximum number of tracks in an event
  
public:
  
  // Constructors and destructors
  HighForestReader();                                   // Default constructor
  HighForestReader(Int_t dataType);                       // Custom constructor
  HighForestReader(const HighForestReader& in);             // Copy constructor
  virtual ~HighForestReader();                          // Destructor
  HighForestReader& operator=(const HighForestReader& obj); // Equal sign operator
  
  // Methods
  void ReadForestFromFile(TFile *inputFile);   // Read the forest from a file
  void BurnForest();                           // Burn the forest  
  void GetEvent(Int_t nEvent);                 // Get the nEventh event from the file
  
  // Getters for leaves in jet tree
  Float_t GetJetPt(Int_t iJet) const;         // Getter for jet pT
  Float_t GetJetPhi(Int_t iJet) const;        // Getter for jet phi
  Float_t GetJetEta(Int_t iJet) const;        // Getter for jet eta
  Float_t GetJetRawPt(Int_t iJet) const;      // Getter for jet raw pT
  Float_t GetJetMaxTrackPt(Int_t iJet) const; // Getter for maximum track pT inside a jet
  
  // Getters for leaves in the track tree
  Float_t GetTrackPt(Int_t iTrack) const;                    // Getter for track pT
  Float_t GetTrackPtError(Int_t iTrack) const;               // Getter for track pT error
  Float_t GetTrackPhi(Int_t iTrack) const;                   // Getter for track phi
  Float_t GetTrackEta(Int_t iTrack) const;                   // Getter for track eta
  Bool_t GetTrackHighPurity(Int_t iTrack) const;             // Getter for the high purity of the track
  Float_t GetTrackVertexDistanceZ(Int_t iTrack) const;       // Getter for track distance from primary vertex in z-direction
  Float_t GetTrackVertexDistanceZError(Int_t iTrack) const;  // Getter for error of track distance from primary vertex in z-direction
  Float_t GetTrackVertexDistanceXY(Int_t iTrack) const;      // Getter for track distance from primary vertex in xy-direction
  Float_t GetTrackVertexDistanceXYError(Int_t iTrack) const; // Getter for error of track distance from primary vertex in xy-direction
  Float_t GetTrackChi2(Int_t iTrack) const;                  // Getter for track chi2 value from reconstruction fit
  Int_t GetNTrackDegreesOfFreedom(Int_t iTrack) const;     // Getter for number of degrees of freedom in reconstruction fit
  Int_t GetNHitsTrackerLayer(Int_t iTrack) const;          // Getter for number of hits in tracker layers
  Int_t GetNHitsTrack(Int_t iTrack) const;                 // Getter for number of hits for the track
  Float_t GetTrackEnergyEcal(Int_t iTrack) const;            // Getter for track energy in ECal
  Float_t GetTrackEnergyHcal(Int_t iTrack) const;            // Getter for track energy in HCal
  
private:
  
  // Methods
  void Initialize();   // Connect the branches to the tree
  
  // Trees in the forest
  TTree *fHeavyIonTree;    // Tree for heavy ion event information
  TTree *fJetTree;         // Tree for jet information
  TTree *fHltTree;         // Tree for high level trigger information
  TTree *fSkimTree;        // Tree for event selection information
  TTree *fTrackTree;       // Tree for tracks  PbPb: anaTrack/trackTree pp: ppTrack/trackTree GenParticles: HiGenParticleAna/hi
  TTree *fParticleFlowCandidateTree;  // Tree for particle flow candidates
  
  // Non-common branches for all types of trees
  TBranch *fnJetsBranch;         // Branch for number of jets in an event
  TBranch *fnTracksBranch;       // Branch for number of tracks
  
  // Leaves for jet tree
  Float_t fJetPtArray[fnMaxJet] = {0};         // pT:s of all the jets in an event
  Float_t fJetPhiArray[fnMaxJet] = {0};        // phis of all the jets in an event
  Float_t fJetEtaArray[fnMaxJet] = {0};        // etas of all the jets in an event
  Float_t fJetRawPtArray[fnMaxJet] = {0};      // raw jet pT for all the jets in an event
  Float_t fJetMaxTrackPtArray[fnMaxJet] = {0}; // maximum track pT inside a jet for all the jets in an event
  
  // Leaves for the track tree
  Float_t fTrackPtArray[fnMaxTrack] = {0};                    // Array for track pT:s
  Float_t fTrackPtErrorArray[fnMaxTrack] = {0};               // Array for track pT errors
  Float_t fTrackPhiArray[fnMaxTrack] = {0};                   // Array for track phis
  Float_t fTrackEtaArray[fnMaxTrack] = {0};                   // Array for track etas
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






















