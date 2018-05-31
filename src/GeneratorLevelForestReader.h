// Reader for jet trees from CMS data
//
//===========================================================
// GeneratorLevelForestReader.h
//
// Author: Jussi Viinikainen
//===========================================================

#ifndef GENERATORLEVELFORESTREADER_H
#define GENERATORLEVELFORESTREADER_H

// Own includes
#include "ForestReader.h"

using namespace std;

class GeneratorLevelForestReader : public ForestReader{
  
private:
  static const Int_t fnMaxJet = 60;        // Maximum number of jets in an event
  
public:
  
  // Constructors and destructors
  GeneratorLevelForestReader();                                       // Default constructor
  GeneratorLevelForestReader(Int_t dataType);                         // Custom constructor
  GeneratorLevelForestReader(const GeneratorLevelForestReader& in);             // Copy constructor
  virtual ~GeneratorLevelForestReader();                              // Destructor
  GeneratorLevelForestReader& operator=(const GeneratorLevelForestReader& obj); // Equal sign operator
  
  // Methods
  void ReadForestFromFile(TFile *inputFile);   // Read the forest from a file
  void BurnForest();                           // Burn the forest  
  void GetEvent(Int_t nEvent);                 // Get the nEventh event from the file
  
  // Getters for leaves in jet tree
  Float_t GetJetPt(Int_t iJet) const;         // Getter for jet pT
  Float_t GetJetPhi(Int_t iJet) const;        // Getter for jet phi
  Float_t GetJetEta(Int_t iJet) const;        // Getter for jet eta
  Float_t GetJetRawPt(Int_t iJet) const;      // Getter for jet raw pT (not relevant for generator jets, just return value that passes cuts)
  Float_t GetJetMaxTrackPt(Int_t iJet) const; // Getter for maximum track pT inside a jet (not relevant for generator jets, just return value that passes cuts)
  
  // Getters for leaves in the track tree relevant for generator level track
  Float_t GetTrackPt(Int_t iTrack) const;          // Getter for track pT
  Float_t GetTrackPhi(Int_t iTrack) const;         // Getter for track phi
  Float_t GetTrackEta(Int_t iTrack) const;         // Getter for track eta
  Int_t GetTrackCharge(Int_t iTrack) const;        // Getter for track charge
  Int_t GetTrackSubevent(Int_t iTrack) const;      // Getter for track subevent index
  
  // Getters for leaves in the track tree that are just there to provide values that pass cuts
  Float_t GetTrackPtError(Int_t iTrack) const;               // Getter for track pT error
  Bool_t GetTrackHighPurity(Int_t iTrack) const;             // Getter for the high purity of the track
  Float_t GetTrackVertexDistanceZ(Int_t iTrack) const;       // Getter for track distance from primary vertex in z-direction
  Float_t GetTrackVertexDistanceZError(Int_t iTrack) const;  // Getter for error of track distance from primary vertex in z-direction
  Float_t GetTrackVertexDistanceXY(Int_t iTrack) const;      // Getter for track distance from primary vertex in xy-direction
  Float_t GetTrackVertexDistanceXYError(Int_t iTrack) const; // Getter for error of track distance from primary vertex in xy-direction
  Float_t GetTrackChi2(Int_t iTrack) const;                  // Getter for track chi2 value from reconstruction fit
  Int_t GetNTrackDegreesOfFreedom(Int_t iTrack) const;       // Getter for number of degrees of freedom in reconstruction fit
  Int_t GetNHitsTrackerLayer(Int_t iTrack) const;            // Getter for number of hits in tracker layers
  Int_t GetNHitsTrack(Int_t iTrack) const;                   // Getter for number of hits for the track
  Float_t GetTrackEnergyEcal(Int_t iTrack) const;            // Getter for track energy in ECal
  Float_t GetTrackEnergyHcal(Int_t iTrack) const;            // Getter for track energy in HCal
  
  // Getters for leaves in the particle flow candidate tree
  Int_t GetParticleFlowCandidateId(Int_t iCandidate) const;      // Getter for particle flow candidate ID
  Float_t GetParticleFlowCandidatePt(Int_t iCandidate) const;    // Getter for particle flow candidate pT
  Float_t GetParticleFlowCandidatePhi(Int_t iCandidate) const;   // Getter for particle flow candidate phi
  Float_t GetParticleFlowCandidateEta(Int_t iCandidate) const;   // Getter for particle flow candidate eta
  Int_t GetNParticleFlowCandidates() const;                      // Getter for number of particle flow candidates in an event
  
private:
  
  // Methods
  void Initialize();       // Connect the branches to the tree

  // Trees in the forest
  TTree *fHeavyIonTree;    // Tree for heavy ion event information
  TTree *fJetTree;         // Tree for jet information
  TTree *fHltTree;         // Tree for HLT information
  TTree *fSkimTree;        // Tree for event selection information
  TTree *fTrackTree;       // Tree for tracks  PbPb: anaTrack/trackTree pp: ppTrack/trackTree GenParticles: HiGenParticleAna/hi
  
  // Leaves for jet tree
  Float_t fJetPtArray[fnMaxJet] = {0};         // pT:s of all the jets in an event
  Float_t fJetPhiArray[fnMaxJet] = {0};        // phis of all the jets in an event
  Float_t fJetEtaArray[fnMaxJet] = {0};        // etas of all the jets in an event
  
  // Leaves for the track tree
  vector<float> *fTrackPtArray;       // Array for track pT:s
  vector<float> *fTrackPhiArray;      // Array for track phis
  vector<float> *fTrackEtaArray;      // Array for track etas
  vector<int> *fTrackChargeArray;     // Array for track charges
  vector<int> *fTrackSubeventArray;   // Array for track subevent indices (0 = PYTHIA, (>0) = HYDJET)
  
};

#endif