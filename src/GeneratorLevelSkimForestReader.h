// Reader for jet trees from CMS data
//
//===========================================================
// GeneratorLevelSkimForestReader.h
//
// Author: Jussi Viinikainen
//===========================================================

#ifndef GENERATORLEVELSKIMFORESTREADER_H
#define GENERATORLEVELSKIMFORESTREADER_H

// Own includes
#include "SkimForestReader.h"

using namespace std;

class GeneratorLevelSkimForestReader : public SkimForestReader{
  
public:
  
  // Constructors and destructors
  GeneratorLevelSkimForestReader();                                                     // Default constructor
  GeneratorLevelSkimForestReader(Int_t dataType, Int_t readMode, Int_t jetType, Bool_t matchJets); // Custom constructor
  GeneratorLevelSkimForestReader(const GeneratorLevelSkimForestReader& in);             // Copy constructor
  virtual ~GeneratorLevelSkimForestReader();                                            // Destructor
  GeneratorLevelSkimForestReader& operator=(const GeneratorLevelSkimForestReader& obj); // Equal sign operator
  
  // Getters for leaves in jet tree
  Float_t GetJetRawPt(Int_t iJet) const;      // Getter for jet raw pT (not relevant for generator jets, just return value that passes cuts)
  Float_t GetJetMaxTrackPt(Int_t iJet) const; // Getter for maximum track pT inside a jet (not relevant for generator jets, just return value that passes cuts)
  
  // Getters for leaves in the track tree relevant for generator level tracks
  Int_t GetTrackCharge(Int_t iTrack) const;                  // Getter for track charge (relevant only for generator level tracks)
  Int_t GetTrackSubevent(Int_t iTrack) const;                // Getter for track subevent index (relevant only for generator level tracks)
  
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
  
  // Check if generator level jet has a matching reconstructed jet
  Bool_t HasMatchingJet(Int_t iJet) const;   // Check if generator level jet has a matching reconstructed jet
  
private:
  
  // Methods
  void Initialize();  // Connect the branches to the tree
  
  // Additional leaves for the track tree
  vector<int> *fTrackChargeArray;     // Array for track charges
  vector<int> *fTrackSubeventArray;   // Array for track subevent indices (0 = PYTHIA, (>0) = HYDJET)
  
};

#endif






















