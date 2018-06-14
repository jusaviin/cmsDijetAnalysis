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
#include <vector>

// Root includes
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>

using namespace std;

class ForestReader{
  
public:
  
  // Possible data types to be read with the reader class
  enum enumDataTypes{kPp, kPbPb, kPpMC, kPbPbMC, kLocalTest, knDataTypes};
  
  // Constructors and destructors
  ForestReader();                                   // Default constructor
  ForestReader(Int_t dataType);                     // Custom constructor
  ForestReader(const ForestReader& in);             // Copy constructor
  virtual ~ForestReader();                          // Destructor
  ForestReader& operator=(const ForestReader& obj); // Equal sign operator
  
  // Methods
  virtual void GetEvent(Int_t nEvent) = 0;                 // Get the nth event in tree
  Int_t GetNEvents() const;                                // Get the number of events
  virtual void ReadForestFromFile(TFile *inputFile) = 0;   // Read the forest from a file
  virtual void BurnForest() = 0;                           // Burn the forest
  
  // Getters for leaves in heavy ion tree
  Float_t GetVz() const;              // Getter for vertex z position
  Float_t GetCentrality() const;      // Getter for centrality
  Int_t GetHiBin() const;             // Getter for CMS hiBin
  Float_t GetPtHat() const;           // Getter for pT hat
  
  // Getters for leaves in jet tree
  virtual Float_t GetJetPt(Int_t iJet) const = 0;         // Getter for jet pT
  virtual Float_t GetJetPhi(Int_t iJet) const = 0;        // Getter for jet phi
  virtual Float_t GetJetEta(Int_t iJet) const = 0;        // Getter for jet eta
  virtual Int_t GetNJets() const;                         // Getter for number of jets
  virtual Float_t GetJetRawPt(Int_t iJet) const = 0;      // Getter for jet raw pT
  virtual Float_t GetJetMaxTrackPt(Int_t iJet) const = 0; // Getter for maximum track pT inside a jet
  
  // Getters for leaves in HLT tree
  Int_t GetCaloJetFilterBit() const;  // Getter for calorimeter jet filter bit
  
  // Getters for leaves in skim tree
  Int_t GetPrimaryVertexFilterBit() const;           // Getter for primary vertex filter bit
  Int_t GetBeamScrapingFilterBit() const;            // Getter got beam scraping filter bit
  Int_t GetHBHENoiseFilterBit() const;               // Getter for HB/HE noise filter bit
  Int_t GetCollisionEventSelectionFilterBit() const; // Getter for collision event selection filter bit
  Int_t GetHfCoincidenceFilterBit() const;           // Getter for hadronic forward coincidence filter bit
  Int_t GetClusterCompatibilityFilterBit() const;    // Getter for cluster compatibility filter bit
  
  // Getters for leaves in the track tree
  virtual Float_t GetTrackPt(Int_t iTrack) const = 0;                    // Getter for track pT
  virtual Float_t GetTrackPtError(Int_t iTrack) const = 0;               // Getter for track pT error
  virtual Float_t GetTrackPhi(Int_t iTrack) const = 0;                   // Getter for track phi
  virtual Float_t GetTrackEta(Int_t iTrack) const = 0;                   // Getter for track eta
  Int_t GetNTracks() const;                                              // Getter for number of tracks
  virtual Bool_t GetTrackHighPurity(Int_t iTrack) const = 0;             // Getter for the high purity of the track
  virtual Float_t GetTrackVertexDistanceZ(Int_t iTrack) const = 0;       // Getter for track distance from primary vertex in z-direction
  virtual Float_t GetTrackVertexDistanceZError(Int_t iTrack) const = 0;  // Getter for error of track distance from primary vertex in z-direction
  virtual Float_t GetTrackVertexDistanceXY(Int_t iTrack) const = 0;      // Getter for track distance from primary vertex in xy-direction
  virtual Float_t GetTrackVertexDistanceXYError(Int_t iTrack) const = 0; // Getter for error of track distance from primary vertex in xy-direction
  virtual Float_t GetTrackChi2(Int_t iTrack) const = 0;                  // Getter for track chi2 value from reconstruction fit
  virtual Int_t GetNTrackDegreesOfFreedom(Int_t iTrack) const = 0;       // Getter for number of degrees of freedom in reconstruction fit
  virtual Int_t GetNHitsTrackerLayer(Int_t iTrack) const = 0;            // Getter for number of hits in tracker layers
  virtual Int_t GetNHitsTrack(Int_t iTrack) const = 0;                   // Getter for number of hits for the track
  virtual Float_t GetTrackEnergyEcal(Int_t iTrack) const = 0;            // Getter for track energy in ECal
  virtual Float_t GetTrackEnergyHcal(Int_t iTrack) const = 0;            // Getter for track energy in HCal
  virtual Int_t GetTrackCharge(Int_t iTrack) const = 0;                  // Getter for track charge (only for generator level tracks)
  virtual Int_t GetTrackSubevent(Int_t iTrack) const = 0;                // Getter for track subevent index (only for generator level tracks)
  virtual Int_t GetTrackMCStatus(Int_t iTrack) const = 0;                // Getter for track MC status (only for generator level tracks)
  
  // Getters for leaves in the particle flow candidate tree
  virtual Int_t GetParticleFlowCandidateId(Int_t iCandidate) const;      // Getter for particle flow candidate ID
  virtual Float_t GetParticleFlowCandidatePt(Int_t iCandidate) const;    // Getter for particle flow candidate pT
  virtual Float_t GetParticleFlowCandidatePhi(Int_t iCandidate) const;   // Getter for particle flow candidate phi
  virtual Float_t GetParticleFlowCandidateEta(Int_t iCandidate) const;   // Getter for particle flow candidate eta
  virtual Int_t GetNParticleFlowCandidates() const;                      // Getter for number of particle flow candidates in an event
  
  // Setter for data type
  void SetDataType(Int_t dataType); // Setter for data type
  
protected:
  
  // Methods
  virtual void Initialize() = 0;  // Connect the branches to the tree
  
  Int_t fDataType;  // Type of data read with the tree. 0 = pp, 1 = PbPb, 2 = ppMC, 3 = PbPbMC, 4 = LocalTest
  
  // Branches for heavy ion tree
  TBranch *fHiVzBranch;            // Branch for vertex z-position
  TBranch *fHiBinBranch;           // Branch for centrality
  TBranch *fPtHatBranch;           // Branch for pT hat
  
  // Branches for jet tree
  TBranch *fJetPtBranch;         // Branch for jet pT
  TBranch *fJetPhiBranch;        // Branch for jet phi
  TBranch *fJetEtaBranch;        // Branch for jet eta
  TBranch *fJetRawPtBranch;      // Branch for raw jet pT
  TBranch *fJetMaxTrackPtBranch; // Maximum pT for a track inside a jet
  
  // Branches for HLT tree
  TBranch *fCaloJetFilterBranch;    // Branch for calo jet filter bit
  
  // Branches for skim tree
  TBranch *fPrimaryVertexBranch;           // Branch for primary vertex filter bit
  TBranch *fBeamScrapingBranch;            // Branch for beam scraping filter bit
  TBranch *fCollisionEventSelectionBranch; // Branch for collision event selection filter bit
  TBranch *fHBHENoiseBranch;               // Branch for HB/HE noise filter bit
  TBranch *fHfCoincidenceBranch;           // Branch for energy recorded in at least 3 HF calorimeter towers
  TBranch *fClusterCompatibilityBranch;    // Branch for cluster compatibility
  
  // Branches for track tree
  TBranch *fTrackPtBranch;                    // Branch for track pT
  TBranch *fTrackPtErrorBranch;               // Branch for track pT error
  TBranch *fTrackPhiBranch;                   // Branch for track phi
  TBranch *fTrackEtaBranch;                   // Branch for track eta
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
  
  // Branched for particle flow candidate tree
  TBranch *fParticleFlowCandidateIdBranch;    // Branch for particle flow candidate ID
  TBranch *fParticleFlowCandidatePtBranch;    // Branch for particle flow candidate pT
  TBranch *fParticleFlowCandidatePhiBranch;   // Branch for particle flow candidate phi
  TBranch *fParticleFlowCandidateEtaBranch;   // Branch for particle flow candidate eta
    
  // Leaves for heavy ion tree
  Float_t fVertexZ;    // Vertex z-position
  Int_t fHiBin;        // HiBin = Centrality percentile * 2
  Float_t fPtHat;      // pT hat
  
  // Leaves for jet tree
  Int_t fnJets;                                // number of jets in an event
  
  // Leaves for the HLT tree
  Int_t fCaloJetFilterBit;    // Filter bit for calorimeter jets
  
  // Leaves for the skim tree
  Int_t fPrimaryVertexFilterBit;           // Filter bit for primary vertex
  Int_t fBeamScrapingFilterBit;            // Filter bit for beam scraping
  Int_t fCollisionEventSelectionFilterBit; // Filter bit for collision event selection
  Int_t fHBHENoiseFilterBit;               // Filter bit for HB/HE noise
  Int_t fHfCoincidenceFilterBit;           // Filter bit for energy recorded in at least 3 HF calorimeter towers
  Int_t fClusterCompatibilityFilterBit;    // Filter bit for cluster compatibility
  
  // Leaves for the track tree
  Int_t fnTracks;  // Number of tracks
  
  // Leaves for the particle flow candidate tree
  vector<int> *fParticleFlowCandidateIdArray;     // Vector for particle flow candidate ID:s
  vector<float> *fParticleFlowCandidatePtArray;   // Vector for particle flow candidate pT:s
  vector<float> *fParticleFlowCandidatePhiArray;  // Vector for particle flow candidate phis
  vector<float> *fParticleFlowCandidateEtaArray;  // Vector for particle flow candidate etas
};

#endif






















