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
#include <TChain.h>
#include <TBranch.h>
#include <TFile.h>

using namespace std;

class ForestReader{
  
protected:
  static const Int_t fnMaxParticleFlowCandidates = 10000;        // Maximum number of particle flow candidates
  static const Int_t fMaxEventPlanes = 30;                       // Maximum number of event planes
  
  /*
   * Key for event planes
   *
  Index     Name   Detector Order hmin1 hmax1 hmin2 hmax2 minpt maxpt nsub mcw    rmate1    rmate2
      0      HFm1        HF     1 -5.00 -3.00  0.00  0.00  0.01 30.00 3sub  no      HFp1   trackp1
      1      HFp1        HF     1  3.00  5.00  0.00  0.00  0.01 30.00 3sub  no      HFm1   trackm1
      2       HF1        HF     1 -5.00 -3.00  3.00  5.00  0.01 30.00 3sub  no   trackm1   trackp1
      3   trackm1   Tracker     1 -2.00 -1.00  0.00  0.00  0.30  3.00 3sub  no      HFm1      HFp1
      4   trackp1   Tracker     1  1.00  2.00  0.00  0.00  0.30  3.00 3sub  no      HFm1      HFp1
      5   Castor1    Castor     1 -6.55 -5.10  0.00  0.00  0.01 50.00 3sub  no      HFp1   trackp1
      6      HFm2        HF     2 -5.00 -3.00  0.00  0.00  0.01 30.00 3sub  no      HFp2 trackmid2
      7      HFp2        HF     2  3.00  5.00  0.00  0.00  0.01 30.00 3sub  no      HFm2 trackmid2
      8       HF2        HF     2 -5.00 -3.00  3.00  5.00  0.01 30.00 3sub  no   trackm2   trackp2
      9 trackmid2   Tracker     2 -0.75  0.75  0.00  0.00  0.30  3.00 3sub  no      HFm2      HFp2
     10   trackm2   Tracker     2 -2.00 -1.00  0.00  0.00  0.30  3.00 3sub  no      HFm2      HFp2
     11   trackp2   Tracker     2  1.00  2.00  0.00  0.00  0.30  3.00 3sub  no      HFm2      HFp2
     12   Castor2    Castor     2 -6.55 -5.10  0.00  0.00  0.01 50.00 3sub  no trackmid2      HFp2
     13      HFm3        HF     3 -5.00 -3.00  0.00  0.00  0.01 30.00 3sub  no      HFp3 trackmid3
     14      HFp3        HF     3  3.00  5.00  0.00  0.00  0.01 30.00 3sub  no      HFm3 trackmid3
     15       HF3        HF     3 -5.00 -3.00  3.00  5.00  0.01 30.00 3sub  no   trackm3   trackp3
     16 trackmid3   Tracker     3 -0.75  0.75  0.00  0.00  0.30  3.00 3sub  no      HFm3      HFp3
     17   trackm3   Tracker     3 -2.00 -1.00  0.00  0.00  0.30  3.00 3sub  no      HFm3      HFp3
     18   trackp3   Tracker     3  1.00  2.00  0.00  0.00  0.30  3.00 3sub  no      HFm3      HFp3
     19      HFm4        HF     4 -5.00 -3.00  0.00  0.00  0.01 30.00 3sub  no      HFp4 trackmid4
     20      HFp4        HF     4  3.00  5.00  0.00  0.00  0.01 30.00 3sub  no      HFm4 trackmid4
     21       HF4        HF     4 -5.00 -3.00  3.00  5.00  0.01 30.00 3sub  no   trackm4   trackp4
     22 trackmid4   Tracker     4 -0.75  0.75  0.00  0.00  0.30  3.00 3sub  no      HFm4      HFp4
     23   trackm4   Tracker     4 -2.00 -1.00  0.00  0.00  0.30  3.00 3sub  no      HFm4      HFp4
     24   trackp4   Tracker     4  1.00  2.00  0.00  0.00  0.30  3.00 3sub  no      HFm4      HFp4
     25    HFm1mc        HF     1 -5.00 -3.00  0.00  0.00  0.01 30.00 3sub yes    HFp1mc trackp1mc
     26    HFp1mc        HF     1  3.00  5.00  0.00  0.00  0.01 30.00 3sub yes    HFm1mc trackm1mc
     27 trackm1mc   Tracker     1 -2.20 -1.40  0.00  0.00  0.30  3.00 3sub yes    HFm1mc    HFp1mc
     28 trackp1mc   Tracker     1  1.40  2.20  0.00  0.00  0.30  3.00 3sub yes    HFm1mc    HFp1mc
     29 Castor1mc    Castor     1 -6.55 -5.10  0.00  0.00  0.01 50.00 3sub yes    HFp1mc trackp1mc
  */
  
public:
  
  // Possible data types to be read with the reader class
  enum enumDataTypes{kPp, kPbPb, kPpMC, kPbPbMC, kLocalTest, knDataTypes};
  
  // Constructors and destructors
  ForestReader();                                          // Default constructor
  ForestReader(Int_t dataType, Int_t readMode, Int_t jetType, Int_t jetAxis, Bool_t matchJets, Bool_t doEventPlane); // Custom constructor
  ForestReader(const ForestReader& in);                    // Copy constructor
  virtual ~ForestReader();                                 // Destructor
  ForestReader& operator=(const ForestReader& obj);        // Equal sign operator
  
  // Methods
  virtual void GetEvent(Int_t nEvent) = 0;                 // Get the nth event in tree
  virtual Int_t GetNEvents() const;                        // Get the number of events
  virtual void ReadForestFromFile(TFile *inputFile) = 0;   // Read the forest from a file
  virtual void ReadForestFromFileList(std::vector<TString> fileList) = 0;   // Read the forest from a file list
  virtual void BurnForest() = 0;                           // Burn the forest
  
  // Getters for leaves in heavy ion tree
  Float_t GetVz() const;              // Getter for vertex z position
  Float_t GetCentrality() const;      // Getter for centrality
  Int_t GetHiBin() const;             // Getter for CMS hiBin
  Float_t GetPtHat() const;           // Getter for pT hat
  Float_t GetEventWeight() const;     // Getter for jet weight in 2018 MC
  Int_t GenNEventPlane() const;       // Getter for the number of event planes
  Float_t GetEventPlaneAngle(Int_t iEventPlane) const;        // Getter for the event plane angle for the i:th event plane
  Float_t GetEventPlaneQ(Int_t iEventPlane) const;            // Getter for the magnitude of the q-vector for the i:th event plane
  Float_t GetEventPlaneMultiplicity(Int_t iEventPlane) const; // Getter for the particle multiplicity in the i:th event plane
  
  // Getters for leaves in jet tree
  virtual Float_t GetJetPt(Int_t iJet) const = 0;         // Getter for jet pT
  virtual Float_t GetJetPhi(Int_t iJet) const = 0;        // Getter for jet phi
  virtual Float_t GetJetEta(Int_t iJet) const = 0;        // Getter for jet eta
  virtual Int_t GetNJets() const;                         // Getter for number of jets
  virtual Float_t GetJetRawPt(Int_t iJet) const = 0;      // Getter for jet raw pT
  virtual Float_t GetJetMaxTrackPt(Int_t iJet) const = 0; // Getter for maximum track pT inside a jet
  
  // Getters for leaves in HLT tree
  Int_t GetCaloJetFilterBit() const;           // Getter for calorimeter jet filter bit
  Int_t GetPrescaledCaloJetFilterBit() const;  // Getter for calorimeter jet filter bit
  
  // Getters for leaves in skim tree
  Int_t GetPrimaryVertexFilterBit() const;           // Getter for primary vertex filter bit
  Int_t GetBeamScrapingFilterBit() const;            // Getter got beam scraping filter bit
  Int_t GetHBHENoiseFilterBit() const;               // Getter for HB/HE noise filter bit
  Int_t GetCollisionEventSelectionFilterBit() const; // Getter for collision event selection filter bit
  Int_t GetHfCoincidenceFilterBit() const;           // Getter for hadronic forward coincidence filter bit
  Int_t GetClusterCompatibilityFilterBit() const;    // Getter for cluster compatibility filter bit
  
  // Specific functions for jet closure plots
  virtual Bool_t HasMatchingJet(Int_t iJet) const = 0; // Check if generator level jet has a matching reconstructed jet
  virtual Float_t GetMatchedPt(Int_t iJet) const = 0;  // Getter for matched jet pT (reco for gen and vice versa)
  virtual Float_t GetMatchedEta(Int_t iJet) const = 0; // Getter for matched jet eta (reco for gen and vice versa)
  virtual Float_t GetMatchedPhi(Int_t iJet) const = 0; // Getter for matched jet phi (reco for gen and vice versa)
  virtual Int_t GetPartonFlavor(Int_t iJet) const = 0; // Parton flavor for the parton initiating the jet
  
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
  
  virtual Int_t GetTrackAlgorithm(Int_t iTrack) const;                   // Getter for track algorithm
  virtual Int_t GetTrackOriginalAlgorithm(Int_t iTrack) const;           // Getter for track original algorithm
  virtual Float_t GetTrackMVA(Int_t iTrack) const;                       // Getter for track MVA
  
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
  
  Int_t fDataType;     // Type of data read with the tree. 0 = pp, 1 = PbPb, 2 = ppMC, 3 = PbPbMC, 4 = LocalTest
  Int_t fReadMode;     // Different forests have different naming conventions. 0 = General forests, 1 = PYTHIA8 forest
  Int_t fJetType;      // Choose the type of jets usedfor analysis. 0 = Calo jets, 1 = PF jets
  Int_t fJetAxis;      // Jet axis used for the jets. 0 = Anti-kT, 1 = Leading particle flow candidate, 2 = WTA
  Bool_t fMatchJets;   // Match generator and reconstructed level jets
  Bool_t fDoEventPlane; // Include event plane branches in the tree
  
  // Branches for heavy ion tree
  TBranch *fHiVzBranch;                   // Branch for vertex z-position
  TBranch *fHiBinBranch;                  // Branch for centrality
  TBranch *fPtHatBranch;                  // Branch for pT hat
  TBranch *fnEventPlaneBranch;            // Branch for the number of event planes
  TBranch *fEventPlaneAngleBranch;        // Branch for the event plane angles
  TBranch *fEventPlaneQBranch;            // Branch for the event plane Q-vector magnitude
  TBranch *fEventPlaneMultiplicityBranch; // Branch for the particle multiplicity in the event plane
  
  // Branches for jet tree
  TBranch *fJetPtBranch;         // Branch for jet pT
  TBranch *fJetPhiBranch;        // Branch for jet phi
  TBranch *fJetEtaBranch;        // Branch for jet eta
  TBranch *fJetRawPtBranch;      // Branch for raw jet pT
  TBranch *fJetMaxTrackPtBranch; // Maximum pT for a track inside a jet
  TBranch *fJetRefPtBranch;      // Branch for reference generator level pT for a reconstructed jet
  TBranch *fJetRefFlavorBranch;  // Branch for flavor for the parton initiating the jet
  TBranch *fJetMatchedPtBranch;  // Branch for the matched jet pT (reco to gen or vice versa)
  TBranch *fJetMatchedEtaBranch; // Branch for the matched jet eta (reco to gen or vice versa)
  TBranch *fJetMatchedPhiBranch; // Branch for the matched jet phi (reco to gen or vice versa)
  TBranch *fEventWeightBranch;     // Branch for jet weight for 2018 MC
  
  // Branches for HLT tree
  TBranch *fCaloJetFilterBranch;         // Branch for calo jet filter bit
  TBranch *fCaloJetFilterPrescaleBranch; // Branch for prescaled calo jet filter bit
  
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
  TBranch *nfParticleFlowCandidateBranch;     // Branch for the number of particle flow candidates, needed for PYTHIA8 forest
  TBranch *fParticleFlowCandidateIdBranch;    // Branch for particle flow candidate ID
  TBranch *fParticleFlowCandidatePtBranch;    // Branch for particle flow candidate pT
  TBranch *fParticleFlowCandidatePhiBranch;   // Branch for particle flow candidate phi
  TBranch *fParticleFlowCandidateEtaBranch;   // Branch for particle flow candidate eta
    
  // Leaves for heavy ion tree
  Float_t fVertexZ;    // Vertex z-position
  Int_t fHiBin;        // HiBin = Centrality percentile * 2
  Float_t fPtHat;      // pT hat
  Int_t fnEventPlane;  // Number of event planes
  Float_t fEventPlaneAngle[fMaxEventPlanes] = {0};          // Event plane angles
  Float_t fEventPlaneQ[fMaxEventPlanes] = {0};              // Event plane q-vector magnitudes
  Float_t fEventPlaneMultiplicity[fMaxEventPlanes] = {0};   // Multiplicity of particles in the event plane
  
  // Leaves for jet tree
  Int_t fnJets;          // number of jets in an event
  Int_t fnMatchedJets;   // number of matched jets in an event
  Float_t fEventWeight;    // jet weight in the 2018 MC tree
  
  // Leaves for the HLT tree
  Int_t fCaloJetFilterBit;         // Filter bit for calorimeter jets
  Int_t fCaloJetFilterBitPrescale; // Prescaled filter bit needed for the minimum bias file in PbPb mixing
  
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
  Int_t fnParticleFlowCandidates;                                      // For PYTHIA8 forest an array is needed instead of vector
  Int_t fParticleFlowCandidateIdArray[fnMaxParticleFlowCandidates];    // For PYTHIA8 forest an array is needed instead of vector
  Float_t fParticleFlowCandidatePtArray[fnMaxParticleFlowCandidates];  // For PYTHIA8 forest an array is needed instead of vector
  Float_t fParticleFlowCandidatePhiArray[fnMaxParticleFlowCandidates]; // For PYTHIA8 forest an array is needed instead of vector
  Float_t fParticleFlowCandidateEtaArray[fnMaxParticleFlowCandidates]; // For PYTHIA8 forest an array is needed instead of vector
  vector<int> *fParticleFlowCandidateIdVector;                         // Vector for particle flow candidate ID:s
  vector<float> *fParticleFlowCandidatePtVector;                       // Vector for particle flow candidate pT:s
  vector<float> *fParticleFlowCandidatePhiVector;                      // Vector for particle flow candidate phis
  vector<float> *fParticleFlowCandidateEtaVector;                      // Vector for particle flow candidate etas
};

#endif
