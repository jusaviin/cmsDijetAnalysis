// Reader for trees needed in event mixing from CMS data
//
//===========================================================
// GeneratorLevelMixingForestReader.h
//
// Author: Jussi Viinikainen
//===========================================================

#ifndef GENERATORLEVELMIXINGFORESTREADER_H
#define GENERATORLEVELMIXINGFORESTREADER_H

// Own includes
#include "ForestReader.h"

using namespace std;

class GeneratorLevelMixingForestReader : public ForestReader{
  
private:
  static const Int_t fnMaxTrack = 60000;   // Maximum number of tracks in an event
  
public:
  
  // Constructors and destructors
  GeneratorLevelMixingForestReader();                                              // Default constructor
  GeneratorLevelMixingForestReader(Int_t dataType, Int_t readMode, Int_t jetType, Int_t jetAxis, Bool_t matchJets, Bool_t doEventPlane, Bool_t minimumBiasMode); // Custom constructor
  GeneratorLevelMixingForestReader(const GeneratorLevelMixingForestReader& in);                    // Copy constructor
  virtual ~GeneratorLevelMixingForestReader();                                     // Destructor
  GeneratorLevelMixingForestReader& operator=(const GeneratorLevelMixingForestReader& obj);        // Equal sign operator
  
  // Methods
  void ReadForestFromFile(TFile *inputFile);   // Read the forest from a file
  void ReadForestFromFileList(std::vector<TString> inputFileList);   // Read the forest from a file
  void BurnForest();                           // Burn the forest  
  void GetEvent(Int_t nEvent);                 // Get the nEventh event from the file
  
  // Getter for number of events
  Int_t GetNEvents() const;                   // Get the number of events
  
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
  Int_t GetNTrackDegreesOfFreedom(Int_t iTrack) const;       // Getter for number of degrees of freedom in reconstruction fit
  Int_t GetNHitsTrackerLayer(Int_t iTrack) const;            // Getter for number of hits in tracker layers
  Int_t GetNHitsTrack(Int_t iTrack) const;                   // Getter for number of hits for the track
  Float_t GetTrackEnergyEcal(Int_t iTrack) const;            // Getter for track energy in ECal
  Float_t GetTrackEnergyHcal(Int_t iTrack) const;            // Getter for track energy in HCal
  Int_t GetTrackCharge(Int_t iTrack) const;                  // Getter for track charge (relevant only for generator level tracks)
  Int_t GetTrackSubevent(Int_t iTrack) const;                // Getter for track subevent index (relevant only for generator level tracks)
  Int_t GetTrackMCStatus(Int_t iTrack) const;                // Getter for track MC status (only for generator level tracks)
  
  // New variables for 2018 data
  Int_t GetTrackAlgorithm(Int_t iTrack) const;               // Getter for track algorithm
  Float_t GetTrackMVA(Int_t iTrack) const;                   // Getter for track MVA
  
  // Check if generator level jet has a matching reconstructed jet
  Bool_t HasMatchingJet(Int_t iJet) const;      // Check if generator level jet has a matching reconstructed jet
  Float_t GetMatchedPt(Int_t iJet) const;       // Getter for matched generator level jet pT
  Float_t GetMatchedEta(Int_t iJet) const;      // Getter for matched generator level jet eta
  Float_t GetMatchedPhi(Int_t iJet) const;      // Getter for matched generator level jet phi
  Int_t GetPartonFlavor(Int_t iJet) const;      // Parton flavor for the parton initiating the jet
  
private:
  
  // Methods
  void Initialize();       // Connect the branches to the tree

  // Trees in the forest
  TChain *fHeavyIonTree;    // Tree for heavy ion event information
  TChain *fSkimTree;        // Tree for event selection information
  TChain *fTrackTree;       // Tree for tracks  PbPb: anaTrack/trackTree pp: ppTrack/trackTree GenParticles: HiGenParticleAna/hi
  
  // Leaves for the track tree
  vector<float> *fTrackPtArray;       // Array for track pT:s
  vector<float> *fTrackPhiArray;      // Array for track phis
  vector<float> *fTrackEtaArray;      // Array for track etas
  vector<int> *fTrackChargeArray;     // Array for track charges
  vector<int> *fTrackSubeventArray;   // Array for track subevent indices (0 = PYTHIA, (>0) = HYDJET)
  
};

#endif
