// Reader for jet trees from CMS data
//
//===========================================================
// SkimForestReader.h
//
// Author: Jussi Viinikainen
//===========================================================

#ifndef SKIMFORESTREADER_H
#define SKIMFORESTREADER_H

// Own includes
#include "ForestReader.h"

using namespace std;

class SkimForestReader : public ForestReader{
  
public:
  
  // Constructors and destructors
  SkimForestReader();                                   // Default constructor
  SkimForestReader(Int_t dataType);                       // Custom constructor
  SkimForestReader(const SkimForestReader& in);             // Copy constructor
  virtual ~SkimForestReader();                          // Destructor
  SkimForestReader& operator=(const SkimForestReader& obj); // Equal sign operator
  
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
  Int_t GetNTrackDegreesOfFreedom(Int_t iTrack) const;       // Getter for number of degrees of freedom in reconstruction fit
  Int_t GetNHitsTrackerLayer(Int_t iTrack) const;            // Getter for number of hits in tracker layers
  Int_t GetNHitsTrack(Int_t iTrack) const;                   // Getter for number of hits for the track
  Float_t GetTrackEnergyEcal(Int_t iTrack) const;            // Getter for track energy in ECal
  Float_t GetTrackEnergyHcal(Int_t iTrack) const;            // Getter for track energy in HCal
  Int_t GetTrackCharge(Int_t iTrack) const;                  // Getter for track charge (relevant only for generator level tracks)
  Int_t GetTrackSubevent(Int_t iTrack) const;                // Getter for track subevent index (relevant only for generator level tracks)
  
private:
  
  // Methods
  void Initialize();  // Connect the branches to the tree
  
  // Trees in the forest
  TTree *fEventTree;    // Tree for heavy ion event information
  
  // Leaves for jet tree
  vector<float> *fJetPtArray;        // pT:s of all the jets in an event
  vector<float> *fJetPhiArray;        // phis of all the jets in an event
  vector<float> *fJetEtaArray;        // etas of all the jets in an event
  vector<float> *fJetRawPtArray;      // raw jet pT for all the jets in an event
  vector<float> *fJetMaxTrackPtArray; // maximum track pT inside a jet for all the jets in an event
  
  // Leaves for the track tree
  vector<float> *fTrackPtArray;                    // vector for track pT:s
  vector<float> *fTrackPtErrorArray;               // vector for track pT errors
  vector<float> *fTrackPhiArray;                   // vector for track phis
  vector<float> *fTrackEtaArray;                   // vector for track etas
  vector<bool> *fHighPurityTrackArray;             // vector for the high purity of tracks
  vector<float> *fTrackVertexDistanceZArray;       // vector for track distance from primary vertex in z-direction
  vector<float> *fTrackVertexDistanceZErrorArray;  // vector for error for track distance from primary vertex in z-direction
  vector<float> *fTrackVertexDistanceXYArray;      // vector for track distance from primary vertex in xy-direction
  vector<float> *fTrackVertexDistanceXYErrorArray; // vector for error for track distance from primary vertex in xy-direction
  vector<float> *fTrackChi2Array;                  // vector for track chi2 value from reconstruction fit
  vector<int> *fnTrackDegreesOfFreedomArray;       // vector for number of degrees of freedom in reconstruction fit
  vector<int> *fnHitsTrackerLayerArray;            // vector for number of hits in tracker layers
  vector<int> *fnHitsTrackArray;                   // vector for number of hits for the track
  vector<float> *fTrackEnergyEcalArray;            // vector for track energy in ECal
  vector<float> *fTrackEnergyHcalArray;            // vector for track energy in HCal
  
};

#endif






















