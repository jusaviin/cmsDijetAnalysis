#ifndef DIJETCARD_H
#define DIJETCARD_H

// C++ includes
#include <iostream>

// Root includes
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>
#include <TVectorT.h>

/*
 * DijetCard class
 *
 * This class reads the ConfigurationCard from the input root file and decodes the
 * necessary information for the analysis.
 */
class DijetCard {
  
private:
  
  TFile *fInputFile;         // Input file from which all the data is read
  TString fCardDirectory;    // Path to the ConfigurationCard directory
  int fDataType;             // Total number of centrality bins in the analysis
  int fMonteCarloType;       // Type of Monte Carlo used for jet-track correlations
  int fAsymmetryBinType;     // Used asymmetry binning (0 = AJ, 1 = xJ)
  TString fDataTypeString;   // Total number of eta gaps in the analysis
  
  void FindDataTypeString(); // Construct a data type string based on information on the card
  void ReadVectors();        // Read the vectors from the file
  
  // Vectors for all the lines inside the card
  TVectorT<float> *fDataTypeVector;                 // Vector for data type
  TVectorT<float> *fMcCorrelationTypeVector;        // Vector for Monte Carlo correlation type
  TVectorT<float> *fMatchJetsVector;                // Vector for jet matching information
  TVectorT<float> *fForestTypeVector;               // Vector for the used forest type
  TVectorT<float> *fReadModeVector;                 // Vector telling if the used file was skim of high forest
  TVectorT<float> *fJetTypeVector;                  // Vector telling if we used particle flow or calorimeter jets
  TVectorT<float> *fJetAxisVector;                  // Vector telling if we used anti-kT or leading PF candidate jet axis
  TVectorT<float> *fJetEtaCutVector;                // Vector for the used eta cut for jets
  TVectorT<float> *fSearchEtaCutVector;             // Vector for the eta region from which the jets are searched in the analysis
  TVectorT<float> *fMaxPtCutVector;                 // Vector for maximum pT accepted for the leading jet
  TVectorT<float> *fMinPtCutVector;                 // Vector for the minimum pT accepted for the leading jet
  TVectorT<float> *fSubleadingPtCutVector;          // Vector for minimum pT allowed for the subleading jet
  TVectorT<float> *fDeltaPhiCutVector;              // Vector for minimum deltaPhi accepted between leading and subleading jet axes
  TVectorT<float> *fMinMaxTrackPtFractionVector;    // Vector for minimum fraction of jet pT taken by the highest pT track in jet
  TVectorT<float> *fMaxMaxTrackPtFractionVector;    // Vector for maximum fraction of jet pT taken by the highest pT track in jet
  TVectorT<float> *fTrackEtaCutVector;              // Vector for the eta cut for tracks
  TVectorT<float> *fMinTrackPtCutVector;            // Vector for the minimum accepted track pT
  TVectorT<float> *fMaxTrackPtRelativeErrorVector;  // Vector for the maximum relative error allowed for track pT
  TVectorT<float> *fVertexMaxDistanceVector;        // Vector for maximum allowed distance of tracks from reconstructed vertex
  TVectorT<float> *fCalorimeterSignalLimitPtVector; // Vector for limit for track pT above which a signal in calorimeters is required
  TVectorT<float> *fHighPtEtFractionVector;         // Vector for the minimum fraction between pT and Et for high pT tracks
  TVectorT<float> *fChi2QualityCutVector;           // Vector for the maximum accepted chi2 for reconstructed tracks
  TVectorT<float> *fMinimumTrackHitsVector;         // Vector for minimum number of hits in tracking for a track
  TVectorT<float> *fSubeventCutVector;              // Vector for cut on subevent index
  TVectorT<float> *fZVertexCutVector;               // Vector for maximum accepted vz in the event
  TVectorT<float> *fLowPtHatCutVector;              // Vector for the minimum accepted pT hat
  TVectorT<float> *fHighPtHatCutVector;             // Vector for the maximum accepted pT hat
  TVectorT<float> *fAsymmetryBinTypeVector;         // Vector for the used asymmetry binning type (AJ or xJ)
  TVectorT<float> *fCentralityBinEdgesVector;       // Vector for the used centrality bin edges
  TVectorT<float> *fTrackPtBinEdgesVector;          // Vector for the used track pT bin edges
  TVectorT<float> *fAsymmetryBinEdgesVector;        // Vector for the used dijet asymmetry bin edges
  TVectorT<float> *fPtHatBinEdgesVector;            // Vector for the used pT hat bin edges
  TVectorT<float> *fDoEventMixingVector;            // Vector for the event mixing flag
  TVectorT<float> *fMixWithPoolVector;              // Vector telling if we do mixing with or without mixing pool
  TVectorT<float> *fNMixedEventsPerDijetVector;     // Vector for the number of events used in the event mixing
  TVectorT<float> *fVzToleranceVector;              // Vector for maximum accepted distance between same event vz and mixed event vz
  TVectorT<float> *fMixingVzBinWidthVector;         // Vector telling vz bin width, in case we use mixing pool
  TVectorT<float> *fMixingHiBinWidthVector;         // Vector telling hibin width, in case we use mixing pool
  TVectorT<float> *fMixingPoolDepthVector;          // Vector telling the mixing pool depth
  
public:
  
  DijetCard(TFile *inFile); // Contructor with input file
  ~DijetCard();             // Destructor
  
  TString GetDataType() const;   // Getter for data type string
  void Write(TDirectory *file);  // Write the contents of the card to a file
  double GetMaxDeltaEta() const; // Get maximum deltaEta possible using jet and track cuts in the card
  
  int GetNAsymmetryBins() const; // Get the number of asymmetry bins
  double GetLowBinBorderAsymmetry(const int iBin) const;  // Get the low border of i:th asymmetry bin
  double GetHighBinBorderAsymmetry(const int iBin) const; // Get the high border of i:th asymmetry bin
  const char *GetAsymmetryBinType(TString latexIt = "") const; // Get a description of the used asymmetry bin type
  
};

#endif
