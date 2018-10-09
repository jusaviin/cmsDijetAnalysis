#ifndef JFFCORRECTOR_H
#define JFFCORRECTOR_H

// Root includes
#include <TFile.h>
#include <TH1.h>

// Own includes
#include "DijetCard.h"
#include "DijetHistogramManager.h"

/*
 * Class for drawing the histograms produced in the dijet analysis
 */
class JffCorrector {

public:
  
  JffCorrector();                                       // Default constructor
  JffCorrector(TFile *inputFile);                       // Constructor
  JffCorrector(TFile *inputFile, TFile* spilloverFile); // Constructor
  JffCorrector(const JffCorrector& in);                 // Copy constructor
  ~JffCorrector();                                      // Destructor
  
  // Setter for input file
  void ReadInputFile(TFile *inputFile); // Setter for input file
  void ReadSpilloverFile(TFile *spilloverFile);  // Setter for spillover file
  
  // Getters JFF correction histograms
  TH1D* GetJetShapeJffCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt) const;  // Jet shape JFF correction histograms
  TH2D* GetDeltaEtaDeltaPhiJffCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt) const;  // DeltaEta-DeltaPhi JFF correction histograms
  TH2D* GetDeltaEtaDeltaPhiSpilloverCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt) const;  // DeltaEta-DeltaPhi spillover correction histograms
  
  // Return information, if correction is ready to be obtained
  bool CorrectionReady();  // True if histograms loaded from file, otherwise false
  bool SpilloverReady();   // True if spillover correction is loaded, otherwise false
  
private:
  
  // Data members
  bool fFileLoaded;                // Flag if the input file has been loaded
  bool fSpilloverLoaded;           // Flag if the spillover file has been loaded

  // JFF correction histograms for jet shape
  TH1D *fhJetShapeCorrection[DijetHistogramManager::knJetTrackCorrelations][DijetHistogramManager::knCentralityBins][DijetHistogramManager::knTrackPtBins];  // JFF correction histograms for jet shape
  TH2D *fhDeltaEtaDeltaPhiCorrection[DijetHistogramManager::knJetTrackCorrelations][DijetHistogramManager::knCentralityBins][DijetHistogramManager::knTrackPtBins];  // JFF correction histograms for deltaEta-deltaPhi histograms
  TH2D *fhDeltaEtaDeltaPhiSpilloverCorrection[DijetHistogramManager::knJetTrackCorrelations][DijetHistogramManager::knCentralityBins][DijetHistogramManager::knTrackPtBins];
  
};

#endif
