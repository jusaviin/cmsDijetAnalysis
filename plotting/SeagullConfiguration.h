#ifndef SEAGULLCONFIGURATION_H
#define SEAGULLCONFIGURATION_H

// Own includes
#include "DijetHistogramManager.h"

/*
 * SeagullConfiguration class
 *
 * Class whose purpose is to provide configuration for seagull correction
 */
class SeagullConfiguration {
  
public:
 
  SeagullConfiguration();   // Contructor
  ~SeagullConfiguration();  // Destructor
  
  int GetSeagullMethod(const int iJetTrack, const int iAsymmetry, const int iCentrality, const int iTrackPt) const;
  int GetSeagullVeto(const int iJetTrack, const int iAsymmetry, const int iCentrality, const int iTrackPt) const;
  
private:
  
  // Methods
  void InitializeArrays();

  // Seagull method arrays for different jet-track correlation types and datasets
  int trackLeadingJetPbPb[DijetHistogramManager::kMaxAsymmetryBins][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins];
  int trackSubleadingJetPbPb[DijetHistogramManager::kMaxAsymmetryBins][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins];
  int trackLeadingJetPbPbMCRecoGen[DijetHistogramManager::kMaxAsymmetryBins][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins];
  int trackLeadingJetPbPbMCRecoGenSubeNon0[DijetHistogramManager::kMaxAsymmetryBins][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins];
  
  // Seagull veto arrays for different jet-track correlation types and datasets
  int vetoTrackLeadingJetPbPb[DijetHistogramManager::kMaxAsymmetryBins][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins];
  int vetoTrackSubleadingJetPbPb[DijetHistogramManager::kMaxAsymmetryBins][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins];
  int vetoTrackLeadingJetPbPbMCRecoGen[DijetHistogramManager::kMaxAsymmetryBins][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins];

};

#endif
