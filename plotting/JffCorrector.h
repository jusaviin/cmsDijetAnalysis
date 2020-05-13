#ifndef JFFCORRECTOR_H
#define JFFCORRECTOR_H

// C++ includes
#include <sstream>
#include <fstream>
#include <string>

// Root includes
#include <TFile.h>
#include <TH1.h>

// Own includes
#include "DijetCard.h"
#include "DijetHistogramManager.h"
#include "DijetMethods.h"

/*
 * Class for drawing the histograms produced in the dijet analysis
 */
class JffCorrector {

public:
  
  enum enumUncertaintySources{kBackgroundFluctuation, kFragmentationBias, kJetEnergyScale, kTriggerEfficiency, kTrackingEfficiency, kResidualTracking, kTrackingDeltaR, kPairAcceptance, kBackgroundSubtraction, kTotal, knUncertaintySources};
  enum enumLongRangeUncertaintySources{kBackgroundGlue, kEtaSide, kEtaRegion, kSameMixed, kVzVariation, kTotalLongRange, knLongRangeUncertaintySources};

  
  JffCorrector();                                       // Default constructor
  JffCorrector(TFile *inputFile);                       // Constructor
  JffCorrector(TFile *inputFile, TFile* spilloverFile); // Constructor
  JffCorrector(TFile *inputFile, TFile* spilloverFile, TFile* trackingFile); // Constructor
  JffCorrector(const JffCorrector& in);                 // Copy constructor
  ~JffCorrector();                                      // Destructor
  
  // Setter for input file
  void ReadInputFile(TFile *inputFile);               // Read the histograms related to JFF correction
  void ReadSpilloverFile(TFile *spilloverFile);       // Read the histograms related to spillover correction
  void ReadSpilloverDeltaRFile(TFile *spilloverFile); // Read the spillover correction histograms as a function of DeltaR
  void ReadTrackDeltaRFile(TFile *trackFile);         // Read the histograms related to residual R-dependent tracking correction
  void ReadSystematicFile(TFile *systematicFile);     // Read the histograms related to systematic uncertainties
  void ReadJetReconstructionBiasFile(const char *fileName);     // Read the correction to v2 due to jet reconstruction bias
  void ReadLongRangeSystematicFile(const char *systematicFile); // Read the histograms related to systematic uncertainties of long range correlations
  
  // Getters for correction histograms
  TH1D* GetJetShapeJffCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry = DijetHistogramManager::kMaxAsymmetryBins) const;  // Jet shape JFF correction histograms
  TH2D* GetDeltaEtaDeltaPhiJffCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry = DijetHistogramManager::kMaxAsymmetryBins) const;  // DeltaEta-DeltaPhi JFF correction histograms
  TH2D* GetDeltaEtaDeltaPhiSpilloverCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry = DijetHistogramManager::kMaxAsymmetryBins) const;  // DeltaEta-DeltaPhi spillover correction histograms
  TH2D* GetDeltaEtaDeltaPhiSpilloverCorrectionAsymmetryScale(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, const int iAsymmetry = DijetHistogramManager::kMaxAsymmetryBins, int usedBin = -1) const;  // DeltaEta-DeltaPhi spillover correction histograms. Use scaled xj integrated distribution to calculate correction in different xj bins
  TH1D* GetJetShapeSpilloverCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry = DijetHistogramManager::kMaxAsymmetryBins) const; // Spillover correction as a function of DeltaR
  TH1D* GetJetShapeSpilloverCorrectionManualTune(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry = DijetHistogramManager::kMaxAsymmetryBins) const; // Manually tuned spillover correction as a function of DeltaR
  TH2D* GetDeltaEtaDeltaPhiTrackDeltaRCorrection(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry = DijetHistogramManager::kMaxAsymmetryBins) const; // DeltaEta-DeltaPhi residual tracking correction histograms
  double GetTrackDeltaRResidualScale(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry = DijetHistogramManager::kMaxAsymmetryBins) const; // DeltaEta-DeltaPhi residual tracking correction histograms
  
  // Getters for corrections for long range correlations
  double GetJetReconstructionBiasCorrection(const int iFlow, const int iCentrality, const int iTrackPt, int iAsymmetry = DijetHistogramManager::kMaxAsymmetryBins) const; // Jet reconstruction bias correction
  
  // Getters related to systematic uncertainties
  TString GetUncertaintyName(const int iUncertainty) const;
  TH1D* GetJetShapeSystematicUncertainty(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry = DijetHistogramManager::kMaxAsymmetryBins, int iUncertainty = kTotal) const;
  TH1D* GetDeltaEtaSystematicUncertainty(const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt, int iAsymmetry = DijetHistogramManager::kMaxAsymmetryBins, int iUncertainty = kTotal) const;
  
  TString GetLongRangeUncertaintyName(const int iUncertainty) const;
  double GetLongRangeSystematicUncertainty(const int iFlow, const int iCentrality, const int iTrackPt, int iAsymmetry = DijetHistogramManager::kMaxAsymmetryBins) const;
  
  // Setters
  void SetUncertaintySmooth(const bool smooth);           // Setter for smoothing the uncertainties
  void SetDeltaEtaSymmetrization(const bool symmetrize);  // Setter for symmetrizing the deltaEta uncertainties
  
  // Return information, if correction is ready to be obtained
  bool CorrectionReady();         // True if histograms loaded from file, otherwise false
  bool SpilloverReady();          // True if spillover correction is loaded, otherwise false
  bool SpilloverDeltaRReady();    // True if spillover histograms as a function of deltaR are loaded
  bool SystematicsReady();        // True if systematic uncertainties are loaded, otherwise false
  bool TrackingCorrectionReady(); // True if residual tracking correction in loaded
  
private:
  
  // Names for different uncertainties
  TString uncertaintyName[knUncertaintySources] = {"backgroundFluctuation", "fragmentationBias", "jetEnergyScale", "triggerEfficiency", "trackingEfficiency", "residualTracking", "trackingDeltaR", "pairAcceptance", "backgroundSubtraction", "total"};
  
  TString longRangeUncertaintyName[knLongRangeUncertaintySources] = {"backgroundShift", "etaSide", "etaRegion", "sameVsMixed", "vzSelection", "total"};
  
  // Data members
  bool fFileLoaded;                // Flag if the input file has been loaded
  int  fJffAsymmetryBins;          // Number of asymmetry bins in the JFF correction file
  int  fJffTrackPtBins;            // Number of track pT bins in the JFF correction file
  bool fSpilloverLoaded;           // Flag if the spillover file has been loaded
  bool fSpilloverDeltaRLoaded;     // Flag if the spillover histograms as a function of DeltaR have been loaded
  int  fSpilloverAsymmetryBins;    // Number of asymmetry bins in the spillover file
  int  fSpilloverTrackPtBins;      // Number of track pT bins in the spillover file
  bool fSystematicErrorLoaded;     // Flag if a systematic error file has been loaded
  int  fSystematicAsymmetryBins;   // Number of asymmetry bins in the systematic uncertainty file
  int  fSystematicTrackPtBins;     // Number of track pT bins in the systematic uncertainty file
  int  fLongRangeAsymmetryBins;    // Number of asymmetry bins in the long range correction file
  bool fTrackingCorrectionLoaded;  // Flag if the tracking correction had been loaded
  int  fTrackingAsymmetryBins;     // Number of asymmetry bins in the residual tracking correction file
  int  fTrackingPtBins;            // Number of track pT bins in the residual tracking correction file
  bool fSmoothUncertainty;         // Smooth the uncertainty histograms
  bool fSymmetrizeDeltaEta;        // Symmetrize the deltaEta uncertainty histograms

  // JFF correction histograms for jet shape
  TH1D *fhJetShapeCorrection[DijetHistogramManager::knJetTrackCorrelations][DijetHistogramManager::kMaxAsymmetryBins+1][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins];  // JFF correction histograms for jet shape
  TH2D *fhDeltaEtaDeltaPhiCorrection[DijetHistogramManager::knJetTrackCorrelations][DijetHistogramManager::kMaxAsymmetryBins+1][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins];  // JFF correction histograms for deltaEta-deltaPhi histograms
  
  // Spillover correction
  TH2D *fhDeltaEtaDeltaPhiSpilloverCorrection[DijetHistogramManager::knJetTrackCorrelations][DijetHistogramManager::kMaxAsymmetryBins+1][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins];
  TH1D *fhJetShapeSpilloverCorrection[DijetHistogramManager::knJetTrackCorrelations][DijetHistogramManager::kMaxAsymmetryBins+1][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins];
  TH1D *fhJetShapeSpilloverCorrectionManualTune[DijetHistogramManager::knJetTrackCorrelations][DijetHistogramManager::kMaxAsymmetryBins+1][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins];
  
  // Residual tracking correction
  TH2D *fhDeltaEtaDeltaPhiTrackingDeltaRCorrection[DijetHistogramManager::knJetTrackCorrelations][DijetHistogramManager::kMaxAsymmetryBins+1][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins];
  TH1D *fhTrackDeltaRResidualScale[DijetHistogramManager::knJetTrackCorrelations][DijetHistogramManager::kMaxAsymmetryBins+1][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins];
  
  // Systematic uncertainty
  TH1D *fhJetShapeUncertainty[DijetHistogramManager::knJetTrackCorrelations][DijetHistogramManager::kMaxAsymmetryBins+1][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins][knUncertaintySources];
  TH1D *fhDeltaEtaUncertainty[DijetHistogramManager::knJetTrackCorrelations][DijetHistogramManager::kMaxAsymmetryBins+1][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins][knUncertaintySources];
  
  // Corrections to Fourier components from jet reconstruction bias
  double fJetReconstructionBiasCorrection[DijetHistogramManager::kMaxAsymmetryBins+1][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins][4] = {{{{0}}}}; // Last bin = Different flow components
  
  // Systematic uncertainty for long range correlations
  double fLongRangeUncertaintyTable[DijetHistogramManager::kMaxAsymmetryBins+1][DijetHistogramManager::kMaxCentralityBins][DijetHistogramManager::kMaxTrackPtBins][4] = {{{{0}}}}; // Last bin = Different flow components
  
};

#endif
