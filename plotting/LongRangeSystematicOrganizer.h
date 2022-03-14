#ifndef LONGRANGESYSTEMATICORGANIZER_H
#define LONGRANGESYSTEMATICORGANIZER_H

// C++ includes
#include <iostream>

// Root includes
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMath.h>

/*
 * Class for drawing the histograms produced in the dijet analysis
 */
class LongRangeSystematicOrganizer {

private:
  static const int knMaxFlow = 4;  // Maximum number of flow components
  static const int knMaxXj = 3;    // Maximum number of xj bins
  
public:
  
  enum enumLongRangeUncertaintySources{kDeltaEtaRegion, kDeltaEtaRegionDihadron, kDeltaEtaSide, kVz, kAngleSmear, kJetCollection, kMCTuning, kMCMethod, kMCFit, kMinBias, kTracking, kQuarkGluonFraction, kJEC, kJER, kAll, knUncertaintySources};
  
  LongRangeSystematicOrganizer();                       // Default constructor
  LongRangeSystematicOrganizer(TFile *inputFile);                       // Constructor
  LongRangeSystematicOrganizer(const LongRangeSystematicOrganizer& in);                 // Copy constructor
  ~LongRangeSystematicOrganizer();                                      // Destructor
  
  // Setter for input file
  void ReadInputFile(TFile *inputFile);               // Read the long range systematic uncertainty file
  
  // Getters for the long range systematic uncertainties
  TString GetLongRangeUncertaintyName(const int iUncertainty, const bool groupedUncertainty = false) const;
  TString GetUncertaintyAxisName(const int iUncertainty, const bool groupedUncertainty = false) const;
  TGraphErrors* GetLongRangeSystematicUncertainty(const int iFlow, const int iUncertainty = kAll, int iAsymmetry = knMaxXj, const bool groupedUncertainty = false) const;
  
  int GetNUncertaintySources(const bool groupedUncertainty = false) const;
  
  // Define how to do histogram grouping
  void SetGroupingStrategy(const int groupingStrategy);
  
  // Adjust the central points for uncertainty graphs according to an input graph
  void AdjustCentralPoints(TGraphErrors* resultGraph[knMaxFlow]);
  
private:
  
  void GroupUncertaintyHistograms(); // Group the individual uncertainty histograms
  
  TString fLongRangeUncertaintyName[knUncertaintySources] = {"deltaEtaRegion", "deltaEtaRegionDihadron", "deltaEtaSide", "vz", "angleSmear", "jetCollection", "mcTuning", "mcMethod", "mcFit", "minBias", "tracking", "quarkGluonFraction", "jetEnergyCorrection", "jetEnergyResolution", "all"};
  TString fUncertaintyAxisName[knUncertaintySources] = {"#Delta#eta region", "#Delta#eta region hh", "#Delta#eta side", "v_{z}", "Angle smear", "Jet collection", "MC C-shift", "MC method", "MC Fit", "MinBias", "Tracking", "q/g fraction", "JEC", "JER", "all"};
  TString fGroupedUncertaintyName[knUncertaintySources] = {"Acceptance correction", "Long range extraction", "vz", "Jet axis resolution", "Jet reconstruction bias", "Dijet bias for dihadron", "Tracking", "JEC", "JER", "Total", "Empty", "Empty", "Empty", "Empty", "Empty"};
  TString fGroupedUncertaintyAxisName[knUncertaintySources] = {"Acceptance", "Long range", "v_{z}", "Angle smear", "Jet reco bias", "MinBias", "Tracking", "JEC", "JER", "Total", "Empty", "Empty", "Empty", "Empty", "Empty"};

  // Systematic uncertainty for long range correlations
  TGraphErrors* fLongRangeUncertaintyGraph[knUncertaintySources][knMaxFlow][knMaxXj+1];
  TGraphErrors* fGroupedLongRangeUncertaintyGraph[knUncertaintySources][knMaxFlow][knMaxXj+1];
  
  // Variables to deal with uncertainty source grouping
  int fnGroupedUncertaintySources;
  int fGroupingStrategy;  // 0 = Group similar uncertainties for paper. 1 = Group to remove old obsolete sources
  
};

#endif
