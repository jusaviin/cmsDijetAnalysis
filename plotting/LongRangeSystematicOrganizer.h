#ifndef LONGRANGESYSTEMATICORGANIZER_H
#define LONGRANGESYSTEMATICORGANIZER_H

// C++ includes
#include <iostream>

// Root includes
#include <TFile.h>
#include <TGraphErrors.h>

/*
 * Class for drawing the histograms produced in the dijet analysis
 */
class LongRangeSystematicOrganizer {

private:
  static const int knMaxFlow = 4;  // Maximum number of flow components
  static const int knMaxXj = 3;    // Maximum number of xj bins
  
public:
  
  enum enumLongRangeUncertaintySources{kDeltaEtaRegion, kDeltaEtaSide, kAngleSmear, kJetCollection, kMinBias, kAll, knUncertaintySources};

  
  LongRangeSystematicOrganizer();                       // Default constructor
  LongRangeSystematicOrganizer(TFile *inputFile);                       // Constructor
  LongRangeSystematicOrganizer(const LongRangeSystematicOrganizer& in);                 // Copy constructor
  ~LongRangeSystematicOrganizer();                                      // Destructor
  
  // Setter for input file
  void ReadInputFile(TFile *inputFile);               // Read the long range systematic uncertainty file
  
  // Getters for the long range systematic uncertainties
  TString GetLongRangeUncertaintyName(const int iUncertainty) const;
  TGraphErrors* GetLongRangeSystematicUncertainty(const int iFlow, const int iUncertainty = kAll, int iAsymmetry = knMaxXj) const;
  
private:
  
  TString fLongRangeUncertaintyName[knUncertaintySources] = {"deltaEtaRegion", "deltaEtaSide", "angleSmear", "jetCollection", "minBias", "all"};

  // Systematic uncertainty for long range correlations
  TGraphErrors * fLongRangeUncertaintyGraph[knUncertaintySources][knMaxFlow][knMaxXj+1];
  
};

#endif
