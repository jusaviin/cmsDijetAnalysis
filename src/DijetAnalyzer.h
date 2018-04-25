// Class for the main analysis algorithms for leading-subleading jet analysis

#ifndef DIJETANALYZER_H
#define DIJETANALYZER_H

// C++ includes
#include <vector>

// Root includes
#include <TString.h>

// Own includes
#include "ConfigurationCard.h"
#include "DijetHistograms.h"
#include "TrkCorr.h"

class DijetAnalyzer{
  
public:
  
  // Constructors and destructor
  DijetAnalyzer(); // Default constructor
  DijetAnalyzer(std::vector<TString> fileNameVector, ConfigurationCard *newCard); // Custom constructor
  DijetAnalyzer(const DijetAnalyzer& in); // Copy constructor
  virtual ~DijetAnalyzer(); // Destructor
  DijetAnalyzer& operator=(const DijetAnalyzer& obj); // Equal sign operator
  
  // Methods
  void RunAnalysis();                      // Run the dijet analysis
  DijetHistograms* GetHistograms() const;  // Getter for histograms
  
private:
  
  std::vector<TString> fFileNames;   // Vector for all the files to loop over
  ConfigurationCard *fCard;          // Configuration card for the analysis
  DijetHistograms *fHistograms;      // Filled histograms
  TrkCorr *fTrackCorrection;         // Track correction class

};

#endif
