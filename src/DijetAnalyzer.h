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

class DijetAnalyzer{
  
public:
  
  // Constructors and destructor
  DijetAnalyzer(); // Default constructor
  DijetAnalyzer(std::vector<TString> fileNameVector, ConfigurationCard *newCard); // Custom constructor
  DijetAnalyzer(const DijetAnalyzer& in); // Copy constructor
  virtual ~DijetAnalyzer(); // Destructor
  DijetAnalyzer& operator=(const DijetAnalyzer& obj); // Equal sign operator
  
  // Methods
  void RunAnalysis(bool runLocal, int debugLevel);  // Run the dijet analysis
  DijetHistograms* GetHistograms();  // Getter for histograms
  
private:
  
  std::vector<TString> fFileNames;
  ConfigurationCard *fCard;
  DijetHistograms *fHistograms;

};

#endif
