// Class for histograms needed in the dijet analysis

#ifndef DIJETHISTOGRAMS_H
#define DIJETHISTOGRAMS_H

// Root includes
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>

// Own includes
#include "ConfigurationCard.h"

class DijetHistograms{
  
public:
  
  // Enumeration for event types to event histogram
  enum enumEventTypes {kAll, kPrimaryVertex, kHBHENoise, kCollisionEventSelection, kBeamScraping, kCaloJet, kVzCut, kDijet, knEventTypes};
  
  // Constructors and destructor
  DijetHistograms(); // Default constructor
  DijetHistograms(ConfigurationCard *newCard); // Custom constructor
  DijetHistograms(const DijetHistograms& in); // Copy constructor
  virtual ~DijetHistograms(); // Destructor
  DijetHistograms& operator=(const DijetHistograms& obj); // Equal sign operator
  
  // Methods
  void CreateHistograms();                   // Create all histograms
  void Write() const;                        // Write the histograms to a file that is opened somewhere else
  void Write(TString outputFileName) const;  // Write the histograms to a file
  void SetCard(ConfigurationCard *newCard);  // Set a new configuration card for the histogram class
  
  // Histograms defined public to allow easier access to them. Should not be abused
  // Notation in comments: l = leading jet, s = subleading jet, uc = uncorrected, ptw = pT weighted
  TH1D *fhVertexZ;      // Vertex z-position
  TH1D *fhEvents;       // Number of events. Bin 1 = all. Bin 2 = good vz bin. Bin 3 = Dijet found.
  TH1D *fhCentrality;   // Centrality information. -0.5 for pp or PYTHIA.
  THnSparse *fhDijet;   // Dijet information. Axes: [l-pT][l-phi][l-eta][s-pT][s-phi][s-eta][dphi][Ajj][cent]
  THnSparse *fhAnyJet;  // Any jet information. Axes: [jet pT][jet phi][jet eta][cent]
  THnSparse *fhTrack;   // Track histogram. Axes: [pT][phi][eta][uc pT][uc phi][uc eta][l-dphi][uc l-dhi][ptw l-dphi][l-deta][uc l-eta][ptw l-deta][s-dphi][uc s-dhi][ptw s-dphi][s-deta][uc s-deta][ptw s-deta][Ajj][cent]
  
private:
  
  ConfigurationCard *fCard;    // Card for binning info
  const TString kEventTypeStrings[knEventTypes] = {"All", "PrimVertex", "HBHENoise", "CollEvtSel", "BeamScrape", "CaloJet", "v_{z} cut", "Dijet"}; // Strings corresponding to event types
  
};

#endif
