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
  TH1D *fhVertexZ;                          // Vertex z-position
  TH1D *fhEvents;                           // Number of events. Bin 1 = all. Bin 2 = good vz bin. Bin 3 = Dijet found.
  TH1D *fhCentrality;                       // Centrality information. -0.5 for pp or PYTHIA.
  THnSparseD *fhDijetDphi;                  // deltaPhi for dijet events  Axes: [dPhi][cent]
  THnSparseD *fhDijetAsymmetry;             // Asymmetry in dijet events  Axes: [Ajj][cent]
  THnSparseD *fhLeadingJetPt;               // pT for the leading jet     Axes: [pT][cent]
  THnSparseD *fhSubleadingJetPt;            // pT for the subleading jet  Axes: [pT][cent]
  THnSparseD *fhAnyJetPt;                   // pT for all jets            Axes: [pT][cent]
  THnSparseD *fhLeadingJetPhi;              // phi for the leading jet    Axes: [phi][cent][pT]
  THnSparseD *fhSubleadingJetPhi;           // phi for the subleading jet Axes: [phi][cent][pT]
  THnSparseD *fhAnyJetPhi;                  // phi for all jets           Axes: [phi][cent][pT]
  THnSparseD *fhLeadingJetEta;              // eta for the leading jet    Axes: [eta][cent][pT]
  THnSparseD *fhSubleadingJetEta;           // eta for the subleading jet Axes: [eta][cent][pT]
  THnSparseD *fhAnyJetEta;                  // eta for all jets           Axes: [eta][cent][pT]
  THnSparseD *fhDijetAsymmetryVsDphi;       // Asymmetry vs. deltaPhi     Axes: [Ajj][dPhi][cent]
  THnSparseD *fhDijetLeadingVsSubleadingPt; // Leading jet pT vs. subleading jet pT Axes: [pT][pT][cent]

  
private:
  
  ConfigurationCard *fCard;    // Card for binning info
  const TString kEventTypeStrings[knEventTypes] = {"All", "PrimVertex", "HBHENoise", "CollEvtSel", "BeamScrape", "CaloJet", "v_{z} cut", "Dijet"}; // Strings corresponding to event types
  
};

#endif
