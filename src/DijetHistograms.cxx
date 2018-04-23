// Histograms needed in the dijet analysis

// C++ includes
#include <assert.h>

// Root includes
#include <TFile.h>
#include <TMath.h>

// Own includes
#include "DijetHistograms.h"

/*
 * Default constructor
 */
DijetHistograms::DijetHistograms() :
  fhVertexZ(0),
  fhEvents(0),
  fhCentrality(0),
  fhDijet(0),
  fhAnyJet(0),
  fCard(0)
{
  // Default constructor
}

/*
 * Custom constructor
 */
DijetHistograms::DijetHistograms(ConfigurationCard *newCard) :
  fhVertexZ(0),
  fhEvents(0),
  fhCentrality(0),
  fhDijet(0),
  fhAnyJet(0),
  fCard(newCard)
{
  // Custom constructor
}

/*
 * Copy constructor
 */
DijetHistograms::DijetHistograms(const DijetHistograms& in) :
  fhVertexZ(in.fhVertexZ),
  fhEvents(in.fhEvents),
  fhCentrality(in.fhCentrality),
  fhDijet(in.fhDijet),
  fhAnyJet(in.fhAnyJet),
  fCard(in.fCard)
{
  // Copy constructor
}

/*
 * Assingment operator
 */
DijetHistograms& DijetHistograms::operator=(const DijetHistograms& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fhVertexZ = in.fhVertexZ;
  fhEvents = in.fhEvents;
  fhCentrality = in.fhCentrality;
  fhDijet = in.fhDijet;
  fhAnyJet = in.fhAnyJet;
  fCard = in.fCard;
  
  return *this;
}

/*
 * Destructor
 */
DijetHistograms::~DijetHistograms(){
  // destructor
  delete fhVertexZ;
  delete fhEvents;
  delete fhCentrality;
  delete fhDijet;
  delete fhAnyJet;
}

/*
 * Set the configuration card used for the histogram class
 */
void DijetHistograms::SetCard(ConfigurationCard *newCard){  
  fCard = newCard;
}

/*
 * Create the necessary histograms
 */
void DijetHistograms::CreateHistograms(){
  
  // Common binning information for histograms
  double minCentrality = -0.75;
  double maxCentrality = 100.25;
  int nCentralityBins = 202;
  
  double minPt = 0;
  double maxPt = 300;
  int nPtBins = 150;
  
  double minPhi = -TMath::Pi();
  double maxPhi = TMath::Pi();
  int nPhiBins = 72;
  
  double minEta = -2.0;
  double maxEta = 2.0;
  int nEtaBins = 40;
  
  double minDeltaPhi = 0;
  double maxDeltaPhi = TMath::Pi();
  int nDeltaPhiBins = 30;
  
  double minAsymmetry = 0;
  double maxAsymmetry = 0.75;
  int nAsymmetryBins = 25;
  
  double minVz = -20;
  double maxVz = 20;
  int nVzBins = 80;
  
  // Arrays for creating THnSparses
  const int dimensionDijet = 9;
  const int dimensionAnyJet = 4;
  int nBinsDijet[dimensionDijet]; int nBinsAnyjet[dimensionAnyJet];
  double lowBinBorderDijet[dimensionDijet]; double lowBinBorderAnyJet[dimensionAnyJet];
  double highBinBorderDijet[dimensionDijet]; double highBinBorderAnyJet[dimensionAnyJet];
  
  // ======== Plain TH1 histograms ========
  
  fhVertexZ = new TH1D("vertexZ","vertexZ",nVzBins,minVz,maxVz); fhVertexZ->Sumw2();
  fhEvents = new TH1D("nEvents","nEvents",knEventTypes,-0.5,knEventTypes-0.5); fhEvents->Sumw2();
  fhCentrality = new TH1D("centrality","centrality",nCentralityBins,minCentrality,maxCentrality); fhCentrality->Sumw2();
  
  // For the event histogram, label each bin corresponding to an event cut
  for(int i = 0; i < knEventTypes; i++){
    fhEvents->GetXaxis()->SetBinLabel(i+1,kEventTypeStrings[i]);
  }
  
  // ======== THnSparse histograms for dijets ========
  
  // Axis 0 for the dijet histogram: leading jet pT
  nBinsDijet[0] = nPtBins;         // nBins for leading jet pT
  lowBinBorderDijet[0] = minPt;    // low bin border for leading jet pT
  highBinBorderDijet[0] = maxPt;   // high bin border for leading jet pT
  
  // Axis 1 for the dijet histogram: leading jet phi
  nBinsDijet[1] = nPhiBins;        // nBins for leading jet phi
  lowBinBorderDijet[1] = minPhi;   // low bin border for leading jet phi
  highBinBorderDijet[1] = maxPhi;  // high bin border for leading jet phi
  
  // Axis 2 for the dijet histogram: leading jet eta
  nBinsDijet[2] = nEtaBins;        // nBins for leading jet eta
  lowBinBorderDijet[2] = minEta;   // low bin border for leading jet eta
  highBinBorderDijet[2] = maxEta;  // high bin border for leading jet eta
  
  // Axis 3 for the dijet histogram: subleading jet pT
  nBinsDijet[3] = nPtBins;         // nBins for subleading jet pT
  lowBinBorderDijet[3] = minPt;    // low bin border for subleading jet pT
  highBinBorderDijet[3] = maxPt;   // high bin border for subleading jet pT
  
  // Axis 4 for the dijet histogram: subleading jet phi
  nBinsDijet[4] = nPhiBins;        // nBins for subleading jet phi
  lowBinBorderDijet[4] = minPhi;   // low bin border for subleading jet phi
  highBinBorderDijet[4] = maxPhi;  // high bin border for subleading jet phi
  
  // Axis 5 for the dijet histogram: subleading jet eta
  nBinsDijet[5] = nEtaBins;        // nBins for subleading jet eta
  lowBinBorderDijet[5] = minEta;   // low bin border for subleading jet eta
  highBinBorderDijet[5] = maxEta;  // high bin border for subleading jet eta
  
  // Axis 6 for the dijet histogram: deltaPhi
  nBinsDijet[6] = nDeltaPhiBins;       // nBins for deltaPhi
  lowBinBorderDijet[6] = minDeltaPhi;  // low bin border for deltaPhi
  highBinBorderDijet[6] = maxDeltaPhi; // high bin border for deltaPhi
  
  // Axis 7 for the dijet histogram: asymmetry
  nBinsDijet[7] = nAsymmetryBins;         // nBins for asymmetry
  lowBinBorderDijet[7] = minAsymmetry;    // low bin border for asymmetry
  highBinBorderDijet[7] = maxAsymmetry;   // high bin border for asymmetry
  
  // Axis 8 for the dijet histogram: centrality
  nBinsDijet[8] = nCentralityBins;       // nBins for centrality
  lowBinBorderDijet[8] = minCentrality;  // low bin border for centrality
  highBinBorderDijet[8] = maxCentrality; // high bin border for centrality
  
  // Create the dijet histogram using the above binning information
  fhDijet = new THnSparseD("dijet","dijet",dimensionDijet,nBinsDijet,lowBinBorderDijet,highBinBorderDijet); fhDijet->Sumw2();
  
  // ======== THnSparse histograms for all jets ========
  
  // Axis 0 for the dijet histogram: leading jet pT
  nBinsAnyjet[0] = nPtBins;         // nBins for any jet pT
  lowBinBorderAnyJet[0] = minPt;    // low bin border for any jet pT
  highBinBorderAnyJet[0] = maxPt;   // high bin border for any jet pT
  
  // Axis 1 for the dijet histogram: leading jet phi
  nBinsAnyjet[1] = nPhiBins;        // nBins for any jet phi
  lowBinBorderAnyJet[1] = minPhi;   // low bin border for any jet phi
  highBinBorderAnyJet[1] = maxPhi;  // high bin border for any jet phi
  
  // Axis 2 for the dijet histogram: leading jet eta
  nBinsAnyjet[2] = nEtaBins;        // nBins for any jet eta
  lowBinBorderAnyJet[2] = minEta;   // low bin border for any jet eta
  highBinBorderAnyJet[2] = maxEta;  // high bin border for any jet eta
  
  // Axis 3 for the dijet histogram: centrality
  nBinsAnyjet[3] = nCentralityBins;       // nBins for centrality
  lowBinBorderAnyJet[3] = minCentrality;  // low bin border for centrality
  highBinBorderAnyJet[3] = maxCentrality; // high bin border for centrality
  
  // Create the histogram for all jets using the above binning information
  fhAnyJet = new THnSparseD("anyJet","anyJet",dimensionAnyJet,nBinsAnyjet,lowBinBorderAnyJet,highBinBorderAnyJet); fhAnyJet->Sumw2();
  
}

/*
 * Write the histograms to file
 */
void DijetHistograms::Write() const{
  
  // Write the histograms to file
  fhVertexZ->Write();
  fhEvents->Write();
  fhCentrality->Write();
  fhDijet->Write();
  fhAnyJet->Write();
  
}

/*
 * Write the histograms to a given file
 */
void DijetHistograms::Write(TString outputFileName) const{
  
  // Define the output file
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  
  // Write the histograms to file
  Write();
  
  // Close the file after everything is written
  outputFile->Close();
}


