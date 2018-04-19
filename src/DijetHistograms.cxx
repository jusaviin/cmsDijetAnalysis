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
  fhDijetDphi(0),
  fhDijetAsymmetry(0),
  fhLeadingJetPt(0),
  fhSubleadingJetPt(0),
  fhAnyJetPt(0),
  fhLeadingJetPhi(0),
  fhSubleadingJetPhi(0),
  fhAnyJetPhi(0),
  fhLeadingJetEta(0),
  fhSubleadingJetEta(0),
  fhAnyJetEta(0),
  fhDijetAsymmetryVsDphi(0),
  fhDijetLeadingVsSubleadingPt(0),
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
  fhDijetDphi(0),
  fhDijetAsymmetry(0),
  fhLeadingJetPt(0),
  fhSubleadingJetPt(0),
  fhAnyJetPt(0),
  fhLeadingJetPhi(0),
  fhSubleadingJetPhi(0),
  fhAnyJetPhi(0),
  fhLeadingJetEta(0),
  fhSubleadingJetEta(0),
  fhAnyJetEta(0),
  fhDijetAsymmetryVsDphi(0),
  fhDijetLeadingVsSubleadingPt(0),
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
  fhDijetDphi(in.fhDijetDphi),
  fhDijetAsymmetry(in.fhDijetAsymmetry),
  fhLeadingJetPt(in.fhLeadingJetPt),
  fhSubleadingJetPt(in.fhSubleadingJetPt),
  fhAnyJetPt(in.fhAnyJetPt),
  fhLeadingJetPhi(in.fhLeadingJetPhi),
  fhSubleadingJetPhi(in.fhSubleadingJetPhi),
  fhAnyJetPhi(in.fhAnyJetPhi),
  fhLeadingJetEta(in.fhLeadingJetEta),
  fhSubleadingJetEta(in.fhSubleadingJetEta),
  fhAnyJetEta(in.fhAnyJetEta),
  fhDijetAsymmetryVsDphi(in.fhDijetAsymmetryVsDphi),
  fhDijetLeadingVsSubleadingPt(in.fhDijetLeadingVsSubleadingPt),
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
  fhDijetDphi = in.fhDijetDphi;
  fhDijetAsymmetry = in.fhDijetAsymmetry;
  fhLeadingJetPt = in.fhLeadingJetPt;
  fhSubleadingJetPt = in.fhSubleadingJetPt;
  fhAnyJetPt = in.fhAnyJetPt;
  fhLeadingJetPhi = in.fhLeadingJetPhi;
  fhSubleadingJetPhi = in.fhSubleadingJetPhi;
  fhAnyJetPhi = in.fhAnyJetPhi;
  fhLeadingJetEta = in.fhLeadingJetEta;
  fhSubleadingJetEta = in.fhSubleadingJetEta;
  fhAnyJetEta = in.fhAnyJetEta;
  fhDijetAsymmetryVsDphi = in.fhDijetAsymmetryVsDphi;
  fhDijetLeadingVsSubleadingPt = in.fhDijetLeadingVsSubleadingPt;
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
  delete fhDijetDphi;
  delete fhDijetAsymmetry;
  delete fhLeadingJetPt;
  delete fhSubleadingJetPt;
  delete fhAnyJetPt;
  delete fhLeadingJetPhi;
  delete fhSubleadingJetPhi;
  delete fhAnyJetPhi;
  delete fhLeadingJetEta;
  delete fhSubleadingJetEta;
  delete fhAnyJetEta;
  delete fhDijetAsymmetryVsDphi;
  delete fhDijetLeadingVsSubleadingPt;
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
  int nBins2D[2]; int nBins3D[3];
  double lowBinBorder2D[2]; double lowBinBorder3D[3];
  double highBinBorder2D[2]; double highBinBorder3D[3];
  
  // ======== Plain TH1 histograms ========
  
  fhVertexZ = new TH1D("vertexZ","vertexZ",nVzBins,minVz,maxVz); fhVertexZ->Sumw2();
  fhEvents = new TH1D("nEvents","nEvents",knEventTypes,-0.5,knEventTypes-0.5); fhEvents->Sumw2();
  fhCentrality = new TH1D("centrality","centrality",nCentralityBins,minCentrality,maxCentrality); fhCentrality->Sumw2();
  
  // For the event histogram, label each bin corresponding to an event cut
  for(int i = 0; i < knEventTypes; i++){
    fhEvents->GetXaxis()->SetBinLabel(i+1,kEventTypeStrings[i]);
  }
  
  // ======== THnSparse histograms with centrality as the second axis =========
  
  // Fill in the centrality information common for all histograms
  nBins2D[1] = nCentralityBins;       // nBins for centrality
  lowBinBorder2D[1] = minCentrality;  // low bin border for centrality
  highBinBorder2D[1] = maxCentrality; // high bin border for centrality
  
  // Create the histograms
  // Dijet deltaPhi
  nBins2D[0] = nDeltaPhiBins;       // nBins for deltaPhi
  lowBinBorder2D[0] = minDeltaPhi;  // low bin border for deltaPhi
  highBinBorder2D[0] = maxDeltaPhi; // high bin border for deltaPhi
  fhDijetDphi = new THnSparseD("dijetDphi","dijetDphi",2,nBins2D,lowBinBorder2D,highBinBorder2D); fhDijetDphi->Sumw2();
  
  // Dijet astymmetry
  nBins2D[0] = nAsymmetryBins;         // nBins for asymmetry
  lowBinBorder2D[0] = minAsymmetry;    // low bin border for asymmetry
  highBinBorder2D[0] = maxAsymmetry;   // high bin border for asymmetry
  fhDijetAsymmetry = new THnSparseD("dijetAsymmetry","dijetAsymmetry",2,nBins2D,lowBinBorder2D,highBinBorder2D); fhDijetAsymmetry->Sumw2();
  
  // Different jet pT histograms
  nBins2D[0] = nPtBins;         // nBins for jet pT
  lowBinBorder2D[0] = minPt;    // low bin border for jet pT
  highBinBorder2D[0] = maxPt;   // high bin border for jet pT
  
  // Leading jet pT
  fhLeadingJetPt = new THnSparseD("leadingJetPt","leadingJetPt",2,nBins2D,lowBinBorder2D,highBinBorder2D); fhLeadingJetPt->Sumw2();
  
  // Subleading jet pT
  fhSubleadingJetPt = new THnSparseD("subleadingJetPt","subleadingJetPt",2,nBins2D,lowBinBorder2D,highBinBorder2D); fhSubleadingJetPt->Sumw2();
  
  // All jets pT
  fhAnyJetPt = new THnSparseD("anyJetPt","anyJetPt",2,nBins2D,lowBinBorder2D,highBinBorder2D); fhAnyJetPt->Sumw2();
  
  // ======= THnSparse histograms with centrality as the second axis and pT as the third axis =======
  
  // Fill in the centrality information common for all histograms
  nBins3D[1] = nCentralityBins;        // nBins for centrality
  lowBinBorder3D[1] = minCentrality;   // low bin border for centrality
  highBinBorder3D[1] = maxCentrality;  // high bin border for centrality
  
  // Fill in the pT information common for all histograms
  nBins3D[2] = nPtBins;        // nBins for jet pT
  lowBinBorder3D[2] = minPt;   // low bin border for jet pT
  highBinBorder3D[2] = maxPt;  // high bin border for jet pT
  
  // Different jet phi histograms
  nBins3D[0] = nPhiBins;        // nBins for jet phi
  lowBinBorder3D[0] = minPhi;   // low bin border for jet phi
  highBinBorder3D[0] = maxPhi;  // high bin border for jet phi
  
  // Leading jet phi
  fhLeadingJetPhi = new THnSparseD("leadingJetPhi","leadingJetPhi",3,nBins3D,lowBinBorder3D,highBinBorder3D); fhLeadingJetPhi->Sumw2();
  
  // Subleading jet phi
  fhSubleadingJetPhi = new THnSparseD("subleadingJetPhi","subleadingJetPhi",3,nBins3D,lowBinBorder3D,highBinBorder3D); fhSubleadingJetPhi->Sumw2();
  
  // All jets phi
  fhAnyJetPhi = new THnSparseD("anyJetPhi","anyJetPhi",3,nBins3D,lowBinBorder3D,highBinBorder3D); fhAnyJetPhi->Sumw2();
  
  // Different jet eta histograms
  nBins3D[0] = nEtaBins;        // nBins for jet eta
  lowBinBorder3D[0] = minEta;   // low bin border for jet eta
  highBinBorder3D[0] = maxEta;  // high bin border for jet eta
  
  // Leading jet phi
  fhLeadingJetEta = new THnSparseD("leadingJetEta","leadingJetEta",3,nBins3D,lowBinBorder3D,highBinBorder3D); fhLeadingJetEta->Sumw2();
  
  // Subleading jet phi
  fhSubleadingJetEta = new THnSparseD("subleadingJetEta","subleadingJetEta",3,nBins3D,lowBinBorder3D,highBinBorder3D); fhSubleadingJetEta->Sumw2();
  
  // All jets phi
  fhAnyJetEta = new THnSparseD("anyJetEta","anyJetEta",3,nBins3D,lowBinBorder3D,highBinBorder3D); fhAnyJetEta->Sumw2();
  
  // ========= THnSparse histograms with centrality as the third axis =========
  
  // Fill in the centrality information common for all histograms
  nBins3D[2] = nCentralityBins;        // nBins for centrality
  lowBinBorder3D[2] = minCentrality;   // low bin border for centrality
  highBinBorder3D[2] = maxCentrality;  // high bin border for centrality
  
  // Dijet asymmetry versus delta phi
  nBins3D[0] = nDeltaPhiBins;         // nBins for deltaPhi
  lowBinBorder3D[0] = minDeltaPhi;    // low bin border for deltaPhi
  highBinBorder3D[0] = maxDeltaPhi;   // high bin border for deltaPhi
  nBins3D[1] = nAsymmetryBins;        // nBins for asymmetry
  lowBinBorder3D[1] = minAsymmetry;   // low bin border for asymmetry
  highBinBorder3D[1] = maxAsymmetry;  // high bin border for asymmetry
  fhDijetAsymmetryVsDphi = new THnSparseD("dijetAsymmetryVsDphi","dijetAsymmetryVsDphi",3,nBins3D,lowBinBorder3D,highBinBorder3D); fhDijetAsymmetryVsDphi->Sumw2();
  
  // Leading jet pT versus subleading jet pT
  nBins3D[0] = nPtBins;          // nBins for leading jet pT
  lowBinBorder3D[0] = minPt;     // low bin border for leading jet pT
  highBinBorder3D[0] = maxPt;    // high bin border for leading jet pT
  nBins3D[1] = nPtBins;          // nBins for subleading jet pT
  lowBinBorder3D[1] = minPt;     // low bin border for subleading jet pT
  highBinBorder3D[1] = maxPt;    // high bin border for subleading jet pT
  fhDijetLeadingVsSubleadingPt = new THnSparseD("dijetLeadingSubleadingPt","dijetLeadingSubleadingPt",3,nBins3D,lowBinBorder3D,highBinBorder3D); fhDijetLeadingVsSubleadingPt->Sumw2();
  
  
}

/*
 * Write the histograms to file
 */
void DijetHistograms::Write() const{
  
  // Write the histograms to file
  fhVertexZ->Write();
  fhEvents->Write();
  fhCentrality->Write();
  fhDijetDphi->Write();
  fhDijetAsymmetry->Write();
  fhLeadingJetPt->Write();
  fhSubleadingJetPt->Write();
  fhAnyJetPt->Write();
  fhLeadingJetPhi->Write();
  fhSubleadingJetPhi->Write();
  fhAnyJetPhi->Write();
  fhLeadingJetEta->Write();
  fhSubleadingJetEta->Write();
  fhAnyJetEta->Write();
  fhDijetAsymmetryVsDphi->Write();
  fhDijetLeadingVsSubleadingPt->Write();
  
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


