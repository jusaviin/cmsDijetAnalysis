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
  fhDijetAsymmetryVsDphi(0),
  fhDijetLeadingVsSubleadingPt(0),
  fCard(0)
{
  // Default constructor
}

/*
 * Default constructor
 */
DijetHistograms::DijetHistograms(ConfigurationCard *newCard) :
  fhVertexZ(0),
  fhEvents(0),
  fhCentrality(0),
  fhDijetDphi(0),
  fhDijetAsymmetry(0),
  fhLeadingJetPt(0),
  fhSubleadingJetPt(0),
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
  
  // Arrays for creating THnSpares
  int nBins2D[2]; int nBins3D[3];
  double lowBinBorder2D[2]; double lowBinBorder3D[3];
  double highBinBorder2D[2]; double highBinBorder3D[3];
  
  // ======== Plain TH1 histograms ========
  
  fhVertexZ = new TH1D("vertexZ","vertexZ",80,-20,20); fhVertexZ->Sumw2();
  fhEvents = new TH1D("nEvents","nEvents",10,0.5,10.5); fhEvents->Sumw2();
  fhCentrality = new TH1D("centrality","centrality",101,-1,100); fhCentrality->Sumw2();
  
  // ======== THnSparse histograms with centrality as the second axis =========
  
  // Fill in the centrality information common for all histograms
  nBins2D[1] = nCentralityBins;
  lowBinBorder2D[1] = minCentrality;
  highBinBorder2D[1] = maxCentrality;
  
  // Create the histograms
  // Dijet deltaPhi
  nBins2D[0] = 30;                  // nBins for deltaPhi
  lowBinBorder2D[0] = 0;            // low bin border for deltaPhi
  highBinBorder2D[0] = TMath::Pi(); // high bin border for deltaPhi
  fhDijetDphi = new THnSparseD("dijetDphi","dijetDphi",2,nBins2D,lowBinBorder2D,highBinBorder2D); fhDijetDphi->Sumw2();
  
  // Dijet astymmetry
  nBins2D[0] = 25;             // nBins for asymmetry
  lowBinBorder2D[0] = 0;       // low bin border for asymmetry
  highBinBorder2D[0] = 0.75;   // high bin border for asymmetry
  fhDijetAsymmetry = new THnSparseD("dijetAsymmetry","dijetAsymmetry",2,nBins2D,lowBinBorder2D,highBinBorder2D); fhDijetAsymmetry->Sumw2();
  
  // Leading jet pT
  nBins2D[0] = 75;          // nBins for leading jet pT
  lowBinBorder2D[0] = 0;    // low bin border for leading jet pT
  highBinBorder2D[0] = 150; // high bin border for leading jet pT
  fhLeadingJetPt = new THnSparseD("leadingJetPt","leadingJetPt",2,nBins2D,lowBinBorder2D,highBinBorder2D); fhLeadingJetPt->Sumw2();
  
  // Subleading jet pT
  nBins2D[0] = 75;          // nBins for subleading jet pT
  lowBinBorder2D[0] = 0;    // low bin border for subleading jet pT
  highBinBorder2D[0] = 150; // high bin border for subleading jet pT
  fhSubleadingJetPt = new THnSparseD("subleadingJetPt","subleadingJetPt",2,nBins2D,lowBinBorder2D,highBinBorder2D); fhSubleadingJetPt->Sumw2();
  
  // ========= THnSparse histograms with centrality as the third axis =========
  
  // Fill in the centrality information common for all histograms
  nBins3D[2] = nCentralityBins;
  lowBinBorder3D[2] = minCentrality;
  highBinBorder3D[2] = maxCentrality;
  
  // Dijet asymmetry versus delta phi
  nBins3D[0] = 30;                  // nBins for deltaPhi
  lowBinBorder3D[0] = 0;            // low bin border for deltaPhi
  highBinBorder3D[0] = TMath::Pi(); // high bin border for deltaPhi
  nBins3D[1] = 25;                  // nBins for asymmetry
  lowBinBorder3D[1] = 0;            // low bin border for asymmetry
  highBinBorder3D[1] = 0.75;        // high bin border for asymmetry
  fhDijetAsymmetryVsDphi = new THnSparseD("dijetAsymmetryVsDphi","dijetAsymmetryVsDphi",3,nBins3D,lowBinBorder3D,highBinBorder3D); fhDijetAsymmetryVsDphi->Sumw2();
  
  // Leading jet pT versus subleading jet pT
  nBins3D[0] = 75;           // nBins for leading jet pT
  lowBinBorder3D[0] = 0;     // low bin border for leading jet pT
  highBinBorder3D[0] = 150;  // high bin border for leading jet pT
  nBins3D[1] = 55;           // nBins for subleading jet pT
  lowBinBorder3D[1] = 0;     // low bin border for subleading jet pT
  highBinBorder3D[1] = 150;  // high bin border for subleading jet pT
  fhDijetLeadingVsSubleadingPt = new THnSparseD("dijetLeadingSubleadingPt","dijetLeadingSubleadingPt",3,nBins3D,lowBinBorder3D,highBinBorder3D); fhDijetLeadingVsSubleadingPt->Sumw2();
  
  
}

/*
 * Write the histograms to file
 */
void DijetHistograms::Write(){
  
  // Write the histograms to file
  fhVertexZ->Write();
  fhEvents->Write();
  fhCentrality->Write();
  fhDijetDphi->Write();
  fhDijetAsymmetry->Write();
  fhLeadingJetPt->Write();
  fhSubleadingJetPt->Write();
  fhDijetAsymmetryVsDphi->Write();
  fhDijetLeadingVsSubleadingPt->Write();
  
}

/*
 * Write the histograms to a given file
 */
void DijetHistograms::Write(TString outputFileName){
  
  // Define the output file
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  
  // Write the histograms to file
  Write();
  
  // Close the file after everything is written
  outputFile->Close();
}


