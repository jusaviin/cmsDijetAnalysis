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
  fhTrack(0),
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
  fhTrack(0),
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
  fhTrack(in.fhTrack),
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
  fhTrack = in.fhTrack;
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
  delete fhTrack;
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
  
  // ======== Common binning information for histograms =========
  
  // Centrality
  double minCentrality = -0.75;   // Minimum centrality bin, is negative since hiBin is -1 for pp
  double maxCentrality = 100.25;  // Maximum centrality bin
  int nCentralityBins = 202;      // Number of centrality bins
  
  // Jet pT
  double minPtJet = 0;     // Minimum jet pT
  double maxPtJet = 300;   // Maximum jet pT
  int nPtBinsJet = 150;    // Number of jet pT bins
  
  //Track pT
  double minPtTrack = 0;   // Minimum track pT
  double maxPtTrack = 30;  // Maximum track pT   (Hallie's analysis = 20)
  int nPtBinsTrack = 600;  // Number of track pT bins (Hallie's analysis = 500)
  
  // Phi
  double minPhi = -TMath::Pi();  // Minimum phi
  double maxPhi = TMath::Pi();   // Maximum phi
  int nPhiBins = 72;             // Number of phi bins
  
  // Eta
  double minEta = -2.5;    // Minimum eta (current eta cut for tracks = 2.4)
  double maxEta = 2.5;     // Maximum eta (current eta cut for tracks = 2.4)
  int nEtaBins = 50;       // Number of eta bins
  
  // DeltaPhi in [0,pi]
  double minDeltaPhi = 0;             // Minimum deltaPhi
  double maxDeltaPhi = TMath::Pi();   // Maximum deltaPhi
  int nDeltaPhiBins = 30;             // Number of deltaPhi bins
  
  // DeltaPhi in [-pi/2,3pi/2]
  double minDeltaPhiTrack = -TMath::Pi()/2.0;    // Minimum deltaPhi for two dimensional plots
  double maxDeltaPhiTrack = 3.0*TMath::Pi()/2.0; // Maximum deltaPhi for two dimensional plots
  int nDeltaPhiBinsTrack = 60;                   // Number of deltaPhi bins for two dimensional plots
  
  // DeltaEta
  double minDeltaEta = -5.0;   // Minimum deltaEta
  double maxDeltaEta = 5.0;    // Maximum deltaEta
  int nDeltaEtaBins = 100;     // Number of deltaEta bins
  
  // Dijet asymmetry
  double minAsymmetry = 0;     // Minimum asymmetry
  double maxAsymmetry = 0.75;  // Maximum asymmetry
  int nAsymmetryBins = 25;     // Number of asymmetry bins
  
  // Vertex z-position
  double minVz = -20;   // Minimum vz
  double maxVz = 20;    // Maximum vz
  int nVzBins = 80;     // Number of vz bins
  
  // Arrays for creating THnSparses
  const int dimensionDijet = 9;
  const int dimensionAnyJet = 4;
  const int dimensionTrack = 20;
  int nBinsDijet[dimensionDijet];
  int nBinsAnyjet[dimensionAnyJet];
  int nBinsTrack[dimensionTrack];
  double lowBinBorderDijet[dimensionDijet];
  double lowBinBorderAnyJet[dimensionAnyJet];
  double lowBinBorderTrack[dimensionTrack];
  double highBinBorderDijet[dimensionDijet];
  double highBinBorderAnyJet[dimensionAnyJet];
  double highBinBorderTrack[dimensionTrack];
  
  // ======== Plain TH1 histograms ========
  
  fhVertexZ = new TH1D("vertexZ","vertexZ",nVzBins,minVz,maxVz); fhVertexZ->Sumw2();
  fhEvents = new TH1D("nEvents","nEvents",knEventTypes,-0.5,knEventTypes-0.5); fhEvents->Sumw2();
  fhCentrality = new TH1D("centrality","centrality",nCentralityBins,minCentrality,maxCentrality); fhCentrality->Sumw2();
  
  // For the event histogram, label each bin corresponding to an event cut
  for(int i = 0; i < knEventTypes; i++){
    fhEvents->GetXaxis()->SetBinLabel(i+1,kEventTypeStrings[i]);
  }
  
  // ======== THnSparse for dijets ========
  
  // Axis 0 for the dijet histogram: leading jet pT
  nBinsDijet[0] = nPtBinsJet;         // nBins for leading jet pT
  lowBinBorderDijet[0] = minPtJet;    // low bin border for leading jet pT
  highBinBorderDijet[0] = maxPtJet;   // high bin border for leading jet pT
  
  // Axis 1 for the dijet histogram: leading jet phi
  nBinsDijet[1] = nPhiBins;        // nBins for leading jet phi
  lowBinBorderDijet[1] = minPhi;   // low bin border for leading jet phi
  highBinBorderDijet[1] = maxPhi;  // high bin border for leading jet phi
  
  // Axis 2 for the dijet histogram: leading jet eta
  nBinsDijet[2] = nEtaBins;        // nBins for leading jet eta
  lowBinBorderDijet[2] = minEta;   // low bin border for leading jet eta
  highBinBorderDijet[2] = maxEta;  // high bin border for leading jet eta
  
  // Axis 3 for the dijet histogram: subleading jet pT
  nBinsDijet[3] = nPtBinsJet;         // nBins for subleading jet pT
  lowBinBorderDijet[3] = minPtJet;    // low bin border for subleading jet pT
  highBinBorderDijet[3] = maxPtJet;   // high bin border for subleading jet pT
  
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
  
  // ======== THnSparse for all jets ========
  
  // Axis 0 for the dijet histogram: leading jet pT
  nBinsAnyjet[0] = nPtBinsJet;         // nBins for any jet pT
  lowBinBorderAnyJet[0] = minPtJet;    // low bin border for any jet pT
  highBinBorderAnyJet[0] = maxPtJet;   // high bin border for any jet pT
  
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
  
  // ======== THnSparse for tracks and jet-track correlations with leading and subleading jets ========
  
  // Axis 0 for the track histogram: track pT
  nBinsTrack[0] = nPtBinsTrack;         // nBins for track pT
  lowBinBorderTrack[0] = minPtTrack;    // low bin border for track pT
  highBinBorderTrack[0] = maxPtTrack;   // high bin border for track pT
  
  // Axis 1 for the track histogram: track phi
  nBinsTrack[1] = nPhiBins;         // nBins for track phi
  lowBinBorderTrack[1] = minPhi;    // low bin border for track phi
  highBinBorderTrack[1] = maxPhi;   // high bin border for track phi
  
  // Axis 2 for the track histogram: track eta
  nBinsTrack[2] = nEtaBins;         // nBins for track eta
  lowBinBorderTrack[2] = minEta;    // low bin border for track eta
  highBinBorderTrack[2] = maxEta;   // high bin border for track eta
  
  // Axis 3 for the track histogram: uncorrected track pT
  nBinsTrack[3] = nPtBinsTrack;         // nBins for uncorrected track pT
  lowBinBorderTrack[3] = minPtTrack;    // low bin border for uncorrected track pT
  highBinBorderTrack[3] = maxPtTrack;   // high bin border for uncorrected track pT
  
  // Axis 4 for the track histogram: uncorrected track phi
  nBinsTrack[4] = nPhiBins;         // nBins for uncorrected track phi
  lowBinBorderTrack[4] = minPhi;    // low bin border for uncorrected track phi
  highBinBorderTrack[4] = maxPhi;   // high bin border for uncorrected track phi
  
  // Axis 5 for the track histogram: uncorrected track eta
  nBinsTrack[5] = nEtaBins;         // nBins for uncorrected track eta
  lowBinBorderTrack[5] = minEta;    // low bin border for uncorrected track eta
  highBinBorderTrack[5] = maxEta;   // high bin border for uncorrected track eta
  
  // Axis 6 for the track histogram: deltaPhi between track and leading jet
  nBinsTrack[6] = nDeltaPhiBinsTrack;       // nBins for deltaPhi between track and leading jet
  lowBinBorderTrack[6] = minDeltaPhiTrack;  // low bin border for deltaPhi between track and leading jet
  highBinBorderTrack[6] = maxDeltaPhiTrack; // high bin border for deltaPhi between track and leading jet
  
  // Axis 7 for the track histogram: deltaPhi between uncorrected track and leading jet
  nBinsTrack[7] = nDeltaPhiBinsTrack;       // nBins for deltaPhi between uncorrected track and leading jet
  lowBinBorderTrack[7] = minDeltaPhiTrack;  // low bin border for deltaPhi between uncorrected track and leading jet
  highBinBorderTrack[7] = maxDeltaPhiTrack; // high bin border for deltaPhi between uncorrected track and leading jet
  
  // Axis 8 for the track histogram: deltaPhi between pT weighted track and leading jet
  nBinsTrack[8] = nDeltaPhiBinsTrack;       // nBins for deltaPhi between pT weighted track and leading jet
  lowBinBorderTrack[8] = minDeltaPhiTrack;  // low bin border for deltaPhi between pT weighted track and leading jet
  highBinBorderTrack[8] = maxDeltaPhiTrack; // high bin border for deltaPhi between pT weighted track and leading jet
  
  // Axis 9 for the track histogram: deltaEta between track and leading jet
  nBinsTrack[9] = nDeltaEtaBins;         // nBins for deltaEta between track and leading jet
  lowBinBorderTrack[9] = minDeltaEta;    // low bin border deltaEta between track and leading jet
  highBinBorderTrack[9] = maxDeltaEta;   // high bin border deltaEta between track and leading jet
  
  // Axis 10 for the track histogram: deltaEta between uncorrected track and leading jet
  nBinsTrack[10] = nDeltaEtaBins;         // nBins for deltaEta between uncorrected track and leading jet
  lowBinBorderTrack[10] = minDeltaEta;    // low bin border for deltaEta between uncorrected track and leading jet
  highBinBorderTrack[10] = maxDeltaEta;   // high bin border for deltaEta between uncorrected track and leading jet
  
  // Axis 11 for the track histogram: deltaEta between pT weighted track and leading jet
  nBinsTrack[11] = nDeltaEtaBins;         // nBins for deltaEta between pT weighted track and leading jet
  lowBinBorderTrack[11] = minDeltaEta;    // low bin border for deltaEta between pT weighted track and leading jet
  highBinBorderTrack[11] = maxDeltaEta;   // high bin border for deltaEta between pT weighted track and leading jet
  
  // Axis 12 for the track histogram: deltaPhi between track and subleading jet
  nBinsTrack[12] = nDeltaPhiBinsTrack;       // nBins for deltaPhi between track and subleading jet
  lowBinBorderTrack[12] = minDeltaPhiTrack;  // low bin border for deltaPhi between track and subleading jet
  highBinBorderTrack[12] = maxDeltaPhiTrack; // high bin border for deltaPhi between track and subleading jet
  
  // Axis 13 for the track histogram: deltaPhi between uncorrected track and subleading jet
  nBinsTrack[13] = nDeltaPhiBinsTrack;       // nBins for deltaPhi between uncorrected track and subleading jet
  lowBinBorderTrack[13] = minDeltaPhiTrack;  // low bin border for deltaPhi between uncorrected track and subleading jet
  highBinBorderTrack[13] = maxDeltaPhiTrack; // high bin border for deltaPhi between uncorrected track and subleading jet
  
  // Axis 14 for the track histogram: deltaPhi between pT weighted track and subleading jet
  nBinsTrack[14] = nDeltaPhiBinsTrack;       // nBins for deltaPhi between pT weighted track and subleading jet
  lowBinBorderTrack[14] = minDeltaPhiTrack;  // low bin border for deltaPhi between pT weighted track and subleading jet
  highBinBorderTrack[14] = maxDeltaPhiTrack; // high bin border for deltaPhi between pT weighted track and subleading jet
  
  // Axis 15 for the track histogram: deltaEta between track and subleading jet
  nBinsTrack[15] = nDeltaEtaBins;         // nBins for deltaEta between track and subleading jet
  lowBinBorderTrack[15] = minDeltaEta;    // low bin border deltaEta between track and subleading jet
  highBinBorderTrack[15] = maxDeltaEta;   // high bin border deltaEta between track and subleading jet
  
  // Axis 16 for the track histogram: deltaEta between uncorrected track and subleading jet
  nBinsTrack[16] = nDeltaEtaBins;         // nBins for deltaEta between uncorrected track and subleading jet
  lowBinBorderTrack[16] = minDeltaEta;    // low bin border for deltaEta between uncorrected track and subleading jet
  highBinBorderTrack[16] = maxDeltaEta;   // high bin border for deltaEta between uncorrected track and subleading jet
  
  // Axis 17 for the track histogram: deltaEta between pT weighted track and subleading jet
  nBinsTrack[17] = nDeltaEtaBins;         // nBins for deltaEta between pT weighted track and subleading jet
  lowBinBorderTrack[17] = minDeltaEta;    // low bin border for deltaEta between pT weighted track and subleading jet
  highBinBorderTrack[17] = maxDeltaEta;   // high bin border for deltaEta between pT weighted track and subleading jet
  
  // Axis 18 for the track histogram: dijet asymmetry
  nBinsTrack[18] = nAsymmetryBins;         // nBins for dijet asymmetry
  lowBinBorderTrack[18] = minAsymmetry;    // low bin border for dijet asymmetry
  highBinBorderTrack[18] = maxAsymmetry;   // high bin border for dijet asymmetry
  
  // Axis 19 for the track histogram: centrality
  nBinsTrack[19] = nCentralityBins;         // nBins for centrality
  lowBinBorderTrack[19] = minCentrality;    // low bin border for centrality
  highBinBorderTrack[19] = maxCentrality;   // high bin border for centrality
  
  // Create the histogram for tracks using the above binning information
  fhTrack = new THnSparseD("track","track",dimensionTrack,nBinsTrack,lowBinBorderTrack,highBinBorderTrack); fhTrack->Sumw2();
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
  fhTrack->Write();
  
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


