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
  fhTrackUncorrected(0),
  fhTrackPtWeighted(0),
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
  fhTrackCuts(0),
  fhCentrality(0),
  fhDijet(0),
  fhAnyJet(0),
  fhTrack(0),
  fhTrackUncorrected(0),
  fhTrackPtWeighted(0),
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
  fhTrackCuts(in.fhTrackCuts),
  fhCentrality(in.fhCentrality),
  fhDijet(in.fhDijet),
  fhAnyJet(in.fhAnyJet),
  fhTrack(in.fhTrack),
  fhTrackUncorrected(in.fhTrackUncorrected),
  fhTrackPtWeighted(in.fhTrackPtWeighted),
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
  fhTrackCuts = in.fhTrackCuts;
  fhCentrality = in.fhCentrality;
  fhDijet = in.fhDijet;
  fhAnyJet = in.fhAnyJet;
  fhTrack = in.fhTrack;
  fhTrackUncorrected = in.fhTrackUncorrected;
  fhTrackPtWeighted = in.fhTrackPtWeighted;
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
  delete fhTrackCuts;
  delete fhCentrality;
  delete fhDijet;
  delete fhAnyJet;
  delete fhTrack;
  delete fhTrackUncorrected;
  delete fhTrackPtWeighted;
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
  const int dimensionTrack = 9;
  const int dimensionTrackPtWeight = 6;
  int nBinsDijet[dimensionDijet];
  int nBinsAnyjet[dimensionAnyJet];
  int nBinsTrack[dimensionTrack];
  int nBinsTrackPtWeight[dimensionTrackPtWeight];
  double lowBinBorderDijet[dimensionDijet];
  double lowBinBorderAnyJet[dimensionAnyJet];
  double lowBinBorderTrack[dimensionTrack];
  double lowBinBorderTrackPtWeight[dimensionTrackPtWeight];
  double highBinBorderDijet[dimensionDijet];
  double highBinBorderAnyJet[dimensionAnyJet];
  double highBinBorderTrack[dimensionTrack];
  double highBinBorderTrackPtWeight[dimensionTrackPtWeight];
  
  // ======== Plain TH1 histograms ========
  
  fhVertexZ = new TH1D("vertexZ","vertexZ",nVzBins,minVz,maxVz); fhVertexZ->Sumw2();
  fhEvents = new TH1D("nEvents","nEvents",knEventTypes,-0.5,knEventTypes-0.5); fhEvents->Sumw2();
  fhTrackCuts = new TH1D("trackCuts","trackCuts",knTrackCuts,-0.5,knTrackCuts-0.5); fhTrackCuts->Sumw2();
  fhCentrality = new TH1D("centrality","centrality",nCentralityBins,minCentrality,maxCentrality); fhCentrality->Sumw2();
  
  // For the event histogram, label each bin corresponding to an event cut
  for(int i = 0; i < knEventTypes; i++){
    fhEvents->GetXaxis()->SetBinLabel(i+1,kEventTypeStrings[i]);
  }
  
  // For the track cut histogram, label each bin corresponding to a track cut
  for(int i = 0; i < knTrackCuts; i++){
    fhTrackCuts->GetXaxis()->SetBinLabel(i+1,kTrackCutStrings[i]);
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
  
  // Axis 3 for the track histogram: deltaPhi between track and leading jet
  nBinsTrack[3] = nDeltaPhiBinsTrack;       // nBins for deltaPhi between track and leading jet
  lowBinBorderTrack[3] = minDeltaPhiTrack;  // low bin border for deltaPhi between track and leading jet
  highBinBorderTrack[3] = maxDeltaPhiTrack; // high bin border for deltaPhi between track and leading jet
  
  // Axis 4 for the track histogram: deltaEta between track and leading jet
  nBinsTrack[4] = nDeltaEtaBins;         // nBins for deltaEta between track and leading jet
  lowBinBorderTrack[4] = minDeltaEta;    // low bin border deltaEta between track and leading jet
  highBinBorderTrack[4] = maxDeltaEta;   // high bin border deltaEta between track and leading jet
  
  // Axis 5 for the track histogram: deltaPhi between track and subleading jet
  nBinsTrack[5] = nDeltaPhiBinsTrack;       // nBins for deltaPhi between track and subleading jet
  lowBinBorderTrack[5] = minDeltaPhiTrack;  // low bin border for deltaPhi between track and subleading jet
  highBinBorderTrack[5] = maxDeltaPhiTrack; // high bin border for deltaPhi between track and subleading jet
  
  // Axis 6 for the track histogram: deltaEta between track and subleading jet
  nBinsTrack[6] = nDeltaEtaBins;         // nBins for deltaEta between track and subleading jet
  lowBinBorderTrack[6] = minDeltaEta;    // low bin border deltaEta between track and subleading jet
  highBinBorderTrack[6] = maxDeltaEta;   // high bin border deltaEta between track and subleading jet
  
  // Axis 7 for the track histogram: dijet asymmetry
  nBinsTrack[7] = nAsymmetryBins;         // nBins for dijet asymmetry
  lowBinBorderTrack[7] = minAsymmetry;    // low bin border for dijet asymmetry
  highBinBorderTrack[7] = maxAsymmetry;   // high bin border for dijet asymmetry
  
  // Axis 8 for the track histogram: centrality
  nBinsTrack[8] = nCentralityBins;         // nBins for centrality
  lowBinBorderTrack[8] = minCentrality;    // low bin border for centrality
  highBinBorderTrack[8] = maxCentrality;   // high bin border for centrality
  
  // Create histograms for tracks and uncorrected tracks using the above binning information
  fhTrack = new THnSparseD("track","track",dimensionTrack,nBinsTrack,lowBinBorderTrack,highBinBorderTrack); fhTrack->Sumw2();
  fhTrackUncorrected = new THnSparseD("trackUncorrected","trackUncorrected",dimensionTrack,nBinsTrack,lowBinBorderTrack,highBinBorderTrack); fhTrackUncorrected->Sumw2();
  
  // ======== THnSparse for pT weighted jet-track correlations with leading and subleading jets ========
  
  // Axis 0 for the track histogram: deltaPhi between track and leading jet
  nBinsTrackPtWeight[0] = nDeltaPhiBinsTrack;       // nBins for deltaPhi between track and leading jet
  lowBinBorderTrackPtWeight[0] = minDeltaPhiTrack;  // low bin border for deltaPhi between track and leading jet
  highBinBorderTrackPtWeight[0] = maxDeltaPhiTrack; // high bin border for deltaPhi between track and leading jet
  
  // Axis 1 for the track histogram: deltaEta between track and leading jet
  nBinsTrackPtWeight[1] = nDeltaEtaBins;         // nBins for deltaEta between track and leading jet
  lowBinBorderTrackPtWeight[1] = minDeltaEta;    // low bin border deltaEta between track and leading jet
  highBinBorderTrackPtWeight[1] = maxDeltaEta;   // high bin border deltaEta between track and leading jet
  
  // Axis 2 for the track histogram: deltaPhi between track and subleading jet
  nBinsTrackPtWeight[2] = nDeltaPhiBinsTrack;       // nBins for deltaPhi between track and subleading jet
  lowBinBorderTrackPtWeight[2] = minDeltaPhiTrack;  // low bin border for deltaPhi between track and subleading jet
  highBinBorderTrackPtWeight[2] = maxDeltaPhiTrack; // high bin border for deltaPhi between track and subleading jet
  
  // Axis 3 for the track histogram: deltaEta between track and subleading jet
  nBinsTrackPtWeight[3] = nDeltaEtaBins;         // nBins for deltaEta between track and subleading jet
  lowBinBorderTrackPtWeight[3] = minDeltaEta;    // low bin border deltaEta between track and subleading jet
  highBinBorderTrackPtWeight[3] = maxDeltaEta;   // high bin border deltaEta between track and subleading jet
  
  // Axis 4 for the track histogram: dijet asymmetry
  nBinsTrackPtWeight[4] = nAsymmetryBins;         // nBins for dijet asymmetry
  lowBinBorderTrackPtWeight[4] = minAsymmetry;    // low bin border for dijet asymmetry
  highBinBorderTrackPtWeight[4] = maxAsymmetry;   // high bin border for dijet asymmetry
  
  // Axis 5 for the track histogram: centrality
  nBinsTrackPtWeight[5] = nCentralityBins;         // nBins for centrality
  lowBinBorderTrackPtWeight[5] = minCentrality;    // low bin border for centrality
  highBinBorderTrackPtWeight[5] = maxCentrality;   // high bin border for centrality
  
  // Create histograms for pT weighted tracks using the above binning information
  fhTrackPtWeighted = new THnSparseD("trackPtWeighted","trackPtWeighted",dimensionTrackPtWeight,nBinsTrackPtWeight,lowBinBorderTrackPtWeight,highBinBorderTrackPtWeight); fhTrackPtWeighted->Sumw2();
}

/*
 * Write the histograms to file
 */
void DijetHistograms::Write() const{
  
  // Write the histograms to file
  fhVertexZ->Write();
  fhEvents->Write();
  fhTrackCuts->Write();
  fhCentrality->Write();
  fhDijet->Write();
  fhAnyJet->Write();
  fhTrack->Write();
  fhTrackUncorrected->Write();
  fhTrackPtWeighted->Write();
  
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


