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
  Double_t minCentrality = -0.75;   // Minimum centrality bin, is negative since hiBin is -1 for pp
  Double_t maxCentrality = 100.25;  // Maximum centrality bin
  Int_t nCentralityBins = 202;      // Number of centrality bins
  
  // Jet pT
  Double_t minPtJet = 0;     // Minimum jet pT
  Double_t maxPtJet = 300;   // Maximum jet pT
  Int_t nPtBinsJet = 150;    // Number of jet pT bins
  
  //Track pT
  Double_t minPtTrack = 0;   // Minimum track pT
  Double_t maxPtTrack = 30;  // Maximum track pT   (Hallie's analysis = 20)
  Int_t nPtBinsTrack = 600;  // Number of track pT bins (Hallie's analysis = 500)
  
  // Phi
  Double_t minPhi = -TMath::Pi();  // Minimum phi
  Double_t maxPhi = TMath::Pi();   // Maximum phi
  Int_t nPhiBins = 72;             // Number of phi bins
  
  // Eta
  Double_t minEta = -2.5;    // Minimum eta (current eta cut for tracks = 2.4)
  Double_t maxEta = 2.5;     // Maximum eta (current eta cut for tracks = 2.4)
  Int_t nEtaBins = 50;       // Number of eta bins
  
  // DeltaPhi in [0,pi]
  Double_t minDeltaPhi = 0;             // Minimum deltaPhi
  Double_t maxDeltaPhi = TMath::Pi();   // Maximum deltaPhi
  Int_t nDeltaPhiBins = 30;             // Number of deltaPhi bins
  
  // DeltaPhi in [-pi/2,3pi/2]
  Double_t minDeltaPhiTrack = -TMath::Pi()/2.0;    // Minimum deltaPhi for two dimensional plots
  Double_t maxDeltaPhiTrack = 3.0*TMath::Pi()/2.0; // Maximum deltaPhi for two dimensional plots
  Int_t nDeltaPhiBinsTrack = 60;                   // Number of deltaPhi bins for two dimensional plots
  
  // DeltaEta
  Double_t minDeltaEta = -5.0;   // Minimum deltaEta
  Double_t maxDeltaEta = 5.0;    // Maximum deltaEta
  Int_t nDeltaEtaBins = 100;     // Number of deltaEta bins
  
  // Dijet asymmetry
  Double_t minAsymmetry = 0;     // Minimum asymmetry
  Double_t maxAsymmetry = 0.75;  // Maximum asymmetry
  Int_t nAsymmetryBins = 25;     // Number of asymmetry bins
  
  // Vertex z-position
  Double_t minVz = -20;   // Minimum vz
  Double_t maxVz = 20;    // Maximum vz
  Int_t nVzBins = 80;     // Number of vz bins
  
  // Arrays for creating THnSparses
  const Int_t dimensionDijet = 9;
  const Int_t dimensionAnyJet = 4;
  const Int_t dimensionTrack = 9;
  const Int_t dimensionTrackPtWeight = 6;
  Int_t nBinsDijet[dimensionDijet];
  Int_t nBinsAnyjet[dimensionAnyJet];
  Int_t nBinsTrack[dimensionTrack];
  Int_t nBinsTrackPtWeight[dimensionTrackPtWeight];
  Double_t lowBinBorderDijet[dimensionDijet];
  Double_t lowBinBorderAnyJet[dimensionAnyJet];
  Double_t lowBinBorderTrack[dimensionTrack];
  Double_t lowBinBorderTrackPtWeight[dimensionTrackPtWeight];
  Double_t highBinBorderDijet[dimensionDijet];
  Double_t highBinBorderAnyJet[dimensionAnyJet];
  Double_t highBinBorderTrack[dimensionTrack];
  Double_t highBinBorderTrackPtWeight[dimensionTrackPtWeight];
  
  // ======== Plain TH1 histograms ========
  
  fhVertexZ = new TH1D("vertexZ","vertexZ",nVzBins,minVz,maxVz); fhVertexZ->Sumw2();
  fhEvents = new TH1D("nEvents","nEvents",knEventTypes,-0.5,knEventTypes-0.5); fhEvents->Sumw2();
  fhTrackCuts = new TH1D("trackCuts","trackCuts",knTrackCuts,-0.5,knTrackCuts-0.5); fhTrackCuts->Sumw2();
  fhCentrality = new TH1D("centrality","centrality",nCentralityBins,minCentrality,maxCentrality); fhCentrality->Sumw2();
  
  // For the event histogram, label each bin corresponding to an event cut
  for(Int_t i = 0; i < knEventTypes; i++){
    fhEvents->GetXaxis()->SetBinLabel(i+1,kEventTypeStrings[i]);
  }
  
  // For the track cut histogram, label each bin corresponding to a track cut
  for(Int_t i = 0; i < knTrackCuts; i++){
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


