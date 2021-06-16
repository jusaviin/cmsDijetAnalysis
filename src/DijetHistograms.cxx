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
  fhVertexZWeighted(0),
  fhVertexZDijet(0),
  fhEvents(0),
  fhTrackCuts(0),
  fhTrackCutsInclusive(0),
  fhCentrality(0),
  fhCentralityWeighted(0),
  fhCentralityDijet(0),
  fhPtHat(0),
  fhPtHatWeighted(0),
  fhLeadingJet(0),
  fhLeadingDijet(0),
  fhSubleadingDijet(0),
  fhDijet(0),
  fhAnyJet(0),
  fhTrack(0),
  fhTrackInclusive(0),
  fhTrackLeadingJet(0),
  fhTrackSubleadingJet(0),
  fhTrackUncorrected(0),
  fhTrackInclusiveUncorrected(0),
  fhTrackLeadingJetUncorrected(0),
  fhTrackSubleadingJetUncorrected(0),
  fhTrackLeadingJetPtWeighted(0),
  fhTrackSubleadingJetPtWeighted(0),
  fhTrackJetInclusive(0),
  fhTrackJetInclusivePtWeighted(0),
  fhJetPtClosure(0),
  fCard(0)
{
  // Default constructor
  for(int i = 0; i < 4; i++){
    fhQvector[i] = NULL;
    fhQvectorNorm[i] = NULL;
    fhEventPlaneMult[i] = NULL;
    fhQvectorDijet[i] = NULL;
    fhQvectorNormDijet[i] = NULL;
    fhEventPlaneMultDijet[i] = NULL;
  }
}

/*
 * Custom constructor
 */
DijetHistograms::DijetHistograms(ConfigurationCard *newCard) :
  fhVertexZ(0),
  fhVertexZWeighted(0),
  fhVertexZDijet(0),
  fhEvents(0),
  fhTrackCuts(0),
  fhTrackCutsInclusive(0),
  fhCentrality(0),
  fhCentralityWeighted(0),
  fhCentralityDijet(0),
  fhPtHat(0),
  fhPtHatWeighted(0),
  fhLeadingJet(0),
  fhLeadingDijet(0),
  fhSubleadingDijet(0),
  fhDijet(0),
  fhAnyJet(0),
  fhTrack(0),
  fhTrackInclusive(0),
  fhTrackLeadingJet(0),
  fhTrackSubleadingJet(0),
  fhTrackUncorrected(0),
  fhTrackInclusiveUncorrected(0),
  fhTrackLeadingJetUncorrected(0),
  fhTrackSubleadingJetUncorrected(0),
  fhTrackLeadingJetPtWeighted(0),
  fhTrackSubleadingJetPtWeighted(0),
  fhTrackJetInclusive(0),
  fhTrackJetInclusivePtWeighted(0),
  fhJetPtClosure(0),
  fCard(newCard)
{
  // Custom constructor
  for(int i = 0; i < 4; i++){
    fhQvector[i] = NULL;
    fhQvectorNorm[i] = NULL;
    fhEventPlaneMult[i] = NULL;
    fhQvectorDijet[i] = NULL;
    fhQvectorNormDijet[i] = NULL;
    fhEventPlaneMultDijet[i] = NULL;
  }
}

/*
 * Copy constructor
 */
DijetHistograms::DijetHistograms(const DijetHistograms& in) :
  fhVertexZ(in.fhVertexZ),
  fhVertexZWeighted(in.fhVertexZWeighted),
  fhVertexZDijet(in.fhVertexZDijet),
  fhEvents(in.fhEvents),
  fhTrackCuts(in.fhTrackCuts),
  fhTrackCutsInclusive(in.fhTrackCutsInclusive),
  fhCentrality(in.fhCentrality),
  fhCentralityWeighted(in.fhCentralityWeighted),
  fhCentralityDijet(in.fhCentralityDijet),
  fhPtHat(in.fhPtHat),
  fhPtHatWeighted(in.fhPtHatWeighted),
  fhLeadingJet(in.fhLeadingJet),
  fhLeadingDijet(in.fhLeadingDijet),
  fhSubleadingDijet(in.fhSubleadingDijet),
  fhDijet(in.fhDijet),
  fhAnyJet(in.fhAnyJet),
  fhTrack(in.fhTrack),
  fhTrackInclusive(in.fhTrackInclusive),
  fhTrackLeadingJet(in.fhTrackLeadingJet),
  fhTrackSubleadingJet(in.fhTrackSubleadingJet),
  fhTrackUncorrected(in.fhTrackUncorrected),
  fhTrackInclusiveUncorrected(in.fhTrackInclusiveUncorrected),
  fhTrackLeadingJetUncorrected(in.fhTrackLeadingJetUncorrected),
  fhTrackSubleadingJetUncorrected(in.fhTrackSubleadingJetUncorrected),
  fhTrackLeadingJetPtWeighted(in.fhTrackLeadingJetPtWeighted),
  fhTrackSubleadingJetPtWeighted(in.fhTrackSubleadingJetPtWeighted),
  fhTrackJetInclusive(in.fhTrackJetInclusive),
  fhTrackJetInclusivePtWeighted(in.fhTrackJetInclusivePtWeighted),
  fhJetPtClosure(in.fhJetPtClosure),
  fCard(in.fCard)
{
  // Copy constructor
  for(int i = 0; i < 4; i++){
    fhQvector[i] = in.fhQvector[i];
    fhQvectorNorm[i] = in.fhQvectorNorm[i];
    fhEventPlaneMult[i] = in.fhEventPlaneMult[i];
    fhQvectorDijet[i] = in.fhQvectorDijet[i];
    fhQvectorNormDijet[i] = in.fhQvectorNormDijet[i];
    fhEventPlaneMultDijet[i] = in.fhEventPlaneMultDijet[i];
  }
}

/*
 * Assingment operator
 */
DijetHistograms& DijetHistograms::operator=(const DijetHistograms& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fhVertexZ = in.fhVertexZ;
  fhVertexZWeighted = in.fhVertexZWeighted;
  fhVertexZDijet = in.fhVertexZDijet;
  fhEvents = in.fhEvents;
  fhTrackCuts = in.fhTrackCuts;
  fhTrackCutsInclusive = in.fhTrackCutsInclusive;
  fhCentrality = in.fhCentrality;
  fhCentralityWeighted = in.fhCentralityWeighted;
  fhCentralityDijet = in.fhCentralityDijet;
  fhPtHat = in.fhPtHat;
  fhPtHatWeighted = in.fhPtHatWeighted;
  fhLeadingJet = in.fhLeadingJet;
  fhLeadingDijet = in.fhLeadingDijet;
  fhSubleadingDijet = in.fhSubleadingDijet;
  fhDijet = in.fhDijet;
  fhAnyJet = in.fhAnyJet;
  fhTrack = in.fhTrack;
  fhTrackInclusive = in.fhTrackInclusive;
  fhTrackLeadingJet = in.fhTrackLeadingJet;
  fhTrackSubleadingJet = in.fhTrackSubleadingJet;
  fhTrackUncorrected = in.fhTrackUncorrected;
  fhTrackInclusiveUncorrected = in.fhTrackInclusiveUncorrected;
  fhTrackLeadingJetUncorrected = in.fhTrackLeadingJetUncorrected;
  fhTrackSubleadingJetUncorrected = in.fhTrackSubleadingJetUncorrected;
  fhTrackLeadingJetPtWeighted = in.fhTrackLeadingJetPtWeighted;
  fhTrackSubleadingJetPtWeighted = in.fhTrackSubleadingJetPtWeighted;
  fhTrackJetInclusive = in.fhTrackJetInclusive;
  fhTrackJetInclusivePtWeighted = in.fhTrackJetInclusivePtWeighted;
  fhJetPtClosure = in.fhJetPtClosure;
  fCard = in.fCard;
  
  for(int i = 0; i < 4; i++){
    fhQvector[i] = in.fhQvector[i];
    fhQvectorNorm[i] = in.fhQvectorNorm[i];
    fhEventPlaneMult[i] = in.fhEventPlaneMult[i];
    fhQvectorDijet[i] = in.fhQvectorDijet[i];
    fhQvectorNormDijet[i] = in.fhQvectorNormDijet[i];
    fhEventPlaneMultDijet[i] = in.fhEventPlaneMultDijet[i];
  }
  
  return *this;
}

/*
 * Destructor
 */
DijetHistograms::~DijetHistograms(){
  // destructor
  delete fhVertexZ;
  delete fhVertexZWeighted;
  delete fhVertexZDijet;
  delete fhEvents;
  delete fhTrackCuts;
  delete fhTrackCutsInclusive;
  delete fhCentrality;
  delete fhCentralityWeighted;
  delete fhCentralityDijet;
  delete fhPtHat;
  delete fhPtHatWeighted;
  delete fhLeadingJet;
  delete fhLeadingDijet;
  delete fhSubleadingDijet;
  delete fhDijet;
  delete fhAnyJet;
  delete fhTrack;
  delete fhTrackInclusive;
  delete fhTrackLeadingJet;
  delete fhTrackSubleadingJet;
  delete fhTrackUncorrected;
  delete fhTrackInclusiveUncorrected;
  delete fhTrackLeadingJetUncorrected;
  delete fhTrackSubleadingJetUncorrected;
  delete fhTrackLeadingJetPtWeighted;
  delete fhTrackSubleadingJetPtWeighted;
  delete fhTrackJetInclusive;
  delete fhTrackJetInclusivePtWeighted;
  delete fhJetPtClosure;
  
  for(int i = 0; i < 4; i++){
    delete fhQvector[i];
    delete fhQvectorNorm[i];
    delete fhEventPlaneMult[i];
    delete fhQvectorDijet[i];
    delete fhQvectorNormDijet[i];
    delete fhEventPlaneMultDijet[i];
  }

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
  const Double_t minCentrality = -0.75;   // Minimum centrality bin, is negative since hiBin is -1 for pp
  const Double_t maxCentrality = 100.25;  // Maximum centrality bin
  const Int_t nCentralityBins = 202;      // Number of centrality bins
  
  // Jet pT
  const Double_t minPtJet = 0;     // Minimum jet pT
  const Double_t maxPtJet = 500;   // Maximum jet pT
  const Int_t nPtBinsJet = 100;    // Number of jet pT bins
  
  //Track pT
  const Double_t minPtTrack = 0;   // Minimum track pT for track histograms
  const Double_t maxPtTrack = 20;  // Maximum track pT for track histograms (Hallie's analysis = 20)
  const Int_t nPtBinsTrack = 500;  // Number of track pT bins for track histograms (Hallie's analysis = 500)
  
  // Phi
  const Double_t minPhi = -TMath::Pi();  // Minimum phi
  const Double_t maxPhi = TMath::Pi();   // Maximum phi
  const Int_t nPhiBins = 64;             // Number of phi bins
  
  // Eta
  const Double_t minEta = -2.5;    // Minimum eta (current eta cut for tracks = 2.4)
  const Double_t maxEta = 2.5;     // Maximum eta (current eta cut for tracks = 2.4)
  const Int_t nEtaBins = 50;       // Number of eta bins
  
  // DeltaPhi in [0,pi]
  const Double_t minDeltaPhi = 0;             // Minimum deltaPhi
  const Double_t maxDeltaPhi = TMath::Pi();   // Maximum deltaPhi
  const Int_t nDeltaPhiBins = 32;             // Number of deltaPhi bins
  
  // DeltaPhi in [-pi/2,3pi/2]
  const Double_t minDeltaPhiJetTrack = -TMath::Pi()/2.0;    // Minimum deltaPhi for jet-track correlations
  const Double_t maxDeltaPhiJetTrack = 3.0*TMath::Pi()/2.0; // Maximum deltaPhi for jet-track correlations
  const Int_t nDeltaPhiBinsJetTrack = 200;                  // Number of deltaPhi bins for jet-track correlations (match the common number in UIC group)
  
  // DeltaEta
  const Double_t minDeltaEtaJetTrack = -5.0;   // Minimum deltaEta for jet-track correlations
  const Double_t maxDeltaEtaJetTrack = 5.0;    // Maximum deltaEta for jet-track correlations
  const Int_t nDeltaEtaBinsJetTrack = 500;     // Number of deltaEta bins for jet-track correlations (match the common number in UIC group)
  
  // Dijet asymmetry AJ
  const Double_t minAsymmetry = 0;     // Minimum asymmetry
  const Double_t maxAsymmetry = 0.75;  // Maximum asymmetry
  const Int_t nAsymmetryBins = 25;     // Number of asymmetry bins
  
  // Dijet asymmetry xJ
  const Double_t minXj = 0;  // Minimum xJ asymmetry
  const Double_t maxXj = 1;  // Maximum xJ asymmetry
  const Int_t nXjBins = 40;
  
  // Vertex z-position
  const Double_t minVz = -20;   // Minimum vz
  const Double_t maxVz = 20;    // Maximum vz
  const Int_t nVzBins = 80;     // Number of vz bins
  
  // pT hat
  const Double_t minPtHat = 0;     // Minimum pT hat
  const Double_t maxPtHat = 460;   // Maximum pT hat
  const Int_t nFinePtHatBins = 230; // Number of fine pT hat bins
  
  // Correlation types
  const Double_t minCorrelationType = -0.5;                   // Correlation type indexing starts from zero
  const Double_t maxCorrelationType = knCorrelationTypes-0.5; // Maximum correlation type index
  const Int_t nCorrelationTypeBins = knCorrelationTypes;      // Make a bin with width of 1 for each correlation type
  
  // Closure type (0 = leading, 1 = subleading, 2 = inclusive)
  const Double_t minClosureType = -0.5;                // Closure type indexing starts from zero
  const Double_t maxClosureType = knClosureTypes-0.5;  // Maximum closure type index
  const Int_t nClosureTypeBins = knClosureTypes;       // Make a bin width of 1 for each closure type
  
  // Generator level pT binning for closure histograms
  const Double_t minClosurePt = 50;                             // Minimum gen jet pT for closure plots
  const Double_t maxClosurePt = 500;                            // Maximum gen jet pT for closure plots
  const Int_t nClosurePtBins = (maxClosurePt-minClosurePt)/10;  // Bin width of 10 for the Gen pT in closure plots
  
  // Particle type for closure plots (0 = quark, 1 = gluon)
  const Double_t minClosureParticleType = -0.5;                        // Closure particle type indexing starts from zero
  const Double_t maxClosureParticleType = knClosureParticleTypes-0.5;  // Maximum closure particle type index
  const Int_t nClosureParticleTypeBins = knClosureParticleTypes;       // Bin width for particle type is 1
  
  // Binning for reco/gen ratio for closure histograms
  const Double_t minClosureRatio = 0;    // Minimum ratio for the closure plots
  const Double_t maxClosureRatio = 2;    // Maximum ratio for the closure plots
  const Int_t nClosureRatioBins = 40;    // Number of closure ratio bins
  
  // Magnitude for event plane q-vector
  const Double_t minQVector = 0;    // Minimum value for the magnitude of q-vector
  const Double_t maxQVector = 5;  // Maximum value for the magnitude of q-vector
  const Int_t nQVectorBins = 4;   // Number of q-vector magnitude bins
  
  Bool_t doEventPlane = (fCard->Get("IncludeEventPlane") == 1);
  
  // Centrality bins for THnSparses (We run into memory issues, if have all the bins)
  const Int_t nWideCentralityBins = fCard->GetNBin("CentralityBinEdges");
  Double_t wideCentralityBins[nWideCentralityBins+1];
  for(Int_t iCentrality = 0; iCentrality < nWideCentralityBins+1; iCentrality++){
    wideCentralityBins[iCentrality] = fCard->Get("CentralityBinEdges",iCentrality);
  }
  
  // Track pT bins for THnSparses (We run into memory issues, if have all the bins)
  const Int_t nWideTrackPtBins = fCard->GetNBin("TrackPtBinEdges");
  Double_t wideTrackPtBins[nWideTrackPtBins+1];
  for(Int_t iTrackPt = 0; iTrackPt < nWideTrackPtBins+1; iTrackPt++){
    wideTrackPtBins[iTrackPt] = fCard->Get("TrackPtBinEdges",iTrackPt);
  }
  
  // Dijet asymmetry bins for THnSparses (We run into memory issues, if have all the bins)
  const Int_t nWideAsymmetryBins = fCard->GetNBin("AsymmetryBinEdges");
  Double_t wideAsymmetryBins[nWideAsymmetryBins+1];
  for(Int_t iAsymmetryBin = 0; iAsymmetryBin < nWideAsymmetryBins+1; iAsymmetryBin++){
    wideAsymmetryBins[iAsymmetryBin] = fCard->Get("AsymmetryBinEdges",iAsymmetryBin);
  }
  
  // Bins for the pT hat histogram
  const Int_t nPtHatBins = fCard->GetNBin("PtHatBinEdges");
  Double_t ptHatBins[nPtHatBins+1];
  for(Int_t iPtHat = 0; iPtHat < nPtHatBins+1; iPtHat++){
    ptHatBins[iPtHat] = fCard->Get("PtHatBinEdges",iPtHat);
  }
  
  // Binning for Qvector
  const Int_t nWideQvectorBins = fCard->GetNBin("QvectorBinEdges");
  Double_t wideQvectorBins[nWideQvectorBins+1];
  for(Int_t iQ = 0; iQ < nWideQvectorBins+1; iQ++){
    wideQvectorBins[iQ] = fCard->Get("QvectorBinEdges",iQ);
  }
  
  // Arrays for creating THnSparses
  Int_t nBins5D[5];
  Int_t nBins6D[6];
  Int_t nBins7D[7];
  Double_t lowBinBorder5D[5];
  Double_t lowBinBorder6D[6];
  Double_t lowBinBorder7D[7];
  Double_t highBinBorder5D[5];
  Double_t highBinBorder6D[6];
  Double_t highBinBorder7D[7];
  
  // ======== Plain TH1 histograms ========
  
  fhVertexZ = new TH1F("vertexZ","vertexZ",nVzBins,minVz,maxVz); fhVertexZ->Sumw2();
  fhVertexZWeighted = new TH1F("vertexZweighted","vertexZweighted",nVzBins,minVz,maxVz); fhVertexZWeighted->Sumw2();
  fhVertexZDijet = new TH1F("vertexZdijet","vertexZdijet",nVzBins,minVz,maxVz); fhVertexZDijet->Sumw2();
  fhEvents = new TH1F("nEvents","nEvents",knEventTypes,-0.5,knEventTypes-0.5); fhEvents->Sumw2();
  fhTrackCuts = new TH1F("trackCuts","trackCuts",knTrackCuts,-0.5,knTrackCuts-0.5); fhTrackCuts->Sumw2();
  fhTrackCutsInclusive = new TH1F("trackCutsInclusive","trackCutsInclusive",knTrackCuts,-0.5,knTrackCuts-0.5); fhTrackCutsInclusive->Sumw2();
  fhCentrality = new TH1F("centrality","centrality",nCentralityBins,minCentrality,maxCentrality); fhCentrality->Sumw2();
  fhCentralityWeighted = new TH1F("centralityWeighted","centralityWeighted",nCentralityBins,minCentrality,maxCentrality); fhCentralityWeighted->Sumw2();
  fhCentralityDijet = new TH1F("centralityDijet","centralityDijet",nCentralityBins,minCentrality,maxCentrality); fhCentralityDijet->Sumw2();
  fhPtHat = new TH1F("pthat","pthat",nPtHatBins,ptHatBins); fhPtHat->Sumw2();
  fhPtHatWeighted = new TH1F("pthatWeighted","pthatWeighted",nFinePtHatBins,minPtHat,maxPtHat); fhPtHatWeighted->Sumw2();
  
  for(int i = 0; i < 4; i++){
    fhQvector[i] = new TH1F(Form("qVector%d",i),Form("qVector%d",i),100,0,200); fhQvector[i]->Sumw2();
    fhQvectorNorm[i] = new TH1F(Form("qVectorNorm%d",i),Form("qVectorNorm%d",i),120,0,6); fhQvectorNorm[i]->Sumw2();
    fhEventPlaneMult[i] = new TH1F(Form("eventPlaneMult%d",i),Form("eventPlaneMult%d",i),100,0,1800); fhEventPlaneMult[i]->Sumw2();
    fhQvectorDijet[i] = new TH1F(Form("qVectorDijet%d",i),Form("qVectorDijet%d",i),100,0,200); fhQvectorDijet[i]->Sumw2();
    fhQvectorNormDijet[i] = new TH1F(Form("qVectorNormDijet%d",i),Form("qVectorNormDijet%d",i),120,0,6); fhQvectorNormDijet[i]->Sumw2();
    fhEventPlaneMultDijet[i] = new TH1F(Form("eventPlaneMultDijet%d",i),Form("eventPlaneMultDijet%d",i),100,0,1800); fhEventPlaneMultDijet[i]->Sumw2();
  }
  
  
  // For the event histogram, label each bin corresponding to an event cut
  for(Int_t i = 0; i < knEventTypes; i++){
    fhEvents->GetXaxis()->SetBinLabel(i+1,kEventTypeStrings[i]);
  }
  
  // If we are using PF jets, change the axis label for that
  if(fCard->Get("JetType") == 1) fhEvents->GetXaxis()->SetBinLabel(kCaloJet+1,"PFJet");
  
  // For the track cut histogram, label each bin corresponding to a track cut
  for(Int_t i = 0; i < knTrackCuts; i++){
    fhTrackCuts->GetXaxis()->SetBinLabel(i+1,kTrackCutStrings[i]);
    fhTrackCutsInclusive->GetXaxis()->SetBinLabel(i+1,kTrackCutStrings[i]);
  }
  
  // ======== THnSparses for leading and subleading jets ========
  
  // Axis 0 for the jet histogram: leading/subleading jet pT
  nBins6D[0] = nPtBinsJet;         // nBins for leading/subleading jet pT
  lowBinBorder6D[0] = minPtJet;    // low bin border for leading/subleading jet pT
  highBinBorder6D[0] = maxPtJet;   // high bin border for leading/subleading jet pT
  
  // Axis 1 for the jet histogram: leading/subleading jet phi
  nBins6D[1] = nPhiBins;        // nBins for leading/subleading jet phi
  lowBinBorder6D[1] = minPhi;   // low bin border for leading/subleading jet phi
  highBinBorder6D[1] = maxPhi;  // high bin border for leading/subleading jet phi
  
  // Axis 2 for the jet histogram: leading/subleading jet eta
  nBins6D[2] = nEtaBins;        // nBins for leading/subleading jet eta
  lowBinBorder6D[2] = minEta;   // low bin border for leading/subleading jet eta
  highBinBorder6D[2] = maxEta;  // high bin border for leading/subleading jet eta
  
  // Axis 3 for the jet histogram: asymmetry
  nBins6D[3] = nWideAsymmetryBins;     // nBins for wide asymmetry bins
  lowBinBorder6D[3] = minAsymmetry;    // low bin border for asymmetry
  highBinBorder6D[3] = maxAsymmetry;   // high bin border for asymmetry
  
  // Axis 4 for the jet histogram: centrality
  nBins6D[4] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorder6D[4] = minCentrality;  // low bin border for centrality
  highBinBorder6D[4] = maxCentrality; // high bin border for centrality
  
  if(doEventPlane){ //
    
    // Axis 5 for the jet histogram: magnitude of event plane q-vector
    nBins6D[5] = nWideQvectorBins;    // nBins for event plane q-vector magnitude
    lowBinBorder6D[5] = minQVector;   // low bin border for event plane q-vector magnitude
    highBinBorder6D[5] = maxQVector;  // high bin border for event plane q-vector magnitude
    
  } else {
  
    // Axis 5 for the jet histogram: jet flavor (quark/gluon)
    nBins6D[5] = nClosureParticleTypeBins;        // nBins for jet flavor
    lowBinBorder6D[5] = minClosureParticleType;   // low bin border for jet flavor
    highBinBorder6D[5] = maxClosureParticleType;  // high bin border for jet flavor
    
  }
  
  // Create the histograms for leading and subleading jets using the above binning information
  fhLeadingDijet = new THnSparseF("leadingJet","leadingJet",6,nBins6D,lowBinBorder6D,highBinBorder6D); fhLeadingDijet->Sumw2();
  fhSubleadingDijet = new THnSparseF("subleadingJet","subleadingJet",6,nBins6D,lowBinBorder6D,highBinBorder6D); fhSubleadingDijet->Sumw2();
  
  // Set custom dijet asymmetry bins for histograms
  fhLeadingDijet->SetBinEdges(3,wideAsymmetryBins);
  fhSubleadingDijet->SetBinEdges(3,wideAsymmetryBins);
  
  // Set custom centrality bins for histograms
  fhLeadingDijet->SetBinEdges(4,wideCentralityBins);
  fhSubleadingDijet->SetBinEdges(4,wideCentralityBins);
  
  if(doEventPlane){ //
    fhLeadingDijet->SetBinEdges(5,wideQvectorBins);
    fhSubleadingDijet->SetBinEdges(5,wideQvectorBins);
  }

  // ========= THnSparse for dijets =========
  
  // Axis 0 for the dijet histogram: leading jet pT
  nBins6D[0] = nPtBinsJet;         // nBins for leading jet pT
  lowBinBorder6D[0] = minPtJet;    // low bin border for leading jet pT
  highBinBorder6D[0] = maxPtJet;   // high bin border for leading jet pT
  
  // Axis 1 for the dijet histogram: subleading jet pT
  nBins6D[1] = nPtBinsJet;         // nBins for subleading jet pT
  lowBinBorder6D[1] = minPtJet;    // low bin border for subleading jet pT
  highBinBorder6D[1] = maxPtJet;   // high bin border for subleading jet pT
  
  // Axis 2 for the dijet histogram: deltaPhi
  nBins6D[2] = nDeltaPhiBins;       // nBins for deltaPhi
  lowBinBorder6D[2] = minDeltaPhi;  // low bin border for deltaPhi
  highBinBorder6D[2] = maxDeltaPhi; // high bin border for deltaPhi
  
  // Axis 3 for the dijet histogram: asymmetry AJ
  nBins6D[3] = nAsymmetryBins;         // nBins for asymmetry
  lowBinBorder6D[3] = minAsymmetry;    // low bin border for asymmetry
  highBinBorder6D[3] = maxAsymmetry;   // high bin border for asymmetry
  
  // Axis 4 for the dijet histogram: centrality
  nBins6D[4] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorder6D[4] = minCentrality;  // low bin border for centrality
  highBinBorder6D[4] = maxCentrality; // high bin border for centrality
  
  // Axis 5 for the dijet histogram: asymmetry xJ
  nBins6D[5] = nXjBins;         // nBins for xJ asymmetry
  lowBinBorder6D[5] = minXj;    // low bin border for xJ asymmetry
  highBinBorder6D[5] = maxXj;   // high bin border for xJ asymmetry
  
  // Create the dijet histogram using the above binning information
  fhDijet = new THnSparseF("dijet","dijet",6,nBins6D,lowBinBorder6D,highBinBorder6D); fhDijet->Sumw2();
  
  // Set custom centrality bins for histograms
  fhDijet->SetBinEdges(4,wideCentralityBins);
  
  // ======== THnSparse for all jets ========
  
  // Axis 0 for the any jet histogram: jet pT
  nBins5D[0] = nPtBinsJet;         // nBins for any jet pT
  lowBinBorder5D[0] = minPtJet;    // low bin border for any jet pT
  highBinBorder5D[0] = maxPtJet;   // high bin border for any jet pT
  
  // Axis 1 for the any jet histogram: jet phi
  nBins5D[1] = nPhiBins;        // nBins for any jet phi
  lowBinBorder5D[1] = minPhi;   // low bin border for any jet phi
  highBinBorder5D[1] = maxPhi;  // high bin border for any jet phi
  
  // Axis 2 for the any jet histogram: jet eta
  nBins5D[2] = nEtaBins;        // nBins for any jet eta
  lowBinBorder5D[2] = minEta;   // low bin border for any jet eta
  highBinBorder5D[2] = maxEta;  // high bin border for any jet eta
  
  // Axis 3 for the any jet histogram: centrality
  nBins5D[3] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorder5D[3] = minCentrality;  // low bin border for centrality
  highBinBorder5D[3] = maxCentrality; // high bin border for centrality
  
  // Axis 4 for the jet histogram: jet flavor (quark/gluon)
  nBins5D[4] = nClosureParticleTypeBins;        // nBins for jet flavor
  lowBinBorder5D[4] = minClosureParticleType;   // low bin border for jet flavor
  highBinBorder5D[4] = maxClosureParticleType;  // high bin border for jet flavor
  
  // Create the histogram for all jets using the above binning information
  fhAnyJet = new THnSparseF("anyJet","anyJet",5,nBins5D,lowBinBorder5D,highBinBorder5D); fhAnyJet->Sumw2();
  fhLeadingJet = new THnSparseF("anyLeadingJet","anyLeadingJet",5,nBins5D,lowBinBorder5D,highBinBorder5D); fhLeadingJet->Sumw2();

  // Set custom centrality bins for histograms
  fhAnyJet->SetBinEdges(3,wideCentralityBins);
  fhLeadingJet->SetBinEdges(3,wideCentralityBins);
  
  // ======== THnSparses for tracks and uncorrected tracks ========
  
  // Axis 0 for the track histogram: track pT
  nBins5D[0] = nPtBinsTrack;         // nBins for track pT
  lowBinBorder5D[0] = minPtTrack;    // low bin border for track pT
  highBinBorder5D[0] = maxPtTrack;   // high bin border for track pT
  
  // Axis 1 for the track histogram: track phi
  nBins5D[1] = nPhiBins;         // nBins for track phi
  lowBinBorder5D[1] = minPhi;    // low bin border for track phi
  highBinBorder5D[1] = maxPhi;   // high bin border for track phi
  
  // Axis 2 for the track histogram: track eta
  nBins5D[2] = nEtaBins;         // nBins for track eta
  lowBinBorder5D[2] = minEta;    // low bin border for track eta
  highBinBorder5D[2] = maxEta;   // high bin border for track eta
  
  // Axis 3 for the track histogram: centrality
  nBins5D[3] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorder5D[3] = minCentrality;  // low bin border for centrality
  highBinBorder5D[3] = maxCentrality; // high bin border for centrality
  
  // Axis 4 for the track histogram: correlation type
  nBins5D[4] = nCorrelationTypeBins;        // nBins for correlation types
  lowBinBorder5D[4] = minCorrelationType;   // low bin border for correlation types
  highBinBorder5D[4] = maxCorrelationType;  // high bin border for correlation types
  
  // Create the histograms for tracks and uncorrected tracks using the above binning information
  fhTrack = new THnSparseF("track","track",5,nBins5D,lowBinBorder5D,highBinBorder5D); fhTrack->Sumw2();
  fhTrackInclusive = new THnSparseF("trackInclusive","trackInclusive",4,nBins5D,lowBinBorder5D,highBinBorder5D); fhTrackInclusive->Sumw2();
  fhTrackUncorrected = new THnSparseF("trackUncorrected","trackUncorrected",5,nBins5D,lowBinBorder5D,highBinBorder5D); fhTrackUncorrected->Sumw2();
  fhTrackInclusiveUncorrected = new THnSparseF("trackInclusiveUncorrected","trackInclusiveUncorrected",4,nBins5D,lowBinBorder5D,highBinBorder5D); fhTrackInclusiveUncorrected->Sumw2();

  // Set custom centrality bins for histograms
  fhTrack->SetBinEdges(3,wideCentralityBins);
  fhTrackInclusive->SetBinEdges(3,wideCentralityBins);
  fhTrackUncorrected->SetBinEdges(3,wideCentralityBins);
  fhTrackInclusiveUncorrected->SetBinEdges(3,wideCentralityBins);
  
  // ======== THnSparses for correlation between tracks and leading or subleading jets ========
  
  // Axis 0 for the track-jet correlation histogram: track pT
  nBins7D[0] = nWideTrackPtBins;     // nBins for wide track pT bins
  lowBinBorder7D[0] = minPtTrack;    // low bin border for track pT
  highBinBorder7D[0] = maxPtTrack;   // high bin border for track pT
  
  // Axis 1 for the track-jet correlation histogram: deltaPhi between track and jet
  nBins7D[1] = nDeltaPhiBinsJetTrack;       // nBins for deltaPhi between track and jet
  lowBinBorder7D[1] = minDeltaPhiJetTrack;  // low bin border for deltaPhi between track and jet
  highBinBorder7D[1] = maxDeltaPhiJetTrack; // high bin border for deltaPhi between track and jet
  
  // Axis 2 for the track-jet correlation histogram: deltaEta between track and jet
  nBins7D[2] = nDeltaEtaBinsJetTrack;         // nBins for deltaEta between track and jet
  lowBinBorder7D[2] = minDeltaEtaJetTrack;    // low bin border deltaEta between track and jet
  highBinBorder7D[2] = maxDeltaEtaJetTrack;   // high bin border deltaEta between track and jet
  
  // Axis 3 for the track-jet correlation histogram: dijet asymmetry
  nBins7D[3] = nWideAsymmetryBins;     // nBins for wide dijet asymmetry bins
  lowBinBorder7D[3] = minAsymmetry;    // low bin border for dijet asymmetry
  highBinBorder7D[3] = maxAsymmetry;   // high bin border for dijet asymmetry
  
  // Axis 4 for the track-jet correlation histogram: centrality
  nBins7D[4] = nWideCentralityBins;     // nBins for centrality
  lowBinBorder7D[4] = minCentrality;    // low bin border for centrality
  highBinBorder7D[4] = maxCentrality;   // high bin border for centrality
  
  // Axis 5 for the track-jet correlation histogram: correlation type
  nBins7D[5] = nCorrelationTypeBins;        // nBins for correlation types
  lowBinBorder7D[5] = minCorrelationType;   // low bin border for correlation types
  highBinBorder7D[5] = maxCorrelationType;  // high bin border for correlation types
  
  if(doEventPlane){ //
    
    // Axis 6 for the track-jet correlation histogram: magnitude of event plane q-vector
    nBins7D[6] = nWideQvectorBins;    // nBins for event plane q-vector magnitude
    lowBinBorder7D[6] = minQVector;   // low bin border for event plane q-vector magnitude
    highBinBorder7D[6] = maxQVector;  // high bin border for event plane q-vector magnitude
    
  } else {
  
    // Axis 6 for the track-jet correlation histogram: jet flavor (quark/gluon)
    nBins7D[6] = nClosureParticleTypeBins;        // nBins for jet flavors
    lowBinBorder7D[6] = minClosureParticleType;   // low bin border for jet flavor
    highBinBorder7D[6] = maxClosureParticleType;  // high bin border for jet flavor
    
  }
  
  // Create histograms for track-jet correlations
  fhTrackLeadingJet = new THnSparseF("trackLeadingJet","trackLeadingJet",7,nBins7D,lowBinBorder7D,highBinBorder7D); fhTrackLeadingJet->Sumw2();
  fhTrackSubleadingJet = new THnSparseF("trackSubleadingJet","trackSubleadingJet",7,nBins7D,lowBinBorder7D,highBinBorder7D); fhTrackSubleadingJet->Sumw2();
  
  // Set custom centrality, asymmetry and track pT bins for histograms
  fhTrackLeadingJet->SetBinEdges(0,wideTrackPtBins);
  fhTrackSubleadingJet->SetBinEdges(0,wideTrackPtBins);
  fhTrackLeadingJet->SetBinEdges(3,wideAsymmetryBins);
  fhTrackSubleadingJet->SetBinEdges(3,wideAsymmetryBins);
  fhTrackLeadingJet->SetBinEdges(4,wideCentralityBins);
  fhTrackSubleadingJet->SetBinEdges(4,wideCentralityBins);
  if(doEventPlane){ //
    fhTrackLeadingJet->SetBinEdges(6,wideQvectorBins);
    fhTrackSubleadingJet->SetBinEdges(6,wideQvectorBins);
  }

  // Create histograms for uncorrected track-jet correlations
  fhTrackLeadingJetUncorrected = new THnSparseF("trackLeadingJetUncorrected","trackLeadingJetUncorrected",7,nBins7D,lowBinBorder7D,highBinBorder7D); fhTrackLeadingJetUncorrected->Sumw2();
  fhTrackSubleadingJetUncorrected = new THnSparseF("trackSubleadingJetUncorrected","trackSubleadingJetUncorrected",7,nBins7D,lowBinBorder7D,highBinBorder7D); fhTrackSubleadingJetUncorrected->Sumw2();
  
  // Set custom centrality, asymmetry and track pT bins for histograms
  fhTrackLeadingJetUncorrected->SetBinEdges(0,wideTrackPtBins);
  fhTrackSubleadingJetUncorrected->SetBinEdges(0,wideTrackPtBins);
  fhTrackLeadingJetUncorrected->SetBinEdges(3,wideAsymmetryBins);
  fhTrackSubleadingJetUncorrected->SetBinEdges(3,wideAsymmetryBins);
  fhTrackLeadingJetUncorrected->SetBinEdges(4,wideCentralityBins);
  fhTrackSubleadingJetUncorrected->SetBinEdges(4,wideCentralityBins);
  if(doEventPlane){ //
    fhTrackLeadingJetUncorrected->SetBinEdges(6,wideQvectorBins);
    fhTrackSubleadingJetUncorrected->SetBinEdges(6,wideQvectorBins);
  }
  
  // Create histograms for pT weighted track-jet correlations
  fhTrackLeadingJetPtWeighted = new THnSparseF("trackLeadingJetPtWeighted","trackLeadingJetPtWeighted",7,nBins7D,lowBinBorder7D,highBinBorder7D); fhTrackLeadingJetPtWeighted->Sumw2();
  fhTrackSubleadingJetPtWeighted = new THnSparseF("trackSubleadingJetPtWeighted","trackSubleadingJetPtWeighted",7,nBins7D,lowBinBorder7D,highBinBorder7D); fhTrackSubleadingJetPtWeighted->Sumw2();
  
  // Set custom centrality, asymmetry and track pT bins for histograms
  fhTrackLeadingJetPtWeighted->SetBinEdges(0,wideTrackPtBins);
  fhTrackSubleadingJetPtWeighted->SetBinEdges(0,wideTrackPtBins);
  fhTrackLeadingJetPtWeighted->SetBinEdges(3,wideAsymmetryBins);
  fhTrackSubleadingJetPtWeighted->SetBinEdges(3,wideAsymmetryBins);
  fhTrackLeadingJetPtWeighted->SetBinEdges(4,wideCentralityBins);
  fhTrackSubleadingJetPtWeighted->SetBinEdges(4,wideCentralityBins);
  if(doEventPlane){ //
    fhTrackLeadingJetPtWeighted->SetBinEdges(6,wideQvectorBins);
    fhTrackSubleadingJetPtWeighted->SetBinEdges(6,wideQvectorBins);
  }
  
  // ======== THnSparses for correlation between tracks and inclusive jets ========
  
  // Axis 0 for the track-inclusive jet correlation histogram: track pT
  nBins6D[0] = nWideTrackPtBins;     // nBins for wide track pT bins
  lowBinBorder6D[0] = minPtTrack;    // low bin border for track pT
  highBinBorder6D[0] = maxPtTrack;   // high bin border for track pT
  
  // Axis 1 for the track-inclusive jet correlation histogram: deltaPhi between track and jet
  nBins6D[1] = nDeltaPhiBinsJetTrack;       // nBins for deltaPhi between track and jet
  lowBinBorder6D[1] = minDeltaPhiJetTrack;  // low bin border for deltaPhi between track and jet
  highBinBorder6D[1] = maxDeltaPhiJetTrack; // high bin border for deltaPhi between track and jet
  
  // Axis 2 for the track-inclusive jet correlation histogram: deltaEta between track and jet
  nBins6D[2] = nDeltaEtaBinsJetTrack;         // nBins for deltaEta between track and jet
  lowBinBorder6D[2] = minDeltaEtaJetTrack;    // low bin border deltaEta between track and jet
  highBinBorder6D[2] = maxDeltaEtaJetTrack;   // high bin border deltaEta between track and jet
  
  // Axis 3 for the track-inclusive jet correlation histogram: centrality
  nBins6D[3] = nWideCentralityBins;     // nBins for centrality
  lowBinBorder6D[3] = minCentrality;    // low bin border for centrality
  highBinBorder6D[3] = maxCentrality;   // high bin border for centrality
  
  // Axis 4 for the track-inclusive jet correlation histogram: correlation type
  nBins6D[4] = nCorrelationTypeBins;        // nBins for correlation types
  lowBinBorder6D[4] = minCorrelationType;   // low bin border for correlation types
  highBinBorder6D[4] = maxCorrelationType;  // high bin border for correlation types
  
  if(doEventPlane){ //
    
    // Axis 5 for the track-inclusive jet correlation histogram: magnitude of event plane q-vector
    nBins6D[5] = nWideQvectorBins;    // nBins for event plane q-vector magnitude
    lowBinBorder6D[5] = minQVector;   // low bin border for event plane q-vector magnitude
    highBinBorder6D[5] = maxQVector;  // high bin border for event plane q-vector magnitude
    
  } else {
  
    // Axis 5 for the track-inclusive jet correlation histogram: jet flavor (quark/gluon)
    nBins6D[5] = nClosureParticleTypeBins;        // nBins for jet flavors
    lowBinBorder6D[5] = minClosureParticleType;   // low bin border for jet flavor
    highBinBorder6D[5] = maxClosureParticleType;  // high bin border for jet flavor
    
  }
  
  // Create histograms for track-inclusive jet correlations
  fhTrackJetInclusive = new THnSparseF("trackJetInclusive","trackJetInclusive",6,nBins6D,lowBinBorder6D,highBinBorder6D); fhTrackJetInclusive->Sumw2();
  fhTrackJetInclusivePtWeighted = new THnSparseF("trackJetInclusivePtWeighted","trackJetInclusivePtWeighted",6,nBins6D,lowBinBorder6D,highBinBorder6D); fhTrackJetInclusivePtWeighted->Sumw2();
  
  // Set custom centrality and track pT bins for histograms
  fhTrackJetInclusive->SetBinEdges(0,wideTrackPtBins);
  fhTrackJetInclusivePtWeighted->SetBinEdges(0,wideTrackPtBins);
  fhTrackJetInclusive->SetBinEdges(3,wideCentralityBins);
  fhTrackJetInclusivePtWeighted->SetBinEdges(3,wideCentralityBins);
  if(doEventPlane){ //
    fhTrackJetInclusive->SetBinEdges(5,wideQvectorBins);
    fhTrackJetInclusivePtWeighted->SetBinEdges(5,wideQvectorBins);
  }
  
  // ======== THnSparses for jet pT closures ========
  
  // Axis 0 for the jet pT closure histogram: closure type (leading/subleading/inclusive)
  nBins7D[0] = nClosureTypeBins;         // nBins for closure types
  lowBinBorder7D[0] = minClosureType;    // low bin border for closure types
  highBinBorder7D[0] = maxClosureType;   // high bin border for closure types
  
  // Axis 1 for the jet pT closure histogram: generator level jet pT
  nBins7D[1] = nClosurePtBins;       // nBins for generator level pT bins in closure plots
  lowBinBorder7D[1] = minClosurePt;  // low bin border generator level pT in closure plots
  highBinBorder7D[1] = maxClosurePt; // high bin border generator level pT in closure plots
  
  // Axis 2 for the jet pT closure histogram: reconstructed jet pT
  nBins7D[2] = nClosurePtBins;       // nBins for reconstructed jet pT bins in closure plots
  lowBinBorder7D[2] = minClosurePt;  // low bin border for reconstructed jet pT in closure plots
  highBinBorder7D[2] = maxClosurePt; // high bin border for reconstructed jet pT in closure plots
  
  // Axis 3 for the jet pT closure histogram: generator level jet eta
  nBins7D[3] = nEtaBins;             // nBins for jet eta
  lowBinBorder7D[3] = minEta;        // low bin border for jet eta
  highBinBorder7D[3] = maxEta;       // high bin border for jet eta
  
  // Axis 4 for the jet pT closure histogram: centrality
  nBins7D[4] = nWideCentralityBins;     // nBins for centrality
  lowBinBorder7D[4] = minCentrality;    // low bin border for centrality
  highBinBorder7D[4] = maxCentrality;   // high bin border for centrality
  
  // Axis 5 for the jet pT closure histogram: ref parton = quark/gluon
  nBins7D[5] = nClosureParticleTypeBins;         // nBins for reference parton
  lowBinBorder7D[5] = minClosureParticleType;    // low bin border for reference parton
  highBinBorder7D[5] = maxClosureParticleType;   // high bin border for reference parton
  
  // Axis 6 for the jet pT closure histogram: reco/gen ratio for closure
  nBins7D[6] = nClosureRatioBins;        // nBins for closure ratio
  lowBinBorder7D[6] = minClosureRatio;   // low bin border for closure ratio
  highBinBorder7D[6] = maxClosureRatio;  // high bin border for closure ratio
  
  // Create histograms for jet pT closure
  fhJetPtClosure = new THnSparseF("jetPtClosure","jetPtClosure",7,nBins7D,lowBinBorder7D,highBinBorder7D); fhJetPtClosure->Sumw2();
  
  // Set custom centrality bins for histograms
  fhJetPtClosure->SetBinEdges(4,wideCentralityBins);
}

/*
 * Write the histograms to file
 */
void DijetHistograms::Write() const{
  
  // Write the histograms to file
  fhVertexZ->Write();
  fhVertexZWeighted->Write();
  fhVertexZDijet->Write();
  fhEvents->Write();
  fhTrackCuts->Write();
  fhTrackCutsInclusive->Write();
  fhCentrality->Write();
  fhCentralityWeighted->Write();
  fhCentralityDijet->Write();
  fhPtHat->Write();
  fhPtHatWeighted->Write();
  fhLeadingJet->Write();
  fhLeadingDijet->Write();
  fhSubleadingDijet->Write();
  fhDijet->Write();
  fhAnyJet->Write();
  fhTrack->Write();
  fhTrackInclusive->Write();
  fhTrackLeadingJet->Write();
  fhTrackSubleadingJet->Write();
  fhTrackUncorrected->Write();
  fhTrackInclusiveUncorrected->Write();
  fhTrackLeadingJetUncorrected->Write();
  fhTrackSubleadingJetUncorrected->Write();
  fhTrackLeadingJetPtWeighted->Write();
  fhTrackSubleadingJetPtWeighted->Write();
  fhTrackJetInclusive->Write();
  fhTrackJetInclusivePtWeighted->Write();
  fhJetPtClosure->Write();
  
  for(int i = 0; i < 4; i++){
    fhQvector[i]->Write();
    fhQvectorNorm[i]->Write();
    fhEventPlaneMult[i]->Write();
    fhQvectorDijet[i]->Write();
    fhQvectorNormDijet[i]->Write();
    fhEventPlaneMultDijet[i]->Write();
  }
  
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
  
  // Delete the outputFile object
  delete outputFile;
}


