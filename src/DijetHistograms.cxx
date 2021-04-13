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
  fhMultiplicity(0),
  fhMultiplicityDijet(0),
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
  fhJetEventPlaneForwardRap(0),
  fhJetEventPlaneMidRap(0),
  fCard(0)
{
  // Default constructor
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
  fhMultiplicity(0),
  fhMultiplicityDijet(0),
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
  fhJetEventPlaneForwardRap(0),
  fhJetEventPlaneMidRap(0),
  fCard(newCard)
{
  // Custom constructor
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
  fhMultiplicity(in.fhMultiplicity),
  fhMultiplicityDijet(in.fhMultiplicityDijet),
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
  fhJetEventPlaneForwardRap(in.fhJetEventPlaneForwardRap),
  fhJetEventPlaneMidRap(in.fhJetEventPlaneMidRap),
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
  fhMultiplicity = in.fhMultiplicity;
  fhMultiplicityDijet = in.fhMultiplicityDijet;
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
  fhJetEventPlaneForwardRap = in.fhJetEventPlaneForwardRap;
  fhJetEventPlaneMidRap = in.fhJetEventPlaneMidRap;
  fCard = in.fCard;
  
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
  delete fhMultiplicity;
  delete fhMultiplicityDijet;
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
  
  // Additional histograms for event plane study
  delete fhJetEventPlaneForwardRap;
  delete fhJetEventPlaneMidRap;
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
  
  // Binning for Q-vector selection
  const Double_t minQvector = 0;         // Minimun normalized Q-vector value
  const Double_t maxQvector = 5;         // Maximun normalized Q-vector value
  const Int_t nBinsQvector = 9;          // Number of Q-vector magnitude bins
  
  // Binning for multiplicity
  const Double_t minMultiplicity = 0;
  const Double_t maxMultiplicity = 4000;
  const Int_t nMultiplicityBins = 400;
  
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
  
//  // Extra variables for jet smearing study
//
//  // Jet pT uncertainty
//  const Double_t minJetPtUncertainty = 0.04;     // Minimum jet pT uncertainty
//  const Double_t maxJetPtUncertainty = 0.14;   // Maximum jet pT uncertainty
//  const Int_t nBinsJetPtUncertainty = 50;     // Number of jet pT uncertainty bins
  
  // Smearing study variables complete!
  
  // Arrays for creating THnSparses
  const Int_t nAxesMultiplicity = 2;
  Int_t nBinsMultiplicity[nAxesMultiplicity];
  Double_t lowBinBorderMultiplicity[nAxesMultiplicity];
  Double_t highBinBorderMultiplicity[nAxesMultiplicity];
  
  const Int_t nAxesJets = 6;   // 6 is nominal, 3 more added for smearing study
  Int_t nBinsJets[nAxesJets];
  Double_t lowBinBorderJets[nAxesJets];
  Double_t highBinBorderJets[nAxesJets];
  
  const Int_t nAxesDijet = 6; // 6 is nominal, 4 more added for xj study
  Int_t nBinsDijet[nAxesDijet];
  Double_t lowBinBorderDijet[nAxesDijet];
  Double_t highBinBorderDijet[nAxesDijet];
  
  const Int_t nAxesAnyJet = 5; // 5 is nominal, 3 more added for smearing study
  Int_t nBinsAnyJet[nAxesAnyJet];
  Double_t lowBinBorderAnyJet[nAxesAnyJet];
  Double_t highBinBorderAnyJet[nAxesAnyJet];
  
  const Int_t nAxesTrack = 5;
  Int_t nBinsTrack[nAxesTrack];
  Double_t lowBinBorderTrack[nAxesTrack];
  Double_t highBinBorderTrack[nAxesTrack];
  
  const Int_t nAxesJetTrack = 7;
  Int_t nBinsJetTrack[nAxesJetTrack];
  Double_t lowBinBorderJetTrack[nAxesJetTrack];
  Double_t highBinBorderJetTrack[nAxesJetTrack];
  
  const Int_t nAxesJetTrackInclusive = 6;
  Int_t nBinsJetTrackInclusive[nAxesJetTrackInclusive];
  Double_t lowBinBorderJetTrackInclusive[nAxesJetTrackInclusive];
  Double_t highBinBorderJetTrackInclusive[nAxesJetTrackInclusive];
  
  const Int_t nAxesJetClosure = 8;  // 7 is nominal, 1 more added for xj study
  Int_t nBinsJetClosure[nAxesJetClosure];
  Double_t lowBinBorderJetClosure[nAxesJetClosure];
  Double_t highBinBorderJetClosure[nAxesJetClosure];
  
  const Int_t nAxesAdditional = 3;
  Int_t nBinsAdditional[nAxesAdditional];
  Double_t lowBinBorderAdditional[nAxesAdditional];
  Double_t highBinBorderAdditional[nAxesAdditional];
  
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
  
  // ======== THnSparses for multiplicity ========
  
  // Axis 0 for the multiplicity histogram: multiplicity
  nBinsMultiplicity[0] = nMultiplicityBins;       // nBins for leading/subleading jet pT
  lowBinBorderMultiplicity[0] = minMultiplicity;  // low bin border for leading/subleading jet pT
  highBinBorderMultiplicity[0] = maxMultiplicity; // high bin border for leading/subleading jet pT
  
  // Axis 1 for the multiplicity histogram: centrality
  nBinsMultiplicity[1] = nWideCentralityBins;     // nBins for wide centrality bins
  lowBinBorderMultiplicity[1] = minCentrality;    // low bin border for centrality
  highBinBorderMultiplicity[1] = maxCentrality;   // high bin border for centrality
  
  // Create the histograms for leading and subleading jets using the above binning information
  fhMultiplicity = new THnSparseF("multiplicity", "multiplicity", nAxesMultiplicity, nBinsMultiplicity, lowBinBorderMultiplicity, highBinBorderMultiplicity); fhMultiplicity->Sumw2();
  fhMultiplicityDijet = new THnSparseF("multiplicityDijet", "multiplicityDijet", nAxesMultiplicity, nBinsMultiplicity, lowBinBorderMultiplicity, highBinBorderMultiplicity); fhMultiplicityDijet->Sumw2();
  
  // Set custom centrality bins for histograms
  fhMultiplicity->SetBinEdges(1,wideCentralityBins);
  fhMultiplicityDijet->SetBinEdges(1,wideCentralityBins);
  
  
  // ======== THnSparses for leading and subleading jets ========
  
  // Axis 0 for the jet histogram: leading/subleading jet pT
  nBinsJets[0] = nPtBinsJet;         // nBins for leading/subleading jet pT
  lowBinBorderJets[0] = minPtJet;    // low bin border for leading/subleading jet pT
  highBinBorderJets[0] = maxPtJet;   // high bin border for leading/subleading jet pT
  
  // Axis 1 for the jet histogram: leading/subleading jet phi
  nBinsJets[1] = nPhiBins;        // nBins for leading/subleading jet phi
  lowBinBorderJets[1] = minPhi;   // low bin border for leading/subleading jet phi
  highBinBorderJets[1] = maxPhi;  // high bin border for leading/subleading jet phi
  
  // Axis 2 for the jet histogram: leading/subleading jet eta
  nBinsJets[2] = nEtaBins;        // nBins for leading/subleading jet eta
  lowBinBorderJets[2] = minEta;   // low bin border for leading/subleading jet eta
  highBinBorderJets[2] = maxEta;  // high bin border for leading/subleading jet eta
  
  // Axis 3 for the jet histogram: asymmetry
  nBinsJets[3] = nWideAsymmetryBins;     // nBins for wide asymmetry bins
  lowBinBorderJets[3] = minXj;           // low bin border for asymmetry
  highBinBorderJets[3] = maxXj;          // high bin border for asymmetry
  
  // Axis 4 for the jet histogram: centrality
  nBinsJets[4] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderJets[4] = minCentrality;  // low bin border for centrality
  highBinBorderJets[4] = maxCentrality; // high bin border for centrality
  
  // Axis 5 for the jet histogram: jet flavor (quark/gluon)
  nBinsJets[5] = nClosureParticleTypeBins;        // nBins for jet flavor
  lowBinBorderJets[5] = minClosureParticleType;   // low bin border for jet flavor
  highBinBorderJets[5] = maxClosureParticleType;  // high bin border for jet flavor
  
//  // Extra axes for smearing study:
//
//  // Axis 6 for the jet histogram: leading/subleading jet pT (smeared)
//  nBinsJets[6] = nPtBinsJet;         // nBins for leading/subleading jet pT (smeared)
//  lowBinBorderJets[6] = minPtJet;    // low bin border for leading/subleading jet pT (smeared)
//  highBinBorderJets[6] = maxPtJet;   // high bin border for leading/subleading jet pT (smeared)
//
//  // Axis 7 for the jet histogram: leading/subleading jet pT minus uncertainty
//  nBinsJets[7] = nBinsJetPtUncertainty;         // nBins for leading/subleading jet pT
//  lowBinBorderJets[7] = minJetPtUncertainty;    // low bin border for leading/subleading jet pT
//  highBinBorderJets[7] = maxJetPtUncertainty;   // high bin border for leading/subleading jet pT
//
//  // Axis 8 for the jet histogram: leading/subleading jet pT plus uncertainty
//  nBinsJets[8] = nBinsJetPtUncertainty;         // nBins for leading/subleading jet pT
//  lowBinBorderJets[8] = minJetPtUncertainty;    // low bin border for leading/subleading jet pT
//  highBinBorderJets[8] = maxJetPtUncertainty;   // high bin border for leading/subleading jet pT
//
//  // Smearing study
  
  
  // Create the histograms for leading and subleading jets using the above binning information
  fhLeadingDijet = new THnSparseF("leadingJet","leadingJet",nAxesJets,nBinsJets,lowBinBorderJets,highBinBorderJets); fhLeadingDijet->Sumw2();
  fhSubleadingDijet = new THnSparseF("subleadingJet","subleadingJet",nAxesJets,nBinsJets,lowBinBorderJets,highBinBorderJets); fhSubleadingDijet->Sumw2();
  
  // Set custom dijet asymmetry bins for histograms
  fhLeadingDijet->SetBinEdges(3,wideAsymmetryBins);
  fhSubleadingDijet->SetBinEdges(3,wideAsymmetryBins);
  
  // Set custom centrality bins for histograms
  fhLeadingDijet->SetBinEdges(4,wideCentralityBins);
  fhSubleadingDijet->SetBinEdges(4,wideCentralityBins);

  // ========= THnSparse for dijets =========
  
  // Axis 0 for the dijet histogram: leading jet pT
  nBinsDijet[0] = nPtBinsJet;         // nBins for leading jet pT
  lowBinBorderDijet[0] = minPtJet;    // low bin border for leading jet pT
  highBinBorderDijet[0] = maxPtJet;   // high bin border for leading jet pT
  
  // Axis 1 for the dijet histogram: subleading jet pT
  nBinsDijet[1] = nPtBinsJet;         // nBins for subleading jet pT
  lowBinBorderDijet[1] = minPtJet;    // low bin border for subleading jet pT
  highBinBorderDijet[1] = maxPtJet;   // high bin border for subleading jet pT
  
  // Axis 2 for the dijet histogram: deltaPhi
  nBinsDijet[2] = nDeltaPhiBins;       // nBins for deltaPhi
  lowBinBorderDijet[2] = minDeltaPhi;  // low bin border for deltaPhi
  highBinBorderDijet[2] = maxDeltaPhi; // high bin border for deltaPhi
  
  // Axis 3 for the dijet histogram: dijet momentum balance xj
  nBinsDijet[3] = nXjBins;         // nBins for xj
  lowBinBorderDijet[3] = minXj;    // low bin border for xj
  highBinBorderDijet[3] = maxXj;   // high bin border for xj
  
  // Axis 4 for the dijet histogram: centrality
  nBinsDijet[4] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderDijet[4] = minCentrality;  // low bin border for centrality
  highBinBorderDijet[4] = maxCentrality; // high bin border for centrality
  
  // Axis 5 for the dijet histogram: matched dijet momentum balance xj
  nBinsDijet[5] = nXjBins;         // nBins for matched xj
  lowBinBorderDijet[5] = minXj;    // low bin border for matched xj
  highBinBorderDijet[5] = maxXj;   // high bin border for matched xj
  
//  // Extra axes for xj study:
//
//  // Axis 6 for the jet histogram: leading jet phi
//  nBinsDijet[6] = nPhiBins;        // nBins for leading jet phi
//  lowBinBorderDijet[6] = minPhi;   // low bin border for leading jet phi
//  highBinBorderDijet[6] = maxPhi;  // high bin border for leading jet phi
//
//  // Axis 7 for the jet histogram: leading jet eta
//  nBinsDijet[7] = nEtaBins;        // nBins for leading jet eta
//  lowBinBorderDijet[7] = minEta;   // low bin border for leading jet eta
//  highBinBorderDijet[7] = maxEta;  // high bin border for leading jet eta
//
//  // Axis 8 for the jet histogram: subleading jet phi
//  nBinsDijet[8] = nPhiBins;        // nBins for subleading jet phi
//  lowBinBorderDijet[8] = minPhi;   // low bin border for subleading jet phi
//  highBinBorderDijet[8] = maxPhi;  // high bin border for subleading jet phi
//
//  // Axis 9 for the jet histogram: subleading jet eta
//  nBinsDijet[9] = nEtaBins;        // nBins for subleading jet eta
//  lowBinBorderDijet[9] = minEta;   // low bin border for subleading jet eta
//  highBinBorderDijet[9] = maxEta;  // high bin border for subleading jet eta
//
//  // xj study
  
  // Create the dijet histogram using the above binning information
  fhDijet = new THnSparseF("dijet","dijet",nAxesDijet,nBinsDijet,lowBinBorderDijet,highBinBorderDijet); fhDijet->Sumw2();
  
  // Set custom centrality bins for histograms
  fhDijet->SetBinEdges(4,wideCentralityBins);
  
  // ======== THnSparse for all jets ========
  
  // Axis 0 for the any jet histogram: jet pT
  nBinsAnyJet[0] = nPtBinsJet;         // nBins for any jet pT
  lowBinBorderAnyJet[0] = minPtJet;    // low bin border for any jet pT
  highBinBorderAnyJet[0] = maxPtJet;   // high bin border for any jet pT
  
  // Axis 1 for the any jet histogram: jet phi
  nBinsAnyJet[1] = nPhiBins;        // nBins for any jet phi
  lowBinBorderAnyJet[1] = minPhi;   // low bin border for any jet phi
  highBinBorderAnyJet[1] = maxPhi;  // high bin border for any jet phi
  
  // Axis 2 for the any jet histogram: jet eta
  nBinsAnyJet[2] = nEtaBins;        // nBins for any jet eta
  lowBinBorderAnyJet[2] = minEta;   // low bin border for any jet eta
  highBinBorderAnyJet[2] = maxEta;  // high bin border for any jet eta
  
  // Axis 3 for the any jet histogram: centrality
  nBinsAnyJet[3] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderAnyJet[3] = minCentrality;  // low bin border for centrality
  highBinBorderAnyJet[3] = maxCentrality; // high bin border for centrality
  
  // Axis 4 for the jet histogram: jet flavor (quark/gluon)
  nBinsAnyJet[4] = nClosureParticleTypeBins;        // nBins for jet flavor
  lowBinBorderAnyJet[4] = minClosureParticleType;   // low bin border for jet flavor
  highBinBorderAnyJet[4] = maxClosureParticleType;  // high bin border for jet flavor
  
//  // Extra axes for smearing study:
//
//  // Axis 5 for the jet histogram: any jet pT (smeared)
//  nBinsAnyJet[5] = nPtBinsJet;         // nBins for any jet pT (smeared)
//  lowBinBorderAnyJet[5] = minPtJet;    // low bin border for any jet pT (smeared)
//  highBinBorderAnyJet[5] = maxPtJet;   // high bin border for any jet pT (smeared)
//
//  // Axis 6 for the jet histogram: any jet pT minus uncertainty
//  nBinsAnyJet[6] = nBinsJetPtUncertainty;         // nBins for any jet pT
//  lowBinBorderAnyJet[6] = minJetPtUncertainty;    // low bin border for any jet pT
//  highBinBorderAnyJet[6] = maxJetPtUncertainty;   // high bin border for any jet pT
//
//  // Axis 7 for the jet histogram: any jet pT plus uncertainty
//  nBinsAnyJet[7] = nBinsJetPtUncertainty;         // nBins for any jet pT
//  lowBinBorderAnyJet[7] = minJetPtUncertainty;    // low bin border for any jet pT
//  highBinBorderAnyJet[7] = maxJetPtUncertainty;   // high bin border for any jet pT
//
//  // Smearing study
  
  // Create the histogram for all jets using the above binning information
  fhAnyJet = new THnSparseF("anyJet","anyJet",nAxesAnyJet,nBinsAnyJet,lowBinBorderAnyJet,highBinBorderAnyJet); fhAnyJet->Sumw2();
  fhLeadingJet = new THnSparseF("anyLeadingJet","anyLeadingJet",nAxesAnyJet,nBinsAnyJet,lowBinBorderAnyJet,highBinBorderAnyJet); fhLeadingJet->Sumw2();

  // Set custom centrality bins for histograms
  fhAnyJet->SetBinEdges(3,wideCentralityBins);
  fhLeadingJet->SetBinEdges(3,wideCentralityBins);
  
  // ======== THnSparses for tracks and uncorrected tracks ========
  
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
  
  // Axis 3 for the track histogram: centrality
  nBinsTrack[3] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderTrack[3] = minCentrality;  // low bin border for centrality
  highBinBorderTrack[3] = maxCentrality; // high bin border for centrality
  
  // Axis 4 for the track histogram: correlation type
  nBinsTrack[4] = nCorrelationTypeBins;        // nBins for correlation types
  lowBinBorderTrack[4] = minCorrelationType;   // low bin border for correlation types
  highBinBorderTrack[4] = maxCorrelationType;  // high bin border for correlation types
  
  // Create the histograms for tracks and uncorrected tracks using the above binning information
  fhTrack = new THnSparseF("track","track",nAxesTrack,nBinsTrack,lowBinBorderTrack,highBinBorderTrack); fhTrack->Sumw2();
  fhTrackInclusive = new THnSparseF("trackInclusive","trackInclusive",nAxesTrack-1,nBinsTrack,lowBinBorderTrack,highBinBorderTrack); fhTrackInclusive->Sumw2();
  fhTrackUncorrected = new THnSparseF("trackUncorrected","trackUncorrected",nAxesTrack,nBinsTrack,lowBinBorderTrack,highBinBorderTrack); fhTrackUncorrected->Sumw2();
  fhTrackInclusiveUncorrected = new THnSparseF("trackInclusiveUncorrected","trackInclusiveUncorrected",nAxesTrack-1,nBinsTrack,lowBinBorderTrack,highBinBorderTrack); fhTrackInclusiveUncorrected->Sumw2();

  // Set custom centrality bins for histograms
  fhTrack->SetBinEdges(3,wideCentralityBins);
  fhTrackInclusive->SetBinEdges(3,wideCentralityBins);
  fhTrackUncorrected->SetBinEdges(3,wideCentralityBins);
  fhTrackInclusiveUncorrected->SetBinEdges(3,wideCentralityBins);
  
  // ======== THnSparses for correlation between tracks and leading or subleading jets ========
  
  // Axis 0 for the track-jet correlation histogram: track pT
  nBinsJetTrack[0] = nWideTrackPtBins;     // nBins for wide track pT bins
  lowBinBorderJetTrack[0] = minPtTrack;    // low bin border for track pT
  highBinBorderJetTrack[0] = maxPtTrack;   // high bin border for track pT
  
  // Axis 1 for the track-jet correlation histogram: deltaPhi between track and jet
  nBinsJetTrack[1] = nDeltaPhiBinsJetTrack;       // nBins for deltaPhi between track and jet
  lowBinBorderJetTrack[1] = minDeltaPhiJetTrack;  // low bin border for deltaPhi between track and jet
  highBinBorderJetTrack[1] = maxDeltaPhiJetTrack; // high bin border for deltaPhi between track and jet
  
  // Axis 2 for the track-jet correlation histogram: deltaEta between track and jet
  nBinsJetTrack[2] = nDeltaEtaBinsJetTrack;         // nBins for deltaEta between track and jet
  lowBinBorderJetTrack[2] = minDeltaEtaJetTrack;    // low bin border deltaEta between track and jet
  highBinBorderJetTrack[2] = maxDeltaEtaJetTrack;   // high bin border deltaEta between track and jet
  
  // Axis 3 for the track-jet correlation histogram: dijet asymmetry
  nBinsJetTrack[3] = nWideAsymmetryBins;     // nBins for wide dijet asymmetry bins
  lowBinBorderJetTrack[3] = minXj;           // low bin border for dijet asymmetry
  highBinBorderJetTrack[3] = maxXj;          // high bin border for dijet asymmetry
  
  // Axis 4 for the track-jet correlation histogram: centrality
  nBinsJetTrack[4] = nWideCentralityBins;     // nBins for centrality
  lowBinBorderJetTrack[4] = minCentrality;    // low bin border for centrality
  highBinBorderJetTrack[4] = maxCentrality;   // high bin border for centrality
  
  // Axis 5 for the track-jet correlation histogram: correlation type
  nBinsJetTrack[5] = nCorrelationTypeBins;        // nBins for correlation types
  lowBinBorderJetTrack[5] = minCorrelationType;   // low bin border for correlation types
  highBinBorderJetTrack[5] = maxCorrelationType;  // high bin border for correlation types
  
  // Axis 6 for the track-jet correlation histogram: jet flavor (quark/gluon)
  nBinsJetTrack[6] = nClosureParticleTypeBins;        // nBins for jet flavors
  lowBinBorderJetTrack[6] = minClosureParticleType;   // low bin border for jet flavor
  highBinBorderJetTrack[6] = maxClosureParticleType;  // high bin border for jet flavor
  
  // Create histograms for track-jet correlations
  fhTrackLeadingJet = new THnSparseF("trackLeadingJet", "trackLeadingJet", nAxesJetTrack, nBinsJetTrack, lowBinBorderJetTrack, highBinBorderJetTrack); fhTrackLeadingJet->Sumw2();
  fhTrackSubleadingJet = new THnSparseF("trackSubleadingJet", "trackSubleadingJet", nAxesJetTrack, nBinsJetTrack, lowBinBorderJetTrack, highBinBorderJetTrack); fhTrackSubleadingJet->Sumw2();
  
  // Set custom centrality, asymmetry and track pT bins for histograms
  fhTrackLeadingJet->SetBinEdges(0,wideTrackPtBins);
  fhTrackSubleadingJet->SetBinEdges(0,wideTrackPtBins);
  fhTrackLeadingJet->SetBinEdges(3,wideAsymmetryBins);
  fhTrackSubleadingJet->SetBinEdges(3,wideAsymmetryBins);
  fhTrackLeadingJet->SetBinEdges(4,wideCentralityBins);
  fhTrackSubleadingJet->SetBinEdges(4,wideCentralityBins);

  // Create histograms for uncorrected track-jet correlations
  fhTrackLeadingJetUncorrected = new THnSparseF("trackLeadingJetUncorrected", "trackLeadingJetUncorrected", nAxesJetTrack, nBinsJetTrack, lowBinBorderJetTrack, highBinBorderJetTrack); fhTrackLeadingJetUncorrected->Sumw2();
  fhTrackSubleadingJetUncorrected = new THnSparseF("trackSubleadingJetUncorrected", "trackSubleadingJetUncorrected", nAxesJetTrack, nBinsJetTrack, lowBinBorderJetTrack, highBinBorderJetTrack); fhTrackSubleadingJetUncorrected->Sumw2();
  
  // Set custom centrality, asymmetry and track pT bins for histograms
  fhTrackLeadingJetUncorrected->SetBinEdges(0,wideTrackPtBins);
  fhTrackSubleadingJetUncorrected->SetBinEdges(0,wideTrackPtBins);
  fhTrackLeadingJetUncorrected->SetBinEdges(3,wideAsymmetryBins);
  fhTrackSubleadingJetUncorrected->SetBinEdges(3,wideAsymmetryBins);
  fhTrackLeadingJetUncorrected->SetBinEdges(4,wideCentralityBins);
  fhTrackSubleadingJetUncorrected->SetBinEdges(4,wideCentralityBins);
  
  // Create histograms for pT weighted track-jet correlations
  fhTrackLeadingJetPtWeighted = new THnSparseF("trackLeadingJetPtWeighted", "trackLeadingJetPtWeighted", nAxesJetTrack, nBinsJetTrack, lowBinBorderJetTrack, highBinBorderJetTrack); fhTrackLeadingJetPtWeighted->Sumw2();
  fhTrackSubleadingJetPtWeighted = new THnSparseF("trackSubleadingJetPtWeighted", "trackSubleadingJetPtWeighted", nAxesJetTrack, nBinsJetTrack, lowBinBorderJetTrack, highBinBorderJetTrack); fhTrackSubleadingJetPtWeighted->Sumw2();
  
  // Set custom centrality, asymmetry and track pT bins for histograms
  fhTrackLeadingJetPtWeighted->SetBinEdges(0,wideTrackPtBins);
  fhTrackSubleadingJetPtWeighted->SetBinEdges(0,wideTrackPtBins);
  fhTrackLeadingJetPtWeighted->SetBinEdges(3,wideAsymmetryBins);
  fhTrackSubleadingJetPtWeighted->SetBinEdges(3,wideAsymmetryBins);
  fhTrackLeadingJetPtWeighted->SetBinEdges(4,wideCentralityBins);
  fhTrackSubleadingJetPtWeighted->SetBinEdges(4,wideCentralityBins);
  
  // ======== THnSparses for correlation between tracks and inclusive jets ========
  
  // Axis 0 for the track-inclusive jet correlation histogram: track pT
  nBinsJetTrackInclusive[0] = nWideTrackPtBins;     // nBins for wide track pT bins
  lowBinBorderJetTrackInclusive[0] = minPtTrack;    // low bin border for track pT
  highBinBorderJetTrackInclusive[0] = maxPtTrack;   // high bin border for track pT
  
  // Axis 1 for the track-inclusive jet correlation histogram: deltaPhi between track and jet
  nBinsJetTrackInclusive[1] = nDeltaPhiBinsJetTrack;       // nBins for deltaPhi between track and jet
  lowBinBorderJetTrackInclusive[1] = minDeltaPhiJetTrack;  // low bin border for deltaPhi between track and jet
  highBinBorderJetTrackInclusive[1] = maxDeltaPhiJetTrack; // high bin border for deltaPhi between track and jet
  
  // Axis 2 for the track-inclusive jet correlation histogram: deltaEta between track and jet
  nBinsJetTrackInclusive[2] = nDeltaEtaBinsJetTrack;         // nBins for deltaEta between track and jet
  lowBinBorderJetTrackInclusive[2] = minDeltaEtaJetTrack;    // low bin border deltaEta between track and jet
  highBinBorderJetTrackInclusive[2] = maxDeltaEtaJetTrack;   // high bin border deltaEta between track and jet
  
  // Axis 3 for the track-inclusive jet correlation histogram: centrality
  nBinsJetTrackInclusive[3] = nWideCentralityBins;     // nBins for centrality
  lowBinBorderJetTrackInclusive[3] = minCentrality;    // low bin border for centrality
  highBinBorderJetTrackInclusive[3] = maxCentrality;   // high bin border for centrality
  
  // Axis 4 for the track-inclusive jet correlation histogram: correlation type
  nBinsJetTrackInclusive[4] = nCorrelationTypeBins;        // nBins for correlation types
  lowBinBorderJetTrackInclusive[4] = minCorrelationType;   // low bin border for correlation types
  highBinBorderJetTrackInclusive[4] = maxCorrelationType;  // high bin border for correlation types
  
  // Axis 5 for the track-inclusive jet correlation histogram: jet flavor (quark/gluon)
  nBinsJetTrackInclusive[5] = nClosureParticleTypeBins;        // nBins for jet flavors
  lowBinBorderJetTrackInclusive[5] = minClosureParticleType;   // low bin border for jet flavor
  highBinBorderJetTrackInclusive[5] = maxClosureParticleType;  // high bin border for jet flavor
  
  // Create histograms for track-inclusive jet correlations
  fhTrackJetInclusive = new THnSparseF("trackJetInclusive", "trackJetInclusive", nAxesJetTrackInclusive, nBinsJetTrackInclusive, lowBinBorderJetTrackInclusive, highBinBorderJetTrackInclusive); fhTrackJetInclusive->Sumw2();
  fhTrackJetInclusivePtWeighted = new THnSparseF("trackJetInclusivePtWeighted", "trackJetInclusivePtWeighted", nAxesJetTrackInclusive, nBinsJetTrackInclusive, lowBinBorderJetTrackInclusive, highBinBorderJetTrackInclusive); fhTrackJetInclusivePtWeighted->Sumw2();
  
  // Set custom centrality and track pT bins for histograms
  fhTrackJetInclusive->SetBinEdges(0,wideTrackPtBins);
  fhTrackJetInclusivePtWeighted->SetBinEdges(0,wideTrackPtBins);
  fhTrackJetInclusive->SetBinEdges(3,wideCentralityBins);
  fhTrackJetInclusivePtWeighted->SetBinEdges(3,wideCentralityBins);
  
  // ======== THnSparses for jet pT closures ========
  
  // Axis 0 for the jet pT closure histogram: closure type (leading/subleading/inclusive)
  nBinsJetClosure[0] = nClosureTypeBins;         // nBins for closure types
  lowBinBorderJetClosure[0] = minClosureType;    // low bin border for closure types
  highBinBorderJetClosure[0] = maxClosureType;   // high bin border for closure types
  
  // Axis 1 for the jet pT closure histogram: generator level jet pT
  nBinsJetClosure[1] = nClosurePtBins;       // nBins for generator level pT bins in closure plots
  lowBinBorderJetClosure[1] = minClosurePt;  // low bin border generator level pT in closure plots
  highBinBorderJetClosure[1] = maxClosurePt; // high bin border generator level pT in closure plots
  
  // Axis 2 for the jet pT closure histogram: reconstructed jet pT
  nBinsJetClosure[2] = nClosurePtBins;       // nBins for reconstructed jet pT bins in closure plots
  lowBinBorderJetClosure[2] = minClosurePt;  // low bin border for reconstructed jet pT in closure plots
  highBinBorderJetClosure[2] = maxClosurePt; // high bin border for reconstructed jet pT in closure plots
  
  // Axis 3 for the jet pT closure histogram: generator level jet eta
  nBinsJetClosure[3] = nEtaBins;             // nBins for jet eta
  lowBinBorderJetClosure[3] = minEta;        // low bin border for jet eta
  highBinBorderJetClosure[3] = maxEta;       // high bin border for jet eta
  
  // Axis 4 for the jet pT closure histogram: centrality
  nBinsJetClosure[4] = nWideCentralityBins;     // nBins for centrality
  lowBinBorderJetClosure[4] = minCentrality;    // low bin border for centrality
  highBinBorderJetClosure[4] = maxCentrality;   // high bin border for centrality
  
  // Axis 5 for the jet pT closure histogram: ref parton = quark/gluon
  nBinsJetClosure[5] = nClosureParticleTypeBins;         // nBins for reference parton
  lowBinBorderJetClosure[5] = minClosureParticleType;    // low bin border for reference parton
  highBinBorderJetClosure[5] = maxClosureParticleType;   // high bin border for reference parton
  
  // Axis 6 for the jet pT closure histogram: reco/gen ratio for closure
  nBinsJetClosure[6] = nClosureRatioBins;        // nBins for closure ratio
  lowBinBorderJetClosure[6] = minClosureRatio;   // low bin border for closure ratio
  highBinBorderJetClosure[6] = maxClosureRatio;  // high bin border for closure ratio
  
  // Extra axes added for xj study
  
  // Axis 7 for the jet pT closure histogram: dijet momentum balance xj
  //nBinsJetClosure[7] = nWideAsymmetryBins;     // nBins for wide dijet asymmetry bins
  //lowBinBorderJetClosure[7] = minXj;           // low bin border for dijet asymmetry
  //highBinBorderJetClosure[7] = maxXj;          // high bin border for dijet asymmetry
  
  // Axis 7 for the jet pT closure histogram: angle with respect to reaction plane
  nBinsJetClosure[7] = nDeltaPhiBinsJetTrack;     // nBins for wide dijet asymmetry bins
  lowBinBorderJetClosure[7] = minDeltaPhiJetTrack;           // low bin border for dijet asymmetry
  highBinBorderJetClosure[7] = maxDeltaPhiJetTrack;          // high bin border for dijet asymmetry
  
  // xj study
  
  // Create histograms for jet pT closure
  fhJetPtClosure = new THnSparseF("jetPtClosure", "jetPtClosure", nAxesJetClosure, nBinsJetClosure, lowBinBorderJetClosure, highBinBorderJetClosure); fhJetPtClosure->Sumw2();
  
  // Set custom centrality bins for histograms
  fhJetPtClosure->SetBinEdges(4,wideCentralityBins);
  
  // Set up wide xj bins for extra axis used in xj study
  //fhJetPtClosure->SetBinEdges(7,wideAsymmetryBins);
  // xj study
  
  // ======== Extra THnSparses for additional event plane study ========
  
  // Axis 1 for the additional histogram: DeltaPhi between jet and event plane angle
  nBinsAdditional[0] = nDeltaPhiBinsJetTrack;         // nBins for deltaPhi between jet and event plane
  lowBinBorderAdditional[0] = minDeltaPhiJetTrack;    // low bin border for deltaPhi between jet and event plane
  highBinBorderAdditional[0] = maxDeltaPhiJetTrack;   // high bin border for deltaPhi between jet and event plane
  
  // Axis 2 for the additional histogram: Normalized Q-vector from mid rapidity event plane
  nBinsAdditional[1] = nBinsQvector;                  // nBins for deltaPhi between jet and event plane
  lowBinBorderAdditional[1] = minQvector;             // low bin border for deltaPhi between jet and event plane
  highBinBorderAdditional[1] = maxQvector;            // high bin border for deltaPhi between jet and event plane
  
  // Axis 3 for the additional histogram: Centrality
  nBinsAdditional[2] = nWideCentralityBins;           // nBins for centrality
  lowBinBorderAdditional[2] = minCentrality;          // low bin border for centrality
  highBinBorderAdditional[2] = maxCentrality;         // high bin border for centrality
  
  // Create additional histograms for event plane study
  fhJetEventPlaneForwardRap = new THnSparseF("jetEventPlaneForwardRap", "jetEventPlaneForwardRap", nAxesAdditional, nBinsAdditional, lowBinBorderAdditional, highBinBorderAdditional); fhJetEventPlaneForwardRap->Sumw2();
  fhJetEventPlaneMidRap = new THnSparseF("jetEventPlaneMidRap", "jetEventPlaneMidRap", nAxesAdditional, nBinsAdditional, lowBinBorderAdditional, highBinBorderAdditional); fhJetEventPlaneMidRap->Sumw2();
  
  // Set custom centrality bins for histograms
  fhJetEventPlaneForwardRap->SetBinEdges(2,wideCentralityBins);
  fhJetEventPlaneMidRap->SetBinEdges(2,wideCentralityBins);
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
  fhMultiplicity->Write();
  fhMultiplicityDijet->Write();
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
  
  // Additional histgorams for the evgent plane study
  fhJetEventPlaneForwardRap->Write();
  fhJetEventPlaneMidRap->Write();
  
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


