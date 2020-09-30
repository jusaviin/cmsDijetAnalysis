// Class for the main analysis algorithms for the leading-subleading jet analysis

// Root includes
#include <TFile.h>
#include <TMath.h>

// Own includes
#include "DijetAnalyzer.h"

using namespace std;

/*
 * Combination of third order polynomial and Gaussian function
 * Functional form: f(x) = a + bx + cx^2 + dx^3 + h/s * e^(-0.5*((x-g)/s)^2)
 */
double jetPolyGauss(double *x, double *par){
  return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]+(par[4]/par[5])*TMath::Exp(-0.5*TMath::Power((x[0]-par[6])/par[5],2));
}

/*
 * Default constructor
 */
DijetAnalyzer::DijetAnalyzer() :
  fFileNames(0),
  fCard(0),
  fHistograms(0),
  fTrackCorrection(),
  fJffCorrection(),
  fVzWeightFunction(0),
  fCentralityWeightFunction(0),
  fPtWeightFunction(0),
  fDijetWeightFunction(0),
  fSmearingFunction(0),
  fTrackEfficiencyCorrector2018(),
  fJetCorrector2018(),
  fJetUncertainty2018(),
  fRng(0),
  fDataType(-1),
  fForestType(0),
  fReadMode(0),
  fJetType(0),
  fMatchJets(false),
  fMatchDijet(false),
  fMatchLeadingJet(false),
  fDebugLevel(0),
  fLocalRun(0),
  fVzWeight(1),
  fCentralityWeight(1),
  fPtHatWeight(1),
  fTotalEventWeight(1),
  fnEventsInMixingFile(0),
  fnMixedEventsPerDijet(0),
  fMixingStartIndex(0),
  fMixingPoolDepth(1),
  fMixingVzBinWidth(1),
  fMixingHiBinWidth(1),
  fRunningMixingIndex(0),
  fMaximumMixingVz(0),
  fMaximumMixingHiBin(0),
  fMixingVzTolerance(0),
  fMixedEventVz(0),
  fMixedEventHiBin(0),
  fJetAxis(0),
  fVzCut(0),
  fMinimumPtHat(0),
  fMaximumPtHat(0),
  fJetEtaCut(0),
  fJetSearchEtaCut(0),
  fJetMaximumPtCut(0),
  fLeadingJetMinPtCut(0),
  fSubleadingJetMinPtCut(0),
  fDeltaPhiCut(0),
  fMinimumMaxTrackPtFraction(0),
  fMaximumMaxTrackPtFraction(0),
  fJetUncertaintyMode(0),
  fTrackEtaCut(0),
  fTrackMinPtCut(0),
  fMaxTrackPtRelativeError(0),
  fMaxTrackDistanceToVertex(0),
  fCalorimeterSignalLimitPt(0),
  fHighPtEtFraction(0),
  fChi2QualityCut(0),
  fMinimumTrackHits(0),
  fSubeventCut(0),
  fMcCorrelationType(0),
  fAsymmetryBinType(0),
  fFillEventInformation(false),
  fFillJetHistograms(false),
  fFillTrackHistograms(false),
  fFillRegularJetTrackCorrelation(false),
  fFillUncorrectedJetTrackCorrelation(false),
  fFillPtWeightedJetTrackCorrelation(false),
  fFillInclusiveJetTrackCorrelation(false),
  fFillJetPtClosure(false),
  fFillDijetJetTrackCorrelation(false)
{
  // Default constructor
  fHistograms = new DijetHistograms();
  fHistograms->CreateHistograms();
  
  // Initialize readers to null
  fJetReader = NULL;
  fTrackReader[DijetHistograms::kSameEvent] = NULL;
  fTrackReader[DijetHistograms::kMixedEvent] = NULL;
  
  // Clear the mixing pool
  for(Int_t iVz = 0; iVz < kMaxMixingVzBins; iVz++){
    for(Int_t iHiBin = 0; iHiBin < kMaxMixingHiBins; iHiBin++){
      fMixingPool[iVz][iHiBin].clear();
    }
  }
}

/*
 * Custom constructor
 */
DijetAnalyzer::DijetAnalyzer(std::vector<TString> fileNameVector, ConfigurationCard *newCard, bool runLocal) :
  fFileNames(fileNameVector),
  fCard(newCard),
  fHistograms(0),
  fJetCorrector2018(),
  fJetUncertainty2018(),
  fForestType(0),
  fJetType(0),
  fMatchJets(false),
  fMatchDijet(false),
  fMatchLeadingJet(false),
  fVzWeight(1),
  fCentralityWeight(1),
  fPtHatWeight(1),
  fTotalEventWeight(1),
  fnEventsInMixingFile(0),
  fMixingStartIndex(0),
  fRunningMixingIndex(0),
  fMixedEventVz(0),
  fMixedEventHiBin(0),
  fMinimumPtHat(0),
  fMaximumPtHat(0),
  fJetEtaCut(0),
  fJetSearchEtaCut(0),
  fJetMaximumPtCut(0),
  fLeadingJetMinPtCut(0),
  fSubleadingJetMinPtCut(0),
  fDeltaPhiCut(0),
  fMinimumMaxTrackPtFraction(0),
  fMaximumMaxTrackPtFraction(0),
  fJetUncertaintyMode(0),
  fTrackEtaCut(0),
  fTrackMinPtCut(0),
  fMaxTrackPtRelativeError(0),
  fMaxTrackDistanceToVertex(0),
  fCalorimeterSignalLimitPt(0),
  fHighPtEtFraction(0),
  fChi2QualityCut(0),
  fMinimumTrackHits(0),
  fSubeventCut(0),
  fMcCorrelationType(0)
{
  // Custom constructor
  fHistograms = new DijetHistograms(fCard);
  fHistograms->CreateHistograms();
  
  // Initialize readers to null
  fJetReader = NULL;
  fTrackReader[DijetHistograms::kSameEvent] = NULL;
  fTrackReader[DijetHistograms::kMixedEvent] = NULL;
  
  // Amount of debugging messages
  fDebugLevel = fCard->Get("DebugLevel");
  
  // Flog for local running or CRAB running
  fLocalRun = runLocal ? 1 : 0;
  
  // Jet axis type
  fJetAxis = fCard->Get("JetAxis");
  
  // vz cut
  fVzCut = fCard->Get("ZVertexCut");          // Event cut vor the z-position of the primary vertex
  
  // Initialize the mixing pool
  for(Int_t iVz = 0; iVz < kMaxMixingVzBins; iVz++){
    for(Int_t iHiBin = 0; iHiBin < kMaxMixingHiBins; iHiBin++){
      fMixingPool[iVz][iHiBin].clear();
    }
  }
  
  // Mixing pool parameters
  fnMixedEventsPerDijet = fCard->Get("NMixedEventsPerDijet");
  fMixingPoolDepth = fCard->Get("MixingPoolDepth");
  fMixingVzBinWidth = fCard->Get("MixingVzBinWidth");
  fMixingHiBinWidth = fCard->Get("MixingHiBinWidth");
  fMixingVzTolerance = fCard->Get("VzTolerance");
  fMaximumMixingVz = FindMixingVzBin(fVzCut);
  fMaximumMixingHiBin = 0;  // This will be changed for heavy ions in the code below
  
  // Asymmetry binning
  fAsymmetryBinType = fCard->Get("AsymmetryBinType");   // 0 = AJ, 1 = xJ
  
  // pT weight function for Pythia to match 2017 MC and data pT spectra
  fPtWeightFunction = new TF1("fPtWeightFunction","pol3",0,500);
  fPtWeightFunction->SetParameters(0.699073,0.00287672,-6.35568e-06,5.02292e-09);
  
  // Weight function derived for leading jet
  fDijetWeightFunction = new TF1("fDijetWeightFunction","pol3",0,500);
  fDijetWeightFunction->SetParameters(0.723161,0.00236126,-3.90984e-06,3.10631e-09);
  
  // Function for smearing the jet pT for systemtic uncertainties
  fSmearingFunction = new TF1("fSmearingFunction","pol4",0,500);
  
  // Weight function derived for subleading jet
  //fDijetWeightFunction = new TF1("fDijetWeightFunction",jetPolyGauss,0,500,7);
  //fDijetWeightFunction->SetParameters(0.683267,2.97273e-03,-6.71989e-06,6.26340e-09,1.81108,16.8066,84.4603);
  
  // Find the correct folder for track correction tables based on data type
  fDataType = fCard->Get("DataType");
  fReadMode = fCard->Get("ReadMode");
  bool ppData = true;
  if(fDataType == ForestReader::kPp || fDataType == ForestReader::kPpMC || fDataType == ForestReader::kLocalTest){
    
    // Track correction for 2015 pp data
    fTrackCorrection = new TrkCorr("trackCorrectionTables/TrkCorr_July22_Iterative_pp_eta2p4/");
    
    // Track correction for 2017 pp data
    fTrackEfficiencyCorrector2018 = new TrkEff2017pp(false, "trackCorrectionTables/pp2017/");
    
    // Common vz weight function used by UIC group for pp MC 2015
    if(fReadMode < 2000){
      fVzWeightFunction = new TF1("fvz","gaus",-15,15);  // Weight function for 2015 MC
      fVzWeightFunction->SetParameter(0,1.10477);
      fVzWeightFunction->SetParameter(1,2.52738);
      fVzWeightFunction->SetParameter(2,1.30296e1);
    } else {
      fVzWeightFunction = new TF1("fvz","pol6",-15,15);  // Weight function for 2017 MC
      fVzWeightFunction->SetParameters(0.973805, 0.00339418, 0.000757544, -1.37331e-06, -2.82953e-07, -3.06778e-10, 3.48615e-09);
    }
    fCentralityWeightFunction = NULL;
    
  } else if (fDataType == ForestReader::kPbPb || fDataType == ForestReader::kPbPbMC){
    
    // Track correction for 2015 PbPb data
    fTrackCorrection = new XiaoTrkCorr("trackCorrectionTables/xiaoCorrection/eta_symmetry_cymbalCorr_FineBin.root");
    
    // Track correction for 2018 PbPb data
    fTrackEfficiencyCorrector2018 = new TrkEff2018PbPb("general", false, "trackCorrectionTables/PbPb2018/");
    
    // Flag for PbPb data
    ppData = false;
    
    // Common vz weight function used by UIC group for PbPb MC
    fVzWeightFunction = new TF1("fvz","pol6",-15,15);
    
    if(fReadMode < 2000){
      fVzWeightFunction->SetParameters(1.18472, -0.132675, 0.00857998, -0.000326085, -1.48786e-06, 4.68665e-07, -7.32942e-09); // Parameters for 2015 MC
    } else {
      fVzWeightFunction->SetParameters(1.00656, -0.0193651, 0.000976851, -1.043e-05, -9.79808e-06, 9.07733e-08, 1.79165e-08); // Parameters for 2018 MC
    }
    
    // Comment centrality weight function used by UIC group for PbPb MC
    
    if(fReadMode < 2000){ // Weight function for 2015 MC
      fCentralityWeightFunction = new TF1("fcent1", "[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4+[7]*exp([5]+[6]*x)", 0, 180); // 2015
      fCentralityWeightFunction->SetParameters(4.40810, -7.75301e-02, 4.91953e-04, -1.34961e-06, 1.44407e-09, -160, 1, 3.68078e-15); // 2015
    } else { // Weight function for 2018 MC
      fCentralityWeightFunction = new TF1("fcent","pol6",0,90); // 2018
      fCentralityWeightFunction->SetParameters(4.64945,-0.201337, 0.00435794,-7.00799e-05,8.18299e-07,-5.52604e-09,1.54472e-11); // 2018
    }
    
    // Set the number of HiBins for mixing pool
    fMaximumMixingHiBin = FindMixingHiBin(200);  // 200 is the maximum number for HiBin in the forests
    
  } else {
    fTrackCorrection = new TrkCorr(""); // Bad data type, no corrections initialized
    fVzWeightFunction = NULL;
    fCentralityWeightFunction = NULL;
  }
  
  // Initialize the class for JFF correction
  fJffCorrection = new JffCorrection(ppData);
  
  // Read which histograms to fill this run
  int filledHistograms = fCard->Get("FilledHistograms");
  std::bitset<knFillTypes> bitChecker(filledHistograms);
  fFillEventInformation = bitChecker.test(kFillEventInformation);
  fFillJetHistograms = bitChecker.test(kFillJets);
  fFillTrackHistograms = bitChecker.test(kFillTracks);
  fFillRegularJetTrackCorrelation = bitChecker.test(kFillRegularJetTrackCorrelation);
  fFillUncorrectedJetTrackCorrelation = bitChecker.test(kFillUncorrectedJetTrackCorrelation);
  fFillPtWeightedJetTrackCorrelation = bitChecker.test(kFillPtWeightedJetTrackCorrelation);
  fFillInclusiveJetTrackCorrelation = bitChecker.test(kFillInclusiveJetTrackCorrelation);
  fFillJetPtClosure = bitChecker.test(kFillJetPtClosure);
  fFillDijetJetTrackCorrelation = (fFillRegularJetTrackCorrelation || fFillPtWeightedJetTrackCorrelation || fFillPtWeightedJetTrackCorrelation || fFillTrackHistograms);
  
  // Do a sanity check for given bin widths for mixing pool
  if(fMaximumMixingHiBin >= kMaxMixingHiBins){
    cout << "Error! There are more than allowed number of hiBins in mixing pool! Please increse bin width in JCard!" << endl;
    assert(0);
  }
  
  if(fMaximumMixingVz >= kMaxMixingVzBins){
    cout << "Error! There are more than allowed number of vz in mixing pool! Please increse bin width in JCard!" << endl;
    assert(0);
  }
  
  // Initialize the random number generator with a random seed
  fRng = new TRandom3();
  fRng->SetSeed(0);
  
}

/*
 * Copy constructor
 */
DijetAnalyzer::DijetAnalyzer(const DijetAnalyzer& in) :
  fJetReader(in.fJetReader),
  fFileNames(in.fFileNames),
  fCard(in.fCard),
  fHistograms(in.fHistograms),
  fTrackCorrection(in.fTrackCorrection),
  fJffCorrection(in.fJffCorrection),
  fVzWeightFunction(in.fVzWeightFunction),
  fCentralityWeightFunction(in.fCentralityWeightFunction),
  fPtWeightFunction(in.fPtWeightFunction),
  fDijetWeightFunction(in.fDijetWeightFunction),
  fSmearingFunction(in.fSmearingFunction),
  fRng(in.fRng),
  fDataType(in.fDataType),
  fForestType(in.fForestType),
  fReadMode(in.fReadMode),
  fJetType(in.fJetType),
  fMatchJets(in.fMatchJets),
  fMatchDijet(in.fMatchDijet),
  fMatchLeadingJet(in.fMatchLeadingJet),
  fDebugLevel(in.fDebugLevel),
  fLocalRun(in.fLocalRun),
  fVzWeight(in.fVzWeight),
  fCentralityWeight(in.fCentralityWeight),
  fPtHatWeight(in.fPtHatWeight),
  fTotalEventWeight(in.fTotalEventWeight),
  fnEventsInMixingFile(in.fnEventsInMixingFile),
  fnMixedEventsPerDijet(in.fnMixedEventsPerDijet),
  fMixingStartIndex(in.fMixingStartIndex),
  fMixingPoolDepth(in.fMixingPoolDepth),
  fMixingVzBinWidth(in.fMixingVzBinWidth),
  fMixingHiBinWidth(in.fMixingHiBinWidth),
  fRunningMixingIndex(in.fRunningMixingIndex),
  fMixingVzTolerance(in.fMixingVzTolerance),
  fMixedEventVz(in.fMixedEventVz),
  fMixedEventHiBin(in.fMixedEventHiBin),
  fJetAxis(in.fJetAxis),
  fVzCut(in.fVzCut),
  fMinimumPtHat(in.fMinimumPtHat),
  fMaximumPtHat(in.fMaximumPtHat),
  fJetEtaCut(in.fJetEtaCut),
  fJetSearchEtaCut(in.fJetSearchEtaCut),
  fJetMaximumPtCut(in.fJetMaximumPtCut),
  fLeadingJetMinPtCut(in.fLeadingJetMinPtCut),
  fSubleadingJetMinPtCut(in.fSubleadingJetMinPtCut),
  fDeltaPhiCut(in.fDeltaPhiCut),
  fMinimumMaxTrackPtFraction(in.fMinimumMaxTrackPtFraction),
  fMaximumMaxTrackPtFraction(in.fMaximumMaxTrackPtFraction),
  fJetUncertaintyMode(in.fJetUncertaintyMode),
  fTrackEtaCut(in.fTrackEtaCut),
  fTrackMinPtCut(in.fTrackMinPtCut),
  fMaxTrackPtRelativeError(in.fMaxTrackPtRelativeError),
  fMaxTrackDistanceToVertex(in.fMaxTrackDistanceToVertex),
  fCalorimeterSignalLimitPt(in.fCalorimeterSignalLimitPt),
  fHighPtEtFraction(in.fHighPtEtFraction),
  fChi2QualityCut(in.fChi2QualityCut),
  fMinimumTrackHits(in.fMinimumTrackHits),
  fSubeventCut(in.fSubeventCut),
  fMcCorrelationType(in.fMcCorrelationType),
  fAsymmetryBinType(in.fAsymmetryBinType)
{
  // Copy constructor
  fTrackReader[DijetHistograms::kSameEvent] = in.fTrackReader[DijetHistograms::kSameEvent];
  fTrackReader[DijetHistograms::kMixedEvent] = in.fTrackReader[DijetHistograms::kMixedEvent];
  
  // Copy the mixing pool
  for(Int_t iVz = 0; iVz < kMaxMixingVzBins; iVz++){
    for(Int_t iHiBin = 0; iHiBin < kMaxMixingHiBins; iHiBin++){
      fMixingPool[iVz][iHiBin] = in.fMixingPool[iVz][iHiBin];
    }
  }
}

/*
 * Assingment operator
 */
DijetAnalyzer& DijetAnalyzer::operator=(const DijetAnalyzer& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fJetReader = in.fJetReader;
  fTrackReader[DijetHistograms::kSameEvent] = in.fTrackReader[DijetHistograms::kSameEvent];
  fTrackReader[DijetHistograms::kMixedEvent] = in.fTrackReader[DijetHistograms::kMixedEvent];
  fFileNames = in.fFileNames;
  fCard = in.fCard;
  fHistograms = in.fHistograms;
  fTrackCorrection = in.fTrackCorrection;
  fJffCorrection = in.fJffCorrection;
  fVzWeightFunction = in.fVzWeightFunction;
  fCentralityWeightFunction = in.fCentralityWeightFunction;
  fPtWeightFunction = in.fPtWeightFunction;
  fDijetWeightFunction = in.fDijetWeightFunction;
  fSmearingFunction = in.fSmearingFunction;
  fRng = in.fRng;
  fDataType = in.fDataType;
  fForestType = in.fForestType;
  fReadMode = in.fReadMode;
  fJetType = in.fJetType;
  fMatchJets = in.fMatchJets;
  fMatchDijet = in.fMatchDijet;
  fMatchLeadingJet = in.fMatchLeadingJet;
  fDebugLevel = in.fDebugLevel;
  fLocalRun = in.fLocalRun;
  fVzWeight = in.fVzWeight;
  fCentralityWeight = in.fCentralityWeight;
  fPtHatWeight = in.fPtHatWeight;
  fTotalEventWeight = in.fTotalEventWeight;
  fnEventsInMixingFile = in.fnEventsInMixingFile;
  fnMixedEventsPerDijet = in.fnMixedEventsPerDijet;
  fMixingStartIndex = in.fMixingStartIndex;
  fMixingPoolDepth = in.fMixingPoolDepth;
  fMixingVzBinWidth = in.fMixingVzBinWidth;
  fMixingHiBinWidth = in.fMixingHiBinWidth;
  fRunningMixingIndex = in.fRunningMixingIndex;
  fMixingVzTolerance = in.fMixingVzTolerance;
  fMixedEventVz = in.fMixedEventVz;
  fMixedEventHiBin = in.fMixedEventHiBin;
  fJetAxis = in.fJetAxis;
  fVzCut = in.fVzCut;
  fMinimumPtHat = in.fMinimumPtHat;
  fMaximumPtHat = in.fMaximumPtHat;
  fJetEtaCut = in.fJetEtaCut;
  fJetSearchEtaCut = in.fJetSearchEtaCut;
  fJetMaximumPtCut = in.fJetMaximumPtCut;
  fLeadingJetMinPtCut = in.fLeadingJetMinPtCut;
  fSubleadingJetMinPtCut = in.fSubleadingJetMinPtCut;
  fDeltaPhiCut = in.fDeltaPhiCut;
  fMinimumMaxTrackPtFraction = in.fMinimumMaxTrackPtFraction;
  fMaximumMaxTrackPtFraction = in.fMaximumMaxTrackPtFraction;
  fJetUncertaintyMode = in.fJetUncertaintyMode;
  fTrackEtaCut = in.fTrackEtaCut;
  fTrackMinPtCut = in.fTrackMinPtCut;
  fMaxTrackPtRelativeError = in.fMaxTrackPtRelativeError;
  fMaxTrackDistanceToVertex = in.fMaxTrackDistanceToVertex;
  fCalorimeterSignalLimitPt = in.fCalorimeterSignalLimitPt;
  fHighPtEtFraction = in.fHighPtEtFraction;
  fChi2QualityCut = in.fChi2QualityCut;
  fMinimumTrackHits = in.fMinimumTrackHits;
  fSubeventCut = in.fSubeventCut;
  fMcCorrelationType = in.fMcCorrelationType;
  fAsymmetryBinType = in.fAsymmetryBinType;
  
  // Copy the mixing pool
  for(Int_t iVz = 0; iVz < kMaxMixingVzBins; iVz++){
    for(Int_t iHiBin = 0; iHiBin < kMaxMixingHiBins; iHiBin++){
      fMixingPool[iVz][iHiBin] = in.fMixingPool[iVz][iHiBin];
    }
  }
  
  return *this;
}

/*
 * Destructor
 */
DijetAnalyzer::~DijetAnalyzer(){
  // destructor
  delete fHistograms;
  delete fTrackCorrection;
  delete fJffCorrection;
  if(fVzWeightFunction) delete fVzWeightFunction;
  if(fTrackEfficiencyCorrector2018) delete fTrackEfficiencyCorrector2018;
  if(fJetCorrector2018) delete fJetCorrector2018;
  if(fJetUncertainty2018) delete fJetUncertainty2018;
  if(fCentralityWeightFunction) delete fCentralityWeightFunction;
  if(fPtWeightFunction) delete fPtWeightFunction;
  if(fDijetWeightFunction) delete fDijetWeightFunction;
  if(fSmearingFunction) delete fSmearingFunction;
  if(fRng) delete fRng;
  if(fJetReader) delete fJetReader;
  if(fTrackReader[DijetHistograms::kSameEvent] && (fMcCorrelationType == kGenReco || fMcCorrelationType == kRecoGen)) delete fTrackReader[DijetHistograms::kSameEvent];
  if(fTrackReader[DijetHistograms::kMixedEvent]) delete fTrackReader[DijetHistograms::kMixedEvent];
}

/*
 * Get the number of particle flow candidates within a range of 0.4 of the given eta-phi angle
 *
 *  Arguments:
 *   Double_t jetPhi = Phi angle of the jet
 *   Double_t jetEta = Eta angle of the jet
 *
 *  return: A tuple with the following information:
 *            0: Number of particle flow candidates within jet cone of R = 0.4 from jet axis
 *            1: Leading particle flow candidate phi
 *            2: Leading particle flow candidate eta
 */
std::tuple<Int_t,Double_t,Double_t> DijetAnalyzer::GetNParticleFlowCandidatesInJet(const Double_t jetPhi, const Double_t jetEta){
  
  // Variables for particle flow candidate properties
  Double_t leadingParticleFlowCandidatePt = 0;
  Double_t leadingParticleFlowCandidatePhi = 0;
  Double_t leadingParticleFlowCandidateEta = 0;
  Double_t particleFlowCandidatePt = 0;
  Double_t particleFlowCandidatePhi = 0;
  Double_t particleFlowCandidateEta = 0;
  Double_t deltaPhi = 0;
  
  // Start counting from zero
  Int_t nParticleFlowCandidatesInThisJet = 0;
  Double_t distanceToThisJet = 0;
  
  // Loop over all particle flow candidates
  for(Int_t iParticleFlowCandidate = 0; iParticleFlowCandidate < fJetReader->GetNParticleFlowCandidates(); iParticleFlowCandidate++){
    if(fJetReader->GetParticleFlowCandidateId(iParticleFlowCandidate) != 1) continue; // Require ID 1 for candidates
    particleFlowCandidatePt = fJetReader->GetParticleFlowCandidatePt(iParticleFlowCandidate);
    if(particleFlowCandidatePt < 2) continue; // Require minimum pT of 2 GeV for candidates
    
    // Specific cuts only for generator level tracks
    if(fMcCorrelationType == kGenReco || fMcCorrelationType == kGenGen){
      if(fJetReader->GetTrackCharge(iParticleFlowCandidate) == 0) continue;                 // Require that the track is charged
      if(!PassSubeventCut(fJetReader->GetTrackSubevent(iParticleFlowCandidate))) continue;  // Require desired subevent
      if(fJetReader->GetTrackMCStatus(iParticleFlowCandidate) != 1) continue;               // Require final state particles
    }
    
    particleFlowCandidatePhi = fJetReader->GetParticleFlowCandidatePhi(iParticleFlowCandidate);
    particleFlowCandidateEta = fJetReader->GetParticleFlowCandidateEta(iParticleFlowCandidate);
    deltaPhi = jetPhi-particleFlowCandidatePhi;
    
    // Transform deltaPhi to interval [-pi,pi]
    while(deltaPhi > (TMath::Pi())){deltaPhi += -2*TMath::Pi();}
    while(deltaPhi < (-1.0*TMath::Pi())){deltaPhi += 2*TMath::Pi();}
    
    distanceToThisJet = TMath::Sqrt(TMath::Power(deltaPhi,2)+TMath::Power(jetEta-particleFlowCandidateEta,2));
    if(distanceToThisJet > 0.4) continue;  // Require the particle to be inside the jet cone of R = 0.4
    nParticleFlowCandidatesInThisJet++;
    
    // Update the information for the leading particle flow candidate in jet
    if(particleFlowCandidatePt > leadingParticleFlowCandidatePt){
      leadingParticleFlowCandidatePt = particleFlowCandidatePt;
      leadingParticleFlowCandidatePhi = particleFlowCandidatePhi;
      leadingParticleFlowCandidateEta = particleFlowCandidateEta;
    }
    
  }
  
  // For generator level jets, set nParticleFlowCandidates to -1 to show that no correction is needed
  if(fMcCorrelationType == kGenReco || fMcCorrelationType == kGenGen) nParticleFlowCandidatesInThisJet = -1;
  
  // Return a tuple of the number of particle flow candidates and information about leading particle flow candidate
  return std::make_tuple(nParticleFlowCandidatesInThisJet,leadingParticleFlowCandidatePhi,leadingParticleFlowCandidateEta);
}

/*
 * Main analysis loop
 */
void DijetAnalyzer::RunAnalysis(){
  
  //****************************************
  //        Event selection cuts
  //****************************************
  
  fMinimumPtHat = fCard->Get("LowPtHatCut");  // Minimum accepted pT hat value
  fMaximumPtHat = fCard->Get("HighPtHatCut"); // Maximum accepted pT hat value
  
  //****************************************
  //        Dijet selection cuts
  //****************************************
  
  fJetEtaCut = fCard->Get("JetEtaCut");           // Eta cut around midrapidity
  fJetSearchEtaCut = fCard->Get("SearchEtaCut");  // Eta cut when searching for a dijet
  fJetMaximumPtCut = fCard->Get("MaxPtCut");      // Maximum pT accepted for leading jet (and tracks)
  fLeadingJetMinPtCut = fCard->Get("MinPtCut");   // Minimum pT cut for leading jet
  fSubleadingJetMinPtCut = fCard->Get("SubleadingPtCut"); // Minimum pT cut for subleading jet
  fDeltaPhiCut = fCard->Get("DeltaPhiCut");       // DeltaPhi cut for the dijet system
  fMinimumMaxTrackPtFraction = fCard->Get("MinMaxTrackPtFraction");  // Cut for jets consisting only from soft particles
  fMaximumMaxTrackPtFraction = fCard->Get("MaxMaxTrackPtFraction");  // Cut for jets consisting only from one high pT particle
  fJetUncertaintyMode = fCard->Get("JetUncertainty");  // Select whether to use nominal jet pT or vary it within uncertainties
  
  //****************************************
  //        Track selection cuts
  //****************************************
  
  fTrackEtaCut = fCard->Get("TrackEtaCut");     // Eta cut around midrapidity
  fTrackMinPtCut = fCard->Get("MinTrackPtCut"); // Minimum pT cut
  fMaxTrackPtRelativeError = fCard->Get("MaxTrackPtRelativeError");   // Maximum relative error for pT
  fMaxTrackDistanceToVertex = fCard->Get("VertexMaxDistance");        // Maximum distance to primary vetrex
  fCalorimeterSignalLimitPt = fCard->Get("CalorimeterSignalLimitPt"); // Require signal in calorimeters for track above this pT
  fHighPtEtFraction = fCard->Get("HighPtEtFraction"); // For high pT tracks, minimum required Et as a fraction of track pT
  fChi2QualityCut = fCard->Get("Chi2QualityCut");     // Quality cut for track reconstruction
  fMinimumTrackHits = fCard->Get("MinimumTrackHits"); // Quality cut for track hits
  
  fSubeventCut = fCard->Get("SubeventCut");     // Required subevent type
  
  //****************************************
  //    Correlation type for Monte Carlo
  //****************************************
  if(fDataType == ForestReader::kPp || fDataType == ForestReader::kPbPb) {
    fMcCorrelationType = -100;
  } else {
    fMcCorrelationType = fCard->Get("McCorrelationType");         // Correlation type for Monte Carlo
  }
  Int_t mixingFileIndex = fCard->Get("MixingFileIndex");        // Select the used mixing file for PbPb MC

  //****************************************
  //       Jet selection and matching
  //****************************************
  fForestType = fCard->Get("ForestType");
  fJetType = fCard->Get("JetType");
  fMatchJets = (fCard->Get("MatchJets") >= 1);
  fMatchDijet = (fCard->Get("MatchJets") == 2 || fCard->Get("MatchJets") == 4);
  fMatchLeadingJet = (fCard->Get("MatchJets") == 3 || fCard->Get("MatchJets") == 5);
  Bool_t reverseMatchVeto = (fCard->Get("MatchJets") == 4 || fCard->Get("MatchJets") == 5);
  Bool_t findMatchedDijet = fCard->Get("MatchJets") == 6;
  
  //*************************************************************
  //    Turn off certain track cuts for generated tracks and pp
  //*************************************************************
  
  if(fMcCorrelationType == kGenGen || fMcCorrelationType == kRecoGen || fReadMode == 2017){
    fCalorimeterSignalLimitPt = 10000;
    fChi2QualityCut = 10000;
    fMinimumTrackHits = 0;
  }
  
  //****************************************
  //            All cuts set!
  //****************************************
  
  //************************************************
  //  Define variables needed in the analysis loop
  //************************************************
  
  // Input files and forest readers for analysis
  TFile *inputFile;
  TFile *copyInputFile; // If we read forest for tracks and jets with different readers, we need to give different file pointers to them
  TFile *mixedEventFile;
  
  // Event variables
  Int_t nEvents = 0;                // Number of events
  Bool_t dijetFound = false;        // Is there a dijet in the event?
  Bool_t twoJetsFound = false;      // Are there two jets in the event?
  Bool_t matchVeto = false;         // Veto the found dijet because it did not match with Gen or Reco
  Double_t vz = 0;                  // Vertex z-position
  Double_t centrality = 0;          // Event centrality
  Int_t hiBin = 0;                  // CMS hiBin (centrality * 2)
  Double_t ptHat = 0;               // pT hat for MC events
  
  // Combining bools to make the code more readable
  Bool_t useDifferentReaderFotJetsAndTracks = (fMcCorrelationType == kRecoGen || fMcCorrelationType == kGenReco); // Use different forest reader for jets and tracks
  
  // Event mixing information
  Bool_t mixEvents = (fCard->Get("DoEventMixing") == 1);  // Do or do not do event mixing
  Bool_t mixWithPool = (fCard->Get("MixWithPool") == 1);  // Select whether to use mixing pool or directly vz and centrality difference in mixing file
  Bool_t onlyMix = (fCard->Get("OnlyMix") == 1); // Only fill mixed event histograms. Option to generate more mixing events and merge them with previous runs including same event histograms without duplicating same event statistics.
  Bool_t doEventPlane = (fCard->Get("IncludeEventPlane") == 1); // Tell if the event plane branches are included in the data files
  std::vector<TString> mixingFiles;  // List of mixing files
  
  // Variables for jets
  Double_t dijetAsymmetry = -99;    // Dijet momuntem balance used for binning
  Double_t dijetXj = -99;           // Dijet momentum balance xj
  Double_t dijetMatchedXj = -99;    // Matched dijet momentum balance
  Double_t leadingJetPt = 0;        // Leading jet pT
  Double_t leadingJetPhi = 0;       // Leading jet phi
  Double_t leadingJetEta = 0;       // Leading jet eta
  Double_t leadingJetFlavor = 0;    // Flavor of the leading jet
  Double_t subleadingJetPt = 0;     // Subleading jet pT
  Double_t subleadingJetPhi = 0;    // Subleading jet phi
  Double_t subleadingJetEta = 0;    // Subleading jet eta
  Double_t subleadingJetFlavor = 0; // Flavor of the subleading jet
  Double_t thirdJetPt = 0;          // Third highest jet pT
  Double_t leadingParticleFlowCandidatePhi = 0;     // Leading particle lofw candidate phi
  Double_t leadingParticleFlowCandidateEta = 0;     // Leading particle flow candidate eta
  Double_t subleadingParticleFlowCandidatePhi = 0;  // Subleading particle flow candidate phi
  Double_t subleadingParticleFlowCandidateEta = 0;  // Subleading particle flow candidate eta
  Double_t swapJetPt = 0;           // Swapping helper variable
  Double_t swapJetPhi = 0;          // Swapping helper variable
  Double_t swapJetEta = 0;          // Swapping helper variable
  Int_t secondHighestIndex = -1;    // Index of the subleading jet in the event
  Int_t highestIndex = -1;          // Index of the leading jet in the event
  Int_t swapIndex = -1;             // Swapping helper variable
  Double_t jetPt = 0;               // pT of the i:th jet in the event
  Double_t jetPtCorrected = 0;      // Jet pT corrected with the JFF correction
  Double_t jetPhi = 0;              // phi of the i:th jet in the event
  Double_t jetEta = 0;              // eta of the i:th jet in the event
  Int_t jetFlavor = 0;              // Flavor of the jet. 0 = Quark jet. 1 = Gluon jet.
  Double_t highestAnyPt = 0;        // Highest pT filled for all jets
  Double_t highestPhi = 0;          // phi of any leading jet
  Double_t highestEta = 0;          // eta of any leading jet
  Double_t dphi = 0;                // deltaPhi for the considered jets
  Double_t leadingJetInfo[4] = {0};       // Array for leading jet pT, phi and eta
  Double_t subleadingJetInfo[4] = {0};    // Array for subleading jet pT, phi and eta
  Double_t inclusiveJetInfo[60][4] = {{0}}; // Array for jet pT, phi and eta for all the jets in the event
  Int_t nJetsInThisEvent = 0;             // Number of jets in the event
  Double_t matchedLeadingJetPt = 1;       // Leading matched jet pT
  Double_t matchedLeadingJetPhi = 0;      // Leading matched jet phi
  Double_t matchedLeadingJetEta = 0;      // Leading matched jet eta
  Double_t matchedSubleadingJetPt = 0;    // Subleading matched jet pT
  Double_t matchedSubleadingJetPhi = 0;   // Subleading matched jet phi
  Double_t matchedSubleadingJetEta = 0;   // Subleading matched jet eta
  Double_t matchedDeltaPhi = 0;           // DeltaPhi between matched leading and subleading jets
  Double_t jetPtWeight = 1;               // Weighting for jet pT
  Double_t triggerEfficiencyWeight = 1;   // Weight for trigger efficiency in pT
  
  // Variables for smearing study
  Double_t smearingFactor = 0;       // Larger of the JEC uncertainties
//  Double_t jetPtSmeared = 0;          // Smeared jet pT
//  Double_t jetPtErrorUp = 0;          // Uncertainty to be added to the jet pT
//  Double_t jetPtErrorDown = 0;        // Uncertainty to be subtracted from the jet pT
//  Double_t highestJetPtSmeared = 0;   // Smeared pT of the highest jet
//  Double_t highestJetPtErrorUp = 0;   // Uncertainty to be added to the jet with highest pT
//  Double_t highestJetPtErrorDown = 0; // Uncertainty to be subtracted from the jet with highest pT
  
  // Variables for jet matching and closure
  Int_t unmatchedCounter = 0;       // Number of jets that fail the matching
  Int_t matchedCounter = 0;         // Number of jets that are matched
  Int_t dijetCounter = 0;           // Counter for dijets before matching
  Int_t matchedDijetCounter = 0;    // Counter for dijets after matching
  Int_t jetSwapCounter = 0;         // Counter for leading and subleading jets that are swapped between reco and gen
  Int_t nonSensicalPartonIndex = 0; // Parton index is -999 even though jets are matched
  Int_t partonFlavor = -999;        // Code for parton flavor in Monte Carlo
  Int_t highThirdJet = 0;           // TODO: Debuggery
  
  // Variables for tracks
  Double_t fillerTrack[4];            // Track histogram filler
  Double_t trackEfficiencyCorrection; // Track efficiency correction
  Int_t nTracks;                      // Number of tracks in an event
  Double_t trackPt = 0;               // Track pT
  Double_t trackEta = 0;              // Track eta
  Double_t trackPhi = 0;              // Track phi
  
  // Variables for particle flow candidates
  Int_t nParticleFlowCandidatesInThisJet = 0;  // Number of particle flow candidates in the current jet
  
  // File name helper variables
  TString currentFile;
  TString currentMixedEventFile;
  
  // Fillers for THnSparses
  const Int_t nFillJet = 5;         // 5 is nominal, 8 used for smearing study
  const Int_t nFillLeadingJet = 6;  // 6 is nominal, 9 used for smearing study
  const Int_t nFillDijet = 6;      // 6 is nominal, 10 used for xj study
  Double_t fillerJet[nFillJet];
  Double_t fillerLeadingJet[nFillLeadingJet];
  Double_t fillerDijet[nFillDijet];
  
  // For 2018 PbPb and 2017 pp data, we need to correct jet pT
  std::string correctionFileRelative[5] = {"jetEnergyCorrections/Spring18_ppRef5TeV_V4_DATA_L2Relative_AK4PF.txt", "jetEnergyCorrections/Autumn18_HI_V6_DATA_L2Relative_AK4PF.txt", "jetEnergyCorrections/Spring18_ppRef5TeV_V4_MC_L2Relative_AK4PF.txt", "jetEnergyCorrections/Autumn18_HI_V6_MC_L2Relative_AK4PF.txt", "jetEnergyCorrections/Autumn18_HI_V4_DATA_L2Relative_AK4PF.txt"};
  std::string correctionFileResidual[5] = {"jetEnergyCorrections/Spring18_ppRef5TeV_V4_DATA_L2L3Residual_AK4PF.txt", "jetEnergyCorrections/Autumn18_HI_V6_DATA_L2L3Residual_AK4PF.txt", "CorrectionNotAppliedPF.txt", "CorrectionNotAppliedPF.txt", "jetEnergyCorrections/Autumn18_HI_V6_DATA_L2L3Residual_AK4PF.txt"};
  std::string uncertaintyFile[5] = {"jetEnergyCorrections/Spring18_ppRef5TeV_V4_DATA_Uncertainty_AK4PF.txt", "jetEnergyCorrections/Autumn18_HI_V6_DATA_Uncertainty_AK4PF.txt", "jetEnergyCorrections/Spring18_ppRef5TeV_V4_MC_Uncertainty_AK4PF.txt", "jetEnergyCorrections/Autumn18_HI_V6_MC_Uncertainty_AK4PF.txt", "jetEnergyCorrections/Autumn18_HI_V6_DATA_Uncertainty_AK4PF.txt"};
  
  // For calo jets, use the correction files for calo jets (otherwise same name, but replace PF with Calo)
  if(fJetType == 0){
    size_t pfIndex = 0;
    pfIndex = correctionFileRelative[fDataType].find("PF", pfIndex);
    correctionFileRelative[fDataType].replace(pfIndex, 2, "Calo");
    pfIndex = 0;
    pfIndex = correctionFileResidual[fDataType].find("PF", pfIndex);
    correctionFileResidual[fDataType].replace(pfIndex, 2, "Calo");
    pfIndex = 0;
    pfIndex = uncertaintyFile[fDataType].find("PF", pfIndex);
    uncertaintyFile[fDataType].replace(pfIndex, 2, "Calo");
    
    correctionFileRelative[1] = "jetEnergyCorrections/75X_dataRun2_v13_L2Relative_AKPu4Calo_stabbed.txt";
    correctionFileRelative[3] = "jetEnergyCorrections/75X_dataRun2_v13_L2Relative_AKPu4Calo_stabbed.txt";
    
  }
    
  vector<string> correctionFiles;
  correctionFiles.push_back(correctionFileRelative[fDataType]);
  if(fJetType > 0 && (fDataType == ForestReader::kPbPb || fDataType == ForestReader::kPp))  correctionFiles.push_back(correctionFileResidual[fDataType]);
  
  fJetCorrector2018 = new JetCorrector(correctionFiles);
  fJetUncertainty2018 = new JetUncertainty(uncertaintyFile[fDataType]);
  
  //************************************************
  //      Find forest readers for data files
  //************************************************
  
  if(fMcCorrelationType == kGenReco || fMcCorrelationType == kGenGen){
    if(fForestType == kSkimForest) {
      fJetReader = new GeneratorLevelSkimForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane);
    } else {
      fJetReader = new GeneratorLevelForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane);
    }
  } else {
    if(fForestType == kSkimForest) {
      fJetReader = new SkimForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane);
    } else {
      fJetReader = new HighForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane);
    }
  }
  
  // Select the reader for tracks based on forest and MC correlation type
  if(fMcCorrelationType == kRecoGen && fForestType == kSkimForest){
    fTrackReader[DijetHistograms::kSameEvent] = new GeneratorLevelSkimForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane);
  } else if(fMcCorrelationType == kRecoGen){
    fTrackReader[DijetHistograms::kSameEvent] = new GeneratorLevelForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane);
  } else if (fMcCorrelationType == kGenReco && fForestType == kSkimForest){
    fTrackReader[DijetHistograms::kSameEvent] = new SkimForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane);
  } else if (fMcCorrelationType == kGenReco){
    fTrackReader[DijetHistograms::kSameEvent] = new HighForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane);
  } else {
    fTrackReader[DijetHistograms::kSameEvent] = fJetReader;
  }
  
  // If mixing events, create ForestReader for that. For PbPb and PbPbMC, the Forest in mixing file is a skim forest
  if(mixEvents){
    if(fDataType == ForestReader::kPbPb){
      if(fReadMode > 2000){
        fTrackReader[DijetHistograms::kMixedEvent] = new MixingForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane); // Reader for 2018 mixing files
      } else {
        fTrackReader[DijetHistograms::kMixedEvent] = new SkimForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane); // Reader for 2015 mixing files
      }
    } else if(fDataType == ForestReader::kPbPbMC){
      if (fMcCorrelationType == kRecoGen || fMcCorrelationType == kGenGen) { // Mixed event reader for generator tracks
        if(fReadMode > 2000){
          fTrackReader[DijetHistograms::kMixedEvent] = new GeneratorLevelMixingForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane); // Reader for 2018 syntax
        } else {
          fTrackReader[DijetHistograms::kMixedEvent] = new GeneratorLevelSkimForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane); // Reader for 2015 syntax
        }
        
      } else {
        if(fReadMode > 2000){
          fTrackReader[DijetHistograms::kMixedEvent] = new MixingForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane); // Reader for 2018 syntax
        } else {
          fTrackReader[DijetHistograms::kMixedEvent] = new SkimForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane); // Reader for 2015 syntax
        }
      }
    } else if (fMcCorrelationType == kRecoGen || fMcCorrelationType == kGenGen) { // Mixed event reader for generator tracks
      if(fForestType == kSkimForest) {
        fTrackReader[DijetHistograms::kMixedEvent] = new GeneratorLevelSkimForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane);
      } else {
        fTrackReader[DijetHistograms::kMixedEvent] = new GeneratorLevelForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane);
      }
    } else {
      if(fForestType == kSkimForest) {
        fTrackReader[DijetHistograms::kMixedEvent] = new SkimForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane);
      } else {
        fTrackReader[DijetHistograms::kMixedEvent] = new HighForestReader(fDataType,fReadMode,fJetType,fJetAxis,fMatchJets,doEventPlane);
      }
    }
    
    // If we are doing mixing, find and open the mixing file
    // PbPb data has different data file for mixing, other data sets use the regular data files for mixing
    // This is done outside the file loop in the cases where mixing is not done from the same file as data
    if(fDataType == ForestReader::kPbPb || fDataType == ForestReader::kPbPbMC){
      
      // 2018 syntax
      if(fReadMode > 2000){
        
        // Define names for mixing file lists for different datasets
        const char* fileListName[2][5] = {
          {"none", "mixingFileList/PbPbData2018_MinBiasFiles.txt", "none", "mixingFileList/PbPbMC2018_MinBiasFiles.txt", "none"}, // Crab
          {"none", "mixingFileList/mixingFilesPbPb.txt", "none", "mixingFileList/mixingFilesPbPbMC.txt", "none"}}; // Local test
        
        // Create a stream to read the input file
        std::string lineInFile;
        std::ifstream mixingFileStream(fileListName[fLocalRun][fDataType]);
        while (std::getline(mixingFileStream,lineInFile)) {
          mixingFiles.push_back(lineInFile);
        }
        fTrackReader[DijetHistograms::kMixedEvent]->ReadForestFromFileList(mixingFiles); // TODO: Testing new reader
      
      } else { // 2015 syntax
        
        if(fDataType == ForestReader::kPbPb){
          const char* pbpbMixingFileNames[] = {"root://xrootd.rcac.purdue.edu//store/user/kumarv/UICWork/Data2015_PbPb_5TeV_MinBiasSkim/Data2015_finalTrkCut_1Mevts.root","root://xrootd.rcac.purdue.edu//store/user/kumarv/UICWork/Data2015_PbPb_5TeV_MinBiasSkim/Data2015_finalTrkCut_1Mevts_copy1.root"};
          currentMixedEventFile = pbpbMixingFileNames[mixingFileIndex];
        } else {
          const char* pbpbMCMixingFileNames[] = {"root://xrootd.rcac.purdue.edu//store/user/wangx/PbPbMC_Py6H_skim_looseTrkCuts_finalJFFs_lowPtGenPartCut_CymbalTune/Pythia6Hydjet_PbPbMC_cymbalTune_mixMerged.root","root://xrootd.rcac.purdue.edu//store/user/jviinika/PbPbMC_Py6H_skim_looseTrkCuts_finalJFFs_lowPtGenPartCut_CymbalTune/Pythia6Hydjet_PbPbMC_cymbalTune_mixMerged.root"};
          currentMixedEventFile = pbpbMCMixingFileNames[mixingFileIndex];
        }
        
        mixedEventFile = TFile::Open(currentMixedEventFile);
        fTrackReader[DijetHistograms::kMixedEvent]->ReadForestFromFile(mixedEventFile);
        
      } // 2015 syntax if
      
      // Read the number of events in the mixing file
      fnEventsInMixingFile = fTrackReader[DijetHistograms::kMixedEvent]->GetNEvents();
      
      // Prapare for event mixing based on the chosen method
      if(mixWithPool){
        CreateMixingPool();
      } else {
        PrepareMixingVectors();
      }
      
    } // Mixing from different file than the data file
    
  } // If for mixing
  
  //************************************************
  //       Main analysis loop over all files
  //************************************************
  
  // Loop over files
  Int_t nFiles = fFileNames.size();
  for(Int_t iFile = 0; iFile < nFiles; iFile++) {
    
    //************************************************
    //              Find and open files
    //************************************************
    
    // Find the filename and open the input file
    currentFile = fFileNames.at(iFile);
    inputFile = TFile::Open(currentFile);
    if(useDifferentReaderFotJetsAndTracks) copyInputFile = TFile::Open(currentFile);

    // Mixing in the case we use the data file also for mixing (true for pp)
    if(mixEvents && fDataType != ForestReader::kPbPb && fDataType != ForestReader::kPbPbMC){
      currentMixedEventFile = fFileNames.at(iFile);
      mixedEventFile = TFile::Open(currentMixedEventFile);
      mixingFiles.push_back(currentMixedEventFile.Data());
    }
    
    // Check that the file exists
    if(!inputFile){
      cout << "Error! Could not find the file: " << currentFile.Data() << endl;
      assert(0);
    }

    // Check that the file is open
    if(!inputFile->IsOpen()){
      cout << "Error! Could not open the file: " << currentFile.Data() << endl;
      assert(0);
    }
    
    // Check that the file is not zombie
    if(inputFile->IsZombie()){
      cout << "Error! The following file is a zombie: " << currentFile.Data() << endl;
      assert(0);
    }
    
    // Check that the mixing file exists if we are reading only single file
    if(mixEvents && ((fDataType != ForestReader::kPbPb && fDataType != ForestReader::kPbPbMC) || fReadMode < 2000)){
      
      if(!mixedEventFile){
        cout << "Error! Could not find the mixing file: " << currentMixedEventFile.Data() << endl;
        assert(0);
      }
      
      // Check that the mixing file is open
      if(!mixedEventFile->IsOpen()){
        cout << "Error! Could not open the mixing file: " << currentMixedEventFile.Data() << endl;
        assert(0);
      }
      
      // Check that the file is not zombie
      if(mixedEventFile->IsZombie()){
        cout << "Error! The following mixing file is a zombie: " << currentMixedEventFile.Data() << endl;
        assert(0);
      }
    }

    // Print the used files
    if(fDebugLevel > 0) cout << "Reading from file: " << currentFile.Data() << endl;
    if(fDebugLevel > 0 && mixEvents) {
      cout << "Mixing files: " << currentMixedEventFile.Data() << endl;
      for(std::vector<TString>::iterator mixIterator = mixingFiles.begin(); mixIterator != mixingFiles.end(); mixIterator++){
        cout << *mixIterator << endl;
      }
      cout << endl;
    }
    
    //************************************************
    //        Read forest and prepare mixing
    //************************************************
    
    // If file is good, read the forest from the file
    fJetReader->ReadForestFromFile(inputFile);  // There might be a memory leak in handling the forest...
    if(useDifferentReaderFotJetsAndTracks) fTrackReader[DijetHistograms::kSameEvent]->ReadForestFromFile(copyInputFile); // If we mix reco and gen, the reader for jets and tracks is different
    nEvents = fJetReader->GetNEvents();
    
    // Read also the forest for event mixing and prepare mixing pool
    if(mixEvents && fDataType != ForestReader::kPbPb && fDataType != ForestReader::kPbPbMC){
      fTrackReader[DijetHistograms::kMixedEvent]->ReadForestFromFile(mixedEventFile);
      fnEventsInMixingFile = fTrackReader[DijetHistograms::kMixedEvent]->GetNEvents();
      
      // Prapare for event mixing based on the chosen method
      if(mixWithPool){
        CreateMixingPool();
      } else {
        PrepareMixingVectors();
      }
    }

    //************************************************
    //         Main event loop for each file
    //************************************************
    
    for(Int_t iEvent = 0; iEvent < nEvents; iEvent++){
      
      //************************************************
      //         Read basic event information
      //************************************************
      
      // Print to console how the analysis is progressing
      if(fDebugLevel > 1 && iEvent % 1000 == 0) cout << "Analyzing event " << iEvent << endl;
      
      // Read the event to memory
      fJetReader->GetEvent(iEvent);
      
      // If track reader is not the same as jet reader, read the event to memory in trackReader
      if(useDifferentReaderFotJetsAndTracks) fTrackReader[DijetHistograms::kSameEvent]->GetEvent(iEvent);

      // Get vz, centrality and pT hat information
      vz = fJetReader->GetVz();
      centrality = fJetReader->GetCentrality();
      hiBin = fJetReader->GetHiBin();
      ptHat = fJetReader->GetPtHat();
      
      // We need to apply pT hat cuts before getting pT hat weight. There might be rare events above the upper
      // limit from which the weights are calculated, which could cause the code to crash.
      if(ptHat < fMinimumPtHat || ptHat >= fMaximumPtHat) continue;
      
      // Get the weighting for the event
      fVzWeight = GetVzWeight(vz);
      fCentralityWeight = GetCentralityWeight(hiBin);
      
      // Event weight is different for 2015 and 2018 MC
      if(fReadMode < 2000){
        fPtHatWeight = GetPtHatWeight(ptHat); // 2015 MC
      } else {
        fPtHatWeight = fJetReader->GetEventWeight(); // 2018 MC
      }
      
      fTotalEventWeight = fVzWeight*fCentralityWeight*fPtHatWeight;
      
      fHistograms->fhEvents->Fill(DijetHistograms::kAll);          // All the events looped over
      
      //  ============================================
      //  ===== Apply all the event quality cuts =====
      //  ============================================
      
      if(!PassEventCuts(fJetReader,fFillEventInformation,DijetHistograms::kSameEvent)) continue;
      
      // Fill the event information histograms for the events that pass the event cuts
      if(fFillEventInformation){
        fHistograms->fhVertexZ->Fill(vz);                            // z vertex distribution from all events
        fHistograms->fhVertexZWeighted->Fill(vz,fVzWeight);          // z-vertex distribution weighted with the weight function
        fHistograms->fhCentrality->Fill(centrality);                 // Centrality filled from all events
        fHistograms->fhCentralityWeighted->Fill(centrality,fCentralityWeight); // Centrality weighted with the centrality weighting function
        fHistograms->fhPtHat->Fill(ptHat);                           // pT hat histogram
        fHistograms->fhPtHatWeighted->Fill(ptHat,fPtHatWeight);      // pT het histogram weighted with corresponding cross section and event number
      }
      
      // Q-vector cut for event plane. Used in a MC study to try to match flow with data.
      if(doEventPlane){
        
        // Variables for event plane
        double eventPlaneQ = 0;            // Magnitude of the event plane Q-vector
        double eventPlaneMultiplicity = 0; // Particle multiplicity in the event plane
        //int iEventPlane = 9; // For event planes, see the big comment in ForestReader.h
        double eventPlaneQx = 0;
        double eventPlaneQy = 0;
        
        // Manual calculation for Q-vector
        // Loop over all track in the event
        nTracks = fTrackReader[DijetHistograms::kSameEvent]->GetNTracks();
        for(Int_t iTrack = 0; iTrack < nTracks; iTrack++){
          
          // Check that all the track cuts are passed
          //if(!PassTrackCuts(iTrack,fHistograms->fhTrackCutsInclusive,DijetHistograms::kSameEvent)) continue;
          
          // Get the track information
          trackPt = fTrackReader[DijetHistograms::kSameEvent]->GetTrackPt(iTrack);
          trackEta = fTrackReader[DijetHistograms::kSameEvent]->GetTrackEta(iTrack);
          trackPhi = fTrackReader[DijetHistograms::kSameEvent]->GetTrackPhi(iTrack);
          //trackEfficiencyCorrection = GetTrackEfficiencyCorrection(DijetHistograms::kSameEvent,iTrack);
          
          if(TMath::Abs(trackEta) > 0.75) continue;
          if(fTrackReader[DijetHistograms::kSameEvent]->GetTrackSubevent(iTrack) == 0) continue;
          if(trackPt > 3) continue;
                    
          eventPlaneQx += TMath::Cos(2*(trackPhi));
          eventPlaneQy += TMath::Sin(2*(trackPhi));
          eventPlaneMultiplicity += 1;
          
        }
        
        if(eventPlaneMultiplicity == 0) eventPlaneMultiplicity += 1;
        eventPlaneQ = TMath::Sqrt(eventPlaneQx*eventPlaneQx + eventPlaneQy*eventPlaneQy);
        
        //eventPlaneQ = fJetReader->GetEventPlaneQ(iEventPlane);  // 8 is second order event plane from both sides of HF
        //eventPlaneMultiplicity = fJetReader->GetEventPlaneMultiplicity(iEventPlane);
        //eventPlaneQ /= TMath::Sqrt(eventPlaneMultiplicity);
                
        if(eventPlaneQ < 3.333) continue;
      }
      
      
      // ======================================
      // ===== Event quality cuts applied =====
      // ======================================
      
      // Reset the variables used in dijet finding
      twoJetsFound = false;
      dijetFound = false;
      highestIndex = -1;
      secondHighestIndex = -1;
      leadingJetPt = 0;
      subleadingJetPt = 0;
      thirdJetPt = 0;
      highestAnyPt = 0;
      highestPhi = 0;
      highestEta = 0;
      nJetsInThisEvent = 0;
      
//      // Extra variables for smearing study
//      highestJetPtSmeared = 0;
//      highestJetPtErrorUp = 0;
//      highestJetPtErrorDown = 0;
      
      //************************************************
      //    Loop over all jets and find leading jet
      //************************************************
      
      // Search for leading jet and fill histograms for all jets within the eta range
      for(Int_t jetIndex = 0; jetIndex < fJetReader->GetNJets(); jetIndex++) {
        jetPt = fJetReader->GetJetPt(jetIndex);
        jetPhi = fJetReader->GetJetPhi(jetIndex);
        jetEta = fJetReader->GetJetEta(jetIndex);
        jetFlavor = 0;
        
        // For data, instead of jet flavor, mark positive vz with 1 and negative with 0
        // This is used in one of the systematic checks for long range correlations
        if((fDataType == ForestReader::kPp || fDataType == ForestReader::kPbPb) && vz > 0) jetFlavor = 1;
        
        //  ========================================
        //  ======== Apply jet quality cuts ========
        //  ========================================
        
        if(TMath::Abs(jetEta) >= fJetSearchEtaCut) continue; // Cut for search eta range
        if(fMinimumMaxTrackPtFraction >= fJetReader->GetJetMaxTrackPt(jetIndex)/fJetReader->GetJetRawPt(jetIndex)) continue; // Cut for jets with only very low pT particles
        if(fMaximumMaxTrackPtFraction <= fJetReader->GetJetMaxTrackPt(jetIndex)/fJetReader->GetJetRawPt(jetIndex)) continue; // Cut for jets where all the pT is taken by one track
        
        // Jet matching between reconstructed and generator level jets
        if(fMatchJets && !fJetReader->HasMatchingJet(jetIndex)) {
          unmatchedCounter++;
          continue;
        }
        
        // Require also reference parton flavor to be quark [-6,-1] U [1,6] or gluon (21)
        // We need to match gen jets to reco to get the parton flavor, but for reco jets it is always available in the forest
        // Here should implement an option if only quark and gluon tagged jets should be allowed in final results!
        if(fMatchJets){
          
          if(fMatchJets) matchedCounter++; // For debugging purposes, count the number of matched jets
          jetFlavor = 0;    // Jet flavor. 0 = Quark jet.
          
          partonFlavor = fJetReader->GetPartonFlavor(jetIndex);
          if(partonFlavor == -999) nonSensicalPartonIndex++;
          if(partonFlavor < -6 || partonFlavor > 21 || (partonFlavor > 6 && partonFlavor < 21) || partonFlavor == 0) continue;
          if(TMath::Abs(partonFlavor) == 21) jetFlavor = 1; // 1 = Gluon jet
          
        }
        
        //  ========================================
        //  ======= Jet quality cuts applied =======
        //  ========================================
        
        // For 2018 data: do a correction for the jet pT
        fJetCorrector2018->SetJetPT(jetPt);
        fJetCorrector2018->SetJetEta(jetEta);
        fJetCorrector2018->SetJetPhi(jetPhi);
        
        fJetUncertainty2018->SetJetPT(jetPt);
        fJetUncertainty2018->SetJetEta(jetEta);
        fJetUncertainty2018->SetJetPhi(jetPhi);
        
        jetPtCorrected = fJetCorrector2018->GetCorrectedPT();
        
//        // Extra code for smearing study
//
//        // Add random smearing of 20 % to the jet pT
//        jetPtSmeared = jetPtCorrected*fRng->Gaus(1,0.2);     // Smearing for 20 % of the jet pT
//
//        // For the uncertainties, calculate the relative uncertainty
//        jetPtErrorUp = fJetUncertainty2018->GetUncertainty().second;
//        jetPtErrorDown = fJetUncertainty2018->GetUncertainty().first;
//
//        // Extra code for smearing study
        
        // Only do the correction for 2018 data and reconstructed Monte Carlo
        if(fReadMode > 2000 && !(fMcCorrelationType == kGenGen || fMcCorrelationType == kGenReco)) {
          jetPt = jetPtCorrected;
          
          // If we are making runs using variation of jet pT within uncertainties, modify the jet pT here
          if(fJetUncertaintyMode == 1) jetPt = jetPt * (1 - fJetUncertainty2018->GetUncertainty().first);
          if(fJetUncertaintyMode == 2) jetPt = jetPt * (1 + fJetUncertainty2018->GetUncertainty().second);
          
          // If we are using smearing scenario, modify the jet pT using gaussian smearing
          if(fJetUncertaintyMode == 3){
            smearingFactor = GetSmearingFactor(jetPt, centrality);
            jetPt = jetPt * fRng->Gaus(1,smearingFactor);
          }
          
        }
        
        // Find leading, subleading and third highest jet pT
        if(jetPt > leadingJetPt){
          thirdJetPt = subleadingJetPt;
          subleadingJetPt = leadingJetPt;
          secondHighestIndex = highestIndex;
          leadingJetPt = jetPt;
          highestIndex = jetIndex;
          leadingJetFlavor = jetFlavor;
        } else if(jetPt > subleadingJetPt){
          thirdJetPt = subleadingJetPt;
          subleadingJetPt = jetPt;
          secondHighestIndex = jetIndex;
        } else if (jetPt > thirdJetPt){
          thirdJetPt = jetPt;
        }
        
        // For jets within the specified eta range, collect any jet histograms and inclusive jet-track correlations
        if(TMath::Abs(jetEta) < fJetEtaCut){
          
          //************************************************
          //       Fill histograms for jet pT closure
          //************************************************
          
          // Only fill is matching is enabled and histograms are selected for filling
          if(fFillJetPtClosure && fMatchJets) FillJetPtClosureHistograms(jetIndex, DijetHistograms::kInclusiveClosure);
          
          //************************************************
          //         Fill histograms for all jets
          //************************************************
          
          // Only fill the any jet histogram if selected
          if(fFillJetHistograms){
            
            // Apply JFF correction for jet pT only if we are using calorimeter jets in 2015 data
            if(fJetType == 0 && fReadMode < 2000){
              std::tie(nParticleFlowCandidatesInThisJet,leadingParticleFlowCandidatePhi,leadingParticleFlowCandidateEta) = GetNParticleFlowCandidatesInJet(jetPhi,jetEta);
              jetPtCorrected = fJffCorrection->GetCorrection(nParticleFlowCandidatesInThisJet,hiBin,jetPt,jetEta);
            } else {
              jetPtCorrected = jetPt;
            }
            
            // Find the pT weight for the jet
            jetPtWeight = GetJetPtWeight(jetPtCorrected);
            
            // Fill the axes in correct order
            fillerJet[0] = jetPtCorrected;          // Axis 0 = any jet pT
            fillerJet[1] = jetPhi;                  // Axis 1 = any jet phi
            fillerJet[2] = jetEta;                  // Axis 2 = any jet eta
            fillerJet[3] = centrality;              // Axis 3 = centrality
            fillerJet[4] = jetFlavor;               // Axis 4 = flavor of the jet
            
//            // Extra axes filled for smearing study
//            fillerJet[5] = jetPtSmeared;            // Axis 5 = Smeared jet pT
//            fillerJet[6] = jetPtErrorUp;            // Axis 6 = Relative uncertainty added to the jet pT
//            fillerJet[7] = jetPtErrorDown;          // Axis 7 = Relative uncertainty subtracted from the jet pT
//            // Extra axes filled for smearing study
            
            fHistograms->fhAnyJet->Fill(fillerJet,fTotalEventWeight*jetPtWeight); // Fill the data point to histogram
            
            // Remember the hishest pT filled to any jet histograms
            if(jetPtCorrected > highestAnyPt){
              highestAnyPt = jetPtCorrected;
              highestPhi = jetPhi;
              highestEta = jetEta;
              
//              // Extra axes filled for the smearing study
//              highestJetPtSmeared = jetPtSmeared;
//              highestJetPtErrorUp = jetPtErrorUp;
//              highestJetPtErrorDown = jetPtErrorDown;
            }
            
          } // Check if we want to fill any jet histograms
          
          //************************************************
          //   Do jet-track correlation for inclusive jets
          //************************************************
          
          // If we are filling the correlation histograms and jets pass the pT cuts, do inclusive jet-track correlations
          if(fFillInclusiveJetTrackCorrelation){
            
            // Apply the JFF correction for inclusive jet pT only if we are using calorimeter jets
            // If we are using leading particle flow candidate axis, we need to determine the direction
            // of the leading particle flow candidate whether we do the correction or not.
            if((fJetType == 0 && fReadMode < 2000) || fJetAxis == 1){
            std::tie(nParticleFlowCandidatesInThisJet,leadingParticleFlowCandidatePhi,leadingParticleFlowCandidateEta) = GetNParticleFlowCandidatesInJet(jetPhi,jetEta);
            }
            if(fJetType == 0 && fReadMode < 2000){
              jetPtCorrected = fJffCorrection->GetCorrection(nParticleFlowCandidatesInThisJet,hiBin,jetPt,jetEta);
            } else {
              jetPtCorrected = jetPt;
            }
            
            // Check that the inclusive jet passes the pT cuts for the jet
            if((jetPtCorrected > fLeadingJetMinPtCut) && (jetPtCorrected < fJetMaximumPtCut)){
              
              // Fill the array with jet information
              inclusiveJetInfo[nJetsInThisEvent][0] = jetPtCorrected;
              inclusiveJetInfo[nJetsInThisEvent][1] = jetPhi;
              inclusiveJetInfo[nJetsInThisEvent][2] = jetEta;
              inclusiveJetInfo[nJetsInThisEvent][3] = jetFlavor;
              
              // If we are using leading particle flow candidate as a probe for jet axis, change the jet info
              if(fJetAxis == 1){
                inclusiveJetInfo[nJetsInThisEvent][1] = leadingParticleFlowCandidatePhi;
                inclusiveJetInfo[nJetsInThisEvent][2] = leadingParticleFlowCandidateEta;
              }
              
              // Correlate inclusive jets with tracks
              if(!onlyMix) CorrelateTracksAndJets(inclusiveJetInfo[nJetsInThisEvent],inclusiveJetInfo[nJetsInThisEvent],DijetHistograms::kSameEvent,true);
              
              // Increase the number of jets we have found in the event
              nJetsInThisEvent++;
              
            } // Jet passes the pT cuts
          } // Check if we fill inclusive jet-track correlation histograms
        } // Eta cut
        
      } // End of search for leading jet loop
      
      //************************************************
      //     Fill histograms for all leading jets
      //************************************************
      
      // Fill histograms for all leading jets
      if(fFillJetHistograms && highestAnyPt > 0){
        
        // Find the jet pT weight
        jetPtWeight = GetJetPtWeight(highestAnyPt);
        
        fillerJet[0] = highestAnyPt;          // Axis 0 = any leading jet pT
        fillerJet[1] = highestPhi;            // Axis 1 = any leading jet phi
        fillerJet[2] = highestEta;            // Axis 2 = any leading jet eta
        fillerJet[3] = centrality;            // Axis 3 = centrality
        fillerJet[4] = leadingJetFlavor;      // Axis 4 = any leading jet flavor
        
//        // Extra axes for the smearing study
//        fillerJet[5] = highestJetPtSmeared;   // Axis 5 = Smeared jet pT
//        fillerJet[6] = highestJetPtErrorUp;   // Axis 6 = Relative uncertainty added to the jet pT
//        fillerJet[7] = highestJetPtErrorDown; // Axis 7 = Relative uncertainty subtracted from the jet pT
//        // Extra axes for the smearing study
        
        fHistograms->fhLeadingJet->Fill(fillerJet,fTotalEventWeight*jetPtWeight);
      }
      
      //************************************************
      //              Apply dijet criteria
      //************************************************
      
      // Check if at least two jets were found
      if(highestIndex > -1 && secondHighestIndex > -1) twoJetsFound = true;
      
      // Only apply the dijet cuts for events with at least two jets
      if(twoJetsFound){
        
        // Read the eta and phi values for leading and subleading jets
        leadingJetPhi = fJetReader->GetJetPhi(highestIndex);
        leadingJetEta = fJetReader->GetJetEta(highestIndex);
        leadingJetFlavor = 0; // 0 = Quark jet
        subleadingJetPhi = fJetReader->GetJetPhi(secondHighestIndex);
        subleadingJetEta = fJetReader->GetJetEta(secondHighestIndex);
        subleadingJetFlavor = 0; // 0 = Quark jet
        
        // If matching jets or using reco jets, find the jet flavor
        if(fMatchJets || fMcCorrelationType == kRecoGen || fMcCorrelationType == kRecoReco){
          if(TMath::Abs(fJetReader->GetPartonFlavor(highestIndex)) == 21) leadingJetFlavor = 1; // 1 = Gluon jet
          if(TMath::Abs(fJetReader->GetPartonFlavor(secondHighestIndex)) == 21) subleadingJetFlavor = 1; // 1 = Gluon jet
        }
        
        // For data, instead of jet flavor, mark positive vz with 1 and negative with 0
        // This is used in one of the systematic checks for long range correlations
        // TODO: Debuggery
        //if((fDataType == ForestReader::kPp || fDataType == ForestReader::kPbPb) && vz > 0){
        //  leadingJetFlavor = 1;
        //  subleadingJetFlavor = 1;
        //}
        
        // Apply the JFF correction for leading and subleading jet pT only if we are using calo jets
        // Determine the directions of the leading particle flow candidates for the leading and subleading
        // jets also if leading particle flow candidate is used to estimate the jet axis
        if((fJetType == 0 && fReadMode < 2000) || fJetAxis == 1){
          std::tie(nParticleFlowCandidatesInThisJet,leadingParticleFlowCandidatePhi,leadingParticleFlowCandidateEta) = GetNParticleFlowCandidatesInJet(leadingJetPhi,leadingJetEta);
          if(fJetType == 0 && fReadMode < 2000) leadingJetPt =   fJffCorrection->GetCorrection(nParticleFlowCandidatesInThisJet,hiBin,leadingJetPt,leadingJetEta);
          std::tie(nParticleFlowCandidatesInThisJet,subleadingParticleFlowCandidatePhi,subleadingParticleFlowCandidateEta) = GetNParticleFlowCandidatesInJet(subleadingJetPhi,subleadingJetEta);
          if(fJetType == 0 && fReadMode < 2000) subleadingJetPt =  fJffCorrection->GetCorrection(nParticleFlowCandidatesInThisJet,hiBin,subleadingJetPt,subleadingJetEta);
        }
        
        // If after the correction subleading jet becomes leading jet and vice versa, swap the leading/subleading info
        if(subleadingJetPt > leadingJetPt){
          swapJetPt = leadingJetPt;   leadingJetPt = subleadingJetPt;    subleadingJetPt = swapJetPt;
          swapJetPhi = leadingJetPhi; leadingJetPhi = subleadingJetPhi;  subleadingJetPhi = swapJetPhi;
          swapJetEta = leadingJetEta; leadingJetEta = subleadingJetEta;  subleadingJetEta = swapJetEta;
          swapIndex = highestIndex;   highestIndex = secondHighestIndex; secondHighestIndex = swapIndex;
          swapJetPhi = leadingParticleFlowCandidatePhi;
          leadingParticleFlowCandidatePhi = subleadingParticleFlowCandidatePhi;
          subleadingParticleFlowCandidatePhi = leadingParticleFlowCandidatePhi;
          swapJetEta = leadingParticleFlowCandidateEta;
          leadingParticleFlowCandidateEta = subleadingParticleFlowCandidateEta;
          subleadingParticleFlowCandidateEta = leadingParticleFlowCandidateEta;
          if(leadingJetFlavor != subleadingJetFlavor){ // Flavor can only be 1 or 0
            leadingJetFlavor = subleadingJetFlavor;
            subleadingJetFlavor = TMath::Abs(leadingJetFlavor - 1);
          }
        }
        
        dijetFound = true;
        dphi =  leadingJetPhi - subleadingJetPhi;
        if(dphi < 0) dphi = -dphi;
        if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
        
        // Apply dijet cuts
        if((leadingJetPt >= fJetMaximumPtCut) ||              // Maximum leading jet pT cut
           (leadingJetPt <= fLeadingJetMinPtCut) ||           // Leading jet minimum pT cut
           (subleadingJetPt <= fSubleadingJetMinPtCut) ||     // Subleading jet minimum pT cut
           (TMath::Abs(leadingJetEta) >= fJetEtaCut) ||       // Leading jet eta cut
           (TMath::Abs(subleadingJetEta) >= fJetEtaCut)||     // Subleading jet eta cut
           (TMath::Abs(dphi) <= fDeltaPhiCut)){               // DeltaPhi cut
          dijetFound = false;
        }
        
        if(dijetFound) dijetCounter++;
      } // End of dijet cuts (two jets found if)
      
      // TODO: Debuggery
      if(thirdJetPt > subleadingJetPt/2 && dijetFound) {
        highThirdJet++;
        
        // For data, mark the events with high pT third jet with 1. 0 means there is no such jet
        if(fDataType == ForestReader::kPbPb || fDataType == ForestReader::kPp){
          leadingJetFlavor = 1;
          subleadingJetFlavor = 1;
        }
      }
      
      //************************************************
      //       Fill histograms for inclusive tracks
      //************************************************
      
      // Inclusive track histograms
      if(fFillTrackHistograms && !onlyMix){
        
        // Loop over all track in the event
        nTracks = fTrackReader[DijetHistograms::kSameEvent]->GetNTracks();
        for(Int_t iTrack = 0; iTrack < nTracks; iTrack++){
          
          // Check that all the track cuts are passed
          if(!PassTrackCuts(iTrack,fHistograms->fhTrackCutsInclusive,DijetHistograms::kSameEvent)) continue;
          
          // Get the efficiency correction
          trackPt = fTrackReader[DijetHistograms::kSameEvent]->GetTrackPt(iTrack);
          trackEta = fTrackReader[DijetHistograms::kSameEvent]->GetTrackEta(iTrack);
          trackPhi = fTrackReader[DijetHistograms::kSameEvent]->GetTrackPhi(iTrack);
          trackEfficiencyCorrection = GetTrackEfficiencyCorrection(DijetHistograms::kSameEvent,iTrack);
          
          // Fill track histograms
          fillerTrack[0] = trackPt;      // Axis 0: Track pT
          fillerTrack[1] = trackPhi;     // Axis 1: Track phi
          fillerTrack[2] = trackEta;     // Axis 2: Track eta
          fillerTrack[3] = centrality;   // Axis 3: Centrality
          fHistograms->fhTrackInclusive->Fill(fillerTrack,trackEfficiencyCorrection*fTotalEventWeight);  // Fill the track histogram
          fHistograms->fhTrackInclusiveUncorrected->Fill(fillerTrack,fTotalEventWeight);                 // Fill the uncorrected track histogram
          
        }
      }
      
      //************************************************
      //       Fill histograms for dijet events
      //************************************************
      
      // If we require to match the dijet, check that the matching jet pairs also fulfill the dijet definition
      matchVeto = false;
      if((dijetFound && (fMatchDijet || fMatchLeadingJet)) || (findMatchedDijet && twoJetsFound)){
        matchedLeadingJetPt = fJetReader->GetMatchedPt(highestIndex);
        matchedLeadingJetPhi = fJetReader->GetMatchedPhi(highestIndex);
        matchedLeadingJetEta = fJetReader->GetMatchedEta(highestIndex);
        matchedSubleadingJetPt = fJetReader->GetMatchedPt(secondHighestIndex);
        matchedSubleadingJetPhi = fJetReader->GetMatchedPhi(secondHighestIndex);
        matchedSubleadingJetEta = fJetReader->GetMatchedEta(secondHighestIndex);
        
        // Check for leading-subleading swapping
        if((matchedLeadingJetPt < matchedSubleadingJetPt) && fMatchDijet) {
          matchVeto = true;
          jetSwapCounter++;
        }
        
        // Swap matched leading and subleading jet values if they are in the wrong order
        if(matchedSubleadingJetPt > matchedLeadingJetPt){
          swapJetPt = matchedLeadingJetPt;   matchedLeadingJetPt = matchedSubleadingJetPt;    matchedSubleadingJetPt = swapJetPt;
          swapJetPhi = matchedLeadingJetPhi; matchedLeadingJetPhi = matchedSubleadingJetPhi;  matchedSubleadingJetPhi = swapJetPhi;
          swapJetEta = matchedLeadingJetEta; matchedLeadingJetEta = matchedSubleadingJetEta;  matchedSubleadingJetEta = swapJetEta;
        }
        
        // Calculate the deltaPhi between mathed jets
        matchedDeltaPhi =  matchedLeadingJetPhi - matchedSubleadingJetPhi;
        if(matchedDeltaPhi < 0) matchedDeltaPhi = -matchedDeltaPhi;
        if(matchedDeltaPhi > TMath::Pi()) matchedDeltaPhi = 2*TMath::Pi() - matchedDeltaPhi;
        
        // Apply dijet cuts for matched jets
        if((matchedLeadingJetPt >= fJetMaximumPtCut) ||                            // Maximum leading jet pT cut
           (matchedLeadingJetPt <= fLeadingJetMinPtCut) ||                         // Leading jet minimum pT cut
           ((matchedSubleadingJetPt <= fSubleadingJetMinPtCut) && fMatchDijet) ||  // Subleading jet minimum pT cut
           (TMath::Abs(matchedLeadingJetEta) >= fJetEtaCut) ||                     // Leading jet eta cut
           ((TMath::Abs(matchedSubleadingJetEta) >= fJetEtaCut) && fMatchDijet) || // Subleading jet eta cut
           ((TMath::Abs(matchedDeltaPhi) <= fDeltaPhiCut) && fMatchDijet)){        // DeltaPhi cut
          matchVeto = true;
        }
        
      } // Matching the jet
      
      // If reverse veto is in place, do the opposite as with regular veto
      if(reverseMatchVeto) matchVeto = !matchVeto;
      
      // If a dijet is found and not vetoed, fill some information to fHistograms
      if((dijetFound && !matchVeto) || (dijetFound && findMatchedDijet)){
                        
        matchedDijetCounter++;
        
        // Histograms for jet pT closure
        if(fFillJetPtClosure && fMatchJets){
          
          // Last argument only used for xj study
          FillJetPtClosureHistograms(highestIndex, DijetHistograms::kLeadingClosure, subleadingJetPt/leadingJetPt);
          FillJetPtClosureHistograms(secondHighestIndex, DijetHistograms::kSubleadingClosure, subleadingJetPt/leadingJetPt);
          
        }
        
        // Dijet event information
        if(fFillEventInformation){
          fHistograms->fhEvents->Fill(DijetHistograms::kDijet);
          fHistograms->fhVertexZDijet->Fill(vz,fTotalEventWeight); // TODO: Total weight here instead of vz
          fHistograms->fhCentralityDijet->Fill(centrality,fTotalEventWeight); // TODO: Total weight here instead of centrality
        }
        
        // Single jet and dijet histograms in dijet events
        if(fFillJetHistograms){
          
          // Calculate the asymmetry
          dijetXj = subleadingJetPt/leadingJetPt;
          
          dijetMatchedXj = matchedSubleadingJetPt/matchedLeadingJetPt;
          dijetAsymmetry = dijetXj;
          
          jetPtWeight = GetDijetWeight(leadingJetPt);
          triggerEfficiencyWeight = GetTriggerEfficiencyWeight(leadingJetPt, centrality);
          
//          // For smearing study, smear the jet pT and find uncertainties for jet pT
//          jetPtSmeared = leadingJetPt*fRng->Gaus(1,0.2);     // Smearing for 20 % of the jet pT
//          fJetUncertainty2018->SetJetPT(leadingJetPt);
//          fJetUncertainty2018->SetJetEta(leadingJetEta);
//          fJetUncertainty2018->SetJetPhi(leadingJetPhi);
//          jetPtErrorUp = fJetUncertainty2018->GetUncertainty().second;
//          jetPtErrorDown = fJetUncertainty2018->GetUncertainty().first;
          
          // Fill the leading jet histogram
          fillerLeadingJet[0] = leadingJetPt;                   // Axis 0: Leading jet pT
          fillerLeadingJet[1] = leadingJetPhi;                  // Axis 1: Leading jet phi
          fillerLeadingJet[2] = leadingJetEta;                  // Axis 2: Leading jet eta
          fillerLeadingJet[3] = dijetAsymmetry;                 // Axis 3: Asymmetry
          fillerLeadingJet[4] = centrality;                     // Axis 4: Centrality
          fillerLeadingJet[5] = leadingJetFlavor;               // Axis 5: Leading jet flavor
          
//          // Extra axes filled for the smearing study
//          fillerLeadingJet[6] = jetPtSmeared;                   // Axis 6: Smeared leading jet pT
//          fillerLeadingJet[7] = jetPtErrorUp;                   // Axis 7: Uncertainty added to the leading jet pT
//          fillerLeadingJet[8] = jetPtErrorDown;                 // Axis 8: Uncertainty subtracted from the leading jet pT
          
          fHistograms->fhLeadingDijet->Fill(fillerLeadingJet,fTotalEventWeight*jetPtWeight*triggerEfficiencyWeight);    // Fill the data point to leading jet histogram
               
//          // For smearing study, smear the jet pT and find uncertainties for jet pT
//          jetPtSmeared = subleadingJetPt*fRng->Gaus(1,0.2);     // Smearing for 20 % of the jet pT
//          fJetUncertainty2018->SetJetPT(subleadingJetPt);
//          fJetUncertainty2018->SetJetEta(subleadingJetEta);
//          fJetUncertainty2018->SetJetPhi(subleadingJetPhi);
//          jetPtErrorUp = fJetUncertainty2018->GetUncertainty().second;
//          jetPtErrorDown = fJetUncertainty2018->GetUncertainty().first;
          
          // Fill the subleading jet histogram
          fillerLeadingJet[0] = subleadingJetPt;                // Axis 0: Subleading jet pT
          fillerLeadingJet[1] = subleadingJetPhi;               // Axis 1: Subleading jet phi
          fillerLeadingJet[2] = subleadingJetEta;               // Axis 2: Subleading jet eta
          fillerLeadingJet[3] = dijetAsymmetry;                 // Axis 3: Asymmetry
          fillerLeadingJet[4] = centrality;                     // Axis 4: Centrality
          fillerLeadingJet[5] = subleadingJetFlavor;            // Axis 5: Subleading jet flavor
          
//          // Extra axes filled for the smearing study
//          fillerLeadingJet[6] = jetPtSmeared;                   // Axis 6: Smeared subleading jet pT
//          fillerLeadingJet[7] = jetPtErrorUp;                   // Axis 7: Uncertainty added to the subleading jet pT
//          fillerLeadingJet[8] = jetPtErrorDown;                 // Axis 8: Uncertainty subtracted from the subleading jet pT
          
          fHistograms->fhSubleadingDijet->Fill(fillerLeadingJet,fTotalEventWeight*jetPtWeight*triggerEfficiencyWeight); // Fill the data point to subleading jet histogram
          
          // Fill the dijet histogram
          fillerDijet[0] = leadingJetPt;                   // Axis 0: Leading jet pT
          fillerDijet[1] = subleadingJetPt;                // Axis 1: Subleading jet pT
          fillerDijet[2] = TMath::Abs(dphi);               // Axis 2: deltaPhi
          fillerDijet[3] = dijetXj;                        // Axis 3: Dijet momentum balance xj
          fillerDijet[4] = centrality;                     // Axis 4: Centrality
          fillerDijet[5] = dijetMatchedXj;                 // Axis 5: Matched dijet momentum balance xj
          
//          // Extra axes for xj study
//          fillerDijet[6] = leadingJetPhi;                  // Axis 1: Leading jet phi
//          fillerDijet[7] = leadingJetEta;                  // Axis 2: Leading jet eta
//          fillerDijet[8] = subleadingJetPhi;               // Axis 1: Subleading jet phi
//          fillerDijet[9] = subleadingJetEta;               // Axis 2: Subleading jet eta
          
          fHistograms->fhDijet->Fill(fillerDijet,fTotalEventWeight*jetPtWeight*triggerEfficiencyWeight);         // Fill the data point to dijet histogram
          
        }
        
        // Fill the arrays with leading and subleading jet information for correlation with tracks
        leadingJetInfo[0] = leadingJetPt;
        leadingJetInfo[1] = leadingJetPhi;
        leadingJetInfo[2] = leadingJetEta;
        leadingJetInfo[3] = leadingJetFlavor;
        subleadingJetInfo[0] = subleadingJetPt;
        subleadingJetInfo[1] = subleadingJetPhi;
        subleadingJetInfo[2] = subleadingJetEta;
        subleadingJetInfo[3] = subleadingJetFlavor;
        
        // If we are using leading particle flow candidate as a probe for jet axis, change the jet info
        if(fJetAxis == 1){
          leadingJetInfo[1] = leadingParticleFlowCandidatePhi;
          leadingJetInfo[2] = leadingParticleFlowCandidateEta;
          subleadingJetInfo[1] = subleadingParticleFlowCandidatePhi;
          subleadingJetInfo[2] = subleadingParticleFlowCandidateEta;
        }
        
        // Do not do the jet-track correlation are not filling the relevant histograms
        if(fFillDijetJetTrackCorrelation){
          
          // Correlate jets with tracks in dijet events
          if(!onlyMix) CorrelateTracksAndJets(leadingJetInfo,subleadingJetInfo,DijetHistograms::kSameEvent);
          
        }
          
      } // Dijet in event
      
      // Do event mixing for tracks and inclusive jets, leading jet and subleading jet using the selected mixing method
      // This is best done for everything at once to avoid loading mixing event several times (slow process)
      // Do not go to mixing if no inclusive jets are found and we are doing inclusive jet mixing or no dijet are found
      // and we are doing dijet mixing
      if(mixEvents && ((fFillInclusiveJetTrackCorrelation && (nJetsInThisEvent > 0)) || (fFillDijetJetTrackCorrelation && (dijetFound && !matchVeto)))){
        if(mixWithPool) {
          MixTracksAndJets(inclusiveJetInfo,leadingJetInfo,subleadingJetInfo,iEvent,nJetsInThisEvent,vz,hiBin,(dijetFound && !matchVeto));
        } else {
          MixTracksAndJetsWithoutPool(inclusiveJetInfo,leadingJetInfo,subleadingJetInfo,iEvent,nJetsInThisEvent,vz,hiBin,(dijetFound && !matchVeto));
        }
      }
      
    } // Event loop
    
    //************************************************
    //      Cleanup at the end of the file loop
    //************************************************
    
    // Close the input files after the event has been read
    inputFile->Close();
    if(useDifferentReaderFotJetsAndTracks) copyInputFile->Close();
    if(mixEvents && fDataType != ForestReader::kPbPb && fDataType != ForestReader::kPbPbMC) mixedEventFile->Close();
    
  } // File loop
  
  // TODO: Debuggery
  cout << "Number of high third jet events: " << highThirdJet << endl;
  cout << "Number of dijets: " << dijetCounter << endl;
  if(fMatchJets && fDebugLevel > 1){
    cout << "A total of " << unmatchedCounter << " jets were not matched!" << endl;
    cout << "A total of " << nonSensicalPartonIndex << " out of " << matchedCounter <<  " matched jets had reference parton flavor -999!" << endl;
    cout << "All dijets: " << dijetCounter << ". Matched dijets: " << matchedDijetCounter << ". Not matched due to leading/subleading swap: " << jetSwapCounter << endl;
  }
  
}

/*
 * Method for all jet-track correlations
 *
 *  const Double_t leadingJetInfo[3] = Array containing leading jet pT, phi and eta
 *  const Double_t subleadingJetInfo[3] = Array containing subleading jet pT, phi and eta
 *  const Int_t correlationType = DijetHistograms::kSameEvent for same event correlations, DijetHistograms::kMixedEvent for mixed event correlations
 *  const Bool_t useInclusiveJets = True: Correlation done for inclusive jets. False: Correlation done for leading and subleading jets
 */
void DijetAnalyzer::CorrelateTracksAndJets(const Double_t leadingJetInfo[4], const Double_t subleadingJetInfo[4], const Int_t correlationType, const Bool_t useInclusiveJets){
  
  // Define a filler for THnSparses
  Double_t fillerJetTrack[7];
  Double_t fillerJetTrackInclusive[6];
  Double_t fillerTrack[5];
  
  // Event information
  Double_t centrality = fTrackReader[correlationType]->GetCentrality();
  Double_t jetPtWeight = 1;
  Double_t triggerEfficiencyWeight = 1;
  
  // Variables for tracks
  Double_t trackPt;       // Track pT
  Double_t trackEta;      // Track eta
  Double_t trackPhi;      // Track phi
  Double_t trackEfficiencyCorrection;  // Efficiency correction for the track
  Double_t deltaPhiTrackLeadingJet;    // DeltaPhi between track and leading jet
  Double_t deltaEtaTrackLeadingJet;    // DeltaEta between track and leading jet
  Double_t deltaPhiTrackSubleadingJet; // DeltaPhi between track and subleading jet
  Double_t deltaEtaTrackSubleadingJet; // DeltaEta between track and subleading jet
  
  // Read the leading jet and subleading jet information from arrays
  Double_t leadingJetPt = leadingJetInfo[0];
  Double_t leadingJetPhi = leadingJetInfo[1];
  Double_t leadingJetEta = leadingJetInfo[2];
  Double_t leadingJetFlavor = leadingJetInfo[3];
  Double_t subleadingJetPt = subleadingJetInfo[0];
  Double_t subleadingJetPhi = subleadingJetInfo[1];
  Double_t subleadingJetEta = subleadingJetInfo[2];
  Double_t subleadingJetFlavor = subleadingJetInfo[3];
  
  // Calculate the dijet asymmetry. Choose AJ or xJ based on selection from JCard.
  Double_t dijetAsymmetry = (fAsymmetryBinType == 0) ? (leadingJetPt - subleadingJetPt)/(leadingJetPt + subleadingJetPt) : subleadingJetPt/leadingJetPt;
  
  // Loop over all track in the event
  Int_t nTracks = fTrackReader[correlationType]->GetNTracks();
  for(Int_t iTrack = 0; iTrack < nTracks; iTrack++){
    
    // Check that all the track cuts are passed
    if(!PassTrackCuts(iTrack,fHistograms->fhTrackCuts,correlationType)) continue;
    
    // Get the most important track information to variables
    trackPt = fTrackReader[correlationType]->GetTrackPt(iTrack);
    trackPhi = fTrackReader[correlationType]->GetTrackPhi(iTrack);
    trackEta = fTrackReader[correlationType]->GetTrackEta(iTrack);
    
    // Get the efficiency correction
    trackEfficiencyCorrection = GetTrackEfficiencyCorrection(correlationType,iTrack);
    
    // Calculate deltaEta and deltaPhi between track and leading and subleading jets
    deltaEtaTrackLeadingJet = leadingJetEta - trackEta;
    deltaPhiTrackLeadingJet = leadingJetPhi - trackPhi;
    deltaEtaTrackSubleadingJet = subleadingJetEta - trackEta;
    deltaPhiTrackSubleadingJet = subleadingJetPhi - trackPhi;
    
    // Transform deltaPhis to interval [-pi/2,3pi/2]
    while(deltaPhiTrackLeadingJet > (1.5*TMath::Pi())){deltaPhiTrackLeadingJet += -2*TMath::Pi();}
    while(deltaPhiTrackSubleadingJet > (1.5*TMath::Pi())){deltaPhiTrackSubleadingJet += -2*TMath::Pi();}
    while(deltaPhiTrackLeadingJet < (-0.5*TMath::Pi())){deltaPhiTrackLeadingJet += 2*TMath::Pi();}
    while(deltaPhiTrackSubleadingJet < (-0.5*TMath::Pi())){deltaPhiTrackSubleadingJet += 2*TMath::Pi();}
    
    // Fill track histograms if selected to do so. Inclusive tracks are filled in separate place
    if(fFillTrackHistograms && !useInclusiveJets){
      fillerTrack[0] = trackPt;                    // Axis 0: Track pT
      fillerTrack[1] = trackPhi;                   // Axis 1: Track phi
      fillerTrack[2] = trackEta;                   // Axis 2: Track eta
      fillerTrack[3] = centrality;                 // Axis 3: Centrality
      fillerTrack[4] = correlationType;            // Axis 4: Correlation type (same or mixed event)
      fHistograms->fhTrack->Fill(fillerTrack,trackEfficiencyCorrection*fTotalEventWeight);  // Fill the track histogram
      fHistograms->fhTrackUncorrected->Fill(fillerTrack,fTotalEventWeight);                 // Fill the uncorrected track histogram
    }
    
    // Fill the selected jet-track correlation histograms
    
    if(useInclusiveJets){
      
      jetPtWeight = GetJetPtWeight(leadingJetPt);
      triggerEfficiencyWeight = GetTriggerEfficiencyWeight(leadingJetPt, centrality);
      
      fillerJetTrackInclusive[0] = trackPt;                    // Axis 0: Track pT
      fillerJetTrackInclusive[1] = deltaPhiTrackLeadingJet;    // Axis 1: DeltaPhi between track and inclusive jet
      fillerJetTrackInclusive[2] = deltaEtaTrackLeadingJet;    // Axis 2: DeltaEta between track and inclusive jet
      fillerJetTrackInclusive[3] = centrality;                 // Axis 3: Centrality
      fillerJetTrackInclusive[4] = correlationType;            // Axis 4: Correlation type (same or mixed event)
      fillerJetTrackInclusive[5] = leadingJetFlavor;           // Axis 5: Jet flavor (quark of gluon)
      
      if(fFillInclusiveJetTrackCorrelation){
        fHistograms->fhTrackJetInclusive->Fill(fillerJetTrackInclusive, trackEfficiencyCorrection*fTotalEventWeight*jetPtWeight*triggerEfficiencyWeight); // Fill the track-inclusive jet correlation histogram
        fHistograms->fhTrackJetInclusivePtWeighted->Fill(fillerJetTrackInclusive, trackEfficiencyCorrection*trackPt*fTotalEventWeight*jetPtWeight*triggerEfficiencyWeight); // Fill the track-inclusive jet correlation histogram
      }
    } else {
      
      jetPtWeight = GetDijetWeight(leadingJetPt);
      triggerEfficiencyWeight = GetTriggerEfficiencyWeight(leadingJetPt, centrality);
      
      // Fill the track-leading jet correlation histograms
      fillerJetTrack[0] = trackPt;                    // Axis 0: Track pT
      fillerJetTrack[1] = deltaPhiTrackLeadingJet;    // Axis 1: DeltaPhi between track and leading jet
      fillerJetTrack[2] = deltaEtaTrackLeadingJet;    // Axis 2: DeltaEta between track and leading jet
      fillerJetTrack[3] = dijetAsymmetry;             // Axis 3: Dijet asymmetry
      fillerJetTrack[4] = centrality;                 // Axis 4: Centrality
      fillerJetTrack[5] = correlationType;            // Axis 5: Correlation type (same or mixed event)
      fillerJetTrack[6] = leadingJetFlavor;           // Axis 6: Leading jet flavor (quark or gluon)
      if(fFillRegularJetTrackCorrelation) fHistograms->fhTrackLeadingJet->Fill(fillerJetTrack, trackEfficiencyCorrection*fTotalEventWeight*jetPtWeight*triggerEfficiencyWeight); // Fill the track-leading jet correlation histogram
      if(fFillUncorrectedJetTrackCorrelation) fHistograms->fhTrackLeadingJetUncorrected->Fill(fillerJetTrack, fTotalEventWeight*jetPtWeight*triggerEfficiencyWeight);                // Fill the uncorrected track-leading jet correlation histogram
      if(fFillPtWeightedJetTrackCorrelation) fHistograms->fhTrackLeadingJetPtWeighted->Fill(fillerJetTrack, trackEfficiencyCorrection*trackPt*fTotalEventWeight*jetPtWeight*triggerEfficiencyWeight); // Fill the pT weighted track-leading jet correlation histogram
            
      // Fill the track-subleading jet correlation histograms
      fillerJetTrack[0] = trackPt;                    // Axis 0: Track pT
      fillerJetTrack[1] = deltaPhiTrackSubleadingJet; // Axis 1: DeltaPhi between track and subleading jet
      fillerJetTrack[2] = deltaEtaTrackSubleadingJet; // Axis 2: DeltaEta between track and subleading jet
      fillerJetTrack[3] = dijetAsymmetry;             // Axis 3: Dijet asymmetry
      fillerJetTrack[4] = centrality;                 // Axis 4: Centrality
      fillerJetTrack[5] = correlationType;            // Axis 5: Correlation type (same or mixed event)
      fillerJetTrack[6] = subleadingJetFlavor;        // Axis 6: Subleading jet flavor (quark or gluon)
      if(fFillRegularJetTrackCorrelation) fHistograms->fhTrackSubleadingJet->Fill(fillerJetTrack, trackEfficiencyCorrection*fTotalEventWeight*jetPtWeight*triggerEfficiencyWeight); // Fill the track-subleading jet correlation histogram
      if(fFillUncorrectedJetTrackCorrelation) fHistograms->fhTrackSubleadingJetUncorrected->Fill(fillerJetTrack, fTotalEventWeight*jetPtWeight*triggerEfficiencyWeight);                // Fill the uncorrected track-subleading jet correlation histogram
      if(fFillPtWeightedJetTrackCorrelation) fHistograms->fhTrackSubleadingJetPtWeighted->Fill(fillerJetTrack, trackEfficiencyCorrection*trackPt*fTotalEventWeight*jetPtWeight*triggerEfficiencyWeight); // Fill the pT weighted track-subleading jet correlation histogram
    }
    
  } // Loop over tracks
}

/*
 * Fill the jet pT closure histograms
 *
 *  const Int_t jetIndex = Index of a jet for which the closure is filled
 *  const Int_t closureType = Leading/subleading/inclusive
 *  const Double_t xj = Dijet momentum balance for the dijet part of which this jet is
 */
void DijetAnalyzer::FillJetPtClosureHistograms(const Int_t jetIndex, const Int_t closureType, const Double_t xj){

  // Define a filler for the closure histogram
  const Int_t nAxesClosure = 8; // 7 nominal, 8 used for xj study
  Double_t fillerClosure[nAxesClosure];
  
  // Find the pT of the matched gen jet and flavor of reference parton
  Float_t matchedGenPt = fJetReader->GetMatchedPt(jetIndex);
  Float_t matchedGenEta = fJetReader->GetMatchedEta(jetIndex);
  Float_t matchedGenPhi = fJetReader->GetMatchedPhi(jetIndex);
  Int_t referencePartonFlavor = fJetReader->GetPartonFlavor(jetIndex);
  Int_t hiBin = fJetReader->GetHiBin();
  
  // Find the centrality of the event and the pT of the reconstructed jet
  Double_t recoPt = fJetReader->GetJetPt(jetIndex);
  Double_t centrality = fJetReader->GetCentrality();
  Double_t jetEta = fJetReader->GetJetEta(jetIndex);
  Double_t jetPhi = fJetReader->GetJetPhi(jetIndex);
  
  // If we are using generator level jets, swap reco and gen variables
  if(fMcCorrelationType == kGenReco || fMcCorrelationType == kGenGen){
    Double_t swapper = matchedGenPt;
    matchedGenPt = recoPt;
    recoPt = swapper;
    swapper = matchedGenEta;
    matchedGenEta = jetEta;
    jetEta = swapper;
    swapper = matchedGenPhi;
    matchedGenPhi = jetPhi;
    jetPhi = swapper;
  }
  
  // Helper variable for smearing study
  Double_t smearingFactor;
  
  // For 2018 data, we need to correct the reconstructed pT with jet energy correction
  fJetCorrector2018->SetJetPT(recoPt);
  fJetCorrector2018->SetJetEta(jetEta);
  fJetCorrector2018->SetJetPhi(jetPhi);
  
  if(fReadMode > 2000){
    recoPt = fJetCorrector2018->GetCorrectedPT();
    
    // If we are using smearing scenario, modify the reconstructed jet pT using gaussian smearing
    if(fJetUncertaintyMode == 3){
      smearingFactor = GetSmearingFactor(recoPt, centrality);
      recoPt = recoPt * fRng->Gaus(1,smearingFactor);
    }
    
  }
  
  // Apply JFF correction for jet pT only if we are using calorimeter jets in 2015 data
  if(fJetType == 0 && fReadMode < 2000){
    int nParticleFlowCandidatesInThisJet;
    double leadingParticleFlowCandidatePhi;
    double leadingParticleFlowCandidateEta;
    std::tie(nParticleFlowCandidatesInThisJet,leadingParticleFlowCandidatePhi,leadingParticleFlowCandidateEta) = GetNParticleFlowCandidatesInJet(jetPhi,jetEta);
    recoPt = fJffCorrection->GetCorrection(nParticleFlowCandidatesInThisJet,hiBin,recoPt,jetEta);
  }
  
  // Define index for parton flavor using algoritm: [-6,-1] U [1,6] -> kQuark, 21 -> kGluon, anything else -> -1
  Int_t referencePartonIndex = -1;
  if(referencePartonFlavor >= -6 && referencePartonFlavor <= 6 && referencePartonFlavor != 0) referencePartonIndex = DijetHistograms::kQuark;
  if(referencePartonFlavor == 21) referencePartonIndex = DijetHistograms::kGluon;
  
  // Fill the different axes for the filler
  fillerClosure[0] = closureType;          // Axis 0: Type of closure (leading/subleading/inclusive)
  fillerClosure[1] = matchedGenPt;         // Axis 1: pT of the matched generator level jet
  fillerClosure[2] = recoPt;               // Axis 2: pT of the matched reconstructed jet
  fillerClosure[3] = jetEta;               // Axis 3: eta of the jet under consideration
  fillerClosure[4] = centrality;           // Axis 4: Centrality of the event
  fillerClosure[5] = referencePartonIndex; // Axis 5: Reference parton type (quark/gluon)
  fillerClosure[6] = recoPt/matchedGenPt;  // Axis 6: Reconstructed level jet to generator level jet pT ratio
  
  // Extra axis used for xj study
  fillerClosure[7] = xj;                   // Axis 7: Dijet momentum balance
  
  // Fill the closure histogram
  fHistograms->fhJetPtClosure->Fill(fillerClosure,fTotalEventWeight);
  
}

/*
 * Do the jet-track correlations with mixed events
 *
 *  const Double_t inclusiveJetInfo[60][4] = Array containing pT, phi and eta for all jets above the leading jet cut in the event
 *  const Double_t leadingJetInfo[4] = Array containing leading jet pT, phi and eta
 *  const Double_t subleadingJetInfo[4] = Array containing subleading jet pT, phi and eta
 *  const Int_t avoidIndex = Index of the current event. Do not use this index for mixing in order not to mix with the same event
 *  const Int_t nJetsInThisEvent = Number of jets above the leading jet cut in this event
 *  const Double_t vz = Vertex z-position for the main event
 *  const Int_t hiBin = HiBin of the main event
 *  const Bool_t dijetInEvent = True: Event has dijet. False: No dijet in this event
 */
void DijetAnalyzer::MixTracksAndJets(const Double_t inclusiveJetInfo[60][4], const Double_t leadingJetInfo[4], const Double_t subleadingJetInfo[4], const Int_t avoidIndex, const Int_t nJetsInThisEvent, const Double_t vz, const Int_t hiBin, const Bool_t dijetInEvent){
  
  // Start mixing from the first event index
  Int_t mixedEventIndex;                        // Index of current event in mixing loop
  Int_t mixedEventPosition;                     // Position of the mixed event in the event pool
  Int_t eventsMixed = 0;                        // Number of events mixed thus far
  Int_t mixingVzBin = FindMixingVzBin(vz);      // Bin index for the vz bin in mixing pool
  Int_t mixingHiBin = FindMixingHiBin(hiBin);   // Bin index for the centrality bin in the mixing pool
  //Bool_t startOfLoop = true;                    // Bool for writing debugging messages is same event are used many times for mixing
  //Int_t nLoops = 0;                             // Counter of how many times same event are used for mixing
  //Int_t firstMixingIndex = 0;                   // Index for the first event used in event mixing
  
  // Initialize the mixed event randomizer
  TRandom3 *mixedEventRandomizer = new TRandom3();  // Randomizer for starting point in the mixed event file
  mixedEventRandomizer->SetSeed(0);
  
  // Continue mixing until we have reached required number of events
  while (eventsMixed < fnMixedEventsPerDijet) {
    
    // Original implementation for the mixing pool using running index to go over pool and a message if events are used more than once
    /*
    // If the running index is outside of the vector range, start over
    if(fRunningMixingIndex >= fMixingPool[mixingVzBin][mixingHiBin].size()) fRunningMixingIndex = 0;
    
    // If we are back at the first index, we are going to the next loop over the vector in mixing pool
    if(fRunningMixingIndex == firstMixingIndex){
      startOfLoop = true;
    }
    
    // At the start of loop over mixing pool, increment loop counter
    if(startOfLoop){
      startOfLoop = false;
      firstMixingIndex = fRunningMixingIndex;
      nLoops++;
      
      // Write a debug message to the console if same events are used for mixing several times
      if(nLoops > 1 && fDebugLevel > 0){
        cout << "Note! In the mixing for event " << avoidIndex << " same mixing events are used for " << nLoops << " times!" << endl;
        cout << "Size of the mixing vector: " << fMixingPool[mixingVzBin][mixingHiBin].size() << endl;
      }
    }*/
    
    // Implementation in Kurt's code: Read a random event to mix with the current event. No check for duplicates.
    
    // Get the index for mixing event from a random position in the mixed event pool
    mixedEventPosition = mixedEventRandomizer->Rndm()*(fMixingPool[mixingVzBin][mixingHiBin].size());
    // Note: use size instead of size-1 because conversion from double to int only checks the integer part of the number
    // and using size-1 makes it impossible to reach the last index in vector (unless random floating point number is exactly 1)
    // If this extremely rare case happens, change the index by one to stay in bounds of the vector.
    if(mixedEventPosition == fMixingPool[mixingVzBin][mixingHiBin].size()) mixedEventPosition--;
    mixedEventIndex = fMixingPool[mixingVzBin][mixingHiBin].at(mixedEventPosition);
    
    // Never mix with an event having the same index to avoid same event correlations
    if(mixedEventIndex == avoidIndex) continue;
    
    // Get the event defined by the index in the mixing pool
    fTrackReader[DijetHistograms::kMixedEvent]->GetEvent(mixedEventIndex);
    
    // Do the correlation with a jet from the current event and track from mixing event
    if(fFillInclusiveJetTrackCorrelation){
      for(int iJet = 0; iJet < nJetsInThisEvent; iJet++){
        CorrelateTracksAndJets(inclusiveJetInfo[iJet],inclusiveJetInfo[iJet],DijetHistograms::kMixedEvent,true);
      }
    }
    
    // Do the correlations with the dijet from current event and tracks from mixing event
    if(fFillDijetJetTrackCorrelation && dijetInEvent){
      CorrelateTracksAndJets(leadingJetInfo,subleadingJetInfo,DijetHistograms::kMixedEvent,false);
    }
    eventsMixed++;
    
  } // While loop for finding events to mix
  
  delete mixedEventRandomizer;
  
}

/*
 * Do the jet-track correlations with mixed events in a poolless way
 *
 *  const Double_t inclusiveJetInfo[60][3] = Array containing pT, phi and eta for all jets above the leading jet cut in the event
 *  const Double_t leadingJetInfo[3] = Array containing leading jet pT, phi and eta
 *  const Double_t subleadingJetInfo[3] = Array containing subleading jet pT, phi and eta
 *  const Int_t avoidIndex = Index of the current event. Do not use this index for mixing in order not to mix with the same event
 *  const Int_t nJetsInThisEvent = Number of jets above the leading jet cut in this event
 *  const Double_t vz = Vertex z-position for the main event
 *  const Int_t hiBin = HiBin of the main event
 *  const Bool_t dijetInEvent = True: Event has dijet. False: No dijet in this event
 */
void DijetAnalyzer::MixTracksAndJetsWithoutPool(const Double_t inclusiveJetInfo[60][4], const Double_t leadingJetInfo[4], const Double_t subleadingJetInfo[4], const Int_t avoidIndex, const Int_t nJetsInThisEvent, const Double_t vz, const Int_t hiBin, const Bool_t dijetInEvent){
  
  // Start mixing from the first event index
  Int_t mixedEventIndex = fMixingStartIndex; // Index of current event in mixing loop
  Int_t eventsMixed = 0;                     // Number of events mixed thus far
  Bool_t allEventsWentThrough = false;       // Flag if we have gone through all the events in mixing file without finding enough events to mix with
  Bool_t checkForDuplicates = false;         // Start checking for duplicates in the second round over the file
  Bool_t skipEvent = false;                  // Flag if we should skip the event we are looking at the moment
  Double_t additionalTolerance = 0;          // Additional tolerance added for vz if there are not enough events in the defined range
  Int_t hiBinTolerance = 0;                  // Additional tolerance added for HiBin if there are not enough events in the defined range
  std::vector<Int_t> mixedEventIndices;      // Vector for holding the event we have already mixed with
  
  // Clear the vector before going into mixing
  mixedEventIndices.clear();
  
  // Continue mixing until we have reached required number of evens
  while (eventsMixed < fnMixedEventsPerDijet) {
    
    // By default, do not skip an event
    skipEvent = false;
    
    // If we have checked all the events but not found enough event to mix with, increase vz and hiBin tolerance
    if(allEventsWentThrough){
      if(fDebugLevel > 0){
        cout << "Could only find " << eventsMixed << " events to mix with event " << avoidIndex << endl;
        cout << "Increasing vz tolerance by 0.2 and hiBin tolerance by 1" << endl;
      }
      
      hiBinTolerance += 1;
      additionalTolerance += 0.2;
      allEventsWentThrough = false;
      checkForDuplicates = true;
    }
    
    // Increment the counter for event index to be mixed with the current event
    mixedEventIndex++;
    
    // If we are out of bounds from the event in data file, reset the counter
    if(mixedEventIndex == fnEventsInMixingFile) {
      mixedEventIndex = -1;
      continue;
    }
    
    // If we come back to the first event index, we have gone through all the events without finding 20 similar events form the file
    if(mixedEventIndex == fMixingStartIndex) allEventsWentThrough = true;
    
    // Do not mix with the same event
    if(mixedEventIndex == avoidIndex) continue;
    
    // Do not mix with events used in the previous iteration over the file
    if(checkForDuplicates){
      for(Int_t iMixedEvent : mixedEventIndices) {
        if(iMixedEvent == mixedEventIndex) skipEvent = true;
      }
    }
    if(skipEvent) continue;
    
    // Match vz and hiBin between the dijet event and mixed event
    if(TMath::Abs(fMixedEventVz.at(mixedEventIndex) - vz) > (fMixingVzTolerance + additionalTolerance)) continue;
    if(TMath::Abs(fMixedEventHiBin.at(mixedEventIndex) - hiBin) > hiBinTolerance + 1e-4) continue;
    
    // If match vz and hiBin, then load the event from the mixed event tree and mark that we have mixed this event
    fTrackReader[DijetHistograms::kMixedEvent]->GetEvent(mixedEventIndex);
    if(CheckForSameEvent(avoidIndex,mixedEventIndex)) continue;  // Check that the mixing file does not have exactly the same event as regular data file
    mixedEventIndices.push_back(mixedEventIndex);
    
    // Do the correlation with a jet from the current event and track from mixing event
    if(fFillInclusiveJetTrackCorrelation){
      for(int iJet = 0; iJet < nJetsInThisEvent; iJet++){
        CorrelateTracksAndJets(inclusiveJetInfo[iJet],inclusiveJetInfo[iJet],DijetHistograms::kMixedEvent,true);
      }
    }
    
    // Do the correlations with the dijet from current event and tracks from mixing event
    if(fFillDijetJetTrackCorrelation && dijetInEvent){
      CorrelateTracksAndJets(leadingJetInfo,subleadingJetInfo,DijetHistograms::kMixedEvent,false);
    }
    
    eventsMixed++;
    
  } // While loop for finding events to mix
  
  // For the next event, start mixing the events from where we were left with in the previous event
  fMixingStartIndex = mixedEventIndex;
  
}

/*
 * Get the proper vz weighting depending on analyzed system
 *
 *  Arguments:
 *   const Double_t vz = Vertex z position for the event
 *
 *   return: Multiplicative correction factor for vz
 */
Double_t DijetAnalyzer::GetVzWeight(const Double_t vz) const{
  if(fDataType == ForestReader::kPp || fDataType == ForestReader::kPbPb) return 1;  // No correction for real data
  if((fDataType == ForestReader::kPpMC && fReadMode < 2000) || fDataType == ForestReader::kLocalTest) return 1.0/fVzWeightFunction->Eval(vz); // Weight for 2015 pp MC
  if(fDataType == ForestReader::kPbPbMC || fDataType == ForestReader::kPpMC) return fVzWeightFunction->Eval(vz); // Weight for 2018 MC
  return -1; // Return crazy value for unknown data types, so user will not miss it
}

/*
 * Get the proper centrality weighting depending on analyzed system
 *
 *  Arguments:
 *   const Int_t hiBin = CMS hiBin
 *
 *   return: Multiplicative correction factor for the given CMS hiBin
 */
Double_t DijetAnalyzer::GetCentralityWeight(const Int_t hiBin) const{
  if(fDataType != ForestReader::kPbPbMC) return 1;
  
  // Different range for centrality weight function for 2015 and 2018.
  if(fReadMode < 2000){
    return (hiBin < 194) ? fCentralityWeightFunction->Eval(hiBin) : 1;  // No weighting for the most peripheral centrality bins 2015
  } else {
    return (hiBin < 194) ? fCentralityWeightFunction->Eval(hiBin/2.0) : 1;  // No weighting for the most peripheral centrality bins 2018
  }
}

/*
 * Get the proper jet pT weighting depending on analyzed system
 *
 *  Arguments:
 *   const Double_t jetPt = Jet pT for the weighted jet
 *
 *   return: Multiplicative correction factor for the jet pT
 */
Double_t DijetAnalyzer::GetJetPtWeight(const Double_t jetPt) const{
  if(fDataType == ForestReader::kPbPb || fDataType == ForestReader::kPp) return 1.0;  // No weight for data
  
  // Only weight 2017 and 2018 MC
  if(fReadMode < 2000){
    return 1.0;  // No weighting for the most peripheral centrality bins 2015
  } else {
    return fPtWeightFunction->Eval(jetPt);  // No weighting for the most peripheral centrality bins 2018
  }
}

/*
 * Get a trigger efficiency weight
 *
 *  The numbers are obtained with centrality binning 0-10, 10-30, 30-50, 50-90.
 *  If different binning is used, you will get wrong weights using this function.
 *
 *  Arguments:
 *   const Double_t jetPt = Jet pT for triggering jet
 *   const Double_t centrality = Centrality of the event
 *
 *   return: Multiplicative correction factor to account for trigger efficiency
 */
Double_t DijetAnalyzer::GetTriggerEfficiencyWeight(const Double_t jetPt, const Double_t centrality) const{
  return 1;  // Trigger efficiency weighting is disabled
  if(fDataType == ForestReader::kPp || fDataType == ForestReader::kPpMC) return 1; // No weight for pp or Pythia.
  if(fDataType == ForestReader::kPbPbMC && fReadMode < 2019) return 1;  // No weight if trigger not used in PbPb MC
  
  // For jets above 250 GeV, the efficiency is 1 within errors.
  if(jetPt >= 250) return 1;
  
  // Efficiency tables for different centralities and jet pT:s
  // In the tables each number represents trigger efficiency in a 5 GeV wide bin starting from 0 up until 500 GeV.
  double triggerEfficiencyTable[4][100] = {
    // Centrality: 0-10 %
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0752386, 0.0869706, 0.111874, 0.142504, 0.183662, 0.235107, 0.286736, 0.348411, 0.405188, 0.464605, 0.528057, 0.586119, 0.643719, 0.697473, 0.747104, 0.791623, 0.8325, 0.865514, 0.896188, 0.920954, 0.939565, 0.95004, 0.963755, 0.966277, 0.971066, 0.977124, 0.982302, 0.979379, 0.981975, 0.984987, 0.985513, 0.983373, 0.987799, 0.988764, 0.990123, 0.987269, 0.990842, 0.991504, 0.991676, 0.993958, 0.993402, 0.990756, 0.985158, 0.990476, 0.993007, 0.994681, 0.994228, 0.995253, 0.991135, 0.995927, 0.993318, 0.997347, 1, 1, 1, 0.996226, 0.995885, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
    // Centrality: 10-30 %
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0776977, 0.0938393, 0.122861, 0.161901, 0.217981, 0.275934, 0.337427, 0.400179, 0.458118, 0.512993, 0.579228, 0.630421, 0.686575, 0.736895, 0.781909, 0.831157, 0.866775, 0.895112, 0.92346, 0.942908, 0.955689, 0.967166, 0.973694, 0.977134, 0.982978, 0.984804, 0.987916, 0.987677, 0.988022, 0.989157, 0.991467, 0.991854, 0.991879, 0.99356, 0.991771, 0.993716, 0.99505, 0.995378, 0.993445, 0.994758, 0.996177, 0.996522, 1, 0.996855, 1, 1, 0.996914, 0.998296, 0.998134, 0.997904, 1, 1, 1, 1, 0.996109, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.985714, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.952381, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
    // Centrality: 30-50 %
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0225389, 0.0357853, 0.0625854, 0.114652, 0.186116, 0.26918, 0.352576, 0.424547, 0.487839, 0.553765, 0.618316, 0.681713, 0.740319, 0.789176, 0.831279, 0.873132, 0.91255, 0.92833, 0.945552, 0.959786, 0.971815, 0.975448, 0.971097, 0.980313, 0.984648, 0.988947, 0.986328, 0.990765, 0.987174, 0.993049, 0.994162, 0.992768, 0.993483, 0.991127, 0.996967, 0.995381, 0.992867, 0.998447, 0.992687, 0.99793, 0.997831, 1, 0.997283, 0.99654, 0.996815, 1, 0.995349, 1, 1, 1, 0.993421, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
    // Centrality: 50-90 %
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00391389, 0.00688976, 0.0192549, 0.0619815, 0.141893, 0.241798, 0.343937, 0.413673, 0.477422, 0.559137, 0.639819, 0.717685, 0.764019, 0.824119, 0.862554, 0.898145, 0.928348, 0.945772, 0.95873, 0.973406, 0.97358, 0.981122, 0.981965, 0.987635, 0.987374, 0.996862, 0.987764, 0.987788, 0.995098, 0.982363, 0.997881, 0.98977, 1, 0.993789, 0.993127, 0.986842, 1, 0.99422, 0.994152, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1}
  };
  
  // Transform jet pT and centrality to bin indices
  Int_t jetPtBin = jetPt / 5;
  Int_t centralityBin = 0;
  if(centrality > fCard->Get("CentralityBinEdges",2)) centralityBin++;
  if(centrality > fCard->Get("CentralityBinEdges",3)) centralityBin++;
  if(centrality > fCard->Get("CentralityBinEdges",4)) centralityBin++;
  
  // Read the trigger efficiency weight from the table
  return 1.0/triggerEfficiencyTable[centralityBin][jetPtBin];
  
}

/*
 * Get a smearing factor corresponding to worsening the smearing resolution in MC by 20 %
 * This is obtained by multiplying the MC smearing resolution by 0.666 and using this as additional
 * smearing for the data. Smearing factor depends on jet pT and centrality.
 *
 *  Arguments:
 *   Double_t jetPt = Jet pT
 *   const Double_t centrality = Centrality of the event
 *
 *  return: Additional smearing factor
 */
Double_t DijetAnalyzer::GetSmearingFactor(Double_t jetPt, const Double_t centrality) {
  
  // For all the jets above 500 GeV, use the resolution for 500 GeV jet
  if(jetPt > 500) jetPt = 500;
  
  // Find the correct centrality bin
  Int_t centralityBin = 0;
  if(centrality > fCard->Get("CentralityBinEdges",2)) centralityBin++;
  if(centrality > fCard->Get("CentralityBinEdges",3)) centralityBin++;
  if(centrality > fCard->Get("CentralityBinEdges",4)) centralityBin++;
  
  // Set the parameters to the smearing function. pp and PbPb have different smearing function parameters
  if(fDataType == ForestReader::kPp || fDataType == ForestReader::kPpMC){
    // Settings for pp
    fSmearingFunction->SetParameters(0.174881, -0.00091979, 3.50064e-06, -6.52541e-09, 4.64199e-12);
    
  } else {
    
    // Parameters for the smearing function
    Double_t resolutionFit[4][5] = {
      {0.451855, -0.00331992, 1.25897e-05, -2.26434e-08, 1.55081e-11},
      {0.366326, -0.00266997, 1.04733e-05, -1.95302e-08, 1.38409e-11},
      {0.268453, -0.00184878, 7.45201e-06, -1.43486e-08, 1.04726e-11},
      {0.202255, -0.00114677, 4.2566e-06, -7.69286e-09, 5.32617e-12}
    };
    
    for(int iParameter = 0; iParameter < 5; iParameter++){
      // Settings for PbPb
      
      fSmearingFunction->SetParameter(iParameter, resolutionFit[centralityBin][iParameter]);
    }
  }
  
  // After the smearing function is set, read the value to return
  return fSmearingFunction->Eval(jetPt)*0.666;
  
}

/*
 * Get the proper leading jet pT weighting depending on analyzed system
 *
 *  Arguments:
 *   const Double_t jetPt = Jet pT for the weighted jet
 *
 *   return: Multiplicative correction factor for the leading jet pT
 */
Double_t DijetAnalyzer::GetDijetWeight(const Double_t jetPt) const{
  if(fDataType == ForestReader::kPbPb || fDataType == ForestReader::kPp) return 1.0;  // No weight for data
  return 1.0; // TODO: Current dijet weight is disabled
  
  // Only weight 2017 and 2018 MC
  if(fReadMode < 2000){
    return 1.0;  // No weighting for the most peripheral centrality bins 2015
  } else {
    return fDijetWeightFunction->Eval(jetPt);  // No weighting for the most peripheral centrality bins 2018
  }
}

/*
 * Get the proper pT hat weighting for MC
 *
 *  Arguments:
 *   const Int_t ptHat = pT hat value in the event
 *
 *   return: Multiplicative correction factor for the given pT hat
 */
Double_t DijetAnalyzer::GetPtHatWeight(const Double_t ptHat) const{
  if(fDataType == ForestReader::kPp || fDataType == ForestReader::kPbPb) return 1; // No correction for real data
  
  // The event are counted for specific bin edges, so those edges will be used for correction
  const Int_t nBins = 9;
  Int_t correctionBinEdges[nBins] = {30, 50, 80, 120, 170, 220, 280, 370, 460};
  
  // Search the correct bin for the given pT hat value
  Int_t currentBin = -1;
  for(Int_t iBin = 0; iBin < nBins; iBin++){
    if(ptHat < correctionBinEdges[iBin]){
      currentBin = iBin;
      break;
    }
  }
  
  // Cross sections for each bin are given in the twiki https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiForest2015
  //                        pT hat =    15       30       50       80       120      170      220      280      370     460
  Double_t crossSections[nBins+1] = {5.335e-1,3.378e-2,3.778e-3,4.423e-4,6.147e-5,1.018e-5,2.477e-6,6.160e-7,1.088e-7,2.537e-8}; // PYTHIA6 tune Z2
  
  // Different cross sections for PYTHIA8
  if(fReadMode == 1 || fReadMode == 2){
    //                               pT hat =    15       30       50       80       120      170      220      280      370     460
    Double_t pythia8CrossSections[nBins+1] = {5.269e-1,3.455e-2,4.068e-3,4.959e-4,7.096e-5,1.223e-5,3.031e-6,7.746e-7,1.410e-7,3.216e-8}; // PYTHIA8 tune CUETP8M1
    for(int iCrossSection = 0; iCrossSection < nBins+1; iCrossSection++){
      crossSections[iCrossSection] = pythia8CrossSections[iCrossSection];
    }
  }
  
  // Number of events for different pT hat bins in the high forest files
  //  pT hat =            15 30     50     80     120    170    220    280    370   460
  //Int_t ppMcEvents[nBins] = {0,444104,322347,383263,468748,447937,259209,234447,39275};  // File list ppMC_Pythia6_forest_5TeV.txt
  Int_t ppMcEvents[nBins] = {0,  0   ,   0  ,480628,481952,449945,259726,234589,39297};  // File list ppMC_Pythia6_wtaForest.txt
  
  // Different number of events for PYTHIA8
  if(fReadMode == 1){
    //  pT hat =               15 30 50 80    120     170    220    280    370     460
    Int_t pythia8Events[nBins] = {0,0,0,1704437,1063981,833592,953778,1083837,183494};  // File list officialPythia8Forest5TeV.txt
    for(int iPtHatBin = 0; iPtHatBin < nBins; iPtHatBin++){
      ppMcEvents[iPtHatBin] = pythia8Events[iPtHatBin];
    }
  }
  
  // Different number of events for PYTHIA8 made by Dhanush
  if(fReadMode == 2){
    //  pT hat =               15 30 50    80     120    170    220    280   370  460
    Int_t pythia8Events2[nBins] = {0,0,175666,305545,258399,189794,196579,54724,9359};  // File list Pythia8_ak4Calo_5TeV.txt
    for(int iPtHatBin = 0; iPtHatBin < nBins; iPtHatBin++){
      ppMcEvents[iPtHatBin] = pythia8Events2[iPtHatBin];
    }
  }
  
  // Event numbers change a bit in skims, since pT files for pT hat bin 30 are cut out. These numbers are good for list mergedSkimPpPythia5TeV.txt
  if(fForestType == kSkimForest){
    //  pT hat =             15 30 50   80     120     170   220    280    370   460
    Int_t skimEvents[nBins] = {0,0,272976,377670,467966,447818,259188,234443,39272};
    for(Int_t i = 0; i < nBins; i++){
      ppMcEvents[i] = skimEvents[i];
    }
  }
  
  // Xiao's skims with WTA axis start from ptHat 80. These numbers are good for list pythia6FileListWTA.txt
  if(fJetAxis == 2 || fReadMode == 3){
    //  pT hat =             15 30 50   80     120     170   220    280    370   460
    Int_t wtaEvents[nBins] = {0,0,   0  ,430463,475122,448889,259493,234524,39285};
    for(Int_t i = 0; i < nBins; i++){
      ppMcEvents[i] = wtaEvents[i];
    }
  }
  
  // Return the weight for pp
  if(fDataType == ForestReader::kPpMC || fDataType == ForestReader::kLocalTest) {
    if(ppMcEvents[currentBin] == 0){ // This should never happen
      cout << "WARNING! No events found for pT hat bin " << currentBin << endl;
      return 0;
    }
    return (crossSections[currentBin]-crossSections[currentBin+1])/ppMcEvents[currentBin];
  }
  
  // Number of events for different pT hat bins in the forest file list PbPbMC_5TeVPythia6+Hydjet_forests.txt
  //  pT hat =             15 30 50       80     120     170     220     280    370    460
    Int_t PbPbMcEvents[nBins] = {0,0,1761973,2792672,2880895,2696425,2910642,786932,130273}; // Events in WTA forest
  //Int_t PbPbMcEvents[nBins] = {0,0,1761973,2776457,2875707,2695801,2898124,783568,129714}; // Events excluding pthat100 files
  //Int_t PbPbMcEvents[nBins] = {0,0,1761973,4554185,3889938,2847904,2935134,793592,131390}; // Events including pthat100 files
  //Int_t PbPbMcEvents[nBins] = {0,0,0,2571563,2850815,2680567,2891375,781744,129417}; // Old Kurt events
  
  // Also for PbPb, there is a small change for event numbers in skimmed files. These numbers are for PbPbMC_Pythia6HydjetCymbal_list.txt
  if(fForestType == kSkimForest){
    //  pT hat =            15 30 50 80    120     170     220     280    370    460
    Int_t skimEvents[nBins] = {0,0,0,2564435,2846560,2671428,2882855,779493,128987};
    for(Int_t i = 0; i < nBins; i++){
      PbPbMcEvents[i] = skimEvents[i];
    }
  }
  
  // Return the weight for PbPb
  if(PbPbMcEvents[currentBin] == 0) { // This should never happen
    cout << "WARNING! No events found for pT hat bin " << currentBin << endl;
    return 0;
  }
  return (crossSections[currentBin]-crossSections[currentBin+1])/PbPbMcEvents[currentBin];
}

/*
 * Check is a track passes the required subevent cut
 *
 *  Arguments:
 *   const Int_t subeventIndex = Subevent index for the track in consideration
 *
 *  return: true if subevent cut is passes, false if not
 */
Bool_t DijetAnalyzer::PassSubeventCut(const Int_t subeventIndex) const{
  if(fSubeventCut == kSubeventAny) return true;
  if((fSubeventCut == kSubeventZero) && (subeventIndex == 0)) return true;
  if((fSubeventCut == kSubeventNonZero) && (subeventIndex > 0)) return true;
  return false;
}

/*
 * Check if the event passes all the track cuts
 *
 *  Arguments:
 *   ForestReader *eventReader = ForestReader containing the event information checked for event cuts
 *   const Bool_t fillHistograms = Flag for filling the event information histograms.
 *   const Int_t correlationType = Index for correlation type (same or mixed event)
 *
 *   return = True if all event cuts are passed, false otherwise
 */
Bool_t DijetAnalyzer::PassEventCuts(ForestReader *eventReader, const Bool_t fillHistograms, const Int_t correlationType){

  // Cut for primary vertex. Only applied for data.
  if(eventReader->GetPrimaryVertexFilterBit() == 0) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(DijetHistograms::kPrimaryVertex);
  
  // Cut for HB/HE noise. Only applied for data.
  if(eventReader->GetHBHENoiseFilterBit() == 0) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(DijetHistograms::kHBHENoise);
  
  // Cut for collision event selection. Only applied for PbPb data.
  if(eventReader->GetCollisionEventSelectionFilterBit() == 0) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(DijetHistograms::kCollisionEventSelection);
  
  // Cut for beam scraping. Only applied for pp data.
  if(eventReader->GetBeamScrapingFilterBit() == 0) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(DijetHistograms::kBeamScraping);
  
  // Cut for energy deposition in at least 3 hadronic forward towers. Only applied for PbPb data.
  if(eventReader->GetHfCoincidenceFilterBit() == 0) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(DijetHistograms::kHfCoincidence);
  
  // Cut for cluster compatibility. Only applied for PbPb data.
  if(eventReader->GetClusterCompatibilityFilterBit() == 0) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(DijetHistograms::kClusterCompatibility);
  
  // Jet trigger requirement.
  if(fDataType != ForestReader::kPbPb && eventReader->GetCaloJetFilterBit() == 0) return false;  // Not PbPb, use unprescaled trigger
  else if(fDataType == ForestReader::kPbPb && correlationType == DijetHistograms::kSameEvent && eventReader->GetCaloJetFilterBit() == 0) return false; // Regular PbPb, use unprescaled trigger. No trigger requirement for mixed events.
  if(fillHistograms) fHistograms->fhEvents->Fill(DijetHistograms::kCaloJet);
  
  // Cut for vertex z-position
  if(TMath::Abs(eventReader->GetVz()) > fVzCut) return false;
  if(fillHistograms) fHistograms->fhEvents->Fill(DijetHistograms::kVzCut);
  
  return true;
  
}

/*
 * Check if a track passes all the track cuts
 *
 *  Arguments:
 *   const Int_t iTrack = Index of the checked track in reader
 *   TH1F *trackCutHistogram = Histogram to which the track cut performance is filled
 *   const Int_t correlationType = Same or mixed event. Histograms filled only for same event
 *
 *   return: True if all track cuts are passed, false otherwise
 */
Bool_t DijetAnalyzer::PassTrackCuts(const Int_t iTrack, TH1F *trackCutHistogram, const Int_t correlationType){
  
  // Only fill the track cut histograms for same event data
  if(correlationType == DijetHistograms::kSameEvent && fFillTrackHistograms) trackCutHistogram->Fill(DijetHistograms::kAllTracks);
  
  // Cuts specific to generator level MC tracks
  if(fTrackReader[correlationType]->GetTrackCharge(iTrack) == 0) return false;  // Require that the track is charged
  if(correlationType == DijetHistograms::kSameEvent && fFillTrackHistograms) trackCutHistogram->Fill(DijetHistograms::kMcCharge);
  
  if(!PassSubeventCut(fTrackReader[correlationType]->GetTrackSubevent(iTrack))) return false;  // Require desired subevent
  if(correlationType == DijetHistograms::kSameEvent && fFillTrackHistograms) trackCutHistogram->Fill(DijetHistograms::kMcSube);
  
  if(fTrackReader[correlationType]->GetTrackMCStatus(iTrack) != 1) return false;  // Require final state particles
  if(correlationType == DijetHistograms::kSameEvent && fFillTrackHistograms) trackCutHistogram->Fill(DijetHistograms::kMcStatus);
  
  Double_t trackPt = fTrackReader[correlationType]->GetTrackPt(iTrack);
  Double_t trackEta = fTrackReader[correlationType]->GetTrackEta(iTrack);
  Double_t trackEt = (fTrackReader[correlationType]->GetTrackEnergyEcal(iTrack)+fTrackReader[correlationType]->GetTrackEnergyHcal(iTrack))/TMath::CosH(trackEta);
  
  //  ==== Apply cuts for tracks and collect information on how much track are cut in each step ====
  
  // Cut for track pT
  if(trackPt <= fTrackMinPtCut) return false;                     // Minimum pT cut
  if(trackPt >= fJetMaximumPtCut) return false;                   // Maximum pT cut (same as for leading jets)
  if(correlationType == DijetHistograms::kSameEvent && fFillTrackHistograms) trackCutHistogram->Fill(DijetHistograms::kPtCuts);
  
  // Cut for track eta
  if(TMath::Abs(trackEta) >= fTrackEtaCut) return false;          // Eta cut
  if(correlationType == DijetHistograms::kSameEvent && fFillTrackHistograms) trackCutHistogram->Fill(DijetHistograms::kEtaCut);
  
  // New cut for 2018 data based on track algorithm and MVA
  //if(fTrackReader[correlationType]->GetTrackOriginalAlgorithm(iTrack) == 14) return false; // Test a cut from Matt TODO TODO TEST TEST
  if(fTrackReader[correlationType]->GetTrackAlgorithm(iTrack) == 6 && fTrackReader[correlationType]->GetTrackMVA(iTrack) < 0.98 && fReadMode > 2017) return false; // Only apply this cut to 2018 PbPb
  if(correlationType == DijetHistograms::kSameEvent && fFillTrackHistograms) trackCutHistogram->Fill(DijetHistograms::kTrackAlgorithm);
  
  // Cut for high purity
  if(!fTrackReader[correlationType]->GetTrackHighPurity(iTrack)) return false;     // High purity cut
  if(correlationType == DijetHistograms::kSameEvent && fFillTrackHistograms) trackCutHistogram->Fill(DijetHistograms::kHighPurity);
  
  // Cut for relative error for track pT
  if(fTrackReader[correlationType]->GetTrackPtError(iTrack)/trackPt >= fMaxTrackPtRelativeError) return false; // Cut for track pT relative error
  if(correlationType == DijetHistograms::kSameEvent && fFillTrackHistograms) trackCutHistogram->Fill(DijetHistograms::kPtError);
  
  // Cut for track distance from primary vertex
  if(TMath::Abs(fTrackReader[correlationType]->GetTrackVertexDistanceZ(iTrack)/fTrackReader[correlationType]->GetTrackVertexDistanceZError(iTrack)) >= fMaxTrackDistanceToVertex) return false; // Mysterious cut about track proximity to vertex in z-direction
  if(TMath::Abs(fTrackReader[correlationType]->GetTrackVertexDistanceXY(iTrack)/fTrackReader[correlationType]->GetTrackVertexDistanceXYError(iTrack)) >= fMaxTrackDistanceToVertex) return false; // Mysterious cut about track proximity to vertex in xy-direction
  if(correlationType == DijetHistograms::kSameEvent && fFillTrackHistograms) trackCutHistogram->Fill(DijetHistograms::kVertexDistance);
  
  // Cut for energy deposition in calorimeters for high pT tracks
  if(!(trackPt < fCalorimeterSignalLimitPt || (trackEt >= fHighPtEtFraction*trackPt))) return false;  // For high pT tracks, require signal also in calorimeters
  if(correlationType == DijetHistograms::kSameEvent && fFillTrackHistograms) trackCutHistogram->Fill(DijetHistograms::kCaloSignal);
  
  // Cuts for track reconstruction quality
  // TODO: Track layer hit number is missing from MC forests. Estimate it from NHits in this case. This should be removed after MC forests are rerun
  if(fReadMode > 2000 && fDataType == ForestReader::kPbPbMC){
    if( fTrackReader[correlationType]->GetTrackChi2(iTrack) / (1.0*fTrackReader[correlationType]->GetNTrackDegreesOfFreedom(iTrack)) / (fTrackReader[correlationType]->GetNHitsTrack(iTrack) / 1.43) >= fChi2QualityCut) return false; // Track reconstruction quality cut
  } else {
    if( fTrackReader[correlationType]->GetTrackChi2(iTrack) / (1.0*fTrackReader[correlationType]->GetNTrackDegreesOfFreedom(iTrack)) / (1.0*fTrackReader[correlationType]->GetNHitsTrackerLayer(iTrack)) >= fChi2QualityCut) return false; // Track reconstruction quality cut
  }
  if(fTrackReader[correlationType]->GetNHitsTrack(iTrack) < fMinimumTrackHits) return false; // Cut for minimum number of hits per track
  if(correlationType == DijetHistograms::kSameEvent && fFillTrackHistograms) trackCutHistogram->Fill(DijetHistograms::kReconstructionQuality);
  
  // If passed all checks, return true
  return true;
}

/*
 * Get the track efficiency correction for a given track
 *
 *  Arguments:
 *   const Int_t correlationType = Same or mixed event. Histograms filled only for same event
 *   const Int_t iTrack = Index of the track for which the efficiency correction is obtained
 *
 *   return: Multiplicative track efficiency correction
 */
Double_t DijetAnalyzer::GetTrackEfficiencyCorrection(const Int_t correlationType, const Int_t iTrack){
  
  // No correction for generator level tracks
  if(fMcCorrelationType == kRecoGen || fMcCorrelationType == kGenGen) return 1;
  
  // Calculate minimum distance of a track to closest jet. This is needed for track efficiency correction
  Float_t trackRMin = 999;   // Initialize the minimum distance to a jet to some very large value
  Float_t trackR = 999;      // Distance of the track to current jet in the loop
  Float_t trackPt = fTrackReader[correlationType]->GetTrackPt(iTrack);    // Track pT
  Float_t trackEta = fTrackReader[correlationType]->GetTrackEta(iTrack);  // Track eta
  Float_t trackPhi = fTrackReader[correlationType]->GetTrackPhi(iTrack);  // Track phi
  Int_t hiBin = fTrackReader[correlationType]->GetHiBin();                // hiBin for 2018 track correction
  
  // Weight factor only for 2017 pp MC as instructed be the tracking group
  double preWeight = 1.0;
  if(fReadMode > 2000 && fDataType == ForestReader::kPpMC) preWeight = 0.979;
  
  // For PbPb2018 and pp2017, there is an efficiency table from which the correction comes
  if(fReadMode > 2000) return preWeight * fTrackEfficiencyCorrector2018->getCorrection(trackPt, trackEta, hiBin);
  
  // The rest gives a proper correction for 2015 data
  
  // For heavy ion correction, track RMin is not used
  if(fDataType != ForestReader::kPbPb && fDataType != ForestReader::kPbPbMC){
    
    // Need distance to nearest jet for the efficiency correction
    for(Int_t iJet = 0; iJet < fTrackReader[correlationType]->GetNJets(); iJet++){
      
      // Require the same jet quality cuts as when searching for dijets
      if(TMath::Abs(fTrackReader[correlationType]->GetJetEta(iJet)) >= fJetSearchEtaCut) continue; // Require jet eta to be in the search range for dijets
      if(fTrackReader[correlationType]->GetJetPt(iJet) <= fSubleadingJetMinPtCut) continue; // Only consider jets that pass the pT cut for subleading jets
      if(fMinimumMaxTrackPtFraction >= fTrackReader[correlationType]->GetJetMaxTrackPt(iJet)/fTrackReader[correlationType]->GetJetRawPt(iJet)) continue; // Cut for jets with only very low pT particles
      if(fMaximumMaxTrackPtFraction <= fTrackReader[correlationType]->GetJetMaxTrackPt(iJet)/fTrackReader[correlationType]->GetJetRawPt(iJet)) continue; // Cut for jets where all the pT is taken by one track
      // if(TMath::Abs(chargedSum[k]/rawpt[k]) < 0.01) continue; // Jet quality cut from readme file. TODO: Check if should be applied
      
      // Note: The ACos(Cos(jetPhi-trackPhi)) structure transforms deltaPhi to interval [0,Pi]
      trackR = TMath::Power(fTrackReader[correlationType]->GetJetEta(iJet)-trackEta,2)+TMath::Power(TMath::ACos(TMath::Cos(fTrackReader[correlationType]->GetJetPhi(iJet)-trackPhi)),2);
      if(trackRMin*trackRMin>trackR) trackRMin=TMath::Sqrt(trackR);
      
    } // Loop for calculating Rmin
    
  } // If for heavy ions
    
  // Find and return the track efficiency correction
  double trackEfficiency = 1;
  trackEfficiency = fTrackCorrection->getTrkCorr(trackPt, trackEta, trackPhi, hiBin, trackRMin);
  
  return trackEfficiency;
  
}

/*
 * Create the mixing pool
 */
void DijetAnalyzer::CreateMixingPool(){
  
  // Print out a debug message
  if(fDebugLevel > 1) cout << "Creating the mixing pool" << endl;
  
  // Start by emptying the pool from possible earlier files used
  for(Int_t iVz = 0; iVz < kMaxMixingVzBins; iVz++){
    for(Int_t iHiBin = 0; iHiBin < kMaxMixingHiBins; iHiBin++){
      fMixingPool[iVz][iHiBin].clear();
    }
  }
  
  // Initialize the mixed event randomizer
  TRandom3 *mixedEventRandomizer = new TRandom3();  // Randomizer for starting point in the mixed event file
  mixedEventRandomizer->SetSeed(0);
  
  // Start reading the file from random position
  Int_t firstMixingEvent = fnEventsInMixingFile*mixedEventRandomizer->Rndm();            // Start mixing from random spot in file
  if(firstMixingEvent == fnEventsInMixingFile) firstMixingEvent--;                       // Move the index to allowed range
  
  // Loop through the file and and collect the events to mixing table
  Double_t mixedEventVz;   // vz in the event
  Int_t mixedEventHiBin;   // Centrality in the event
  Int_t binForVz;          // Bin given to vz
  Int_t binForCentrality;  // Bin given to centrality
  Int_t iCurrentEvent;     // Index of the event
  
  for(Int_t iMixedEvent = 0; iMixedEvent < fnEventsInMixingFile; iMixedEvent++){
    
    // Read the events from file starting from random position
    iCurrentEvent = iMixedEvent + firstMixingEvent;
    if(iCurrentEvent >= fnEventsInMixingFile) iCurrentEvent -= fnEventsInMixingFile;
    
    // Read vz and HiBin for the current event
    fTrackReader[DijetHistograms::kMixedEvent]->GetEvent(iCurrentEvent);
    if(!PassEventCuts(fTrackReader[DijetHistograms::kMixedEvent],false,DijetHistograms::kMixedEvent)) continue;
    mixedEventVz = fTrackReader[DijetHistograms::kMixedEvent]->GetVz();
    mixedEventHiBin = fTrackReader[DijetHistograms::kMixedEvent]->GetHiBin();
    
    // Determine the bins for vz and centrality
    binForVz = FindMixingVzBin(mixedEventVz);
    binForCentrality = FindMixingHiBin(mixedEventHiBin);
    
    // If we are under the required mixing pool depth, add the event to pool
    if(fMixingPool[binForVz][binForCentrality].size() < fMixingPoolDepth) fMixingPool[binForVz][binForCentrality].push_back(iCurrentEvent);

  }
  
  // Delete the randomizer before returning
  delete mixedEventRandomizer;
  
  // Validate that this pool can be used for mixing
  ValidateMixingPool();
  
}

/*
 * Check that there are events in each of the vz and centrality bins in mixing pool
 */
void DijetAnalyzer::ValidateMixingPool(){
  
  Int_t offset;

  // Loop over all the mixing vectors and check that they have events
  for(Int_t iVz = 0; iVz <= fMaximumMixingVz; iVz++){
    for(Int_t iHiBin = 0; iHiBin <= fMaximumMixingHiBin; iHiBin++){
      offset = 1;
      while(fMixingPool[iVz][iHiBin].size() <= 1){  // We must have at least 2 events in each bin
        
        // Print a debug message about the events in the mixing bin
        if(fDebugLevel > 0) {
          cout << "No events in mixing pool for bin vz: " << iVz << " hibin: " << iHiBin << endl;
          cout << "Filling this bin from adjacent bins with distance " << offset << endl;
        }
        
        // TODO: Validation with offset to both vz and hiBin
        
        // While there are no events in the bin, check as many adjacent bins as needed to find some events
        if(iVz - offset >= 0){
          if(fMixingPool[iVz-offset][iHiBin].size() > 0){
            for(Int_t iEvent = 0; iEvent < fMixingPool[iVz-offset][iHiBin].size(); iEvent++){
              fMixingPool[iVz][iHiBin].push_back(fMixingPool[iVz-offset][iHiBin].at(iEvent));
            }
          }
        }
        if(iVz + offset <= fMaximumMixingVz){
          if(fMixingPool[iVz+offset][iHiBin].size() > 0){
            for(Int_t iEvent = 0; iEvent < fMixingPool[iVz+offset][iHiBin].size(); iEvent++){
              if(fMixingPool[iVz][iHiBin].size() < fMixingPoolDepth) fMixingPool[iVz][iHiBin].push_back(fMixingPool[iVz+offset][iHiBin].at(iEvent));
            }
          }
        }
        if(iHiBin - offset >= 0){
          if(fMixingPool[iVz][iHiBin-offset].size() > 0){
            for(Int_t iEvent = 0; iEvent < fMixingPool[iVz][iHiBin-offset].size(); iEvent++){
              if(fMixingPool[iVz][iHiBin].size() < fMixingPoolDepth) fMixingPool[iVz][iHiBin].push_back(fMixingPool[iVz][iHiBin-offset].at(iEvent));
            }
          }
        }
        if(iHiBin + offset <= fMaximumMixingHiBin){
          if(fMixingPool[iVz][iHiBin+offset].size() > 0){
            for(Int_t iEvent = 0; iEvent < fMixingPool[iVz][iHiBin+offset].size(); iEvent++){
              if(fMixingPool[iVz][iHiBin].size() < fMixingPoolDepth) fMixingPool[iVz][iHiBin].push_back(fMixingPool[iVz][iHiBin+offset].at(iEvent));
            }
          }
        }
        
        offset++;
      } // While loop for empty mixing pool size
    } // Centrality loop
  } // vz loop
}

/*
 * In case we are doing mixing without the pool, prepare vectors for vz and hiBin from the mixing file for faster event matching
 */
void DijetAnalyzer::PrepareMixingVectors(){
  
  // Print out debug message
  if(fDebugLevel > 1) cout << "Preparing for poolless mixing" << endl;
  
  // Initialize the mixed event randomizer
  TRandom3 *mixedEventRandomizer = new TRandom3();  // Randomizer for starting point in the mixed event file
  mixedEventRandomizer->SetSeed(0);
  
  // Start reading the file from a random point
  fMixingStartIndex = fnEventsInMixingFile*mixedEventRandomizer->Rndm();  // Start mixing from random spot in file
  if(fMixingStartIndex == fnEventsInMixingFile) fMixingStartIndex--;       // Move the index to allowed range
  
  // Read vz and hiBin from each event in event mixing file to memory.
  // This way we avoid loading different mixed events in a loop several times
  fMixedEventVz.clear();     // Clear the vectors for any possible
  fMixedEventHiBin.clear();  // contents they might have
  for(Int_t iMixedEvent = 0; iMixedEvent < fnEventsInMixingFile; iMixedEvent++){
    fTrackReader[DijetHistograms::kMixedEvent]->GetEvent(iMixedEvent);
    if(PassEventCuts(fTrackReader[DijetHistograms::kMixedEvent],false,DijetHistograms::kMixedEvent)){
      fMixedEventVz.push_back(fTrackReader[DijetHistograms::kMixedEvent]->GetVz());
      fMixedEventHiBin.push_back(fTrackReader[DijetHistograms::kMixedEvent]->GetHiBin());
    } else { // If event cuts not passed, input values such that events will never be mixed with these
      fMixedEventVz.push_back(100);
      fMixedEventHiBin.push_back(1000);
    }
  }
  
  // Delete the randomizer before returning
  delete mixedEventRandomizer;
}

/*
 *  Find a bin for the input vz
 */
Int_t DijetAnalyzer::FindMixingVzBin(const Double_t vz) const{
  if(TMath::Abs(vz-fVzCut) < 1e-3) return ((fVzCut*2)/fMixingVzBinWidth) - 1;  // If we are exactly at the upper limit, put the result to last available bin
  return (vz + fVzCut)/fMixingVzBinWidth;
}

/*
 *  Find a bin for the input hiBin. There are 200 hiBins
 */
Int_t DijetAnalyzer::FindMixingHiBin(const Int_t hiBin) const{
  if(hiBin < 0) return 0; // For pp, return 0 as hiBin
  if(hiBin == 200 && fMixingHiBinWidth > 1) return (hiBin/fMixingHiBinWidth) - 1; // If we are exactly at the upper limit, put the result to last available bin (to avoid the lowest statistics bin to be alone)
  return hiBin/fMixingHiBinWidth;
}

/*
 * Getter for dijet histograms
 */
DijetHistograms* DijetAnalyzer::GetHistograms() const{
  return fHistograms;
}

/*
 * Check if the mixed event is the same as regular event or not
 *
 *  return: True is events are the same, false if not
 */
Bool_t DijetAnalyzer::CheckForSameEvent(const Int_t sameEventIndex, const Int_t mixedEventIndex) const{
  
  // First, check if the events have the same number of tracks
  const Int_t nTracksSameEvent = fTrackReader[DijetHistograms::kSameEvent]->GetNTracks();
  const Int_t nTracksMixedEvent = fTrackReader[DijetHistograms::kMixedEvent]->GetNTracks();
  
  // If there is different number of tracks, events must be different
  if(nTracksSameEvent != nTracksMixedEvent) return false;
  
  // Loop over the tracks and check whether they have the same pT ot not
  Double_t sameEventTrackPt;
  Double_t mixedEventTrackPt;
  for(Int_t iTrack = 0; iTrack < nTracksSameEvent; iTrack++){
    sameEventTrackPt = fTrackReader[DijetHistograms::kSameEvent]->GetTrackPt(iTrack);
    mixedEventTrackPt = fTrackReader[DijetHistograms::kMixedEvent]->GetTrackPt(iTrack);
    
    // If tracks with same index have different pT, the events must be different
    if(TMath::Abs(sameEventTrackPt-mixedEventTrackPt) > 0.001) return false;
  }
  
  // If all the tracks have exactly the same pT, the events must be the same
  cout << "Attention!" << endl;
  cout << "Attention!" << endl;
  cout << "Attention!" << endl;
  cout << "Events " << sameEventIndex << " in the same event file and " << mixedEventIndex << " in the mixed event file are the same!" << endl;
  cout << "Attention!" << endl;
  cout << "Attention!" << endl;
  cout << "Attention!" << endl;
  return true;
}
