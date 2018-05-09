#ifndef DIJETDRAWER_H
#define DIJETDRAWER_H

// Root includes
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>

// Own includes
#include "DijetCard.h"
#include "JDrawer.h"

/*
 * Class for drawing the histograms produced in the dijet analysis
 */
class DijetDrawer {
 
  // First, define the dimensions for histogram arrays
private:
  static const int knCentralityBins = 4;   // Number of centrality bins
  static const int knTrackPtBins = 6;      // Number of track pT bins
  static const int knCorrelationTypes = 3; // Number of correlation type bins (same event/mixed event/corrected same event)
  static const int knDeltaPhiBins = 4;     // Number of delta phi slices (whole phi/near side/away side/between peaks)
  
public:
  
  /*
   * Constructor for DijetDrawer
   */
  DijetDrawer(TFile *inputFile){
    
    // Connect the inputFile to the class and read the collisions system information
    fInputFile = inputFile;
    fCard = new DijetCard(inputFile);
    fCollisionSystem = fCard->GetDataType();
    
    // Set default values for flags and drawing settings
    fDrawEventInformation = false;
    fDrawDijetHistograms = false;
    fDrawLeadingJetHistograms = false;
    fDrawSubleadingJetHistograms = false;
    fDrawAnyJetHistograms = false;
    fDrawTracks = false;
    fDrawUncorrectedTracks = false;
    fDrawTrackLeadingJetCorrelations = false;
    fDrawUncorrectedTrackLeadingJetCorrelations = false;
    fDrawPtWeightedTrackLeadingJetCorrelations = false;
    fDrawTrackSubleadingJetCorrelations = false;
    fDrawUncorrectedTrackSubleadingJetCorrelations = false;
    fDrawPtWeightedTrackSubleadingJetCorrelations = false;
    
    // Draw mixed event histograms for selected jet-track corraletion histograms
    fDrawSameEvent = false;
    fDrawMixedEvent = false;
    fDrawCorrected = false;
    
    // Choose if you want to write the figures to pdf file
    fSaveFigures = false;
    fFigureFormat = "pdf";
    
    // Logarithmic scales for figures for pT distributions
    fLogPt = true;          // pT distributions
    fLogCorrelation = true; // track-jet deltaPhi-deltaEta distributions
    
    // Plotting style for 2D and 3D plots
    fColorPalette = kRainBow;
    fStyle2D = "colz";
    fStyle3D = "surf1";
    
    // Default drawn centrality bins
    fFirstDrawnCentralityBin = 0;
    fLastDrawnCentralityBin = knCentralityBins-1;
    
    // Default binning for centrality
    for(int iCentrality = 0; iCentrality < knCentralityBins + 1; iCentrality++){
      fCentralityBinIndices[iCentrality] = iCentrality+1;
    }
    
    // Default binning for track pT
    for(int iTrackPt = 0; iTrackPt < knTrackPtBins + 1; iTrackPt++){
      fTrackPtBinIndices[iTrackPt] = iTrackPt + 1;
    }
    
    // Default binning for deltaPhi
    for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
      fLowDeltaPhiBinIndices[iDeltaPhi] = iDeltaPhi+1;
      fHighDeltaPhiBinIndices[iDeltaPhi] = iDeltaPhi+2;
    }
  }
  
  /*
   * Destructor for DijetDrawer
   */
  ~DijetDrawer(){
    // Destructor
    delete fCard;
  }
  
  /*
   * Set up centrality bin indices according to provided bin borders
   */
  void SetCentralityBins(double *binBorders){
    SetBinIndices(knCentralityBins,fCentralityBinIndices,binBorders,4);
  }
  
  /*
   * Set up track pT bin indices according to provided bin borders
   */
  void SetTrackPtBins(double *binBorders){
    SetBinIndices(knTrackPtBins,fTrackPtBinIndices,binBorders,0);
  }
  
  /*
   * Set up deltaPhi bin indices according to provided bin borders
   */
  void SetDeltaPhiBins(double *lowBinBorders, double *highBinBorders){
    SetBinIndices(knDeltaPhiBins,fLowDeltaPhiBinIndices,fHighDeltaPhiBinIndices,lowBinBorders,highBinBorders,1);
  }
  
  /*
   * Loads the histograms from the inputfile
   */
  void LoadHistograms(){
    
    // Load the event information histograms
    if(fDrawEventInformation){
      fhVertexZ = (TH1D*) fInputFile->Get("vertexZ");                 // Vertex z position
      fhEvents = (TH1D*) fInputFile->Get("nEvents");                  // Number of events surviving different event cuts
      fhTrackCuts = (TH1D*) fInputFile->Get("trackCuts");             // Number of tracks surviving different track cuts
      fhCentrality = (TH1D*) fInputFile->Get("centrality");           // Centrality in all events
      fhCentralityDijet = (TH1D*) fInputFile->Get("centralityDijet"); // Centrality in dijet events
    }
    
    /* Load leading jet histograms
     *
     * THnSparse for leading jets:
     *
     *   Histogram name        Axis index       Content of axis
     * ----------------------------------------------------------
     *     leadingJet            Axis 0         Leading jet pT
     *     leadingJet            Axis 1         Leading jet phi
     *     leadingJet            Axis 2         Leading jet eta
     *     leadingJet            Axis 3         Dijet asymmetry
     *     leadingJet            Axis 4           Centrality
     */
    if(fDrawLeadingJetHistograms){
      LoadSingleJetHistograms(fhLeadingJetPt,fhLeadingJetPhi,fhLeadingJetEta,fhLeadingJetEtaPhi,"leadingJet",4);
    }
  }
  
  // Setter for drawing event information
  void SetDrawEventInformation(bool drawOrNot){
    fDrawEventInformation = drawOrNot;
  }
  
  // Setter for drawing dijet histograms
  void SetDrawDijetHistograms(bool drawOrNot){
    fDrawDijetHistograms = drawOrNot;
  }
  
  // Setter for drawing leading jet histograms
  void SetDrawLeadingJetHistograms(bool drawOrNot){
    fDrawLeadingJetHistograms = drawOrNot;
  }
  
  // Setter for drawing subleading jet histograms
  void SetDrawSubleadingJetHistograms(bool drawOrNot){
    fDrawSubleadingJetHistograms = drawOrNot;
  }
  
  // Setter for drawing all jet histograms
  void SetDrawAnyJetHistograms(bool drawOrNot){
    fDrawAnyJetHistograms = drawOrNot;
  }
  
  // Setter for drawing jet histograms
  void SetDrawAllJets(bool drawLeading, bool drawSubleading, bool drawAny){
    SetDrawLeadingJetHistograms(drawLeading);
    SetDrawSubleadingJetHistograms(drawSubleading);
    SetDrawAnyJetHistograms(drawAny);
  }
  
  // Setter for drawing tracks
  void SetDrawTracks(bool drawOrNot){
    fDrawTracks = drawOrNot;
  }
  
  // Setter for drawing uncorrected tracks
  void SetDrawTracksUncorrected(bool drawOrNot){
    fDrawUncorrectedTracks = drawOrNot;
  }
  
  // Setter for drawing track histograms
  void SetDrawAllTracks(bool drawTracks, bool drawUncorrected){
    SetDrawTracks(drawTracks);
    SetDrawTracksUncorrected(drawUncorrected);
  }
  
  // Setter for drawing leading jet-track correlations
  void SetDrawTrackLeadingJetCorrelations(bool drawOrNot){
    fDrawTrackLeadingJetCorrelations = drawOrNot;
  }
  
  // Setter for drawing uncorrected leading jet-track correlations
  void SetDrawTrackLeadingJetCorrelationsUncorrected(bool drawOrNot){
    fDrawUncorrectedTrackLeadingJetCorrelations = drawOrNot;
  }
  
  // Setter for drawing pT weighted leading jet-track correlations
  void SetDrawTrackLeadingJetCorrelationsPtWeighted(bool drawOrNot){
    fDrawPtWeightedTrackLeadingJetCorrelations = drawOrNot;
  }
  
  // Setter for drawing all correlations related to tracks and leading jets
  void SetDrawAllTrackLeadingJetCorrelations(bool drawLeading, bool drawUncorrected, bool drawPtWeighted){
    SetDrawTrackLeadingJetCorrelations(drawLeading);
    SetDrawTrackLeadingJetCorrelationsUncorrected(drawUncorrected);
    SetDrawTrackLeadingJetCorrelationsPtWeighted(drawPtWeighted);
  }
  
  // Setter for drawing subleading jet-track correlations
  void SetDrawTrackSubleadingJetCorrelations(bool drawOrNot){
    fDrawTrackSubleadingJetCorrelations = drawOrNot;
  }
  
  // Setter for drawing uncorrected subleading jet-track correlations
  void SetDrawTrackSubleadingJetCorrelationsUncorrected(bool drawOrNot){
    fDrawUncorrectedTrackSubleadingJetCorrelations = drawOrNot;
  }
  
  // Setter for drawing pT weighted subleading jet-track correlations
  void SetDrawTrackSubleadingJetCorrelationsPtWeighted(bool drawOrNot){
    fDrawPtWeightedTrackSubleadingJetCorrelations = drawOrNot;
  }
  
  // Setter for drawing all correlations related to tracks and subleading jets
  void SetDrawAllTrackSubleadingJetCorrelations(bool drawSubleading, bool drawUncorrected, bool drawPtWeighted){
    SetDrawTrackSubleadingJetCorrelations(drawSubleading);
    SetDrawTrackSubleadingJetCorrelationsUncorrected(drawUncorrected);
    SetDrawTrackSubleadingJetCorrelationsPtWeighted(drawPtWeighted);
  }
  
  // Setter for drawing same event correlation distributions
  void SetDrawSameEvent(bool drawOrNot){
    fDrawSameEvent = drawOrNot;
  }
  
  // Setter for drawing mixed event correlation distributions
  void SetDrawMixedEvent(bool drawOrNot){
    fDrawMixedEvent = drawOrNot;
  }
  
  // Setter for drawing corrected correlation distributions
  void SetDrawCorrectedCorrelations(bool drawOrNot){
    fDrawCorrected = drawOrNot;
  }
  
  // Setter for drawing different correlation types
  void SetDrawCorrelationTypes(bool sameEvent, bool mixedEvent, bool corrected){
    SetDrawSameEvent(sameEvent);
    SetDrawMixedEvent(mixedEvent);
    SetDrawCorrectedCorrelations(corrected);
  }
  
  // Setter for saving the figures to a file
  void SetSaveFigures(bool saveOrNot, const char *format = "pdf"){
    fSaveFigures = saveOrNot;
    fFigureFormat = format;
  }
  
  // Setter for logarithmic pT axis
  void SetLogPt(bool isLog){
    fLogPt = isLog;
  }
  
  // Setter for logarithmic z axis for correlation plots
  void SetLogCorrelation(bool isLog){
    fLogCorrelation = isLog;
  }
  
  // Setter for logarithmix axes
  void SetLogAxes(bool pt, bool correlation){
    SetLogPt(pt);
    SetLogCorrelation(correlation);
  }

  // Setter for color palette
  void SetColorPalette(int color){
    fColorPalette = color;
  }
  
  // Setter for 2D drawing style
  void SetDrawingStyle2D(const char* style){
    fStyle2D = style;
  }
  
  // Setter for 3D drawing style
  void SetDrawingStyle3D(const char* style){
    fStyle3D = style;
  }
  
  // Setter for 2D drawing style
  void SetDrawingStyles(int color, const char* style2D, const char* style3D){
    SetColorPalette(color);
    SetDrawingStyle2D(style2D);
    SetDrawingStyle3D(style3D);
  }
  
  // Setter for drawn centrality bins
  void SetCentralityBinRange(int first, int last){
    fFirstDrawnCentralityBin = first;
    fLastDrawnCentralityBin = last;
    
    // Sanity check for drawn centrality bins
    if(fFirstDrawnCentralityBin < 0) fFirstDrawnCentralityBin = 0;
    if(fLastDrawnCentralityBin < fFirstDrawnCentralityBin) fLastDrawnCentralityBin = fFirstDrawnCentralityBin;
    if(fLastDrawnCentralityBin > knCentralityBins-1) fLastDrawnCentralityBin = knCentralityBins-1;
  }
  
private:
  
  // Data members
  TFile *fInputFile;  // File from which the histograms are read
  DijetCard *fCard;   // Card inside the data file for binning, cut collision system etc. information
  TString fCollisionSystem;  // String for collision system (pp,PbPb,pp MC,PbPb MC,localTest)
  
  // ==============================================
  // ======== Flags for histograms to draw ========
  // ==============================================
  
  bool fDrawEventInformation;
  bool fDrawDijetHistograms;
  bool fDrawLeadingJetHistograms;
  bool fDrawSubleadingJetHistograms;
  bool fDrawAnyJetHistograms;
  bool fDrawTracks;
  bool fDrawUncorrectedTracks;
  bool fDrawTrackLeadingJetCorrelations;
  bool fDrawUncorrectedTrackLeadingJetCorrelations;
  bool fDrawPtWeightedTrackLeadingJetCorrelations;
  bool fDrawTrackSubleadingJetCorrelations;
  bool fDrawUncorrectedTrackSubleadingJetCorrelations;
  bool fDrawPtWeightedTrackSubleadingJetCorrelations;
  
  // ==============================================
  // ============== Drawing settings ==============
  // ==============================================
  
  // Draw mixed event histograms for selected jet-track corraletion histograms
  bool fDrawSameEvent;
  bool fDrawMixedEvent;
  bool fDrawCorrected;
  
  // Choose if you want to write the figures to pdf file
  bool fSaveFigures;
  const char* fFigureFormat;
  
  // Logarithmic scales for figures for pT distributions
  bool fLogPt;          // pT distributions
  bool fLogCorrelation; // track-jet deltaPhi-deltaEta distributions
  
  // Plotting style for 2D and 3D plots
  int fColorPalette;
  const char* fStyle2D;
  const char* fStyle3D;
  
  // Drawn centrality bins
  int fFirstDrawnCentralityBin;
  int fLastDrawnCentralityBin;
  
  // =============================================
  // ============ Binning information ============
  // =============================================
  int fCentralityBinIndices[knCentralityBins+1];
  int fTrackPtBinIndices[knTrackPtBins+1];
  int fLowDeltaPhiBinIndices[knDeltaPhiBins];
  int fHighDeltaPhiBinIndices[knDeltaPhiBins];
  
  // =============================================
  // ===== Histograms for the dijet analysis =====
  // =============================================
  
  // Vertex z position
  TH1D *fhVertexZ;
  
  // Number of events surviving different event cuts
  TH1D *fhEvents;
  
  // Number of tracks surviving different track cuts
  TH1D *fhTrackCuts;
  
  // Centrality of all and dijet events
  TH1D *fhCentrality;
  TH1D *fhCentralityDijet;
  
  // Histograms for leading jets
  TH1D *fhLeadingJetPt[knCentralityBins];               // Leading jet pT histograms
  TH1D *fhLeadingJetPhi[knCentralityBins];              // Leading jet phi histograms
  TH1D *fhLeadingJetEta[knCentralityBins];              // Leading jet eta histograms
  TH2D *fhLeadingJetEtaPhi[knCentralityBins];           // 2D eta-phi histogram for leading jet
  
  // Histograms for subleading jets
  TH1D *fhSubleadingJetPt[knCentralityBins];            // Subleading jet pT histograms
  TH1D *fhSubleadingJetPhi[knCentralityBins];           // Subleading jet phi histograms
  TH1D *fhSubleadingJetEta[knCentralityBins];           // Subleading jet eta histograms
  TH2D *fhSubleadingJetEtaPhi[knCentralityBins];        // 2D eta-phi histogram for subleading jet
  
  // Histograms for all jets
  TH1D *fhAnyJetPt[knCentralityBins] ;                  // Any jet pT histograms
  TH1D *fhAnyJetPhi[knCentralityBins];                  // Any jet phi histograms
  TH1D *fhAnyJetEta[knCentralityBins];                  // Any jet eta histograms
  TH2D *fhAnyJetEtaPhi[knCentralityBins];               // 2D eta-phi histogram for all jets
  
  // Histograms for dijets
  TH1D *fhDijetDphi[knCentralityBins];                  // Dijet deltaPhi histograms
  TH1D *fhDijetAsymmetry[knCentralityBins];             // Dijet asymmetry histograms
  TH2D *fhDijetLeadingVsSubleadingPt[knCentralityBins]; // Leading versus subleading jet pT 2D histograms
  
  // Histograms for tracks in dijet events
  TH1D *fhTrackPt[knCorrelationTypes][knCentralityBins] ;                   // Track pT histograms
  TH1D *fhTrackPhi[knCorrelationTypes][knCentralityBins];                   // Track phi histograms
  TH1D *fhTrackEta[knCorrelationTypes][knCentralityBins];                   // Track eta histograms
  TH2D *fhTrackEtaPhi[knCorrelationTypes][knCentralityBins];                // 2D eta-phi histogram for track
  
  // Histograms for uncorrected tracks in dijet events
  TH1D *fhTrackPtUncorrected[knCorrelationTypes][knCentralityBins] ;        // Uncorrected track pT histograms
  TH1D *fhTrackPhiUncorrected[knCorrelationTypes][knCentralityBins];        // Uncorrected track phi histograms
  TH1D *fhTrackEtaUncorrected[knCorrelationTypes][knCentralityBins];        // Uncorrected track eta histograms
  TH2D *fhTrackEtaPhiUncorrected[knCorrelationTypes][knCentralityBins];     // 2D eta-phi histogram for uncorrected tracks
  
  // Histograms for track-leading jet correlations
  TH1D *fhTrackLeadingJetDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins];                // DeltaPhi between track and leading jet
  TH1D *fhTrackLeadingJetDeltaEta[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins]; // DeltaEta between track and leading jet
  TH2D *fhTrackLeadingJetDeltaEtaDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins];        // DeltaEta and deltaPhi between track and leading jet
  
  // Histograms for uncorrected track-leading jet correlations
  TH1D *fhTrackLeadingJetDeltaPhiUncorrected[knCorrelationTypes][knCentralityBins][knTrackPtBins];                // Uncorrected deltaPhi between track and leading jet
  TH1D *fhTrackLeadingJetDeltaEtaUncorrected[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins]; // Uncorrected deltaEta between track and leading jet
  TH2D *fhTrackLeadingJetDeltaEtaDeltaPhiUncorrected[knCorrelationTypes][knCentralityBins][knTrackPtBins];        // Uncorrected deltaEta and deltaPhi between track and leading jet
  
  // Histograms for pT weighted track-leading jet correlations
  TH1D *fhTrackLeadingJetDeltaPhiPtWeighted[knCorrelationTypes][knCentralityBins][knTrackPtBins];                // pT weighted deltaPhi between track and leading jet
  TH1D *fhTrackLeadingJetDeltaEtaPtWeighted[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins]; // pT weighted deltaEta between track and leading jet
  TH2D *fhTrackLeadingJetDeltaEtaDeltaPhiPtWeighted[knCorrelationTypes][knCentralityBins][knTrackPtBins];        // pT weighted deltaEta and deltaPhi between track and leading jet
  
  // Histograms for track-subleading jet correlations
  TH1D *fhTrackSubleadingJetDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins];                // DeltaPhi between track and subleading jet
  TH1D *fhTrackSubleadingJetDeltaEta[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins]; // DeltaEta between track and subleading jet
  TH2D *fhTrackSubleadingJetDeltaEtaDeltaPhi[knCorrelationTypes][knCentralityBins][knTrackPtBins];        // DeltaEta and deltaPhi between track and subleading jet
  
  // Histograms for uncorrected track-subleading jet correlations
  TH1D *fhTrackSubleadingJetDeltaPhiUncorrected[knCorrelationTypes][knCentralityBins][knTrackPtBins];                // Uncorrected deltaPhi between track and subleading jet
  TH1D *fhTrackSubleadingJetDeltaEtaUncorrected[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins]; // Uncorrected deltaEta between track and subleading jet
  TH2D *fhTrackSubleadingJetDeltaEtaDeltaPhiUncorrected[knCorrelationTypes][knCentralityBins][knTrackPtBins];        // Uncorrected deltaEta and deltaPhi between track and subleading jet
  
  // Histograms for pT weighted track-subleading jet correlations
  TH1D *fhTrackSubleadingJetDeltaPhiPtWeighted[knCorrelationTypes][knCentralityBins][knTrackPtBins];                // pT weighted deltaPhi between track and subleading jet
  TH1D *fhTrackSubleadingJetDeltaEtaPtWeighted[knCorrelationTypes][knCentralityBins][knTrackPtBins][knDeltaPhiBins]; // pT weighted deltaEta between track and subleading jet
  TH2D *fhTrackSubleadingJetDeltaEtaDeltaPhiPtWeighted[knCorrelationTypes][knCentralityBins][knTrackPtBins];        // pT weighted deltaEta and deltaPhi between track and subleading jet
  
  /*
   * Read the bin indices for given bin borders
   *
   *  Arguments:
   *   const int nBins = Number of bins for the indices
   *   int *binIndices = Array of integers to be filled with bin index information read from the file
   *   const double *binBorders = Array for bin borders that are searched from the file
   *   const int iAxis = Index of the axis used for reading bin indices
   */
  void SetBinIndices(const int nBins, int *binIndices, const double *binBorders, const int iAxis){
    TH1D* hBinner = FindHistogram(fInputFile,"trackLeadingJet",iAxis,0,0,0);
    for(int iBin = 0; iBin < nBins+1; iBin++){
      binIndices[iBin] = hBinner->GetXaxis()->FindBin(binBorders[iBin]);
    }
  }
  
  /*
   * Read the bin indices for given bin borders
   *
   *  Arguments:
   *   const int nBins = Number of bins for the indices
   *   int *lowBinIndices = Array of integers to be filled with bin low edge index information read from the file
   *   int *highBinIndices = Array of integers to be filled with bin high edge index information read from the file
   *   const double *lowBinBorders = Array for low bin borders that are searched from the file
   *   const double *highBinBorders = Array for high bin borders that are searched from the file
   *   const int iAxis = Index of the axis used for reading bin indices
   */
  void SetBinIndices(const int nBins, int *lowBinIndices, int *highBinIndices, const double *lowBinBorders, const double *highBinBorders, const int iAxis){
    TH1D* hBinner = FindHistogram(fInputFile,"trackLeadingJet",iAxis,0,0,0);
    for(int iBin = 0; iBin < nBins; iBin++){
      lowBinIndices[iBin] = hBinner->GetXaxis()->FindBin(lowBinBorders[iBin]);
      highBinIndices[iBin] = hBinner->GetXaxis()->FindBin(highBinBorders[iBin]);
    }
  }
  
  /*
   * Extract a 2D histogram from a given centrality bin from THnSparseD
   *
   *  Arguments:
   *   TFile *inputFile = Inputfile containing the THnSparse to be read
   *   const char *name = Name of the THnSparse that is read
   *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
   *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
   *   int nAxes = Number of axes that are restained for the projection
   *   int *restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
   *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
   *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
   */
  TH2D* FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex){
    THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
    for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
    char newName[100];
    sprintf(newName,"%s%d",histogramArray->GetName(),lowBinIndex[0]);
    TH2D *projectedHistogram = (TH2D*) histogramArray->Projection(yAxis,xAxis);
    projectedHistogram->SetName(newName);
    return projectedHistogram;
  }
  
  /*
   * Extract a 2D histogram from a given centrality bin from THnSparseD
   *
   *  Arguments:
   *   TFile *inputFile = Inputfile containing the THnSparse to be read
   *   const char *name = Name of the THnSparse that is read
   *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
   *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
   *   int restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
   *   int lowBinIndex = Index of the lowest considered bin in the restriction axis
   *   int highBinIndex = Index of the highest considered bin in the restriction axis
   *   int restrictionAxis2 = Index for the axis in THnSparse that is used as a second restriction for the projection
   *   int lowBinIndex2 = Index of the lowest considered bin in the second restriction axis
   *   int highBinIndex2 = Index of the highest considered bin in the second restriction axis
   */
  TH2D* FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0){
    THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
    histogramArray->GetAxis(restrictionAxis)->SetRange(lowBinIndex,highBinIndex);
    if(highBinIndex2 > 0 && lowBinIndex2 > 0) histogramArray->GetAxis(restrictionAxis2)->SetRange(lowBinIndex2,highBinIndex2);
    char newName[100];
    sprintf(newName,"%s%d",histogramArray->GetName(),lowBinIndex);
    TH2D *projectedHistogram = (TH2D*) histogramArray->Projection(yAxis,xAxis);
    projectedHistogram->SetName(newName);
    return projectedHistogram;
  }
  
  /*
   * Extract a 2D histogram from a given centrality bin from THnSparseD
   *
   *  Arguments:
   *   TFile *inputFile = Inputfile containing the THnSparse to be read
   *   const char *name = Name of the THnSparse that is read
   *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH1D
   *   int nAxes = Number of axes that are restained for the projection
   *   int *restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
   *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
   *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
   */
  TH1D* FindHistogram(TFile *inputFile, const char *name, int xAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex){
    THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
    for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
    char newName[100];
    sprintf(newName,"%s%d",histogramArray->GetName(),lowBinIndex[0]);
    TH1D *projectedHistogram = (TH1D*) histogramArray->Projection(xAxis);
    projectedHistogram->SetName(newName);
    return projectedHistogram;
  }
  
  /*
   * Extract a histogram from a given centrality bin from THnSparseD
   *
   *  Arguments:
   *   TFile *inputFile = Inputfile containing the THnSparse to be read
   *   const char *name = Name of the THnSparse that is read
   *   int xAxis = Index for the axis in THnSparse that is projected to TH1D
   *   int restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
   *   int lowBinIndex = Index of the lowest considered bin in the restriction axis
   *   int highBinIndex = Index of the highest considered bin in the restriction axis
   *   int restrictionAxis2 = Index for the axis in THnSparse that is used as a second restriction for the projection
   *   int lowBinIndex2 = Index of the lowest considered bin in the second restriction axis
   *   int highBinIndex2 = Index of the highest considered bin in the second restriction axis
   */
  TH1D* FindHistogram(TFile *inputFile, const char *name, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0){
    THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
    histogramArray->GetAxis(restrictionAxis)->SetRange(lowBinIndex,highBinIndex);
    if(highBinIndex2 > 0 && lowBinIndex2 > 0) histogramArray->GetAxis(restrictionAxis2)->SetRange(lowBinIndex2,highBinIndex2);
    char newName[100];
    sprintf(newName,"%s%d",histogramArray->GetName(),lowBinIndex);
    TH1D *projectedHistogram = (TH1D*) histogramArray->Projection(xAxis);
    projectedHistogram->SetName(newName);
    return projectedHistogram;
  }
  
  /*
   * Loader for single jet histograms
   *
   *  Arguments:
   *    TH1D *hJetPt[knCentralityBins] = Array of jet pT histograms
   *    TH1D *hJetPhi[knCentralityBins] = Array of jet phi histograms
   *    TH1D *hJetEta[knCentralityBins] = Array of jet eta histograms
   *    TH2D *hJetEtaPhi[knCentralityBins] = Array of jet eta-phi histograms
   *    const char* name = Name of the histogram in the input file
   *    const int iCentralityAxis = Index of centrality axis in THnSparse
   *
   * THnSparse for single jets:
   *
   *   Histogram name: leadingJet/subleadingJet/anyJet
   *
   *     Axis index       Content of axis         Exception
   * ----------------------------------------------------------
   *       Axis 0         Leading jet pT
   *       Axis 1         Leading jet phi
   *       Axis 2         Leading jet eta
   *       Axis 3         Dijet asymmetry    (for anyJet: Centrality)
   *       Axis 4           Centrality       (for anyJet: Nothing)
   */
   */
  void LoadSingleJetHistograms(TH1D *hJetPt[knCentralityBins], TH1D* hJetPhi[knCentralityBins], TH1D* hJetEta[knCentralityBins], TH2D* hJetEtaPhi[knCentralityBins], const char* name, const int iCentralityAxis){
    
    // Define helper variables
    int duplicateRemoverCentrality = -1;
    int lowerCentralityBin = 0;
    int higherCentralityBin = 0;
    
    for(int iCentralityBin = fFirstDrawnCentralityBin; iCentralityBin <= fLastDrawnCentralityBin; iCentralityBin++){
      
      // Select the bin indices
      if(iCentralityBin == fLastDrawnCentralityBin) duplicateRemoverCentrality = 0;
      lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
      higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
      
      hJetPt[iCentralityBin] = FindHistogram(fInputFile,name,0,iCentralityAxis,lowerCentralityBin,higherCentralityBin);
      hJetPhi[iCentralityBin] = FindHistogram(fInputFile,name,1,iCentralityAxis,lowerCentralityBin,higherCentralityBin);
      hJetEta[iCentralityBin] = FindHistogram(fInputFile,name,2,iCentralityAxis,lowerCentralityBin,higherCentralityBin);
      hJetEtaPhi[iCentralityBin] = FindHistogram2D(fInputFile,name,1,2,iCentralityAxis,lowerCentralityBin,higherCentralityBin);
    }
  }
  
};

#endif
