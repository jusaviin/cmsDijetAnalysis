/*
 * Implementation of DijetHistogramManager
 */

// Root includes
#include <THnSparse.h>
#include <TPad.h>

// Own includes
#include "DijetHistogramManager.h"

/*
 * Constructor
 */
DijetHistogramManager::DijetHistogramManager(TFile *inputFile) :
  fInputFile(inputFile),
  fLoadEventInformation(false),
  fLoadDijetHistograms(false),
  fFirstLoadedCentralityBin(0),
  fLastLoadedCentralityBin(knCentralityBins-1),
  fFirstLoadedTrackPtBin(0),
  fLastLoadedTrackPtBin(knTrackPtBins-1)
{
  // Read card from inputfile and collision system from card
  fCard = new DijetCard(inputFile);
  TString collisionSystem = fCard->GetDataType();
  
  // Make a string for collision system based on information on the card
  fSystemAndEnergy = Form("%s 5.02 TeV",collisionSystem.Data());
  fCompactSystemAndEnergy = fSystemAndEnergy;
  fCompactSystemAndEnergy.ReplaceAll(" ","");
  fCompactSystemAndEnergy.ReplaceAll(".","v");
  
  // Create a new DijetMethods
  fMethods = new DijetMethods();
  
  // Do not draw anything by default
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    fLoadJetTrackCorrelations[iJetTrack] = false;
  }
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    fLoadTracks[iTrackType] = false;
  }
  for(int iJetType = 0; iJetType < knSingleJetCategories; iJetType++){
    fLoadSingleJets[iJetType] = false;
  }
  
  // Default binning for centrality
  for(int iCentrality = 0; iCentrality < knCentralityBins + 1; iCentrality++){
    fCentralityBinIndices[iCentrality] = iCentrality+1;
    fCentralityBinBorders[iCentrality] = 0;
  }
  
  // Default binning for track pT
  for(int iTrackPt = 0; iTrackPt < knTrackPtBins + 1; iTrackPt++){
    fTrackPtBinIndices[iTrackPt] = iTrackPt + 1;
    fFineTrackPtBinIndices[iTrackPt] = iTrackPt + 1;
    fTrackPtBinBorders[iTrackPt] = 0;
  }
  
  // Default binning for deltaPhi
  for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
    fLowDeltaPhiBinIndices[iDeltaPhi] = iDeltaPhi+1;
    fHighDeltaPhiBinIndices[iDeltaPhi] = iDeltaPhi+2;
    fDeltaPhiString[iDeltaPhi] = "";
    fCompactDeltaPhiString[iDeltaPhi] = "";
  }
  
  // Initialize all the other histograms to null
  fhVertexZ = NULL;         // Vertex z position
  fhEvents = NULL;          // Number of events surviving different event cuts
  fhTrackCuts = NULL;       // Number of tracks surviving different track cuts
  fhCentrality = NULL;      // Centrality of all events
  fhCentralityDijet = NULL; // Centrality of dijet events
  
  // Centrality loop
  for(int iCentrality = 0; iCentrality < knCentralityBins; iCentrality++){
    fhDijetDphi[iCentrality] = NULL;                  // Dijet deltaPhi histograms
    fhDijetAsymmetry[iCentrality] = NULL;             // Dijet asymmetry histograms
    fhDijetLeadingVsSubleadingPt[iCentrality] = NULL; // Leading versus subleading jet pT 2D histograms
    
    // Single jet category loop
    for(int iJetCategory = 0; iJetCategory < knSingleJetCategories; iJetCategory++){
      fhJetPt[iJetCategory][iCentrality] = NULL;      // Jet pT histograms
      fhJetPhi[iJetCategory][iCentrality] = NULL;     // Jet phi histograms
      fhJetEta[iJetCategory][iCentrality] = NULL;     // Jet eta histograms
      fhJetEtaPhi[iJetCategory][iCentrality] = NULL;  // 2D eta-phi histogram for jets
    } // Single jet categories loop
    
    // Event correlation type loop
    for(int iCorrelationType = 0; iCorrelationType < knCorrelationTypes; iCorrelationType++){
      
      // Loop over track categories
      for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
        fhTrackPt[iTrackType][iCorrelationType][iCentrality] = NULL;   // Track pT histograms
        
        // Loop over track pT bins
        for(int iTrackPt = 0; iTrackPt < knTrackPtBins + 1; iTrackPt++){
          fhTrackPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt] = NULL;    // Track phi histograms
          fhTrackEta[iTrackType][iCorrelationType][iCentrality][iTrackPt] = NULL;    // Track eta histograms
          fhTrackEtaPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt] = NULL; // 2D eta-phi histogram for track
        } // Track pT loop
        
      } // Track category loop
      
      // Loop over jet-track correlation types
      for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
        
        // Loop over track pT bins
        for(int iTrackPt = 0; iTrackPt < knTrackPtBins; iTrackPt++){
          fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt] = NULL;         // DeltaPhi between jet and track
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt] = NULL; // DeltaEta and deltaPhi between jet and track
          
          // Loop over deltaPhi bins
          for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
            fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iCentrality][iTrackPt][iDeltaPhi] = NULL; // DeltaEta between jet and track
          } // DeltaPhi loop
        } // Track pT loop
      } // Jet-track correlation type loop
    } // Event correlation type loop
    
    // Jet shape histograms
    for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
      for(int iTrackPt = 0; iTrackPt < knTrackPtBins; iTrackPt++){
        for(int iJetShape = 0; iJetShape < knJetShapeTypes; iJetShape++){
          fhJetShape[iJetShape][iJetTrack][iCentrality][iTrackPt] = NULL;
        } // Jet shape type loop
      } // Track pT loop
    } // Jet-track correlation type loop
  } // Centrality loop
}

/*
 * Copy constructor
 */
DijetHistogramManager::DijetHistogramManager(const DijetHistogramManager& in) :
  fInputFile(in.fInputFile),
  fCard(in.fCard),
  fSystemAndEnergy(in.fSystemAndEnergy),
  fCompactSystemAndEnergy(in.fCompactSystemAndEnergy),
  fMethods(in.fMethods),
  fLoadEventInformation(in.fLoadEventInformation),
  fLoadDijetHistograms(in.fLoadDijetHistograms),
  fFirstLoadedCentralityBin(in.fFirstLoadedCentralityBin),
  fLastLoadedCentralityBin(in.fLastLoadedCentralityBin),
  fFirstLoadedTrackPtBin(in.fFirstLoadedTrackPtBin),
  fLastLoadedTrackPtBin(in.fLastLoadedTrackPtBin),
  fhVertexZ(in.fhVertexZ),
  fhEvents(in.fhEvents),
  fhTrackCuts(in.fhTrackCuts),
  fhCentrality(in.fhCentrality),
  fhCentralityDijet(in.fhCentralityDijet)
{
  // Copy constructor
  
  // Copy all values
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    fLoadJetTrackCorrelations[iJetTrack] = in.fLoadJetTrackCorrelations[iJetTrack];
  }
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    fLoadTracks[iTrackType] = in.fLoadTracks[iTrackType];
  }
  for(int iJetType = 0; iJetType < knSingleJetCategories; iJetType++){
    fLoadSingleJets[iJetType] = in.fLoadSingleJets[iJetType];
  }
  
  // Copy binning for centrality
  for(int iCentrality = 0; iCentrality < knCentralityBins + 1; iCentrality++){
    fCentralityBinIndices[iCentrality] = in.fCentralityBinIndices[iCentrality];
    fCentralityBinBorders[iCentrality] = in.fCentralityBinBorders[iCentrality];
  }
  
  // Copy binning for track pT
  for(int iTrackPt = 0; iTrackPt < knTrackPtBins + 1; iTrackPt++){
    fTrackPtBinIndices[iTrackPt] = in.fTrackPtBinIndices[iTrackPt];
    fFineTrackPtBinIndices[iTrackPt] = in.fFineTrackPtBinIndices[iTrackPt];
    fTrackPtBinBorders[iTrackPt] = in.fTrackPtBinBorders[iTrackPt];
  }
  
  // Copy binning for deltaPhi
  for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
    fLowDeltaPhiBinIndices[iDeltaPhi] = in.fLowDeltaPhiBinIndices[iDeltaPhi];
    fHighDeltaPhiBinIndices[iDeltaPhi] = in.fHighDeltaPhiBinIndices[iDeltaPhi];
    fDeltaPhiString[iDeltaPhi] = in.fDeltaPhiString[iDeltaPhi];
    fCompactDeltaPhiString[iDeltaPhi] = in.fCompactDeltaPhiString[iDeltaPhi];
  }
  
  // Centrality loop
  for(int iCentrality = 0; iCentrality < knCentralityBins; iCentrality++){
    fhDijetDphi[iCentrality] = in.fhDijetDphi[iCentrality];                                   // Dijet deltaPhi histograms
    fhDijetAsymmetry[iCentrality] = in.fhDijetAsymmetry[iCentrality];                         // Dijet asymmetry histograms
    fhDijetLeadingVsSubleadingPt[iCentrality] = in.fhDijetLeadingVsSubleadingPt[iCentrality]; // Leading versus subleading jet pT 2D histograms
    
    // Single jet category loop
    for(int iJetCategory = 0; iJetCategory < knSingleJetCategories; iJetCategory++){
      fhJetPt[iJetCategory][iCentrality] = in.fhJetPt[iJetCategory][iCentrality];         // Jet pT histograms
      fhJetPhi[iJetCategory][iCentrality] = in.fhJetPhi[iJetCategory][iCentrality];       // Jet phi histograms
      fhJetEta[iJetCategory][iCentrality] = in.fhJetEta[iJetCategory][iCentrality];       // Jet eta histograms
      fhJetEtaPhi[iJetCategory][iCentrality] = in.fhJetEtaPhi[iJetCategory][iCentrality]; // 2D eta-phi histogram for jets
    } // Single jet categories loop
    
    // Event correlation type loop
    for(int iCorrelationType = 0; iCorrelationType < knCorrelationTypes; iCorrelationType++){
      
      // Loop over track categories
      for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
        fhTrackPt[iTrackType][iCorrelationType][iCentrality] = in.fhTrackPt[iTrackType][iCorrelationType][iCentrality];   // Track pT histograms
        
        // Loop over track pT bins
        for(int iTrackPt = 0; iTrackPt < knTrackPtBins + 1; iTrackPt++){
          fhTrackPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt] = in.fhTrackPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt];    // Track phi histograms
          fhTrackEta[iTrackType][iCorrelationType][iCentrality][iTrackPt] = in.fhTrackEta[iTrackType][iCorrelationType][iCentrality][iTrackPt];    // Track eta histograms
          fhTrackEtaPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt] = in.fhTrackEtaPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt]; // 2D eta-phi histogram for track
        } // Track pT loop
        
      } // Track category loop
      
      // Loop over jet-track correlation types
      for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
        
        // Loop over track pT bins
        for(int iTrackPt = 0; iTrackPt < knTrackPtBins; iTrackPt++){
          fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt] = in.fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt];         // DeltaPhi between jet and track
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt] = in.fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt]; // DeltaEta and deltaPhi between jet and track
          
          // Loop over deltaPhi bins
          for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
            fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iCentrality][iTrackPt][iDeltaPhi] = in.fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iCentrality][iTrackPt][iDeltaPhi]; // DeltaEta between jet and track
          } // DeltaPhi loop
        } // Track pT loop
      } // Jet-track correlation type loop
    } // Event correlation type loop
    
    // Jet shape histograms
    for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
      for(int iTrackPt = 0; iTrackPt < knTrackPtBins; iTrackPt++){
        for(int iJetShape = 0; iJetShape < knJetShapeTypes; iJetShape++){
          fhJetShape[iJetShape][iJetTrack][iCentrality][iTrackPt] = in.fhJetShape[iJetShape][iJetTrack][iCentrality][iTrackPt];
        } // Jet shape type loop
      } // Track pT loop
    } // Jet-track correlation type loop
  } // Centrality loop
}

/*
 * Destructor
 */
DijetHistogramManager::~DijetHistogramManager(){
  delete fCard;
  delete fMethods;
}

/*
 * Apply the mixed event correction to all jet-track correlation histograms that are selected for analysis
 * After that subtract the background form the mixed event corrected distributions
 */
void DijetHistogramManager::ProcessHistograms(){
  DoMixedEventCorrection();  // Mixed event correction needs to be done first, as we need the corrected histograms for the background subtraction
  SubtractBackgroundAndCalculateJetShape(); // Subtract the background and take projections of processed two-dimensional histograms. After that, calculate jet shape
}

/*
 * Apply mixed event correction to all jet-track correlation histograms that are selected for analysis
 */
void DijetHistogramManager::DoMixedEventCorrection(){
  
  /*
   * Because background subtraction always needs information about leading and subleading jets, only loop over half the array
   * and do mixed event correction and at the same time for corresponding leading and subleading jet track correlation histograms.
   * Note that the array checks whether the correlation histogram is loaded instead of drawn, because the loop is only over
   * the leading jet indices and these are loaded but not drawn in case only subleading jet track correlation histograms are
   * selected to be drawn.
   */
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations/2; iJetTrack++){
    if(!fLoadJetTrackCorrelations[iJetTrack]) continue; // Only correct the histograms that are selected for analysis
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
        
        // Do the mixed event correction for leading jet-track correlation histogram
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kCorrected][iCentralityBin][iTrackPtBin] = fMethods->MixedEventCorrect(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kSameEvent][iCentralityBin][iTrackPtBin],fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kMixedEvent][iCentralityBin][iTrackPtBin],fhJetTrackDeltaEtaDeltaPhi[iJetTrack+knJetTrackCorrelations/2][kMixedEvent][iCentralityBin][iTrackPtBin]);
        
        // Do the mixed event correction for subleading jet-track correlation histogram
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack+knJetTrackCorrelations/2][kCorrected][iCentralityBin][iTrackPtBin] = fMethods->MixedEventCorrect(fhJetTrackDeltaEtaDeltaPhi[iJetTrack+knJetTrackCorrelations/2][kSameEvent][iCentralityBin][iTrackPtBin],fhJetTrackDeltaEtaDeltaPhi[iJetTrack+knJetTrackCorrelations/2][kMixedEvent][iCentralityBin][iTrackPtBin],fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kMixedEvent][iCentralityBin][iTrackPtBin]);
        
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation category loop
}

/*
 * Subtract the background and take projections of processed two-dimensional histograms
 */
void DijetHistogramManager::SubtractBackgroundAndCalculateJetShape(){
  
  // Helper variables to make the code more readable
  char histogramName[200];
  int nBins;
  int connectedIndex;
  
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    if(!fLoadJetTrackCorrelations[iJetTrack]) continue; // Only correct the histograms that are selected for analysis
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){

        // Get the subleading/leading jet index connected to the currect leading/subleading correlation type
        connectedIndex = GetConnectedIndex(iJetTrack);
        
        // Subtract the background from the mixed event corrected histogram
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackgroundSubtracted][iCentralityBin][iTrackPtBin] = fMethods->SubtractBackground(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kCorrected][iCentralityBin][iTrackPtBin],fhJetTrackDeltaEtaDeltaPhi[connectedIndex][kCorrected][iCentralityBin][iTrackPtBin]);
        
        // Get also the background and background overlap region for QA purposes
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackground][iCentralityBin][iTrackPtBin] = fMethods->GetBackground();
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackgroundOverlap][iCentralityBin][iTrackPtBin] = fMethods->GetBackgroundOverlap();
        
        // Calculate the jet shape from the background subtracted histogram
        fhJetShape[kJetShape][iJetTrack][iCentralityBin][iTrackPtBin] = fMethods->GetJetShape(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackgroundSubtracted][iCentralityBin][iTrackPtBin]);
        
        // Get the number of two-dimensional histogram bins used for each deltaR bin in the jet shape histogram
        fhJetShape[kJetShapeBinCount][iJetTrack][iCentralityBin][iTrackPtBin] = fMethods->GetJetShapeCounts();
        
        // Get the mapping histogram of Rbins to deltaPhi-deltaEta bins
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kJetShapeBinMap][iCentralityBin][iTrackPtBin] = fMethods->GetJetShapeBinMap();
        
        // Project the deltaPhi and deltaEta histograms from the processed two-dimensional histograms
        for(int iCorrelationType = kCorrected; iCorrelationType < knCorrelationTypes; iCorrelationType++){
          
          sprintf(histogramName,"%sDeltaPhiProjection",fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->GetName());
          nBins = fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->GetNbinsY();
          fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin] = (TH1D*) fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->ProjectionX(histogramName,1,nBins)->Clone();  // Exclude underflow and overflow bins by specifying range
          
          for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
            sprintf(histogramName,"%sDeltaEtaProjection%d",fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->GetName(),iDeltaPhi);
            fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaPhi] = (TH1D*) fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->ProjectionY(histogramName,fLowDeltaPhiBinIndices[iDeltaPhi],fHighDeltaPhiBinIndices[iDeltaPhi])->Clone();
          } // DeltaPhi loop
        } // Correlation type loop
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation category loop
}

/*
 * Get the index that of the same type of leading/subleading jet-track correlation as the
 * given subleading/leading jet-track correlation index.
 */
int DijetHistogramManager::GetConnectedIndex(const int jetTrackIndex) const{
  int connectedIndex = jetTrackIndex + knJetTrackCorrelations/2;
  if(connectedIndex >= knJetTrackCorrelations) connectedIndex -= knJetTrackCorrelations;
  return connectedIndex;
}

/*
 * Load all the selected histograms from the inputfile
 */
void DijetHistogramManager::LoadHistograms(){
  
  // Load the event information histograms
  if(fLoadEventInformation){
    fhVertexZ = (TH1D*) fInputFile->Get("vertexZ");                 // Vertex z position
    fhEvents = (TH1D*) fInputFile->Get("nEvents");                  // Number of events surviving different event cuts
    fhTrackCuts = (TH1D*) fInputFile->Get("trackCuts");             // Number of tracks surviving different track cuts
    fhCentrality = (TH1D*) fInputFile->Get("centrality");           // Centrality in all events
    fhCentralityDijet = (TH1D*) fInputFile->Get("centralityDijet"); // Centrality in dijet events
  }
  
  // Load single jet histograms
  LoadSingleJetHistograms();
  
  // Load dijet histograms
  LoadDijetHistograms();
  
  // Load track histograms
  LoadTrackHistograms();
  
  // Load all track jet correlation histograms
  LoadJetTrackCorrelationHistograms();
  
}

/*
 * Loader for single jet histograms
 *
 * THnSparse for single jets:
 *
 *   Histogram name: leadingJet/subleadingJet/anyJet
 *
 *     Axis index       Content of axis         Exception
 * ----------------------------------------------------------
 *       Axis 0             Jet pT
 *       Axis 1             Jet phi
 *       Axis 2             Jet eta
 *       Axis 3         Dijet asymmetry    (for anyJet: Centrality)
 *       Axis 4           Centrality       (for anyJet: Nothing)
 */
void DijetHistogramManager::LoadSingleJetHistograms(){
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  
  int centralityIndex[] = {4,4,3}; // TODO: Change the main analysis file such that anyJet has dummy axis for asymmetry to simplify loading here
  
  for(int iJetCategory = 0; iJetCategory < knSingleJetCategories; iJetCategory++){
    if(!fLoadSingleJets[iJetCategory]) continue;  // Only load the selected histograms
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      
      // Select the bin indices
      if(iCentralityBin == fLastLoadedCentralityBin) duplicateRemoverCentrality = 0;
      lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
      higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
      
      fhJetPt[iJetCategory][iCentralityBin] = FindHistogram(fInputFile,fSingleJetHistogramName[iJetCategory],0,centralityIndex[iJetCategory],lowerCentralityBin,higherCentralityBin);
      fhJetPhi[iJetCategory][iCentralityBin] = FindHistogram(fInputFile,fSingleJetHistogramName[iJetCategory],1,centralityIndex[iJetCategory],lowerCentralityBin,higherCentralityBin);
      fhJetEta[iJetCategory][iCentralityBin] = FindHistogram(fInputFile,fSingleJetHistogramName[iJetCategory],2,centralityIndex[iJetCategory],lowerCentralityBin,higherCentralityBin);
      fhJetEtaPhi[iJetCategory][iCentralityBin] = FindHistogram2D(fInputFile,fSingleJetHistogramName[iJetCategory],1,2,centralityIndex[iJetCategory],lowerCentralityBin,higherCentralityBin);
    } // Loop over centrality bins
  } // Loop over single jet categories
}

/*
 * Loader for dijet histograms
 *
 *  Arguments:
 *    TH1D *hDeltaPhi[knCentralityBins] = Array of dijet deltaPhi histogams
 *    TH1D *hAsymmetry[knCentralityBins] = Array of dijet asymmetry histograms
 *    TH2D *hLeadingSubleadingPt[knCentralityBins] = Array of leading jet pT vs. subleading jet pT histograms
 *    const char* name = Name of the histogram in the input file
 *    const int iCentralityAxis = Index of centrality axis in THnSparse
 *
 * THnSparse for dijets:
 *
 *   Histogram name        Axis index       Content of axis
 * ----------------------------------------------------------
 *        dijet              Axis 0         Leading jet pT
 *        dijet              Axis 1        Subleading jet pT
 *        dijet              Axis 2         Dijet deltaPhi
 *        dijet              Axis 3         Dijet asymmetry
 *        dijet              Axis 4           Centrality
 */
void DijetHistogramManager::LoadDijetHistograms(){
  
  if(!fLoadDijetHistograms) return; // Do not load the histograms if they are not selected for drawing
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  
  for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
    
    // Select the bin indices
    if(iCentralityBin == fLastLoadedCentralityBin) duplicateRemoverCentrality = 0;
    lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
    higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
    
    fhDijetDphi[iCentralityBin] = FindHistogram(fInputFile,"dijet",2,4,lowerCentralityBin,higherCentralityBin);
    fhDijetAsymmetry[iCentralityBin] = FindHistogram(fInputFile,"dijet",3,4,lowerCentralityBin,higherCentralityBin);
    fhDijetLeadingVsSubleadingPt[iCentralityBin] = FindHistogram2D(fInputFile,"dijet",0,1,4,lowerCentralityBin,higherCentralityBin);
  }
}

/*
 * Loader for track histograms
 *
 * THnSparse for tracks:
 *
 *   Histogram name: track/trackUncorrected
 *
 *     Axis index       Content of axis
 * ----------------------------------------
 *       Axis 0            Track pT
 *       Axis 1            Track phi
 *       Axis 2            Track eta
 *       Axis 3            Centrality
 *       Axis 4         Correlation type
 */
void DijetHistogramManager::LoadTrackHistograms(){
  
  // Define arrays to help find the histograms
  int axisIndices[3] = {0};
  int lowLimits[3] = {0};
  int highLimits[3] = {0};
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int duplicateRemoverTrackPt = -1;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  
  // Loop over all track histograms
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    if(!fLoadTracks[iTrackType]) continue;  // Only load the selected track types
    for(int iCorrelationType = 0; iCorrelationType <= kMixedEvent; iCorrelationType++){  // Data file contains only same and mixed event distributions
      for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
        
        // Select the bin indices
        if(iCentralityBin == fLastLoadedCentralityBin) duplicateRemoverCentrality = 0;
        lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
        higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
        
        // Setup axes with restrictions, (3 = centrality, 4 = correlation type)
        axisIndices[0] = 3; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;
        axisIndices[1] = 4; lowLimits[1] = iCorrelationType+1; highLimits[1] = iCorrelationType+1;
        
        fhTrackPt[iTrackType][iCorrelationType][iCentralityBin] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],0,2,axisIndices,lowLimits,highLimits);
        fhTrackPhi[iTrackType][iCorrelationType][iCentralityBin][knTrackPtBins] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],1,2,axisIndices,lowLimits,highLimits);
        fhTrackEta[iTrackType][iCorrelationType][iCentralityBin][knTrackPtBins] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],2,2,axisIndices,lowLimits,highLimits);
        fhTrackEtaPhi[iTrackType][iCorrelationType][iCentralityBin][knTrackPtBins] = FindHistogram2D(fInputFile,fTrackHistogramNames[iTrackType],1,2,2,axisIndices,lowLimits,highLimits);
        
        for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
          
          // Select the bin indices for track pT
          lowerTrackPtBin = fFineTrackPtBinIndices[iTrackPtBin];
          higherTrackPtBin = fFineTrackPtBinIndices[iTrackPtBin+1]+duplicateRemoverTrackPt;
          
          // Add restriction for pT axis (0)
          axisIndices[2] = 0; lowLimits[2] = lowerTrackPtBin; highLimits[2] = higherTrackPtBin;
          
          // Read the angle histograms in track pT bins
          fhTrackPhi[iTrackType][iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],1,3,axisIndices,lowLimits,highLimits);
          fhTrackEta[iTrackType][iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],2,3,axisIndices,lowLimits,highLimits);
          fhTrackEtaPhi[iTrackType][iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram2D(fInputFile,fTrackHistogramNames[iTrackType],1,2,3,axisIndices,lowLimits,highLimits);
          
        } // Track pT loop
      } // Centrality loop
    } // Correlation type loop
  } // Track category loop
}

/*
 * Loader for track jet correlation histograms
 *
 * THnSparse for tracks:
 *
 *   Histogram name: trackLeadingJet/trackLeadingJetUncorrected/trackLeadingJetPtWeighted
                     trackSubleadingJet/trackSubleadingJetUncorrected/trackSubleadingJetPtWeighted
 *
 *     Axis index                  Content of axis
 * -----------------------------------------------------------
 *       Axis 0                        Track pT
 *       Axis 1             DeltaPhi between track and jet
 *       Axis 2             DeltaEta between track and jet
 *       Axis 3                    Dijet asymmetry
 *       Axis 4                       Centrality
 *       Axis 5                    Correlation type
 */
void DijetHistogramManager::LoadJetTrackCorrelationHistograms(){
  
  // Define arrays to help find the histograms
  int axisIndices[4] = {0};
  int lowLimits[4] = {0};
  int highLimits[4] = {0};
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int duplicateRemoverTrackPt = -1;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  
  // Load all the histograms from the files
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    if(!fLoadJetTrackCorrelations[iJetTrack]) continue; // Only load categories of correlation that are selected
    for(int iCorrelationType = 0; iCorrelationType <= kMixedEvent; iCorrelationType++){ // Data file contains only same and mixed event distributions
      for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
        
        // Select the bin indices
        if(iCentralityBin == fLastLoadedCentralityBin) duplicateRemoverCentrality = 0;
        lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
        higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
        
        for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
          
          // Select the bin indices for track pT
          lowerTrackPtBin = fTrackPtBinIndices[iTrackPtBin];
          higherTrackPtBin = fTrackPtBinIndices[iTrackPtBin+1]+duplicateRemoverTrackPt;
          
          // Setup the axes with restrictions, that are common for all jet-track correlation histograms
          axisIndices[0] = 5; lowLimits[0] = iCorrelationType+1; highLimits[0] = iCorrelationType+1;   // Same/mixed event
          axisIndices[1] = 4; lowLimits[1] = lowerCentralityBin; highLimits[1] = higherCentralityBin;  // Centrality
          axisIndices[2] = 0; lowLimits[2] = lowerTrackPtBin;    highLimits[2] = higherTrackPtBin;     // Track pT
          
          fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram(fInputFile,fJetTrackHistogramNames[iJetTrack],1,3,axisIndices,lowLimits,highLimits);
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram2D(fInputFile,fJetTrackHistogramNames[iJetTrack],1,2,3,axisIndices,lowLimits,highLimits);
          
          // DeltaPhi binning for deltaEta histogram
          for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
            axisIndices[3] = 1; lowLimits[3] = fLowDeltaPhiBinIndices[iDeltaPhi]; highLimits[3] = fHighDeltaPhiBinIndices[iDeltaPhi];  // DeltaPhi
            fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaPhi] = FindHistogram(fInputFile,fJetTrackHistogramNames[iJetTrack],2,4,axisIndices,lowLimits,highLimits);
          } // DeltaPhi loop
        } // Track pT loop
      } // Centrality loop
    } // Correlation type loop
  } // Jet-track correlation category loop
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
 *   int *axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
 */
TH2D* DijetHistogramManager::FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex){
  
  // Read the histogram with the given name from the file
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for eeach histogram that is read from the file
  char newName[200];
  sprintf(newName,"%s",histogramArray->GetName());
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    sprintf(newName,"%s_%d=%d-%d",newName,axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]);
  }
  
  // Project out the histogram and give it the created unique name
  TH2D *projectedHistogram = (TH2D*) histogramArray->Projection(yAxis,xAxis);
  projectedHistogram->SetName(newName);
  
  // Apply bin width normalization to the projected histogram
  projectedHistogram->Scale(1.0,"width");
  
  // Return the projected histogram
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
TH2D* DijetHistogramManager::FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2, int lowBinIndex2, int highBinIndex2){
  int restrictionAxes[2] = {restrictionAxis,restrictionAxis2};
  int lowBinIndices[2] = {lowBinIndex,lowBinIndex2};
  int highBinIndices[2] = {highBinIndex,highBinIndex2};
  int nAxes = 2;
  if(highBinIndex2 == 0 && lowBinIndex2 == 0) nAxes = 1;
  return FindHistogram2D(inputFile,name,xAxis,yAxis,nAxes,restrictionAxes,lowBinIndices,highBinIndices);
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH1D
 *   int nAxes = Number of axes that are restained for the projection
 *   int *axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
 */
TH1D* DijetHistogramManager::FindHistogram(TFile *inputFile, const char *name, int xAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex){
  
  // Read the histogram with the given name from the file
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for eeach histogram that is read from the file
  char newName[200];
  sprintf(newName,"%s",histogramArray->GetName());
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    sprintf(newName,"%s_%d=%d-%d",newName,axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]);
  }
  
  // Project out the histogram and give it the created unique name
  TH1D *projectedHistogram = (TH1D*) histogramArray->Projection(xAxis);
  projectedHistogram->SetName(newName);
  
  // Apply bin width normalization to the projected histogram
  projectedHistogram->Scale(1.0,"width");
  
  // Return the projected histogram
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
TH1D* DijetHistogramManager::FindHistogram(TFile *inputFile, const char *name, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2, int lowBinIndex2, int highBinIndex2){
  int restrictionAxes[2] = {restrictionAxis,restrictionAxis2};
  int lowBinIndices[2] = {lowBinIndex,lowBinIndex2};
  int highBinIndices[2] = {highBinIndex,highBinIndex2};
  int nAxes = 2;
  if(highBinIndex2 == 0 && lowBinIndex2 == 0) nAxes = 1;
  return FindHistogram(inputFile,name,xAxis,nAxes,restrictionAxes,lowBinIndices,highBinIndices);
}

/*
 * Read the bin indices for given bin borders
 *
 *  Arguments:
 *   const int nBins = Number of bins for the indices
 *   int *binIndices = Array of integers to be filled with bin index information read from the file
 *   const double *binBorders = Array for bin borders that are searched from the file
 *   const int iAxis = Index of the axis used for reading bin indices
 */
void DijetHistogramManager::SetBinIndices(const int nBins, double *copyBinBorders, int *binIndices, const double *binBorders, const int iAxis){
  TH1D* hBinner = FindHistogram(fInputFile,"trackLeadingJet",iAxis,0,0,0);
  for(int iBin = 0; iBin < nBins+1; iBin++){
    copyBinBorders[iBin] = binBorders[iBin];
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
void DijetHistogramManager::SetBinIndices(const int nBins, int *lowBinIndices, int *highBinIndices, const double *lowBinBorders, const double *highBinBorders, const int iAxis){
  TH1D* hBinner = FindHistogram(fInputFile,"trackLeadingJet",iAxis,0,0,0);
  for(int iBin = 0; iBin < nBins; iBin++){
    lowBinIndices[iBin] = hBinner->GetXaxis()->FindBin(lowBinBorders[iBin]);
    highBinIndices[iBin] = hBinner->GetXaxis()->FindBin(highBinBorders[iBin]);
  }
}

/*
 * Set up centrality bin indices according to provided bin borders
 */
void DijetHistogramManager::SetCentralityBins(double *binBorders){
  SetBinIndices(knCentralityBins,fCentralityBinBorders,fCentralityBinIndices,binBorders,4);
}

/*
 * Set up track pT bin indices according to provided bin borders
 */
void DijetHistogramManager::SetTrackPtBins(double *binBorders){
  SetBinIndices(knTrackPtBins,fTrackPtBinBorders,fTrackPtBinIndices,binBorders,0);
  
  // The track histograms have finer pT binning, so we need to use different bin indices for them
  TH1D* hTrackPtBinner = FindHistogram(fInputFile,"track",0,0,0,0);
  for(int iTrackPt = 0; iTrackPt < knTrackPtBins+1; iTrackPt++){
    fFineTrackPtBinIndices[iTrackPt] = hTrackPtBinner->GetXaxis()->FindBin(binBorders[iTrackPt]);
  }
}

/*
 * Set up deltaPhi bin indices according to provided bin borders
 */
void DijetHistogramManager::SetDeltaPhiBins(double *lowBinBorders, double *highBinBorders, TString deltaPhiStrings[knDeltaPhiBins], TString compactDeltaPhiStrings[knDeltaPhiBins]){
  SetBinIndices(knDeltaPhiBins,fLowDeltaPhiBinIndices,fHighDeltaPhiBinIndices,lowBinBorders,highBinBorders,1);
  for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
    fDeltaPhiString[iDeltaPhi] = deltaPhiStrings[iDeltaPhi];
    fCompactDeltaPhiString[iDeltaPhi] = compactDeltaPhiStrings[iDeltaPhi];
  }
}

// Setter for drawing event information
void DijetHistogramManager::SetLoadEventInformation(const bool loadOrNot){
  fLoadEventInformation = loadOrNot;
}

// Setter for drawing dijet histograms
void DijetHistogramManager::SetLoadDijetHistograms(const bool loadOrNot){
  fLoadDijetHistograms = loadOrNot;
}

// Setter for drawing leading jet histograms
void DijetHistogramManager::SetLoadLeadingJetHistograms(const bool loadOrNot){
  fLoadSingleJets[kLeadingJet] = loadOrNot;
}

// Setter for drawing subleading jet histograms
void DijetHistogramManager::SetLoadSubleadingJetHistograms(const bool loadOrNot){
  fLoadSingleJets[kSubleadingJet] = loadOrNot;
}

// Setter for drawing all jet histograms
void DijetHistogramManager::SetLoadAnyJetHistograms(const bool loadOrNot){
  fLoadSingleJets[kAnyJet] = loadOrNot;
}

// Setter for drawing jet histograms
void DijetHistogramManager::SetLoadAllJets(const bool drawLeading, const bool drawSubleading, const bool drawAny){
  SetLoadLeadingJetHistograms(drawLeading);
  SetLoadSubleadingJetHistograms(drawSubleading);
  SetLoadAnyJetHistograms(drawAny);
}

// Setter for drawing tracks
void DijetHistogramManager::SetLoadTracks(const bool loadOrNot){
  fLoadTracks[kTrack] = loadOrNot;
}

// Setter for drawing uncorrected tracks
void DijetHistogramManager::SetLoadTracksUncorrected(const bool loadOrNot){
  fLoadTracks[kUncorrectedTrack] = loadOrNot;
}

// Setter for drawing track histograms
void DijetHistogramManager::SetLoadAllTracks(const bool drawTracks, const bool drawUncorrected){
  SetLoadTracks(drawTracks);
  SetLoadTracksUncorrected(drawUncorrected);
}

/*
 * Setter for drawing jet-track correlations.
 *
 * The method sets the drawing of the histogram defined by the primaryIndex.
 * The background subtraction histograms are conntructed using buth leading and subleading jet histograms.
 * Thus when drawing the leading/subleading jet histograms, we need to also load the subleading/leading jet
 * histograms to be able to perform the background subtraction.
 *
 *  Arguments:
 *   const bool loadOrNot = Flag whether these type or correlation should be drawn or not
 *   const int primaryIndex = Index of the primary jet-track correlation type
 *   const int connectedIndex = Index of the type connected to the primary type in background subtraction
 */
void DijetHistogramManager::SetLoadJetTrackCorrelations(const bool loadOrNot, const int primaryIndex, const int connectedIndex){
  
  // If we are loading the connected index, do not disable this since it is needed in background subtraction
  if(fLoadJetTrackCorrelations[connectedIndex]){
    fLoadJetTrackCorrelations[primaryIndex] = true;
  } else {
    fLoadJetTrackCorrelations[primaryIndex] = loadOrNot;
    fLoadJetTrackCorrelations[connectedIndex] = loadOrNot;
  }
  
}

// Setter for drawing leading jet-track correlations
void DijetHistogramManager::SetLoadTrackLeadingJetCorrelations(const bool loadOrNot){
  SetLoadJetTrackCorrelations(loadOrNot,kTrackLeadingJet,kTrackSubleadingJet);
}

// Setter for drawing uncorrected leading jet-track correlations
void DijetHistogramManager::SetLoadTrackLeadingJetCorrelationsUncorrected(const bool loadOrNot){
  SetLoadJetTrackCorrelations(loadOrNot,kUncorrectedTrackLeadingJet,kUncorrectedTrackSubleadingJet);
}

// Setter for drawing pT weighted leading jet-track correlations
void DijetHistogramManager::SetLoadTrackLeadingJetCorrelationsPtWeighted(const bool loadOrNot){
  SetLoadJetTrackCorrelations(loadOrNot,kPtWeightedTrackLeadingJet,kPtWeightedTrackSubleadingJet);
}

// Setter for drawing all correlations related to tracks and leading jets
void DijetHistogramManager::SetLoadAllTrackLeadingJetCorrelations(const bool drawLeading, const bool drawUncorrected, const bool drawPtWeighted){
  SetLoadTrackLeadingJetCorrelations(drawLeading);
  SetLoadTrackLeadingJetCorrelationsUncorrected(drawUncorrected);
  SetLoadTrackLeadingJetCorrelationsPtWeighted(drawPtWeighted);
}

// Setter for drawing subleading jet-track correlations
void DijetHistogramManager::SetLoadTrackSubleadingJetCorrelations(const bool loadOrNot){
  SetLoadJetTrackCorrelations(loadOrNot,kTrackSubleadingJet,kTrackLeadingJet);
}

// Setter for drawing uncorrected subleading jet-track correlations
void DijetHistogramManager::SetLoadTrackSubleadingJetCorrelationsUncorrected(const bool loadOrNot){
  SetLoadJetTrackCorrelations(loadOrNot,kUncorrectedTrackSubleadingJet,kUncorrectedTrackLeadingJet);
}

// Setter for drawing pT weighted subleading jet-track correlations
void DijetHistogramManager::SetLoadTrackSubleadingJetCorrelationsPtWeighted(const bool loadOrNot){
  SetLoadJetTrackCorrelations(loadOrNot,kPtWeightedTrackSubleadingJet,kPtWeightedTrackLeadingJet);
}

// Setter for drawing all correlations related to tracks and subleading jets
void DijetHistogramManager::SetLoadAllTrackSubleadingJetCorrelations(const bool drawSubleading, const bool drawUncorrected, const bool drawPtWeighted){
  SetLoadTrackSubleadingJetCorrelations(drawSubleading);
  SetLoadTrackSubleadingJetCorrelationsUncorrected(drawUncorrected);
  SetLoadTrackSubleadingJetCorrelationsPtWeighted(drawPtWeighted);
}

// Setter for drawn centrality bins
void DijetHistogramManager::SetCentralityBinRange(const int first, const int last){
  fFirstLoadedCentralityBin = first;
  fLastLoadedCentralityBin = last;
  
  // Sanity check for drawn centrality bins
  BinSanityCheck(knCentralityBins,fFirstLoadedCentralityBin,fLastLoadedCentralityBin);
}

// Setter for drawn track pT bins
void DijetHistogramManager::SetTrackPtBinRange(const int first, const int last){
  fFirstLoadedTrackPtBin = first;
  fLastLoadedTrackPtBin = last;
  
  // Sanity check for drawn track pT bins
  BinSanityCheck(knTrackPtBins,fFirstLoadedTrackPtBin,fLastLoadedTrackPtBin);
}

// Sanity check for set bins
void DijetHistogramManager::BinSanityCheck(const int nBins, int first, int last) const{
  if(first < 0) first = 0;
  if(last < first) last = first;
  if(last > nBins-1) last = nBins-1;
}

// Sanity check for input bin index
int DijetHistogramManager::BinIndexCheck(const int nBins, const int binIndex) const{
  if(binIndex < 0) return 0;
  if(binIndex > nBins-1) return nBins-1;
  return binIndex;
}

// Setter for used DijetMethods
void DijetHistogramManager::SetDijetMethods(DijetMethods* newMethods){
  fMethods = newMethods;
}

// Getter for the number of centrality bins
int DijetHistogramManager::GetNCentralityBins() const{
  return knCentralityBins;
}

// Getter for the number of track pT bins
int DijetHistogramManager::GetNTrackPtBins() const{
  return knTrackPtBins;
}

// Getter for correlation type string
TString DijetHistogramManager::GetCorrelationTypeString(int iCorrelationType) const{
  iCorrelationType = BinIndexCheck(knCorrelationTypes,iCorrelationType);
  return fCorrelationTypeString[iCorrelationType];
}

// Getter for compact correlation type string
TString DijetHistogramManager::GetCompactCorrelationTypeString(int iCorrelationType) const{
  iCorrelationType = BinIndexCheck(knCorrelationTypes,iCorrelationType);
  return fCompactCorrelationTypeString[iCorrelationType];
}

// Getter for deltaPhi string
TString DijetHistogramManager::GetDeltaPhiString(int iDeltaPhiRegion) const{
  iDeltaPhiRegion = BinIndexCheck(knDeltaPhiBins,iDeltaPhiRegion);
  return fDeltaPhiString[iDeltaPhiRegion];
}

// Getter for compact deltaPhi string
TString DijetHistogramManager::GetCompactDeltaPhiString(int iDeltaPhiRegion) const{
  iDeltaPhiRegion = BinIndexCheck(knDeltaPhiBins,iDeltaPhiRegion);
  return fCompactDeltaPhiString[iDeltaPhiRegion];
}

// Getter for jet-track correlation histogram name
const char* DijetHistogramManager::GetJetTrackHistogramName(int iJetTrackCorrelation) const{
  iJetTrackCorrelation = BinIndexCheck(knJetTrackCorrelations,iJetTrackCorrelation);
  return fJetTrackHistogramNames[iJetTrackCorrelation];
}

// Getter for name suitable for x-axis in a given jet-track correlation histogram
const char* DijetHistogramManager::GetJetTrackAxisName(int iJetTrackCorrelation) const{
  iJetTrackCorrelation = BinIndexCheck(knJetTrackCorrelations,iJetTrackCorrelation);
  return fJetTrackAxisNames[iJetTrackCorrelation];
}

// Getter for track histogram name
const char* DijetHistogramManager::GetTrackHistogramName(int iTrackType) const{
  iTrackType = BinIndexCheck(knTrackCategories,iTrackType);
  return fTrackHistogramNames[iTrackType];
}

// Getter for name suitable for x-axis in a given track histogram
const char* DijetHistogramManager::GetTrackAxisName(int iTrackType) const{
  iTrackType = BinIndexCheck(knTrackCategories,iTrackType);
  return fTrackAxisNames[iTrackType];
}

// Getter for single jet histogram name
const char* DijetHistogramManager::GetSingleJetHistogramName(int iJetType) const{
  iJetType = BinIndexCheck(knSingleJetCategories,iJetType);
  return fSingleJetHistogramName[iJetType];
}

// Getter for name suitable for x-axis in a given single jet histogram
const char* DijetHistogramManager::GetSingleJetAxisName(int iJetType) const{
  iJetType = BinIndexCheck(knSingleJetCategories,iJetType);
  return fSingleJetAxisNames[iJetType];
}

// Getter for jet shape histogram name
const char* DijetHistogramManager::GetJetShapeHistogramName(int iJetShapeType) const{
  iJetShapeType = BinIndexCheck(knJetShapeTypes,iJetShapeType);
  return fJetShapeHistogramName[iJetShapeType];
}

// Getter for name suitable for x-axis in a given jet shape histogram
const char* DijetHistogramManager::GetJetShapeAxisName(int iJetShapeType) const{
  iJetShapeType = BinIndexCheck(knJetShapeTypes,iJetShapeType);
  return fJetShapeYAxisNames[iJetShapeType];
}

// Getter for collision system
TString DijetHistogramManager::GetSystem() const{
  return fCard->GetDataType();
}

// Getter for i:th centrality bin border
double DijetHistogramManager::GetCentralityBinBorder(const int iCentrality) const{
  return fCentralityBinBorders[iCentrality];
}

// Getter for i:th track pT bin border
double DijetHistogramManager::GetTrackPtBinBorder(const int iTrackPt) const{
  return fTrackPtBinBorders[iTrackPt];
}

// Getters for event information histograms

// Getter for z-vertex histogram
TH1D* DijetHistogramManager::GetHistogramVertexZ() const{
  return fhVertexZ;
}

// Getter for histogram for number of events surviving different event cuts
TH1D* DijetHistogramManager::GetHistogramEvents() const{
  return fhEvents;
}

// Getter for histogram for number of tracks surviving different track cuts
TH1D* DijetHistogramManager::GetHistogramTrackCuts() const{
  return fhTrackCuts;
}

// Getter for centrality histogram in all events
TH1D* DijetHistogramManager::GetHistogramCentrality() const{
  return fhCentrality;
}

// Getter for centrality histogram in dijet events
TH1D* DijetHistogramManager::GetHistogramCentralityDijet() const{
  return fhCentralityDijet;
}

// Getters for single jet histograms

// Getter for jet pT histograms
TH1D* DijetHistogramManager::GetHistogramJetPt(const int iJetType, const int iCentrality) const{
  return fhJetPt[iJetType][iCentrality];
}

// Getter for jet phi histograms
TH1D* DijetHistogramManager::GetHistogramJetPhi(const int iJetType, const int iCentrality) const{
  return fhJetPhi[iJetType][iCentrality];
}

// Getter for jet eta histograms
TH1D* DijetHistogramManager::GetHistogramJetEta(const int iJetType, const int iCentrality) const{
  return fhJetEta[iJetType][iCentrality];
}

// Getter for 2D eta-phi histogram for jets
TH2D* DijetHistogramManager::GetHistogramJetEtaPhi(const int iJetType, const int iCentrality) const{
  return fhJetEtaPhi[iJetType][iCentrality];
}

// Getters for dijet histograms

// Getter for dijet deltaPhi histograms
TH1D* DijetHistogramManager::GetHistogramDijetDeltaPhi(const int iCentrality) const{
  return fhDijetDphi[iCentrality];
}

// Getter for dijet asymmetry histograms
TH1D* DijetHistogramManager::GetHistogramDijetAsymmetry(const int iCentrality) const{
  return fhDijetAsymmetry[iCentrality];
}

// Getter for leading versus subleading jet pT 2D histograms
TH2D* DijetHistogramManager::GetHistogramDijetLeadingVsSubleadingPt(const int iCentrality) const{
  return fhDijetLeadingVsSubleadingPt[iCentrality];
}

// Getters for histograms for tracks in dijet events

// Getter for track pT histograms
TH1D* DijetHistogramManager::GetHistogramTrackPt(const int iTrackType, const int iCorrelationType, const int iCentrality) const{
  return fhTrackPt[iTrackType][iCorrelationType][iCentrality];
}

// Getter for track phi histograms
TH1D* DijetHistogramManager::GetHistogramTrackPhi(const int iTrackType, const int iCorrelationType, const int iCentrality, const int iTrackPt) const{
  return fhTrackPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt];
}

// Getter for track eta histograms
TH1D* DijetHistogramManager::GetHistogramTrackEta(const int iTrackType, const int iCorrelationType, const int iCentrality, const int iTrackPt) const{
  return fhTrackEta[iTrackType][iCorrelationType][iCentrality][iTrackPt];
}

// Getter for 2D eta-phi histogram for track
TH2D* DijetHistogramManager::GetHistogramTrackEtaPhi(const int iTrackType, const int iCorrelationType, const int iCentrality, const int iTrackPt) const{
  return fhTrackEtaPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt];
}

// Getters for track-leading jet correlation histograms

// Getter for deltaPhi between jet and track
TH1D* DijetHistogramManager::GetHistogramJetTrackDeltaPhi(const int iJetTrackCorrelation, const int iCorrelationType, const int iCentrality, const int iTrackPt) const{
  return fhJetTrackDeltaPhi[iJetTrackCorrelation][iCorrelationType][iCentrality][iTrackPt];
}

// Getter for deltaEta between jet and track
TH1D* DijetHistogramManager::GetHistogramJetTrackDeltaEta(const int iJetTrackCorrelation, const int iCorrelationType, const int iCentrality, const int iTrackPt, const int iDeltaPhiRegion) const{
  return fhJetTrackDeltaEta[iJetTrackCorrelation][iCorrelationType][iCentrality][iTrackPt][iDeltaPhiRegion];
}

// Getter for deltaEta and deltaPhi between jet and track
TH2D* DijetHistogramManager::GetHistogramJetTrackDeltaEtaDeltaPhi(const int iJetTrackCorrelation, const int iCorrelationType, const int iCentrality, const int iTrackPt) const{
  return fhJetTrackDeltaEtaDeltaPhi[iJetTrackCorrelation][iCorrelationType][iCentrality][iTrackPt];
}

// Getters for jet shape histograms
TH1D* DijetHistogramManager::GetHistogramJetShape(const int iJetShapeType, const int iJetTrackCorrelation, const int iCentrality, const int iTrackPt) const{
  return fhJetShape[iJetShapeType][iJetTrackCorrelation][iCentrality][iTrackPt];
}

/*
 * Getter for any one-dimensional histogram based on input string
 *
 *  Arguments:
 *   TString name = Name corresponding to loaded histogram
 *   int bin1 = First bin index for the histogram
 *   int bin2 = Second bin index for the histogram
 *   int bin3 = Third bin index for the histogram
 *   int bin4 = Fourth bin index for the histogram
 *   int bin5 = Fifth bin index for the histogram
 *
 *   return: Histogram corresponding to given name and bins
 */
TH1D* DijetHistogramManager::GetOneDimensionalHistogram(TString name, int bin1, int bin2, int bin3, int bin4, int bin5) const{
  if(name.EqualTo("vertexz",TString::kIgnoreCase) || name.EqualTo("fhvertexz",TString::kIgnoreCase)) return GetHistogramVertexZ();
  if(name.EqualTo("events",TString::kIgnoreCase) || name.EqualTo("fhevents",TString::kIgnoreCase)) return GetHistogramEvents();
  if(name.EqualTo("trackcuts",TString::kIgnoreCase) || name.EqualTo("fhtrackcuts",TString::kIgnoreCase)) return GetHistogramTrackCuts();
  if(name.EqualTo("centrality",TString::kIgnoreCase) || name.EqualTo("fhcentrality",TString::kIgnoreCase)) return GetHistogramCentrality();
  if(name.EqualTo("centralitydijet",TString::kIgnoreCase) || name.EqualTo("fhcentralitydijet",TString::kIgnoreCase)) return GetHistogramCentralityDijet();
  if(name.EqualTo("jetpt",TString::kIgnoreCase) || name.EqualTo("fhjetpt",TString::kIgnoreCase)) return GetHistogramJetPt(bin1,bin2);
  if(name.EqualTo("jetphi",TString::kIgnoreCase) || name.EqualTo("fhjetphi",TString::kIgnoreCase)) return GetHistogramJetPhi(bin1,bin2);
  if(name.EqualTo("jeteta",TString::kIgnoreCase) || name.EqualTo("fhjeteta",TString::kIgnoreCase)) return GetHistogramJetEta(bin1,bin2);
  if(name.EqualTo("dijetdeltaphi",TString::kIgnoreCase) || name.EqualTo("dijetdphi",TString::kIgnoreCase) || name.EqualTo("fhdijetdphi",TString::kIgnoreCase)) return GetHistogramDijetDeltaPhi(bin1);
  if(name.EqualTo("dijetasymmetry",TString::kIgnoreCase) || name.EqualTo("fhdijetasymmetry",TString::kIgnoreCase)) return GetHistogramDijetAsymmetry(bin1);
  if(name.EqualTo("trackpt",TString::kIgnoreCase) || name.EqualTo("fhtrackpt",TString::kIgnoreCase)) return GetHistogramTrackPt(bin1,bin2,bin3);
  if(name.EqualTo("trackphi",TString::kIgnoreCase) || name.EqualTo("fhtrackphi",TString::kIgnoreCase)) return GetHistogramTrackPhi(bin1,bin2,bin3,bin4);
  if(name.EqualTo("tracketa",TString::kIgnoreCase) || name.EqualTo("fhtracketa",TString::kIgnoreCase)) return GetHistogramTrackEta(bin1,bin2,bin3,bin4);
  if(name.EqualTo("jettrackdeltaphi",TString::kIgnoreCase) || name.EqualTo("fhjettrackdeltaphi",TString::kIgnoreCase)) return GetHistogramJetTrackDeltaPhi(bin1,bin2,bin3,bin4);
  if(name.EqualTo("jettrackdeltaeta",TString::kIgnoreCase) || name.EqualTo("fhjettrackdeltaeta",TString::kIgnoreCase)) return GetHistogramJetTrackDeltaEta(bin1,bin2,bin3,bin4,bin5);
  if(name.EqualTo("jetshape",TString::kIgnoreCase) || name.EqualTo("fhjetshape",TString::kIgnoreCase)) return GetHistogramJetShape(bin1,bin2,bin3,bin4);
  return NULL;
}

/*
 * Getter for any two-dimensional histogram based on input string
 *
 *  Arguments:
 *   TString name = Name corresponding to loaded histogram
 *   int bin1 = First bin index for the histogram
 *   int bin2 = Second bin index for the histogram
 *   int bin3 = Third bin index for the histogram
 *   int bin4 = Fourth bin index for the histogram
 *
 *   return: Histogram corresponding to given name and bins
 */
TH2D* DijetHistogramManager::GetTwoDimensionalHistogram(TString name, int bin1, int bin2, int bin3, int bin4) const{
  if(name.EqualTo("jetetaphi",TString::kIgnoreCase) || name.EqualTo("fhjetetaphi",TString::kIgnoreCase)) return GetHistogramJetEtaPhi(bin1,bin2);
  if(name.EqualTo("dijetleadingvssubleadingpt",TString::kIgnoreCase) || name.EqualTo("fhdijetleadingvssubleadingpt",TString::kIgnoreCase)) return GetHistogramDijetLeadingVsSubleadingPt(bin1);
  if(name.EqualTo("tracketaphi",TString::kIgnoreCase) || name.EqualTo("fhtracketaphi",TString::kIgnoreCase)) return GetHistogramTrackEtaPhi(bin1,bin2,bin3,bin4);
  if(name.EqualTo("jettrackdeltaetadeltaphi",TString::kIgnoreCase) || name.EqualTo("fhjettrackdeltaetadeltaphi",TString::kIgnoreCase)) return GetHistogramJetTrackDeltaEtaDeltaPhi(bin1,bin2,bin3,bin4);
  return NULL;
}
