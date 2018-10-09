/*
 * Implementation of DijetHistogramManager
 */

// Root includes
#include <THnSparse.h>
#include <TPad.h>

// Own includes
#include "DijetHistogramManager.h"
#include "JffCorrector.h"

/*
 * Default constructor
 */
DijetHistogramManager::DijetHistogramManager() :
  fInputFile(NULL),
  fCard(NULL),
  fSystemAndEnergy(""),
  fCompactSystemAndEnergy(""),
  fApplyJffCorrection(false),
  fApplySpilloverCorrection(false),
  fApplySeagullCorrection(false),
  fLoadEventInformation(false),
  fLoadDijetHistograms(false),
  fLoad2DHistograms(false),
  fFirstLoadedCentralityBin(0),
  fLastLoadedCentralityBin(knCentralityBins-1),
  fFirstLoadedTrackPtBin(0),
  fLastLoadedTrackPtBin(knTrackPtBins-1)
{
  
  // Create a new DijetMethods and JffCorrector
  fMethods = new DijetMethods();
  fJffCorrectionFinder = new JffCorrector();
  
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
    fLowDeltaPhiBinBorders[iDeltaPhi] = 0;
    fHighDeltaPhiBinBorders[iDeltaPhi] = 0;
    fDeltaPhiString[iDeltaPhi] = "";
    fCompactDeltaPhiString[iDeltaPhi] = "";
  }
  
  // Initialize all the other histograms to null
  fhVertexZ = NULL;            // Vertex z position
  fhVertexZWeighted = NULL;    // Weighted vertex z-position (only meaningfull for MC)
  fhEvents = NULL;             // Number of events surviving different event cuts
  fhTrackCuts = NULL;          // Number of tracks surviving different track cuts
  fhTrackCutsInclusive = NULL; // Number of inclusive tracks surviving different track cuts
  fhCentrality = NULL;         // Centrality of all events
  fhCentralityWeighted = NULL; // Weighted centrality distribution in all events (only meaningful for MC)
  fhCentralityDijet = NULL;    // Centrality of dijet events
  fhPtHat = NULL;              // pT hat for MC events (only meaningful for MC)
  fhPtHatWeighted = NULL;      // Weighted pT hat distribution (only meaningful for MC)
  
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
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt] = NULL;      // DeltaEta and deltaPhi between jet and track
          
          // Loop over deltaEta bins
          for(int iDeltaEta = 0; iDeltaEta < knDeltaEtaBins; iDeltaEta++){
            fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt][iDeltaEta] = NULL; // DeltaPhi between jet and track
          }
          
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
 * Constructor
 */
DijetHistogramManager::DijetHistogramManager(TFile *inputFile) :
  DijetHistogramManager()
{
  fInputFile = inputFile;
  
  // Read card from inputfile and collision system from card
  fCard = new DijetCard(inputFile);
  TString collisionSystem = fCard->GetDataType();
  
  // Make a string for collision system based on information on the card
  fSystemAndEnergy = Form("%s 5.02 TeV",collisionSystem.Data());
  fCompactSystemAndEnergy = fSystemAndEnergy;
  fCompactSystemAndEnergy.ReplaceAll(" ","");
  fCompactSystemAndEnergy.ReplaceAll(".","v");
  
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
  fJffCorrectionFinder(in.fJffCorrectionFinder),
  fApplyJffCorrection(in.fApplyJffCorrection),
  fApplySpilloverCorrection(in.fApplySpilloverCorrection),
  fApplySeagullCorrection(in.fApplySeagullCorrection),
  fLoadEventInformation(in.fLoadEventInformation),
  fLoadDijetHistograms(in.fLoadDijetHistograms),
  fLoad2DHistograms(in.fLoad2DHistograms),
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
    fLowDeltaPhiBinBorders[iDeltaPhi] = in.fLowDeltaPhiBinBorders[iDeltaPhi];
    fHighDeltaPhiBinBorders[iDeltaPhi] = in.fHighDeltaPhiBinBorders[iDeltaPhi];
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
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt] = in.fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt]; // DeltaEta and deltaPhi between jet and track
          
          // Loop over deltaEta bins
          for(int iDeltaEta = 0; iDeltaEta < knDeltaEtaBins; iDeltaEta++){
            fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt][iDeltaEta] = in.fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentrality][iTrackPt][iDeltaEta];         // DeltaPhi between jet and track
          }
          
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
  delete fJffCorrectionFinder;
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
  
  // Helper variables
  int connectedIndex;
  double scalingFactor;
  
  // Loop over all jet-track correlation types and apply the mixed event correction
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    if(!fLoadJetTrackCorrelations[iJetTrack]) continue; // Only correct the histograms that are selected for analysis
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      
      // Scaling factor is number of dijets for leading and subleading and number of jets for inclusive correlations
      if(iJetTrack < kTrackInclusiveJet){
        scalingFactor = 1.0/GetPtIntegral(iCentralityBin);
      } else {
        scalingFactor = 1.0/GetInclusiveJetPtIntegral(iCentralityBin);
      }
      
      for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){

        connectedIndex = GetConnectedIndex(iJetTrack);
        
        // Do the mixed event correction for the current jet-track correlation histogram
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kCorrected][iCentralityBin][iTrackPtBin] = fMethods->MixedEventCorrect(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kSameEvent][iCentralityBin][iTrackPtBin],fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kMixedEvent][iCentralityBin][iTrackPtBin],fhJetTrackDeltaEtaDeltaPhi[connectedIndex][kMixedEvent][iCentralityBin][iTrackPtBin]);
        
        // Scale the histograms with the number of jets/dijets
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kCorrected][iCentralityBin][iTrackPtBin]->Scale(scalingFactor);
        
        // Apply the seagull correction after the mixed event correction
        if(fApplySeagullCorrection){
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kCorrected][iCentralityBin][iTrackPtBin] = fMethods->DoSeagullCorrection(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kCorrected][iCentralityBin][iTrackPtBin]);
        }

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
  int nProjectedBins;
  bool isInclusive;
  TH2D *correctionHistogram;
  
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    if(!fLoadJetTrackCorrelations[iJetTrack]) continue; // Only correct the histograms that are selected for analysis
    isInclusive = (iJetTrack >= kTrackInclusiveJet); // Set the flag for inclusive jet-track correlations
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){

        // Get the subleading/leading jet index connected to the currect leading/subleading correlation type
        connectedIndex = GetConnectedIndex(iJetTrack);
        
        // Subtract the background from the mixed event corrected histogram
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackgroundSubtracted][iCentralityBin][iTrackPtBin] = fMethods->SubtractBackground(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kCorrected][iCentralityBin][iTrackPtBin],fhJetTrackDeltaEtaDeltaPhi[connectedIndex][kCorrected][iCentralityBin][iTrackPtBin],isInclusive);
        
        // Get also the background and background overlap region for QA purposes
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackground][iCentralityBin][iTrackPtBin] = fMethods->GetBackground();
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackgroundOverlap][iCentralityBin][iTrackPtBin] = fMethods->GetBackgroundOverlap();
        
        // Apply the JFF correction to the background subtracted deltaEta-deltaPhi distribution
        if(fApplyJffCorrection && fJffCorrectionFinder->CorrectionReady()){
          correctionHistogram = fJffCorrectionFinder->GetDeltaEtaDeltaPhiJffCorrection(iJetTrack,iCentralityBin,iTrackPtBin);
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackgroundSubtracted][iCentralityBin][iTrackPtBin]->Add(correctionHistogram,-1);
        }
        
        // Apply the spillover correction to the background subtracted deltaEta-deltaPhi distribution
        if(fApplySpilloverCorrection && fJffCorrectionFinder->SpilloverReady()){
          correctionHistogram = fJffCorrectionFinder->GetDeltaEtaDeltaPhiSpilloverCorrection(iJetTrack,iCentralityBin,iTrackPtBin);
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackgroundSubtracted][iCentralityBin][iTrackPtBin]->Add(correctionHistogram,-1);
        }
        
        // Calculate the jet shape from the background subtracted histogram
        fhJetShape[kJetShape][iJetTrack][iCentralityBin][iTrackPtBin] = fMethods->GetJetShape(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackgroundSubtracted][iCentralityBin][iTrackPtBin]);
        
        // Get the number of two-dimensional histogram bins used for each deltaR bin in the jet shape histogram
        fhJetShape[kJetShapeBinCount][iJetTrack][iCentralityBin][iTrackPtBin] = fMethods->GetJetShapeCounts();
        
        // Get the mapping histogram of Rbins to deltaPhi-deltaEta bins
        fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kJetShapeBinMap][iCentralityBin][iTrackPtBin] = fMethods->GetJetShapeBinMap();
        
        // Project the deltaPhi and deltaEta histograms from the processed two-dimensional histograms
        for(int iCorrelationType = kCorrected; iCorrelationType < knCorrelationTypes; iCorrelationType++){
          
          // DeltaPhi histogram over whole eta
          sprintf(histogramName,"%sDeltaPhiProjection%d",fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->GetName(),kWholeEta);
          nBins = fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->GetNbinsY();
          fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin][kWholeEta] = (TH1D*) fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->ProjectionX(histogramName,1,nBins)->Clone();  // Exclude underflow and overflow bins by specifying range
          fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin][kWholeEta]->Scale(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->GetYaxis()->GetBinWidth(1));  // For correct normalization, need to divide out deltaEta bin width of the two-dimensional histogram
          
          // DeltaPhi histogram over signal eta region
          sprintf(histogramName,"%sDeltaPhiProjection%d",fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->GetName(),kSignalEtaRegion);
          fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin][kSignalEtaRegion] = (TH1D*)fMethods->ProjectSignalDeltaPhi(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin])->Clone();
          
          // DeltaPhi histogram over background eta region
          sprintf(histogramName,"%sDeltaPhiProjection%d",fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->GetName(),kBackgroundEtaRegion);
          fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin][kBackgroundEtaRegion] = (TH1D*)fMethods->ProjectBackgroundDeltaPhi(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin])->Clone();
          
          for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
            sprintf(histogramName,"%sDeltaEtaProjection%d",fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->GetName(),iDeltaPhi);
            fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaPhi] = (TH1D*) fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->ProjectionY(histogramName,fLowDeltaPhiBinIndices[iDeltaPhi],fHighDeltaPhiBinIndices[iDeltaPhi])->Clone();
            
            // To retain the normalization, we must scale the histograms with the number of bins projected over and by the width of deltaPhi bin
            nProjectedBins = fHighDeltaPhiBinIndices[iDeltaPhi] - fLowDeltaPhiBinIndices[iDeltaPhi] + 1;
            //fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaPhi]->Scale(1.0/nProjectedBins);
            fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaPhi]->Scale(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->GetXaxis()->GetBinWidth(1));
            
          } // DeltaPhi loop
        } // Correlation type loop
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation category loop
}

/*
 *  Apply the JFF correction to relevant histograms
 *
 *  Arguments:
 *   JffCorrector *jffCorrectionFinder = Class holding the JFF correction histograms
 */
void DijetHistogramManager::ApplyJffCorrection(JffCorrector *jffCorrectionFinder){
  
  // Helper histogram for reading the JFF correction
  TH1D *correctionHistogram;
  
  // Loop over all the histogram and apply the JFF correction
  double scalingFactor;
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    if(!fLoadJetTrackCorrelations[iJetTrack]) continue; // Only scale the histograms that are selected for analysis
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      
      // Scaling factor is number of dijets for leading and subleading and number of jets for inclusive correlations
      if(iJetTrack < kTrackInclusiveJet){
        scalingFactor = 1.0/GetPtIntegral(iCentralityBin);
      } else {
        scalingFactor = 1.0/GetInclusiveJetPtIntegral(iCentralityBin);
      }
      
      for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
        correctionHistogram = jffCorrectionFinder->GetJetShapeJffCorrection(iJetTrack,iCentralityBin,iTrackPtBin);
        fhJetShape[kJetShape][iJetTrack][iCentralityBin][iTrackPtBin]->Scale(scalingFactor);  // Need to scale with the number of dijets/all jets since the correction is also normalized to the number of dijets/all jets
        fhJetShape[kJetShape][iJetTrack][iCentralityBin][iTrackPtBin]->Add(correctionHistogram,-1);
      } // Track pT loop
    } // Centrality loop
  } // Jet-track correlation category loop
}

/*
 * Normalize the jet shape histograms such that the pT integrated result is unity in the range deltaR < 1
 */
void DijetHistogramManager::NormalizeJetShape(){
  
  // Helper variables for doing the normalization
  TH1D *jetShapeSum;
  double jetShapeIntegral;
  
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    if(!fLoadJetTrackCorrelations[iJetTrack]) continue; // Only scale the histograms that are selected for analysis
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      jetShapeSum = (TH1D*)fhJetShape[kJetShape][iJetTrack][iCentralityBin][fFirstLoadedTrackPtBin]->Clone(Form("jetShapeSum%d%d",iJetTrack,iCentralityBin));
      
      // First, sum all pT bins together
      for(int iTrackPtBin = fFirstLoadedTrackPtBin+1; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
        jetShapeSum->Add(fhJetShape[kJetShape][iJetTrack][iCentralityBin][iTrackPtBin]);
      } // Track pT loop
      
      // Then calculate the integral for deltaR < 1
      jetShapeIntegral = jetShapeSum->Integral(1,jetShapeSum->FindBin(0.99),"width");
      
      // Finally, normalize each pT bin with the intagral
      for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
        fhJetShape[kJetShape][iJetTrack][iCentralityBin][iTrackPtBin]->Scale(1.0/jetShapeIntegral);
      } // Track pT loop
      
    } // Centrality loop
  } // Jet-track correlation category loop
  
}

/*
 * Get the index that of the same type of leading/subleading jet-track correlation as the
 * given subleading/leading jet-track correlation index.
 */
int DijetHistogramManager::GetConnectedIndex(const int jetTrackIndex) const{
  if(jetTrackIndex < kTrackSubleadingJet) return jetTrackIndex + 3;  // Three leading jet histograms are listed before subleading jet histograms
  if(jetTrackIndex < kTrackInclusiveJet) return jetTrackIndex - 3;   // Three subleading jet histograms are listed after leading jet histograms
  return jetTrackIndex;  // There is no connected index for inclusive jet, so return the index itself
}

/*
 * Load all the selected histograms from the inputfile
 */
void DijetHistogramManager::LoadHistograms(){
  
  // Always load the number of events histogram
  fhEvents = (TH1D*) fInputFile->Get("nEvents");                           // Number of events surviving different event cuts
  
  // Load the event information histograms
  if(fLoadEventInformation){
    fhVertexZ = (TH1D*) fInputFile->Get("vertexZ");                        // Vertex z position
    fhVertexZWeighted = (TH1D*) fInputFile->Get("vertexZweighted");        // MC weighted vertex z position
    fhTrackCuts = (TH1D*) fInputFile->Get("trackCuts");                    // Number of tracks surviving different track cuts
    fhTrackCutsInclusive = (TH1D*) fInputFile->Get("trackCutsInclusive");  // Number of inclusive tracks surviving different track cuts
    fhCentrality = (TH1D*) fInputFile->Get("centrality");                  // Centrality in all events
    fhCentralityWeighted = (TH1D*) fInputFile->Get("centralityWeighted");  // MC weighted centrality in all events
    fhCentralityDijet = (TH1D*) fInputFile->Get("centralityDijet");        // Centrality in dijet events
    fhPtHat = (TH1D*) fInputFile->Get("pthat");                            // pT hat for MC events
    fhPtHatWeighted = (TH1D*) fInputFile->Get("pthatWeighted");            // Weighted pT hat for MC events
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
    
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      
      // Select the bin indices
      if(iCentralityBin == fLastLoadedCentralityBin) {
        duplicateRemoverCentrality = 0;
      } else {
        duplicateRemoverCentrality = -1;
      }
      lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
      higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
      
      // Always load single jet pT histograms
      fhJetPt[iJetCategory][iCentralityBin] = FindHistogram(fInputFile,fSingleJetHistogramName[iJetCategory],0,centralityIndex[iJetCategory],lowerCentralityBin,higherCentralityBin);
      
      if(!fLoadSingleJets[iJetCategory]) continue;  // Only load the remaining single jet histograms is selected
      
      fhJetPhi[iJetCategory][iCentralityBin] = FindHistogram(fInputFile,fSingleJetHistogramName[iJetCategory],1,centralityIndex[iJetCategory],lowerCentralityBin,higherCentralityBin);
      fhJetEta[iJetCategory][iCentralityBin] = FindHistogram(fInputFile,fSingleJetHistogramName[iJetCategory],2,centralityIndex[iJetCategory],lowerCentralityBin,higherCentralityBin);
      if(fLoad2DHistograms) fhJetEtaPhi[iJetCategory][iCentralityBin] = FindHistogram2D(fInputFile,fSingleJetHistogramName[iJetCategory],1,2,centralityIndex[iJetCategory],lowerCentralityBin,higherCentralityBin);
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
  
  if(!fLoadDijetHistograms) return; // Do not load the histograms if they are not selected for loading
  
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
    if(fLoad2DHistograms) fhDijetLeadingVsSubleadingPt[iCentralityBin] = FindHistogram2D(fInputFile,"dijet",0,1,4,lowerCentralityBin,higherCentralityBin);
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
  int axisAdder = 0;
  
  // Loop over all track histograms
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    if(!fLoadTracks[iTrackType]) continue;  // Only load the selected track types
    for(int iCorrelationType = 0; iCorrelationType <= kMixedEvent; iCorrelationType++){  // Data file contains only same and mixed event distributions
      if(iTrackType > kUncorrectedTrack && iCorrelationType == kMixedEvent) continue; // No mixed event histograms for inclusive tracks
      if(iTrackType > kUncorrectedTrack){
        axisAdder = -1;
      } else {
        axisAdder = 0;
      }
      for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
        
        // Select the bin indices
        if(iCentralityBin == fLastLoadedCentralityBin) {
          duplicateRemoverCentrality = 0;
        } else {
          duplicateRemoverCentrality = -1;
        }
        lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
        higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
        
        // Setup axes with restrictions, (3 = centrality, 4 = correlation type)
        axisIndices[0] = 3; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;
        axisIndices[1] = 4; lowLimits[1] = iCorrelationType+1; highLimits[1] = iCorrelationType+1;
        
        fhTrackPt[iTrackType][iCorrelationType][iCentralityBin] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],0,2+axisAdder,axisIndices,lowLimits,highLimits);
        fhTrackPhi[iTrackType][iCorrelationType][iCentralityBin][knTrackPtBins] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],1,2+axisAdder,axisIndices,lowLimits,highLimits);
        fhTrackEta[iTrackType][iCorrelationType][iCentralityBin][knTrackPtBins] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],2,2+axisAdder,axisIndices,lowLimits,highLimits);
        if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCorrelationType][iCentralityBin][knTrackPtBins] = FindHistogram2D(fInputFile,fTrackHistogramNames[iTrackType],1,2,2+axisAdder,axisIndices,lowLimits,highLimits);
        
        for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
          
          // Select the bin indices for track pT
          lowerTrackPtBin = fFineTrackPtBinIndices[iTrackPtBin];
          higherTrackPtBin = fFineTrackPtBinIndices[iTrackPtBin+1]+duplicateRemoverTrackPt;
          
          // Add restriction for pT axis (0)
          axisIndices[2+axisAdder] = 0; lowLimits[2+axisAdder] = lowerTrackPtBin; highLimits[2+axisAdder] = higherTrackPtBin;
          
          // Read the angle histograms in track pT bins
          fhTrackPhi[iTrackType][iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],1,3+axisAdder,axisIndices,lowLimits,highLimits);
          fhTrackEta[iTrackType][iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],2,3+axisAdder,axisIndices,lowLimits,highLimits);
          if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram2D(fInputFile,fTrackHistogramNames[iTrackType],1,2,3+axisAdder,axisIndices,lowLimits,highLimits);
          
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
        if(iCentralityBin == fLastLoadedCentralityBin) {
          duplicateRemoverCentrality = 0;
        } else {
          duplicateRemoverCentrality = -1;
        }
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
          
          // Indexing is different for inclusive jet-track correlation, as they do not have dijet asymmetry
          if(iJetTrack >= kTrackInclusiveJet){
            axisIndices[0] = 4;
            axisIndices[1] = 3;
          }
          
          fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin][kWholeEta] = FindHistogram(fInputFile,fJetTrackHistogramNames[iJetTrack],1,3,axisIndices,lowLimits,highLimits);
          if(fLoad2DHistograms) fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin] = FindHistogram2D(fInputFile,fJetTrackHistogramNames[iJetTrack],1,2,3,axisIndices,lowLimits,highLimits);
          
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
  
  // Create a unique name for each histogram that is read from the file
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
 * Write all the loaded histograms into a file
 *
 *  const char* fileName = Name of the file to which the histograms are written
 *  const char* fileOption = Option given to the file when it is loaded
 */
void DijetHistogramManager::Write(const char* fileName, const char* fileOption){
  
  // Create the output file
  TFile *outputFile = new TFile(fileName,fileOption);
  
  // Helper variable for renaming the saved histograms
  char histogramNamer[200];
  
  // Write the event information histograms to the output file
  if(fLoadEventInformation){
    fhEvents->Write();             // Number of events surviving different event cuts
    fhVertexZ->Write();            // Vertex z position
    fhVertexZWeighted->Write();    // MC weighted vertex z position
    fhTrackCuts->Write();          // Number of tracks surviving different track cuts
    fhTrackCutsInclusive->Write(); // Number of inclusive tracks surviving different track cuts
    fhCentrality->Write();         // Centrality in all events
    fhCentralityWeighted->Write(); // MC weighted centrality in all events
    fhCentralityDijet->Write();    // Centrality in dijet events
    fhPtHat->Write();              // pT hat for MC events (only meaningful for MC)
    fhPtHatWeighted->Write();      // Weighted pT hat distribution (only meaningful for MC)
  }
 
  // Write the single jet histograms to the output file
  for(int iJetCategory = 0; iJetCategory < knSingleJetCategories; iJetCategory++){
    if(!fLoadSingleJets[iJetCategory]) continue;  // Only write the loaded the selected histograms
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory(fSingleJetHistogramName[iJetCategory])) gDirectory->mkdir(fSingleJetHistogramName[iJetCategory]);
    gDirectory->cd(fSingleJetHistogramName[iJetCategory]);
    
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      
      // Single jet pT
      sprintf(histogramNamer,"%sPt_C%d",fSingleJetHistogramName[iJetCategory],iCentralityBin);
      fhJetPt[iJetCategory][iCentralityBin]->Write(histogramNamer);
      
      // Single jet phi
      sprintf(histogramNamer,"%sPhi_C%d",fSingleJetHistogramName[iJetCategory],iCentralityBin);
      fhJetPhi[iJetCategory][iCentralityBin]->Write(histogramNamer);
      
      // Single jet eta
      sprintf(histogramNamer,"%sEta_C%d",fSingleJetHistogramName[iJetCategory],iCentralityBin);
      fhJetEta[iJetCategory][iCentralityBin]->Write(histogramNamer);
      
      //Single jet eta-phi
      sprintf(histogramNamer,"%sEtaPhi_C%d",fSingleJetHistogramName[iJetCategory],iCentralityBin);
      if(fLoad2DHistograms) fhJetEtaPhi[iJetCategory][iCentralityBin]->Write(histogramNamer);
    } // Loop over centrality bins
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Loop over single jet categories
  
  // Write the dijet histograms to the output file
  if(fLoadDijetHistograms){
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory("dijet")) gDirectory->mkdir("dijet");
    gDirectory->cd("dijet");
    
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      
      // Dijet deltaPhi
      sprintf(histogramNamer,"dijetDeltaPhi_C%d",iCentralityBin);
      fhDijetDphi[iCentralityBin]->Write(histogramNamer);
      
      // Dijet asymmetry
      sprintf(histogramNamer,"dijetAsymmetry_C%d",iCentralityBin);
      fhDijetAsymmetry[iCentralityBin]->Write(histogramNamer);
      
      // Leading jet pT vs. subleading jet pT
      sprintf(histogramNamer,"leadingVsSubleadingPt_C%d",iCentralityBin);
      if(fLoad2DHistograms) fhDijetLeadingVsSubleadingPt[iCentralityBin]->Write(histogramNamer);
    }
    
    // Return back to main directory
    gDirectory->cd("../");
  }
  
  // Write the track histograms to the output file
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    if(!fLoadTracks[iTrackType]) continue;  // Only write the loaded track types
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory(fTrackHistogramNames[iTrackType])) gDirectory->mkdir(fTrackHistogramNames[iTrackType]);
    gDirectory->cd(fTrackHistogramNames[iTrackType]);
    
    for(int iCorrelationType = 0; iCorrelationType <= kMixedEvent; iCorrelationType++){  // Tracks have only same and mixed event distributions
      if(iTrackType > kUncorrectedTrack && iCorrelationType == kMixedEvent) continue; // No mixed event histograms for inclusive tracks

      for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
        
        // Track pT
        sprintf(histogramNamer,"%sPt%s_C%d",fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin);
        fhTrackPt[iTrackType][iCorrelationType][iCentralityBin]->Write(histogramNamer);
        
        // pT integrated track phi
        sprintf(histogramNamer,"%sPhi%s_C%dT%d",fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,knTrackPtBins);
        fhTrackPhi[iTrackType][iCorrelationType][iCentralityBin][knTrackPtBins]->Write(histogramNamer);
        
        // pT integrated track eta
        sprintf(histogramNamer,"%sEta%s_C%dT%d",fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,knTrackPtBins);
        fhTrackEta[iTrackType][iCorrelationType][iCentralityBin][knTrackPtBins]->Write(histogramNamer);
        
        // pT integrated track eta-phi
        sprintf(histogramNamer,"%sEtaPhi%s_C%dT%d",fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,knTrackPtBins);
        if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCorrelationType][iCentralityBin][knTrackPtBins]->Write(histogramNamer);
        
        for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
          
          // Track phi in track pT bins
          sprintf(histogramNamer,"%sPhi%s_C%dT%d",fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,iTrackPtBin);
          fhTrackPhi[iTrackType][iCorrelationType][iCentralityBin][iTrackPtBin]->Write(histogramNamer);
          
          // Track eta in track pT bins
          sprintf(histogramNamer,"%sEta%s_C%dT%d",fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,iTrackPtBin);
          fhTrackEta[iTrackType][iCorrelationType][iCentralityBin][iTrackPtBin]->Write(histogramNamer);
          
          // Track eta-phi in track pT bins
          sprintf(histogramNamer,"%sEtaPhi%s_C%dT%d",fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,iTrackPtBin);
          if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCorrelationType][iCentralityBin][iTrackPtBin]->Write(histogramNamer);
          
        } // Track pT loop
      } // Centrality loop
    } // Correlation type loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Track category loop
  
  // Write the jet-track correlation histograms to the output file
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    if(!fLoadJetTrackCorrelations[iJetTrack]) continue; // Only write categories of correlation that are loaded
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory(fJetTrackHistogramNames[iJetTrack])) gDirectory->mkdir(fJetTrackHistogramNames[iJetTrack]);
    gDirectory->cd(fJetTrackHistogramNames[iJetTrack]);
    
    // Loop over correlation types
    for(int iCorrelationType = 0; iCorrelationType < knCorrelationTypes; iCorrelationType++){
      
      // Loop over centrality bins
      for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
        
        // Loop over track pT bins
        for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
          
          // Jet-track deltaPhi
          for(int iDeltaEta = 0; iDeltaEta < knDeltaEtaBins; iDeltaEta++){
           
            if(iDeltaEta > kWholeEta && iCorrelationType < kCorrected) continue; // DeltaEta slicing not implemented for same and mixed event
            sprintf(histogramNamer,"%sDeltaPhi%s%s_C%dT%d",fJetTrackHistogramNames[iJetTrack],fCompactCorrelationTypeString[iCorrelationType].Data(),fCompactDeltaEtaString[iDeltaEta],iCentralityBin,iTrackPtBin);
          fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaEta]->Write(histogramNamer);
          }
          
          // Jet-track deltaEtaDeltaPhi
          sprintf(histogramNamer,"%sDeltaEtaDeltaPhi%s_C%dT%d",fJetTrackHistogramNames[iJetTrack],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,iTrackPtBin);
          if(fLoad2DHistograms) fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin]->Write(histogramNamer);
          
          // DeltaPhi binning for deltaEta histogram
          for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
            sprintf(histogramNamer,"%sDeltaEta%s%s_C%dT%d",fJetTrackHistogramNames[iJetTrack],fCompactCorrelationTypeString[iCorrelationType].Data(),fCompactDeltaPhiString[iDeltaPhi].Data(),iCentralityBin,iTrackPtBin);
            fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaPhi]->Write(histogramNamer);
          } // DeltaPhi loop
        } // Track pT loop
      } // Centrality loop
    } // Correlation type loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Jet-track correlation category loop
  
  // Write the jet shape histograms to the output file
  for(int iJetShape = 0; iJetShape < knJetShapeTypes; iJetShape++){
    
    // Loop over jet-track correlation categories
    for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
      if(!fLoadJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
      
      // Create a directory for the histograms if it does not already exist
      sprintf(histogramNamer,"%s_%s",fJetShapeHistogramName[iJetShape],fJetTrackHistogramNames[iJetTrack]);
      if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
      gDirectory->cd(histogramNamer);
      
      // Loop over centrality
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        
        // Loop over track pT bins
        for(int iTrackPt = fFirstLoadedTrackPtBin; iTrackPt <= fLastLoadedTrackPtBin; iTrackPt++){
          sprintf(histogramNamer,"%s_%s_C%dT%d",fJetShapeHistogramName[iJetShape],fJetTrackHistogramNames[iJetTrack],iCentrality,iTrackPt);
          fhJetShape[iJetShape][iJetTrack][iCentrality][iTrackPt]->Write(histogramNamer);
          
        } // Track pT loop
      } // Centrality loop
      
      // Return back to main directory
      gDirectory->cd("../");
      
    } // Jet-track correlation category loop
    
  } // Jet shape type loop
  
  // Write the card to the output file if it is not already written
  if(!gDirectory->GetDirectory("JCard")) fCard->Write(outputFile);
  
  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;
  
}

/*
 * Load the selected histograms from a file containing readily processed histograms
 */
void DijetHistogramManager::LoadProcessedHistograms(){
  
  // Helper variable for finding names of loaded histograms
  char histogramNamer[200];
  
  // Always load the number of events histogram
  fhEvents = (TH1D*) fInputFile->Get("nEvents");                           // Number of events surviving different event cuts
  
  // Load the event information histograms
  if(fLoadEventInformation){
    fhVertexZ = (TH1D*) fInputFile->Get("vertexZ");                        // Vertex z position
    fhVertexZWeighted = (TH1D*) fInputFile->Get("vertexZweighted");        // MC weighted vertex z position
    fhTrackCuts = (TH1D*) fInputFile->Get("trackCuts");                    // Number of tracks surviving different track cuts
    fhTrackCutsInclusive = (TH1D*) fInputFile->Get("trackCutsInclusive");  // Number of inclusive tracks surviving different track cuts
    fhCentrality = (TH1D*) fInputFile->Get("centrality");                  // Centrality in all events
    fhCentralityWeighted = (TH1D*) fInputFile->Get("centralityWeighted");  // MC weighted centrality in all events
    fhCentralityDijet = (TH1D*) fInputFile->Get("centralityDijet");        // Centrality in dijet events
    fhPtHat = (TH1D*) fInputFile->Get("pthat");                            // pT hat for MC events
    fhPtHatWeighted = (TH1D*) fInputFile->Get("pthatWeighted");            // Weighted pT hat for MC events
  }
  
  // Load the single jet histograms from the input file
  for(int iJetCategory = 0; iJetCategory < knSingleJetCategories; iJetCategory++){
    
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      
      // Always load single jet pT histograms
      sprintf(histogramNamer,"%s/%sPt_C%d",fSingleJetHistogramName[iJetCategory],fSingleJetHistogramName[iJetCategory],iCentralityBin);
      fhJetPt[iJetCategory][iCentralityBin] = (TH1D*) fInputFile->Get(histogramNamer);
      
      if(!fLoadSingleJets[iJetCategory]) continue;  // Only load the loaded the selected histograms
      
      // Single jet phi
      sprintf(histogramNamer,"%s/%sPhi_C%d",fSingleJetHistogramName[iJetCategory],fSingleJetHistogramName[iJetCategory],iCentralityBin);
      fhJetPhi[iJetCategory][iCentralityBin] = (TH1D*) fInputFile->Get(histogramNamer);
      
      // Single jet eta
      sprintf(histogramNamer,"%s/%sEta_C%d",fSingleJetHistogramName[iJetCategory],fSingleJetHistogramName[iJetCategory],iCentralityBin);
      fhJetEta[iJetCategory][iCentralityBin] = (TH1D*) fInputFile->Get(histogramNamer);
      
      //Single jet eta-phi
      sprintf(histogramNamer,"%s/%sEtaPhi_C%d",fSingleJetHistogramName[iJetCategory],fSingleJetHistogramName[iJetCategory],iCentralityBin);
      if(fLoad2DHistograms) fhJetEtaPhi[iJetCategory][iCentralityBin] = (TH2D*) fInputFile->Get(histogramNamer);
    } // Loop over centrality bins
    
  } // Loop over single jet categories
  
  // Load the dijet histograms from the input file
  if(fLoadDijetHistograms){
    
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      
      // Dijet deltaPhi
      sprintf(histogramNamer,"dijet/dijetDeltaPhi_C%d",iCentralityBin);
      fhDijetDphi[iCentralityBin] = (TH1D*) fInputFile->Get(histogramNamer);
      
      // Dijet asymmetry
      sprintf(histogramNamer,"dijet/dijetAsymmetry_C%d",iCentralityBin);
      fhDijetAsymmetry[iCentralityBin] = (TH1D*) fInputFile->Get(histogramNamer);
      
      // Leading jet pT vs. subleading jet pT
      sprintf(histogramNamer,"dijet/leadingVsSubleadingPt_C%d",iCentralityBin);
      if(fLoad2DHistograms) fhDijetLeadingVsSubleadingPt[iCentralityBin] = (TH2D*) fInputFile->Get(histogramNamer);
    }

  }
  
  // Load the track histograms from the input file
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    if(!fLoadTracks[iTrackType]) continue;  // Only load the selected track types
    
    for(int iCorrelationType = 0; iCorrelationType <= kMixedEvent; iCorrelationType++){  // Tracks have only same and mixed event distributions
      if(iTrackType > kUncorrectedTrack && iCorrelationType == kMixedEvent) continue; // No mixed event histograms for inclusive tracks
      
      for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
        
        // Track pT
        sprintf(histogramNamer,"%s/%sPt%s_C%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin);
        fhTrackPt[iTrackType][iCorrelationType][iCentralityBin] = (TH1D*) fInputFile->Get(histogramNamer);
        
        // pT integrated track phi
        sprintf(histogramNamer,"%s/%sPhi%s_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,knTrackPtBins);
        fhTrackPhi[iTrackType][iCorrelationType][iCentralityBin][knTrackPtBins] = (TH1D*) fInputFile->Get(histogramNamer);
        
        // pT integrated track eta
        sprintf(histogramNamer,"%s/%sEta%s_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,knTrackPtBins);
        fhTrackEta[iTrackType][iCorrelationType][iCentralityBin][knTrackPtBins] = (TH1D*) fInputFile->Get(histogramNamer);
        
        // pT integrated track eta-phi
        sprintf(histogramNamer,"%s/%sEtaPhi%s_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,knTrackPtBins);
        if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCorrelationType][iCentralityBin][knTrackPtBins] = (TH2D*) fInputFile->Get(histogramNamer);
        
        for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
          
          // Track phi in track pT bins
          sprintf(histogramNamer,"%s/%sPhi%s_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,iTrackPtBin);
          fhTrackPhi[iTrackType][iCorrelationType][iCentralityBin][iTrackPtBin] = (TH1D*) fInputFile->Get(histogramNamer);
          
          // Track eta in track pT bins
          sprintf(histogramNamer,"%s/%sEta%s_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,iTrackPtBin);
          fhTrackEta[iTrackType][iCorrelationType][iCentralityBin][iTrackPtBin] = (TH1D*) fInputFile->Get(histogramNamer);
          
          // Track eta-phi in track pT bins
          sprintf(histogramNamer,"%s/%sEtaPhi%s_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,iTrackPtBin);
          if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCorrelationType][iCentralityBin][iTrackPtBin] = (TH2D*) fInputFile->Get(histogramNamer);
          
        } // Track pT loop
      } // Centrality loop
    } // Correlation type loop
  } // Track category loop
  
  // Load the jet-track correlation histograms from the input file
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    if(!fLoadJetTrackCorrelations[iJetTrack]) continue; // Only load categories of correlation that are selected
    
    // Loop over correlation types
    for(int iCorrelationType = 0; iCorrelationType < knCorrelationTypes; iCorrelationType++){
      
      // Loop over centrality bins
      for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
        
        // Loop over track pT bins
        for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
          
          // Jet-track deltaPhi
          for(int iDeltaEta = 0; iDeltaEta < knDeltaEtaBins; iDeltaEta++){
            
            if(iDeltaEta > kWholeEta && iCorrelationType < kCorrected) continue; // DeltaEta slicing not implemented for same and mixed event
            
            sprintf(histogramNamer,"%s/%sDeltaPhi%s%s_C%dT%d",fJetTrackHistogramNames[iJetTrack],fJetTrackHistogramNames[iJetTrack],fCompactCorrelationTypeString[iCorrelationType].Data(),fCompactDeltaEtaString[iDeltaEta],iCentralityBin,iTrackPtBin);
            fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaEta] = (TH1D*) fInputFile->Get(histogramNamer);
          }
          
          // Jet-track deltaEtaDeltaPhi
          sprintf(histogramNamer,"%s/%sDeltaEtaDeltaPhi%s_C%dT%d",fJetTrackHistogramNames[iJetTrack],fJetTrackHistogramNames[iJetTrack],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,iTrackPtBin);
          if(fLoad2DHistograms) fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin] = (TH2D*) fInputFile->Get(histogramNamer);
          
          // DeltaPhi binning for deltaEta histogram
          for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
            sprintf(histogramNamer,"%s/%sDeltaEta%s%s_C%dT%d",fJetTrackHistogramNames[iJetTrack],fJetTrackHistogramNames[iJetTrack],fCompactCorrelationTypeString[iCorrelationType].Data(),fCompactDeltaPhiString[iDeltaPhi].Data(),iCentralityBin,iTrackPtBin);
            fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iCentralityBin][iTrackPtBin][iDeltaPhi] = (TH1D*) fInputFile->Get(histogramNamer);
          } // DeltaPhi loop
        } // Track pT loop
      } // Centrality loop
    } // Correlation type loop
  } // Jet-track correlation category loop
  
  // Write the jet shape histograms to the output file
  for(int iJetShape = 0; iJetShape < knJetShapeTypes; iJetShape++){
    
    // Loop over jet-track correlation categories
    for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
      if(!fLoadJetTrackCorrelations[iJetTrack]) continue;  // Only load the selected categories
      
      // Loop over centrality
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        
        // Loop over track pT bins
        for(int iTrackPt = fFirstLoadedTrackPtBin; iTrackPt <= fLastLoadedTrackPtBin; iTrackPt++){
          sprintf(histogramNamer,"%s_%s/%s_%s_C%dT%d",fJetShapeHistogramName[iJetShape],fJetTrackHistogramNames[iJetTrack],fJetShapeHistogramName[iJetShape],fJetTrackHistogramNames[iJetTrack],iCentrality,iTrackPt);
          fhJetShape[iJetShape][iJetTrack][iCentrality][iTrackPt] = (TH1D*) fInputFile->Get(histogramNamer);
          
        } // Track pT loop
      } // Centrality loop
    } // Jet-track correlation category loop
  } // Jet shape type loop

}

/*
 * Read the bin indices for given bin borders
 *
 *  Arguments:
 *   const int nBins = Number of bins for the indices
 *   int *binIndices = Array of integers to be filled with bin index information read from the file
 *   const double *binBorders = Array for bin borders that are searched from the file
 *   const int iAxis = Index of the axis used for reading bin indices
 *   const bool setIndices = true: Set both bin indices and bin borders, false: Set only bin borders
 */
void DijetHistogramManager::SetBinIndices(const int nBins, double *copyBinBorders, int *binIndices, const double *binBorders, const int iAxis, const bool setIndices){
  TH1D* hBinner;
  if(setIndices) hBinner = FindHistogram(fInputFile,"trackLeadingJet",iAxis,0,0,0);
  for(int iBin = 0; iBin < nBins+1; iBin++){
    copyBinBorders[iBin] = binBorders[iBin];
    if(setIndices) binIndices[iBin] = hBinner->GetXaxis()->FindBin(binBorders[iBin]);
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
 * Set up centrality bin borders and indices according to provided bin borders
 *
 *
 */
void DijetHistogramManager::SetCentralityBins(const double *binBorders, const bool setIndices){
  SetBinIndices(knCentralityBins,fCentralityBinBorders,fCentralityBinIndices,binBorders,4,setIndices);
}

/*
 * Set up track pT bin indices according to provided bin borders
 */
void DijetHistogramManager::SetTrackPtBins(const double *binBorders, const bool setIndices){
  SetBinIndices(knTrackPtBins,fTrackPtBinBorders,fTrackPtBinIndices,binBorders,0,setIndices);
  
  // The track histograms have finer pT binning, so we need to use different bin indices for them
  if(setIndices){
    TH1D* hTrackPtBinner = FindHistogram(fInputFile,"track",0,0,0,0);
    for(int iTrackPt = 0; iTrackPt < knTrackPtBins+1; iTrackPt++){
      fFineTrackPtBinIndices[iTrackPt] = hTrackPtBinner->GetXaxis()->FindBin(binBorders[iTrackPt]);
    }
  }
}

/*
 * Set up deltaPhi bin indices according to provided bin borders
 */
void DijetHistogramManager::SetDeltaPhiBins(const double *lowBinBorders, const double *highBinBorders, TString deltaPhiStrings[knDeltaPhiBins], TString compactDeltaPhiStrings[knDeltaPhiBins], const bool setIndices){
  if(setIndices) SetBinIndices(knDeltaPhiBins,fLowDeltaPhiBinIndices,fHighDeltaPhiBinIndices,lowBinBorders,highBinBorders,1);
  for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
    fDeltaPhiString[iDeltaPhi] = deltaPhiStrings[iDeltaPhi];
    fCompactDeltaPhiString[iDeltaPhi] = compactDeltaPhiStrings[iDeltaPhi];
    fLowDeltaPhiBinBorders[iDeltaPhi] = lowBinBorders[iDeltaPhi];
    fHighDeltaPhiBinBorders[iDeltaPhi] = highBinBorders[iDeltaPhi];
  }
}

// Setter for loading event information
void DijetHistogramManager::SetLoadEventInformation(const bool loadOrNot){
  fLoadEventInformation = loadOrNot;
}

// Setter for loading dijet histograms
void DijetHistogramManager::SetLoadDijetHistograms(const bool loadOrNot){
  fLoadDijetHistograms = loadOrNot;
}

// Setter for loading leading jet histograms
void DijetHistogramManager::SetLoadLeadingJetHistograms(const bool loadOrNot){
  fLoadSingleJets[kLeadingJet] = loadOrNot;
}

// Setter for loading subleading jet histograms
void DijetHistogramManager::SetLoadSubleadingJetHistograms(const bool loadOrNot){
  fLoadSingleJets[kSubleadingJet] = loadOrNot;
}

// Setter for loading all jet histograms
void DijetHistogramManager::SetLoadAnyJetHistograms(const bool loadOrNot){
  fLoadSingleJets[kAnyJet] = loadOrNot;
}

// Setter for loading jet histograms
void DijetHistogramManager::SetLoadAllJets(const bool drawLeading, const bool drawSubleading, const bool drawAny){
  SetLoadLeadingJetHistograms(drawLeading);
  SetLoadSubleadingJetHistograms(drawSubleading);
  SetLoadAnyJetHistograms(drawAny);
}

// Setter for loading tracks
void DijetHistogramManager::SetLoadTracks(const bool loadOrNot){
  fLoadTracks[kTrack] = loadOrNot;
}

// Setter for loading uncorrected tracks
void DijetHistogramManager::SetLoadTracksUncorrected(const bool loadOrNot){
  fLoadTracks[kUncorrectedTrack] = loadOrNot;
}

// Setter for loading all track histograms
void DijetHistogramManager::SetLoadAllTracks(const bool drawTracks, const bool drawUncorrected){
  SetLoadTracks(drawTracks);
  SetLoadTracksUncorrected(drawUncorrected);
}

// Setter for loading inclusive tracks
void DijetHistogramManager::SetLoadInclusiveTracks(const bool loadOrNot){
  fLoadTracks[kInclusiveTrack] = loadOrNot;
}

// Setter for loading uncorrected inclusive tracks
void DijetHistogramManager::SetLoadInclusiveTracksUncorrected(const bool loadOrNot){
  fLoadTracks[kUncorrectedInclusiveTrack] = loadOrNot;
}

// Setter for loading all inclusive track histograms
void DijetHistogramManager::SetLoadAllInclusiveTracks(const bool drawTracks, const bool drawUncorrected){
  SetLoadInclusiveTracks(drawTracks);
  SetLoadInclusiveTracksUncorrected(drawUncorrected);
}

/*
 * Setter for loading jet-track correlations.
 *
 * The method sets the loading of the histogram defined by the primaryIndex.
 * The background subtraction histograms are conntructed using buth leading and subleading jet histograms.
 * Thus when loading the leading/subleading jet histograms, we need to also load the subleading/leading jet
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

// Setter for loading leading jet-track correlations
void DijetHistogramManager::SetLoadTrackLeadingJetCorrelations(const bool loadOrNot){
  SetLoadJetTrackCorrelations(loadOrNot,kTrackLeadingJet,kTrackSubleadingJet);
}

// Setter for loading uncorrected leading jet-track correlations
void DijetHistogramManager::SetLoadTrackLeadingJetCorrelationsUncorrected(const bool loadOrNot){
  SetLoadJetTrackCorrelations(loadOrNot,kUncorrectedTrackLeadingJet,kUncorrectedTrackSubleadingJet);
}

// Setter for loading pT weighted leading jet-track correlations
void DijetHistogramManager::SetLoadTrackLeadingJetCorrelationsPtWeighted(const bool loadOrNot){
  SetLoadJetTrackCorrelations(loadOrNot,kPtWeightedTrackLeadingJet,kPtWeightedTrackSubleadingJet);
}

// Setter for loading all correlations related to tracks and leading jets
void DijetHistogramManager::SetLoadAllTrackLeadingJetCorrelations(const bool drawLeading, const bool drawUncorrected, const bool drawPtWeighted){
  SetLoadTrackLeadingJetCorrelations(drawLeading);
  SetLoadTrackLeadingJetCorrelationsUncorrected(drawUncorrected);
  SetLoadTrackLeadingJetCorrelationsPtWeighted(drawPtWeighted);
}

// Setter for loading subleading jet-track correlations
void DijetHistogramManager::SetLoadTrackSubleadingJetCorrelations(const bool loadOrNot){
  SetLoadJetTrackCorrelations(loadOrNot,kTrackSubleadingJet,kTrackLeadingJet);
}

// Setter for loading uncorrected subleading jet-track correlations
void DijetHistogramManager::SetLoadTrackSubleadingJetCorrelationsUncorrected(const bool loadOrNot){
  SetLoadJetTrackCorrelations(loadOrNot,kUncorrectedTrackSubleadingJet,kUncorrectedTrackLeadingJet);
}

// Setter for loading pT weighted subleading jet-track correlations
void DijetHistogramManager::SetLoadTrackSubleadingJetCorrelationsPtWeighted(const bool loadOrNot){
  SetLoadJetTrackCorrelations(loadOrNot,kPtWeightedTrackSubleadingJet,kPtWeightedTrackLeadingJet);
}

// Setter for loading all correlations related to tracks and subleading jets
void DijetHistogramManager::SetLoadAllTrackSubleadingJetCorrelations(const bool drawSubleading, const bool drawUncorrected, const bool drawPtWeighted){
  SetLoadTrackSubleadingJetCorrelations(drawSubleading);
  SetLoadTrackSubleadingJetCorrelationsUncorrected(drawUncorrected);
  SetLoadTrackSubleadingJetCorrelationsPtWeighted(drawPtWeighted);
}

// Setter for loading inclusive jet-track correlations
void DijetHistogramManager::SetLoadTrackInclusiveJetCorrelations(const bool loadOrNot){
  fLoadJetTrackCorrelations[kTrackInclusiveJet] = loadOrNot;
}

// Setter for leading pT weighted inclusive jet-track correlations
void DijetHistogramManager::SetLoadTrackInclusiveJetCorrelationsPtWeighted(const bool loadOrNot){
  fLoadJetTrackCorrelations[kPtWeightedTrackInclusiveJet] = loadOrNot;
}

// Setter for loading all correlations related to tracks and inclusive jets
void DijetHistogramManager::SetLoadAllTrackInclusiveJetCorrelations(const bool loadInclusive, const bool loadPtWeighted){
  SetLoadTrackInclusiveJetCorrelations(loadInclusive);
  SetLoadTrackInclusiveJetCorrelationsPtWeighted(loadPtWeighted);
}

 // Setter for loading two-dimensional histograms
void DijetHistogramManager::SetLoad2DHistograms(const bool loadOrNot){
  fLoad2DHistograms = loadOrNot;
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

// Setter for JFF correction
void DijetHistogramManager::SetJffCorrection(TFile *jffFile, const bool applyCorrection){
  fApplyJffCorrection = applyCorrection;
  if(fApplyJffCorrection){
    fJffCorrectionFinder->ReadInputFile(jffFile);
  }
}

// Setter for spillover correction
void DijetHistogramManager::SetSpilloverCorrection(TFile *spilloverFile, const bool applyCorrection){
  fApplySpilloverCorrection = applyCorrection;
  if(fApplySpilloverCorrection){
    fJffCorrectionFinder->ReadSpilloverFile(spilloverFile);
  }
}

// Setter for seagull correction flag
void DijetHistogramManager::SetSeagullCorrection(const bool applyCorrection){
  fApplySeagullCorrection = applyCorrection;
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

// Getter for i:th low deltaPhi border
double DijetHistogramManager::GetDeltaPhiBorderLow(const int iDeltaPhi) const{
  return fLowDeltaPhiBinBorders[iDeltaPhi];
}

// Getter for i:th high deltaPhi border
double DijetHistogramManager::GetDeltaPhiBorderHigh(const int iDeltaPhi) const{
  return fHighDeltaPhiBinBorders[iDeltaPhi];
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
TH1D* DijetHistogramManager::GetHistogramJetTrackDeltaPhi(const int iJetTrackCorrelation, const int iCorrelationType, const int iCentrality, const int iTrackPt, const int iDeltaEta) const{
  return fhJetTrackDeltaPhi[iJetTrackCorrelation][iCorrelationType][iCentrality][iTrackPt][iDeltaEta];
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
  if(name.EqualTo("jettrackdeltaphi",TString::kIgnoreCase) || name.EqualTo("fhjettrackdeltaphi",TString::kIgnoreCase)) return GetHistogramJetTrackDeltaPhi(bin1,bin2,bin3,bin4,bin5);
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

// Get the first loaded centrality bin
int DijetHistogramManager::GetFirstCentralityBin() const{
  return fFirstLoadedCentralityBin;
}

// Get the last loaded centrality bin
int DijetHistogramManager::GetLastCentralityBin() const{
  return fLastLoadedCentralityBin;
}

// Get the first loaded track pT bin
int DijetHistogramManager::GetFirstTrackPtBin() const{
  return fFirstLoadedTrackPtBin;
}

// Get the last loaded track pT bin
int DijetHistogramManager::GetLastTrackPtBin() const{
  return fLastLoadedTrackPtBin;
}

// Getter for the number of events passing the cuts
int DijetHistogramManager::GetNEvents() const{
  return fhEvents->GetBinContent(fhEvents->FindBin(DijetHistograms::kVzCut));
}

// Getter for the number of dijets
int DijetHistogramManager::GetNDijets() const{
  return fhEvents->GetBinContent(fhEvents->FindBin(DijetHistograms::kDijet));
}

// Getter for integral over leading jet pT
double DijetHistogramManager::GetPtIntegral(int iCentrality) const{
  return fhJetPt[kLeadingJet][iCentrality]->Integral("width");
}

// Getter for integral over inclusive jet pT over 120 GeV
double DijetHistogramManager::GetInclusiveJetPtIntegral(int iCentrality) const{
  return fhJetPt[kAnyJet][iCentrality]->Integral(fhJetPt[kAnyJet][iCentrality]->FindBin(120),fhJetPt[kAnyJet][iCentrality]->GetNbinsX(),"width");
}
