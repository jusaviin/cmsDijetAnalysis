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
  fLoadJetPtClosureHistograms(false),
  fFirstLoadedCentralityBin(0),
  fLastLoadedCentralityBin(1),
  fFirstLoadedTrackPtBin(0),
  fLastLoadedTrackPtBin(1),
  fProcessAsymmetryBins(false),
  fPreprocess(false),
  fnCentralityBins(kMaxCentralityBins),
  fnTrackPtBins(kMaxTrackPtBins),
  fnAsymmetryBins(kMaxAsymmetryBins),
  fAvoidMixingPeak(false),
  fImproviseMixing(false),
  fDefaultMixingDeltaEtaFitRange(0.2)
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
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins + 1; iCentrality++){
    fCentralityBinIndices[iCentrality] = iCentrality+1;
    fCentralityBinBorders[iCentrality] = 0;
  }
  
  // Default binning for track pT
  for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins + 1; iTrackPt++){
    fTrackPtBinIndices[iTrackPt] = iTrackPt + 1;
    fFineTrackPtBinIndices[iTrackPt] = iTrackPt + 1;
    fTrackPtBinBorders[iTrackPt] = 0;
  }
  
  // Default binning for deltaPhi
  TString defaultDeltaPhiString[] = {""," Near side", " Away side", " Between peaks"};
  TString defaultCompactDeltaPhiString[] = {"", "_NearSide", "_AwaySide", "_BetweenPeaks"};
  for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
    fLowDeltaPhiBinIndices[iDeltaPhi] = iDeltaPhi+1;
    fHighDeltaPhiBinIndices[iDeltaPhi] = iDeltaPhi+2;
    fLowDeltaPhiBinBorders[iDeltaPhi] = 0;
    fHighDeltaPhiBinBorders[iDeltaPhi] = 0;
    fDeltaPhiString[iDeltaPhi] = defaultDeltaPhiString[iDeltaPhi];
    fCompactDeltaPhiString[iDeltaPhi] = defaultCompactDeltaPhiString[iDeltaPhi];
  }
  
  // Default binning for jet pT
  int defaultJetIndex[] = {25,35,45,101,0,0,0,0,0,0};
  double defaultJetPtBorder[] = {120,170,220,5020,0,0,0,0,0,0};
  for(int iJetPt = 0; iJetPt < knJetPtBins + 1; iJetPt++){
    fJetPtBinIndices[iJetPt] = defaultJetIndex[iJetPt];
    fJetPtBinBorders[iJetPt] = defaultJetPtBorder[iJetPt];
  }
  
  // Bin naming for saved asymmetry bins
  for(int iAsymmetry = 0; iAsymmetry < kMaxAsymmetryBins; iAsymmetry++){
    fAsymmetryBinName[iAsymmetry] = Form("A%d",iAsymmetry);
  }
  fAsymmetryBinName[kMaxAsymmetryBins] = "";
  
  // Initialize all the other histograms to null
  fhVertexZ = NULL;            // Vertex z position
  fhVertexZWeighted = NULL;    // Weighted vertex z-position (only meaningfull for MC)
  fhVertexZDijet = NULL;       // Vertex z position in dijet events
  fhEvents = NULL;             // Number of events surviving different event cuts
  fhTrackCuts = NULL;          // Number of tracks surviving different track cuts
  fhTrackCutsInclusive = NULL; // Number of inclusive tracks surviving different track cuts
  fhCentrality = NULL;         // Centrality of all events
  fhCentralityWeighted = NULL; // Weighted centrality distribution in all events (only meaningful for MC)
  fhCentralityDijet = NULL;    // Centrality of dijet events
  fhPtHat = NULL;              // pT hat for MC events (only meaningful for MC)
  fhPtHatWeighted = NULL;      // Weighted pT hat distribution (only meaningful for MC)
  
  // Centrality loop
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    fhDijetDphi[iCentrality] = NULL;                  // Dijet deltaPhi histograms
    fhDijetLeadingVsSubleadingPt[iCentrality] = NULL; // Leading versus subleading jet pT 2D histograms
    
    for(int iJetPt = 0; iJetPt < knJetPtBins + 1; iJetPt++){
      fhDijetAsymmetry[iCentrality][iJetPt] = NULL;   // Dijet asymmetry AJ histograms
      fhDijetXj[iCentrality][iJetPt] = NULL;          // Dijet asymmetry xJ histograms
    }
    
    // Single jet category loop
    for(int iJetCategory = 0; iJetCategory < knSingleJetCategories; iJetCategory++){
      for(int iAsymmetry = 0; iAsymmetry < kMaxAsymmetryBins+1; iAsymmetry++){
        fhJetPt[iJetCategory][iCentrality][iAsymmetry] = NULL;      // Jet pT histograms
        fhJetPhi[iJetCategory][iCentrality][iAsymmetry] = NULL;     // Jet phi histograms
        fhJetEta[iJetCategory][iCentrality][iAsymmetry] = NULL;     // Jet eta histograms
        fhJetEtaPhi[iJetCategory][iCentrality][iAsymmetry] = NULL;  // 2D eta-phi histogram for jets
      }
    } // Single jet categories loop
    
    // Event correlation type loop
    for(int iCorrelationType = 0; iCorrelationType < knCorrelationTypes; iCorrelationType++){
      
      // Loop over track categories
      for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
        fhTrackPt[iTrackType][iCorrelationType][iCentrality] = NULL;   // Track pT histograms
        
        // Loop over track pT bins
        for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins + 1; iTrackPt++){
          fhTrackPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt] = NULL;    // Track phi histograms
          fhTrackEta[iTrackType][iCorrelationType][iCentrality][iTrackPt] = NULL;    // Track eta histograms
          fhTrackEtaPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt] = NULL; // 2D eta-phi histogram for track
        } // Track pT loop
        
      } // Track category loop
      
      // Loop over jet-track correlation types
      for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
        
        // Loop over asymmetry bins
        for(int iAsymmetry = 0; iAsymmetry < kMaxAsymmetryBins+1; iAsymmetry++){
          
          // Loop over track pT bins
          for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins; iTrackPt++){
            fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentrality][iTrackPt] = NULL;      // DeltaEta and deltaPhi between jet and track
            
            // Loop over deltaEta bins
            for(int iDeltaEta = 0; iDeltaEta < knDeltaEtaBins; iDeltaEta++){
              fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentrality][iTrackPt][iDeltaEta] = NULL; // DeltaPhi between jet and track
            }
            
            // Loop over deltaPhi bins
            for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
              fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iAsymmetry][iCentrality][iTrackPt][iDeltaPhi] = NULL; // DeltaEta between jet and track
            } // DeltaPhi loop
          } // Track pT loop
        } // Asymmetry loop
      } // Jet-track correlation type loop
    } // Event correlation type loop
    
    // Jet shape and seagull correction histograms
    for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins; iTrackPt++){
        for(int iAsymmetry = 0; iAsymmetry < kMaxAsymmetryBins+1; iAsymmetry++){
          for(int iJetShape = 0; iJetShape < knJetShapeTypes; iJetShape++){
            fhJetShape[iJetShape][iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          } // Jet shape type loop
          
          // Seagull background eta histogram and fit
          fhSeagullDeltaEta[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          fSeagullFit[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = NULL;
          
        } // Asymmetry loop
      } // Track pT loop
    } // Jet-track correlation type loop
    
    // Jet pT closure histograms
    for(int iClosureType = 0; iClosureType < DijetHistograms::knClosureTypes; iClosureType++){
      for(int iGenJetPt = 0; iGenJetPt < knGenJetPtBins; iGenJetPt++){
        for(int iClosureParticle = 0; iClosureParticle < DijetHistograms::knClosureParticleTypes+1; iClosureParticle++){
          fhJetPtClosure[iClosureType][iGenJetPt][iCentrality][iClosureParticle] = NULL;
        } // Closure particle loop
      } // Gen jet pT loop
    } // Closure type loop
  } // Centrality loop
}

/*
 * Constructor
 */
DijetHistogramManager::DijetHistogramManager(TFile *inputFile) :
  DijetHistogramManager()
{
  fInputFile = inputFile;
  
  // Read card from inputfile
  fCard = new DijetCard(inputFile);
  
  // Initialize values using the information in card
  InitializeFromCard();
  
}

/*
 * Constructor
 */
DijetHistogramManager::DijetHistogramManager(TFile *inputFile, DijetCard *card) :
  DijetHistogramManager()
{
  fInputFile = inputFile;
  
  // Initialize values using the information in card
  fCard = card;
  InitializeFromCard();
  
}

/*
 * Initialize several member variables from DijetCard
 */
void DijetHistogramManager::InitializeFromCard(){
  
  // Read the collision system from the card
  TString collisionSystem = fCard->GetDataType();
  
  // Make a string for collision system based on information on the card
  fSystemAndEnergy = Form("%s 5.02 TeV",collisionSystem.Data());
  fCompactSystemAndEnergy = fSystemAndEnergy;
  fCompactSystemAndEnergy.ReplaceAll(" ","");
  fCompactSystemAndEnergy.ReplaceAll(".","v");
  
  // Read the centrality and track pT bins from the card
  fnCentralityBins = fCard->GetNCentralityBins();
  fnTrackPtBins = fCard->GetNTrackPtBins();
  
  for(int iCentrality = 0; iCentrality <= fnCentralityBins; iCentrality++){
    fCentralityBinBorders[iCentrality] = fCard->GetLowBinBorderCentrality(iCentrality);
  }
  fLastLoadedCentralityBin = fnCentralityBins-1;
  
  for(int iTrackPt = 0; iTrackPt <= fnTrackPtBins; iTrackPt++){
    fTrackPtBinBorders[iTrackPt] = fCard->GetLowBinBorderTrackPt(iTrackPt);
  }
  fLastLoadedTrackPtBin = fnTrackPtBins-1;
  
  // Remove centrality selection from pp data and local testing
  if(collisionSystem.Contains("pp") || collisionSystem.Contains("localTest")){
    fLastLoadedCentralityBin = 0;
    fCentralityBinBorders[0] = -0.5;
  }
  
  // Read the number of asymmetry bins from the card and fix the naming for the last bin
  fnAsymmetryBins = fCard->GetNAsymmetryBins();
  fAsymmetryBinName[fnAsymmetryBins] = "";
  
  // Read the deltaPhi bin borders from the card
  if(knDeltaPhiBins == fCard->GetNDeltaPhiBins()){
    for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
      fLowDeltaPhiBinBorders[iDeltaPhi] = fCard->GetLowBinBorderDeltaPhi(iDeltaPhi);
      fHighDeltaPhiBinBorders[iDeltaPhi] = fCard->GetHighBinBorderDeltaPhi(iDeltaPhi);
    }
  }
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
  fLoadJetPtClosureHistograms(in.fLoadJetPtClosureHistograms),
  fFirstLoadedCentralityBin(in.fFirstLoadedCentralityBin),
  fLastLoadedCentralityBin(in.fLastLoadedCentralityBin),
  fFirstLoadedTrackPtBin(in.fFirstLoadedTrackPtBin),
  fLastLoadedTrackPtBin(in.fLastLoadedTrackPtBin),
  fProcessAsymmetryBins(in.fProcessAsymmetryBins),
  fPreprocess(in.fPreprocess),
  fnAsymmetryBins(in.fnAsymmetryBins),
  fAvoidMixingPeak(in.fAvoidMixingPeak),
  fImproviseMixing(in.fImproviseMixing),
  fDefaultMixingDeltaEtaFitRange(in.fDefaultMixingDeltaEtaFitRange),
  fhVertexZ(in.fhVertexZ),
  fhVertexZWeighted(in.fhVertexZWeighted),
  fhVertexZDijet(in.fhVertexZDijet),
  fhEvents(in.fhEvents),
  fhTrackCuts(in.fhTrackCuts),
  fhCentrality(in.fhCentrality),
  fhCentralityWeighted(in.fhCentralityWeighted),
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
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins + 1; iCentrality++){
    fCentralityBinIndices[iCentrality] = in.fCentralityBinIndices[iCentrality];
    fCentralityBinBorders[iCentrality] = in.fCentralityBinBorders[iCentrality];
  }
  
  // Copy binning for track pT
  for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins + 1; iTrackPt++){
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
  
  // Copy binning for jet pT
  for(int iJetPt = 0; iJetPt < knJetPtBins+1; iJetPt++){
    fJetPtBinIndices[iJetPt] = in.fJetPtBinIndices[iJetPt];
    fJetPtBinBorders[iJetPt] = in.fJetPtBinBorders[iJetPt];
  }
  
  // Copy the bin naming for asymmetry bins
  for(int iAsymmetry = 0; iAsymmetry <= kMaxAsymmetryBins; iAsymmetry++){
    fAsymmetryBinName[iAsymmetry] = in.fAsymmetryBinName[iAsymmetry];
  }
  
  // Centrality loop
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    fhDijetDphi[iCentrality] = in.fhDijetDphi[iCentrality];                                   // Dijet deltaPhi histograms
    fhDijetLeadingVsSubleadingPt[iCentrality] = in.fhDijetLeadingVsSubleadingPt[iCentrality]; // Leading versus subleading jet pT 2D histograms
    
    // Jet pT loop
    for(int iJetPt = 0; iJetPt < knJetPtBins + 1; iJetPt++){
      fhDijetAsymmetry[iCentrality][iJetPt] = in.fhDijetAsymmetry[iCentrality][iJetPt];       // Dijet asymmetry AJ histograms
      fhDijetXj[iCentrality][iJetPt] = in.fhDijetXj[iCentrality][iJetPt];                     // Dijet asymmetry xJ histograms
    }
    
    // Single jet category loop
    for(int iJetCategory = 0; iJetCategory < knSingleJetCategories; iJetCategory++){
      for(int iAsymmetry = 0; iAsymmetry < kMaxAsymmetryBins+1; iAsymmetry++){
        fhJetPt[iJetCategory][iCentrality][iAsymmetry] = in.fhJetPt[iJetCategory][iCentrality][iAsymmetry];         // Jet pT histograms
        fhJetPhi[iJetCategory][iCentrality][iAsymmetry] = in.fhJetPhi[iJetCategory][iCentrality][iAsymmetry];       // Jet phi histograms
        fhJetEta[iJetCategory][iCentrality][iAsymmetry] = in.fhJetEta[iJetCategory][iCentrality][iAsymmetry];       // Jet eta histograms
        fhJetEtaPhi[iJetCategory][iCentrality][iAsymmetry] = in.fhJetEtaPhi[iJetCategory][iCentrality][iAsymmetry]; // 2D eta-phi histogram for jets
      } // Asymmetry loop
    } // Single jet categories loop
    
    // Event correlation type loop
    for(int iCorrelationType = 0; iCorrelationType < knCorrelationTypes; iCorrelationType++){
      
      // Loop over track categories
      for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
        fhTrackPt[iTrackType][iCorrelationType][iCentrality] = in.fhTrackPt[iTrackType][iCorrelationType][iCentrality];   // Track pT histograms
        
        // Loop over track pT bins
        for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins + 1; iTrackPt++){
          fhTrackPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt] = in.fhTrackPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt];    // Track phi histograms
          fhTrackEta[iTrackType][iCorrelationType][iCentrality][iTrackPt] = in.fhTrackEta[iTrackType][iCorrelationType][iCentrality][iTrackPt];    // Track eta histograms
          fhTrackEtaPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt] = in.fhTrackEtaPhi[iTrackType][iCorrelationType][iCentrality][iTrackPt]; // 2D eta-phi histogram for track
        } // Track pT loop
        
      } // Track category loop
      
      // Loop over jet-track correlation types
      for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
        
        // Loop over asymmetry bins
        for(int iAsymmetry = 0; iAsymmetry < kMaxAsymmetryBins+1; iAsymmetry++){
          
          // Loop over track pT bins
          for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins; iTrackPt++){
            fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentrality][iTrackPt] = in.fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentrality][iTrackPt]; // DeltaEta and deltaPhi between jet and track
            
            // Loop over deltaEta bins
            for(int iDeltaEta = 0; iDeltaEta < knDeltaEtaBins; iDeltaEta++){
              fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentrality][iTrackPt][iDeltaEta] = in.fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentrality][iTrackPt][iDeltaEta];         // DeltaPhi between jet and track
            }
            
            // Loop over deltaPhi bins
            for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
              fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iAsymmetry][iCentrality][iTrackPt][iDeltaPhi] = in.fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iAsymmetry][iCentrality][iTrackPt][iDeltaPhi]; // DeltaEta between jet and track
            } // DeltaPhi loop
          } // Track pT loop
        } // Asymmetry loop
      } // Jet-track correlation type loop
    } // Event correlation type loop
    
    // Jet shape and seagull correction histograms
    for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins; iTrackPt++){
        for(int iAsymmetry = 0; iAsymmetry < kMaxAsymmetryBins+1; iAsymmetry++){
          for(int iJetShape = 0; iJetShape < knJetShapeTypes; iJetShape++){
            fhJetShape[iJetShape][iJetTrack][iAsymmetry][iCentrality][iTrackPt] = in.fhJetShape[iJetShape][iJetTrack][iAsymmetry][iCentrality][iTrackPt];
          } // Jet shape type loop
          
          
          // Seagull background eta histogram and fit
          fhSeagullDeltaEta[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = in.fhSeagullDeltaEta[iJetTrack][iAsymmetry][iCentrality][iTrackPt];
          fSeagullFit[iJetTrack][iAsymmetry][iCentrality][iTrackPt] = in.fSeagullFit[iJetTrack][iAsymmetry][iCentrality][iTrackPt];
        } // Asymmetry loop
        
      } // Track pT loop
    } // Jet-track correlation type loop
    
    // Jet pT closure histograms
    for(int iClosureType = 0; iClosureType < DijetHistograms::knClosureTypes; iClosureType++){
      for(int iGenJetPt = 0; iGenJetPt < knGenJetPtBins; iGenJetPt++){
        for(int iClosureParticle = 0; iClosureParticle < DijetHistograms::knClosureParticleTypes+1; iClosureParticle++){
          fhJetPtClosure[iClosureType][iGenJetPt][iCentrality][iClosureParticle] = in.fhJetPtClosure[iClosureType][iGenJetPt][iCentrality][iClosureParticle];
        } // Closure particle loop
      } // Gen jet pT loop
    } // Closure type loop
    
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
  
  // If preprocessing flag is set, do nothing if asked to process histograms
  if(!fPreprocess){
    DoMixedEventCorrection();  // Mixed event correction needs to be done first, as we need the corrected histograms for the background subtraction
    SubtractBackgroundAndCalculateJetShape(); // Subtract the background and take projections of processed two-dimensional histograms. After that, calculate jet shape
  }
}

/*
 * Apply mixed event correction to all jet-track correlation histograms that are selected for analysis
 */
void DijetHistogramManager::DoMixedEventCorrection(){
  
  // Helper variables
  int connectedIndex;
  double scalingFactor;
  bool mixingPeakVisible;  // Some bins have a peak around xero in the mixed event distribution due to holes in acceptance
  // In these cases we should not smoothen the mixing to avoid having peak contribution in flat areas
  TH2D *correctionHistogram;
  int seagullMethod = 0;
  
  // Loop over all jet-track correlation types and apply the mixed event correction
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    if(!fLoadJetTrackCorrelations[iJetTrack]) continue; // Only correct the histograms that are selected for analysis
    for(int iAsymmetry = 0; iAsymmetry <= fnAsymmetryBins; iAsymmetry++){
      if((iJetTrack >= kTrackInclusiveJet || !fProcessAsymmetryBins) && iAsymmetry != fnAsymmetryBins) continue; // No asymmetry bins for inclusive jet-track
      
      for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
        
        // Scaling factor is number of dijets for leading and subleading and number of jets for inclusive correlations
        if(iJetTrack < kTrackInclusiveJet){
          scalingFactor = 1.0/GetPtIntegral(iCentralityBin,iAsymmetry);
        } else {
          scalingFactor = 1.0/GetInclusiveJetPtIntegral(iCentralityBin);
        }
        
        for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
          
          connectedIndex = GetConnectedIndex(iJetTrack);
          
          // The peak in the mixed event distribution can be seen the low pT bins ,especially in the central PbPb events
          // The detector efficiency improves for higher pT tracks and the peak in event mixing disappers
          mixingPeakVisible = false;
          if(((iTrackPtBin < 4 && iCentralityBin < 2) || (iTrackPtBin == 0 && iCentralityBin < 3)) && fAvoidMixingPeak) mixingPeakVisible = true;
          
          // For the peripheral bin, the maximum of the distribution due to hole in acceptance is in slightly negative deltaEta
          if(iTrackPtBin == 0 && iCentralityBin == 3 && fAvoidMixingPeak){
            fMethods->SetMixedEventFitRegion(-0.7,-0.35);
          }
          
          // Do the mixed event correction for the current jet-track correlation histogram
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kCorrected][iAsymmetry][iCentralityBin][iTrackPtBin] = fMethods->MixedEventCorrect(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kSameEvent][iAsymmetry][iCentralityBin][iTrackPtBin],fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kMixedEvent][iAsymmetry][iCentralityBin][iTrackPtBin],fhJetTrackDeltaEtaDeltaPhi[connectedIndex][kMixedEvent][iAsymmetry][iCentralityBin][iTrackPtBin],mixingPeakVisible);
          
          // Reset the fit region for event mixing to default
          if(iTrackPtBin == 0 && iCentralityBin == 3 && fAvoidMixingPeak){
            fMethods->SetMixedEventFitRegion(fDefaultMixingDeltaEtaFitRange);
          }
          
          // Remember the normalized mixed events distribution
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kMixedEventNormalized][iAsymmetry][iCentralityBin][iTrackPtBin] = fMethods->GetNormalizedMixedEvent();
          
          // Scale the histograms with the number of jets/dijets
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kCorrected][iAsymmetry][iCentralityBin][iTrackPtBin]->Scale(scalingFactor);
          
          // Apply the seagull correction after the mixed event correction
          if(fApplySeagullCorrection){
            seagullMethod = 0; // TODO: This works for MC, need to check for data after data is produced!!
            if(iJetTrack < kTrackSubleadingJet || iJetTrack > kPtWeightedTrackSubleadingJet){
              if(iCentralityBin == 0){
                if(iTrackPtBin > 0 && iTrackPtBin < 5) seagullMethod = 1;
              }
            }
            fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kCorrected][iAsymmetry][iCentralityBin][iTrackPtBin] = fMethods->DoSeagullCorrection(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kCorrected][iAsymmetry][iCentralityBin][iTrackPtBin],seagullMethod);
            
            // Get the used background eta histogram and fitted function for QA purposes
            fhSeagullDeltaEta[iJetTrack][iAsymmetry][iCentralityBin][iTrackPtBin] = (TH1D*)fMethods->GetBackgroundEta()->Clone();
            fSeagullFit[iJetTrack][iAsymmetry][iCentralityBin][iTrackPtBin] = (TF1*)fMethods->GetSeagullFit()->Clone();
          }
          
          // Apply the spillover correction to the mixed event corrected deltaEta-deltaPhi distribution
          if(fApplySpilloverCorrection && fJffCorrectionFinder->SpilloverReady()){
            if(iTrackPtBin < 5){ // Do not apply spillover correction to the highest pT bin
              if(iJetTrack < kTrackSubleadingJet || iJetTrack > kPtWeightedTrackSubleadingJet){ // Do not apply spillover correction for subleading jets
                correctionHistogram = fJffCorrectionFinder->GetDeltaEtaDeltaPhiSpilloverCorrection(iJetTrack,iCentralityBin,iTrackPtBin); // TODO: Maybe add asymmetry also to spillover, if it matters there
                fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kCorrected][iAsymmetry][iCentralityBin][iTrackPtBin]->Add(correctionHistogram,-1);
              } // If for subleading jets
            } // If for track pT
          }
          
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry loop
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
  int lowDeltaPhiBin, highDeltaPhiBin;
  bool isInclusive;
  double binError;
  double errorScale;
  TH2D *correctionHistogram;
  
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    if(!fLoadJetTrackCorrelations[iJetTrack]) continue; // Only correct the histograms that are selected for analysis
    isInclusive = (iJetTrack >= kTrackInclusiveJet); // Set the flag for inclusive jet-track correlations
    for(int iAsymmetry = 0; iAsymmetry <= fnAsymmetryBins; iAsymmetry++){
      if((iJetTrack >= kTrackInclusiveJet || !fProcessAsymmetryBins) && iAsymmetry != fnAsymmetryBins) continue; // No asymmetry bins for inclusive jet-track
      
      for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
        for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
          
          // Get the subleading/leading jet index connected to the currect leading/subleading correlation type
          connectedIndex = GetConnectedIndex(iJetTrack);
          
          // Subtract the background from the mixed event corrected histogram
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackgroundSubtracted][iAsymmetry][iCentralityBin][iTrackPtBin] = fMethods->SubtractBackground(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kCorrected][iAsymmetry][iCentralityBin][iTrackPtBin],fhJetTrackDeltaEtaDeltaPhi[connectedIndex][kCorrected][iAsymmetry][iCentralityBin][iTrackPtBin],fCard->GetMaxDeltaEta(),isInclusive);
          
          // Get also the background and background overlap region for QA purposes
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackground][iAsymmetry][iCentralityBin][iTrackPtBin] = fMethods->GetBackground();
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackgroundOverlap][iAsymmetry][iCentralityBin][iTrackPtBin] = fMethods->GetBackgroundOverlap();
          
          // Apply the JFF correction to the background subtracted deltaEta-deltaPhi distribution
          if(fApplyJffCorrection && fJffCorrectionFinder->CorrectionReady()){
            correctionHistogram = fJffCorrectionFinder->GetDeltaEtaDeltaPhiJffCorrection(iJetTrack,iCentralityBin,iTrackPtBin); // TODO: Add asymmetry to JFF correction
            fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackgroundSubtracted][iAsymmetry][iCentralityBin][iTrackPtBin]->Add(correctionHistogram,-1);
          }
          
          // Calculate the jet shape from the background subtracted histogram
          fhJetShape[kJetShape][iJetTrack][iAsymmetry][iCentralityBin][iTrackPtBin] = fMethods->GetJetShape(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kBackgroundSubtracted][iAsymmetry][iCentralityBin][iTrackPtBin]);
          
          // Get the number of two-dimensional histogram bins used for each deltaR bin in the jet shape histogram
          fhJetShape[kJetShapeBinCount][iJetTrack][iAsymmetry][iCentralityBin][iTrackPtBin] = fMethods->GetJetShapeCounts();
          
          // Get the mapping histogram of Rbins to deltaPhi-deltaEta bins
          fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kJetShapeBinMap][iAsymmetry][iCentralityBin][iTrackPtBin] = fMethods->GetJetShapeBinMap();
          
          // Project the deltaPhi and deltaEta histograms from the processed two-dimensional histograms
          for(int iCorrelationType = kMixedEvent; iCorrelationType < knCorrelationTypes; iCorrelationType++){
            
            // Only do the projections if mixed event distribution is improvised. Otherwise the projections are done already.
            if(iCorrelationType == kMixedEvent && !fImproviseMixing) continue;
            
            // DeltaPhi histogram over whole eta
            sprintf(histogramName,"%sDeltaPhiProjection%d",fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin]->GetName(),kWholeEta);
            nBins = fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin]->GetNbinsY();
            fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][kWholeEta] = (TH1D*) fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin]->ProjectionX(histogramName,1,nBins)->Clone();  // Exclude underflow and overflow bins by specifying range
            
            fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][kWholeEta]->Scale(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin]->GetYaxis()->GetBinWidth(1));  // For correct normalization, need to divide out deltaEta bin width of the two-dimensional histogram
            
            /*
             * Do error scaling and fourier fit for the background deltaPhi distribution
             *
             * Error scaling is needed because information only from 1.5 < |deltaEta| < 2.5 is used to determine the background.
             * When doing projection, root by default scaled the histogram with square root of bins projected over
             * The whole range of the analysis is |deltaEta| < 4, which means that four times the bins used to determine
             * the background are used in normalixing the error. We can fix this be scaling the error up by sqrt(4) = 2.
             * This number is provided by the DijetMethods, since if we change the limits described above, the number is
             * automatically adjusted for the new limits inside DijetMethods.
             */
            if(iCorrelationType == kBackground){
              
              for(int iBin = 1; iBin < fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][kWholeEta]->GetNbinsX(); iBin++){
                binError = fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][kWholeEta]->GetBinError(iBin);
                errorScale = fMethods->GetBackgroundErrorScalingFactor();
                fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][kWholeEta]->SetBinError(iBin,binError*errorScale);
              }
              fMethods->FourierFit(fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][kWholeEta],knFittedFlowComponents);
            }
            
            // DeltaPhi histogram over signal eta region
            sprintf(histogramName,"%sDeltaPhiProjection%d",fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin]->GetName(),kSignalEtaRegion);
            fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][kSignalEtaRegion] = (TH1D*)fMethods->ProjectSignalDeltaPhi(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin])->Clone();
            
            // DeltaPhi histogram over background eta region
            sprintf(histogramName,"%sDeltaPhiProjection%d",fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin]->GetName(),kBackgroundEtaRegion);
            fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][kBackgroundEtaRegion] = (TH1D*)fMethods->ProjectBackgroundDeltaPhi(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin])->Clone();
              
            for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
              
              lowDeltaPhiBin = fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin]->GetXaxis()->FindBin(fLowDeltaPhiBinBorders[iDeltaPhi]);
              highDeltaPhiBin = fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin]->GetXaxis()->FindBin(fHighDeltaPhiBinBorders[iDeltaPhi]);
              sprintf(histogramName,"%sDeltaEtaProjection%d",fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin]->GetName(),iDeltaPhi);
              fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][iDeltaPhi] = (TH1D*) fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin]->ProjectionY(histogramName,lowDeltaPhiBin,highDeltaPhiBin)->Clone();
              
              // To retain the normalization, we must scale the histograms with the number of bins projected over and by the width of deltaPhi bin
              nProjectedBins = highDeltaPhiBin - lowDeltaPhiBin + 1;
              if(iDeltaPhi > kWholePhi) fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][iDeltaPhi]->Scale(1.0/nProjectedBins);
              fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][iDeltaPhi]->Scale(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin]->GetXaxis()->GetBinWidth(1));
              
            } // DeltaPhi loop
          } // Correlation type loop
        } // Track pT loop
      } // Centrality loop
    } // Asymmetry bin loop
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
      
      for(int iAsymmetry = 0; iAsymmetry <= fnAsymmetryBins; iAsymmetry++){
        if((iJetTrack >= kTrackInclusiveJet || !fProcessAsymmetryBins) && iAsymmetry != fnAsymmetryBins) continue; // No asymmetry bins for inclusive jet-track
        for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
          correctionHistogram = jffCorrectionFinder->GetJetShapeJffCorrection(iJetTrack,iCentralityBin,iTrackPtBin);
          fhJetShape[kJetShape][iJetTrack][iAsymmetry][iCentralityBin][iTrackPtBin]->Scale(scalingFactor);  // Need to scale with the number of dijets/all jets since the correction is also normalized to the number of dijets/all jets
          fhJetShape[kJetShape][iJetTrack][iAsymmetry][iCentralityBin][iTrackPtBin]->Add(correctionHistogram,-1);
        } // Track pT loop
      } // Asymmetry loop
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
    for(int iAsymmetry = 0; iAsymmetry < fnAsymmetryBins; iAsymmetry++){
      if((iJetTrack >= kTrackInclusiveJet || !fProcessAsymmetryBins) && iAsymmetry != fnAsymmetryBins) continue; // No asymmetry bins for inclusive jet-track
      
      for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
        jetShapeSum = (TH1D*)fhJetShape[kJetShape][iJetTrack][iAsymmetry][iCentralityBin][fFirstLoadedTrackPtBin]->Clone(Form("jetShapeSum%d%d%d",iJetTrack,iAsymmetry,iCentralityBin));
        
        // First, sum all pT bins together
        for(int iTrackPtBin = fFirstLoadedTrackPtBin+1; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
          jetShapeSum->Add(fhJetShape[kJetShape][iJetTrack][iAsymmetry][iCentralityBin][iTrackPtBin]);
        } // Track pT loop
        
        // Then calculate the integral for deltaR < 1
        jetShapeIntegral = jetShapeSum->Integral(1,jetShapeSum->FindBin(0.99),"width");
        
        // Finally, normalize each pT bin with the intagral
        for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
          fhJetShape[kJetShape][iJetTrack][iAsymmetry][iCentralityBin][iTrackPtBin]->Scale(1.0/jetShapeIntegral);
        } // Track pT loop
        
      } // Centrality loop
    } // Asymmetry loop
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
    fhVertexZDijet = (TH1D*) fInputFile->Get("vertexZdijet");              // Vertex z position in dijet events
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
  
  // Load jet pT closure histograms
  LoadJetPtClosureHistograms();
  
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
  
  // Define arrays to help find the histograms
  int axisIndices[2] = {0};
  int lowLimits[2] = {0};
  int highLimits[2] = {0};
  
  int centralityIndex[] = {4,4,3,3}; // For jet histograms without dijet requirement, there are only three axes.
  int nAxesArray[] = {2,2,1,1};      // Number of constraining axes. No asymmetry for histograms without dijet requirement
  int nAxes = 2;                     // Number of constraining axes for this iteration
  
  // Bin indices for pT cut 120 histograms
  int lowPtBinIndex[1];
  int highPtBinIndex[1];
  double lowPtBinLimit[1] = {120};
  double highPtBinLimit[1] = {5020};
  
  for(int iJetCategory = 0; iJetCategory < knSingleJetCategories; iJetCategory++){
    
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      
      // Select the bin indices
      if(iCentralityBin == fnCentralityBins-1) {
        duplicateRemoverCentrality = 0;
      } else {
        duplicateRemoverCentrality = -1;
      }
      lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
      higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
      
      for(int iAsymmetry = 0; iAsymmetry < fnAsymmetryBins+1; iAsymmetry++){
        
        if(nAxesArray[iJetCategory] == 1 && (iAsymmetry > 0 && iAsymmetry < fnAsymmetryBins)) continue; // Only fill histograms without asymmetry selection for inclusive jet histograms. Fill pT > 120 for asymmetry bin 0.
        
        axisIndices[0] = centralityIndex[iJetCategory]; axisIndices[1] = 3;
        lowLimits[0] = lowerCentralityBin; lowLimits[1] = iAsymmetry+1;
        highLimits[0] = higherCentralityBin; highLimits[1] = iAsymmetry+1;
        
        // For the last index of array, disable asymmetry contraint
        nAxes = nAxesArray[iJetCategory];
        if(iAsymmetry == fnAsymmetryBins) nAxes = 1;
        
        // Fill histograms without dijet asymmetry information with pT cut 120 to asymmetry bin 0
        if(nAxes == 1 && iAsymmetry == 0){
          nAxes = 2;           // Set the number of constraining axes to 2
          axisIndices[1] = 0;  // Set the second constraining axis to be the pT axis
          SetBinIndices("anyJet",1,lowPtBinIndex,highPtBinIndex,lowPtBinLimit,highPtBinLimit,0);  // Find the bin indices for pT
          lowLimits[1] = lowPtBinIndex[0];   // Set the obtained index for low pT bin
          highLimits[1] = highPtBinIndex[0]; // Set the obtained index for high pT bin
        }
        
        // Always load single jet pT histograms
        fhJetPt[iJetCategory][iCentralityBin][iAsymmetry] = FindHistogram(fInputFile,fSingleJetHistogramName[iJetCategory],0,nAxes,axisIndices,lowLimits,highLimits);
        
        if(!fLoadSingleJets[iJetCategory]) continue;  // Only load the remaining single jet histograms is selected
        
        fhJetPhi[iJetCategory][iCentralityBin][iAsymmetry] = FindHistogram(fInputFile,fSingleJetHistogramName[iJetCategory],1,nAxes,axisIndices,lowLimits,highLimits);
        fhJetEta[iJetCategory][iCentralityBin][iAsymmetry] = FindHistogram(fInputFile,fSingleJetHistogramName[iJetCategory],2,nAxes,axisIndices,lowLimits,highLimits);
        if(fLoad2DHistograms) fhJetEtaPhi[iJetCategory][iCentralityBin][iAsymmetry] = FindHistogram2D(fInputFile,fSingleJetHistogramName[iJetCategory],1,2,nAxes,axisIndices,lowLimits,highLimits);
      } // Loop over dijet asymmetry bins
    } // Loop over centrality bins
  } // Loop over single jet categories
}

/*
 * Loader for dijet histograms
 *
 * THnSparse for dijets:
 *
 *   Histogram name        Axis index       Content of axis
 * ----------------------------------------------------------
 *        dijet              Axis 0         Leading jet pT
 *        dijet              Axis 1        Subleading jet pT
 *        dijet              Axis 2         Dijet deltaPhi
 *        dijet              Axis 3        Dijet asymmetry AJ
 *        dijet              Axis 4           Centrality
 *        dijet              Axis 5        Dijet asymmetry xJ
 */
void DijetHistogramManager::LoadDijetHistograms(){
  
  if(!fLoadDijetHistograms) return; // Do not load the histograms if they are not selected for loading
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int duplicateRemoverJetPt = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int lowerJetPtBin = 0;
  int higherJetPtBin = 0;
  int axisIndices[2] = {0};
  int lowLimits[2] = {0};
  int highLimits[2] = {0};
  
  for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
    
    // Select the centrality bin indices
    if(iCentralityBin == fnCentralityBins-1) duplicateRemoverCentrality = 0;
    lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
    higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
    
    fhDijetDphi[iCentralityBin] = FindHistogram(fInputFile,"dijet",2,4,lowerCentralityBin,higherCentralityBin);
    fhDijetAsymmetry[iCentralityBin][knJetPtBins] = FindHistogram(fInputFile,"dijet",3,4,lowerCentralityBin,higherCentralityBin);
    fhDijetXj[iCentralityBin][knJetPtBins] = FindHistogram(fInputFile,"dijet",5,4,lowerCentralityBin,higherCentralityBin);
    if(fLoad2DHistograms) fhDijetLeadingVsSubleadingPt[iCentralityBin] = FindHistogram2D(fInputFile,"dijet",0,1,4,lowerCentralityBin,higherCentralityBin);

    // Load the asymmetry histograms in jet pT bins
    for(int iJetPt = 0; iJetPt < knJetPtBins; iJetPt++){

      // Select the jet pT bin indices
      duplicateRemoverJetPt = (iJetPt == knJetPtBins - 1) ? 0 : -1;  // Set duplicate remover to 0 for the last jet pT bin
      lowerJetPtBin = fJetPtBinIndices[iJetPt];
      higherJetPtBin = fJetPtBinIndices[iJetPt+1]+duplicateRemoverJetPt;
      
      // Setup axes with restrictions (0 = leading jet pT, 4 = centrality)
      axisIndices[0] = 0; axisIndices[1] = 4;
      lowLimits[0] = lowerJetPtBin; lowLimits[1] = lowerCentralityBin;
      highLimits[0] = higherJetPtBin; highLimits[1] = higherCentralityBin;
      
      // Find the histograms with restrictions
      fhDijetAsymmetry[iCentralityBin][iJetPt] = FindHistogram(fInputFile,"dijet",3,2,axisIndices,lowLimits,highLimits);
      fhDijetXj[iCentralityBin][iJetPt] = FindHistogram(fInputFile,"dijet",5,2,axisIndices,lowLimits,highLimits);
    }
    
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
        if(iCentralityBin == fnCentralityBins-1) {
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
        fhTrackPhi[iTrackType][iCorrelationType][iCentralityBin][fnTrackPtBins] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],1,2+axisAdder,axisIndices,lowLimits,highLimits);
        fhTrackEta[iTrackType][iCorrelationType][iCentralityBin][fnTrackPtBins] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],2,2+axisAdder,axisIndices,lowLimits,highLimits);
        if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCorrelationType][iCentralityBin][fnTrackPtBins] = FindHistogram2D(fInputFile,fTrackHistogramNames[iTrackType],1,2,2+axisAdder,axisIndices,lowLimits,highLimits);
        
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
 * Loader for jet-track correlation histograms
 *
 * THnSparses for jet-track correlations:
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
  int axisIndices[5] = {0};
  int lowLimits[5] = {0};
  int highLimits[5] = {0};
  int nRestrictionAxes;
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int duplicateRemoverTrackPt = -1;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  
  // Load all the histograms from the file
  for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
    if(!fLoadJetTrackCorrelations[iJetTrack]) continue; // Only load categories of correlation that are selected
    for(int iAsymmetry = 0; iAsymmetry <= fnAsymmetryBins; iAsymmetry++){
      if((iJetTrack >= kTrackInclusiveJet || !fProcessAsymmetryBins) && iAsymmetry != fnAsymmetryBins) continue; // No asymmetry bins for inclusive jet-track
      for(int iCorrelationType = 0; iCorrelationType <= kMixedEvent; iCorrelationType++){ // Data file contains only same and mixed event distributions
        for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
          
          // Select the bin indices
          if(iCentralityBin == fnCentralityBins-1) {
            duplicateRemoverCentrality = 0;
          } else {
            duplicateRemoverCentrality = -1;
          }
          lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
          higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
          
          for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
            
            // If specified, improvise mixed event distributions from deltaPhi side band region instead of reading them from file
            // Mixing can only be improvised if the two-dimensional same event distribution is previously loaded
            if(fImproviseMixing && iCorrelationType == kMixedEvent && fLoad2DHistograms){
              fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kMixedEvent][iAsymmetry][iCentralityBin][iTrackPtBin] = fMethods->ImproviseMixedEvent(fhJetTrackDeltaEtaDeltaPhi[iJetTrack][kSameEvent][iAsymmetry][iCentralityBin][iTrackPtBin]);
              continue;
            }
            
            // Select the bin indices for track pT
            lowerTrackPtBin = fTrackPtBinIndices[iTrackPtBin];
            higherTrackPtBin = fTrackPtBinIndices[iTrackPtBin+1]+duplicateRemoverTrackPt;
            
            // Setup the axes with restrictions, that are common for all jet-track correlation histograms
            axisIndices[0] = 5; lowLimits[0] = iCorrelationType+1; highLimits[0] = iCorrelationType+1;   // Same/mixed event
            axisIndices[1] = 4; lowLimits[1] = lowerCentralityBin; highLimits[1] = higherCentralityBin;  // Centrality
            axisIndices[2] = 0; lowLimits[2] = lowerTrackPtBin;    highLimits[2] = higherTrackPtBin;     // Track pT
            axisIndices[3] = 3; lowLimits[3] = iAsymmetry+1;       highLimits[3] = iAsymmetry+1;         // Asymmetry
            
            nRestrictionAxes = 4;
            
            // Remove the asymmetry restriction for the last iAsymmetry bin. That will be asymmetry inclusive.
            if(iAsymmetry == fnAsymmetryBins) nRestrictionAxes = 3;
            
            // Indexing is different for inclusive jet-track correlation, as they do not have dijet asymmetry
            if(iJetTrack >= kTrackInclusiveJet){
              axisIndices[0] = 4;
              axisIndices[1] = 3;
              nRestrictionAxes = 3;
            }
            
            fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][kWholeEta] = FindHistogram(fInputFile,fJetTrackHistogramNames[iJetTrack],1,nRestrictionAxes,axisIndices,lowLimits,highLimits);
            if(fLoad2DHistograms) fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin] = FindHistogram2D(fInputFile,fJetTrackHistogramNames[iJetTrack],1,2,nRestrictionAxes,axisIndices,lowLimits,highLimits);
            
            // DeltaPhi binning for deltaEta histogram
            for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
              axisIndices[3] = 1; lowLimits[3] = fLowDeltaPhiBinIndices[iDeltaPhi]; highLimits[3] = fHighDeltaPhiBinIndices[iDeltaPhi]; // DeltaPhi
              axisIndices[4] = 3; lowLimits[4] = iAsymmetry+1;                      highLimits[4] = iAsymmetry+1;                       // Asymmetry
              nRestrictionAxes = 5;
              
              // No asymmetry restriction for inclusive jets or the last asymmetry bin
              if(iJetTrack >= kTrackInclusiveJet || iAsymmetry == fnAsymmetryBins) nRestrictionAxes = 4;
              
              fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][iDeltaPhi] = FindHistogram(fInputFile,fJetTrackHistogramNames[iJetTrack],2,nRestrictionAxes,axisIndices,lowLimits,highLimits);
            } // DeltaPhi loop
          } // Track pT loop
        } // Centrality loop
      } // Asymmetry loop
    } // Correlation type loop
  } // Jet-track correlation category loop
}

/*
 * Loader for jet pT closure histograms
 *
 * THnSparse for closure histograms:
 *
 *   Histogram name: jetPtClosure
 *
 *     Axis index                  Content of axis
 * -----------------------------------------------------------
 *       Axis 0             Leading / subleading / inclusive
 *       Axis 1              Matched generator level jet pT
 *       Axis 2                       Centrality
 *       Axis 3                      Quark / gluon
 *       Axis 4             Matched reco to gen jet pT ratio
 */
void DijetHistogramManager::LoadJetPtClosureHistograms(){
  
  if(!fLoadJetPtClosureHistograms) return; // Do not load the histograms if they are not selected for loading
  
  // Define arrays to help find the histograms
  int axisIndices[4] = {0};
  int lowLimits[4] = {0};
  int highLimits[4] = {0};
  int nRestrictionAxes = 4;
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  
  // Load all the histograms from the file
  for(int iClosureType = 0; iClosureType < DijetHistograms::knClosureTypes; iClosureType++){
    for(int iGenJetPt = 0; iGenJetPt < knGenJetPtBins; iGenJetPt++){
      for(int iClosureParticle = 0; iClosureParticle < DijetHistograms::knClosureParticleTypes+1; iClosureParticle++){
        for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
          
          // Select the bin indices
          if(iCentralityBin == fnCentralityBins-1) {
            duplicateRemoverCentrality = 0;
          } else {
            duplicateRemoverCentrality = -1;
          }
          lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
          higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
          
          // Setup the axes with restrictions
          nRestrictionAxes = 4;
          axisIndices[0] = 0; lowLimits[0] = iClosureType+1; highLimits[0] = iClosureType+1;   // Leading/subleading/inclusive
          axisIndices[1] = 1; lowLimits[1] = iGenJetPt+1;    highLimits[1] = iGenJetPt+1;      // Gen jet pT
          axisIndices[2] = 2; lowLimits[2] = lowerCentralityBin; highLimits[2] = higherCentralityBin; // Centrality
          axisIndices[3] = 3; lowLimits[3] = iClosureParticle+1; highLimits[3] = iClosureParticle+1;  // Qoark/gluon
          
          // For the last closure particle bin no restrictions for quark/gluon jets
          if(iClosureParticle == DijetHistograms::knClosureParticleTypes) nRestrictionAxes = 3;
          
        fhJetPtClosure[iClosureType][iGenJetPt][iCentralityBin][iClosureParticle] = FindHistogram(fInputFile,"jetPtClosure",4,nRestrictionAxes,axisIndices,lowLimits,highLimits);
          
        } // Centrality loop
      } // Closure particle loop
    } // Gen jet pT loop
  } // Closure type loop
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
  
  // If cannot find histogram, inform that it could not be found and return null
  if(histogramArray == nullptr){
    cout << "Could not find " << name << ". Skipping loading this histogram." << endl;
    return NULL;
  }
  
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
 * Extract a histogram with given restrictions on other axes in THnSparse
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
  
  // If cannot find histogram, inform that it could not be found and return null
  if(histogramArray == nullptr){
    cout << "Could not find " << name << ". Skipping loading this histogram." << endl;
    return NULL;
  }
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for each histogram that is read from the file
  char newName[200];
  sprintf(newName,"%s",histogramArray->GetName());
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    sprintf(newName,"%s_%d=%d-%d",newName,axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]);
  }

  // Project out the histogram and give it the created unique name
  TH1D *projectedHistogram = NULL;
  
  // Check that we are not trying to project a non-existing axis
  if(xAxis < histogramArray->GetNdimensions()){
    projectedHistogram =(TH1D*) histogramArray->Projection(xAxis);
    projectedHistogram->SetName(newName);
  
    // Apply bin width normalization to the projected histogram
    projectedHistogram->Scale(1.0,"width");
  }
  
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
    if(fhVertexZDijet) fhVertexZDijet->Write();  // Vertex z position in dijet events
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
    if(!fLoadSingleJets[iJetCategory]) continue;  // Only write the selected histograms
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory(fSingleJetHistogramName[iJetCategory])) gDirectory->mkdir(fSingleJetHistogramName[iJetCategory]);
    gDirectory->cd(fSingleJetHistogramName[iJetCategory]);
    
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      
      // Check that the histograms are actually there before trying to save them.
      if(fhJetPt[iJetCategory][iCentralityBin][fnAsymmetryBins] == NULL) {
        cout << "Could not find histograms of type " << fSingleJetHistogramName[iJetCategory] << " to write. Will skip writing these." << endl;
        continue;
      }
      
      for(int iAsymmetry = 0; iAsymmetry <= fnAsymmetryBins; iAsymmetry++){
        
        // Single jet pT
        sprintf(histogramNamer,"%sPt_C%d%s",fSingleJetHistogramName[iJetCategory],iCentralityBin,fAsymmetryBinName[iAsymmetry].Data());
        if(fhJetPt[iJetCategory][iCentralityBin][iAsymmetry]) fhJetPt[iJetCategory][iCentralityBin][iAsymmetry]->Write(histogramNamer);
        
        // Single jet phi
        sprintf(histogramNamer,"%sPhi_C%d%s",fSingleJetHistogramName[iJetCategory],iCentralityBin,fAsymmetryBinName[iAsymmetry].Data());
        if(fhJetPhi[iJetCategory][iCentralityBin][iAsymmetry]) fhJetPhi[iJetCategory][iCentralityBin][iAsymmetry]->Write(histogramNamer);
        
        // Single jet eta
        sprintf(histogramNamer,"%sEta_C%d%s",fSingleJetHistogramName[iJetCategory],iCentralityBin,fAsymmetryBinName[iAsymmetry].Data());
        if(fhJetEta[iJetCategory][iCentralityBin][iAsymmetry]) fhJetEta[iJetCategory][iCentralityBin][iAsymmetry]->Write(histogramNamer);
        
        //Single jet eta-phi
        sprintf(histogramNamer,"%sEtaPhi_C%d%s",fSingleJetHistogramName[iJetCategory],iCentralityBin,fAsymmetryBinName[iAsymmetry].Data());
        if(fLoad2DHistograms && fhJetEtaPhi[iJetCategory][iCentralityBin][iAsymmetry]) fhJetEtaPhi[iJetCategory][iCentralityBin][iAsymmetry]->Write(histogramNamer);
      } // Loop over asymmetry bins
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
      fhDijetAsymmetry[iCentralityBin][knJetPtBins]->Write(histogramNamer);
      
      // xJ histograms are newer addition. Do not crash the code if they are not there for older files
      if(fhDijetXj[iCentralityBin][knJetPtBins] != NULL){
        sprintf(histogramNamer,"dijetXj_C%d",iCentralityBin);
        fhDijetXj[iCentralityBin][knJetPtBins]->Write(histogramNamer);
      }
      
      for(int iJetPt = 0; iJetPt < knJetPtBins; iJetPt++){
        
        // Dijet asymmetry AJ in jet pT bins
        sprintf(histogramNamer,"dijetAsymmetry_C%dT%d",iCentralityBin,iJetPt);
        fhDijetAsymmetry[iCentralityBin][iJetPt]->Write(histogramNamer);
        
        // Dijet asymmetry xJ in jet pT bins
        if(fhDijetXj[iCentralityBin][iJetPt] != NULL){
          sprintf(histogramNamer,"dijetXj_C%dT%d",iCentralityBin,iJetPt);
          fhDijetXj[iCentralityBin][iJetPt]->Write(histogramNamer);
        }
      }
      
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
        sprintf(histogramNamer,"%sPhi%s_C%dT%d",fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,fnTrackPtBins);
        fhTrackPhi[iTrackType][iCorrelationType][iCentralityBin][fnTrackPtBins]->Write(histogramNamer);
        
        // pT integrated track eta
        sprintf(histogramNamer,"%sEta%s_C%dT%d",fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,fnTrackPtBins);
        fhTrackEta[iTrackType][iCorrelationType][iCentralityBin][fnTrackPtBins]->Write(histogramNamer);
        
        // pT integrated track eta-phi
        sprintf(histogramNamer,"%sEtaPhi%s_C%dT%d",fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,fnTrackPtBins);
        if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCorrelationType][iCentralityBin][fnTrackPtBins]->Write(histogramNamer);
        
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
      
      // Only same and mixed event histograms are filled for preprocessing. Skip the rest.
      if(iCorrelationType > kMixedEvent && fPreprocess) continue;
      
      // Loop over asymmetry bins
      for(int iAsymmetry = 0; iAsymmetry <= fnAsymmetryBins; iAsymmetry++){
        if((iJetTrack >= kTrackInclusiveJet || !fProcessAsymmetryBins) && iAsymmetry != fnAsymmetryBins) continue; // No asymmetry bins for inclusive jet-track
        
        // Loop over centrality bins
        for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
          
          // Loop over track pT bins
          for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
            
            // Jet-track deltaPhi
            for(int iDeltaEta = 0; iDeltaEta < knDeltaEtaBins; iDeltaEta++){
              
              if(iDeltaEta > kWholeEta && iCorrelationType < kMixedEventNormalized) continue; // DeltaEta slicing not implemented for same and mixed event
              sprintf(histogramNamer,"%sDeltaPhi%s%s_%sC%dT%d",fJetTrackHistogramNames[iJetTrack], fCompactCorrelationTypeString[iCorrelationType].Data(),fCompactDeltaEtaString[iDeltaEta], fAsymmetryBinName[iAsymmetry].Data(),iCentralityBin,iTrackPtBin);
              if(fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][iDeltaEta] != NULL) fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][iDeltaEta]->Write(histogramNamer);
            }
            
            // Jet-track deltaEtaDeltaPhi
            sprintf(histogramNamer,"%sDeltaEtaDeltaPhi%s_%sC%dT%d",fJetTrackHistogramNames[iJetTrack], fCompactCorrelationTypeString[iCorrelationType].Data(),fAsymmetryBinName[iAsymmetry].Data(),iCentralityBin,iTrackPtBin);
            if(fLoad2DHistograms && fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin]) fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin]->Write(histogramNamer);
            
            // DeltaPhi binning for deltaEta histogram
            for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
              sprintf(histogramNamer,"%sDeltaEta%s%s_%sC%dT%d",fJetTrackHistogramNames[iJetTrack], fCompactCorrelationTypeString[iCorrelationType].Data(),fCompactDeltaPhiString[iDeltaPhi].Data(), fAsymmetryBinName[iAsymmetry].Data(),iCentralityBin,iTrackPtBin);
              if(fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][iDeltaPhi]) fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][iDeltaPhi]->Write(histogramNamer);
            } // DeltaPhi loop
          } // Track pT loop
        } // Centrality loop
      } // Asymmetry loop
    } // Correlation type loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Jet-track correlation category loop
  
  // Jet shape histograms are not producess in preprocessing
  if(!fPreprocess){
    
    // Write the jet shape histograms to the output file
    for(int iJetShape = 0; iJetShape < knJetShapeTypes; iJetShape++){
      
      // Loop over jet-track correlation categories
      for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
        if(!fLoadJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
        
        // Create a directory for the histograms if it does not already exist
        sprintf(histogramNamer,"%s_%s",fJetShapeHistogramName[iJetShape],fJetTrackHistogramNames[iJetTrack]);
        if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
        gDirectory->cd(histogramNamer);
        
        // Loop over asymmetry
        for(int iAsymmetry = 0; iAsymmetry <= fnAsymmetryBins; iAsymmetry++){
          if((iJetTrack >= kTrackInclusiveJet || !fProcessAsymmetryBins) && iAsymmetry != fnAsymmetryBins) continue; // No asymmetry bins for inclusive jet-track
          
          // Loop over centrality
          for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
            
            // Loop over track pT bins
            for(int iTrackPt = fFirstLoadedTrackPtBin; iTrackPt <= fLastLoadedTrackPtBin; iTrackPt++){
              sprintf(histogramNamer,"%s_%s_%sC%dT%d",fJetShapeHistogramName[iJetShape], fJetTrackHistogramNames[iJetTrack],fAsymmetryBinName[iAsymmetry].Data(),iCentrality,iTrackPt);
              if(fhJetShape[iJetShape][iJetTrack][iAsymmetry][iCentrality][iTrackPt]) fhJetShape[iJetShape][iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Write(histogramNamer);
              
            } // Track pT loop
          } // Centrality loop
        } // Asymmetry loop
        
        // Return back to main directory
        gDirectory->cd("../");
        
      } // Jet-track correlation category loop
      
    } // Jet shape type loop
    
  } // if for preprocess
  
  // Write the jet pT closure histograms to a file
  if(fLoadJetPtClosureHistograms){
    
    // Loop over closure types (leading/subleading/inclusive)
    for(int iClosureType = 0; iClosureType < DijetHistograms::knClosureTypes; iClosureType++){
      
      // Create a directory for the histograms if it does not already exist
      sprintf(histogramNamer,"jetPtClosure_%s",fSingleJetHistogramName[iClosureType]);
      if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
      gDirectory->cd(histogramNamer);
      
      // Centrality loop
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        
        // Loop over closure particles (quark/gluon/no selection)
        for(int iClosureParticle = 0; iClosureParticle < DijetHistograms::knClosureParticleTypes+1; iClosureParticle++){
          
          // Loop over generator level jet pT bins
          for(int iGenJetPt = 0; iGenJetPt < knGenJetPtBins; iGenJetPt++){
            sprintf(histogramNamer,"jetPtClosure_%s%s_C%dT%d",fSingleJetHistogramName[iClosureType],fClosureParticleName[iClosureParticle],iCentrality,iGenJetPt);
            fhJetPtClosure[iClosureType][iGenJetPt][iCentrality][iClosureParticle]->Write(histogramNamer);
            
          } // Generator level jet pT loop
        } // Closure particle type (quark/gluon) loop
      } // Centrality loop
      
      // Return back to main directory
      gDirectory->cd("../");
    } // Closure type (leading/subleading/inclusive) loop
    
  } // Writing jet pT closure histograms
  
  // Write the card to the output file if it is not already written
  if(!gDirectory->GetDirectory("JCard")) fCard->Write(outputFile);
  
  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;
  
  // Write seagull QA histograms only if the flag for seagull correction is set and we are not only doing preprocessing
  if(fApplySeagullCorrection && !fPreprocess){
    
    // Write the QA histograms to a separate QA file
    TString qaFileName = fileName;
    qaFileName.ReplaceAll(".root","_QA.root");
    TFile *qaFile = new TFile(qaFileName,fileOption);
    
    // Loop over jet-track correlation categories
    for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
      if(!fLoadJetTrackCorrelations[iJetTrack]) continue;  // Only draw the selected categories
      
      // Create a directory for the histograms if it does not already exist
      sprintf(histogramNamer,"seagullDeltaEta_%s",fJetTrackHistogramNames[iJetTrack]);
      if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
      gDirectory->cd(histogramNamer);
      
      // Loop over asymmetry
      for(int iAsymmetry = 0; iAsymmetry <= fnAsymmetryBins; iAsymmetry++){
        if((iJetTrack >= kTrackInclusiveJet || !fProcessAsymmetryBins) && iAsymmetry != fnAsymmetryBins) continue; // No asymmetry bins for inclusive jet-track
        
        // Loop over centrality
        for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
          
          // Loop over track pT bins
          for(int iTrackPt = fFirstLoadedTrackPtBin; iTrackPt <= fLastLoadedTrackPtBin; iTrackPt++){
            sprintf(histogramNamer,"seagullDeltaEta_%s_%sC%dT%d",fJetTrackHistogramNames[iJetTrack], fAsymmetryBinName[iAsymmetry].Data(),iCentrality,iTrackPt);
            if(fhSeagullDeltaEta[iJetTrack][iAsymmetry][iCentrality][iTrackPt]) fhSeagullDeltaEta[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Write(histogramNamer);
            
          } // Track pT loop
        } // Centrality loop
      } // Asymmetry loop
      
      // Return back to main directory
      gDirectory->cd("../");
      
      // Create a directory for the histograms if it does not already exist
      sprintf(histogramNamer,"seagullFit_%s",fJetTrackHistogramNames[iJetTrack]);
      if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
      gDirectory->cd(histogramNamer);
      
      // Loop over asymmetry
      for(int iAsymmetry = 0; iAsymmetry <= fnAsymmetryBins; iAsymmetry++){
        if((iJetTrack >= kTrackInclusiveJet || !fProcessAsymmetryBins) && iAsymmetry != fnAsymmetryBins) continue; // No asymmetry bins for inclusive jet-track
        
        // Loop over centrality
        for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
          
          // Loop over track pT bins
          for(int iTrackPt = fFirstLoadedTrackPtBin; iTrackPt <= fLastLoadedTrackPtBin; iTrackPt++){
            sprintf(histogramNamer,"seagullFit_%s_%sC%dT%d",fJetTrackHistogramNames[iJetTrack],fAsymmetryBinName[iAsymmetry].Data(),iCentrality,iTrackPt);
            if(fSeagullFit[iJetTrack][iAsymmetry][iCentrality][iTrackPt]) fSeagullFit[iJetTrack][iAsymmetry][iCentrality][iTrackPt]->Write(histogramNamer);
            
          } // Track pT loop
        } // Centrality loop
      } // Asymmetry loop
      
      // Return back to main directory
      gDirectory->cd("../");
      
    } // Jet-track correlation category loop
    
    
    // Close the QA file after everything is written
    qaFile->Close();
    
    // Delete the qaFile object
    delete qaFile;
  } // Seagull correction
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
    fhVertexZDijet = (TH1D*) fInputFile->Get("vertexZdijet");              // Vertex z position in dijet events
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
      
      for(int iAsymmetry = 0; iAsymmetry <= fnAsymmetryBins; iAsymmetry++){
        
        // Always load single jet pT histograms
        sprintf(histogramNamer,"%s/%sPt_C%d%s",fSingleJetHistogramName[iJetCategory],fSingleJetHistogramName[iJetCategory], iCentralityBin,fAsymmetryBinName[iAsymmetry].Data());
        fhJetPt[iJetCategory][iCentralityBin][iAsymmetry] = (TH1D*) fInputFile->Get(histogramNamer);
        
        if(!fLoadSingleJets[iJetCategory]) continue;  // Only load the loaded the selected histograms
        
        // Single jet phi
        sprintf(histogramNamer,"%s/%sPhi_C%d%s",fSingleJetHistogramName[iJetCategory],fSingleJetHistogramName[iJetCategory], iCentralityBin,fAsymmetryBinName[iAsymmetry].Data());
        fhJetPhi[iJetCategory][iCentralityBin][iAsymmetry] = (TH1D*) fInputFile->Get(histogramNamer);
        
        // Single jet eta
        sprintf(histogramNamer,"%s/%sEta_C%d%s",fSingleJetHistogramName[iJetCategory],fSingleJetHistogramName[iJetCategory], iCentralityBin,fAsymmetryBinName[iAsymmetry].Data());
        fhJetEta[iJetCategory][iCentralityBin][iAsymmetry] = (TH1D*) fInputFile->Get(histogramNamer);
        
        // Single jet eta-phi
        sprintf(histogramNamer,"%s/%sEtaPhi_C%d%s",fSingleJetHistogramName[iJetCategory],fSingleJetHistogramName[iJetCategory], iCentralityBin,fAsymmetryBinName[iAsymmetry].Data());
        if(fLoad2DHistograms) fhJetEtaPhi[iJetCategory][iCentralityBin][iAsymmetry] = (TH2D*) fInputFile->Get(histogramNamer);
      }
    } // Loop over centrality bins
    
  } // Loop over single jet categories
  
  // Load the dijet histograms from the input file
  if(fLoadDijetHistograms){
    
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      
      // Dijet deltaPhi
      sprintf(histogramNamer,"dijet/dijetDeltaPhi_C%d",iCentralityBin);
      fhDijetDphi[iCentralityBin] = (TH1D*) fInputFile->Get(histogramNamer);
      
      // Dijet asymmetry AJ
      sprintf(histogramNamer,"dijet/dijetAsymmetry_C%d",iCentralityBin);
      fhDijetAsymmetry[iCentralityBin][knJetPtBins] = (TH1D*) fInputFile->Get(histogramNamer);
      
      // Dijet asymmetry xJ
      sprintf(histogramNamer,"dijet/dijetXj_C%d",iCentralityBin);
      fhDijetXj[iCentralityBin][knJetPtBins] = (TH1D*) fInputFile->Get(histogramNamer);
      
      // Asymmetries in jet pT bins
      for(int iJetPt = 0; iJetPt < knJetPtBins; iJetPt++){
        
        // Dijet asymmetry AJ in jet pT bins
        sprintf(histogramNamer,"dijet/dijetAsymmetry_C%dT%d",iCentralityBin,iJetPt);
        fhDijetAsymmetry[iCentralityBin][iJetPt] = (TH1D*) fInputFile->Get(histogramNamer);
        
        // Dijet asymmetry xJ in jet pT bins
        sprintf(histogramNamer,"dijet/dijetXj_C%dT%d",iCentralityBin,iJetPt);
        fhDijetXj[iCentralityBin][iJetPt] = (TH1D*) fInputFile->Get(histogramNamer);
        
      }
      
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
        sprintf(histogramNamer,"%s/%sPhi%s_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,fnTrackPtBins);
        fhTrackPhi[iTrackType][iCorrelationType][iCentralityBin][fnTrackPtBins] = (TH1D*) fInputFile->Get(histogramNamer);
        
        // pT integrated track eta
        sprintf(histogramNamer,"%s/%sEta%s_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,fnTrackPtBins);
        fhTrackEta[iTrackType][iCorrelationType][iCentralityBin][fnTrackPtBins] = (TH1D*) fInputFile->Get(histogramNamer);
        
        // pT integrated track eta-phi
        sprintf(histogramNamer,"%s/%sEtaPhi%s_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],fCompactCorrelationTypeString[iCorrelationType].Data(),iCentralityBin,fnTrackPtBins);
        if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCorrelationType][iCentralityBin][fnTrackPtBins] = (TH2D*) fInputFile->Get(histogramNamer);
        
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
      
      // Loop over asymmetry bins
      for(int iAsymmetry = 0; iAsymmetry <= fnAsymmetryBins; iAsymmetry++){
        if((iJetTrack >= kTrackInclusiveJet || !fProcessAsymmetryBins) && iAsymmetry != fnAsymmetryBins) continue; // No asymmetry bins for inclusive jet-track
        
        // Loop over centrality bins
        for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
          
          // Loop over track pT bins
          for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
            
            // Jet-track deltaPhi
            for(int iDeltaEta = 0; iDeltaEta < knDeltaEtaBins; iDeltaEta++){
              
              if(iDeltaEta > kWholeEta && iCorrelationType < kMixedEventNormalized) continue; // DeltaEta slicing not implemented for same and mixed event
              sprintf(histogramNamer,"%s/%sDeltaPhi%s%s_%sC%dT%d",fJetTrackHistogramNames[iJetTrack],fJetTrackHistogramNames[iJetTrack], fCompactCorrelationTypeString[iCorrelationType].Data(),fCompactDeltaEtaString[iDeltaEta], fAsymmetryBinName[iAsymmetry].Data(),iCentralityBin,iTrackPtBin);
              fhJetTrackDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][iDeltaEta] = (TH1D*) fInputFile->Get(histogramNamer);
            }
            
            // Jet-track deltaEtaDeltaPhi
            sprintf(histogramNamer,"%s/%sDeltaEtaDeltaPhi%s_%sC%dT%d",fJetTrackHistogramNames[iJetTrack],fJetTrackHistogramNames[iJetTrack], fCompactCorrelationTypeString[iCorrelationType].Data(),fAsymmetryBinName[iAsymmetry].Data(),iCentralityBin,iTrackPtBin);
            if(fLoad2DHistograms) fhJetTrackDeltaEtaDeltaPhi[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin] = (TH2D*) fInputFile->Get(histogramNamer);
            
            // DeltaPhi binning for deltaEta histogram
            for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
              sprintf(histogramNamer,"%s/%sDeltaEta%s%s_%sC%dT%d",fJetTrackHistogramNames[iJetTrack],fJetTrackHistogramNames[iJetTrack], fCompactCorrelationTypeString[iCorrelationType].Data(),fCompactDeltaPhiString[iDeltaPhi].Data(), fAsymmetryBinName[iAsymmetry].Data(),iCentralityBin,iTrackPtBin);
              fhJetTrackDeltaEta[iJetTrack][iCorrelationType][iAsymmetry][iCentralityBin][iTrackPtBin][iDeltaPhi] = (TH1D*) fInputFile->Get(histogramNamer);
            } // DeltaPhi loop
          } // Track pT loop
        } // Centrality loop
      } // Asymmetry loop
    } // Correlation type loop
  } // Jet-track correlation category loop
  
  // Load the jet shape histograms from the input file
  for(int iJetShape = 0; iJetShape < knJetShapeTypes; iJetShape++){
    
    // Loop over jet-track correlation categories
    for(int iJetTrack = 0; iJetTrack < knJetTrackCorrelations; iJetTrack++){
      if(!fLoadJetTrackCorrelations[iJetTrack]) continue;  // Only load the selected categories
      
      // Loop ovar asymmetry
      for(int iAsymmetry = 0; iAsymmetry <= fnAsymmetryBins; iAsymmetry++){
        if((iJetTrack >= kTrackInclusiveJet || !fProcessAsymmetryBins) && iAsymmetry != fnAsymmetryBins) continue; // No asymmetry bins for inclusive jet-track
        
        // Loop over centrality
        for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
          
          // Loop over track pT bins
          for(int iTrackPt = fFirstLoadedTrackPtBin; iTrackPt <= fLastLoadedTrackPtBin; iTrackPt++){
            sprintf(histogramNamer,"%s_%s/%s_%s_%sC%dT%d",fJetShapeHistogramName[iJetShape],fJetTrackHistogramNames[iJetTrack], fJetShapeHistogramName[iJetShape],fJetTrackHistogramNames[iJetTrack],fAsymmetryBinName[iAsymmetry].Data(),iCentrality,iTrackPt);
            fhJetShape[iJetShape][iJetTrack][iAsymmetry][iCentrality][iTrackPt] = (TH1D*) fInputFile->Get(histogramNamer);
            
          } // Track pT loop
        } // Centrality loop
      } // Asymmetry loop
    } // Jet-track correlation category loop
  } // Jet shape type loop
  
  // Load the jet pT closure histograms from a processed file
  if(fLoadJetPtClosureHistograms){
    
    // Loop over closure types (leading/subleading/inclusive)
    for(int iClosureType = 0; iClosureType < DijetHistograms::knClosureTypes; iClosureType++){
      
      // Centrality loop
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        
        // Loop over closure particles (quark/gluon/no selection)
        for(int iClosureParticle = 0; iClosureParticle < DijetHistograms::knClosureParticleTypes+1; iClosureParticle++){
          
          // Loop over generator level jet pT bins
          for(int iGenJetPt = 0; iGenJetPt < knGenJetPtBins; iGenJetPt++){
            sprintf(histogramNamer,"jetPtClosure_%s/jetPtClosure_%s%s_C%dT%d",fSingleJetHistogramName[iClosureType],fSingleJetHistogramName[iClosureType],fClosureParticleName[iClosureParticle],iCentrality,iGenJetPt);
            fhJetPtClosure[iClosureType][iGenJetPt][iCentrality][iClosureParticle] = (TH1D*) fInputFile->Get(histogramNamer);
            
          } // Generator level jet pT loop
        } // Closure particle type (quark/gluon) loop
      } // Centrality loop
      
      // Return back to main directory
      gDirectory->cd("../");
    } // Closure type (leading/subleading/inclusive) loop
    
  } // Opening jet pT closure histograms

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
void DijetHistogramManager::SetBinIndices(const int nBins, int *binIndices, const double *binBorders, const int iAxis, const bool setIndices){
  if(!setIndices) return;
  TH1D* hBinner;
  hBinner = FindHistogram(fInputFile,"trackLeadingJet",iAxis,0,0,0);
  for(int iBin = 0; iBin < nBins+1; iBin++){
    binIndices[iBin] = hBinner->GetXaxis()->FindBin(binBorders[iBin]);
  }
}

/*
 * Read the bin indices for given bin borders
 *
 *  Arguments:
 *   const int nBins = Number of bins for the indices
 *   double *copyBinBorders = Array to which a copy of bin borders is made
 *   int *binIndices = Array of integers to be filled with bin index information read from the file
 *   const double *binBorders = Array for bin borders that are searched from the file
 *   const int iAxis = Index of the axis used for reading bin indices
 *   const bool setIndices = true: Set both bin indices and bin borders, false: Set only bin borders
 */
void DijetHistogramManager::SetBinBordersAndIndices(const int nBins, double *copyBinBorders, int *binIndices, const double *binBorders, const int iAxis, const bool setIndices){
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
 *   const char* histogramName = Name of the histogram used to find bin indices
 *   const int nBins = Number of bins for the indices
 *   int *lowBinIndices = Array of integers to be filled with bin low edge index information read from the file
 *   int *highBinIndices = Array of integers to be filled with bin high edge index information read from the file
 *   const double *lowBinBorders = Array for low bin borders that are searched from the file
 *   const double *highBinBorders = Array for high bin borders that are searched from the file
 *   const int iAxis = Index of the axis used for reading bin indices
 */
void DijetHistogramManager::SetBinIndices(const char* histogramName, const int nBins, int *lowBinIndices, int *highBinIndices, const double *lowBinBorders, const double *highBinBorders, const int iAxis){
  TH1D* hBinner = FindHistogram(fInputFile,histogramName,iAxis,0,0,0);
  for(int iBin = 0; iBin < nBins; iBin++){
    lowBinIndices[iBin] = hBinner->GetXaxis()->FindBin(lowBinBorders[iBin]);
    highBinIndices[iBin] = hBinner->GetXaxis()->FindBin(highBinBorders[iBin]);
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
  SetBinIndices("trackLeadingJet",nBins,lowBinIndices,highBinIndices,lowBinBorders,highBinBorders,iAxis);
}

/*
 * Set up centrality bin borders and indices according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given centrality bins
 *  const double *binBorders = New bin borders for centrality
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void DijetHistogramManager::SetCentralityBins(const bool readBinsFromFile, const int nBins, const double *binBorders, const bool setIndices){
  
  // If bins are read from file, do not use the given bin borders
  if(readBinsFromFile){
    SetBinIndices(fnCentralityBins, fCentralityBinIndices, fCentralityBinBorders, 4, setIndices);
  } else { // If the given bin borders are use, update the number of bins and bin borders according to input
    if(nBins <= kMaxCentralityBins){
      fnCentralityBins = nBins;
      SetBinBordersAndIndices(fnCentralityBins,fCentralityBinBorders,fCentralityBinIndices,binBorders,4,setIndices);
    } else {
      cout << "Error! Too many centrality bins given. Maximum number is " << kMaxCentralityBins << ". Will not set bins." << endl;
    }
  }
}

/*
 * Set up track pT bin borders and indices according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given track pT bins
 *  const double *binBorders = New bin borders for track pT
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void DijetHistogramManager::SetTrackPtBins(const bool readBinsFromFile, const int nBins, const double *binBorders, const bool setIndices){
  
  // If bins are read from file, do not use the given bin borders
  if(readBinsFromFile){
    SetBinIndices(fnTrackPtBins, fTrackPtBinIndices, fTrackPtBinBorders, 0, setIndices);
  } else { // If the given bin borders are use, update the number of bins and bin borders according to input
    if(nBins <= kMaxTrackPtBins){
      fnTrackPtBins = nBins;
      SetBinBordersAndIndices(fnTrackPtBins,fTrackPtBinBorders,fTrackPtBinIndices,binBorders,0,setIndices);
    } else {
      cout << "Error! Too many track pT bins given. Maximum number is " << kMaxTrackPtBins << ". Will not set bins." << endl;
    }
  }
  
  // The track histograms have finer pT binning, so we need to use different bin indices for them
  if(setIndices){
    TH1D* hTrackPtBinner = FindHistogram(fInputFile,"track",0,0,0,0);
    for(int iTrackPt = 0; iTrackPt < fnTrackPtBins+1; iTrackPt++){
      fFineTrackPtBinIndices[iTrackPt] = hTrackPtBinner->GetXaxis()->FindBin(fTrackPtBinBorders[iTrackPt]);
    }
  }
}

/*
 * Set up deltaPhi bin indices according to provided bin borders
 */
void DijetHistogramManager::SetDeltaPhiBins(const bool readBinsFromFile, const double *lowBinBorders, const double *highBinBorders, TString deltaPhiStrings[knDeltaPhiBins], TString compactDeltaPhiStrings[knDeltaPhiBins], const bool setIndices){
  if(setIndices) SetBinIndices(knDeltaPhiBins,fLowDeltaPhiBinIndices,fHighDeltaPhiBinIndices,lowBinBorders,highBinBorders,1);
  if(!readBinsFromFile){
    for(int iDeltaPhi = 0; iDeltaPhi < knDeltaPhiBins; iDeltaPhi++){
      fDeltaPhiString[iDeltaPhi] = deltaPhiStrings[iDeltaPhi];
      fCompactDeltaPhiString[iDeltaPhi] = compactDeltaPhiStrings[iDeltaPhi];
      fLowDeltaPhiBinBorders[iDeltaPhi] = lowBinBorders[iDeltaPhi];
      fHighDeltaPhiBinBorders[iDeltaPhi] = highBinBorders[iDeltaPhi];
    }
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

// Setter for loading all leading jet histograms
void DijetHistogramManager::SetLoadAnyLeadingJetHistograms(const bool loadOrNot){
  fLoadSingleJets[kAnyLeadingJet] = loadOrNot;
}

// Setter for loading jet histograms
void DijetHistogramManager::SetLoadAllJets(const bool drawLeading, const bool drawSubleading, const bool drawAny, const bool drawAnyLeading){
  SetLoadLeadingJetHistograms(drawLeading);
  SetLoadSubleadingJetHistograms(drawSubleading);
  SetLoadAnyJetHistograms(drawAny);
  SetLoadAnyLeadingJetHistograms(drawAnyLeading);
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

// Setter for loading jet pT closure histograms
void DijetHistogramManager::SetLoadJetPtClosureHistograms(const bool loadOrNot){
  fLoadJetPtClosureHistograms = loadOrNot;
}

// Setter for drawn centrality bins
void DijetHistogramManager::SetCentralityBinRange(const int first, const int last){
  fFirstLoadedCentralityBin = first;
  fLastLoadedCentralityBin = last;
  
  // Sanity check for drawn centrality bins
  BinSanityCheck(fnCentralityBins,fFirstLoadedCentralityBin,fLastLoadedCentralityBin);
}

// Setter for drawn track pT bins
void DijetHistogramManager::SetTrackPtBinRange(const int first, const int last){
  fFirstLoadedTrackPtBin = first;
  fLastLoadedTrackPtBin = last;
  
  // Sanity check for drawn track pT bins
  BinSanityCheck(fnTrackPtBins,fFirstLoadedTrackPtBin,fLastLoadedTrackPtBin);
}

// Setter for processing asymmetry bins
void DijetHistogramManager::SetAsymmetryProcessing(const bool processAsymmetry){
  fProcessAsymmetryBins = processAsymmetry;
}

// Setter for preprocessing (only load and write same and mixed event)
void DijetHistogramManager::SetPreprocess(const bool preprocess){
  fPreprocess = preprocess;
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
  return fnCentralityBins;
}

// Getter for the number of track pT bins
int DijetHistogramManager::GetNTrackPtBins() const{
  return fnTrackPtBins;
}

// Getter for the number of jet pT bins
int DijetHistogramManager::GetNJetPtBins() const{
  return knJetPtBins;
}

// Getter for the number of dijet asymmetry AJ bins
int DijetHistogramManager::GetNAsymmetryBins() const{
  return fnAsymmetryBins;
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

// Getter for i:th track pT bin border
double DijetHistogramManager::GetJetPtBinBorder(const int iJetPt) const{
  return fJetPtBinBorders[iJetPt];
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

// Getter for z-vertex histogram
TH1D* DijetHistogramManager::GetHistogramVertexZWeighted() const{
  return fhVertexZWeighted;
}

// Getter for z-vertex histogram
TH1D* DijetHistogramManager::GetHistogramVertexZDijet() const{
  return fhVertexZDijet;
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
TH1D* DijetHistogramManager::GetHistogramJetPt(const int iJetType, const int iCentrality, int iAsymmetry) const{
  if(iAsymmetry == kMaxAsymmetryBins) iAsymmetry = fnAsymmetryBins;
  return fhJetPt[iJetType][iCentrality][iAsymmetry];
}

// Getter for jet phi histograms
TH1D* DijetHistogramManager::GetHistogramJetPhi(const int iJetType, const int iCentrality, int iAsymmetry) const{
  if(iAsymmetry == kMaxAsymmetryBins) iAsymmetry = fnAsymmetryBins;
  return fhJetPhi[iJetType][iCentrality][iAsymmetry];
}

// Getter for jet eta histograms
TH1D* DijetHistogramManager::GetHistogramJetEta(const int iJetType, const int iCentrality, int iAsymmetry) const{
  if(iAsymmetry == kMaxAsymmetryBins) iAsymmetry = fnAsymmetryBins;
  return fhJetEta[iJetType][iCentrality][iAsymmetry];
}

// Getter for 2D eta-phi histogram for jets
TH2D* DijetHistogramManager::GetHistogramJetEtaPhi(const int iJetType, const int iCentrality, int iAsymmetry) const{
  if(iAsymmetry == kMaxAsymmetryBins) iAsymmetry = fnAsymmetryBins;
  return fhJetEtaPhi[iJetType][iCentrality][iAsymmetry];
}

// Getters for dijet histograms

// Getter for dijet deltaPhi histograms
TH1D* DijetHistogramManager::GetHistogramDijetDeltaPhi(const int iCentrality) const{
  return fhDijetDphi[iCentrality];
}

// Getter for dijet asymmetry AJ histograms
TH1D* DijetHistogramManager::GetHistogramDijetAsymmetry(const int iCentrality, const int iJetPt) const{
  return fhDijetAsymmetry[iCentrality][iJetPt];
}

// Getter for dijet asymmetry xJ histograms
TH1D* DijetHistogramManager::GetHistogramDijetXj(const int iCentrality, const int iJetPt) const{
  return fhDijetXj[iCentrality][iJetPt];
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
TH1D* DijetHistogramManager::GetHistogramJetTrackDeltaPhi(const int iJetTrackCorrelation, const int iCorrelationType, int iAsymmetry, const int iCentrality, const int iTrackPt, const int iDeltaEta) const{
  if(iAsymmetry == kMaxAsymmetryBins) iAsymmetry = fnAsymmetryBins;
  return fhJetTrackDeltaPhi[iJetTrackCorrelation][iCorrelationType][iAsymmetry][iCentrality][iTrackPt][iDeltaEta];
}

// Getter for deltaEta between jet and track
TH1D* DijetHistogramManager::GetHistogramJetTrackDeltaEta(const int iJetTrackCorrelation, const int iCorrelationType, int iAsymmetry, const int iCentrality, const int iTrackPt, const int iDeltaPhiRegion) const{
  if(iAsymmetry == kMaxAsymmetryBins) iAsymmetry = fnAsymmetryBins;
  return fhJetTrackDeltaEta[iJetTrackCorrelation][iCorrelationType][iAsymmetry][iCentrality][iTrackPt][iDeltaPhiRegion];
}

// Getter for deltaEta and deltaPhi between jet and track
TH2D* DijetHistogramManager::GetHistogramJetTrackDeltaEtaDeltaPhi(const int iJetTrackCorrelation, const int iCorrelationType, int iAsymmetry, const int iCentrality, const int iTrackPt) const{
  if(iAsymmetry == kMaxAsymmetryBins) iAsymmetry = fnAsymmetryBins;
  return fhJetTrackDeltaEtaDeltaPhi[iJetTrackCorrelation][iCorrelationType][iAsymmetry][iCentrality][iTrackPt];
}

// Getters for jet shape histograms
TH1D* DijetHistogramManager::GetHistogramJetShape(const int iJetShapeType, const int iJetTrackCorrelation, int iAsymmetry, const int iCentrality, const int iTrackPt) const{
  if(iAsymmetry == kMaxAsymmetryBins) iAsymmetry = fnAsymmetryBins;
  return fhJetShape[iJetShapeType][iJetTrackCorrelation][iAsymmetry][iCentrality][iTrackPt];
}

// Getter for jet pT closure histograms
TH1D* DijetHistogramManager::GetHistogramJetPtClosure(const int iClosureType, const int iGenPtBin, const int iCentrality, const int iClosureParticle) const{
  return fhJetPtClosure[iClosureType][iGenPtBin][iCentrality][iClosureParticle];
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
 *   int bin6 = Sixth bin index for the histogram
 *
 *   return: Histogram corresponding to given name and bins
 */
TH1D* DijetHistogramManager::GetOneDimensionalHistogram(TString name, int bin1, int bin2, int bin3, int bin4, int bin5, int bin6) const{
  if(name.EqualTo("vertexz",TString::kIgnoreCase) || name.EqualTo("fhvertexz",TString::kIgnoreCase)) return GetHistogramVertexZ();
  if(name.EqualTo("events",TString::kIgnoreCase) || name.EqualTo("fhevents",TString::kIgnoreCase)) return GetHistogramEvents();
  if(name.EqualTo("trackcuts",TString::kIgnoreCase) || name.EqualTo("fhtrackcuts",TString::kIgnoreCase)) return GetHistogramTrackCuts();
  if(name.EqualTo("centrality",TString::kIgnoreCase) || name.EqualTo("fhcentrality",TString::kIgnoreCase)) return GetHistogramCentrality();
  if(name.EqualTo("centralitydijet",TString::kIgnoreCase) || name.EqualTo("fhcentralitydijet",TString::kIgnoreCase)) return GetHistogramCentralityDijet();
  if(name.EqualTo("jetpt",TString::kIgnoreCase) || name.EqualTo("fhjetpt",TString::kIgnoreCase)) return GetHistogramJetPt(bin1,bin2,bin3);
  if(name.EqualTo("jetphi",TString::kIgnoreCase) || name.EqualTo("fhjetphi",TString::kIgnoreCase)) return GetHistogramJetPhi(bin1,bin2,bin3);
  if(name.EqualTo("jeteta",TString::kIgnoreCase) || name.EqualTo("fhjeteta",TString::kIgnoreCase)) return GetHistogramJetEta(bin1,bin2,bin3);
  if(name.EqualTo("dijetdeltaphi",TString::kIgnoreCase) || name.EqualTo("dijetdphi",TString::kIgnoreCase) || name.EqualTo("fhdijetdphi",TString::kIgnoreCase)) return GetHistogramDijetDeltaPhi(bin1);
  if(name.EqualTo("dijetasymmetry",TString::kIgnoreCase) || name.EqualTo("fhdijetasymmetry",TString::kIgnoreCase)) return GetHistogramDijetAsymmetry(bin1);
  if(name.EqualTo("dijetxj",TString::kIgnoreCase) || name.EqualTo("fhdijetxj",TString::kIgnoreCase)) return GetHistogramDijetXj(bin1);
  if(name.EqualTo("trackpt",TString::kIgnoreCase) || name.EqualTo("fhtrackpt",TString::kIgnoreCase)) return GetHistogramTrackPt(bin1,bin2,bin3);
  if(name.EqualTo("trackphi",TString::kIgnoreCase) || name.EqualTo("fhtrackphi",TString::kIgnoreCase)) return GetHistogramTrackPhi(bin1,bin2,bin3,bin4);
  if(name.EqualTo("tracketa",TString::kIgnoreCase) || name.EqualTo("fhtracketa",TString::kIgnoreCase)) return GetHistogramTrackEta(bin1,bin2,bin3,bin4);
  if(name.EqualTo("jettrackdeltaphi",TString::kIgnoreCase) || name.EqualTo("fhjettrackdeltaphi",TString::kIgnoreCase)) return GetHistogramJetTrackDeltaPhi(bin1,bin2,bin3,bin4,bin5,bin6);
  if(name.EqualTo("jettrackdeltaeta",TString::kIgnoreCase) || name.EqualTo("fhjettrackdeltaeta",TString::kIgnoreCase)) return GetHistogramJetTrackDeltaEta(bin1,bin2,bin3,bin4,bin5,bin6);
  if(name.EqualTo("jetshape",TString::kIgnoreCase) || name.EqualTo("fhjetshape",TString::kIgnoreCase)) return GetHistogramJetShape(bin1,bin2,bin3,bin4,bin5);
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
 *   int bin5 = Fifth bin index for the histogram
 *
 *   return: Histogram corresponding to given name and bins
 */
TH2D* DijetHistogramManager::GetTwoDimensionalHistogram(TString name, int bin1, int bin2, int bin3, int bin4, int bin5) const{
  if(name.EqualTo("jetetaphi",TString::kIgnoreCase) || name.EqualTo("fhjetetaphi",TString::kIgnoreCase)) return GetHistogramJetEtaPhi(bin1,bin2);
  if(name.EqualTo("dijetleadingvssubleadingpt",TString::kIgnoreCase) || name.EqualTo("fhdijetleadingvssubleadingpt",TString::kIgnoreCase)) return GetHistogramDijetLeadingVsSubleadingPt(bin1);
  if(name.EqualTo("tracketaphi",TString::kIgnoreCase) || name.EqualTo("fhtracketaphi",TString::kIgnoreCase)) return GetHistogramTrackEtaPhi(bin1,bin2,bin3,bin4);
  if(name.EqualTo("jettrackdeltaetadeltaphi",TString::kIgnoreCase) || name.EqualTo("fhjettrackdeltaetadeltaphi",TString::kIgnoreCase)) return GetHistogramJetTrackDeltaEtaDeltaPhi(bin1,bin2,bin3,bin4,bin5);
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

// Getter for integral over leading jet pT. Include the overflow bin in the integral. There is no jet pT limit
// in for the correlation but the histogram has limited range, so to get all the jets into normalization we
// need to include the overflow bin in the intagral.
double DijetHistogramManager::GetPtIntegral(const int iCentrality, int iAsymmetry) const{
  if(iAsymmetry == kMaxAsymmetryBins) iAsymmetry = fnAsymmetryBins;
  return fhJetPt[kLeadingJet][iCentrality][iAsymmetry]->Integral(1,fhJetPt[kLeadingJet][iCentrality][iAsymmetry]->GetNbinsX()+1,"width");
}

// Getter for integral over all leading jets with pT > 120 GeV in a given centrality bin. Include the overflow bin in the integral
double DijetHistogramManager::GetAnyLeadingJetPtIntegral(const int iCentrality) const{
  return fhJetPt[kAnyLeadingJet][iCentrality][fnAsymmetryBins]->Integral(fhJetPt[kAnyLeadingJet][iCentrality][fnAsymmetryBins]->FindBin(120),fhJetPt[kAnyLeadingJet][iCentrality][fnAsymmetryBins]->GetNbinsX()+1,"width");
}

/*
 * Getter for integral over inclusive jet pT over 120 GeV. Include the overflow bin in the integral
 *
 *  const int iCentrality = Centrality bin
 *  const double minPt = Jet pT above which the integral is calculated
 */
double DijetHistogramManager::GetInclusiveJetPtIntegral(const int iCentrality, const double minPt) const{
  return fhJetPt[kAnyJet][iCentrality][fnAsymmetryBins]->Integral(fhJetPt[kAnyJet][iCentrality][fnAsymmetryBins]->FindBin(minPt+0.001), fhJetPt[kAnyJet][iCentrality][fnAsymmetryBins]->GetNbinsX()+1,"width");
}

// Setter for avoiding possible peaks in mixed event distribution
void DijetHistogramManager::SetAvoidMixingPeak(const bool avoid){
  fAvoidMixingPeak = avoid;
}

// Setter for avoiding possible peaks in mixed event distribution
void DijetHistogramManager::SetImproviseMixing(const bool improvise){
  fImproviseMixing = improvise;
}

// Default fit range used to normalize the mixed event
void DijetHistogramManager::SetDefaultMixingDeltaEtaFitRange(const double fitRange){
  fDefaultMixingDeltaEtaFitRange = fitRange;
}
