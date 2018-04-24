// Class for the main analysis algorithms for the leading-subleading jet analysis

// Root includes
#include <TFile.h>
#include <TMath.h>

// Own includes
#include "DijetAnalyzer.h"
#include "ForestReader.h"

using namespace std;

/*
 * Default constructor
 */
DijetAnalyzer::DijetAnalyzer() :
  fFileNames(0),
  fCard(0),
  fHistograms(0)
{
  // Default constructor
  fHistograms = new DijetHistograms();
  fHistograms->CreateHistograms();
}

/*
 * Custom constructor
 */
DijetAnalyzer::DijetAnalyzer(std::vector<TString> fileNameVector, ConfigurationCard *newCard) :
  fFileNames(fileNameVector),
  fCard(newCard),
  fHistograms(0)
{
  // Custom constructor
  fHistograms = new DijetHistograms();
  fHistograms->CreateHistograms();
}

/*
 * Copy constructor
 */
DijetAnalyzer::DijetAnalyzer(const DijetAnalyzer& in) :
  fFileNames(in.fFileNames),
  fCard(in.fCard),
  fHistograms(in.fHistograms)
{
  // Copy constructor
}

/*
 * Assingment operator
 */
DijetAnalyzer& DijetAnalyzer::operator=(const DijetAnalyzer& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fFileNames = in.fFileNames;
  fCard = in.fCard;
  fHistograms = in.fHistograms;
  
  return *this;
}

/*
 * Destructor
 */
DijetAnalyzer::~DijetAnalyzer(){
  // destructor
  delete fHistograms;
}

/*
 * Main analysis loop
 */
void DijetAnalyzer::RunAnalysis(){
  
  TFile *inputFile;
  ForestReader *treeReader = new ForestReader(fCard->Get("DataType"));
  
  //****************************************
  //        Event selection cuts
  //****************************************
  
  const double vzCut = fCard->Get("ZVertexCut");
  
  //****************************************
  //        Dijet selection cuts
  //****************************************
  
  const double jetEtaCut = fCard->Get("JetEtaCut");
  const double jetSearchEtaCut = fCard->Get("SearchEtaCut");
  const double jetMaximumPtCut = fCard->Get("MaxPtCut");
  const double leadingJetMinPtCut = fCard->Get("MinPtCut");
  const double subleadingJetMinPtCut = fCard->Get("SubleadingPtCut");
  const double deltaPhiCut = fCard->Get("DeltaPhiCut");
  const double minimumMaxTrackPtFraction = fCard->Get("MinMaxTrackPtFraction");
  const double maximumMaxTrackPtFraction = fCard->Get("MaxMaxTrackPtFraction");
  
  //****************************************
  //        Traxk selection cuts
  //****************************************
  
  const double trackEtaCut = fCard->Get("TrackEtaCut");
  const double trackMinPtCut = fCard->Get("MinTrackPtCut");
  
  //****************************************
  //            All cuts set!
  //****************************************

  // Event variables
  int nEvents = 0;
  bool dijetFound = false;
  bool twoJetsFound = false;
  double vz = 0;
  double centrality = 0;
  
  // Variables for jets
  double Aj = -99;
  double leadingPt = 0;
  double subleadingPt = 0;
  int secondHighestIndex = -1;
  int highestIndex = -1;
  double jetPt = 0;
  double dphi = 0;
  
  // Variables for tracks
  double trackPt;       // Track pT
  double trackEta;      // Track eta
  double trackPhi;      // Track phi
  double trackRMin;     // Minimum distance to a jet
  
  // File name helper variables
  TString currentFile;
  
  // Fillers for THnSparses
  double fillerAnyJet[4];
  double fillerDijet[9];
  
  // Amount of debugging messages
  int debugLevel = fCard->Get("DebugLevel");
  
  // Loop over files
  for(int iFile = 0; iFile < (int) fFileNames.size(); iFile++) {
    
    // Find the filename
    currentFile = fFileNames.at(iFile);
    
    // Open the file and check that everything goes fine
    inputFile = TFile::Open(currentFile);

    // Check that the file exists
    if(!inputFile){
      cout << "Warning! Could not open the file: " << currentFile.Data() << endl;
      continue;
    }
    
    // Check that the file is not zombie
    if(inputFile->IsZombie()){
      cout << "Warning! The following file is a zombie: " << currentFile.Data() << endl;
      continue;
    }
    
    // Debug message, if wanted
    if(debugLevel > 0) cout << "Reading from file: " << currentFile.Data() << endl;
    
    // If file is good, read the forest from the file
    treeReader->ReadForestFromFile(inputFile);
    nEvents = treeReader->GetNEvents();
    
    // Event loop
    for(int iEvent = 0; iEvent < nEvents; iEvent++){
      
      // Read the event to memory
      treeReader->GetEvent(iEvent);
      
      // Get vz and centrality information from all events
      vz = treeReader->GetVz();
      centrality = treeReader->GetCentrality();
      fHistograms->fhVertexZ->Fill(vz);                   // z vertex distribution from all events
      fHistograms->fhEvents->Fill(DijetHistograms::kAll); // All the events looped over
      fHistograms->fhCentrality->Fill(centrality);        // Centrality filled from all events
      
      //  ===== Apply all the event quality cuts =====
      
      // Cut for primary vertex. Only applied for data.
      if(treeReader->GetPrimaryVertexFilterBit() == 0) continue;
      fHistograms->fhEvents->Fill(DijetHistograms::kPrimaryVertex);
      
      // Cut for HB/HE noise. Only applied for data.
      if(treeReader->GetHBHENoiseFilterBit() == 0) continue;
      fHistograms->fhEvents->Fill(DijetHistograms::kHBHENoise);
      
      // Cut for collision event selection. Only applied for PbPb data.
      if(treeReader->GetCollisionEventSelectionFilterBit() == 0) continue;
      fHistograms->fhEvents->Fill(DijetHistograms::kCollisionEventSelection);
      
      // Cut for beam scraping. Only applied for pp data.
      if(treeReader->GetBeamScrapingFilterBit() == 0) continue;
      fHistograms->fhEvents->Fill(DijetHistograms::kBeamScraping);
      
      // Cut for calirimeter jet quality. Only applied for data.
      if(treeReader->GetCaloJetFilterBit() == 0) continue;
      fHistograms->fhEvents->Fill(DijetHistograms::kCaloJet);
      
      // Cut for vertex z-position
      if(TMath::Abs(vz) > vzCut) continue;
      fHistograms->fhEvents->Fill(DijetHistograms::kVzCut);
      
      // ===== Event quality cuts applied =====
      
      // Reset the variables used in dijet finding
      twoJetsFound = false;
      dijetFound = false;
      highestIndex = -1;
      secondHighestIndex = -1;
      leadingPt = 0;
      subleadingPt = 0;
      
      // Search for leading jet and fill histograms for all jets within the eta vut
      for(int jetIndex = 0; jetIndex < treeReader->GetNJets(); jetIndex++) {
        jetPt = treeReader->GetJetPt(jetIndex);
        if(TMath::Abs(treeReader->GetJetEta(jetIndex)) >= jetSearchEtaCut) continue; // Cut for search eta range
        if(minimumMaxTrackPtFraction > treeReader->GetJetMaxTrackPt()/treeReader->GetJetRawPt()) continue; // Cut for jets with only very low pT particle
        if(maximumMaxTrackPtFraction < treeReader->GetJetMaxTrackPt()/treeReader->GetJetRawPt()) continue; // Cut for jets where all the pT is taken by one track
        
        // Fill the histogram for all jets within eta range
        if(TMath::Abs(treeReader->GetJetEta(jetIndex)) < jetEtaCut){
          
          // Fill the axes in correct order
          fillerAnyJet[0] = jetPt;                             // Axis 0 = any jet pT
          fillerAnyJet[1] = treeReader->GetJetPhi(jetIndex);   // Axis 1 = any jet phi
          fillerAnyJet[2] = treeReader->GetJetEta(jetIndex);   // Axis 2 = any jet eta
          fillerAnyJet[3] = centrality;                        // Axis 3 = centrality
          fHistograms->fhAnyJet->Fill(fillerAnyJet);           // Fill the data point to histogram
          
        }
        
        if(jetPt <= leadingJetMinPtCut) continue; // Minimum leading jet pT cut
        if(jetPt > leadingPt){
          leadingPt = jetPt;
          highestIndex = jetIndex;
        }
      } // End of search for leading jet loop
      
      // Search for subleading jet
      for(int jetIndex = 0 ; jetIndex < treeReader->GetNJets(); jetIndex++){
        jetPt = treeReader->GetJetPt(jetIndex);
        if(jetIndex == highestIndex) continue; // Do not consider leading particle
        if(TMath::Abs(treeReader->GetJetEta(jetIndex)) >= jetSearchEtaCut) continue; // Cut for search eta range
        if(minimumMaxTrackPtFraction > treeReader->GetJetMaxTrackPt()/treeReader->GetJetRawPt()) continue; // Cut for jets with only very low pT particle
        if(maximumMaxTrackPtFraction < treeReader->GetJetMaxTrackPt()/treeReader->GetJetRawPt()) continue; // Cut for jets where all the pT is taken by one track
        if(jetPt <= subleadingJetMinPtCut) continue; // Minimum subleading jet pT cut
        if(jetPt > subleadingPt){
          subleadingPt = jetPt;
          secondHighestIndex = jetIndex;
        }
      }  //End of subleading jet search
      
      // Check if at least two jets were found
      if(highestIndex > -1 && secondHighestIndex > -1) twoJetsFound = true;
      
      // Only apply the dijet cuts for events with at least two jets
      if(twoJetsFound){
        
        dijetFound = true;
        dphi =  treeReader->GetJetPhi(highestIndex) - treeReader->GetJetPhi(secondHighestIndex);
        if(dphi < 0) dphi = -dphi;
        if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
        
        
        if((treeReader->GetJetPt(highestIndex) >= jetMaximumPtCut) ||              // Maximum leading jet pT cut
           (TMath::Abs(treeReader->GetJetEta(highestIndex)) >= jetEtaCut) ||       // Leading jet eta cut
           (TMath::Abs(treeReader->GetJetEta(secondHighestIndex)) >= jetEtaCut)||  // Subleading jet eta cut
           (TMath::Abs(dphi) <= deltaPhiCut)){                                     // DeltaPhi cut
          dijetFound = false;
        }
      } // End of dijet cuts
      
      // If a dijet is found, fill some information to fHistograms
      if(dijetFound){
        fHistograms->fhEvents->Fill(DijetHistograms::kDijet);
        
        // Calculate the asymmetry
        Aj = (treeReader->GetJetPt(highestIndex) - treeReader->GetJetPt(secondHighestIndex))/(treeReader->GetJetPt(highestIndex) + treeReader->GetJetPt(secondHighestIndex));
        
        // Fill the dijet histogram axes in correct order
        fillerDijet[0] = treeReader->GetJetPt(highestIndex);        // Axis 0: Leading jet pT
        fillerDijet[1] = treeReader->GetJetPhi(highestIndex);       // Axis 1: Leading jet phi
        fillerDijet[2] = treeReader->GetJetEta(highestIndex);       // Axis 2: Leading jet eta
        fillerDijet[3] = treeReader->GetJetPt(secondHighestIndex);  // Axis 3: Subleading jet pT
        fillerDijet[4] = treeReader->GetJetPhi(secondHighestIndex); // Axis 4: Subleading jet phi
        fillerDijet[5] = treeReader->GetJetEta(secondHighestIndex); // Axis 5: Subleading jet eta
        fillerDijet[6] = fabs(dphi);                                // Axis 6: deltaPhi
        fillerDijet[7] = Aj;                                        // Axis 7: Asymmetry
        fillerDijet[8] = centrality;                                // Axis 8: Centrality
        fHistograms->fhDijet->Fill(fillerDijet);                    // Fill the data point to dijet histogram

      }
      
      // Correlate jets with tracks
      int nTracks = treeReader->GetNTracks();
      for(int iTrack = i; iTrack <nTracks; iTrack++){
        trackEta = treeReader->GetTrackEta(iTrack);
        
        // Apply cuts for tracks
        if(fabs(trackEta) > trackEtaCut) continue;                // Eta cut
        if(treeReader->GetTrackHighPurity(iTrack) != 1) continue; // High purity cut
        
        trackPt = treeReader->GetTrackPt(iTrack);
        trackPhi = treeReader->GetTrackPhi(iTrack);
        
      } // Loop over tracks
      
    } // Event loop
    
    // Close the input file after the event has been read
    inputFile->Close();
    
  } // File loop
  
  delete treeReader;  // Delete the created tree reader after the analysis is done
}

/*
 * Getter for dijet histograms
 */
DijetHistograms* DijetAnalyzer::GetHistograms() const{
  return fHistograms;
}
