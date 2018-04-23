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
  
  const double etacut = fCard->Get("EtaCut");
  const double searchetacut = fCard->Get("SearchEtaCut");
  const double pTmaxcut = fCard->Get("MaxPtCut");
  const double pTmincut = fCard->Get("MinPtCut");
  const double leadingjetcut = fCard->Get("MinPtCut");
  const double subleadingjetcut = fCard->Get("SubleadingPtCut");
  const double dphicut = fCard->Get("DeltaPhiCut");
  
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
        if(TMath::Abs(treeReader->GetJetEta(jetIndex)) >= searchetacut) continue;
        
        // Fill the histogram for all jets within eta range
        if(TMath::Abs(treeReader->GetJetEta(jetIndex)) < etacut){
          
          // Fill the axes in correct order
          fillerAnyJet[0] = jetPt;                             // Axis 0 = any jet pT
          fillerAnyJet[1] = treeReader->GetJetPhi(jetIndex);   // Axis 1 = any jet phi
          fillerAnyJet[2] = treeReader->GetJetEta(jetIndex);   // Axis 2 = any jet eta
          fillerAnyJet[3] = centrality;                        // Axis 3 = centrality
          fHistograms->fhAnyJet->Fill(fillerAnyJet);           // Fill the data point to histogram
          
        }
        
        if(jetPt <= leadingjetcut) continue;
        if(jetPt > leadingPt){
          leadingPt = jetPt;
          highestIndex = jetIndex;
        }
      } // End of search for leading jet loop
      
      // Search for subleading jet
      for(int jetIndex = 0 ; jetIndex < treeReader->GetNJets(); jetIndex++){
        jetPt = treeReader->GetJetPt(jetIndex);
        if(jetIndex == highestIndex) continue;
        if(TMath::Abs(treeReader->GetJetEta(jetIndex)) >= searchetacut) continue;
        if(jetPt <= subleadingjetcut) continue;
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
        
        
        if((treeReader->GetJetPt(highestIndex) >= pTmaxcut) ||
           (treeReader->GetJetPt(highestIndex) <= pTmincut) ||
           (TMath::Abs(treeReader->GetJetEta(highestIndex)) >= etacut) ||
           (TMath::Abs(treeReader->GetJetEta(secondHighestIndex)) >= etacut)||
           (TMath::Abs(dphi) <= dphicut)){
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
