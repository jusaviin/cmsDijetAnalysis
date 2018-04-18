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
 *
 *  Arguments:
 *    int debugLevel = Amount of debug messages printed out. 0 = none, 1 = some, 2 = all.
 */
void DijetAnalyzer::RunAnalysis(int debugLevel){
  
  TFile *inputFile;
  ForestReader *treeReader = new ForestReader();
  
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
  const double pTmincut = fCard->Get("MinPtCut");      // Original 120
  const double leadingjetcut = fCard->Get("MinPtCut");  // Original 120
  const double subleadingjetcut = fCard->Get("SubleadingPtCut"); // Original 50
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
  double filler2D[2];
  double filler3D[3];
  
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
      fHistograms->fhVertexZ->Fill(vz);     // z vertex distribution from all events
      fHistograms->fhEvents->Fill(1);       // All the events looped over
      fHistograms->fhCentrality->Fill(centrality);   // Centrality filled from all events
      
      // Apply the vz cut
      if(TMath::Abs(vz) > vzCut) continue;
      fHistograms->fhEvents->Fill(2);       // Events with z vertex within accepteble range
      
      // Reset dijet flags
      twoJetsFound = false;
      dijetFound = false;
      
      // Search for leading jet
      for(int jetIndex = 0; jetIndex < treeReader->GetNJets(); jetIndex++) {
        jetPt = treeReader->GetJetPt(jetIndex);
        if(TMath::Abs(treeReader->GetJetEta(jetIndex)) >= searchetacut) continue;
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
           (TMath::Abs(treeReader->GetJetEta(highestIndex)) >= etacut)||
           (TMath::Abs(dphi) <= dphicut)){
          dijetFound = false;
        }
      } // End of dijet cuts
      
      // If a dijet is found, fill some information to fHistograms
      if(dijetFound){
        fHistograms->fhEvents->Fill(3);
        
        // Put the centralities to fillers
        filler2D[1] = centrality;
        filler3D[2] = centrality;
        
        // Fill dijet deltaPhi histogram
        filler2D[0] = fabs(dphi);
        fHistograms->fhDijetDphi->Fill(filler2D);
        
        // Calculate the asymmetry
        Aj = (treeReader->GetJetPt(highestIndex) - treeReader->GetJetPt(secondHighestIndex))/(treeReader->GetJetPt(highestIndex) + treeReader->GetJetPt(secondHighestIndex));
        
        // Fill the asymmetry histogram
        filler2D[0] = Aj;
        fHistograms->fhDijetAsymmetry->Fill(filler2D);
        
        // Fill the aymmetry versus delta phi histogram
        filler3D[0] = fabs(dphi);
        filler3D[1] = Aj;
        fHistograms->fhDijetAsymmetryVsDphi->Fill(filler3D);
        
        // Fill the leading jet pT histogram
        filler2D[0] = treeReader->GetJetPt(highestIndex);
        fHistograms->fhLeadingJetPt->Fill(filler2D);
        
        // Fill the subleading jet pT histogram
        filler2D[0] = treeReader->GetJetPt(secondHighestIndex);
        fHistograms->fhSubleadingJetPt->Fill(filler2D);
        
        // Fill the histogram for leading jet pT versus subleading jet pT
        filler3D[0] = treeReader->GetJetPt(highestIndex);
        filler3D[1] = treeReader->GetJetPt(secondHighestIndex);
        fHistograms->fhDijetLeadingVsSubleadingPt->Fill(filler3D);
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
DijetHistograms* DijetAnalyzer::GetHistograms(){
  return fHistograms;
}
