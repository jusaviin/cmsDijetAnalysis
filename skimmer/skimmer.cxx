// C++ includes
#include <iostream>
#include <fstream>
#include <vector>

// Root includes
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"

using namespace std;

/*
 * File list reader
 *
 *  Arguments:
 *    std::vector<TString> &fileNameVector = Vector filled with filenames found in the file
 *    TString fileNameFile = Text file containing one analysis file name in each line
 *    int debug = Level of debug messages shown
 *    int locationIndex = Where to find analysis files: 0 = Purdue EOS, 1 = CERN EOS, 2 = Use xrootd to find the data
 *    bool runLocal = True: Local run mode. False: Crab run mode
 */
void ReadFileList(std::vector<TString> &fileNameVector, TString fileNameFile, int debug, int locationIndex, bool runLocal)
{
  
  // Possible location for the input files
  const char *fileLocation[] = {"root://xrootd.rcac.purdue.edu/","root://eoscms.cern.ch/","root://cmsxrootd.fnal.gov/"};
  
  // Set up the file names file for reading
  ifstream file_stream(fileNameFile);
  std::string line;
  fileNameVector.clear();
  if( debug > 0 ) std::cout << "Open file " << fileNameFile.Data() << " to extract files to run over" << std::endl;
  
  // Open the file names file for reading
  if( file_stream.is_open() ) {
    if( debug > 0) std::cout << "Opened " << fileNameFile.Data() << " for reading" << std::endl;
    int lineNumber = 0;
    
    // Loop over the lines in the file
    while( !file_stream.eof() ) {
      getline(file_stream, line);
      if( debug > 0) std::cout << lineNumber << ": " << line << std::endl;
      TString lineString(line);
      
      // Put all non-empty lines to file names vector
      if( lineString.CompareTo("", TString::kExact) != 0 ) {
        
        if(runLocal){
          // For local running, it is assumed that the file name is directly the centents of the line
          fileNameVector.push_back(lineString);
          
        } else {
          // For crab running, the line will have format ["file1", "file2", ... , "fileN"]
          TObjArray *fileNameArray = lineString.Tokenize(" ");  // Tokenize the string from every ' ' character
          int numberOfFiles = fileNameArray->GetEntries();
          TObjString *currentFileNameObject;
          TString currentFileName;
          for(int i = 0; i < numberOfFiles; i++){   // Loop over all the files in the array
            currentFileNameObject = (TObjString *)fileNameArray->At(i);
            currentFileName = currentFileNameObject->String();
            
            // Strip unwanted characters
            currentFileName.Remove(TString::kBoth,'['); // Remove possible parantheses
            currentFileName.Remove(TString::kBoth,']'); // Remove possible parantheses
            currentFileName.Remove(TString::kBoth,','); // Remove commas
            currentFileName.Remove(TString::kBoth,'"'); // Remove quotation marks
            
            // After stripping characters not belonging to the file name, we can add the file to list
            currentFileName.Prepend(fileLocation[locationIndex]);  // If not running locally, we need to give xrootd path
            fileNameVector.push_back(currentFileName);
          }
        }
        
      } // Empty line if
      
      
      lineNumber++;
    } // Loop over lines in the file
    
    // If cannot read the file, give error and end program
  } else {
    std::cout << "Error, could not open " << fileNameFile.Data() << " for reading" << std::endl;
    assert(0);
  }
}

/*
 *  Convert string to boolean value
 */
bool checkBool(string str) {
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  std::istringstream is(str);
  bool b;
  is >> std::boolalpha >> b;
  return b;
}

int main(int argc,char *argv[]){
  
  //==== Read arguments =====
  if ( argc<5 ) {
    cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    cout<<"+ Usage of the macro: " << endl;
    cout<<"+  "<<argv[0]<<" [fileNameFile] [isMC] [outputFileName] [fileLocation] <runLocal>"<<endl;
    cout<<"+  fileNameFile: Text file containing the list of files used in the analysis. For crab analysis a job id should be given here." <<endl;
    cout<<"+  isMC: 1 for MC, 0 for data." <<endl;
    cout<<"+  outputFileName: .root file to which the histograms are written." <<endl;
    cout<<"+  fileLocation: Where to find analysis files: 0 = Purdue EOS, 1 = CERN EOS, 2 = Use xrootd to find the data." << endl;
    cout<<"+  runLocal: True: Search input files from local machine. False (default): Search input files from grid with xrootd." << endl;
    cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    cout << endl << endl;
    exit(1);
  }
  
  // First, check if we are supposed to run locally or on crab
  bool runLocal = false;
  if(argc >= 6) runLocal = checkBool(argv[5]);
  
  // Find the file list name depending on if we run locally or on crab
  TString fileNameFile;
  if(runLocal){
    fileNameFile = argv[1];
  } else {
    fileNameFile = Form("job_input_file_list_%d.txt",atoi(argv[1]));
  }
  
  // Read the other command line arguments
  int isMC = atoi(argv[2]);
  TString outputFileName = argv[3];
  const int fileSearchIndex = atoi(argv[4]);
  
  // Read the file names used for skimming to a vector
  std::vector<TString> fileNameVector;
  fileNameVector.clear();
  ReadFileList(fileNameVector,fileNameFile,1,fileSearchIndex,runLocal);
  
  // Maximum size of arrays
  const Int_t nMaxJet = 250;        // Maximum number of jets in an event
  const Int_t nMaxTrack = 60000;    // Maximum number of tracks in an event
  
  //jetTree[0] = (TTree*)inputFile->Get("ak4CaloJetAnalyzer/t");
  //jetTree[1] = (TTree*)inputFile->Get("ak4PFJetAnalyzer/t");
  
  // Define trees to be read from the files
  const int nJetTrees = 4; // 4 For PbPb, 2 for pp
  TChain *heavyIonTree = new TChain("hiEvtAnalyzer/HiTree");
  TChain *hltTree = new TChain("hltanalysis/HltTree");
  TChain *skimTree = new TChain("skimanalysis/HltTree");
  TChain *jetTree[nJetTrees];
  jetTree[0] = new TChain("akPu4CaloJetAnalyzer/t");     // For pp: ak4CaloJetAnalyzer/t
  jetTree[1] = new TChain("akCs4PFJetAnalyzer/t");       // For pp: ak4PFJetAnalyzer/t
  jetTree[2] = new TChain("akPu4PFJetAnalyzer/t");       // For pp: remove
  jetTree[3] = new TChain("akFlowPuCs4PFJetAnalyzer/t"); // For pp: remove
  TChain *trackTree = new TChain("ppTrack/trackTree");
  TChain *genTrackTree = new TChain("HiGenParticleAna/hi");
  TChain *particleFlowCandidateTree = new TChain("pfcandAnalyzer/pfTree");
  
  // All the branches and leaves come in arrat of two, one for input and one for output
  
  // Branches for heavy ion tree
  TBranch *hiVzBranch;            // Branch for vertex z-position
  TBranch *hiBinBranch;           // Branch for centrality
  TBranch *ptHatBranch;           // Branch for pT hat
  TBranch *eventWeightBranch;     // Branch for jet weight for 2018 MC
  
  // Leaves for heavy ion tree
  Float_t vertexZ;      // Vertex z-position
  Int_t hiBin;          // HiBin = Centrality percentile * 2
  Float_t ptHat;        // pT hat
  Float_t eventWeight;  // jet weight in the 2018 MC tree
  
  // Branches for HLT tree
  TBranch *caloJetFilterBranch80;         // Branch for calo jet 80 filter bit
  TBranch *caloJetFilterBranch100;        // Branch for calo jet 100 filter bit
  
  // Leaves for the HLT tree
  Int_t caloJetFilterBit80;          // Filter bit for calorimeter jets 80
  Int_t caloJetFilterBit100;         // Filter bit for calorimeter jets 100
  
  // Branches for skim tree
  TBranch *primaryVertexBranch;             // Branch for primary vertex filter bit
  TBranch *beamScrapingBranch;              // Branch for beam scraping filter bit
  TBranch *collisionEventSelectionBranch;   // Branch for collision event selection filter bit
  TBranch *collisionEventSelectionBranchV2; // Branch for collision event selection filter bit v2
  TBranch *hBHENoiseBranch;                 // Branch for HB/HE noise filter bit
  TBranch *hfCoincidenceBranch3Th3;         // Branch for energy recorded in at least 3 HF calorimeter towers
  TBranch *hfCoincidenceBranch2Th4;         // Branch for energy recorded in at least 2 HF calorimeter towers
  TBranch *clusterCompatibilityBranch;      // Branch for cluster compatibility
  
  // Leaves for the skim tree
  Int_t primaryVertexFilterBit;             // Filter bit for primary vertex
  Int_t beamScrapingFilterBit;              // Filter bit for beam scraping
  Int_t collisionEventSelectionFilterBit;   // Filter bit for collision event selection
  Int_t collisionEventSelectionFilterBitV2; // Filter bit for collision event selection v2
  Int_t hBHENoiseFilterBit;                 // Filter bit for HB/HE noise
  Int_t hfCoincidenceFilterBit3Th3;         // Filter bit for energy recorded in at least 3 HF calorimeter towers
  Int_t hfCoincidenceFilterBit2Th4;         // Filter bit for energy recorded in at least 2 HF calorimeter towers
  Int_t clusterCompatibilityFilterBit;      // Filter bit for cluster compatibility
  
  // Branches for jet tree
  TBranch *nJetsBranch[nJetTrees];         // Branch for number of jets in an event
  TBranch *nGenJetsBranch[nJetTrees];      // Branch for the number of generator level jets in an event
  TBranch *jetPtBranch[nJetTrees];         // Branch for jet pT
  TBranch *jetPhiBranch[nJetTrees];        // Branch for jet phi
  TBranch *jetPhiBranchWTA[nJetTrees];     // Branch for jet phi with WTA axis
  TBranch *jetEtaBranch[nJetTrees];        // Branch for jet eta
  TBranch *jetEtaBranchWTA[nJetTrees];     // Branch for jet eta with WTA axis
  TBranch *jetRawPtBranch[nJetTrees];      // Branch for raw jet pT
  TBranch *jetMaxTrackPtBranch[nJetTrees]; // Maximum pT for a track inside a jet
  TBranch *jetRefPtBranch[nJetTrees];      // Branch for reference generator level pT for a reconstructed jet
  TBranch *jetRefFlavorBranch[nJetTrees];  // Branch for flavor for the parton initiating the jet
  TBranch *genJetPtBranch[nJetTrees];      // Branch for the generator level jet pT
  TBranch *genJetEtaBranch[nJetTrees];     // Branch for the generetor level jet eta
  TBranch *genJetEtaBranchWTA[nJetTrees];  // Branch for the generetor level jet eta with WTA axis
  TBranch *genJetPhiBranch[nJetTrees];     // Branch for the generator level jet phi
  TBranch *genJetPhiBranchWTA[nJetTrees];  // Branch for the generator level jet phi with WTA axis
  
  // Leaves for jet tree
  Int_t nJets[nJetTrees];                                 // number of jets in an event
  Int_t nGenJets[nJetTrees];                              // number of generator level jets in an event
  Float_t jetPtArray[nJetTrees][nMaxJet] = {{0}};         // pT:s of all the jets in an event
  Float_t jetPhiArray[nJetTrees][nMaxJet] = {{0}};        // phis of all the jets in an event
  Float_t jetPhiArrayWTA[nJetTrees][nMaxJet] = {{0}};     // phis of all the jets in an event  with WTA axis
  Float_t jetEtaArray[nJetTrees][nMaxJet] = {{0}};        // etas of all the jets in an event
  Float_t jetEtaArrayWTA[nJetTrees][nMaxJet] = {{0}};     // etas of all the jets in an event  with WTA axis
  Float_t jetRawPtArray[nJetTrees][nMaxJet] = {{0}};      // raw jet pT for all the jets in an event
  Float_t jetMaxTrackPtArray[nJetTrees][nMaxJet] = {{0}}; // maximum track pT inside a jet for all the jets in an event
  Float_t jetRefPtArray[nJetTrees][nMaxJet] = {{0}};      // reference generator level pT for a reconstructed jet
  Int_t jetRefFlavorArray[nJetTrees][nMaxJet] = {{0}};    // flavor for initiating parton for the reference gen jet
  Float_t genJetPtArray[nJetTrees][nMaxJet] = {{0}};      // pT:s of all the generator level jets in an event
  Float_t genJetPhiArray[nJetTrees][nMaxJet] = {{0}};     // phis of all the generator level jets in an event
  Float_t genJetPhiArrayWTA[nJetTrees][nMaxJet] = {{0}};  // phis of all the generator level jets in an event with WTA axis
  Float_t genJetEtaArray[nJetTrees][nMaxJet] = {{0}};     // etas of all the generator level jets in an event
  Float_t genJetEtaArrayWTA[nJetTrees][nMaxJet] = {{0}};  // etas of all the generator level jets in an event with WTA axis
  
  // Branches for track tree
  TBranch *nTracksBranch;                    // Branch for number of tracks
  TBranch *trackPtBranch;                    // Branch for track pT
  TBranch *trackPtErrorBranch;               // Branch for track pT error
  TBranch *trackPhiBranch;                   // Branch for track phi
  TBranch *trackEtaBranch;                   // Branch for track eta
  TBranch *trackHighPurityBranch;            // Branch for high purity of the track
  TBranch *trackVertexDistanceZBranch;       // Branch for track distance from primary vertex in z-direction
  TBranch *trackVertexDistanceZErrorBranch;  // Branch for error for track distance from primary vertex in z-direction
  TBranch *trackVertexDistanceXYBranch;      // Branch for track distance from primary vertex in xy-direction
  TBranch *trackVertexDistanceXYErrorBranch; // Branch for error for track distance from primary vertex in xy-direction
  TBranch *trackChi2Branch;                  // Branch for track chi2 value from reconstruction fit
  TBranch *nTrackDegreesOfFreedomBranch;     // Branch for number of degrees of freedom in reconstruction fit
  TBranch *nHitsTrackerLayerBranch;          // Branch for number of hits in tracker layers
  TBranch *nHitsTrackBranch;                 // Branch for number of hits for the track
  TBranch *trackEnergyEcalBranch;            // Branch for track energy in ECal
  TBranch *trackEnergyHcalBranch;            // Branch for track energy in HCal
  TBranch *trackAlgorithmBranch;             // Branch for track algorithm
  TBranch *trackOriginalAlgorithmBranch;     // Branch for track original algorithm
  TBranch *trackMVABranch;                   // Branch for track MVA
  TBranch *trackChargeBranch;                   // Branch for track MVA
  
  // Leaves for the track tree
  Int_t nTracks;                                            // Number of tracks
  Float_t trackPtArray[nMaxTrack] = {0};                    // Array for track pT:s
  Float_t trackPtErrorArray[nMaxTrack] = {0};               // Array for track pT errors
  Float_t trackPhiArray[nMaxTrack] = {0};                   // Array for track phis
  Float_t trackEtaArray[nMaxTrack] = {0};                   // Array for track etas
  Bool_t trackHighPurityArray[nMaxTrack] = {0};             // Array for the high purity of tracks
  Float_t trackVertexDistanceZArray[nMaxTrack] = {0};       // Array for track distance from primary vertex in z-direction
  Float_t trackVertexDistanceZErrorArray[nMaxTrack] = {0};  // Array for error for track distance from primary vertex in z-direction
  Float_t trackVertexDistanceXYArray[nMaxTrack] = {0};      // Array for track distance from primary vertex in xy-direction
  Float_t trackVertexDistanceXYErrorArray[nMaxTrack] = {0}; // Array for error for track distance from primary vertex in xy-direction
  Float_t trackChi2Array[nMaxTrack] = {0};                  // Array for track chi2 value from reconstruction fit
  UChar_t nTrackDegreesOfFreedomArray[nMaxTrack] = {0};     // Array for number of degrees of freedom in reconstruction fit
  UChar_t nHitsTrackerLayerArray[nMaxTrack] = {0};          // Array for number of hits in tracker layers
  UChar_t nHitsTrackArray[nMaxTrack] = {0};                 // Array for number of hits for the track
  Float_t trackEnergyEcalArray[nMaxTrack] = {0};            // Array for track energy in ECal
  Float_t trackEnergyHcalArray[nMaxTrack] = {0};            // Array for track energy in HCal
  UChar_t trackAlgorithmArray[nMaxTrack] = {0};             // Array for track algorithm
  UChar_t trackOriginalAlgorithmArray[nMaxTrack] = {0};     // Array for track original algorithm
  Float_t trackMVAArray[nMaxTrack] = {0};                   // Array for track MVA
  Int_t trackChargeArray[nMaxTrack] = {0};                  // Array for track charge
  
  // Branches for generator level track tree
  TBranch *genTrackPtBranch;         // Branch for generator level track pT:s
  TBranch *genTrackPhiBranch;        // Branch for generator level track phis
  TBranch *genTrackEtaBranch;        // Branch for generator level track etas
  TBranch *genTrackChargeBranch;     // Branch for generator level track charges
  TBranch *genTrackSubeventBranch;   // Branch for generator level track subevent indices (0 = PYTHIA, (>0) = HYDJET)
  
  // Leaves for generator level track tree
  vector<float> *genTrackPtArray;       // Array for generator level track pT:s
  vector<float> *genTrackPhiArray;      // Array for generator level track phis
  vector<float> *genTrackEtaArray;      // Array for generator level track etas
  vector<int> *genTrackChargeArray;     // Array for generator level track charges
  vector<int> *genTrackSubeventArray;   // Array for generator level track subevent indices (0 = PYTHIA, (>0) = HYDJET)
  
  // Branches for particle flow candidate ID tree
  TBranch *particleFlowCandidateIdBranch;    // Branch for particle flow candidate ID
  TBranch *particleFlowCandidatePtBranch;    // Branch for particle flow candidate pT
  TBranch *particleFlowCandidatePhiBranch;   // Branch for particle flow candidate phi
  TBranch *particleFlowCandidateEtaBranch;   // Branch for particle flow candidate eta
  
  // Leaves for particle flow candidate tree
  vector<int> *particleFlowCandidateIdVector;       // Vector for particle flow candidate ID:s
  vector<float> *particleFlowCandidatePtVector;     // Vector for particle flow candidate pT:s
  vector<float> *particleFlowCandidatePhiVector;    // Vector for particle flow candidate phis
  vector<float> *particleFlowCandidateEtaVector;    // Vector for particle flow candidate etas
  
  // Add all the files to the chain
  for(std::vector<TString>::iterator listIterator = fileNameVector.begin(); listIterator != fileNameVector.end(); listIterator++){
    
    cout << "Adding file " << *listIterator << " to the chains" << endl;
    
    heavyIonTree->Add(*listIterator);
    hltTree->Add(*listIterator);
    skimTree->Add(*listIterator);
    trackTree->Add(*listIterator);
    genTrackTree->Add(*listIterator);
    particleFlowCandidateTree->Add(*listIterator);
    
    for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
      jetTree[iJetType]->Add(*listIterator);
    }
    
  }

  
  // ========================================== //
  // Read all the branches from the input trees //
  // ========================================== //
  
  // Connect the branches of the heavy ion tree
  heavyIonTree->SetBranchStatus("*",0);
  heavyIonTree->SetBranchStatus("vz",1);
  heavyIonTree->SetBranchAddress("vz",&vertexZ,&hiVzBranch);
  heavyIonTree->SetBranchStatus("hiBin",1);
  heavyIonTree->SetBranchAddress("hiBin",&hiBin,&hiBinBranch);
  
  // ptHat and event weight only for MC
  if(isMC){
    heavyIonTree->SetBranchStatus("pthat",1);
    heavyIonTree->SetBranchAddress("pthat",&ptHat,&ptHatBranch);
    heavyIonTree->SetBranchStatus("weight",1);
    heavyIonTree->SetBranchAddress("weight",&eventWeight,&eventWeightBranch);
  }
  
  // Connect the branches to the HLT tree
  hltTree->SetBranchStatus("*",0);
  
  // PbPb syntax
  hltTree->SetBranchStatus("HLT_HIPuAK4CaloJet80Eta5p1_v1",1); // 2018 syntax
  hltTree->SetBranchAddress("HLT_HIPuAK4CaloJet80Eta5p1_v1",&caloJetFilterBit80,&caloJetFilterBranch80); // 2018 syntax
  hltTree->SetBranchStatus("HLT_HIPuAK4CaloJet100Eta5p1_v1",1); // 2018 syntax
  hltTree->SetBranchAddress("HLT_HIPuAK4CaloJet100Eta5p1_v1",&caloJetFilterBit100,&caloJetFilterBranch100); // 2018 syntax
  
  // pp syntax
  // hltTree->SetBranchStatus("HLT_HIAK4CaloJet80_v1",1); // 2017 syntax
  // hltTree->SetBranchAddress("HLT_HIAK4CaloJet80_v1",&fCaloJetFilterBit,&fCaloJetFilterBranch);
  
  // Connect the branches to the skim tree
  skimTree->SetBranchStatus("*",0);

  skimTree->SetBranchStatus("pprimaryVertexFilter",1);
  skimTree->SetBranchAddress("pprimaryVertexFilter",&primaryVertexFilterBit,&primaryVertexBranch);
  skimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
  skimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&hBHENoiseFilterBit,&hBHENoiseBranch);
  
  skimTree->SetBranchStatus("collisionEventSelectionAOD",1);
  skimTree->SetBranchAddress("collisionEventSelectionAOD", &collisionEventSelectionFilterBit, &collisionEventSelectionBranch);
  skimTree->SetBranchStatus("collisionEventSelectionAODv2",1);
  skimTree->SetBranchAddress("collisionEventSelectionAODv2", &collisionEventSelectionFilterBitV2, &collisionEventSelectionBranchV2);
  skimTree->SetBranchStatus("phfCoincFilter3Th3",1);
  skimTree->SetBranchAddress("phfCoincFilter3Th3", &hfCoincidenceFilterBit3Th3, &hfCoincidenceBranch3Th3);
  skimTree->SetBranchStatus("phfCoincFilter2Th4",1);
  skimTree->SetBranchAddress("phfCoincFilter2Th4", &hfCoincidenceFilterBit2Th4, &hfCoincidenceBranch2Th4);
  
  skimTree->SetBranchStatus("pclusterCompatibilityFilter",1);
  skimTree->SetBranchAddress("pclusterCompatibilityFilter",&clusterCompatibilityFilterBit,&clusterCompatibilityBranch);
  
  skimTree->SetBranchStatus("pBeamScrapingFilter",1);
  skimTree->SetBranchAddress("pBeamScrapingFilter",&beamScrapingFilterBit,&beamScrapingBranch);
  
  // pp syntax
  // skimTree->SetBranchStatus("pPAprimaryVertexFilter",1);
  // skimTree->SetBranchAddress("pPAprimaryVertexFilter",&primaryVertexFilterBit,&primaryVertexBranch);
  
  // Same branch names for all jet collections
  for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
    
    // Connect the branches to the jet tree
    jetTree[iJetType]->SetBranchStatus("*",0);
    jetTree[iJetType]->SetBranchStatus("jtpt",1);
    jetTree[iJetType]->SetBranchAddress("jtpt",&jetPtArray[iJetType],&jetPtBranch[iJetType]);
    
    // Jet eta with E-scheme and WTA axes
    jetTree[iJetType]->SetBranchStatus("jtphi",1);
    jetTree[iJetType]->SetBranchAddress("jtphi",&jetPhiArray[iJetType],&jetPhiBranch[iJetType]);
    jetTree[iJetType]->SetBranchStatus("WTAphi",1);
    jetTree[iJetType]->SetBranchAddress("WTAphi",&jetPhiArrayWTA[iJetType],&jetPhiBranchWTA[iJetType]);
    
    // Jet phi with E-scheme and WTA axes
    jetTree[iJetType]->SetBranchStatus("jteta",1);
    jetTree[iJetType]->SetBranchAddress("jteta",&jetEtaArray[iJetType],&jetEtaBranch[iJetType]);
    jetTree[iJetType]->SetBranchStatus("WTAeta",1);
    jetTree[iJetType]->SetBranchAddress("WTAeta",&jetEtaArrayWTA[iJetType],&jetEtaBranchWTA[iJetType]);
    
    jetTree[iJetType]->SetBranchStatus("nref",1);
    jetTree[iJetType]->SetBranchAddress("nref",&nJets[iJetType],&nJetsBranch[iJetType]);
    jetTree[iJetType]->SetBranchStatus("rawpt",1);
    jetTree[iJetType]->SetBranchAddress("rawpt",&jetRawPtArray[iJetType],&jetRawPtBranch[iJetType]);
    jetTree[iJetType]->SetBranchStatus("trackMax",1);
    jetTree[iJetType]->SetBranchAddress("trackMax",&jetMaxTrackPtArray[iJetType],&jetMaxTrackPtBranch[iJetType]);
    
    // If we are looking at Monte Carlo, connect the reference pT and parton arrays
    if(isMC){
      jetTree[iJetType]->SetBranchStatus("refpt",1);
      jetTree[iJetType]->SetBranchAddress("refpt",&jetRefPtArray[iJetType],&jetRefPtBranch[iJetType]);
      jetTree[iJetType]->SetBranchStatus("refparton_flavor",1);
      jetTree[iJetType]->SetBranchAddress("refparton_flavor",&jetRefFlavorArray[iJetType],&jetRefFlavorBranch[iJetType]);
      jetTree[iJetType]->SetBranchStatus("genpt",1);
      jetTree[iJetType]->SetBranchAddress("genpt",&genJetPtArray[iJetType],&genJetPtBranch[iJetType]);
      
      // Gen jet phi for e-scheme and WTA axes
      jetTree[iJetType]->SetBranchStatus("genphi",1);
      jetTree[iJetType]->SetBranchAddress("genphi",&genJetPhiArray[iJetType],&genJetPhiBranch[iJetType]);
      jetTree[iJetType]->SetBranchStatus("WTAgenphi",1);
      jetTree[iJetType]->SetBranchAddress("WTAgenphi",&genJetPhiArrayWTA[iJetType],&genJetPhiBranchWTA[iJetType]);
      
      // Gen jet eta for e-scheme and WTA axes
      jetTree[iJetType]->SetBranchStatus("geneta",1);
      jetTree[iJetType]->SetBranchAddress("geneta",&genJetEtaArray[iJetType],&genJetEtaBranch[iJetType]);
      jetTree[iJetType]->SetBranchStatus("WTAgeneta",1);
      jetTree[iJetType]->SetBranchAddress("WTAgeneta",&genJetEtaArrayWTA[iJetType],&genJetEtaBranchWTA[iJetType]);
      
      jetTree[iJetType]->SetBranchStatus("ngen",1);
      jetTree[iJetType]->SetBranchAddress("ngen",&nGenJets[iJetType],&nGenJetsBranch[iJetType]);
    }
    
  } // Loop over different jet collections
  
  // Connect the branches to the track tree
  trackTree->SetBranchStatus("*",0);
  
  trackTree->SetBranchStatus("trkPt",1);
  trackTree->SetBranchAddress("trkPt",&trackPtArray,&trackPtBranch);
  trackTree->SetBranchStatus("trkPtError",1);
  trackTree->SetBranchAddress("trkPtError",&trackPtErrorArray,&trackPtErrorBranch);
  trackTree->SetBranchStatus("trkPhi",1);
  trackTree->SetBranchAddress("trkPhi",&trackPhiArray,&trackPhiBranch);
  trackTree->SetBranchStatus("trkEta",1);
  trackTree->SetBranchAddress("trkEta",&trackEtaArray,&trackEtaBranch);
  trackTree->SetBranchStatus("nTrk",1);
  trackTree->SetBranchAddress("nTrk",&nTracks,&nTracksBranch);
  trackTree->SetBranchStatus("highPurity",1);
  trackTree->SetBranchAddress("highPurity",&trackHighPurityArray,&trackHighPurityBranch);
  trackTree->SetBranchStatus("trkDz1",1);
  trackTree->SetBranchAddress("trkDz1",&trackVertexDistanceZArray,&trackVertexDistanceZBranch);
  trackTree->SetBranchStatus("trkDzError1",1);
  trackTree->SetBranchAddress("trkDzError1",&trackVertexDistanceZErrorArray,&trackVertexDistanceZErrorBranch);
  trackTree->SetBranchStatus("trkDxy1",1);
  trackTree->SetBranchAddress("trkDxy1",&trackVertexDistanceXYArray,&trackVertexDistanceXYBranch);
  trackTree->SetBranchStatus("trkDxyError1",1);
  trackTree->SetBranchAddress("trkDxyError1",&trackVertexDistanceXYErrorArray,&trackVertexDistanceXYErrorBranch);
  trackTree->SetBranchStatus("trkChi2",1);
  trackTree->SetBranchAddress("trkChi2",&trackChi2Array,&trackChi2Branch);
  trackTree->SetBranchStatus("trkNdof",1);
  trackTree->SetBranchAddress("trkNdof",&nTrackDegreesOfFreedomArray,&nTrackDegreesOfFreedomBranch);
  trackTree->SetBranchStatus("trkNlayer",1);
  trackTree->SetBranchAddress("trkNlayer",&nHitsTrackerLayerArray,&nHitsTrackerLayerBranch);
  trackTree->SetBranchStatus("trkNHit",1);
  trackTree->SetBranchAddress("trkNHit",&nHitsTrackArray,&nHitsTrackBranch);
  trackTree->SetBranchStatus("pfEcal",1);
  trackTree->SetBranchAddress("pfEcal",&trackEnergyEcalArray,&trackEnergyEcalBranch);
  trackTree->SetBranchStatus("pfHcal",1);
  trackTree->SetBranchAddress("pfHcal",&trackEnergyHcalArray,&trackEnergyHcalBranch);
  
  // Additional information needed for 2018 track cuts
  trackTree->SetBranchStatus("trkAlgo",1);
  trackTree->SetBranchAddress("trkAlgo",&trackAlgorithmArray,&trackAlgorithmBranch);
  trackTree->SetBranchStatus("trkOriginalAlgo",1);
  trackTree->SetBranchAddress("trkOriginalAlgo",&trackOriginalAlgorithmArray,&trackOriginalAlgorithmBranch);
  trackTree->SetBranchStatus("trkMVA",1);
  trackTree->SetBranchAddress("trkMVA",&trackMVAArray,&trackMVABranch);
  
  // Additional saved branches
  trackTree->SetBranchStatus("trkCharge",1);
  trackTree->SetBranchAddress("trkCharge",&trackChargeArray,&trackChargeBranch);
  
  // Generator level tracks only in Monte Carlo
  if(isMC){
    
    // Connect the branches to generator level track tree
    genTrackTree->SetBranchStatus("*",0);
    
    genTrackTree->SetBranchStatus("pt",1);
    genTrackTree->SetBranchAddress("pt",&genTrackPtArray,&genTrackPtBranch);
    genTrackTree->SetBranchStatus("phi",1);
    genTrackTree->SetBranchAddress("phi",&genTrackPhiArray,&genTrackPhiBranch);
    genTrackTree->SetBranchStatus("eta",1);
    genTrackTree->SetBranchAddress("eta",&genTrackEtaArray,&genTrackEtaBranch);
    genTrackTree->SetBranchStatus("chg",1);
    genTrackTree->SetBranchAddress("chg",&genTrackChargeArray,&genTrackChargeBranch);
    genTrackTree->SetBranchStatus("sube",1);
    genTrackTree->SetBranchAddress("sube",&genTrackSubeventArray,&genTrackSubeventBranch);
    
  }
  
  // Connect the branches to the particle flow candidate tree
  particleFlowCandidateTree->SetBranchStatus("*",0);
  
  particleFlowCandidateTree->SetBranchStatus("pfId",1);
  particleFlowCandidateTree->SetBranchAddress("pfId",&particleFlowCandidateIdVector,&particleFlowCandidateIdBranch);
  particleFlowCandidateTree->SetBranchStatus("pfPt",1);
  particleFlowCandidateTree->SetBranchAddress("pfPt",&particleFlowCandidatePtVector,&particleFlowCandidatePtBranch);
  particleFlowCandidateTree->SetBranchStatus("pfPhi",1);
  particleFlowCandidateTree->SetBranchAddress("pfPhi",&particleFlowCandidatePhiVector,&particleFlowCandidatePhiBranch);
  particleFlowCandidateTree->SetBranchStatus("pfEta",1);
  particleFlowCandidateTree->SetBranchAddress("pfEta",&particleFlowCandidateEtaVector,&particleFlowCandidateEtaBranch);
  
  
  // ========================================== //
  //           Define output trees
  // ========================================== //
  
  // Copy the heavy ion tree to the output
  TTree *heavyIonTreeOutput = new TTree("HiTree","");
  
  // Connect the branches of the heavy ion tree
  heavyIonTreeOutput->Branch("vz",&vertexZ,"vz/F");
  heavyIonTreeOutput->Branch("hiBin",&hiBin,"hiBin/I");
  
  // ptHat and event weight only for MC
  if(isMC){
    heavyIonTreeOutput->Branch("pthat",&ptHat,"pthat/F");
    heavyIonTreeOutput->Branch("weight",&eventWeight,"weight/F");
  }
  
  // Copy the HLT tree to the output
  TTree *hltTreeOutput = new TTree("HltTree","");
  
  // Connect the branches of the HLT tree
  hltTreeOutput->Branch("HLT_HIPuAK4CaloJet80Eta5p1_v1",&caloJetFilterBit80,"HLT_HIPuAK4CaloJet80Eta5p1_v1/I");
  hltTreeOutput->Branch("HLT_HIPuAK4CaloJet100Eta5p1_v1",&caloJetFilterBit100,"HLT_HIPuAK4CaloJet100Eta5p1_v1/I");
  
  // Copy the skim tree to the output
  TTree *skimTreeOutput = new TTree("HltTree","");
  
  skimTreeOutput->Branch("pprimaryVertexFilter",&primaryVertexFilterBit,"pprimaryVertexFilter/I");
  skimTreeOutput->Branch("HBHENoiseFilterResultRun2Loose",&hBHENoiseFilterBit,"HBHENoiseFilterResultRun2Loose/I");
  skimTreeOutput->Branch("collisionEventSelectionAOD", &collisionEventSelectionFilterBit, "collisionEventSelectionAOD/I");
  skimTreeOutput->Branch("collisionEventSelectionAODv2", &collisionEventSelectionFilterBitV2, "collisionEventSelectionAODv2/I");
  skimTreeOutput->Branch("phfCoincFilter3Th3", &hfCoincidenceFilterBit3Th3, "phfCoincFilter3Th3/I");
  skimTreeOutput->Branch("phfCoincFilter2Th4", &hfCoincidenceFilterBit2Th4, "phfCoincFilter2Th4/I");
  skimTreeOutput->Branch("pclusterCompatibilityFilter",&clusterCompatibilityFilterBit,"pclusterCompatibilityFilter/I");
  skimTreeOutput->Branch("pBeamScrapingFilter",&beamScrapingFilterBit,"pBeamScrapingFilter/I");
       
  // Copy the jet trees to the output
  TTree *jetTreeOutput[nJetTrees];
  
  for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
    
    jetTreeOutput[iJetType] = new TTree("t","");
    
    jetTreeOutput[iJetType]->Branch("nref",&nJets[iJetType],"nref/I");
    jetTreeOutput[iJetType]->Branch("jtpt",&jetPtArray[iJetType],"jtpt[nref]/F");
    
    // Jet eta with E-scheme and WTA axes
    jetTreeOutput[iJetType]->Branch("jtphi",&jetPhiArray[iJetType],"jtphi[nref]/F");
    jetTreeOutput[iJetType]->Branch("WTAphi",&jetPhiArrayWTA[iJetType],"WTAphi[nref]/F");
    
    // Jet phi with E-scheme and WTA axes
    jetTreeOutput[iJetType]->Branch("jteta",&jetEtaArray[iJetType],"jteta[nref]/F");
    jetTreeOutput[iJetType]->Branch("WTAeta",&jetEtaArrayWTA[iJetType],"WTAeta[nref]/F");
    
    jetTreeOutput[iJetType]->Branch("rawpt",&jetRawPtArray[iJetType],"rawpt[nref]/F");
    jetTreeOutput[iJetType]->Branch("trackMax",&jetMaxTrackPtArray[iJetType],"trackMax[nref]/F");
    
    // If we are looking at Monte Carlo, connect the reference pT and parton arrays
    if(isMC){
      jetTreeOutput[iJetType]->Branch("refpt",&jetRefPtArray[iJetType],"refpt[nref]/F");
      jetTreeOutput[iJetType]->Branch("refparton_flavor",&jetRefFlavorArray[iJetType],"refparton_flavor[nref]/I");
      
      jetTreeOutput[iJetType]->Branch("ngen",&nGenJets[iJetType],"ngen/I");
      
      jetTreeOutput[iJetType]->Branch("genpt",&genJetPtArray[iJetType],"genpt[ngen]/F");
      
      // Gen jet phi for e-scheme and WTA axes
      jetTreeOutput[iJetType]->Branch("genphi",&genJetPhiArray[iJetType],"genphi[ngen]/F");
      jetTreeOutput[iJetType]->Branch("WTAgenphi",&genJetPhiArrayWTA[iJetType],"WTAgenphi[ngen]/F");
      
      // Gen jet eta for e-scheme and WTA axes
      jetTreeOutput[iJetType]->Branch("geneta",&genJetEtaArray[iJetType],"geneta[ngen]/F");
      jetTreeOutput[iJetType]->Branch("WTAgeneta",&genJetEtaArrayWTA[iJetType],"WTAgeneta[ngen]/F");
      
    } // Branches only for MC

  } // Jet type loop
  
  // Copy the track trees to the output
  TTree *trackTreeOutput = new TTree("trackTree","");
  
  Int_t nTracksOutput;                                       // Number of tracks
  Float_t trackPtOutput[nMaxTrack] = {0};                    // Array for track pT:s
  Float_t trackPtErrorOutput[nMaxTrack] = {0};               // Array for track pT errors
  Float_t trackPhiOutput[nMaxTrack] = {0};                   // Array for track phis
  Float_t trackEtaOutput[nMaxTrack] = {0};                   // Array for track etas
  Bool_t trackHighPurityOutput[nMaxTrack] = {0};             // Array for the high purity of tracks
  Float_t trackVertexDistanceZOutput[nMaxTrack] = {0};       // Array for track distance from primary vertex in z-direction
  Float_t trackVertexDistanceZErrorOutput[nMaxTrack] = {0};  // Array for error for track distance from primary vertex in z-direction
  Float_t trackVertexDistanceXYOutput[nMaxTrack] = {0};      // Array for track distance from primary vertex in xy-direction
  Float_t trackVertexDistanceXYErrorOutput[nMaxTrack] = {0}; // Array for error for track distance from primary vertex in xy-direction
  Float_t trackChi2Output[nMaxTrack] = {0};                  // Array for track chi2 value from reconstruction fit
  UChar_t nTrackDegreesOfFreedomOutput[nMaxTrack] = {0};     // Array for number of degrees of freedom in reconstruction fit
  UChar_t nHitsTrackerLayerOutput[nMaxTrack] = {0};          // Array for number of hits in tracker layers
  UChar_t nHitsTrackOutput[nMaxTrack] = {0};                 // Array for number of hits for the track
  Float_t trackEnergyEcalOutput[nMaxTrack] = {0};            // Array for track energy in ECal
  Float_t trackEnergyHcalOutput[nMaxTrack] = {0};            // Array for track energy in HCal
  UChar_t trackAlgorithmOutput[nMaxTrack] = {0};             // Array for track algorithm
  UChar_t trackOriginalAlgorithmOutput[nMaxTrack] = {0};     // Array for track original algorithm
  Float_t trackMVAOutput[nMaxTrack] = {0};                   // Array for track MVA
  Int_t trackChargeOutput[nMaxTrack] = {0};                  // Array for track charge
  
  trackTreeOutput->Branch("nTrk",&nTracksOutput,"nTrk/I");
  trackTreeOutput->Branch("trkPt",&trackPtOutput,"trkPt[nTrk]/F");
  trackTreeOutput->Branch("trkPtError",&trackPtErrorOutput,"trkPtError[nTrk]/F");
  trackTreeOutput->Branch("trkPhi",&trackPhiOutput,"trkPhi[nTrk]/F");
  trackTreeOutput->Branch("trkEta",&trackEtaOutput,"trkEta[nTrk]/F");
  
  trackTreeOutput->Branch("highPurity",&trackHighPurityOutput,"highPurity[nTrk]/O");
  trackTreeOutput->Branch("trkDz1",&trackVertexDistanceZOutput,"trkDz1[nTrk]/F");
  trackTreeOutput->Branch("trkDzError1",&trackVertexDistanceZErrorOutput,"trkDzError1[nTrk]/F");
  trackTreeOutput->Branch("trkDxy1",&trackVertexDistanceXYOutput,"trkDxy1[nTrk]/F");
  trackTreeOutput->Branch("trkDxyError1",&trackVertexDistanceXYErrorOutput,"trkDxyError1[nTrk]/F");
  trackTreeOutput->Branch("trkChi2",&trackChi2Output,"trkChi2[nTrk]/F");
  trackTreeOutput->Branch("trkNdof",&nTrackDegreesOfFreedomOutput,"trkNdof[nTrk]/b");
  trackTreeOutput->Branch("trkNlayer",&nHitsTrackerLayerOutput,"trkNlayer[nTrk]/b");
  trackTreeOutput->Branch("trkNHit",&nHitsTrackOutput,"trkNHit[nTrk]/b");
  trackTreeOutput->Branch("pfEcal",&trackEnergyEcalOutput,"pfEcal[nTrk]/F");
  trackTreeOutput->Branch("pfHcal",&trackEnergyHcalOutput,"pfHcal[nTrk]/F");
  
  // Additional information needed for 2018 track cuts
  trackTreeOutput->Branch("trkAlgo",&trackAlgorithmOutput,"trkAlgo[nTrk]/b");
  trackTreeOutput->Branch("trkOriginalAlgo",&trackOriginalAlgorithmOutput,"trkOriginalAlgo[nTrk]/b");
  trackTreeOutput->Branch("trkMVA",&trackMVAOutput,"trkMVA[nTrk]/F");
  
  // Additional branches
  trackTreeOutput->Branch("trkCharge",&trackChargeOutput,"trkCharge[nTrk]/I");
  
  // Generator level tracks only in Monte Carlo
  TTree *genTrackTreeOutput = new TTree("hi","");
  
  std::vector<float> *genTrackPtVector = new std::vector<float>(); genTrackPtVector->clear();
  std::vector<float> *genTrackPhiVector = new std::vector<float>(); genTrackPhiVector->clear();
  std::vector<float> *genTrackEtaVector = new std::vector<float>(); genTrackEtaVector->clear();
  std::vector<int> *genTrackChargeVector = new std::vector<int>(); genTrackChargeVector->clear();
  std::vector<int> *genTrackSubeventVector = new std::vector<int>(); genTrackSubeventVector->clear();
  
  // Connect the branches to generator level track tree
  if(isMC){
    genTrackTreeOutput->Branch("pt","vector<float>", &genTrackPtVector);
    genTrackTreeOutput->Branch("phi","vector<float>", &genTrackPhiVector);
    genTrackTreeOutput->Branch("eta","vector<float>", &genTrackEtaVector);
    genTrackTreeOutput->Branch("chg","vector<int>", &genTrackChargeVector);
    genTrackTreeOutput->Branch("sube","vector<int>", &genTrackSubeventVector);
  }
  
  // Copy the particle flow candidate tree to the output
  TTree *particleFlowCandidateTreeOutput = new TTree("pfTree","");
  
  std::vector<int> *particleFlowCandidateIdOutputVector = new std::vector<int>(); particleFlowCandidateIdOutputVector->clear();
  std::vector<float> *particleFlowCandidatePtOutputVector = new std::vector<float>(); particleFlowCandidatePtOutputVector->clear();
  std::vector<float> *particleFlowCandidatePhiOutputVector = new std::vector<float>(); particleFlowCandidatePhiOutputVector->clear();
  std::vector<float> *particleFlowCandidateEtaOutputVector = new std::vector<float>(); particleFlowCandidateEtaOutputVector->clear();
  
  particleFlowCandidateTreeOutput->Branch("pfId","vector<int>",&particleFlowCandidateIdOutputVector);
  particleFlowCandidateTreeOutput->Branch("pfPt","vector<float>",&particleFlowCandidatePtOutputVector);
  particleFlowCandidateTreeOutput->Branch("pfPhi","vector<float>",&particleFlowCandidatePhiOutputVector);
  particleFlowCandidateTreeOutput->Branch("pfEta","vector<float>",&particleFlowCandidateEtaOutputVector);
  
  // ========================================== //
  //          Loop over all events              //
  // ========================================== //
  
  int nEvents = heavyIonTree->GetEntries();
  cout << "There are " << nEvents << " events" << endl;
  
  bool passTrackCuts;
  int iTrackOutput;
  
  for(int iEvent = 0; iEvent < nEvents; iEvent++) {
    
    if( iEvent % 1000 == 0 )  std::cout << "iEvent: " << iEvent <<  " of " << nEvents << std::endl;
    
    // ========================================== //
    //        Read the event to input trees
    // ========================================== //
    
    heavyIonTree->GetEntry(iEvent);
    hltTree->GetEntry(iEvent);
    skimTree->GetEntry(iEvent);
    trackTree->GetEntry(iEvent);
    if(isMC) genTrackTree->GetEntry(iEvent);
    particleFlowCandidateTree->GetEntry(iEvent);
    
    for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
      jetTree[iJetType]->GetEntry(iEvent);
    }
    
    heavyIonTreeOutput->Fill();
    hltTreeOutput->Fill();
    skimTreeOutput->Fill();
    
    
    for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
      jetTreeOutput[iJetType]->Fill();
    }
    
    
    // Reco track loop
    nTracksOutput = nTracks;
    iTrackOutput = 0;
    for(int iTrack = 0; iTrack < nTracks; iTrack++){
      
      passTrackCuts = true;
      
      // Do basic track cuts for the reconstructed tracks
      if(trackHighPurityArray[iTrack] != 1) passTrackCuts = false;
      
      if(trackPtErrorArray[iTrack]/trackPtArray[iTrack] > 0.1) passTrackCuts = false;
      if(TMath::Abs(trackVertexDistanceZArray[iTrack]/trackVertexDistanceZErrorArray[iTrack]) > 3) passTrackCuts = false;
      if(TMath::Abs(trackVertexDistanceXYArray[iTrack]/trackVertexDistanceXYErrorArray[iTrack]) > 3) passTrackCuts = false;
      
      if(TMath::Abs(trackEtaArray[iTrack]) >= 2.4) passTrackCuts = false;  //acceptance of the tracker
      
      if(trackPtArray[iTrack] < 0.7) passTrackCuts = false;   // Minimum track pT
      if(trackPtArray[iTrack] > 300 ) passTrackCuts = false;  // Maximum track pT
      
      if(passTrackCuts){
        trackPtOutput[iTrackOutput] = trackPtArray[iTrack];
        trackPhiOutput[iTrackOutput] = trackPhiArray[iTrack];
        trackEtaOutput[iTrackOutput] = trackEtaArray[iTrack];
        trackHighPurityOutput[iTrackOutput] = trackHighPurityArray[iTrack];
        trackVertexDistanceZOutput[iTrackOutput] = trackVertexDistanceZArray[iTrack];
        trackVertexDistanceZErrorOutput[iTrackOutput] = trackVertexDistanceZErrorArray[iTrack];
        trackVertexDistanceXYOutput[iTrackOutput] = trackVertexDistanceXYArray[iTrack];
        trackVertexDistanceXYErrorOutput[iTrackOutput] = trackVertexDistanceXYErrorArray[iTrack];
        trackChi2Output[iTrackOutput] = trackChi2Output[iTrack];
        nTrackDegreesOfFreedomOutput[iTrackOutput] = nTrackDegreesOfFreedomArray[iTrack];
        nHitsTrackerLayerOutput[iTrackOutput] = nHitsTrackerLayerArray[iTrack];
        nHitsTrackOutput[iTrackOutput] = nHitsTrackArray[iTrack];
        trackEnergyEcalOutput[iTrackOutput] = trackEnergyEcalArray[iTrack];
        trackEnergyHcalOutput[iTrackOutput] = trackEnergyHcalArray[iTrack];
        trackAlgorithmOutput[iTrackOutput] = trackAlgorithmArray[iTrack];
        trackOriginalAlgorithmOutput[iTrackOutput] = trackOriginalAlgorithmArray[iTrack];
        trackMVAOutput[iTrackOutput] = trackMVAArray[iTrack];
        trackChargeOutput[iTrackOutput] = trackChargeArray[iTrack];
        iTrackOutput++;
      } else {
        nTracksOutput--;
      }
      
    }
    
    trackTreeOutput->Fill();
    
    // Gen track loop
    if(isMC){
      for(int iTrack = 0; iTrack < genTrackPtArray->size(); iTrack++){
        
        // Cut away low pT tracks and tracks with eta outside of tracker acceptance
        if(TMath::Abs(genTrackEtaArray->at(iTrack)) >= 2.4) continue; //acceptance of the tracker
        
        if(genTrackPtArray->at(iTrack) < 0.7) continue;   // Minimum track pT
        if(genTrackPtArray->at(iTrack) > 300 ) continue;  // Maximum track pT
        
        // Fill the output vectors with gen particles surviving the cuts
        genTrackPtVector->push_back(genTrackPtArray->at(iTrack));
        genTrackPhiVector->push_back(genTrackPhiArray->at(iTrack));
        genTrackEtaVector->push_back(genTrackEtaArray->at(iTrack));
        genTrackChargeVector->push_back(genTrackChargeArray->at(iTrack));
        genTrackSubeventVector->push_back(genTrackSubeventArray->at(iTrack));
      }
      
      genTrackTreeOutput->Fill();
      
    } // Filling gen tracks for MC
    
    // Particle flow candidate loop
    for(int pfi = 0; pfi < particleFlowCandidateIdVector->size(); pfi++) {
      
      particleFlowCandidateIdOutputVector->push_back(particleFlowCandidateIdVector->at(pfi));
      particleFlowCandidatePtOutputVector->push_back(particleFlowCandidatePtVector->at(pfi));
      particleFlowCandidatePhiOutputVector->push_back(particleFlowCandidatePhiVector->at(pfi));
      particleFlowCandidateEtaOutputVector->push_back(particleFlowCandidateEtaVector->at(pfi));
      
    } // particle flow candidate loop
    
    particleFlowCandidateTreeOutput->Fill();
    
    // Clear the vectors before the next event! Otherwise all the tracks pile up cumulatively
    if(isMC){
      genTrackPtVector->clear();
      genTrackPhiVector->clear();
      genTrackEtaVector->clear();
      genTrackChargeVector->clear();
      genTrackSubeventVector->clear();
    }
    
    particleFlowCandidateIdOutputVector->clear();
    particleFlowCandidatePtOutputVector->clear();
    particleFlowCandidatePhiOutputVector->clear();
    particleFlowCandidateEtaOutputVector->clear();
    
  } // Event loop
  
  // Write the skimmed trees to the output file
  
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  
  gDirectory->mkdir("hiEvtAnalyzer");
  gDirectory->cd("hiEvtAnalyzer");
  
  heavyIonTreeOutput->Write();
  
  gDirectory->cd("../");
  gDirectory->mkdir("hltanalysis");
  gDirectory->cd("hltanalysis");
  
  hltTreeOutput->Write();
  
  gDirectory->cd("../");
  gDirectory->mkdir("skimanalysis");
  gDirectory->cd("skimanalysis");
  
  skimTreeOutput->Write();
  
  gDirectory->cd("../");
  
  const char *jetDirectories[] = {"akPu4CaloJetAnalyzer","akCs4PFJetAnalyzer","akPu4PFJetAnalyzer","akFlowPuCs4PFJetAnalyzer"};
  
  for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
    
    gDirectory->mkdir(jetDirectories[iJetType]);
    gDirectory->cd(jetDirectories[iJetType]);
    
    jetTreeOutput[iJetType]->Write();
    
    gDirectory->cd("../");
    
  } // Loop over jet types
  
  gDirectory->mkdir("ppTrack");
  gDirectory->cd("ppTrack");
  
  trackTreeOutput->Write();
  
  gDirectory->cd("../");
  gDirectory->mkdir("HiGenParticleAna");
  gDirectory->cd("HiGenParticleAna");
  
  genTrackTreeOutput->Write();
  
  gDirectory->cd("../");
  gDirectory->mkdir("pfcandAnalyzer");
  gDirectory->cd("pfcandAnalyzer");

  particleFlowCandidateTreeOutput->Write();
  
  gDirectory->cd("../");
  
  outputFile->Close();
  
  return 0;
  
}
