// C++ includes
#include <iostream>   // Input/output stream. Needed for cout.
#include <fstream>    // File stream for intup/output to/from files
#include <stdlib.h>   // Standard utility libraries
#include <assert.h>   // Standard c++ debugging tool. Terminates the program if expression given evaluates to 0.
#include <vector>     // C++ vector class
#include <sstream>    // Libraries for checking boolean input
#include <string>     // Libraries for checking boolean input
#include <iomanip>    // Libraries for checking boolean input
#include <algorithm>  // Libraries for checking boolean input
#include <cctype>     // Libraries for checking boolean input

// Includes from Root
//#include <TROOT.h>    // Not sure if needed...
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>

// Own includes
#include "src/DijetAnalyzer.h"
#include "src/ConfigurationCard.h"
#include "src/DijetHistograms.h"

using namespace std;

/*
 * File list reader
 *
 *  Arguments:
 *    std::vector<TString> &fileNameVector = Vector filled with filenames found in the file
 *    TString fileNameFile = Text file containing one analysis file name in each line
 *    int debug = Level of debug messages shown
 *    bool runLocal = True: Local run mode. False: Crab run mode
 */
void ReadFileList(std::vector<TString> &fileNameVector, TString fileNameFile, int debug, bool runLocal)
{
  
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
            currentFileName.Prepend("root://cmsxrootd.fnal.gov//");  // If not running locally, we need to give xrootd path
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

/*
 *  Main program
 *
 *  Command line arguments:
 *  argv[1] = List of files to be analyzed, given in text file
 *  argv[2] = Debugging level. 0 = none, 1 = some, 2 = all.
 */
int main(int argc, char **argv) {
  
  //==== Read arguments =====
  //TROOT root("nanoDST","nanoDST analysis");  // I do not really know what this does.
  if ( argc<5 ) {
    cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    cout<<"+ Usage of the macro: " << endl;
    cout<<"+  "<<argv[0]<<" [fileNameFile] [configurationCard] [outputFileName] <debugLevel> <runLocal>"<<endl;
    cout<<"+  fileNameFile: Text file containing the list of files used in the analysis. For crab analysis a job id should be given here." <<endl;
    cout<<"+  configurationCard: Card file with binning and cut information for the analysis." <<endl;
    cout<<"+  outputFileName: .root file to which the histograms are written." <<endl;
    cout<<"+  debugLevel: Amount of debug messages shown. 0 = none, 1 = some, 2 = all." <<endl;
    cout<<"+  runLocal: True: Search input files from local machine. False (default): Search input files from root." << endl;
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
  const char *cardName = argv[2];
  TString outputFileName = argv[3];
  int debugLevel = 0;
  if(argc >= 5) debugLevel = atoi(argv[4]);
  
  // Read the card
  ConfigurationCard *dijetCard = new ConfigurationCard(cardName);
  if(debugLevel > 0){
    dijetCard->PrintOut();
    cout << endl;
  }
  
  // Read the file names used for the analysis to a vector
  std::vector<TString> fileNameVector;
  fileNameVector.clear();
  ReadFileList(fileNameVector,fileNameFile,debugLevel,runLocal);
  
  // Variable for histograms in the analysis
  DijetHistograms *histograms;
  
  // Run the analysis over the list of files
  DijetAnalyzer *jetAnalysis = new DijetAnalyzer(fileNameVector, dijetCard);
  jetAnalysis->RunAnalysis(debugLevel);
  histograms = jetAnalysis->GetHistograms();
  
  // Write the histograms and card to file
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  histograms->Write();
  dijetCard->WriteCard(outputFile);
  outputFile->Close();
  
}

