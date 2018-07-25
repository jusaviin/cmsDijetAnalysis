/*
 * Header file for mixed event lookout table
 */

#ifndef MIXEDEVENTLOOKOUTTABLE_H
#define MIXEDEVENTLOOKOUTTABLE_H

#include <TString.h>
#include "ForestReader.h"

class MixedEventLookoutTable{
  
public:
  
  MixedEventLookoutTable(int dataType); // Default constructor
  virtual ~MixedEventLookoutTable(); // Destructor
  
  TString GetMixingFileName(TString dataFileName); // Get a file to be used for mixing
  
private:
  
  static const int kMaxFiles = 1600;        // Maximum number of files in arrays
  TString fDataFileArray[kMaxFiles] = {""};        // Array for data file names
  TString fMixedEventFileArray[kMaxFiles] = {""};  // Mixed event files for a given data file
  
  void SetArraysPbPbMC();   // Set the arrays for PbPb MC
  void SetArraysPythia8();  // Set the arrays for Pythia8 forest
  
};

#endif
