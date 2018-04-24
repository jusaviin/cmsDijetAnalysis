#ifndef INPUTSETTINGS
#define INPUTSETTINGS

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

class TrkSettings {
  public:
  //TrkSettings that will be set when a job number is specified
   int job;
   int nSkip;
   double ptMin, ptMax, centPUMin, centPUMax;

  //TrkSettings that should be the same across all processes
   std::string jobName;
   bool reuseSkim;
   bool checkClosure;
   int nPb, vz_window, nMC;
   bool doEff, doFake, doMult, doSecondary;

   bool doPthat, doVtx, doCentPU;
   std::vector<std::string> MCFiles;
   std::string DataFile;
   std::vector<double> pthatBins, pthatCrossSection;
   std::vector<int> multiRecoBins;

   int nPtBinCoarse, ptBinFine, etaBinFine, phiBinFine, centPUBinFine, jetBinFine;
   std::vector<double> ptBinCoarse;
   std::vector<int> nCentPUBinCoarse;
   std::vector<std::vector<double> > centPUBinCoarse, eventSkip;

   int highPurityDef;
   bool doCaloMatch, doTrackTriggerCuts, doOtherCut, doSplit;

   int nStep, fullIterations, terminateStep;
   std::vector<int> stepOrder;
   std::string jetDefinition;
   
   std::string trackTreeName;

   TrkSettings(std::string = "TrkCorrInputFile.txt");
};

#endif
