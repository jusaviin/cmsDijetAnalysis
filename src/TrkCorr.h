#ifndef TRKCORR
#define TRKCORR

#include "TMath.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TrkSettings.h"
#include <iostream>
#include <vector>

class TrkCorr{
  public:
    double getTrkCorr(float pt, float eta, float phi, float hiBin, float rmin=99, float jtpt=0, int correction=0);
    TrkCorr(std::string inputDirectory = "trkCorrections/");
    ~TrkCorr();    

  private:
    int nFiles;
    int nSteps;

    std::vector<std::vector<TH1D*> > eff;
    std::vector<std::vector<TH1D*> > fake;
    std::vector<std::vector<TH2D*> > eff2;
    std::vector<std::vector<TH2D*> > fake2;
    std::vector<TH2D*> secondary;
    std::vector<TH1D*> multiple;

    TrkSettings * s;
};

#endif
