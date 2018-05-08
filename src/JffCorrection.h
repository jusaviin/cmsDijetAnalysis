#ifndef JFFCORRECTION_H
#define JFFCORRECTION_H

/*
 *  Class for JFF correction for jet pT. Code from Dhanush Hangal.
 */

#include "TFile.h"
#include "TF1.h"
#include <iostream>

class JffCorrection{
  
public:
  JffCorrection(bool ispp);
  double GetCorrection(bool ispp, int nCScand, int hiBin, double jtpt, double jteta);
  
private:
  static const int nCentBins = 4;
  TF1 *f_param_a0[nCentBins];
  TF1 *f_param_a1[nCentBins];
  TF1 *flat_corr[nCentBins];
  TF1 *f_pol1[nCentBins];
  double* centBins;
  TFile *fin;
  double corrpt;
};

#endif
