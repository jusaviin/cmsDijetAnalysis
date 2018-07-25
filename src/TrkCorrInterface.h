#ifndef TRKCORRINTERFACE_H
#define TRKCORRINTERFACE_H

class TrkCorrInterface{
  
public:
  
  virtual double getTrkCorr(float pt, float eta, float phi, float hiBin, float rmin = 99, float jtpt = 0, int correction = 0) = 0;
  virtual ~TrkCorrInterface();
  
};

#endif
