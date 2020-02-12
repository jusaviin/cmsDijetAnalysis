/*
 * Implementation of the SpilloverFluctuationCleaner class
 */

// Own includes
#include "SpilloverFluctuationCleaner.h"

/*
 * Contructor
 */
SpilloverFluctuationCleaner::SpilloverFluctuationCleaner()
{

}

/*
 * Destructor
 */
SpilloverFluctuationCleaner::~SpilloverFluctuationCleaner(){

}

/*
 * Clean the fluctuations from a given bin
 */
void SpilloverFluctuationCleaner::CleanSpilloverFluctuationDeltaR(TH1* jetShapeHistogram, const int iJetTrack, const int iAsymmetry, const int iCentrality, const int iTrackPt) const{
  
  // Currently only implemented for pt weighted leading jet
  if(iJetTrack != DijetHistogramManager::kPtWeightedTrackLeadingJet) return;
  
  // Cleaning configuration for bin 0.0 < xj < 0.6, C = 0-10 %, 1 < pT < 2 GeV
  if(iAsymmetry == 0 && iCentrality == 0 && iTrackPt == 1){
    jetShapeHistogram->SetBinContent(11,8.2853939*0.52);   // Bin 11: 0.5 < DeltaR < 0.6
    jetShapeHistogram->SetBinContent(12,6.9338037*0.52);   // Bin 12: 0.6 < DeltaR < 0.7
    
  // Cleaning configuration for bin 0.0 < xj < 0.6, C = 0-10 %, 2 < pT < 3 GeV
  } else if(iAsymmetry == 0 && iCentrality == 0 && iTrackPt == 2){
    jetShapeHistogram->SetBinContent(12,2.9801714*0.56);   // Bin 12: 0.6 < DeltaR < 0.7
    
  // Cleaning configuration for bin 0.0 < xj < 0.6, C = 0-10 %, 3 < pT < 4 GeV
  } else if(iAsymmetry == 0 && iCentrality == 0 && iTrackPt == 3){
    jetShapeHistogram->SetBinContent(14,0.40682701*0.4);   // Bin 14: 0.8 < DeltaR < 1.0
  
  // Cleaning configuration for bin 0.0 < xj < 0.6, C = 0-10 %, 4 < pT < 8 GeV
  } else if(iAsymmetry == 0 && iCentrality == 0 && iTrackPt == 4){
    jetShapeHistogram->SetBinContent(8,4.8304080*0.58);    // Bin 8:  0.35 < DeltaR < 0.4
    jetShapeHistogram->SetBinContent(10,1.8283094*0.56);   // Bin 10: 0.45 < DeltaR < 0.5
    jetShapeHistogram->SetBinContent(12,0.83188098*0.56);  // Bin 12:  0.6 < DeltaR < 0.7
    jetShapeHistogram->SetBinContent(13,0.45755407*0.6);   // Bin 13:  0.7 < DeltaR < 0.8
    
  // Cleaning configuration for bin 0.0 < xj < 0.6, C = 10-30 %, 1 < pT < 2 GeV
  } else if(iAsymmetry == 0 && iCentrality == 1 && iTrackPt == 1){
    jetShapeHistogram->SetBinContent(12,5.2977733*0.68);   // Bin 12: 0.6 < DeltaR < 0.7
    jetShapeHistogram->SetBinContent(13,4.1776487*0.68);   // Bin 13: 0.7 < DeltaR < 0.8
    
  // Cleaning configuration for bin 0.0 < xj < 0.6, C = 10-30 %, 2 < pT < 3 GeV
  } else if(iAsymmetry == 0 && iCentrality == 1 && iTrackPt == 2){
    jetShapeHistogram->SetBinContent(9,6.0743246*0.66);    // Bin 9:  0.4 < DeltaR < 0.45
    jetShapeHistogram->SetBinContent(11,3.2986801*0.66);   // Bin 11: 0.5 < DeltaR < 0.6
   
  // Cleaning configuration for bin 0.0 < xj < 0.6, C = 10-30 %, 8 < pT < 12 GeV
  } else if(iAsymmetry == 0 && iCentrality == 1 && iTrackPt == 5){
    jetShapeHistogram->SetBinContent(10,0.61743523*0.88);  // Bin 10: 0.45 < DeltaR < 0.5
   
  // Cleaning configuration for bin 0.0 < xj < 0.6, C = 30-50 %, 2 < pT < 3 GeV
  } else if(iAsymmetry == 0 && iCentrality == 2 && iTrackPt == 2){
    jetShapeHistogram->SetBinContent(10,3.0561156*0.77);   // Bin 10: 0.45 < DeltaR < 0.5
    jetShapeHistogram->SetBinContent(13,0.82120308*0.82);  // Bin 13:  0.7 < DeltaR < 0.8
    
  // Cleaning configuration for bin 0.0 < xj < 0.6, C = 30-50 %, 3 < pT < 4 GeV
  } else if(iAsymmetry == 0 && iCentrality == 2 && iTrackPt == 3){
    jetShapeHistogram->SetBinContent(10,1.3586914*0.74);   // Bin 10: 0.45 < DeltaR < 0.5
    jetShapeHistogram->SetBinContent(11,0.80764494*0.74);  // Bin 11:  0.5 < DeltaR < 0.6
   
  // Cleaning configuration for bin 0.6 < xj < 0.8, C = 0-10 %, 0.7 < pT < 1 GeV
  } else if(iAsymmetry == 1 && iCentrality == 0 && iTrackPt == 0){
    jetShapeHistogram->SetBinContent(6,4.4646498*0.54);    // Bin 6: 0.25 < DeltaR < 0.3
    jetShapeHistogram->SetBinContent(14,1.8661903*0.54);   // Bin 14: 0.8 < DeltaR < 1.0
    
  // Cleaning configuration for bin 0.6 < xj < 0.8, C = 0-10 %, 1 < pT < 2 GeV
  } else if(iAsymmetry == 1 && iCentrality == 0 && iTrackPt == 1){
    jetShapeHistogram->SetBinContent(6,17.647416*0.65);    // Bin 6: 0.25 < DeltaR < 0.3
    
  // Cleaning configuration for bin 0.6 < xj < 0.8, C = 0-10 %, 3 < pT < 4 GeV
  } else if(iAsymmetry == 1 && iCentrality == 0 && iTrackPt == 3){
    jetShapeHistogram->SetBinContent(11,1.4693197*0.82);   // Bin 11: 0.5 < DeltaR < 0.6
    jetShapeHistogram->SetBinContent(13,0.88926778*0.82);  // Bin 13: 0.7 < DeltaR < 0.8
    
  // Cleaning configuration for bin 0.6 < xj < 0.8, C = 0-10 %, 4 < pT < 8 GeV
  } else if(iAsymmetry == 1 && iCentrality == 0 && iTrackPt == 4){
    jetShapeHistogram->SetBinContent(12,0.74733002*0.82);  // Bin 12: 0.6 < DeltaR < 0.7
    jetShapeHistogram->SetBinContent(13,0.42728863*0.86);  // Bin 13: 0.7 < DeltaR < 0.8
    
  // Cleaning configuration for bin 0.6 < xj < 0.8, C = 0-10 %, 8 < pT < 12 GeV
  } else if(iAsymmetry == 1 && iCentrality == 0 && iTrackPt == 5){
    jetShapeHistogram->SetBinContent(11,0.33611123*1);     // Bin 11: 0.5 < DeltaR < 0.6
    
  // Cleaning configuration for bin 0.6 < xj < 0.8, C = 10-30 %, 1 < pT < 2 GeV
  } else if(iAsymmetry == 1 && iCentrality == 1 && iTrackPt == 1){
    jetShapeHistogram->SetBinContent(13,4.3375146*0.86);   // Bin 13: 0.7 < DeltaR < 0.8
    
  // Cleaning configuration for bin 0.6 < xj < 0.8, C = 30-50 %, 0.7 < pT < 1 GeV
  } else if(iAsymmetry == 1 && iCentrality == 2 && iTrackPt == 0){
    jetShapeHistogram->SetBinContent(11,1.4457849*0.89);   // Bin 11: 0.5 < DeltaR < 0.6
    
  // Cleaning configuration for bin 0.8 < xj < 1.0, C = 0-10 %, 0.7 < pT < 1 GeV
  } else if(iAsymmetry == 2 && iCentrality == 0 && iTrackPt == 0){
    jetShapeHistogram->SetBinContent(9,3.2980074*0.72);    // Bin 9: 0.4 < DeltaR < 0.45
    jetShapeHistogram->SetBinContent(12,2.1015206*0.8);    // Bin 12: 0.6 < DeltaR < 0.7
    
  // Cleaning configuration for bin 0.8 < xj < 1.0, C = 0-10 %, 1 < pT < 2 GeV
  } else if(iAsymmetry == 2 && iCentrality == 0 && iTrackPt == 1){
    jetShapeHistogram->SetBinContent(11,9.3967489*0.79);   // Bin 11: 0.5 < DeltaR < 0.6
    jetShapeHistogram->SetBinContent(13,6.1780839*0.86);   // Bin 13: 0.7 < DeltaR < 0.8
    
  // Cleaning configuration for bin 0.8 < xj < 1.0, C = 0-10 %, 2 < pT < 3 GeV
  } else if(iAsymmetry == 2 && iCentrality == 0 && iTrackPt == 2){
    jetShapeHistogram->SetBinContent(14,1.6302581*1);      // Bin 14: 0.8 < DeltaR < 1.0
    
  // Cleaning configuration for bin 0.8 < xj < 1.0, C = 0-10 %, 3 < pT < 4 GeV
  } else if(iAsymmetry == 2 && iCentrality == 0 && iTrackPt == 3){
    jetShapeHistogram->SetBinContent(14,0.48591667*1);     // Bin 14: 0.8 < DeltaR < 1.0
    
  // Cleaning configuration for bin 0.8 < xj < 1.0, C = 0-10 %, 4 < pT < 8 GeV
  } else if(iAsymmetry == 2 && iCentrality == 0 && iTrackPt == 4){
    jetShapeHistogram->SetBinContent(12,0.73699003*0.9);   // Bin 12: 0.6 < DeltaR < 0.7
    jetShapeHistogram->SetBinContent(13,0.57238499*0.94);  // Bin 13: 0.7 < DeltaR < 0.8
    
  // Cleaning configuration for bin 0.8 < xj < 1.0, C = 0-10 %, 8 < pT < 12 GeV
  } else if(iAsymmetry == 2 && iCentrality == 0 && iTrackPt == 5){
    jetShapeHistogram->SetBinContent(7,2.6180083*0.95);    // Bin 7: 0.3 < DeltaR < 0.35
    
  // Cleaning configuration for bin 0.8 < xj < 1.0, C = 10-30 %, 0.7 < pT < 1 GeV
  } else if(iAsymmetry == 2 && iCentrality == 1 && iTrackPt == 0){
    jetShapeHistogram->SetBinContent(12,1.9804747*0.85);   // Bin 12: 0.6 < DeltaR < 0.7
    jetShapeHistogram->SetBinContent(13,1.6187124*0.85);   // Bin 13: 0.7 < DeltaR < 0.8
    jetShapeHistogram->SetBinContent(14,1.2260020*0.85);   // Bin 14: 0.8 < DeltaR < 1.0
    
  // Cleaning configuration for bin 0.8 < xj < 1.0, C = 10-30 %, 3 < pT < 4 GeV
  } else if(iAsymmetry == 2 && iCentrality == 1 && iTrackPt == 3){
    jetShapeHistogram->SetBinContent(11,1.3141761*0.85);   // Bin 11: 0.5 < DeltaR < 0.6
    jetShapeHistogram->SetBinContent(13,0.69913071*0.94);  // Bin 13: 0.7 < DeltaR < 0.8
    
  // Cleaning configuration for bin 0.8 < xj < 1.0, C = 10-30 %, 4 < pT < 8 GeV
  } else if(iAsymmetry == 2 && iCentrality == 1 && iTrackPt == 4){
    jetShapeHistogram->SetBinContent(12,0.89740300*1);     // Bin 12: 0.6 < DeltaR < 0.7
    
  // Cleaning configuration for bin xj integrted, C = 0-10 %, 1 < pT < 2 GeV
  } else if(iAsymmetry == 3 && iCentrality == 0 && iTrackPt == 1){
    jetShapeHistogram->SetBinContent(11,8.5798017*0.64);   // Bin 11: 0.5 < DeltaR < 0.6
    jetShapeHistogram->SetBinContent(12,6.9298858*0.64);   // Bin 12: 0.6 < DeltaR < 0.7
    
  // Cleaning configuration for bin xj integrted, C = 0-10 %, 3 < pT < 4 GeV
  } else if(iAsymmetry == 3 && iCentrality == 0 && iTrackPt == 3){
    jetShapeHistogram->SetBinContent(11,1.4697244*0.67);   // Bin 11: 0.5 < DeltaR < 0.6
    jetShapeHistogram->SetBinContent(13,0.85793671*0.71);  // Bin 13: 0.7 < DeltaR < 0.8
    
  // Cleaning configuration for bin xj integrted, C = 0-10 %, 4 < pT < 8 GeV
  } else if(iAsymmetry == 3 && iCentrality == 0 && iTrackPt == 4){
    jetShapeHistogram->SetBinContent(12,0.85350667*0.84);  // Bin 12: 0.6 < DeltaR < 0.7
    jetShapeHistogram->SetBinContent(13,0.55984141*0.92);  // Bin 13: 0.7 < DeltaR < 0.8
    
  // Cleaning configuration for bin xj integrted, C = 0-10 %, 8 < pT < 12 GeV
  } else if(iAsymmetry == 3 && iCentrality == 0 && iTrackPt == 5){
    jetShapeHistogram->SetBinContent(11,0.31237633*1);     // Bin 11: 0.5 < DeltaR < 0.6
    
  // Cleaning configuration for bin xj integrted, C = 10-30 %, 1 < pT < 2 GeV
  } else if(iAsymmetry == 3 && iCentrality == 1 && iTrackPt == 1){
    jetShapeHistogram->SetBinContent(13,4.4053483*0.83);   // Bin 13: 0.7 < DeltaR < 0.8
    
  // Cleaning configuration for bin xj integrted, C = 10-30 %, 2 < pT < 3 GeV
  } else if(iAsymmetry == 3 && iCentrality == 1 && iTrackPt == 2){
    jetShapeHistogram->SetBinContent(11,3.3867984*0.77);   // Bin 11: 0.5 < DeltaR < 0.6
    
  // Cleaning configuration for bin xj integrted, C = 10-30 %, 8 < pT < 12 GeV
  } else if(iAsymmetry == 3 && iCentrality == 1 && iTrackPt == 5){
    jetShapeHistogram->SetBinContent(10,0.66608731*0.9);   // Bin 10: 0.45 < DeltaR < 0.5
    
  // Cleaning configuration for bin xj integrted, C = 30-50 %, 0.7 < pT < 1 GeV
  } else if(iAsymmetry == 3 && iCentrality == 2 && iTrackPt == 0){
    jetShapeHistogram->SetBinContent(9,2.1125765*0.77);    // Bin 9: 0.4 < DeltaR < 0.45
    jetShapeHistogram->SetBinContent(10,1.9024753*0.79);   // Bin 10: 0.45 < DeltaR < 0.5
    
  }

}
