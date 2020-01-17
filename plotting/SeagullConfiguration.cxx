/*
 * Implementation of the SeagullConfiguration class
 */

// Own includes
#include "SeagullConfiguration.h"

/*
 * Contructor
 */
SeagullConfiguration::SeagullConfiguration()
{
  for(int iAsymmetry = 0; iAsymmetry < DijetHistogramManager::kMaxAsymmetryBins; iAsymmetry++){
    for(int iCentrality = 0; iCentrality < DijetHistogramManager::kMaxCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < DijetHistogramManager::kMaxTrackPtBins; iTrackPt++){
        trackLeadingJetPbPb[iAsymmetry][iCentrality][iTrackPt] = 0;
        trackSubleadingJetPbPb[iAsymmetry][iCentrality][iTrackPt] = 0;
        vetoTrackLeadingJetPbPb[iAsymmetry][iCentrality][iTrackPt] = 0;
        vetoTrackSubleadingJetPbPb[iAsymmetry][iCentrality][iTrackPt] = 0;
      } // Track pT loop
    } // Centrality loop
  } // Asymmetry loop
  
  InitializeArrays();
}

/*
 * Destructor
 */
SeagullConfiguration::~SeagullConfiguration(){

}

/*
 * Initialization for seagull method arrays. Set all bins that need nonzero value
 */
void SeagullConfiguration::InitializeArrays(){
  
  // ======================================================== //
  //                                                          //
  //  M     M   EEEEEE   TTTTTTT   H    H    OOOOO    DDDD    //
  //  MMM MMM   E           T      H    H   O     O   D   D   //
  //  M  M  M   EEEEEE      T      HHHHHH   O     O   D    D  //
  //  M     M   E           T      H    H   O     O   D   D   //
  //  M     M   EEEEEE      T      H    H    OOOOO    DDDD    //
  //                                                          //
  // ======================================================== //
  
  // Defaul seagull method has index 0. The meaning of different indices are explained below
  //
  //  0: Assume no dip in the middle and fit second order polynomial
  //  1: Assume dip in the middle and symmetrize background deltaEta. Fit with exponential function
  //  2: Same as 0, but use first order polynomial instead of second order
  //  3: Same as 1, but restrict fit to region 0-2
  //  4: Assume dip and falling tails of distribution and symmetrize deltaEta. Fit exponential function together with second order polynomial
  //  5: Same as 1, but do not symmetrize background deltaEta
  //  6: Same as 4, but do not symmetrize background deltaEta
  
  // ===========================================
  // == Track-leading jet histograms for PbPb ==
  // ===========================================
  
  // xj integrated
  trackLeadingJetPbPb[3][0][0] = 6;   // C = 0-10,  0.7 < pT < 1 GeV
  trackLeadingJetPbPb[3][0][1] = 6;   // C = 0-10,    1 < pT < 2 GeV
  trackLeadingJetPbPb[3][0][2] = 5;   // C = 0-10,    2 < pT < 3 GeV
  trackLeadingJetPbPb[3][0][3] = 1;   // C = 0-10,    3 < pT < 4 GeV
  trackLeadingJetPbPb[3][0][4] = 1;   // C = 0-10,    4 < pT < 8 GeV
  
  trackLeadingJetPbPb[3][1][0] = 6;   // C = 10-30, 0.7 < pT < 1 GeV
  trackLeadingJetPbPb[3][1][1] = 6;   // C = 10-30,   1 < pT < 2 GeV
  trackLeadingJetPbPb[3][1][2] = 5;   // C = 10-30,   2 < pT < 3 GeV
  trackLeadingJetPbPb[3][1][3] = 1;   // C = 10-30,   3 < pT < 4 GeV
  trackLeadingJetPbPb[3][1][4] = 1;   // C = 10-30,   4 < pT < 8 GeV
  
  trackLeadingJetPbPb[3][2][0] = 1;   // C = 30-50, 0.7 < pT < 1 GeV
  trackLeadingJetPbPb[3][2][1] = 1;   // C = 30-50,   1 < pT < 2 GeV
  trackLeadingJetPbPb[3][2][2] = 1;   // C = 30-50,   2 < pT < 3 GeV

  // 0 < xj < 0.6
  trackLeadingJetPbPb[0][0][0] = 1;   // C = 0-10,  0.7 < pT < 1 GeV
  trackLeadingJetPbPb[0][0][1] = 6;   // C = 0-10,    1 < pT < 2 GeV
  trackLeadingJetPbPb[0][0][2] = 1;   // C = 0-10,    2 < pT < 3 GeV
  trackLeadingJetPbPb[0][0][3] = 4;   // C = 0-10,    3 < pT < 4 GeV
  trackLeadingJetPbPb[0][0][4] = 1;   // C = 0-10,    4 < pT < 8 GeV
  
  trackLeadingJetPbPb[0][1][0] = 1;   // C = 10-30, 0.7 < pT < 1 GeV
  trackLeadingJetPbPb[0][1][1] = 1;   // C = 10-30,   1 < pT < 2 GeV
  trackLeadingJetPbPb[0][1][2] = 1;   // C = 10-30,   2 < pT < 3 GeV
  trackLeadingJetPbPb[0][1][3] = 1;   // C = 10-30,   3 < pT < 4 GeV
  trackLeadingJetPbPb[0][1][4] = 1;   // C = 10-30,   4 < pT < 8 GeV
  
  trackLeadingJetPbPb[0][2][1] = 1;   // C = 30-50,   1 < pT < 2 GeV
    
  // 0.6 < xj < 0.8
  trackLeadingJetPbPb[1][0][0] = 1;   // C = 0-10,  0.7 < pT < 1 GeV
  trackLeadingJetPbPb[1][0][1] = 6;   // C = 0-10,    1 < pT < 2 GeV
  trackLeadingJetPbPb[1][0][2] = 1;   // C = 0-10,    2 < pT < 3 GeV
  trackLeadingJetPbPb[1][0][3] = 1;   // C = 0-10,    3 < pT < 4 GeV
  
  trackLeadingJetPbPb[1][1][0] = 1;   // C = 10-30, 0.7 < pT < 1 GeV
  trackLeadingJetPbPb[1][1][1] = 5;   // C = 10-30,   1 < pT < 2 GeV
  trackLeadingJetPbPb[1][1][2] = 1;   // C = 10-30,   2 < pT < 3 GeV
  trackLeadingJetPbPb[1][1][3] = 1;   // C = 10-30,   3 < pT < 4 GeV
  
  trackLeadingJetPbPb[1][2][1] = 1;   // C = 30-50,   1 < pT < 2 GeV
  
  // 0.8 < xj < 1
  trackLeadingJetPbPb[2][0][1] = 6;   // C = 0-10,    1 < pT < 2 GeV
  trackLeadingJetPbPb[2][0][2] = 6;   // C = 0-10,    2 < pT < 3 GeV
  
  trackLeadingJetPbPb[2][1][1] = 6;   // C = 10-30,   1 < pT < 2 GeV
  trackLeadingJetPbPb[2][1][2] = 6;   // C = 10-30,   2 < pT < 3 GeV
  
  trackLeadingJetPbPb[2][2][1] = 1;   // C = 30-50,   1 < pT < 2 GeV  (Note: Do not use this for escheme)

  
  // ==============================================
  // == Track-subleading jet histograms for PbPb ==
  // ==============================================
  
  // xj integrated
  trackSubleadingJetPbPb[3][0][0] = 6;   // C = 0-10,  0.7 < pT < 1 GeV
  trackSubleadingJetPbPb[3][0][1] = 6;   // C = 0-10,    1 < pT < 2 GeV
  trackSubleadingJetPbPb[3][0][2] = 6;   // C = 0-10,    2 < pT < 3 GeV
  
  trackSubleadingJetPbPb[3][1][0] = 6;   // C = 10-30, 0.7 < pT < 1 GeV
  trackSubleadingJetPbPb[3][1][1] = 6;   // C = 10-30,   1 < pT < 2 GeV
  trackSubleadingJetPbPb[3][1][2] = 5;   // C = 10-30,   2 < pT < 3 GeV
  
  trackSubleadingJetPbPb[3][2][0] = 1;   // C = 30-50, 0.7 < pT < 1 GeV
  
  // 0 < xj < 0.6
  trackSubleadingJetPbPb[0][0][1] = 6;   // C = 0-10,    1 < pT < 2 GeV
  
  // 0.6 < xj < 0.8
  trackSubleadingJetPbPb[1][0][0] = 1;   // C = 0-10,  0.7 < pT < 1 GeV
  trackSubleadingJetPbPb[1][0][1] = 1;   // C = 0-10,    1 < pT < 2 GeV
  trackSubleadingJetPbPb[1][0][2] = 5;   // C = 0-10,    2 < pT < 3 GeV
  
  trackSubleadingJetPbPb[1][1][1] = 6;   // C = 10-30,   1 < pT < 2 GeV
  
  trackSubleadingJetPbPb[1][2][0] = 1;   // C = 30-50, 0.7 < pT < 2 GeV
  
  // 0.8 < xj < 1
  trackSubleadingJetPbPb[2][0][1] = 6;   // C = 0-10,    1 < pT < 2 GeV
  
  trackSubleadingJetPbPb[2][1][1] = 6;   // C = 10-30,   1 < pT < 2 GeV
  
  // ========================================== //
  //                                            //
  //  V       V   EEEEEEE   TTTTTTT    OOOOO    //
  //   V     V    E            T      O     O   //
  //    V   V     EEEEEEE      T      O     O   //
  //     V V      E            T      O     O   //
  //      V       EEEEEEE      T       OOOOO    //
  //                                            //
  // ========================================== //
  
  // Option to define veto for seagull correction in each bin. By default no veto is issued (0). The options are:
  //
  // 0: Regular correction => Do seagull correction if better chi2/ndf than for flat line
  // 1: Force correction => Do not compare chi2 to flat line
  // 2: Skip correction => No seagull correction applied
  
  // ===========================================
  // == Track-leading jet histograms for PbPb ==
  // ===========================================
  
  // xj integrated
  vetoTrackLeadingJetPbPb[3][0][5] = 2;   // C = 0-10,    8 < pT < 12 GeV
  vetoTrackLeadingJetPbPb[3][0][6] = 2;   // C = 0-10,   12 < pT < 300 GeV
  
  vetoTrackLeadingJetPbPb[3][1][5] = 2;   // C = 10-30,   8 < pT < 12 GeV
  vetoTrackLeadingJetPbPb[3][1][6] = 2;   // C = 10-30,  12 < pT < 300 GeV
  
  vetoTrackLeadingJetPbPb[3][2][5] = 2;   // C = 30-50,   8 < pT < 12 GeV
  vetoTrackLeadingJetPbPb[3][2][6] = 2;   // C = 30-50,  12 < pT < 300 GeV
  
  vetoTrackLeadingJetPbPb[3][3][5] = 2;   // C = 50-90,   8 < pT < 12 GeV
  vetoTrackLeadingJetPbPb[3][3][6] = 2;   // C = 50-90,  12 < pT < 300 GeV
  
  // 0 < xj < 0.6
  vetoTrackLeadingJetPbPb[0][0][5] = 2;   // C = 0-10,    8 < pT < 12 GeV
  vetoTrackLeadingJetPbPb[0][0][6] = 2;   // C = 0-10,   12 < pT < 300 GeV
  
  vetoTrackLeadingJetPbPb[0][1][5] = 2;   // C = 10-30,   8 < pT < 12 GeV
  vetoTrackLeadingJetPbPb[0][1][6] = 2;   // C = 10-30,  12 < pT < 300 GeV
  
  vetoTrackLeadingJetPbPb[0][2][2] = 2;   // C = 30-50,   2 < pT < 3 GeV
  vetoTrackLeadingJetPbPb[0][2][4] = 2;   // C = 30-50,   4 < pT < 8 GeV
  vetoTrackLeadingJetPbPb[0][2][5] = 2;   // C = 30-50,   8 < pT < 12 GeV
  vetoTrackLeadingJetPbPb[0][2][6] = 2;   // C = 30-50,  12 < pT < 300 GeV
  
  vetoTrackLeadingJetPbPb[0][3][2] = 2;   // C = 50-90,   2 < pT < 3 GeV
  vetoTrackLeadingJetPbPb[0][3][3] = 1;   // C = 50-90,   3 < pT < 4 GeV
  vetoTrackLeadingJetPbPb[0][3][4] = 2;   // C = 50-90,   4 < pT < 8 GeV
  vetoTrackLeadingJetPbPb[0][3][5] = 2;   // C = 50-90,   8 < pT < 12 GeV
  vetoTrackLeadingJetPbPb[0][3][6] = 2;   // C = 50-90,  12 < pT < 300 GeV
  
  // 0.6 < xj < 0.8
  vetoTrackLeadingJetPbPb[1][0][4] = 2;   // C = 0-10,    4 < pT < 8 GeV
  vetoTrackLeadingJetPbPb[1][0][5] = 2;   // C = 0-10,    8 < pT < 12 GeV
  vetoTrackLeadingJetPbPb[1][0][6] = 2;   // C = 0-10,   12 < pT < 300 GeV
  
  vetoTrackLeadingJetPbPb[1][1][4] = 2;   // C = 10-30,   4 < pT < 8 GeV
  vetoTrackLeadingJetPbPb[1][1][5] = 2;   // C = 10-30,   8 < pT < 12 GeV
  vetoTrackLeadingJetPbPb[1][1][6] = 2;   // C = 10-30,  12 < pT < 300 GeV
  
  vetoTrackLeadingJetPbPb[1][2][4] = 2;   // C = 30-50,   4 < pT < 8 GeV
  vetoTrackLeadingJetPbPb[1][2][5] = 2;   // C = 30-50,   8 < pT < 12 GeV
  vetoTrackLeadingJetPbPb[1][2][6] = 2;   // C = 30-50,  12 < pT < 300 GeV
  
  vetoTrackLeadingJetPbPb[1][3][5] = 2;   // C = 50-90,   8 < pT < 12 GeV
  vetoTrackLeadingJetPbPb[1][3][6] = 2;   // C = 50-90,  12 < pT < 300 GeV
  
  // 0.8 < xj < 1.0
  vetoTrackLeadingJetPbPb[2][0][5] = 2;   // C = 0-10,    8 < pT < 12 GeV
  vetoTrackLeadingJetPbPb[2][0][6] = 2;   // C = 0-10,   12 < pT < 300 GeV
  
  vetoTrackLeadingJetPbPb[2][1][5] = 2;   // C = 10-30,   8 < pT < 12 GeV
  vetoTrackLeadingJetPbPb[2][1][6] = 2;   // C = 10-30,  12 < pT < 300 GeV
  
  vetoTrackLeadingJetPbPb[2][2][5] = 2;   // C = 30-50,   8 < pT < 12 GeV
  vetoTrackLeadingJetPbPb[2][2][6] = 2;   // C = 30-50,  12 < pT < 300 GeV
  
  vetoTrackLeadingJetPbPb[2][3][4] = 2;   // C = 50-90,   4 < pT < 8 GeV
  vetoTrackLeadingJetPbPb[2][3][5] = 2;   // C = 50-90,   8 < pT < 12 GeV
  vetoTrackLeadingJetPbPb[2][3][6] = 2;   // C = 50-90,  12 < pT < 300 GeV
  
  // ==============================================
  // == Track-subleading jet histograms for PbPb ==
  // ==============================================
  
  // xj integrated
  vetoTrackSubleadingJetPbPb[3][0][5] = 2;   // C = 0-10,    8 < pT < 12 GeV
  vetoTrackSubleadingJetPbPb[3][0][6] = 2;   // C = 0-10,   12 < pT < 300 GeV
  
  vetoTrackSubleadingJetPbPb[3][1][5] = 2;   // C = 10-30,   8 < pT < 12 GeV
  vetoTrackSubleadingJetPbPb[3][1][6] = 2;   // C = 10-30,  12 < pT < 300 GeV
  
  vetoTrackSubleadingJetPbPb[3][2][5] = 2;   // C = 30-50,   8 < pT < 12 GeV
  vetoTrackSubleadingJetPbPb[3][2][6] = 2;   // C = 30-50,  12 < pT < 300 GeV
  
  vetoTrackSubleadingJetPbPb[3][3][5] = 2;   // C = 50-90,   8 < pT < 12 GeV
  vetoTrackSubleadingJetPbPb[3][3][6] = 2;   // C = 50-90,  12 < pT < 300 GeV
  
  // 0 < xj < 0.6
  vetoTrackSubleadingJetPbPb[0][0][5] = 2;   // C = 0-10,    8 < pT < 12 GeV
  vetoTrackSubleadingJetPbPb[0][0][6] = 2;   // C = 0-10,   12 < pT < 300 GeV
  
  vetoTrackSubleadingJetPbPb[0][1][5] = 2;   // C = 10-30,   8 < pT < 12 GeV
  vetoTrackSubleadingJetPbPb[0][1][6] = 2;   // C = 10-30,  12 < pT < 300 GeV
  
  vetoTrackSubleadingJetPbPb[0][2][3] = 1;   // C = 30-50,   3 < pT < 4 GeV
  vetoTrackSubleadingJetPbPb[0][2][4] = 2;   // C = 30-50,   4 < pT < 8 GeV
  vetoTrackSubleadingJetPbPb[0][2][5] = 2;   // C = 30-50,   8 < pT < 12 GeV
  vetoTrackSubleadingJetPbPb[0][2][6] = 2;   // C = 30-50,  12 < pT < 300 GeV
  
  vetoTrackSubleadingJetPbPb[0][3][1] = 2;   // C = 50-90,   1 < pT < 2 GeV
  vetoTrackSubleadingJetPbPb[0][3][2] = 1;   // C = 50-90,   2 < pT < 3 GeV
  vetoTrackSubleadingJetPbPb[0][3][5] = 2;   // C = 50-90,   8 < pT < 12 GeV
  vetoTrackSubleadingJetPbPb[0][3][6] = 2;   // C = 50-90,  12 < pT < 300 GeV
  
  // 0.6 < xj < 0.8
  vetoTrackSubleadingJetPbPb[1][0][1] = 1;   // C = 0-10,    1 < pT < 2 GeV
  vetoTrackSubleadingJetPbPb[1][0][2] = 1;   // C = 0-10,    2 < pT < 3 GeV
  vetoTrackSubleadingJetPbPb[1][0][5] = 2;   // C = 0-10,    8 < pT < 12 GeV
  vetoTrackSubleadingJetPbPb[1][0][6] = 2;   // C = 0-10,   12 < pT < 300 GeV
  
  vetoTrackSubleadingJetPbPb[1][1][3] = 2;   // C = 10-30,   3 < pT < 4 GeV
  vetoTrackSubleadingJetPbPb[1][1][5] = 2;   // C = 10-30,   8 < pT < 12 GeV
  vetoTrackSubleadingJetPbPb[1][1][6] = 2;   // C = 10-30,  12 < pT < 300 GeV
  
  vetoTrackSubleadingJetPbPb[1][2][1] = 2;   // C = 30-50,   1 < pT < 2 GeV
  vetoTrackSubleadingJetPbPb[1][2][3] = 2;   // C = 30-50,   3 < pT < 4 GeV
  vetoTrackSubleadingJetPbPb[1][2][5] = 2;   // C = 30-50,   8 < pT < 12 GeV
  vetoTrackSubleadingJetPbPb[1][2][6] = 2;   // C = 30-50,  12 < pT < 300 GeV
  
  vetoTrackSubleadingJetPbPb[1][3][3] = 2;   // C = 50-90,   3 < pT < 4 GeV
  vetoTrackSubleadingJetPbPb[1][3][4] = 2;   // C = 50-90,   4 < pT < 8 GeV
  vetoTrackSubleadingJetPbPb[1][3][5] = 2;   // C = 50-90,   8 < pT < 12 GeV
  vetoTrackSubleadingJetPbPb[1][3][6] = 2;   // C = 50-90,  12 < pT < 300 GeV
  
  // 0.8 < xj < 1.0
  vetoTrackSubleadingJetPbPb[2][0][5] = 2;   // C = 0-10,    8 < pT < 12 GeV
  vetoTrackSubleadingJetPbPb[2][0][6] = 2;   // C = 0-10,   12 < pT < 300 GeV
  
  vetoTrackSubleadingJetPbPb[2][1][0] = 2;   // C = 10-30, 0.7 < pT < 1 GeV
  vetoTrackSubleadingJetPbPb[2][1][3] = 2;   // C = 10-30,   3 < pT < 4 GeV
  vetoTrackSubleadingJetPbPb[2][1][4] = 2;   // C = 10-30,   4 < pT < 8 GeV
  vetoTrackSubleadingJetPbPb[2][1][5] = 2;   // C = 10-30,   8 < pT < 12 GeV
  vetoTrackSubleadingJetPbPb[2][1][6] = 2;   // C = 10-30,  12 < pT < 300 GeV
  
  vetoTrackSubleadingJetPbPb[2][2][0] = 2;   // C = 30-50, 0.7 < pT < 1 GeV
  vetoTrackSubleadingJetPbPb[2][2][3] = 1;   // C = 30-50,   3 < pT < 4 GeV
  vetoTrackSubleadingJetPbPb[2][2][4] = 2;   // C = 30-50,   4 < pT < 8 GeV
  vetoTrackSubleadingJetPbPb[2][2][5] = 2;   // C = 30-50,   8 < pT < 12 GeV
  vetoTrackSubleadingJetPbPb[2][2][6] = 2;   // C = 30-50,  12 < pT < 300 GeV
  
  vetoTrackSubleadingJetPbPb[2][3][0] = 2;   // C = 50-90, 0.7 < pT < 1 GeV
  vetoTrackSubleadingJetPbPb[2][3][4] = 2;   // C = 50-90,   4 < pT < 8 GeV
  vetoTrackSubleadingJetPbPb[2][3][5] = 2;   // C = 50-90,   8 < pT < 12 GeV
  vetoTrackSubleadingJetPbPb[2][3][6] = 2;   // C = 50-90,  12 < pT < 300 GeV
  
}

/*
 * Get the seagull method for the given bin
 *
 *  const int iJetTrack = Index for jet track correlation
 *  const int iAsymmetry = Index for asymmetry bin
 *  const int iCentrality = Index for the centrality bin
 *  const int iTrackPt = Index for the track pT bin
 */
int SeagullConfiguration::GetSeagullMethod(const int iJetTrack, const int iAsymmetry, const int iCentrality, const int iTrackPt) const{
  
  // Read the value from correct table based on the provided jet track index
  switch (iJetTrack) {
    case DijetHistogramManager::kTrackLeadingJet:
      return trackLeadingJetPbPb[iAsymmetry][iCentrality][iTrackPt];
      break;
      
    case DijetHistogramManager::kPtWeightedTrackLeadingJet:
      return trackLeadingJetPbPb[iAsymmetry][iCentrality][iTrackPt];
      break;
      
    case DijetHistogramManager::kTrackSubleadingJet:
      return trackSubleadingJetPbPb[iAsymmetry][iCentrality][iTrackPt];
      break;
      
    case DijetHistogramManager::kPtWeightedTrackSubleadingJet:
      return trackSubleadingJetPbPb[iAsymmetry][iCentrality][iTrackPt];
      break;
      
    case DijetHistogramManager::kTrackInclusiveJet:
      return trackLeadingJetPbPb[iAsymmetry][iCentrality][iTrackPt];
      break;
      
    case DijetHistogramManager::kPtWeightedTrackInclusiveJet:
      return trackLeadingJetPbPb[iAsymmetry][iCentrality][iTrackPt];
      break;
      
    default:
      return 0;
      break;
  }
  
  // This should never happen, but does not hurt to have.
  return 0;

}

/*
 * Get the seagull veto flag for the given bin
 *
 *  const int iJetTrack = Index for jet track correlation
 *  const int iAsymmetry = Index for asymmetry bin
 *  const int iCentrality = Index for the centrality bin
 *  const int iTrackPt = Index for the track pT bin
 */
int SeagullConfiguration::GetSeagullVeto(const int iJetTrack, const int iAsymmetry, const int iCentrality, const int iTrackPt) const{
  
  // Read the value from correct table based on the provided jet track index
  switch (iJetTrack) {
    case DijetHistogramManager::kTrackLeadingJet:
      return vetoTrackLeadingJetPbPb[iAsymmetry][iCentrality][iTrackPt];
      break;
      
    case DijetHistogramManager::kPtWeightedTrackLeadingJet:
      return vetoTrackLeadingJetPbPb[iAsymmetry][iCentrality][iTrackPt];
      break;
      
    case DijetHistogramManager::kTrackSubleadingJet:
      return vetoTrackSubleadingJetPbPb[iAsymmetry][iCentrality][iTrackPt];
      break;
      
    case DijetHistogramManager::kPtWeightedTrackSubleadingJet:
      return vetoTrackSubleadingJetPbPb[iAsymmetry][iCentrality][iTrackPt];
      break;
      
    case DijetHistogramManager::kTrackInclusiveJet:
      return vetoTrackLeadingJetPbPb[iAsymmetry][iCentrality][iTrackPt];
      break;
      
    case DijetHistogramManager::kPtWeightedTrackInclusiveJet:
      return vetoTrackLeadingJetPbPb[iAsymmetry][iCentrality][iTrackPt];
      break;
      
    default:
      return 0;
      break;
  }
  
  // This should never happen, but does not hurt to have.
  return 0;

}
