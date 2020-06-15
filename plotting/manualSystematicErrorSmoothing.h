/*
 * Configuration to manually smoothen systematic uncertainties in some fluctuating bins 
 */

/*
 *  Systematic uncertainty smoothening for jet resolution
 *
 *  Files used to derive the smoothing factors with compareDijetHistograms.C:
 *   dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_allCorrections_processed_2020-02-17.root
 *   dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_20eveMix_20pSmear_wtaAxis_allCorrectionsUnsmeared_processed_2020-05-13.root
 *
 *  Arguments:
 *   int iJetTrack = Index for the jet track correlation
 *   int iAsymmetry = Index for the dijet momentum balance bin
 *   int iCentrality = Index for the centrality bin
 *   int iTrackPt = Index for the track pT bin
 *   int iBin = Index for the bin in the histogram
 *
 *  return: Factor with which the uncertainty should be divided to iron out fluctuating bins
 *
 */
double getSmoothingJetResolution(int iJetTrack, int iAsymmetry, int iCentrality, int iTrackPt, int iBin){

  const int nJetTrack = 8;
  const int nAsymmetry = 4;
  const int nCentrality = 4;
  const int nTrackPt = 7;
  const int nBin = 15;

  double resolutionTable[nJetTrack][nAsymmetry][nCentrality][nTrackPt][nBin];

  // Initialize the resolution table to one
  for(int iJetTrack = 0; iJetTrack < nJetTrack; iJetTrack++){
    for(int iAsymmetry = 0; iAsymmetry < nAsymmetry; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentrality; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPt; iTrackPt++){
          for(int iBin = 0; iBin < nBin; iBin++){
            resolutionTable[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iBin] = 1;
          }
        }
      }
    }
  }
  
  // Manual configuration here
  
  // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  // ================================================================ //
  //                       Leading jet shape
  // ================================================================ //
  
  // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  // 0.0 < xj < 0.6 | Centrality: 0-10 | 0.7 < pT < 1 GeV
  resolutionTable[2][0][0][0][10] = 3;
  resolutionTable[2][0][0][0][11] = 2;
  resolutionTable[2][0][0][0][13] = 2;
  resolutionTable[2][0][0][0][14] = 4;
  
  // 0.0 < xj < 0.6 | Centrality: 0-10 | 1 < pT < 2 GeV
  resolutionTable[2][0][0][1][11] = 2;
  resolutionTable[2][0][0][1][14] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 0-10 | 2 < pT < 3 GeV
  resolutionTable[2][0][0][2][12] = 5;
  
  // 0.0 < xj < 0.6 | Centrality: 0-10 | 3 < pT < 4 GeV
  resolutionTable[2][0][0][3][10] = 2;
  resolutionTable[2][0][0][3][11] = 3;
  resolutionTable[2][0][0][3][13] = 3;
  resolutionTable[2][0][0][3][14] = 50;
  
  // 0.0 < xj < 0.6 | Centrality: 30-50 | 0.7 < pT < 1 GeV
  resolutionTable[2][0][2][0][8] = 2;
  resolutionTable[2][0][2][0][10] = 3;
  resolutionTable[2][0][2][0][12] = 3;
  resolutionTable[2][0][2][0][14] = 4;
  
  // 0.0 < xj < 0.6 | Centrality: 30-50 | 1 < pT < 2 GeV
  resolutionTable[2][0][2][1][11] = 2;
  resolutionTable[2][0][2][1][12] = 3;
  resolutionTable[2][0][2][1][13] = 2;
  resolutionTable[2][0][2][1][14] = 4;
  
  // 0.0 < xj < 0.6 | Centrality: 30-50 | 2 < pT < 3 GeV
  resolutionTable[2][0][2][2][6] = 2;
  resolutionTable[2][0][2][2][7] = 2;
  resolutionTable[2][0][2][2][11] = 2;
  resolutionTable[2][0][2][2][13] = 3;
  
  // 0.0 < xj < 0.6 | Centrality: 30-50 | 3 < pT < 4 GeV
  resolutionTable[2][0][2][3][6] = 2;
  resolutionTable[2][0][2][3][10] = 2;
  resolutionTable[2][0][2][3][12] = 4;
  resolutionTable[2][0][2][3][13] = 2;
  resolutionTable[2][0][2][3][14] = 4;
  
  // 0.0 < xj < 0.6 | Centrality: 30-50 | 4 < pT < 8 GeV
  resolutionTable[2][0][2][4][13] = 3;
  resolutionTable[2][0][2][4][14] = 3;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 0.7 < pT < 1 GeV
  resolutionTable[2][0][3][0][13] = 4;
  resolutionTable[2][0][3][0][14] = 4 ;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 1 < pT < 2 GeV
  resolutionTable[2][0][3][1][8] = 2;
  resolutionTable[2][0][3][1][10] = 2;
  resolutionTable[2][0][3][1][11] = 2;
  resolutionTable[2][0][3][1][13] = 2;
  resolutionTable[2][0][3][1][14] = 4;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 2 < pT < 3 GeV
  resolutionTable[2][0][3][2][6] = 2;
  resolutionTable[2][0][3][2][12] = 3;
  resolutionTable[2][0][3][2][13] = 2;
  resolutionTable[2][0][3][2][14] = 3;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 3 < pT < 4 GeV
  resolutionTable[2][0][3][3][10] = 2;
  resolutionTable[2][0][3][3][12] = 2;
  resolutionTable[2][0][3][3][13] = 3;
  resolutionTable[2][0][3][3][14] = 3;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 4 < pT < 8 GeV
  resolutionTable[2][0][3][4][5] = 2;
  resolutionTable[2][0][3][4][7] = 2;
  resolutionTable[2][0][3][4][12] = 2;
  resolutionTable[2][0][3][4][14] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 8 < pT < 12 GeV
  resolutionTable[2][0][3][5][5] = 2;
  resolutionTable[2][0][3][5][9] = 3;
  resolutionTable[2][0][3][5][10] = 3;
  resolutionTable[2][0][3][5][12] = 2;
  resolutionTable[2][0][3][5][14] = 8;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 12 < pT < 300 GeV
  resolutionTable[2][0][3][6][5] = 3;
  resolutionTable[2][0][3][6][7] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 0-10 | 0.7 < pT < 1 GeV
  resolutionTable[2][1][0][0][13] = 4;
  resolutionTable[2][1][0][0][14] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 0-10 | 1 < pT < 2 GeV
  resolutionTable[2][1][0][1][14] = 4;
  
  // 0.6 < xj < 0.8 | Centrality: 0-10 | 2 < pT < 3 GeV
  resolutionTable[2][1][0][2][7] = 2;
  resolutionTable[2][1][0][2][8] = 2;
  resolutionTable[2][1][0][2][12] = 4;
  resolutionTable[2][1][0][2][13] = 3;
  resolutionTable[2][1][0][2][14] = 3;
  
  // 0.6 < xj < 0.8 | Centrality: 0-10 | 4 < pT < 8 GeV
  resolutionTable[2][1][0][4][6] = 2;
  resolutionTable[2][1][0][4][9] = 2;
  resolutionTable[2][1][0][4][12] = 2;
  resolutionTable[2][1][0][4][13] = 4;
  
  // 0.6 < xj < 0.8 | Centrality: 10-30 | 0.7 < pT < 1 GeV
  resolutionTable[2][1][1][0][2] = 2;
  resolutionTable[2][1][1][0][6] = 2;
  resolutionTable[2][1][1][0][8] = 2;
  resolutionTable[2][1][1][0][10] = 2;
  resolutionTable[2][1][1][0][12] = 2;
  resolutionTable[2][1][1][0][13] = 2;
  resolutionTable[2][1][1][0][14] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 10-30 | 2 < pT < 3 GeV
  resolutionTable[2][1][1][2][10] = 2;
  resolutionTable[2][1][1][2][11] = 2;
  resolutionTable[2][1][1][2][12] = 2;
  resolutionTable[2][1][1][2][13] = 2;
  resolutionTable[2][1][1][2][14] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 10-30 | 3 < pT < 4 GeV
  resolutionTable[2][1][1][3][10] = 2;
  resolutionTable[2][1][1][3][11] = 2;
  resolutionTable[2][1][1][3][12] = 5;
  resolutionTable[2][1][1][3][13] = 5;
  resolutionTable[2][1][1][3][14] = 5;
  
  // 0.6 < xj < 0.8 | Centrality: 30-50 | 0.7 < pT < 1 GeV
  resolutionTable[2][1][2][0][3] = 2;
  resolutionTable[2][1][2][0][8] = 2;
  resolutionTable[2][1][2][0][10] = 3;
  resolutionTable[2][1][2][0][11] = 3;
  resolutionTable[2][1][2][0][12] = 2;
  resolutionTable[2][1][2][0][13] = 2;
  resolutionTable[2][1][2][0][14] = 4;
  
  // 0.6 < xj < 0.8 | Centrality: 30-50 | 1 < pT < 2 GeV
  resolutionTable[2][1][2][1][10] = 2;
  resolutionTable[2][1][2][1][13] = 4;
  resolutionTable[2][1][2][1][14] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 30-50 | 2 < pT < 3 GeV
  resolutionTable[2][1][2][2][10] = 4;
  resolutionTable[2][1][2][2][13] = 2;
  resolutionTable[2][1][2][2][14] = 8;
  
  // 0.6 < xj < 0.8 | Centrality: 30-50 | 3 < pT < 4 GeV
  resolutionTable[2][1][2][3][10] = 4;
  resolutionTable[2][1][2][3][12] = 4;
  resolutionTable[2][1][2][3][14] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 30-50 | 4 < pT < 8 GeV
  resolutionTable[2][1][2][4][14] = 4;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 0.7 < pT < 1 GeV
  resolutionTable[2][1][3][0][2] = 3;
  resolutionTable[2][1][3][0][5] = 2;
  resolutionTable[2][1][3][0][6] = 2;
  resolutionTable[2][1][3][0][7] = 2;
  resolutionTable[2][1][3][0][9] = 2;
  resolutionTable[2][1][3][0][10] = 2;
  resolutionTable[2][1][3][0][11] = 2;
  resolutionTable[2][1][3][0][12] = 3;
  resolutionTable[2][1][3][0][14] = 4;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 1 < pT < 2 GeV
  resolutionTable[2][1][3][1][8] = 3;
  resolutionTable[2][1][3][1][9] = 3;
  resolutionTable[2][1][3][1][10] = 3;
  resolutionTable[2][1][3][1][11] = 3;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 3 < pT < 4 GeV
  resolutionTable[2][1][3][3][9] = 2;
  resolutionTable[2][1][3][3][10] = 2;
  resolutionTable[2][1][3][3][11] = 2;
  resolutionTable[2][1][3][3][12] = 2;
  resolutionTable[2][1][3][3][13] = 2;
  resolutionTable[2][1][3][3][14] = 4;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 4 < pT < 8 GeV
  resolutionTable[2][1][3][4][13] = 3;
  resolutionTable[2][1][3][4][14] = 4;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 8 < pT < 12 GeV
  resolutionTable[2][1][3][4][7] = 2;
  resolutionTable[2][1][3][4][8] = 2;
  resolutionTable[2][1][3][4][9] = 3;
  resolutionTable[2][1][3][4][10] = 3;
  resolutionTable[2][1][3][4][11] = 3;
  resolutionTable[2][1][3][4][12] = 3;
  resolutionTable[2][1][3][4][13] = 2;
  resolutionTable[2][1][3][4][14] = 8;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 1 < pT < 2 GeV
  resolutionTable[2][2][0][1][6] = 2;
  resolutionTable[2][2][0][1][12] = 2;
  resolutionTable[2][2][0][1][13] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 2 < pT < 3 GeV
  resolutionTable[2][2][0][2][6] = 2;
  resolutionTable[2][2][0][2][13] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 3 < pT < 4 GeV
  resolutionTable[2][2][0][3][9] = 2;
  resolutionTable[2][2][0][3][10] = 4;
  resolutionTable[2][2][0][3][13] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 4 < pT < 8 GeV
  resolutionTable[2][2][0][4][13] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 8 < pT < 12 GeV
  resolutionTable[2][2][0][5][6] = 2;
  resolutionTable[2][2][0][5][14] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 12 < pT < 300 GeV
  resolutionTable[2][2][0][6][12] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 10-30 | 0.7 < pT < 1 GeV
  resolutionTable[2][2][1][0][5] = 2;
  resolutionTable[2][2][1][0][7] = 2;
  resolutionTable[2][2][1][0][12] = 2;
  resolutionTable[2][2][1][0][14] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 10-30 | 1 < pT < 2 GeV
  resolutionTable[2][2][1][1][11] = 2;
  resolutionTable[2][2][1][1][13] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 10-30 | 2 < pT < 3 GeV
  resolutionTable[2][2][1][2][9] = 2;
  resolutionTable[2][2][1][2][12] = 4;
  resolutionTable[2][2][1][2][13] = 3;
  resolutionTable[2][2][1][2][14] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 10-30 | 3 < pT < 4 GeV
  resolutionTable[2][2][1][3][7] = 2;
  resolutionTable[2][2][1][3][11] = 4;
  resolutionTable[2][2][1][3][13] = 10;
  resolutionTable[2][2][1][3][14] = 10;
  
  // 0.8 < xj < 1.0 | Centrality: 10-30 | 8 < pT < 12 GeV
  resolutionTable[2][2][1][5][11] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 0.7 < pT < 1 GeV
  resolutionTable[2][2][2][0][1] = 3;
  resolutionTable[2][2][2][0][3] = 3;
  resolutionTable[2][2][2][0][6] = 2;
  resolutionTable[2][2][2][0][7] = 2;
  resolutionTable[2][2][2][0][10] = 2;
  resolutionTable[2][2][2][0][11] = 2;
  resolutionTable[2][2][2][0][12] = 2;
  resolutionTable[2][2][2][0][13] = 6;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 1 < pT < 2 GeV
  resolutionTable[2][2][2][1][1] = 2;
  resolutionTable[2][2][2][1][9] = 2;
  resolutionTable[2][2][2][1][12] = 2;
  resolutionTable[2][2][2][1][13] = 2;
  resolutionTable[2][2][2][1][14] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 2 < pT < 3 GeV
  resolutionTable[2][2][2][2][1] = 2;
  resolutionTable[2][2][2][2][2] = 2;
  resolutionTable[2][2][2][2][12] = 2;
  resolutionTable[2][2][2][2][13] = 2;
  resolutionTable[2][2][2][2][14] = 8;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 3 < pT < 4 GeV
  resolutionTable[2][2][2][3][8] = 2;
  resolutionTable[2][2][2][3][9] = 2;
  resolutionTable[2][2][2][3][10] = 2;
  resolutionTable[2][2][2][3][11] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 4 < pT < 8 GeV
  resolutionTable[2][2][2][4][11] = 2;
  resolutionTable[2][2][2][4][12] = 2;
  resolutionTable[2][2][2][4][13] = 2;
  resolutionTable[2][2][2][4][14] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 8 < pT < 12 GeV
  resolutionTable[2][2][2][5][11] = 2;
  resolutionTable[2][2][2][5][12] = 3;
  resolutionTable[2][2][2][5][13] = 4;
  resolutionTable[2][2][2][5][14] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 0.7 < pT < 1 GeV
  resolutionTable[2][2][3][0][3] = 3;
  resolutionTable[2][2][3][0][5] = 2;
  resolutionTable[2][2][3][0][6] = 2;
  resolutionTable[2][2][3][0][8] = 2;
  resolutionTable[2][2][3][0][13] = 3;
  resolutionTable[2][2][3][0][14] = 5;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 1 < pT < 2 GeV
  resolutionTable[2][2][3][1][9] = 2;
  resolutionTable[2][2][3][1][10] = 2;
  resolutionTable[2][2][3][1][11] = 2;
  resolutionTable[2][2][3][1][12] = 2;
  resolutionTable[2][2][3][1][14] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 2 < pT < 3 GeV
  resolutionTable[2][2][3][2][10] = 3;
  resolutionTable[2][2][3][2][13] = 2;
  resolutionTable[2][2][3][2][14] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 3 < pT < 4 GeV
  resolutionTable[2][2][3][3][3] = 3;
  resolutionTable[2][2][3][3][9] = 2;
  resolutionTable[2][2][3][3][13] = 6;
  resolutionTable[2][2][3][3][14] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 4 < pT < 8 GeV
  resolutionTable[2][2][3][4][7] = 2;
  resolutionTable[2][2][3][4][8] = 2;
  resolutionTable[2][2][3][4][10] = 2;
  resolutionTable[2][2][3][4][13] = 3;
  
  // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  // ================================================================ //
  //                     Subleading jet shape
  // ================================================================ //
  
  // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  // xj integrated | Centrality 10-30 | 2 < pT < 3 GeV
  resolutionTable[5][3][1][2][12] = 2;
  resolutionTable[5][3][1][2][14] = 3;
  
  // xj integrated | Centrality 30-50 | 0.7 < pT < 1 GeV
  resolutionTable[5][3][2][0][12] = 2;
  resolutionTable[5][3][2][0][13] = 2;
  resolutionTable[5][3][2][0][14] = 2;
  
  // xj integrated | Centrality 30-50 | 2 < pT < 3 GeV
  resolutionTable[5][3][2][2][12] = 2;
  resolutionTable[5][3][2][2][14] = 4;
  
  // xj integrated | Centrality 30-50 | 4 < pT < 8 GeV
  resolutionTable[5][3][2][4][12] = 2;
  resolutionTable[5][3][2][4][14] = 2;
  
  // xj integrated | Centrality 50-90 | 1 < pT < 2 GeV
  resolutionTable[5][3][2][1][11] = 2;
  resolutionTable[5][3][2][1][12] = 2;
  resolutionTable[5][3][2][1][13] = 2;
  resolutionTable[5][3][2][1][14] = 2;
  
  // xj integrated | Centrality 50-90 | 3 < pT < 4 GeV
  resolutionTable[5][3][2][3][13] = 2;
  resolutionTable[5][3][2][3][14] = 4;
  
  // 0.0 < xj < 0.6 | Centrality: 0-10 | 1 < pT < 2 GeV
  resolutionTable[5][0][0][1][9] = 2;
  resolutionTable[5][0][0][1][10] = 2;
  resolutionTable[5][0][0][1][11] = 2;
  resolutionTable[5][0][0][1][13] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 0-10 | 2 < pT < 3 GeV
  resolutionTable[2][0][0][2][14] = 3;
  
  // 0.0 < xj < 0.6 | Centrality: 30-50 | 0.7 < pT < 1 GeV
  resolutionTable[5][0][2][0][6] = 2;
  resolutionTable[5][0][2][0][7] = 2;
  resolutionTable[5][0][2][0][8] = 2;
  resolutionTable[5][0][2][0][9] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 30-50 | 1 < pT < 2 GeV
  resolutionTable[5][0][2][1][13] = 2;
  resolutionTable[5][0][2][1][14] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 30-50 | 2 < pT < 3 GeV
  resolutionTable[5][0][2][2][9] = 2;
  resolutionTable[5][0][2][2][10] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 0.7 < pT < 1 GeV
  resolutionTable[5][0][3][0][5] = 2;
  resolutionTable[5][0][3][0][14] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 2 < pT < 3 GeV
  resolutionTable[5][0][3][2][7] = 2;
  resolutionTable[5][0][3][2][8] = 2;
  resolutionTable[5][0][3][2][9] = 2;
  resolutionTable[5][0][3][2][14] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 3 < pT < 4 GeV
  resolutionTable[5][0][3][3][12] = 2;
  resolutionTable[5][0][3][3][14] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 4 < pT < 8 GeV
  resolutionTable[5][0][3][4][8] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 8 < pT < 12 GeV
  resolutionTable[5][0][3][5][3] = 2;
  resolutionTable[5][0][3][5][6] = 4;
  
  // 0.6 < xj < 0.8 | Centrality: 0-10 | 0.7 < pT < 1 GeV
  resolutionTable[5][1][0][0][11] = 2;
  resolutionTable[5][1][0][0][13] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 0-10 | 1 < pT < 2 GeV
  resolutionTable[5][1][0][1][7] = 2;
  resolutionTable[5][1][0][1][8] = 2;
  resolutionTable[5][1][0][1][9] = 3;
  resolutionTable[5][1][0][1][12] = 3;
  
  // 0.6 < xj < 0.8 | Centrality: 0-10 | 2 < pT < 3 GeV
  resolutionTable[5][1][0][2][2] = 2;
  resolutionTable[5][1][0][2][3] = 2;
  resolutionTable[5][1][0][2][4] = 2;
  resolutionTable[5][1][0][2][5] = 2;
  resolutionTable[5][1][0][2][6] = 2;
  resolutionTable[5][1][0][2][7] = 2;
  resolutionTable[5][1][0][2][8] = 2;
  resolutionTable[5][1][0][2][9] = 2;
  resolutionTable[5][1][0][2][10] = 2;
  resolutionTable[5][1][0][2][11] = 3;
  resolutionTable[5][1][0][2][12] = 3;
  resolutionTable[5][1][0][2][14] = 3;
  
  // 0.6 < xj < 0.8 | Centrality: 0-10 | 3 < pT < 4 GeV
  resolutionTable[5][1][0][3][12] = 3;
  resolutionTable[5][1][0][3][14] = 3;
  
  // 0.6 < xj < 0.8 | Centrality: 0-10 | 4 < pT < 8 GeV
  resolutionTable[5][1][0][4][8] = 3;
  
  // 0.6 < xj < 0.8 | Centrality: 10-30 | 0.7 < pT < 1 GeV
  resolutionTable[5][1][1][0][9] = 2;
  resolutionTable[5][1][1][0][10] = 2;
  resolutionTable[5][1][1][0][12] = 3;
  resolutionTable[5][1][1][0][14] = 3;
  
  // 0.6 < xj < 0.8 | Centrality: 10-30 | 1 < pT < 2 GeV
  resolutionTable[5][1][1][1][11] = 2;
  resolutionTable[5][1][1][1][13] = 2;
  resolutionTable[5][1][1][1][14] = 4;
  
  // 0.6 < xj < 0.8 | Centrality: 10-30 | 2 < pT < 3 GeV
  resolutionTable[5][1][1][2][9] = 2;
  resolutionTable[5][1][1][2][10] = 2;
  resolutionTable[5][1][1][2][13] = 3;
  
  // 0.6 < xj < 0.8 | Centrality: 10-30 | 3 < pT < 4 GeV
  resolutionTable[5][1][1][3][13] = 3;
  
  // 0.6 < xj < 0.8 | Centrality: 30-50 | 2 < pT < 3 GeV
  resolutionTable[5][1][2][2][10] = 2;
  resolutionTable[5][1][2][2][11] = 2;
  resolutionTable[5][1][2][2][12] = 3;
  resolutionTable[5][1][2][2][14] = 3;
  
  // 0.6 < xj < 0.8 | Centrality: 30-50 | 3 < pT < 4 GeV
  resolutionTable[5][1][2][3][13] = 2;
  resolutionTable[5][1][2][3][14] = 3;
  
  // 0.6 < xj < 0.8 | Centrality: 30-50 | 4 < pT < 8 GeV
  resolutionTable[5][1][2][4][12] = 3;
  resolutionTable[5][1][2][4][13] = 3;
  resolutionTable[5][1][2][4][14] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 30-50 | 8 < pT < 12 GeV
  resolutionTable[5][1][2][5][12] = 2;
  resolutionTable[5][1][2][5][13] = 2;
  resolutionTable[5][1][2][5][14] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 0.7 < pT < 1 GeV
  resolutionTable[5][1][3][0][1] = 3;
  resolutionTable[5][1][3][0][4] = 2;
  resolutionTable[5][1][3][0][7] = 2;
  resolutionTable[5][1][3][0][10] = 3;
  resolutionTable[5][1][3][0][13] = 4;
  resolutionTable[5][1][3][0][14] = 3;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 1 < pT < 2 GeV
  resolutionTable[5][1][3][1][5] = 2;
  resolutionTable[5][1][3][1][6] = 2;
  resolutionTable[5][1][3][1][13] = 4;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 2 < pT < 3 GeV
  resolutionTable[5][1][3][2][5] = 2;
  resolutionTable[5][1][3][2][6] = 2;
  resolutionTable[5][1][3][2][8] = 3;
  resolutionTable[5][1][3][2][14] = 4;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 3 < pT < 4 GeV
  resolutionTable[5][1][3][3][9] = 1.5;
  resolutionTable[5][1][3][3][13] = 3;
  resolutionTable[5][1][3][3][14] = 3;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 4 < pT < 8 GeV
  resolutionTable[5][1][3][4][7] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 1 < pT < 2 GeV
  resolutionTable[5][2][0][1][1] = 3;
  resolutionTable[5][2][0][1][2] = 3;
  resolutionTable[5][2][0][1][3] = 3;
  resolutionTable[5][2][0][1][4] = 3;
  resolutionTable[5][2][0][1][5] = 3;
  resolutionTable[5][2][0][1][6] = 3;
  resolutionTable[5][2][0][1][7] = 3;
  resolutionTable[5][2][0][1][8] = 3;
  resolutionTable[5][2][0][1][9] = 3;
  resolutionTable[5][2][0][1][10] = 3;
  resolutionTable[5][2][0][1][11] = 3;
  resolutionTable[5][2][0][1][12] = 4;
  resolutionTable[5][2][0][1][13] = 4;
  resolutionTable[5][2][0][1][14] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 2 < pT < 3 GeV
  resolutionTable[5][2][0][2][13] = 4;
  resolutionTable[5][2][0][2][14] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 3 < pT < 4 GeV
  resolutionTable[5][2][0][3][9] = 3;
  resolutionTable[5][2][0][3][13] = 2;
  resolutionTable[5][2][0][3][14] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 4 < pT < 8 GeV
  resolutionTable[5][2][0][4][12] = 3;
  resolutionTable[5][2][0][4][13] = 3;
  resolutionTable[5][2][0][4][14] = 3;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 8 < pT < 12 GeV
  resolutionTable[5][2][0][5][10] = 4;
  resolutionTable[5][2][0][5][12] = 3;
  resolutionTable[5][2][0][5][13] = 3;
  resolutionTable[5][2][0][5][14] = 3;
  
  // 0.8 < xj < 1.0 | Centrality: 10-30 | 0.7 < pT < 1 GeV
  resolutionTable[5][2][1][0][9] = 2;
  resolutionTable[5][2][1][0][10] = 2;
  resolutionTable[5][2][1][0][11] = 2;
  resolutionTable[5][2][1][0][12] = 2;
  resolutionTable[5][2][1][0][13] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 10-30 | 1 < pT < 2 GeV
  resolutionTable[5][2][1][1][6] = 2;
  resolutionTable[5][2][1][1][7] = 2;
  resolutionTable[5][2][1][1][9] = 2;
  resolutionTable[5][2][1][1][13] = 2;
  resolutionTable[5][2][1][1][14] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 10-30 | 2 < pT < 3 GeV
  resolutionTable[5][2][1][2][10] = 2;
  resolutionTable[5][2][1][2][11] = 2;
  resolutionTable[5][2][1][2][13] = 3;
  
  // 0.8 < xj < 1.0 | Centrality: 10-30 | 3 < pT < 4 GeV
  resolutionTable[5][2][1][3][8] = 2;
  resolutionTable[5][2][1][3][10] = 2;
  resolutionTable[5][2][1][3][13] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 10-30 | 4 < pT < 8 GeV
  resolutionTable[5][2][1][4][12] = 2;
  resolutionTable[5][2][1][4][13] = 2;
  resolutionTable[5][2][1][4][14] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 10-30 | 8 < pT < 12 GeV
  resolutionTable[5][2][1][5][9] = 3;
  resolutionTable[5][2][1][5][10] = 3;
  resolutionTable[5][2][1][5][12] = 3;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 0.7 < pT < 1 GeV
  resolutionTable[5][2][2][0][8] = 2;
  resolutionTable[5][2][2][0][12] = 2;
  resolutionTable[5][2][2][0][14] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 1 < pT < 2 GeV
  resolutionTable[5][2][2][1][11] = 2;
  resolutionTable[5][2][2][1][14] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 2 < pT < 3 GeV
  resolutionTable[5][2][2][2][2] = 2;
  resolutionTable[5][2][2][2][8] = 2;
  resolutionTable[5][2][2][2][13] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 3 < pT < 4 GeV
  resolutionTable[5][2][2][3][8] = 2;
  resolutionTable[5][2][2][3][9] = 3;
  resolutionTable[5][2][2][3][10] = 2;
  resolutionTable[5][2][2][3][14] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 4 < pT < 8 GeV
  resolutionTable[5][2][2][4][7] = 3;
  resolutionTable[5][2][2][4][8] = 2;
  resolutionTable[5][2][2][4][11] = 2;
  resolutionTable[5][2][2][4][13] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 8 < pT < 12 GeV
  resolutionTable[5][2][2][5][5] = 2;
  resolutionTable[5][2][2][5][7] = 2;
  resolutionTable[5][2][2][5][9] = 4;
  resolutionTable[5][2][2][5][12] = 3;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 12 < pT < 300 GeV
  resolutionTable[5][2][2][6][6] = 4;
  resolutionTable[5][2][2][6][7] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 0.7 < pT < 1 GeV
  resolutionTable[5][2][3][0][1] = 3;
  resolutionTable[5][2][3][0][4] = 2;
  resolutionTable[5][2][3][0][5] = 2;
  resolutionTable[5][2][3][0][9] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 1 < pT < 2 GeV
  resolutionTable[5][2][3][1][3] = 2;
  resolutionTable[5][2][3][1][4] = 2;
  resolutionTable[5][2][3][1][7] = 2;
  resolutionTable[5][2][3][1][9] = 2;
  resolutionTable[5][2][3][1][11] = 2;
  resolutionTable[5][2][3][1][12] = 2;
  resolutionTable[5][2][3][1][14] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 2 < pT < 3 GeV
  resolutionTable[5][2][3][2][12] = 2;
  resolutionTable[5][2][3][2][14] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 3 < pT < 4 GeV
  resolutionTable[5][2][3][3][9] = 2;
  resolutionTable[5][2][3][3][10] = 2;
  resolutionTable[5][2][3][3][11] = 2;
  resolutionTable[5][2][3][3][14] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 4 < pT < 8 GeV
  resolutionTable[5][2][3][4][10] = 2;
  resolutionTable[5][2][3][4][12] = 2;
  resolutionTable[5][2][3][4][14] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 8 < pT < 12 GeV
  resolutionTable[5][2][3][5][3] = 2;
  resolutionTable[5][2][3][5][5] = 2;
  resolutionTable[5][2][3][5][7] = 2;
  resolutionTable[5][2][3][5][9] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 12 < pT < 300 GeV
  resolutionTable[5][2][3][6][7] = 3;
  
  // Manual configuration ready, return the needed factor
  
  return resolutionTable[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iBin];

}

/*
 *  Systematic uncertainty smoothening for jet resolution
 *
 *  Files used to derive the smoothing factors with compareDijetHistograms.C:
 *   dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_100trig_JECv6_xjBins_wtaAxis_subleadingJffTuning_allCorrections_processed_2020-02-17.root
 *   dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6minus_wtaAxis_allCorrections_finalTuning_onlyFinalResults_processed_2019-12-16.root
 *   dijetPbPb2018_akFlowPuCs4PFJets_noUncOrInc_25eveMix_xjBins_100trig_JECv6plus_wtaAxis_allCorrections_finalTuning_onlyFinalResults_processed_2019-12-29.root
 *
 *  Arguments:
 *   int iJetTrack = Index for the jet track correlation
 *   int iAsymmetry = Index for the dijet momentum balance bin
 *   int iCentrality = Index for the centrality bin
 *   int iTrackPt = Index for the track pT bin
 *   int iBin = Index for the bin in the histogram
 *
 *  return: Factor with which the uncertainty should be divided to iron out fluctuating bins
 *
 */
double getSmoothingJetEnergyScale(int iJetTrack, int iAsymmetry, int iCentrality, int iTrackPt, int iBin){

  const int nJetTrack = 8;
  const int nAsymmetry = 4;
  const int nCentrality = 4;
  const int nTrackPt = 7;
  const int nBin = 15;

  double scaleTable[nJetTrack][nAsymmetry][nCentrality][nTrackPt][nBin];

  // Initialize the resolution table to one
  for(int iJetTrack = 0; iJetTrack < nJetTrack; iJetTrack++){
    for(int iAsymmetry = 0; iAsymmetry < nAsymmetry; iAsymmetry++){
      for(int iCentrality = 0; iCentrality < nCentrality; iCentrality++){
        for(int iTrackPt = 0; iTrackPt < nTrackPt; iTrackPt++){
          for(int iBin = 0; iBin < nBin; iBin++){
            scaleTable[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iBin] = 1;
          }
        }
      }
    }
  }
  
  // Manual configuration here
  
  // ========================
  //    Leading jet shape
  // ========================
  
  // xj integrated | Centrality: 0-10 | 3 < pT < 4 GeV
  scaleTable[2][3][0][3][11] = 2;
  scaleTable[2][3][0][3][12] = 2;
  scaleTable[2][3][0][3][13] = 2;
  scaleTable[2][3][0][3][14] = 4;
  
  // xj integrated | Centrality: 10-30 | 2 < pT < 3 GeV
  scaleTable[2][3][1][2][13] = 2;
  scaleTable[2][3][1][2][14] = 3;
  
  // xj integrated | Centrality: 10-30 | 3 < pT < 4 GeV
  scaleTable[2][3][1][3][1] = 1.6;
  scaleTable[2][3][1][3][2] = 1.6;
  scaleTable[2][3][1][3][3] = 1.6;
  scaleTable[2][3][1][3][4] = 1.6;
  scaleTable[2][3][1][3][5] = 1.3;
  scaleTable[2][3][1][3][6] = 1.3;
  scaleTable[2][3][1][3][7] = 1.3;
  scaleTable[2][3][1][3][8] = 1.3;
  scaleTable[2][3][1][3][9] = 1.3;
  scaleTable[2][3][1][3][11] = 2;
  scaleTable[2][3][1][3][13] = 2;
  scaleTable[2][3][1][3][14] = 2;
  
  // xj integrated | Centrality: 30-50 | 0.7 < pT < 1 GeV
  scaleTable[2][3][2][0][9] = 2;
  scaleTable[2][3][2][0][11] = 2;
  scaleTable[2][3][2][0][12] = 2.5;
  scaleTable[2][3][2][0][13] = 2.5;
  scaleTable[2][3][2][0][14] = 2.5;
  
  // xj integrated | Centrality: 30-50 | 2 < pT < 3 GeV
  scaleTable[2][3][2][2][1] = 2;
  scaleTable[2][3][2][2][2] = 2.4;
  scaleTable[2][3][2][2][3] = 2.4;
  scaleTable[2][3][2][2][4] = 2.4;
  scaleTable[2][3][2][2][5] = 2.1;
  scaleTable[2][3][2][2][6] = 2;
  scaleTable[2][3][2][2][7] = 2;
  scaleTable[2][3][2][2][8] = 2;
  scaleTable[2][3][2][2][9] = 2;
  scaleTable[2][3][2][2][10] = 2;
  scaleTable[2][3][2][2][11] = 2.5;
  scaleTable[2][3][2][2][12] = 3;
  scaleTable[2][3][2][2][13] = 5;
  scaleTable[2][3][2][2][14] = 5;
  
  // xj integrated | Centrality: 30-50 | 3 < pT < 4 GeV
  scaleTable[2][3][2][3][14] = 2;
  
  // xj integrated | Centrality: 50-90 | 0.7 < pT < 1 GeV
  scaleTable[2][3][3][0][13] = 2;
  scaleTable[2][3][3][0][14] = 2;
  
  // xj integrated | Centrality: 50-90 | 1 < pT < 2 GeV
  scaleTable[2][3][3][1][13] = 2;
  scaleTable[2][3][3][1][14] = 3;
  
  // 0.0 < xj < 0.6 | Centrality: 10-30 | 0.7 < pT < 1 GeV
  scaleTable[2][0][1][0][10] = 2;
  scaleTable[2][0][1][0][11] = 2.5;
  scaleTable[2][0][1][0][12] = 3.5;
  scaleTable[2][0][1][0][13] = 3.5;
  scaleTable[2][0][1][0][14] = 3.5;
  
  // 0.0 < xj < 0.6 | Centrality: 10-30 | 2 < pT < 3 GeV
  scaleTable[2][0][1][2][13] = 2;
  scaleTable[2][0][1][2][14] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 30-50 | 0.7 < pT < 1 GeV
  scaleTable[2][0][2][0][10] = 2;
  scaleTable[2][0][2][0][12] = 3;
  scaleTable[2][0][2][0][13] = 2;
  scaleTable[2][0][2][0][14] = 4;
  
  // 0.0 < xj < 0.6 | Centrality: 30-50 | 1 < pT < 2 GeV
  scaleTable[2][0][2][1][14] = 3;
  
  // 0.0 < xj < 0.6 | Centrality: 30-50 | 2 < pT < 3 GeV
  scaleTable[2][0][2][1][13] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 0.7 < pT < 1 GeV
  scaleTable[2][0][3][0][10] = 2;
  scaleTable[2][0][3][0][11] = 3;
  scaleTable[2][0][3][0][12] = 4;
  scaleTable[2][0][3][0][13] = 4;
  scaleTable[2][0][3][0][14] = 10;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 1 < pT < 2 GeV
  scaleTable[2][0][3][1][14] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 2 < pT < 3 GeV
  scaleTable[2][0][3][2][9] = 2;
  scaleTable[2][0][3][2][14] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 0-10 | 0.7 < pT < 1 GeV
  scaleTable[2][1][0][0][12] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 0-10 | 1 < pT < 2 GeV
  scaleTable[2][1][0][1][1] = 2;
  scaleTable[2][1][0][1][2] = 2;
  scaleTable[2][1][0][1][3] = 2;
  scaleTable[2][1][0][1][4] = 2;
  scaleTable[2][1][0][1][5] = 2;
  scaleTable[2][1][0][1][6] = 2;
  scaleTable[2][1][0][1][7] = 2;
  scaleTable[2][1][0][1][8] = 2;
  scaleTable[2][1][0][1][9] = 2;
  scaleTable[2][1][0][1][10] = 2;
  scaleTable[2][1][0][1][11] = 2;
  scaleTable[2][1][0][1][12] = 2;
  scaleTable[2][1][0][1][13] = 3;
  scaleTable[2][1][0][1][14] = 4;
  
  // 0.6 < xj < 0.8 | Centrality: 0-10 | 3 < pT < 4 GeV
  scaleTable[2][1][0][3][14] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 10-30 | 2 < pT < 3 GeV
  scaleTable[2][1][1][2][13] = 2;
  scaleTable[2][1][1][2][14] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 10-30 | 3 < pT < 4 GeV
  scaleTable[2][1][1][3][13] = 3;
  scaleTable[2][1][1][3][14] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 30-50 | 0.7 < pT < 1 GeV
  scaleTable[2][1][2][0][11] = 1.5;
  scaleTable[2][1][2][0][13] = 2;
  scaleTable[2][1][2][0][14] = 3;
  
  // 0.6 < xj < 0.8 | Centrality: 30-50 | 2 < pT < 3 GeV
  scaleTable[2][1][2][2][14] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 0.7 < pT < 1 GeV
  scaleTable[2][1][3][0][12] = 2;
  scaleTable[2][1][3][0][13] = 5;
  scaleTable[2][1][3][0][14] = 6;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 1 < pT < 2 GeV
  scaleTable[2][1][3][1][13] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 2 < pT < 3 GeV
  scaleTable[2][1][3][2][13] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 0.7 < pT < 1 GeV
  scaleTable[2][2][0][0][13] = 2;
  scaleTable[2][2][0][0][14] = 1.5;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 1 < pT < 2 GeV
  scaleTable[2][2][0][1][13] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 3 < pT < 4 GeV
  scaleTable[2][2][0][3][13] = 1.5;
  
  // 0.8 < xj < 1.0 | Centrality: 10-30 | 0.7 < pT < 1 GeV
  scaleTable[2][2][1][0][10] = 1.5;
  scaleTable[2][2][1][0][12] = 1.5;
  scaleTable[2][2][1][0][13] = 1.5;
  scaleTable[2][2][1][0][14] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 10-30 | 2 < pT < 3 GeV
  scaleTable[2][2][1][2][4] = 1.5;
  scaleTable[2][2][1][2][5] = 1.5;
  scaleTable[2][2][1][2][6] = 1.5;
  scaleTable[2][2][1][2][7] = 2;
  scaleTable[2][2][1][2][8] = 2;
  scaleTable[2][2][1][2][9] = 2;
  scaleTable[2][2][1][2][10] = 2;
  scaleTable[2][2][1][2][11] = 2;
  scaleTable[2][2][1][2][12] = 2;
  scaleTable[2][2][1][2][13] = 2;
  scaleTable[2][2][1][2][14] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 10-30 | 3 < pT < 4 GeV
  scaleTable[2][2][1][2][3] = 1.5;
  scaleTable[2][2][1][2][7] = 1.5;
  scaleTable[2][2][1][2][9] = 1.5;
  scaleTable[2][2][1][2][10] = 1.5;
  scaleTable[2][2][1][2][11] = 3;
  scaleTable[2][2][1][2][12] = 1.5;
  scaleTable[2][2][1][2][13] = 2.5;
  scaleTable[2][2][1][2][14] = 3;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 1 < pT < 2 GeV
  scaleTable[2][2][2][1][1] = 2;
  scaleTable[2][2][2][1][2] = 2;
  scaleTable[2][2][2][1][3] = 2;
  scaleTable[2][2][2][1][4] = 2;
  scaleTable[2][2][2][1][5] = 2;
  scaleTable[2][2][2][1][6] = 2;
  scaleTable[2][2][2][1][7] = 2;
  scaleTable[2][2][2][1][8] = 2;
  scaleTable[2][2][2][1][9] = 2;
  scaleTable[2][2][2][1][10] = 2;
  scaleTable[2][2][2][1][11] = 2;
  scaleTable[2][2][2][1][12] = 2;
  scaleTable[2][2][2][1][13] = 3;
  scaleTable[2][2][2][1][14] = 3.5;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 2 < pT < 3 GeV
  scaleTable[2][2][2][2][13] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 3 < pT < 4 GeV
  scaleTable[2][2][2][3][14] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 0.7 < pT < 1 GeV
  scaleTable[2][2][3][0][13] = 1.5;
  scaleTable[2][2][3][0][14] = 1.5;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 1 < pT < 2 GeV
  scaleTable[2][2][3][1][11] = 2;
  scaleTable[2][2][3][1][13] = 2;
  scaleTable[2][2][3][1][14] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 3 < pT < 4 GeV
  scaleTable[2][2][3][3][14] = 1.5;
  
  // ========================
  //   Subleading jet shape
  // ========================
  
  // xj integrated | Centrality: 0-10 | 0.7 < pT < 1 GeV
  scaleTable[5][3][0][0][5] = 0.6;
  scaleTable[5][3][0][0][6] = 0.6;
  scaleTable[5][3][0][0][12] = 1.5;
  scaleTable[5][3][0][0][14] = 2;
  
  // xj integrated | Centrality: 0-10 | 1 < pT < 2 GeV
  scaleTable[5][3][0][1][6] = 1.3;
  scaleTable[5][3][0][1][7] = 1.3;
  scaleTable[5][3][0][1][10] = 1.5;
  scaleTable[5][3][0][1][11] = 1.5;
  scaleTable[5][3][0][1][12] = 2;
  scaleTable[5][3][0][1][14] = 2;
  
  // xj integrated | Centrality: 0-10 | 2 < pT < 3 GeV
  scaleTable[5][3][0][2][12] = 1.5;
  scaleTable[5][3][0][2][13] = 2;
  
  // xj integrated | Centrality: 10-30 | 0.7 < pT < 1 GeV
  scaleTable[5][3][1][0][12] = 1.5;
  scaleTable[5][3][1][0][13] = 2;
  scaleTable[5][3][1][0][14] = 2.5;
  
  // xj integrated | Centrality: 10-30 | 1 < pT < 2 GeV
  scaleTable[5][3][1][1][13] = 1.5;
  scaleTable[5][3][1][1][14] = 1.5;
  
  // xj integrated | Centrality: 10-30 | 3 < pT < 4 GeV
  scaleTable[5][3][1][3][14] = 2;
  
  // xj integrated | Centrality: 30-50 | 0.7 < pT < 1 GeV
  scaleTable[5][3][2][0][13] = 1.5;
  scaleTable[5][3][2][0][14] = 2;
  
  // xj integrated | Centrality: 50-90 | 0.7 < pT < 1 GeV
  scaleTable[5][3][3][0][12] = 1.5;
  scaleTable[5][3][3][0][13] = 1.5;
  scaleTable[5][3][3][0][14] = 2;
  
  // xj integrated | Centrality: 50-90 | 1 < pT < 2 GeV
  scaleTable[5][3][3][1][13] = 2;
  scaleTable[5][3][3][1][14] = 1.5;
  
  // xj integrated | Centrality: 50-90 | 3 < pT < 4 GeV
  scaleTable[5][3][3][3][13] = 1.5;
  scaleTable[5][3][3][3][14] = 3;
  
  // 0.0 < xj < 0.6 | Centrality: 0-10 | 0.7 < pT < 1 GeV
  scaleTable[5][0][0][0][14] = 3;
  
  // 0.0 < xj < 0.6 | Centrality: 0-10 | 1 < pT < 2 GeV
  scaleTable[5][0][0][1][12] = 3;
  scaleTable[5][0][0][1][14] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 0-10 | 2 < pT < 3 GeV
  scaleTable[5][0][0][2][13] = 3;
  scaleTable[5][0][0][2][14] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 10-30 | 2 < pT < 3 GeV
  scaleTable[5][0][1][2][13] = 1.5;
  scaleTable[5][0][1][2][14] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 30-50 | 1 < pT < 2 GeV
  scaleTable[5][0][2][1][14] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 30-50 | 2 < pT < 3 GeV
  scaleTable[5][0][2][2][9] = 2;
  scaleTable[5][0][2][2][14] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 0.7 < pT < 1 GeV
  scaleTable[5][0][3][0][13] = 2;
  scaleTable[5][0][3][0][14] = 2;
  
  // 0.0 < xj < 0.6 | Centrality: 50-90 | 2 < pT < 3 GeV
  scaleTable[5][0][3][2][6] = 1.5;
  scaleTable[5][0][3][2][10] = 1.5;
  scaleTable[5][0][3][2][13] = 1.5;
  scaleTable[5][0][3][2][14] = 1.5;
  
  // 0.6 < xj < 0.8 | Centrality: 0-10 | 0.7 < pT < 1 GeV
  scaleTable[5][1][0][0][1] = 2;
  scaleTable[5][1][0][0][5] = 2;
  scaleTable[5][1][0][0][7] = 2;
  scaleTable[5][1][0][0][8] = 2;
  scaleTable[5][1][0][0][10] = 2;
  scaleTable[5][1][0][0][11] = 2;
  scaleTable[5][1][0][0][12] = 2;
  scaleTable[5][1][0][0][13] = 2;
  scaleTable[5][1][0][0][14] = 4;
  
  // 0.6 < xj < 0.8 | Centrality: 0-10 | 1 < pT < 2 GeV
  scaleTable[5][1][0][1][3] = 1.5;
  scaleTable[5][1][0][1][4] = 1.5;
  scaleTable[5][1][0][1][5] = 1.5;
  scaleTable[5][1][0][1][6] = 1.5;
  scaleTable[5][1][0][1][7] = 1.5;
  scaleTable[5][1][0][1][8] = 1.5;
  scaleTable[5][1][0][1][9] = 1.5;
  scaleTable[5][1][0][1][10] = 2;
  scaleTable[5][1][0][1][11] = 2;
  scaleTable[5][1][0][1][12] = 2;
  scaleTable[5][1][0][1][13] = 2;
  scaleTable[5][1][0][1][14] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 0-10 | 2 < pT < 3 GeV
  scaleTable[5][1][0][1][13] = 3;
  scaleTable[5][1][0][1][14] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 10-30 | 1 < pT < 2 GeV
  scaleTable[5][1][1][1][12] = 2;
  scaleTable[5][1][1][1][13] = 2;
  scaleTable[5][1][1][1][14] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 30-50 | 0.7 < pT < 1 GeV
  scaleTable[5][1][2][0][12] = 1.5;
  scaleTable[5][1][2][0][13] = 2;
  scaleTable[5][1][2][0][14] = 1.5;
  
  // 0.6 < xj < 0.8 | Centrality: 30-50 | 2 < pT < 3 GeV
  scaleTable[5][1][2][2][9] = 2;
  scaleTable[5][1][2][2][12] = 2;
  scaleTable[5][1][2][2][14] = 3;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 0.7 < pT < 1 GeV
  scaleTable[5][1][3][0][11] = 2;
  scaleTable[5][1][3][0][12] = 3;
  scaleTable[5][1][3][0][13] = 3;
  scaleTable[5][1][3][0][14] = 4;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 1 < pT < 2 GeV
  scaleTable[5][1][3][1][13] = 2;
  scaleTable[5][1][3][1][14] = 2;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 2 < pT < 3 GeV
  scaleTable[5][1][3][2][12] = 2;
  scaleTable[5][1][3][2][13] = 2.5;
  scaleTable[5][1][3][2][14] = 3;
  
  // 0.6 < xj < 0.8 | Centrality: 50-90 | 3 < pT < 4 GeV
  scaleTable[5][1][3][3][14] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 0.7 < pT < 1 GeV
  scaleTable[5][2][0][0][10] = 1.5;
  scaleTable[5][2][0][0][13] = 1.5;
  scaleTable[5][2][0][0][14] = 1.5;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 1 < pT < 2 GeV
  scaleTable[5][2][0][1][2] = 1.5;
  scaleTable[5][2][0][1][3] = 1.5;
  scaleTable[5][2][0][1][4] = 1.5;
  scaleTable[5][2][0][1][5] = 1.5;
  scaleTable[5][2][0][1][6] = 1.5;
  scaleTable[5][2][0][1][7] = 1.5;
  scaleTable[5][2][0][1][8] = 1.5;
  scaleTable[5][2][0][1][9] = 1.5;
  scaleTable[5][2][0][1][10] = 1.5;
  scaleTable[5][2][0][1][11] = 1.5;
  
  // 0.8 < xj < 1.0 | Centrality: 0-10 | 3 < pT < 4 GeV
  scaleTable[5][2][0][3][12] = 1.5;
  scaleTable[5][2][0][3][14] = 1.5;
  
  // 0.8 < xj < 1.0 | Centrality: 10-30 | 0.7 < pT < 1 GeV
  scaleTable[5][2][1][0][2] = 3;
  
  // 0.8 < xj < 1.0 | Centrality: 10-30 | 1 < pT < 2 GeV
  scaleTable[5][2][1][1][14] = 2;
  
  // 0.8 < xj < 1.0 | Centrality: 30-50 | 2 < pT < 3 GeV
  scaleTable[5][2][2][2][12] = 1.5;
  scaleTable[5][2][2][2][14] = 3;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 0.7 < pT < 1 GeV
  scaleTable[5][2][3][0][10] = 1.5;
  scaleTable[5][2][3][0][12] = 2;
  scaleTable[5][2][3][0][14] = 4;
  
  // 0.8 < xj < 1.0 | Centrality: 50-90 | 3 < pT < 4 GeV
  scaleTable[5][2][3][3][14] = 2;
  
  return scaleTable[iJetTrack][iAsymmetry][iCentrality][iTrackPt][iBin];

}
