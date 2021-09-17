#include "JDrawer.h"

/*
 * Extract a histogram with given restrictions on other axes in THnSparse
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH1D
 *   int nAxes = Number of axes that are restained for the projection
 *   int *axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
 */
TH1D* findHistogram(TFile *inputFile, const char *name, int xAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex){
  
  // Read the histogram with the given name from the file
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  
  // If cannot find histogram, inform that it could not be found and return null
  if(histogramArray == nullptr){
    cout << "Could not find " << name << ". Skipping loading this histogram." << endl;
    return NULL;
  }
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for each histogram that is read from the file
  char newName[200];
  sprintf(newName,"%s",histogramArray->GetName());
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    sprintf(newName,"%s_%d=%d-%d",newName,axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]);
  }

  // Project out the histogram and give it the created unique name
  TH1D *projectedHistogram = NULL;
  
  // Check that we are not trying to project a non-existing axis
  if(xAxis < histogramArray->GetNdimensions()){
    projectedHistogram = (TH1D*) histogramArray->Projection(xAxis);
    projectedHistogram->SetName(newName);
  
    // Apply bin width normalization to the projected histogram
    projectedHistogram->Scale(1.0,"width");
  }
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
 *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
 *   int nAxes = Number of axes that are restained for the projection
 *   int *axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
 */
TH2D* findHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex){
  
  // Read the histogram with the given name from the file
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  
  // If cannot find histogram, inform that it could not be found and return null
  if(histogramArray == nullptr){
    cout << "Could not find " << name << ". Skipping loading this histogram." << endl;
    return NULL;
  }
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for eeach histogram that is read from the file
  char newName[200];
  sprintf(newName,"%s",histogramArray->GetName());
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    sprintf(newName,"%s_%d=%d-%d",newName,axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]);
  }
  
  // Project out the histogram and give it the created unique name
  TH2D *projectedHistogram = (TH2D*) histogramArray->Projection(yAxis,xAxis);
  projectedHistogram->SetName(newName);
  
  // Apply bin width normalization to the projected histogram
  projectedHistogram->Scale(1.0,"width");
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Explore the jets which are badly matching between reco and gen xj
 *
 * Project the things out of THnSparses here instead of using histogram manager so that
 * several different projections can be done quickly and easily.
 * Also as these are non-standard histograms, no need to include them in histogram manager.
 */
void xjSmearingExplorer(){

  // Open the input file
  TString fileName = "data/PbPbMC2018_RecoReco_akFlowJet_onlyJets_5pCentShift_recoDijet_xjMatrixWithMap_2020-04-16.root";
  TFile *inputFile = TFile::Open(fileName);

  // Define the axes that need restrictions

  /*
   * THnSparse for dijets:
   *
   *   Histogram name        Axis index       Content of axis
   * ----------------------------------------------------------
   *        dijet              Axis 0         Leading jet pT
   *        dijet              Axis 1        Subleading jet pT
   *        dijet              Axis 2         Dijet deltaPhi
   *        dijet              Axis 3     Dijet momentum balance xj
   *        dijet              Axis 4           Centrality
   *        dijet              Axis 5     xj from matched (reco/gen) dijet
   *        dijet              Axis 6         Leading jet phi
   *        dijet              Axis 7         Leading jet eta
   *        dijet              Axis 8        Subleading jet phi
   *        dijet              Axis 9        Subleading jet eta
   */

  // For test at first, get the xj map in the most central bin
  int axisIndices[3] = {0};
  int lowLimits[3] = {0};
  int highLimits[3] = {0};

  axisIndices[0] = 4;
  lowLimits[0] = 2;
  highLimits[0] = 2;

  TH1D *leadingPt = findHistogram(inputFile, "dijet", 0, 1, axisIndices, lowLimits, highLimits);
  TH1D *leadingPhi = findHistogram(inputFile, "dijet", 6, 1, axisIndices, lowLimits, highLimits);
  TH1D *subleadingPt = findHistogram(inputFile, "dijet", 1, 1, axisIndices, lowLimits, highLimits);
  TH1D *deltaPhi = findHistogram(inputFile, "dijet", 2, 1, axisIndices, lowLimits, highLimits);
  TH2D *centralSmear = findHistogram2D(inputFile, "dijet", 3, 5, 1, axisIndices, lowLimits, highLimits);

  // Next, add cuts on the reconstructed and generator level xj
  double epsilon = 0.001;

  // =========
  // xj matrix requiring leading jet pT 120-150 GeV
  // =========

  axisIndices[1] = 0;
  lowLimits[1] = leadingPt->GetXaxis()->FindBin(120 + epsilon);
  highLimits[1] = leadingPt->GetXaxis()->FindBin(150 - epsilon);

  TH2D *lowLeadingPtCentralXjMatrix = findHistogram2D(inputFile, "dijet", 3, 5, 2, axisIndices, lowLimits, highLimits);

  // =========
  // xj matrix requiring leading jet pT > 150 GeV
  // =========

  axisIndices[1] = 0;
  lowLimits[1] = leadingPt->GetXaxis()->FindBin(150 + epsilon);
  highLimits[1] = leadingPt->GetXaxis()->FindBin(500 - epsilon);

  TH2D *highLeadingPtCentralXjMatrix = findHistogram2D(inputFile, "dijet", 3, 5, 2, axisIndices, lowLimits, highLimits);

  // =========
  // xj matrix requiring subleading jet pT 50-80 GeV
  // =========

  axisIndices[1] = 1;
  lowLimits[1] = subleadingPt->GetXaxis()->FindBin(50 + epsilon);
  highLimits[1] = subleadingPt->GetXaxis()->FindBin(80 - epsilon);

  TH2D *lowSubleadingPtCentralXjMatrix = findHistogram2D(inputFile, "dijet", 3, 5, 2, axisIndices, lowLimits, highLimits);

  // =========
  // xj matrix requiring subleading jet pT > 80 GeV
  // =========

  axisIndices[1] = 1;
  lowLimits[1] = subleadingPt->GetXaxis()->FindBin(80 + epsilon);
  highLimits[1] = subleadingPt->GetXaxis()->FindBin(500 - epsilon);

  TH2D *highSubleadingPtCentralXjMatrix = findHistogram2D(inputFile, "dijet", 3, 5, 2, axisIndices, lowLimits, highLimits);

  // =========
  // xj matrix requiring dijet deltaPhi < 2.95
  // =========

  axisIndices[1] = 2;
  lowLimits[1] = deltaPhi->GetXaxis()->FindBin(0 + epsilon);
  highLimits[1] = deltaPhi->GetXaxis()->FindBin(2.93 - epsilon);

  TH2D *lowDeltaPhiCentralXjMatrix = findHistogram2D(inputFile, "dijet", 3, 5, 2, axisIndices, lowLimits, highLimits);

  // =========
  // xj matrix requiring dijet deltaPhi > 2.95
  // =========

  axisIndices[1] = 2;
  lowLimits[1] = deltaPhi->GetXaxis()->FindBin(1.98 + epsilon);
  highLimits[1] = deltaPhi->GetXaxis()->FindBin(3.12 - epsilon);

  TH2D *highDeltaPhiCentralXjMatrix = findHistogram2D(inputFile, "dijet", 3, 5, 2, axisIndices, lowLimits, highLimits);

  // =========
  // xj matrix requiring -2 < leading jet phi < 3
  // =========

  axisIndices[1] = 6;
  lowLimits[1] = leadingPhi->GetXaxis()->FindBin(-2 + epsilon);
  highLimits[1] = leadingPhi->GetXaxis()->FindBin(3 - epsilon);

  TH2D *goodLeadingPhiCentralXjMatrix = findHistogram2D(inputFile, "dijet", 3, 5, 2, axisIndices, lowLimits, highLimits);

  // =========
  // xj matrix requiring -3 < leading jet phi < -2
  // =========

  axisIndices[1] = 6;
  lowLimits[1] = leadingPhi->GetXaxis()->FindBin(-3 + epsilon);
  highLimits[1] = leadingPhi->GetXaxis()->FindBin(-2 - epsilon);

  TH2D *badLeadingPhiCentralXjMatrix = findHistogram2D(inputFile, "dijet", 3, 5, 2, axisIndices, lowLimits, highLimits);

  // =========
  // Leading jet eta-phi for low reconstructed and high generator level xj
  // =========

  // Require reconstructed xj in the bin 0.0 < xj < 0.6
  axisIndices[1] = 3;
  lowLimits[1] = centralSmear->GetXaxis()->FindBin(0.0 + epsilon);
  highLimits[1] = centralSmear->GetXaxis()->FindBin(0.6 - epsilon);

  // Require generator level xj in the interval 0.7 < xj < 1.0
  axisIndices[2] = 5;
  lowLimits[2] = centralSmear->GetYaxis()->FindBin(0.7 + epsilon);
  highLimits[2] = centralSmear->GetYaxis()->FindBin(1.0 - epsilon);

  TH2D *leadingMap = findHistogram2D(inputFile, "dijet", 6, 7, 3, axisIndices, lowLimits, highLimits);

  // =========
  // Leading jet eta-phi for low reconstructed and low generator level xj
  // =========

  // Require reconstructed xj in the bin 0.0 < xj < 0.6
  axisIndices[1] = 3;
  lowLimits[1] = centralSmear->GetXaxis()->FindBin(0.0 + epsilon);
  highLimits[1] = centralSmear->GetXaxis()->FindBin(0.6 - epsilon);

  // Require generator level xj in the interval 0.0 < xj < 0.6
  axisIndices[2] = 5;
  lowLimits[2] = centralSmear->GetYaxis()->FindBin(0.0 + epsilon);
  highLimits[2] = centralSmear->GetYaxis()->FindBin(0.6 - epsilon);
  TH2D *leadingMapMatched = findHistogram2D(inputFile, "dijet", 6, 7, 3, axisIndices, lowLimits, highLimits);

  
  // Can print out mean pT from the histograms
  double lowBinBorder = 0;
  double highBinBorder = 0.6;
  int lowBorder = goodLeadingPhiCentralXjMatrix->GetXaxis()->FindBin(lowBinBorder+0.001);
  int highBorder = goodLeadingPhiCentralXjMatrix->GetXaxis()->FindBin(highBinBorder-0.001);
  
  TH1D *projectionX = goodLeadingPhiCentralXjMatrix->ProjectionX("goodProjectionX");
  TH1D *projectionY = goodLeadingPhiCentralXjMatrix->ProjectionY("goodProjectionY", lowBorder, highBorder);
  
  projectionX->GetXaxis()->SetRangeUser(lowBinBorder, highBinBorder);
  
  cout << endl;
  cout << "Good leading phi central bin" << endl;
  cout << "Bin " << lowBinBorder << " < xj < " <<highBinBorder << endl;
  cout << "Mean reconstructed xj = " << projectionX->GetMean() << " +- " << projectionX->GetMeanError() << endl;
  cout << "Mean generator level xj = " << projectionY->GetMean() << " +- " << projectionY->GetMeanError() << endl;
  cout << endl;
  
  projectionX = badLeadingPhiCentralXjMatrix->ProjectionX("badProjectionX");
  projectionY = badLeadingPhiCentralXjMatrix->ProjectionY("badProjectionY", lowBorder, highBorder);
  
  projectionX->GetXaxis()->SetRangeUser(lowBinBorder, highBinBorder);
  
  cout << endl;
  cout << "Bad leading phi central bin" << endl;
  cout << "Bin " << lowBinBorder << " < xj < " <<highBinBorder << endl;
  cout << "Mean reconstructed xj = " << projectionX->GetMean() << " +- " << projectionX->GetMeanError() << endl;
  cout << "Mean generator level xj = " << projectionY->GetMean() << " +- " << projectionY->GetMeanError() << endl;
  cout << endl;
  
  projectionX = centralSmear->ProjectionX("defaultProjectionX");
  projectionY = centralSmear->ProjectionY("defaultProjectionY", lowBorder, highBorder);
  
  projectionX->GetXaxis()->SetRangeUser(lowBinBorder+epsilon, highBinBorder-epsilon);
  
  cout << endl;
  cout << "Default leading phi central bin" << endl;
  cout << "Bin " << lowBinBorder << " < xj < " <<highBinBorder << endl;
  cout << "Mean reconstructed xj = " << projectionX->GetMean() << " +- " << projectionX->GetMeanError() << endl;
  cout << "Mean generator level xj = " << projectionY->GetMean() << " +- " << projectionY->GetMeanError() << endl;
  cout << endl;
  
  TLine *diagonalLine = new TLine(0,0,1,1);
  diagonalLine->SetLineStyle(2);

  gStyle->SetPalette(kRainBow);
  JDrawer *drawer = new JDrawer();

  // Change the right margin better suited for 2D-drawing
  drawer->SetRightMargin(0.14);
    
  // For xj matrix, use logarithmic scale to better see the structures
  drawer->SetLogZ(true);
    
  // Make the plot square
  drawer->SetCanvasSize(700,600);  // 600 for range 0-1, 1250 for range 0-2

  drawer->DrawHistogram(centralSmear, "Reconstructed x_{j}", "Generator level x_{j}", "C = 0-10", "colz");
  diagonalLine->Draw();

  // xj matrices in different jet pT regions
  drawer->DrawHistogram(lowLeadingPtCentralXjMatrix, "Reconstructed x_{j}", "Generator level x_{j}", "Leading jet pT 120-150 GeV, C = 0-10", "colz");
  diagonalLine->Draw();
  drawer->DrawHistogram(highLeadingPtCentralXjMatrix, "Reconstructed x_{j}", "Generator level x_{j}", "Leading jet pT > 150 GeV, C = 0-10", "colz");
  diagonalLine->Draw();
  drawer->DrawHistogram(lowSubleadingPtCentralXjMatrix, "Reconstructed x_{j}", "Generator level x_{j}", "Subleading jet pT 50-80 GeV, C = 0-10", "colz");
  diagonalLine->Draw();
  drawer->DrawHistogram(highSubleadingPtCentralXjMatrix, "Reconstructed x_{j}", "Generator level x_{j}", "Subleading jet pT > 80 GeV, C = 0-10", "colz");
  diagonalLine->Draw();
  drawer->DrawHistogram(lowDeltaPhiCentralXjMatrix, "Reconstructed x_{j}", "Generator level x_{j}", "#Delta#phi < 2.95, C = 0-10", "colz");
  diagonalLine->Draw();
  drawer->DrawHistogram(highDeltaPhiCentralXjMatrix, "Reconstructed x_{j}", "Generator level x_{j}", "#Delta#phi > 2.95, C = 0-10", "colz");
  diagonalLine->Draw();
  drawer->DrawHistogram(goodLeadingPhiCentralXjMatrix, "Reconstructed x_{j}", "Generator level x_{j}", "-2 < Leading jet #phi < 3, C = 0-10", "colz");
  diagonalLine->Draw();
  drawer->DrawHistogram(badLeadingPhiCentralXjMatrix, "Reconstructed x_{j}", "Generator level x_{j}", "-3 < Leading jet #phi < -2, C = 0-10", "colz");
  diagonalLine->Draw();

  drawer->SetCanvasSize(800,600);  // 600 for range 0-1, 1250 for range 0-2
  drawer->SetTopMargin(0.1);

  // Leading jet eta-phi maps for different xj regions
  drawer->DrawHistogram(leadingMap, "Reco leading jet #phi", "Reco leading jet #eta", "0.0 < x_{j}^{reco} < 0.6, 0.7 < x_{j}^{gen} < 1.0, C = 0-10", "colz");
  drawer->DrawHistogram(leadingMapMatched, "Reco leading jet #phi", "Reco leading jet #eta", "0.0 < x_{j}^{reco} < 0.6, 0.0 < x_{j}^{gen} < 0.6, C = 0-10", "colz");
}
