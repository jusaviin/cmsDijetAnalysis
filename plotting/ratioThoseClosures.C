#include "JDrawer.h"

void ratioThoseClosures(){

  // Define file names
  const int nFiles = 4;
  TString fileNames[] = {"closureNominal.root", "closureSmear10.root", "closureSmear20.root", "closureSmear30.root"};
  
  // Read the files
  TFile* closureFile[nFiles];
  for(int iFile = 0; iFile < nFiles; iFile++){
    closureFile[iFile] = TFile::Open(fileNames[iFile]);
  }
  
  // Read the histograms from the files
  const int nCentralityBins = 4;
  TH1D* closureSigmaHistogram[nFiles][nCentralityBins];
  TH1D* closureSigmaRatio[nFiles-1][nCentralityBins];
  
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      closureSigmaHistogram[iFile][iCentrality] = (TH1D*) closureFile[iFile]->Get(Form("jetPtClosureSigma_Corr0_Cent%d_xj3_Part2", iCentrality));
      
      // Calculate the ratios
      if(iFile > 0){
        closureSigmaRatio[iFile-1][iCentrality] = (TH1D*) closureSigmaHistogram[iFile][iCentrality]->Clone(Form("ratio%d%d",iFile,iCentrality));
        closureSigmaRatio[iFile-1][iCentrality]->Divide(closureSigmaHistogram[0][iCentrality]);
      }
    }
  }
  
  // Draw the ratios
  JDrawer *drawer = new JDrawer();
  for(int iFile = 0; iFile < nFiles-1; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      closureSigmaRatio[iFile][iCentrality]->GetYaxis()->SetRangeUser(0.6,1.4);
      drawer->DrawHistogram(closureSigmaRatio[iFile][iCentrality], "Gen pT", "Ratio", Form("iFile: %d, Cent: %d", iFile, iCentrality));
    }
  }
  
}
