#include "DijetMethods.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"

/*
 * Convert Xiao's analysis file to a format that can be used with this analysis code
 */
void convertXiaoToThisAnalysis(){
  
  // Open the file to be converted
  TString directoryName = "data/xiao/";
  TString fileName = "bjtc_data_newjes_36p";
  TFile *xiaoFile = TFile::Open(directoryName+fileName+".root");
  
  // We need to fake JCard information in order to use Xiao's data in this analysis framework
  TString fakeCardFileName = "data/dijetPbPb2018_akFlowPuCs4PfJets_jet80trigger_onlyInclusive_JECv5b_processed_2019-09-10.root";
  TFile* fakeCardFile = TFile::Open(fakeCardFileName);
  DijetCard *fakeCard = new DijetCard(fakeCardFile);
  
  // Binning information
  const int nCentralityBins = 3;
  const int nTrackPtBins = 7;
  
  // Histograms to be transferred
  TH1D *jetPt[nCentralityBins];
  TH1D *jetPhi[nCentralityBins];
  TH1D *jetEta[nCentralityBins];
  TH2D *jetTrackSameEvent[nCentralityBins][nTrackPtBins];
  TH2D *jetTrackMixedEvent[nCentralityBins][nTrackPtBins];
  
  // The two dimensional histogram need to be rotated to be used with this analysis
  TH2D *jetTrackSameEventRotated[nCentralityBins][nTrackPtBins];
  TH2D *jetTrackMixedEventRotated[nCentralityBins][nTrackPtBins];
  
  // Read the histograms from the file
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    jetPt[iCentrality] = (TH1D*) xiaoFile->Get(Form("jetQASets/incl_RecoLevel_pt_C%d", iCentrality));
    jetPhi[iCentrality] = (TH1D*) xiaoFile->Get(Form("jetQASets/incl_RecoLevel_phi_C%d", iCentrality));
    jetEta[iCentrality] = (TH1D*) xiaoFile->Get(Form("jetQASets/incl_RecoLevel_eta_C%d", iCentrality));
    
    jetTrackSameEvent[iCentrality][0] = (TH2D*) xiaoFile->Get(Form("incl_RecoJet_RecoTrk/incl_RecoJet_RecoTrk_P0_C%d", iCentrality));
    jetTrackSameEvent[iCentrality][0]->SetName(Form("incl_RecoJet_RecoTrk_P0_C%d_underflow", iCentrality));
    
    jetTrackMixedEvent[iCentrality][0] = (TH2D*) xiaoFile->Get(Form("incl_RecoJet_RecoTrk/incl_RecoJet_RecoTrk_mixing_P0_C%d", iCentrality));
    jetTrackMixedEvent[iCentrality][0]->SetName(Form("incl_RecoJet_RecoTrk_mixing_P0_C%d_underflow", iCentrality));
    
    for(int iTrackPt = 1; iTrackPt < nTrackPtBins; iTrackPt++){
      jetTrackSameEvent[iCentrality][iTrackPt] = (TH2D*) xiaoFile->Get(Form("incl_RecoJet_RecoTrk/incl_RecoJet_RecoTrk_P%d_C%d", iTrackPt-1, iCentrality));
      
      jetTrackMixedEvent[iCentrality][iTrackPt] = (TH2D*) xiaoFile->Get(Form("incl_RecoJet_RecoTrk/incl_RecoJet_RecoTrk_mixing_P%d_C%d", iTrackPt-1, iCentrality));
    } // Track pT loop
    
  } // Centrality loop
  
  // Same and mixed event histogram need to be rotated, because of different convention of deltaPhi and deltaEta axes
  DijetMethods *methods = new DijetMethods();
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      
      jetTrackSameEventRotated[iCentrality][iTrackPt] = methods->RotateHistogram(jetTrackSameEvent[iCentrality][iTrackPt]);
      jetTrackMixedEventRotated[iCentrality][iTrackPt] = methods->RotateHistogram(jetTrackMixedEvent[iCentrality][iTrackPt]);
      
    } // Track pT loop
  } // Centrality loop
  
  // Save the histogram to a file in a format that is readable to this analysis
  
  // Create the output file
  TFile *outputFile = new TFile(directoryName+fileName+"_preprocessed.root","UPDATE");
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    // Write the jet histogram for leading, subleading and any jets for compatibility reasons
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory("leadingJet")) gDirectory->mkdir("leadingJet");
    gDirectory->cd("leadingJet");
    
    jetPt[iCentrality]->Write(Form("leadingJetPt_C%d", iCentrality),TObject::kOverwrite);
    
    jetPhi[iCentrality]->Write(Form("leadingJetPhi_C%d", iCentrality),TObject::kOverwrite);
    
    jetEta[iCentrality]->Write(Form("leadingJetEta_C%d", iCentrality),TObject::kOverwrite);
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory("subleadingJet")) gDirectory->mkdir("subleadingJet");
    gDirectory->cd("subleadingJet");
    
    jetPt[iCentrality]->Write(Form("subleadingJetPt_C%d", iCentrality),TObject::kOverwrite);
    
    jetPhi[iCentrality]->Write(Form("subleadingJetPhi_C%d", iCentrality),TObject::kOverwrite);
    
    jetEta[iCentrality]->Write(Form("subleadingJetEta_C%d", iCentrality),TObject::kOverwrite);
    
    // Return back to main directory
    gDirectory->cd("../");
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory("anyJet")) gDirectory->mkdir("anyJet");
    gDirectory->cd("anyJet");
    
    jetPt[iCentrality]->Write(Form("anyJetPt_C%dA0", iCentrality),TObject::kOverwrite);
    jetPt[iCentrality]->Write(Form("anyJetPt_C%d", iCentrality),TObject::kOverwrite);
    
    jetPhi[iCentrality]->Write(Form("anyJetPhi_C%dA0", iCentrality),TObject::kOverwrite);
    jetPhi[iCentrality]->Write(Form("anyJetPhi_C%d", iCentrality),TObject::kOverwrite);
    
    jetEta[iCentrality]->Write(Form("anyJetEta_C%dA0", iCentrality),TObject::kOverwrite);
    jetEta[iCentrality]->Write(Form("anyJetEta_C%d", iCentrality),TObject::kOverwrite);
    
    // Return back to main directory
    gDirectory->cd("../");
    
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      
      // Write the jet track correlation histograms
      
      // Create a directory for the histograms if it does not already exist
      if(!gDirectory->GetDirectory("trackJetInclusive")) gDirectory->mkdir("trackJetInclusive");
      gDirectory->cd("trackJetInclusive");
      
      jetTrackSameEventRotated[iCentrality][iTrackPt]->Write(Form("trackJetInclusiveDeltaEtaDeltaPhi_SameEvent_C%dT%d", iCentrality, iTrackPt),TObject::kOverwrite);
      jetTrackMixedEventRotated[iCentrality][iTrackPt]->Write(Form("trackJetInclusiveDeltaEtaDeltaPhi_MixedEvent_C%dT%d", iCentrality, iTrackPt),TObject::kOverwrite);
      
      // Return back to main directory
      gDirectory->cd("../");
      
      
    } // Track pT loop
  } // Centrality loop
  
  // Write the card to the output file if it is not already written
  if(!gDirectory->GetDirectory("JCard")) fakeCard->Write(outputFile);
  
  // Close the file after everything is written
  outputFile->Close();
}
