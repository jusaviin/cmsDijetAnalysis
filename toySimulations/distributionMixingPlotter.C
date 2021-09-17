#include "../plotting/JDrawer.h"

/*
 * Plotter for figures produced by the toy simulation featuring a hole in the acceptance for jets and tracks
 * or only for the tracks.
 */
void distributionMixingPlotter(const char *inputFileName = "mixingFromDistribution10M.root"){
  
  // Open the data file and read the jet phi from same events
  TFile *inputFile = TFile::Open(inputFileName);
  
  // Figure saving options
  bool saveFigures = false;         // Flag to determine whather or not save the figures

  // Read the two-dimensional histograms from the input file
  TH2D *jetTrackDeltaPhiDeltaEta = (TH2D*) inputFile->Get("jetDeltaPhiDeltaEta");
  TH2D *jetTrackDeltaPhiDeltaEtaSame = (TH2D*) inputFile->Get("jetDeltaPhiDeltaEtaSame");
  
  // Project deltaEta and deltaPhi distributions out of the jet-track correlation
  TH1D *jetTrackDeltaPhi = jetTrackDeltaPhiDeltaEta->ProjectionX("jetTrackDeltaPhi");
  TH1D *jetTrackDeltaEta = jetTrackDeltaPhiDeltaEta->ProjectionY("jetTrackDeltaEta");
  TH1D *jetTrackDeltaPhiSame = jetTrackDeltaPhiDeltaEtaSame->ProjectionX("jetTrackDeltaPhiSame");
  TH1D *jetTrackDeltaEtaSame = jetTrackDeltaPhiDeltaEtaSame->ProjectionY("jetTrackDeltaEtaSame");
  
  // Project deltaPhi distributions in deltaEta bins
  const int nEtaBins = 10;
  const int nPhiBins = 11;
  double etaBinBorders[] = {-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0};
  double phiBinBorders[] = {-1,-0.8,-0.6,-0.41,-0.22,-0.04,0.04,0.22,0.41,0.6,1.2,1.8};
  int etaBin1, etaBin2, phiBin1, phiBin2;
  TH1D *jetTrackDeltaPhiEtaBin[nEtaBins];
  TH1D *jetTrackDeltaPhiSameEtaBin[nEtaBins];
  TH1D *jetTrackDeltaEtaPhiBin[nPhiBins];
  TH1D *jetTrackDeltaEtaSamePhiBin[nPhiBins];
  TH1D *jetTrackDeltaEtaDifference[nPhiBins-1];
  char namer[100];
  
  for(int iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++){
    
    // Find the defined eta and phi bins from the 2D-histogram
    etaBin1 = jetTrackDeltaPhiDeltaEta->GetYaxis()->FindBin(etaBinBorders[iEtaBin]);
    etaBin2 = jetTrackDeltaPhiDeltaEta->GetYaxis()->FindBin(etaBinBorders[iEtaBin+1]);
    
    // Project the deltaPhi in the given deltaEta window
    sprintf(namer,"jetTrackDeltaPhiEtaBin%d",iEtaBin);
    jetTrackDeltaPhiEtaBin[iEtaBin] = jetTrackDeltaPhiDeltaEta->ProjectionX(namer,etaBin1,etaBin2);
    
    sprintf(namer,"jetTrackDeltaPhiSameEtaBin%d",iEtaBin);
    jetTrackDeltaPhiSameEtaBin[iEtaBin] = jetTrackDeltaPhiDeltaEtaSame->ProjectionX(namer,etaBin1,etaBin2);
    
    
  }
  
  for(int iPhiBin = 0; iPhiBin < nPhiBins; iPhiBin++){
    
    // Find the defined eta and phi bins from the 2D-histogram
    phiBin1 = jetTrackDeltaPhiDeltaEta->GetXaxis()->FindBin(phiBinBorders[iPhiBin]);
    phiBin2 = jetTrackDeltaPhiDeltaEta->GetXaxis()->FindBin(phiBinBorders[iPhiBin+1]);
    
    // Project the deltaEta in the given deltaPhi window
    sprintf(namer,"jetTrackDeltaEtaPhiBin%d",iPhiBin);
    jetTrackDeltaEtaPhiBin[iPhiBin] = jetTrackDeltaPhiDeltaEta->ProjectionY(namer,phiBin1,phiBin2);
    
    // Project the deltaEta in the given deltaPhi window
    sprintf(namer,"jetTrackDeltaEtaSamePhiBin%d",iPhiBin);
    jetTrackDeltaEtaSamePhiBin[iPhiBin] = jetTrackDeltaPhiDeltaEtaSame->ProjectionY(namer,phiBin1,phiBin2);
    
    // Subtrack deltaEta histogram with large deltaPhi from the other deltaPhi slices
    if(iPhiBin > 0){
      jetTrackDeltaEtaDifference[iPhiBin-1] = (TH1D*) jetTrackDeltaEtaPhiBin[iPhiBin]->Clone();
      jetTrackDeltaEtaDifference[iPhiBin-1]->Add(jetTrackDeltaEtaPhiBin[0],-1);
    }
    
  }
  
  // Draw all the distributions
  JDrawer *drawer = new JDrawer();
  
  // Change the right margin better suited for 2D-drawing
  drawer->SetRightMargin(0.1);
  
  // Jet-track correlation including hole for jets and tracks
  drawer->DrawHistogram(jetTrackDeltaPhiDeltaEta,"#Delta#phi","#Delta#eta","Jet-track correlation, mixed","surf1");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/distributionMixingJetTrackDeltaPhiDeltaEta.pdf");
  
  // Jet-track correlation including hole for tracks but no hole for jets
  drawer->DrawHistogram(jetTrackDeltaPhiDeltaEtaSame,"#Delta#phi","#Delta#eta","Jets-track correlation, same","colz");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/distributionMixingJetTrackDeltaPhiDeltaEtaSame.pdf");
  
  
  // Change right margin back to 1D-drawing
  drawer->SetRightMargin(0.06);
  
  // DeltaPhi projection of the jet-track correlation with a hole in the acceptance for jets and tracks
  drawer->DrawHistogram(jetTrackDeltaPhi,"#Delta#phi","#frac{dN}{d#Delta#phi}","Jet-track #Delta#phi, mixing");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationJetTrackDeltaPhi.pdf");
  
  // DeltaEta projection of the jet-track correlation with a hole in the acceptance for jets and tracks
  drawer->DrawHistogram(jetTrackDeltaEta,"#Delta#eta","#frac{dN}{d#Delta#eta}","Jet-track #Delta#eta, mixing");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationJetTrackDeltaEta.pdf");
  
  // DeltaPhi projection of the jet-track correlation with a hole in the acceptance only for tracks
  drawer->DrawHistogram(jetTrackDeltaPhiSame,"#Delta#phi","#frac{dN}{d#Delta#phi}","Jet-track #Delta#phi, same");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationJetTrackDeltaPhiNoJetHole.pdf");
  
  // DeltaEta projection of the jet-track correlation with a hole in the acceptance only for tracks
  drawer->DrawHistogram(jetTrackDeltaEtaSame,"#Delta#eta","#frac{dN}{d#Delta#eta}","Jet-track #Delta#eta, same");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationJetTrackDeltaEtaNoJetHole.pdf");
  
  // DeltaPhi projections in DeltaEta bins
  for(int iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++){
    
    // Hole for tracks and jets
    sprintf(namer,"Hole for track and jets, %.2f < #Delta#eta < %.2f",etaBinBorders[iEtaBin],etaBinBorders[iEtaBin+1]);
    drawer->DrawHistogram(jetTrackDeltaPhiEtaBin[iEtaBin],"#Delta#phi","#frac{dN}{d#Delta#phi}",namer);
    
    if(saveFigures){
      sprintf(namer,"figures/toySimulationJetTrackDeltaPhiEtaBin%d.pdf",iEtaBin);
      gPad->GetCanvas()->SaveAs(namer);
    }
    
    // Hole only for tracks
    /*
    sprintf(namer,"Hole only for tracks, %.2f < #Delta#eta < %.2f",etaBinBorders[iEtaBin],etaBinBorders[iEtaBin+1]);
    drawer->DrawHistogram(jetTrackDeltaPhiNoJetHoleEtaBin[iEtaBin],"#Delta#phi","#frac{dN}{d#Delta#phi}",namer);
    
    if(saveFigures){
      sprintf(namer,"figures/toySimulationJetTrackDeltaPhiNoJetHoleEtaBin%d.pdf",iEtaBin);
      gPad->GetCanvas()->SaveAs(namer);
    }
    */
    
  } // DeltaPhi projections in DeltaEta bins
  
  // DeltaEta projections in DeltaPhi bins
  for(int iPhiBin = 0; iPhiBin < nPhiBins; iPhiBin++){
    
    // Hole for tracks and jets
    sprintf(namer,"Mixing, %.2f < #Delta#phi < %.2f",phiBinBorders[iPhiBin],phiBinBorders[iPhiBin+1]);
    drawer->DrawHistogram(jetTrackDeltaEtaPhiBin[iPhiBin],"#Delta#eta","#frac{dN}{d#Delta#eta}",namer);
    
    if(saveFigures){
      sprintf(namer,"figures/toySimulationJetTrackDeltaEtaPhiBin%d.pdf",iPhiBin);
      gPad->GetCanvas()->SaveAs(namer);
    }
    
    // Hole for tracks and jets
    sprintf(namer,"Same, %.2f < #Delta#phi < %.2f",phiBinBorders[iPhiBin],phiBinBorders[iPhiBin+1]);
    drawer->DrawHistogram(jetTrackDeltaEtaSamePhiBin[iPhiBin],"#Delta#eta","#frac{dN}{d#Delta#eta}",namer);
    
    if(saveFigures){
      sprintf(namer,"figures/toySimulationJetTrackDeltaEtaPhiBin%d.pdf",iPhiBin);
      gPad->GetCanvas()->SaveAs(namer);
    }
    
    if(iPhiBin < nPhiBins-1){
      
      // Hole for tracks and jets
      sprintf(namer,"#Delta#eta difference, (%.2f < #Delta#phi < %.2f) - (%.2f < #Delta#phi < %.2f)",phiBinBorders[iPhiBin+1],phiBinBorders[iPhiBin+2],phiBinBorders[0],phiBinBorders[1]);
      drawer->DrawHistogram(jetTrackDeltaEtaDifference[iPhiBin],"#Delta#eta","#frac{dN}{d#Delta#eta}",namer);
      
      if(saveFigures){
        sprintf(namer,"figures/toySimulationJetTrackDeltaEtaDifference%d.pdf",iPhiBin);
        gPad->GetCanvas()->SaveAs(namer);
      }
      
    }
  } // DeltaEta projections in DeltaPhi bins
}
