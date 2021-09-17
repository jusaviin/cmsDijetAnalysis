#include "plotting/JDrawer.h"

/*
 * Plotter for figures produced by the toy simulation featuring a hole in the acceptance for jets and tracks
 * or only for the tracks.
 */
void acceptanceHolePlotter(const char *inputFileName = "toyHoleInAcceptance10M.root"){
  
  // Open the data file and read the jet phi from same events
  TFile *inputFile = TFile::Open(inputFileName);
  
  // Figure saving options
  bool saveFigures = false;         // Flag to determine whather or not save the figures

  // Read the two-dimensional histograms from the input file
  TH2D *jetTrackDeltaPhiDeltaEta = (TH2D*) inputFile->Get("jetDeltaPhiDeltaEta");
  TH2D *jetTrackDeltaPhiDeltaEtaNoJetHole = (TH2D*) inputFile->Get("jetDeltaPhiDeltaEtaNoJetHole");
  TH2D *jetTrackDeltaPhiDeltaEtaNoHole = (TH2D*) inputFile->Get("jetDeltaPhiDeltaEtaNoHole");
  TH2D *jetPhiEta = (TH2D*) inputFile->Get("jetPhiEta");
  TH2D *jetPhiEtaNoHole = (TH2D*) inputFile->Get("jetPhiEtaNoHole");
  TH2D *trackPhiEta = (TH2D*) inputFile->Get("trackPhiEta");
  TH2D *trackPhiEtaNoHole = (TH2D*) inputFile->Get("trackPhiEtaNoHole");
  
  // Project deltaEta and deltaPhi distributions out of the jet-track correlation
  TH1D *jetTrackDeltaPhi = jetTrackDeltaPhiDeltaEta->ProjectionX("jetTrackDeltaPhi");
  TH1D *jetTrackDeltaEta = jetTrackDeltaPhiDeltaEta->ProjectionY("jetTrackDeltaEta");
  TH1D *jetTrackDeltaPhiNoJetHole = jetTrackDeltaPhiDeltaEtaNoJetHole->ProjectionX("jetTrackDeltaPhiNoJetHole");
  TH1D *jetTrackDeltaEtaNoJetHole = jetTrackDeltaPhiDeltaEtaNoJetHole->ProjectionY("jetTrackDeltaEtaNoJetHole");
  TH1D *jetTrackDeltaPhiNoHole = jetTrackDeltaPhiDeltaEtaNoHole->ProjectionX("jetTrackDeltaPhiNoHole");
  TH1D *jetTrackDeltaEtaNoHole = jetTrackDeltaPhiDeltaEtaNoHole->ProjectionY("jetTrackDeltaEtaNoHole");
  
  // Project deltaPhi distributions in deltaEta bins
  const int nEtaBins = 10;
  const int nPhiBins = 11;
  double etaBinBorders[] = {-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0};
  double phiBinBorders[] = {-1,-0.8,-0.6,-0.41,-0.22,-0.04,0.04,0.22,0.41,0.6,0.8,1.0};
  int etaBin1, etaBin2, phiBin1, phiBin2;
  TH1D *jetTrackDeltaPhiEtaBin[nEtaBins];
  TH1D *jetTrackDeltaPhiNoJetHoleEtaBin[nEtaBins];
  TH1D *jetTrackDeltaPhiNoHoleEtaBin[nEtaBins];
  TH1D *jetTrackDeltaEtaPhiBin[nPhiBins];
  TH1D *jetTrackDeltaEtaDifference[nPhiBins-1];
  char namer[100];
  
  for(int iEtaBin = 0; iEtaBin < nEtaBins; iEtaBin++){
    
    // Find the defined eta and phi bins from the 2D-histogram
    etaBin1 = jetTrackDeltaPhiDeltaEta->GetYaxis()->FindBin(etaBinBorders[iEtaBin]);
    etaBin2 = jetTrackDeltaPhiDeltaEta->GetYaxis()->FindBin(etaBinBorders[iEtaBin+1]);
    
    // Project the deltaPhi in the given deltaEta window
    sprintf(namer,"jetTrackDeltaPhiEtaBin%d",iEtaBin);
    jetTrackDeltaPhiEtaBin[iEtaBin] = jetTrackDeltaPhiDeltaEta->ProjectionX(namer,etaBin1,etaBin2);
    
    sprintf(namer,"jetTrackDeltaPhiNoJetHoleEtaBin%d",iEtaBin);
    jetTrackDeltaPhiNoJetHoleEtaBin[iEtaBin] = jetTrackDeltaPhiDeltaEtaNoJetHole->ProjectionX(namer,etaBin1,etaBin2);
    
    sprintf(namer,"jetTrackDeltaPhiNoHoleEtaBin%d",iEtaBin);
    jetTrackDeltaPhiNoHoleEtaBin[iEtaBin] = jetTrackDeltaPhiDeltaEtaNoHole->ProjectionX(namer,etaBin1,etaBin2);
    
  }
  
  for(int iPhiBin = 0; iPhiBin < nPhiBins; iPhiBin++){
    
    // Find the defined eta and phi bins from the 2D-histogram
    phiBin1 = jetTrackDeltaPhiDeltaEta->GetXaxis()->FindBin(phiBinBorders[iPhiBin]);
    phiBin2 = jetTrackDeltaPhiDeltaEta->GetXaxis()->FindBin(phiBinBorders[iPhiBin+1]);
    
    // Project the deltaEta in the given deltaPhi window
    sprintf(namer,"jetTrackDeltaEtaPhiBin%d",iPhiBin);
    jetTrackDeltaEtaPhiBin[iPhiBin] = jetTrackDeltaPhiDeltaEta->ProjectionY(namer,phiBin1,phiBin2);
    
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
  drawer->DrawHistogram(jetTrackDeltaPhiDeltaEta,"#Delta#phi","#Delta#eta","Jet-track correlation, hole for jets and tracks","surf1");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationJetTrackDeltaPhiDeltaEta.pdf");
  
  // Jet-track correlation including hole for tracks but no hole for jets
  drawer->DrawHistogram(jetTrackDeltaPhiDeltaEtaNoJetHole,"#Delta#phi","#Delta#eta","Jets-track correlation, no hole for jets","colz");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationJetTrackDeltaPhiDeltaEtaNoJetHole.pdf");
  
  // Jet-track correlation with no holes in acceptance
  drawer->DrawHistogram(jetTrackDeltaPhiDeltaEtaNoHole,"#Delta#phi","#Delta#eta","Jets-track correlation, no hole for jets or tracks","colz");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationJetTrackDeltaPhiDeltaEtaNoHole.pdf");
  
  // Jet eta-phi distribution with hole in acceptance
  drawer->DrawHistogram(jetPhiEta,"#phi","#eta","Jets","colz");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationJetPhiEta.pdf");
  
  // Jet eta-phi distribution without a hole in acceptance
  drawer->DrawHistogram(jetPhiEtaNoHole,"#phi","#eta","Jets, no hole","colz");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationJetPhiEtaNoHole.pdf");
  
  // Track eta-phi distribution with a hole in the acceptance
  drawer->DrawHistogram(trackPhiEta,"#phi","#eta","Tracks","colz");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationTrackPhiEta.pdf");
  
  // Track eta-phi distribution without a hole in the acceptance
  drawer->DrawHistogram(trackPhiEtaNoHole,"#phi","#eta","Tracks without hole","colz");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationTrackPhiEtaNoHole.pdf");
  
  // Change right margin back to 1D-drawing
  drawer->SetRightMargin(0.06);
  
  // DeltaPhi projection of the jet-track correlation with a hole in the acceptance for jets and tracks
  drawer->DrawHistogram(jetTrackDeltaPhi,"#Delta#phi","#frac{dN}{d#Delta#phi}","Jet-track #Delta#phi, hole for tracks and jets");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationJetTrackDeltaPhi.pdf");
  
  // DeltaEta projection of the jet-track correlation with a hole in the acceptance for jets and tracks
  drawer->DrawHistogram(jetTrackDeltaEta,"#Delta#eta","#frac{dN}{d#Delta#eta}","Jet-track #Delta#eta, hole for tracks and jets");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationJetTrackDeltaEta.pdf");
  
  // DeltaPhi projection of the jet-track correlation with a hole in the acceptance only for tracks
  drawer->DrawHistogram(jetTrackDeltaPhiNoJetHole,"#Delta#phi","#frac{dN}{d#Delta#phi}","Jet-track #Delta#phi, no hole for jets");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationJetTrackDeltaPhiNoJetHole.pdf");
  
  // DeltaEta projection of the jet-track correlation with a hole in the acceptance only for tracks
  drawer->DrawHistogram(jetTrackDeltaEtaNoJetHole,"#Delta#eta","#frac{dN}{d#Delta#eta}","Jet-track #Delta#eta, no hole for jets");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationJetTrackDeltaEtaNoJetHole.pdf");
  
  // DeltaPhi projection of the jet-track correlation with full acceptance
  drawer->DrawHistogram(jetTrackDeltaPhiNoHole,"#Delta#phi","#frac{dN}{d#Delta#phi}","Jet-track #Delta#phi, no holes");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationJetTrackDeltaPhiNoHole.pdf");
  
  // DeltaEta projection of the jet-track correlation with full acceptance
  drawer->DrawHistogram(jetTrackDeltaEtaNoHole,"#Delta#eta","#frac{dN}{d#Delta#eta}","Jet-track #Delta#eta, no holes");
  if(saveFigures) gPad->GetCanvas()->SaveAs("figures/toySimulationJetTrackDeltaEtaNoHole.pdf");
  
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
    sprintf(namer,"Hole for track and jets, %.2f < #Delta#phi < %.2f",phiBinBorders[iPhiBin],phiBinBorders[iPhiBin+1]);
    drawer->DrawHistogram(jetTrackDeltaEtaPhiBin[iPhiBin],"#Delta#eta","#frac{dN}{d#Delta#eta}",namer);
    
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
