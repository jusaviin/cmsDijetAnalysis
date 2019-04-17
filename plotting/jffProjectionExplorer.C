#include "DijetCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetMethods.h"
#include "JDrawer.h"

/*
 * Compare jet shape corrections
 */
void jffProjectionExplorer(){
  
  // Open the files for the JFF corrections
  const int nFiles = 1;
  TFile *standardFile = TFile::Open("newPpTest4.root");
  TFile *rebinFile = TFile::Open("newPbPbTest.root");
  TFile *fitFile = TFile::Open("newPbPbTestSameEvent2.root");
  TFile *files[] = {standardFile,rebinFile,fitFile};
  TFile *inclusiveAnalysisJffPbPb = TFile::Open("data/JffResidual_nominal_Pb_HIN_16_020.root");
  
  // Styling option for different files
  int fileColors[] = {kRed, kBlue, kGreen+3, kMagenta, kCyan};
  const char* fileLegend[] = {"PF Jets","PF jets","sameEvent2","Other method","Another method"};
  
  // Figure saving options
  bool saveFigures = true;         // Flag to determine whather or not save the figures
  
  // Configuration
  const int nCentralityBins = 4;
  const int nTrackPtBins = 6;
  const int nJetTrack = 6;
  double trackPtBorders[] = {0.7,1,2,3,4,8,300};
  double centralityBinBorders[] = {0,10,30,50,100};
  double minXrange[] = {-1.5,-1.5,0};
  
  const int deltaEtaRebin = 4;
  const int deltaPhiRebin = 2;
  
  // Choose correlation types to draw
  bool drawLeading = false;
  bool drawLeadingPtWeighted = false;
  bool drawSubleading = false;
  bool drawSubleadingPtWeighted = false;
  bool drawInclusive = true;
  bool drawInclusivePtWeighted = false;
  
  // Draw old inclusive
  bool drawOldInclusive = false;
  
  // Choose angles to draw
  bool drawDeltaEta = false;
  bool drawDeltaPhi = false;
  bool drawDeltaR = true;
  
  bool drawJetTrack[nJetTrack] = {drawLeading,drawLeadingPtWeighted,drawSubleading,drawSubleadingPtWeighted,drawInclusive,drawInclusivePtWeighted};
  bool drawAngle[] = {drawDeltaEta,drawDeltaPhi,drawDeltaR};
  
  // Jff correction histograms
  // First bin = Leading, pT weighted leading, subeading, pT weighted subleading, inclusive, pT weighted inclusive
  // Second bin = No rebin, rebin2
  TH2D *jffCorrectionDeltaEtaDeltaPhi[nJetTrack][nFiles][nCentralityBins][nTrackPtBins];
  TH1D *jffProjection[3][nJetTrack][nFiles][nCentralityBins][nTrackPtBins]; // First bin deltaEta/deltaPhi/deltaR
  TH1D *jffCorrectionInclusive[nCentralityBins][nTrackPtBins];
  
  // All the variables needed for figure drawing
  JDrawer *drawer = new JDrawer();
  TLegend *legend;
  TH1D* ratio;
  const char *headers[] = {"Track-leading jet","p_{T}w track-leading jet","Track-subleading jet","p_{T}w track-subleading jet","Track-inclusive jet","p_{T}w track-inclusive jet"};
  const char *saveName[] = {"trackLeadingJet","trackLeadingJetPtWeighted","trackSubleadingJet","trackSubleadingJetPtWeighted","trackJetInclusive","trackJetInclusivePtWeighted"};
  const char *deltaLatex[3] = {"#Delta#eta","#Delta#phi","#DeltaR"};
  const char *deltaName[3] = {"DeltaEta","DeltaPhi","DeltaR"};
  
  // Check if we are using pp or PbPb data
  DijetCard *card = new DijetCard(standardFile);
  TString collisionSystem = card->GetDataType();
  bool ppData = (collisionSystem.Contains("pp") || collisionSystem.Contains("localTest"));
  
  // Create a new DijetMethods to easily project eta and phi from the two-dimensional histogram
  DijetMethods *projector = new DijetMethods();
  
  // Read the JFF correction histograms from files
  char namer[200];
  char centralityString[50];
  char compactCentralityString[50];
  int lowBin, highBin;
  for(int iFile = 0; iFile < nFiles; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      if(ppData && (iCentrality > 0)) continue; // No centrality selection for pp
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        
        for(int iJetTrack = 0; iJetTrack < nJetTrack; iJetTrack++){
          if(!drawJetTrack[iJetTrack]) continue;
          
          // Read the corrections for different jet-track correlations
          sprintf(namer,"%sDeltaEtaDeltaPhi/jffCorrection_%sDeltaEtaDeltaPhi_C%dT%d", saveName[iJetTrack], saveName[iJetTrack], iCentrality, iTrackPt);
          jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iFile][iCentrality][iTrackPt] = (TH2D*) files[iFile]->Get(namer);
          cout << "Integral of correction iFile: " << iFile << " Centrality: " << " Track pT: " << " jetTrack: " << iJetTrack << " is: " << jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iFile][iCentrality][iTrackPt]->Integral("width") << endl;
          
          // Project the deltaEta for region -1.5 < deltaPhi < 1.5
          jffProjection[0][iJetTrack][iFile][iCentrality][iTrackPt] = projector->ProjectRegionDeltaEta(jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iFile][iCentrality][iTrackPt],-1.5,1.5,Form("DeltaEtaProjection%d",iFile));
          jffProjection[0][iJetTrack][iFile][iCentrality][iTrackPt]->Rebin(deltaEtaRebin);
          jffProjection[0][iJetTrack][iFile][iCentrality][iTrackPt]->Scale(1.0*projector->GetNBinsProjectedOver()/deltaEtaRebin);
          
          // Project the deltaPhi for region -1.5 < deltaEta < 1.5
          jffProjection[1][iJetTrack][iFile][iCentrality][iTrackPt] = projector->ProjectRegionDeltaPhi(jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iFile][iCentrality][iTrackPt],-1.5,1.5,Form("DeltaPhiProjection%d",iFile));
          jffProjection[1][iJetTrack][iFile][iCentrality][iTrackPt]->Rebin(deltaPhiRebin);
          jffProjection[1][iJetTrack][iFile][iCentrality][iTrackPt]->Scale(1.0*projector->GetNBinsProjectedOver()/deltaPhiRebin);
          
          // Get the deltaR distribution from the DijetMethods
          jffProjection[2][iJetTrack][iFile][iCentrality][iTrackPt] = projector->GetJetShape(jffCorrectionDeltaEtaDeltaPhi[iJetTrack][iFile][iCentrality][iTrackPt]);
          jffProjection[2][iJetTrack][iFile][iCentrality][iTrackPt]->SetName(Form("jetShapeCorrection%d%d%d%d", iJetTrack, iFile, iCentrality, iTrackPt));
          
        }
        
      } // Track pT loop
    } // Centrality loop
  } // File loop
  
  // Read the corrections used in the inclusive jet-track paper
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      jffCorrectionInclusive[iCentrality][iTrackPt] = (TH1D*) inclusiveAnalysisJffPbPb->Get(Form("drResidualJff_Pb_%d_%d",iTrackPt,iCentrality));
    }
  }
  
  // Draw all the projections
  for(int iJetTrack = 0; iJetTrack < nJetTrack; iJetTrack++){
    if(!drawJetTrack[iJetTrack]) continue;
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      if(ppData && (iCentrality > 0)) continue; // No centrality selection for pp
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        for(int iDeltaAngle = 0; iDeltaAngle < 3; iDeltaAngle++){
          if(!drawAngle[iDeltaAngle]) continue;
          
          // First, draw the standard JFF projection for deltaEta/deltaPhi
          jffProjection[iDeltaAngle][iJetTrack][0][iCentrality][iTrackPt]->SetLineColor(fileColors[0]);
          jffProjection[iDeltaAngle][iJetTrack][0][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(minXrange[iDeltaAngle],1);
          drawer->DrawHistogram(jffProjection[iDeltaAngle][iJetTrack][0][iCentrality][iTrackPt],deltaLatex[iDeltaAngle],Form("1/N_{dijet} dN/d%s",deltaLatex[iDeltaAngle])," ");
          
          // Next, any other JFF projections for deltaEta/deltaPhi
          for(int iFile = 1; iFile < nFiles; iFile++){
            
            // Next, draw the rebinned projection for deltaEta/deltaPhi
            jffProjection[iDeltaAngle][iJetTrack][iFile][iCentrality][iTrackPt]->SetLineColor(fileColors[iFile]);
            jffProjection[iDeltaAngle][iJetTrack][iFile][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(minXrange[iDeltaAngle],1.5);
            jffProjection[iDeltaAngle][iJetTrack][iFile][iCentrality][iTrackPt]->Draw("same");
            
          }
          
          if(drawOldInclusive && iDeltaAngle == 2 && iJetTrack == 4){
            
            jffCorrectionInclusive[iCentrality][iTrackPt]->SetLineColor(fileColors[nFiles]);
            jffCorrectionInclusive[iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(minXrange[iDeltaAngle],1.5);
            jffCorrectionInclusive[iCentrality][iTrackPt]->Draw("same");
            
          }
          
          // Draw a legend to the canvas
          if(ppData){
            sprintf(centralityString,"Pythia6");
            sprintf(compactCentralityString,"");
          } else {
            sprintf(centralityString,"P+H %.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
            sprintf(compactCentralityString,"_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
          }
          
          sprintf(namer,"%.1f < p_{T} < %.1f GeV",trackPtBorders[iTrackPt],trackPtBorders[iTrackPt+1]);
          legend = new TLegend(0.58,0.17,0.88,0.57);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->SetHeader(headers[iJetTrack]);
          legend->AddEntry((TObject*)0,centralityString,"");
          legend->AddEntry((TObject*)0,namer,"");
          for(int iFile = nFiles-1; iFile >= 0; iFile--){
            legend->AddEntry(jffProjection[iDeltaAngle][iJetTrack][iFile][iCentrality][iTrackPt],fileLegend[iFile],"l");
          }
          if(drawOldInclusive && iDeltaAngle == 2 && iJetTrack == 4){
            legend->AddEntry(jffCorrectionInclusive[iCentrality][iTrackPt],"HIN-16-020","l");
          }
          legend->Draw();
          
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/jffProjection_%s%sComparison%s_pT=%.0f-%.0f.pdf", saveName[iJetTrack], deltaName[iDeltaAngle], compactCentralityString, trackPtBorders[iTrackPt], trackPtBorders[iTrackPt+1]));
          }
          
        } // DeltaPhi or DeltaEta loop
      } // Track pT loop
    } // Centrality loop
  } // Jet-track loop
  
}
