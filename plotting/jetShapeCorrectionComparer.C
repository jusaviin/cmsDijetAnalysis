#include "DijetCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

/*
 * Compare jet shape corrections
 */
void jetShapeCorrectionComparer(){
  
  // Open the files for the JFF corrections
  TFile *standardFile = TFile::Open("newPpTest2.root");
  TFile *wtaFile = TFile::Open("ppTestWTA.root");
  TFile *files[2] = {standardFile,wtaFile};
  
  // Figure saving options
  bool saveFigures = true;         // Flag to determine whather or not save the figures
  bool drawCorrection = true;      // Draw the comparison of JFF corrections (RecoGen - GenGen) between two files
  bool drawRatio = true;          // Draw the comparison of JFF ratios (RecoGen/GenGen) between two files
  
  // Configuration
  const int nCentralityBins = 4;
  const int nTrackPtBins = 6;
  double trackPtBorders[] = {0.7,1,2,3,4,8,300};
  double centralityBinBorders[] = {0,10,30,50,100};
  
  // Jff correction histograms
  // First bin = Leading, pT weighted leading, subeading, pT weighted subleading
  // Second bin = Standard, WTA
  TH1D *jffCorrection[4][2][nCentralityBins][nTrackPtBins];
  
  // All the variables needed for figure drawing
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  TLegend *legend;
  TH1D* ratio;
  const char *headers[4] = {"Track-leading jet","p_{T}w track-leading jet","Track-subleading jet","p_{T}w track-subleading jet"};
  const char *saveName[4] = {"trackLeadingJet","trackLeadingJetPtWeighted","trackSubleadingJet","trackSubleadingJetPtWeighted"};
  const char *comparisonType[2] = {"Correction","Ratio"};
  const char *axisName[2] = {"RecoGen - GenGen","RecoGen/GenGen"};
  double ratioZoomLow[2] = {-2,0.5};
  double ratioZoomHigh[2] = {4,1.5};
  
  // Check if we are using pp or PbPb data
  DijetCard *card = new DijetCard(standardFile);
  TString collisionSystem = card->GetDataType();
  bool ppData = (collisionSystem.Contains("pp") || collisionSystem.Contains("localTest"));
  
  // Read the JFF correction histograms from files
  char namer[200];
  char centralityString[50];
  char compactCentralityString[50];
  bool notThisType[2] = {!drawCorrection,!drawRatio};
  for(int iType = 0; iType < 2; iType++){
    if(notThisType[iType]) continue; // Only draw the selected types
    for(int iFile = 0; iFile < 2; iFile++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        if(ppData && (iCentrality > 0)) continue; // No centrality selection for pp
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          // Read the correction for leading jet
          sprintf(namer,"JetShape_trackLeadingJet/jff%s_JetShape_trackLeadingJet_C%dT%d",comparisonType[iType],iCentrality,iTrackPt);
          jffCorrection[0][iFile][iCentrality][iTrackPt] = (TH1D*) files[iFile]->Get(namer);
          
          // Read the correction for pT weighted leading jet
          sprintf(namer,"JetShape_trackLeadingJetPtWeighted/jff%s_JetShape_trackLeadingJetPtWeighted_C%dT%d",comparisonType[iType],iCentrality,iTrackPt);
          jffCorrection[1][iFile][iCentrality][iTrackPt] = (TH1D*) files[iFile]->Get(namer);
          
          // Read the correction for subleading jet
          sprintf(namer,"JetShape_trackSubleadingJet/jff%s_JetShape_trackSubleadingJet_C%dT%d",comparisonType[iType],iCentrality,iTrackPt);
          jffCorrection[2][iFile][iCentrality][iTrackPt] = (TH1D*) files[iFile]->Get(namer);
          
          // Read the correction for pT weighted subleading jet
          sprintf(namer,"JetShape_trackSubleadingJetPtWeighted/jff%s_JetShape_trackSubleadingJetPtWeighted_C%dT%d",comparisonType[iType],iCentrality,iTrackPt);
          jffCorrection[3][iFile][iCentrality][iTrackPt] = (TH1D*) files[iFile]->Get(namer);
          
        } // Track pT loop
      } // Centrality loop
    } // File loop
    
    // Draw all the corrections
    for(int iJetTrack = 0; iJetTrack < 4; iJetTrack++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        if(ppData && (iCentrality > 0)) continue; // No centrality selection for pp
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
          
          // First, draw the standard leading JFF correction
          drawer->CreateSplitCanvas();
          jffCorrection[iJetTrack][0][iCentrality][iTrackPt]->SetLineColor(kRed);
          jffCorrection[iJetTrack][0][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
          drawer->DrawHistogramToUpperPad(jffCorrection[iJetTrack][0][iCentrality][iTrackPt],"#DeltaR",axisName[iType]," ");
          
          // Next, draw the WTA leading JFF correction
          jffCorrection[iJetTrack][1][iCentrality][iTrackPt]->SetLineColor(kBlue);
          jffCorrection[iJetTrack][1][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
          jffCorrection[iJetTrack][1][iCentrality][iTrackPt]->Draw("same");
          
          // Draw a legend to the canvas
          if(ppData){
            sprintf(centralityString,"Pythia6");
            sprintf(compactCentralityString,"");
          } else {
            sprintf(centralityString,"P+H %.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
            sprintf(compactCentralityString,"_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
          }
          sprintf(namer,"%.1f < p_{T} < %.1f GeV",trackPtBorders[iTrackPt],trackPtBorders[iTrackPt+1]);
          legend = new TLegend(0.50,0.02,0.85,0.35);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->SetHeader(headers[iJetTrack]);
          legend->AddEntry((TObject*)0,centralityString,"");
          legend->AddEntry((TObject*)0,namer,"");
          legend->AddEntry(jffCorrection[iJetTrack][0][iCentrality][iTrackPt],"Standard","l");
          legend->AddEntry(jffCorrection[iJetTrack][1][iCentrality][iTrackPt],"WTA","l");
          legend->Draw();
          
          // Normalize the smeared jet pT to have the same yield as the reco jet pT
          sprintf(namer,"ratio%d%d",iJetTrack,iTrackPt);
          ratio = (TH1D*) jffCorrection[iJetTrack][1][iCentrality][iTrackPt]->Clone(namer);
          ratio->Divide(jffCorrection[iJetTrack][0][iCentrality][iTrackPt]);
          ratio->GetYaxis()->SetRangeUser(ratioZoomLow[iType],ratioZoomHigh[iType]);
          drawer->DrawHistogramToLowerPad(ratio,"#DeltaR","WTA/standard"," ");
          
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/jetShape%sComparison_%s%s_pT=%.0f-%.0f.pdf",comparisonType[iType],saveName[iJetTrack],compactCentralityString,trackPtBorders[iTrackPt],trackPtBorders[iTrackPt+1]));
          }
          
        } // Track pT loop
      } // Centrality loop
    } // Jet-track loop
  } // Type loop (correction or ratio)
}
