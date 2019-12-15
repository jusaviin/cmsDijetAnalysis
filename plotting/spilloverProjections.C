#include "DijetHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "DijetCard.h"
#include "DijetMethods.h"
#include "JffCorrector.h"
#include "JDrawer.h"

/*
 * Macro for comparing spillover corrections from different sources
 */
void spilloverProjections(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString spilloverFileName = "corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncorr_improvisedMixing_wtaAxis_centShift5_jet100trigger_JECv6_2019-10-04.root";
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncorr_xjBins_improvisedMixing_wtaAxis_jet100trigger_JECv6_2019-09-26.root
  // corrections/spilloverCorrection_PbPbMC_akFlowPuCsPfJets_xjBins_noUncorr_improviseMixing_wta_cutFluctuation_preliminary_2019-08-16.root
  TString spilloverComparisonFileName = "corrections/spilloverCorrection_PbPbMC_akPu4CaloJets_xjBins_noUncorr_improvisedMixing_wtaAxis_JECv5b_preliminary_2019-09-08.root";
  // corrections/spilloverCorrection_PbPbMC_akPu4CaloJets_xjBins_noUncorr_improvisedMixing_wtaAxis_JECv5b_preliminary_2019-09-08.root
  // corrections/spilloverCorrection_akFlowPuCs4PFJet_noUncOrInc_improvisedMixing_symmetrized_looseCut_tightForSubleading_wtaAxis_JECv6_2019-10-23.root
  TString inclusiveUncertaintyFileName = "uncertainties/inclusiveAnalysis/js_AllSource_syst_err.root";

  bool saveFigures = false;
  
  // Open the input files
  TFile *spilloverFile[2];
  spilloverFile[0] = TFile::Open(spilloverFileName);
  spilloverFile[1] = TFile::Open(spilloverComparisonFileName);
  
  TFile *inclusiveUncertaintyFile = TFile::Open(inclusiveUncertaintyFileName);
  
  const int nCentralityBins = 4;
  const int nTrackPtBins = 5;
  const int nAsymmetryBins = 3;
  double centralityBinBorders[] = {0,10,30,50,90};  // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double xjBinBorders[] = {0,0.6,0.8,1}; // Bin borders for xj
  int iJetTrack = DijetHistogramManager::kPtWeightedTrackLeadingJet; // kTrackLeadingJet kPtWeightedTrackLeadingJet
  
  // Define arrays for the histograms
  TH2D *spilloverHistogram[2][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  
  // Projections from two dimensional histograms
  TH1D *spilloverEtaProjection[2][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *spilloverPhiProjection[2][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];
  TH1D *spilloverRProjection[2][nAsymmetryBins+1][nCentralityBins][nTrackPtBins];

  // Inclusive uncertainty for spillover (50 % of the correction)
  TH1D *spilloverUncertaintyInclusive[nCentralityBins][nTrackPtBins];
  TH1D *spilloverInclusiveRatio[nCentralityBins][nTrackPtBins];
  
  // Define a DijetMethods to get the spillover correction as a function of DeltaR
  DijetMethods *calculator = new DijetMethods();
  
  // Define JFF corrector to easily read the spillover histograms from the files
  JffCorrector *spilloverReader[2];
  for(int iFile = 0; iFile < 2; iFile++){
    spilloverReader[iFile] = new JffCorrector();
    spilloverReader[iFile]->ReadSpilloverFile(spilloverFile[iFile]);
  }
  
  // String needed to read the inclusive spillover corrections from the file
  TString inclusiveSources[] = {"rel_spill_err","jff_err","rel_bkg_err","rel_rela_err"};
  TString inclusiveCentrality[] = {"Cent0_Cent10","Cent10_Cent30","Cent30_Cent50","Cent50_Cent100"};
  TString inclusiveTrackPt[] = {"TrkPt07_TrkPt1","TrkPt1_TrkPt2","TrkPt2_TrkPt3","TrkPt3_TrkPt4","TrkPt4_TrkPt8","TrkPt8_TrkPt12","TrkPt12_TrkPt16","TrkPt16_TrkPt20","TrkPt20_TrkPt300"};
  
  // Read the histograms from spillover file and do projections
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      for(int iFile = 0; iFile < 2; iFile++){
        spilloverHistogram[iFile][nAsymmetryBins][iCentrality][iTrackPt] = spilloverReader[iFile]->GetDeltaEtaDeltaPhiSpilloverCorrection(iJetTrack, iCentrality, iTrackPt);
        
        spilloverEtaProjection[iFile][nAsymmetryBins][iCentrality][iTrackPt] = calculator->ProjectRegionDeltaEta(spilloverHistogram[iFile][nAsymmetryBins][iCentrality][iTrackPt], -1, 1, Form("etaProjection%d%d%d", iFile, iCentrality, iTrackPt));
        spilloverEtaProjection[iFile][nAsymmetryBins][iCentrality][iTrackPt]->Scale(calculator->GetNBinsProjectedOver());
        
        spilloverPhiProjection[iFile][nAsymmetryBins][iCentrality][iTrackPt] = calculator->ProjectRegionDeltaPhi(spilloverHistogram[iFile][nAsymmetryBins][iCentrality][iTrackPt], 0, 1, Form("phiProjection%d%d%d", iFile, iCentrality, iTrackPt));
        spilloverPhiProjection[iFile][nAsymmetryBins][iCentrality][iTrackPt]->Scale(calculator->GetNBinsProjectedOver());
        
        spilloverRProjection[iFile][nAsymmetryBins][iCentrality][iTrackPt] = calculator->GetJetShape(spilloverHistogram[iFile][nAsymmetryBins][iCentrality][iTrackPt]);
        spilloverRProjection[iFile][nAsymmetryBins][iCentrality][iTrackPt]->SetName(Form("jetShape%d%d%d", iFile, iCentrality, iTrackPt));
        
      } // File loop
      
      spilloverUncertaintyInclusive[iCentrality][iTrackPt] = (TH1D*) inclusiveUncertaintyFile->Get(Form("dR_Syst_PbPb_%s_%s_rel_spill_err", inclusiveCentrality[iCentrality].Data(), inclusiveTrackPt[iTrackPt].Data()));
      
    } // Track pT loop
  } // Centrality loop
  
  
  // Overlay the projections from different files
  JDrawer *drawer = new JDrawer();
  TLegend *legend;
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  TH1D *helperHistogram;
  
  cout << "Integral sanity check: " << endl;
  int bin1 = spilloverEtaProjection[0][nAsymmetryBins][0][1]->GetXaxis()->FindBin(-0.99);
  int bin2 = spilloverEtaProjection[0][nAsymmetryBins][0][1]->GetXaxis()->FindBin(0.99);
  cout << "Integral of deltaEta: " << spilloverEtaProjection[0][nAsymmetryBins][0][1]->Integral(bin1,bin2,"width") << endl;
  bin1 = spilloverPhiProjection[0][nAsymmetryBins][0][1]->GetXaxis()->FindBin(-0.99);
  bin2 = spilloverPhiProjection[0][nAsymmetryBins][0][1]->GetXaxis()->FindBin(0.99);
  cout << "Integral of deltaPhi: " << spilloverPhiProjection[0][nAsymmetryBins][0][1]->Integral(bin1,bin2,"width") << endl;
  bin1 = spilloverRProjection[0][nAsymmetryBins][0][1]->GetXaxis()->FindBin(0.01);
  bin2 = spilloverRProjection[0][nAsymmetryBins][0][1]->GetXaxis()->FindBin(0.99);
  cout << "Integral of deltaR: " << spilloverRProjection[0][nAsymmetryBins][0][1]->Integral(bin1,bin2,"width") << endl;
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    
    if(iCentrality == nCentralityBins){
      centralityString = "pp";
      compactCentralityString = "_pp";
    } else {
      centralityString = Form("C = %.0f-%.0f %%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
      compactCentralityString = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
    }
    
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      
      trackPtString = Form("Track pT: %.1f-%.1f GeV", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
      compactTrackPtString = Form("_pT=%.1f-%.1f", trackPtBinBorders[iTrackPt], trackPtBinBorders[iTrackPt+1]);
      compactTrackPtString.ReplaceAll(".","v");
      
      // ===================================
      // == Draw the deltaEta projections ==
      // ===================================
      
      spilloverEtaProjection[0][nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(kBlack);
      spilloverEtaProjection[1][nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(kRed);
      
      spilloverEtaProjection[0][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-1,1);
      drawer->DrawHistogram(spilloverEtaProjection[0][nAsymmetryBins][iCentrality][iTrackPt], "#Delta#eta", "#frac{dN}{d#Delta#eta}", " ");
      spilloverEtaProjection[1][nAsymmetryBins][iCentrality][iTrackPt]->Draw("same");
      
      legend = new TLegend(0.7,0.65,0.9,0.92);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0,centralityString.Data(),"");
      legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
      legend->AddEntry(spilloverEtaProjection[0][nAsymmetryBins][iCentrality][iTrackPt], "Flow PF jet", "l");
      legend->AddEntry(spilloverEtaProjection[1][nAsymmetryBins][iCentrality][iTrackPt], "Calo jet", "l");
      legend->Draw();
      
      // ===================================
      // == Draw the deltaPhi projections ==
      // ===================================
      
      spilloverPhiProjection[0][nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(kBlack);
      spilloverPhiProjection[1][nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(kRed);
      
      spilloverPhiProjection[0][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(-1,1);
      drawer->DrawHistogram(spilloverPhiProjection[0][nAsymmetryBins][iCentrality][iTrackPt], "#Delta#varphi", "#frac{dN}{d#Delta#varphi}", " ");
      spilloverPhiProjection[1][nAsymmetryBins][iCentrality][iTrackPt]->Draw("same");
      
      legend = new TLegend(0.7,0.65,0.9,0.92);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0,centralityString.Data(),"");
      legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
      legend->AddEntry(spilloverPhiProjection[0][nAsymmetryBins][iCentrality][iTrackPt], "Flow PF jet", "l");
      legend->AddEntry(spilloverPhiProjection[1][nAsymmetryBins][iCentrality][iTrackPt], "Calo jet", "l");
      legend->Draw();
      
      // =================================
      // == Draw the deltaR projections ==
      // =================================
      
      spilloverRProjection[0][nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(kBlack);
      spilloverRProjection[1][nAsymmetryBins][iCentrality][iTrackPt]->SetLineColor(kRed);
      spilloverUncertaintyInclusive[iCentrality][iTrackPt]->SetLineColor(kBlue);
      
      // Scale the histograms by a fraction of the correction
      //spilloverRProjection[0][nAsymmetryBins][iCentrality][iTrackPt]->Scale(0.15);
      //spilloverRProjection[1][nAsymmetryBins][iCentrality][iTrackPt]->Scale(0.5);
      
      spilloverRProjection[0][nAsymmetryBins][iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
      drawer->DrawHistogram(spilloverRProjection[0][nAsymmetryBins][iCentrality][iTrackPt], "#Deltar", "#frac{dN}{d#Deltar}", " ");
      spilloverRProjection[1][nAsymmetryBins][iCentrality][iTrackPt]->Draw("same");
      spilloverUncertaintyInclusive[iCentrality][iTrackPt]->Draw("same");
      
      legend = new TLegend(0.7,0.65,0.9,0.92);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0,centralityString.Data(),"");
      legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
      legend->AddEntry(spilloverRProjection[0][nAsymmetryBins][iCentrality][iTrackPt], "Flow PF jet", "l");
      legend->AddEntry(spilloverRProjection[1][nAsymmetryBins][iCentrality][iTrackPt], "Calo jet", "l");
      legend->AddEntry(spilloverUncertaintyInclusive[iCentrality][iTrackPt], "Inclusive", "l");
      legend->Draw();
      
      // =============================================
      // == Take a ratio between inclusive and calo ==
      // =============================================
      
      spilloverInclusiveRatio[iCentrality][iTrackPt] = (TH1D*) spilloverRProjection[0][nAsymmetryBins][iCentrality][iTrackPt]->Clone(Form("spilloverRatio%d%d", iCentrality, iTrackPt));
      //helperHistogram = (TH1D*) spilloverUncertaintyInclusive[iCentrality][iTrackPt]->Clone(Form("spilloverRatio%d%d", iCentrality, iTrackPt));
      //for(int iBin = 1; iBin <= helperHistogram->GetNbinsX(); iBin++){
      //  helperHistogram->SetBinContent(iBin, spilloverRProjection[1][nAsymmetryBins][iCentrality][iTrackPt]->GetBinContent(iBin));
      //  helperHistogram->SetBinError(iBin, spilloverRProjection[1][nAsymmetryBins][iCentrality][iTrackPt]->GetBinError(iBin));
      //}
      //spilloverInclusiveRatio[iCentrality][iTrackPt]->Divide(helperHistogram);
      spilloverInclusiveRatio[iCentrality][iTrackPt]->Divide(spilloverRProjection[1][nAsymmetryBins][iCentrality][iTrackPt]);
      spilloverInclusiveRatio[iCentrality][iTrackPt]->GetXaxis()->SetRangeUser(0,1);
      drawer->DrawHistogram(spilloverInclusiveRatio[iCentrality][iTrackPt], "#Deltar", "Ratio", " ");
      
      legend = new TLegend(0.7,0.65,0.9,0.92);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0,centralityString.Data(),"");
      legend->AddEntry((TObject*) 0,trackPtString.Data(),"");
      legend->AddEntry(spilloverInclusiveRatio[iCentrality][iTrackPt], "PF/calo", "l");
      legend->Draw();
      
      // Save the figures into a file
      //if(saveFigures){
      //  gPad->GetCanvas()->SaveAs(Form("figures/spilloverComparisonCentralityShift_C=%.0f-%.0f.pdf", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]));
      //}
      
    } // Track pT loop
  } // Centrality loop
    
}
