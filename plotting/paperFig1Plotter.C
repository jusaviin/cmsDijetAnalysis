
#include "stackHistXiao.h"
#include "xCanvas.h"

const int ncent = 5;
const int nasym = 4;
const int ntrack = 7;

TH1D**** loadHist_cent_trk_asym(TString name1, TFile *f){
  TH1D**** hist =new TH1D***[ncent];
  for(int i =0 ;i<ncent ; ++i){
    hist[i] = new TH1D**[ntrack];
    for(int j=0; j<ntrack; j++){
      hist[i][j]= new TH1D*[nasym];
      for(int l =0 ;l<nasym ; ++l){
        hist[i][j][l] = (TH1D*) f->Get(name1+Form("_%d_%d_%d", i,j,l));
      }
    }
  }
  return hist;
}

TH1D*** loadHist_cent_asym(TString name, TFile *f){
  TH1D*** hist =new TH1D**[ncent];
  for(int i =0 ;i<ncent ; ++i){
    hist[i] = new TH1D*[nasym];
    for(int l =0 ;l<nasym ; ++l){
      hist[i][l] = (TH1D*) f->Get(name+Form("_%d_%d", i,l));
    }
  }
  return hist;
}

std::pair<TH1D***, TH1D***>* loadHistPair_cent_asym(TString name1, TString name2, TFile *f){
  std::pair<TH1D***,TH1D***> * p = new std::pair<TH1D***,TH1D***>();
  p->first  = loadHist_cent_asym(name1, f);
  p->second = loadHist_cent_asym(name2, f);
  return p;
}

stackHist** makeStack(TString name, TH1D**** js, int iasym){
  stackHist **jstk= new stackHist*[ncent];
  for(int i =0 ;i<ncent ; ++i){
    jstk[i] = new stackHist(name+Form("_%d_%d",i,iasym));
    for(int l =0 ;l<ntrack ; ++l){
      jstk[i]->addHist((TH1*) js[i][l][iasym]);
    }
  }
  return jstk;
}

void drawStack(stackHist *jstk, bool isconner = 0){
  jstk->setRange(.5, 3000, "y");
  jstk->setRange(0.0, 0.99, "x");
  jstk->drawStack();
  jstk->hst->GetYaxis()->SetNdivisions(505);
  jstk->hst->GetYaxis()->SetLabelSize(0.09);
  jstk->hst->GetYaxis()->SetTitleOffset(1.2);
  jstk->hst->GetYaxis()->SetTitleSize(0.1);
  jstk->hst->GetYaxis()->SetTitleOffset(0.8);
  jstk->hst->GetYaxis()->SetTitleSize(0.14);
  jstk->hst->GetYaxis()->SetTitle("P(#Deltar)");
  
  jstk->hst->GetXaxis()->SetTitle("#Deltar");
  jstk->hst->GetXaxis()->SetTitleOffset(0.6);
  jstk->hst->GetXaxis()->SetTitleSize(0.14);
  jstk->hst->GetXaxis()->SetLabelSize(0.11);
  jstk->hst->GetXaxis()->SetLabelOffset(-0.01);
  if(isconner){
    jstk->hst->GetYaxis()->SetTitleSize(0.12);
    jstk->hst->GetYaxis()->SetTitleOffset(0.94);
    jstk->hst->GetYaxis()->SetLabelSize(0.08);
    jstk->hst->GetXaxis()->SetTitleOffset(0.79);
    jstk->hst->GetXaxis()->SetTitleSize(0.105);
    jstk->hst->GetXaxis()->SetLabelSize(0.085);
    jstk->hst->GetXaxis()->SetLabelOffset(0.014);
  }
  jstk->hst->Draw();
}

void drawTotal(std::pair<TH1D***,TH1D***>* p, int i, int j){
  auto h = (p->first)[i][j];
  auto e = (p->second)[i][j];
  e->SetFillColorAlpha(kGray+3, .4);
  e->SetMarkerColor(kBlack);
  e->SetMarkerSize(1.4);
  e->SetMarkerStyle(20);
  e->Draw("same e2");
}

void drawHist(TH1D* hist, TH1D* syst, int isconner = 0){
  syst->SetAxisRange(0, .99, "X");
  syst->SetAxisRange(0, 5, "Y");
  hist->SetMarkerSize(.8);
  hist->SetMarkerStyle(24);
  hist->SetMarkerColor(kBlack);
  hist->SetLineColor(kBlack);
  hist->SetTitle("");
  syst->SetFillColorAlpha(kGray+2, .4);
  syst->SetFillStyle(1001);
  syst->SetMarkerSize(1.);
  syst->SetMarkerStyle(4);
  syst->SetMarkerColor(kBlack);
  syst->SetTitle("");
  if(isconner){
  }else {
  }
  syst->Draw("same e2");
  hist->Draw("same");
}

void paperFig1Plotter(TString saveName, TH1D**** ld_js, std::pair<TH1D***,TH1D***>* ld_sum, TH1D**** sub_js, std::pair<TH1D***,TH1D***>*sub_sum, std::vector<TString>& txt){
  // txt[0] will be the title of the figure
  stackHist **jstk_ld  = makeStack("st_ld", ld_js , nasym-1);
  stackHist **jstk_sub = makeStack("st_ld", sub_js, nasym-1);
  TH1D *ratio[ncent];
  TH1D *ratio_err[ncent];
  for(int i=0; i<ncent; ++i){
    ratio[i] = (TH1D*)(sub_sum->first)[i][nasym-1]->Clone(Form("ratio_%d",i));
    ratio[i] -> Divide((ld_sum->first)[i][nasym-1]);
    ratio_err[i] =(TH1D*)(sub_sum->second)[i][nasym-1]->Clone(Form("ratio_err_%d",i));
    ratio_err[i] ->Divide((ld_sum->second)[i][nasym-1]);
  }
  
  TBox *box = new TBox();
  TLatex *mainTitle = new TLatex();
  TLine *line = new TLine();
  const char *pbpbLabel = "PbPb";
  const char *ppLabel = "pp";
  TString cent_lab[4] = {"0-10%", "10-30%", "30-50%", "50-90%"};
  auxi_canvas *bigCanvas = new auxi_canvas("fig1", "", 2500, 1550);
  bigCanvas->SetMargin(0.06, 0.01, 0.08, 0.3); // Margin order: Left, Right, Bottom, Top
  bigCanvas->divide(2,5);
  
  for(int i= 0; i< ncent; i++){
    
    //cout << "In the loop iAsymmetry: " << iAsymmetry << " iCentrality: " << iCentrality << endl;
    
    bigCanvas->CD(5-i);
    gPad->SetLogy();
    drawStack(jstk_ld[i]);
    drawTotal(ld_sum, i, 3);
    
    bigCanvas->CD(10-i);
    gPad->SetLogy();
    drawStack(jstk_sub[i], i==ncent-1);
    drawTotal(sub_sum, i, 3);
    
    //				bigCanvas->CD(15-i);
    //				drawHist(ratio[i], ratio_err[i]);
  }
  for(int i=0; i<ncent; ++i){
    if(i == ncent-1 ){
      bigCanvas->CD(5-i);
      mainTitle->SetTextFont(22);
      mainTitle->SetTextSize(0.12);
      mainTitle->DrawLatexNDC(0.35, 0.88, ppLabel);
      mainTitle->SetTextSize(0.10);
      mainTitle->DrawLatexNDC(0.55, 0.88, "Leading jet");
      bigCanvas->CD(10-i);
      mainTitle->SetTextSize(0.09);
      mainTitle->DrawLatexNDC(0.49, 0.88, "Subleading jet");
    }
    else {
      bigCanvas->CD(5-i);
      mainTitle->SetTextFont(22);
      mainTitle->SetTextSize(0.11);
      mainTitle->DrawLatexNDC(0.19, 0.9, pbpbLabel);
      mainTitle->DrawLatexNDC(0.19, 0.80, cent_lab[i]);
    }
  }
  bigCanvas->cd(0);
  mainTitle->SetTextFont(62);
  mainTitle->SetTextSize(0.06);
  mainTitle->DrawLatexNDC(.12,  .958, "CMS");
  mainTitle->SetTextFont(42);
  mainTitle->SetTextSize(0.055);
  mainTitle->DrawLatexNDC(.21, 0.958, "Particle transverse momentum distributions in jets");
  mainTitle->SetTextSize(0.05);
  mainTitle->DrawLatexNDC(0.26, 0.91, "5.02 TeV   pp 320 pb^{-1}   PbPb 1.7 nb^{-1}");
  mainTitle->SetTextSize(0.042);
  mainTitle->DrawLatexNDC(0.15, 0.86, "anti-k_{T} R = 0.4, |#eta_{jet}| < 1.6, p_{T,1} > 120 GeV, p_{T,2} > 50 GeV, #Delta#phi_{1,2} > #frac{5#pi}{6}");
  TLegend* ptLegend1 = new TLegend(0.04, 0.72, 0.27, 0.82);
  TLegend* ptLegend2 = new TLegend(0.28, 0.72, 0.51, 0.82);
  TLegend* ptLegend3 = new TLegend(0.53 ,0.72, 0.76, 0.82);
  TLegend* ptLegend4 = new TLegend(0.77, 0.72, 1.00, 0.82);
  ptLegend1->SetLineColor(kWhite);
  ptLegend1->SetFillColor(kWhite);
  ptLegend2->SetLineColor(kWhite);
  ptLegend2->SetFillColor(kWhite);
  ptLegend3->SetLineColor(kWhite);
  ptLegend3->SetFillColor(kWhite);
  ptLegend4->SetLineColor(kWhite);
  ptLegend4->SetFillColor(kWhite);
  
  mainTitle->SetTextFont(42);
  mainTitle->SetTextSize(0.035);
  TLegend* lt4 = new TLegend(0.25, 0.85, 0.85, 0.95);
  lt4->SetTextSize(0.06);
  lt4->SetLineColor(kWhite);
  lt4->SetFillColor(kWhite);
  
  bigCanvas->cd(0);
  mainTitle->SetTextFont(42);
  mainTitle->SetTextSize(0.038);
  box->SetFillColor(kWhite);
  float boxwidth = 0.02, boxhight = 0.034;
  float x= 0.236, y=0.0427, offsetx = 0.005, offsety = 0.005;
  box->DrawBox(x, y, x+boxwidth, y+boxhight); mainTitle->DrawLatex(x+offsetx, y+offsety, "0");
  //      box->SetFillColor(kAzure+6);
  x= 0.423, box->DrawBox(x, y, x+boxwidth, y+boxhight); mainTitle->DrawLatex(x+offsetx, y+offsety, "0");
  x= 0.608, box->DrawBox(x, y, x+boxwidth, y+boxhight); mainTitle->DrawLatex(x+offsetx, y+offsety, "0");
  x= 0.795, box->DrawBox(x, y, x+boxwidth, y+boxhight); mainTitle->DrawLatex(x+offsetx, y+offsety, "0");
  x= 0.978,box->DrawBox(x, y, x+boxwidth, y+boxhight); mainTitle->DrawLatex(x+offsetx, y+offsety, "1");
  
  ptLegend1->SetTextSize(0.03);
  ptLegend2->SetTextSize(0.03);
  ptLegend3->SetTextSize(0.03);
  ptLegend4->SetTextSize(0.03);
  ptLegend1->AddEntry(jstk_ld[0]->hist_trunk.at(0), "0.7 < p_{T}^{trk}< 1 GeV","f");
  ptLegend1->AddEntry(jstk_ld[0]->hist_trunk.at(1), "1 < p_{T}^{trk}< 2 GeV","f");
  
  ptLegend2->AddEntry(jstk_ld[0]->hist_trunk.at(2), "2 < p_{T}^{trk}< 3 GeV","f");
  ptLegend2->AddEntry(jstk_ld[0]->hist_trunk.at(3), "3 < p_{T}^{trk}< 4 GeV","f");
  
  ptLegend3->AddEntry(jstk_ld[0]->hist_trunk.at(4), "4 < p_{T}^{trk}< 8 GeV","f");
  ptLegend3->AddEntry(jstk_ld[0]->hist_trunk.at(5), "8 < p_{T}^{trk}< 12 GeV","f");
  
  ptLegend4->AddEntry(jstk_ld[0]->hist_trunk.at(6), "12 < p_{T}^{trk}< 300 GeV","f");
  ptLegend4->AddEntry((ld_sum->second)[4][3], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
  
  ptLegend1->Draw();
  ptLegend2->Draw();
  ptLegend3->Draw();
  ptLegend4->Draw();
  
  bigCanvas->SaveAs(saveName);
  //        mainTitle->DrawLatexNDC(jetShapeTitlePosition[iJetTrack/3]-0.08, 0.97, jetShapeTitle[iJetTrack/3]);
}
