
#include "paperFig1Plotter.C"

void drawRatio(TH1D** ratio, TH1D** err, int i, Color_t color,  bool isconner = 0){
  auto hist = ratio[i];
  auto syst = err[i];
  syst->SetAxisRange(0, .99, "X");
  syst->SetAxisRange(0.05, 3.9, "Y");
  hist->SetMarkerSize(1.5);
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(kWhite);
  hist->SetLineColor(color);
  hist->SetTitle("");
  hist->SetStats(0);
  syst->SetLineColor(color);
  syst->SetFillColorAlpha(color, .3);
  syst->SetFillStyle(1001);
  syst->SetMarkerSize(1.8);
  syst->SetMarkerStyle(4);
  syst->SetMarkerColor(color);
  syst->GetYaxis()->SetNdivisions(505);
  syst->GetXaxis()->SetNdivisions(505);
  syst->GetXaxis()->SetTitle("#Deltar");
  syst->GetXaxis()->CenterTitle();
  syst->GetYaxis()->CenterTitle();
  syst->GetYaxis()->SetTitle("P(#Deltar)_{PbPb}/P(#Deltar)_{pp}");
  syst->SetStats(0);
  syst->SetTitle("");
  if(isconner){
    syst->GetXaxis()->SetTitleSize(0.11);
    syst->GetXaxis()->SetTitleOffset(0.78);
    syst->GetYaxis()->SetTitleSize(0.08);
    syst->GetYaxis()->SetTitleOffset(1);
    syst->GetYaxis()->SetLabelSize(0.08);
    syst->GetXaxis()->SetLabelSize(0.08);
    syst->GetYaxis()->SetLabelOffset(0.012);
    syst->GetXaxis()->SetLabelOffset(0.006);
    cout<<"connor!"<<endl;
  }else {
    syst->GetXaxis()->SetTitleSize(0.13);
    syst->GetXaxis()->SetTitleOffset(0.65);
    syst->GetYaxis()->SetTitleSize(0.08);
    syst->GetYaxis()->SetTitleSize(0.1);
    syst->GetYaxis()->SetTitleOffset(0.8);
    syst->GetYaxis()->SetLabelSize(0.1);
    syst->GetXaxis()->SetLabelSize(0.1);
    syst->GetYaxis()->SetLabelOffset(0.01);
    syst->GetXaxis()->SetLabelOffset(-0.009);
  }
  syst->Draw("same e2");
  hist->Draw("same");
}

void paperFig2Plotter(TString saveName, std::pair<TH1D***,TH1D***>* ld_sum, std::pair<TH1D***,TH1D***>*sub_sum){
  // txt[0] will be the title of the figure
  TH1D *ld_ratio[ncent];
  TH1D *ld_ratio_err[ncent];
  TH1D *sub_ratio[ncent];
  TH1D *sub_ratio_err[ncent];
  for(int i=0; i<ncent-1; ++i){
    ld_ratio[i] =  (TH1D*)(ld_sum->first)[i][nasym-1]->Clone(Form("ld_ratio_%d",i));
    ld_ratio[i] -> Divide((ld_sum->first)[ncent-1][nasym-1]);
    ld_ratio_err[i] =(TH1D*) (ld_sum->second)[i][nasym-1]->Clone(Form("ld_ratio_err_%d",i));
    ld_ratio_err[i] ->Divide((ld_sum->second)[ncent-1][nasym-1]);
    sub_ratio[i]  = (TH1D*)(sub_sum->first)[i][nasym-1]->Clone(Form("sub_ratio_%d",i));
    sub_ratio[i] -> Divide((sub_sum->first)[ncent-1][nasym-1]);
    sub_ratio_err[i]  =(TH1D*)(sub_sum->second)[i][nasym-1]->Clone(Form("sub_ratio_err_%d",i));
    sub_ratio_err[i] ->Divide((sub_sum->second)[ncent-1][nasym-1]);
  }
  
  TBox *box = new TBox();
  TLatex *mainTitle = new TLatex();
  TLine *line = new TLine();
  line->SetLineStyle(2);
  TString pbpbLabel = "PbPb";
  TString ppLabel = "pp";
  TString cent_lab[4] = {"0-10%", "10-30%", "30-50%", "50-90%"};
  auxi_canvas *bigCanvas = new auxi_canvas("fig2", "", 2000, 1100);
  bigCanvas->SetMargin(0.06, 0.01, 0.1, 0.1); // Margin order: Left, Right, Bottom, Top
  bigCanvas->divide(2,4);
  
  for(int i= 0; i< ncent-1; i++){
    
    //cout << "In the loop iAsymmetry: " << iAsymmetry << " iCentrality: " << iCentrality << endl;
    bigCanvas->CD(4-i);
    drawRatio(ld_ratio, ld_ratio_err, i, kAzure+2);
    line->DrawLine(0,1 ,1 ,1);
    
    bigCanvas->CD(8-i);
    drawRatio(sub_ratio, sub_ratio_err, i, kRed+1, i == 3);
    line->DrawLine(0,1 ,1 ,1);
    
  }
  mainTitle->SetTextFont(22);
  for(int i=0; i<ncent-1; ++i){
    if(i == ncent-2 ){
      bigCanvas->CD(4-i);
      mainTitle->SetTextSize(0.10);
      //						mainTitle->DrawLatexNDC(0.5, 0.88, "Leading jet");
      mainTitle->DrawLatexNDC(0.3, 0.9, "Cent: "+cent_lab[i]);
      bigCanvas->CD(8-i);
      mainTitle->SetTextSize(0.08);
      //						mainTitle->DrawLatexNDC(0.49, 0.88, "Subleading jet");
      
    }
    else {
      bigCanvas->CD(4-i);
      mainTitle->SetTextSize(0.10);
      mainTitle->DrawLatexNDC(0.1, 0.9, "Cent: "+cent_lab[i]);
    }
  }
  TLegend* ptLegend1 = new TLegend(0.3, 0.45, 0.77, 0.85);
  ptLegend1->SetLineColor(kWhite);
  ptLegend1->AddEntry(ld_ratio_err[0], "Leading jet", "lpf");
  ptLegend1->AddEntry(sub_ratio_err[0], "Subleading jet", "lpf");
  bigCanvas->CD(1);
  ptLegend1->Draw();
  
  bigCanvas->cd(0);
  mainTitle->SetTextFont(42);
  mainTitle->SetTextSize(0.045);
  box->SetFillColor(kWhite);
  float boxwidth = 0.02, boxhight = 0.034;
  float offsetx = 0.005, offsety = 0.005; 
  float x= 0.28, y=0.06; box->DrawBox(x, y, x+boxwidth, y+boxhight); mainTitle->DrawLatex(x+offsetx, y+offsety, "0");
  x= 0.515, box->DrawBox(x, y, x+boxwidth, y+boxhight); mainTitle->DrawLatex(x+offsetx, y+offsety, "0");
  x= 0.747, box->DrawBox(x, y, x+boxwidth, y+boxhight); mainTitle->DrawLatex(x+offsetx, y+offsety, "0");
  x= 0.978,box->DrawBox(x, y, x+boxwidth, y+boxhight); mainTitle->DrawLatex(x+offsetx, y+offsety, "1");
  
  
  mainTitle->SetTextFont(62);
  mainTitle->SetTextSize(0.06);
  mainTitle->DrawLatexNDC(.01,  .93, "CMS");
  mainTitle->SetTextFont(42);
  mainTitle->SetTextSize(0.045);
  mainTitle->DrawLatexNDC(.086, 0.93, "Particle transverse momentum distribution ratios");
  mainTitle->SetTextSize(0.035);
  mainTitle->DrawLatexNDC(0.64, 0.96, "5.02 TeV   pp 320 pb^{-1}   PbPb 1.7 nb^{-1}");
  mainTitle->SetTextSize(0.027);
  mainTitle->DrawLatexNDC(0.58, 0.92, "anti-k_{T} R = 0.4, |#eta_{jet}| < 1.6, p_{T,1} > 120 GeV, p_{T,2} > 50 GeV, #Delta#phi_{1,2} > #frac{5#pi}{6}");
  
  
  
  bigCanvas->SaveAs(saveName);
  //        mainTitle->DrawLatexNDC(jetShapeTitlePosition[iJetTrack/3]-0.08, 0.97, jetShapeTitle[iJetTrack/3]);
}
