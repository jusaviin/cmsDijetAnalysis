/*
 * Class for convanien drawing af stacked histograms
 * Code written by Xiao Wang
 */

#ifndef STACKHIST_H
#define STACKHIST_H

#include "THStack.h"
#include "TCanvas.h"
#include "vector"
#include "iostream"

class stackHist{
	public:
		TString xtitle;
		TString ytitle;
		std::vector<TH1*> hist_trunk;
		std::vector<TH1*> hist_trunk_up;
		std::vector<TH1*> hist_trunk_down;
		std::vector<Color_t>* color_v;
		TH1* hist_sum;
		THStack* hst;
		THStack* hst_up;
		THStack* hst_down;
		TString st_name;
		int lineW;
		float xmin, xmax, ymin, ymax;
		float xmind, xmaxd, ymind, ymaxd;
		float xtitleSize, ytitleSize;
		float xlabelSize, ylabelSize;
		//TLegend* tl_stack;
		bool doDiff;
		stackHist(TString name, TString title ="")
		{
			st_name = name;
			hst = new THStack(name, title);	
			hst_up = new THStack(name+"_up", title);	
			hst_down = new THStack(name+"_down", title);	
			hist_sum = NULL;
			color_v = NULL;
			xmin =0; ymin =0; xmax =1; ymax =1;
			xtitle = "";
			ytitle = "";
			lineW = 1;
			doDiff=false;
			xtitleSize=0.07;
			ytitleSize=0.07;
			xlabelSize=0.07;
			ylabelSize=0.07;
		}

		~stackHist();

		void addHist(TH1* h);
		void addHist(TH1** h ,int num, float *weight = NULL);
		void addDiff(TH1** h ,int num, float *weight = NULL);
		TH1* sumOver();
		TH1* drawStack(TString opt ="", TString addOpt = "hist");
		TH1* drawDiff(TString opt ="", TString addOpt = "hist");
		void stackConfig(THStack* hstt);
		void setTitleSize(float, TString axis="x");
		void setLabelSize(float, TString axis="x");
		void setRange(float, float, TString ax = "x");
		void setFillColor();
		void drawSum(TString opt = "same");
		TLegend* makeLegend(TString* tl_txt, float x1=0.1,float y1=0.05,float x2=0.9,float y2=0.9,\
				bool isDiff = false, size_t n=0);
};

void stackHist::addHist(TH1* h){
	int nn = hist_trunk.size();
	hist_trunk.push_back((TH1*) h->Clone(Form("_%d", nn)));
}
void stackHist::addHist(TH1** h ,int num, float *weight){
	TH1** hh = h;
	for(int i=0;i<num;i++){
		TString stemp = st_name+Form("_%d",i);
		hist_trunk.push_back((TH1*)hh[i]->Clone(stemp));
		if(weight != NULL) hist_trunk.back()->Scale(weight[i]);
	}
}

void stackHist::addDiff(TH1** h ,int num, float *weight){
	doDiff= true;
	TH1** hh = h;
	TH1* tem_up;
	TH1* tem_down;
	TString stemp;
	for(int i=0;i<num;i++){
		stemp = st_name+Form("_up_%d",i);
		tem_up= (TH1*)hh[i]->Clone(stemp);
		stemp = st_name+Form("_down_%d",i);
		tem_down= (TH1*)hh[i]->Clone(stemp);
		if(weight!=NULL) {
			tem_down->Scale(weight[i]); 
			tem_up  ->Scale(weight[i]);
		}
		for(int j=1;j<hh[i]->GetNbinsX()+1;j++){
			if( hh[i]->GetBinContent(j) <0){
				tem_up->SetBinContent(j,0);
				tem_up->SetBinError(j,0);
			}
			else {
				tem_down->SetBinContent(j,0);
				tem_down->SetBinError(j,0);
			}
		}
		stemp = st_name+Form("_up_%d",i);
//		hist_trunk_up.push_back((TH1*)tem_up->Clone(stemp));
		hist_trunk_up.push_back((TH1*)tem_up);
		stemp = st_name+Form("_down_%d",i);
//		hist_trunk_down.push_back((TH1*)tem_down->Clone(stemp));
		hist_trunk_down.push_back((TH1*)tem_down);
		if(weight != NULL) {
			hist_trunk_up.back()->Scale(weight[i]);
			hist_trunk_down.back()->Scale(weight[i]);
		}
//		tem_up->Delete();
//		tem_down->Delete();
	}
}

TH1* stackHist::sumOver(){
	TString stemp = st_name+"_sum";
	hist_sum = (TH1*)hist_trunk.front()->Clone(stemp);
	for(int i=1;i<hist_trunk.size();i++){
		hist_sum->Add(hist_trunk.at(i));
	}
	hist_sum->SetLineWidth(1);
	return hist_sum;
}

TH1* stackHist::drawStack(TString opt, TString addOpt ){
	setFillColor();
    if(opt == ""){
      for(auto it = hist_trunk.begin(); it !=hist_trunk.end(); ++it){
        hst->Add((TH1*)*it, addOpt);
        (*it)->SetLineWidth(lineW);
      }
    }
    else if(opt == "r"){
      for(auto it = hist_trunk.rbegin(); it !=hist_trunk.rend(); ++it){
        hst->Add((TH1*)*it, addOpt);
        (*it)->SetLineWidth(lineW);
      }
    }
  	hst->Draw();
    stackConfig(hst);
    hst->GetXaxis()->SetRangeUser(xmin, xmax);
    hst->SetMinimum(ymin);
    hst->SetMaximum(ymax);
    hst->Draw();
	return hst->GetHistogram();
}

TH1* stackHist::drawDiff(TString opt, TString addOpt){
	setFillColor();
	if(opt == ""){
		for(auto it = hist_trunk_up.begin(); it !=hist_trunk_up.end(); ++it){
			hst_up->Add((TH1*)*it, addOpt);
			(*it)->SetLineWidth(lineW);
		}
		for(auto it = hist_trunk_down.begin(); it !=hist_trunk_down.end(); ++it){
			hst_down->Add((TH1*)*it, addOpt);
			(*it)->SetLineWidth(lineW);
		}
	}
	else if(opt == "r"){
		for(auto it = hist_trunk_up.rbegin(); it !=hist_trunk_up.rend(); ++it){
			hst_up->Add((TH1*)*it, addOpt);
			(*it)->SetLineWidth(lineW);
		}
		for(auto it = hist_trunk_down.rbegin(); it !=hist_trunk_down.rend(); ++it){
			hst_down->Add((TH1*)*it, addOpt);
			(*it)->SetLineWidth(lineW);
		}
	}
	hst_up->Draw();
	stackConfig(hst_up);
	hst_down->Draw("same");
	hst_up->GetXaxis()->SetRangeUser(xmind, xmaxd);
	hst_up->SetMinimum(ymind);
	hst_up->SetMaximum(ymaxd);
	hst_up->Draw();
	hst_down->Draw("same");
	return hst_up->GetHistogram();
}

void stackHist::stackConfig(THStack* hstt){
	hstt->GetYaxis()->SetLabelSize(ylabelSize);
	hstt->GetXaxis()->SetLabelSize(xlabelSize);
	hstt->GetYaxis()->SetTitleSize(ytitleSize);
	hstt->GetXaxis()->SetTitleSize(xtitleSize);
	hstt->GetYaxis()->SetTitleOffset(0.8);
	hstt->GetXaxis()->SetTitleOffset(0.7);
	hstt->GetXaxis()->SetNdivisions(505);
	hstt->GetXaxis()->SetTickLength(0.03);
	hstt->GetYaxis()->SetTickLength(0.028);
	hstt->GetYaxis()->CenterTitle(true);
	hstt->GetXaxis()->CenterTitle(true);
	hstt->GetXaxis()->SetTitle(xtitle);
	hstt->GetYaxis()->SetTitle(ytitle);
}

void stackHist::setFillColor(){
	if( color_v ==NULL){
		color_v = new vector<Color_t> (0);
		color_v->push_back(kBlue-9);
		color_v->push_back(kYellow-9);
		color_v->push_back(kOrange+1);
		color_v->push_back(kViolet-5);
		color_v->push_back(kGreen+3);
		color_v->push_back(kRed);
		color_v->push_back(kRed+1);
		color_v->push_back(kRed+2);
		color_v->push_back(kRed+3);
		color_v->push_back(kRed+4);
	}
	if(color_v->size()<hist_trunk.size()){
		std::cout<<"color type isn't enough for the # of hists"<<std::endl;
		return;	
	}
	else {
		for(int i=0;i<hist_trunk.size();i++){
			hist_trunk.at(i)->SetFillColor(color_v->at(i));
		}
		if(doDiff){
			for(int i=0;i<hist_trunk_up.size();i++){
				hist_trunk_up.at(i)->SetFillColor(color_v->at(i));
				hist_trunk_down.at(i)->SetFillColor(color_v->at(i));
			}
		}
	}
	return;
}

void stackHist::setRange(float a, float b, TString ax){
	if(ax == "x"){
		xmin=a;
		xmax=b;
	}
	else if( ax=="y"){
		ymin=a;
		ymax=b;
	}
	else if( ax=="xd"){
		xmind=a;
		xmaxd=b;
	}
	else if( ax=="yd"){
		ymind=a;
		ymaxd=b;
	}
}

TLegend* stackHist::makeLegend(TString* tl_txt, float x1, float y1, float x2, float y2, bool isDiff, size_t n){
	TLegend* tl_stack = new TLegend(x1, y1, x2,y2);
	size_t nn;
	if( n!=0) nn = n;
	if( isDiff){
		nn = hist_trunk_up.size();
		for(size_t i=0;i<nn;i++){
			tl_stack->AddEntry(hist_trunk_up.at(i), tl_txt[i], "f");
		}
	}
	else {
		nn = hist_trunk.size();
		for(size_t i=0;i<nn;i++){
			tl_stack->AddEntry(hist_trunk.at(i), tl_txt[i], "f");
		}
	}
	tl_stack->SetLineColorAlpha(kWhite,0);
	tl_stack->SetFillColorAlpha(kWhite,0);
	tl_stack->SetTextSize(0.065);
	return tl_stack;
}

void stackHist::setTitleSize(float size, TString axis){
	if(axis == "x") xtitleSize = size;
	else ytitleSize = size;
}
void stackHist::setLabelSize(float size, TString axis){
	if(axis == "x") xlabelSize = size;
	else ylabelSize = size;
}
void stackHist::drawSum(TString opt){
	hist_sum = sumOver();
	hist_sum->SetMarkerStyle(24);
	hist_sum->SetMarkerSize(0.5);
	hist_sum->SetLineWidth(1);
	hist_sum->Draw(opt);
}

stackHist::~stackHist(){
	delete hst;
	delete hst_up;
	delete hst_down;
	hist_trunk.clear();
	hist_trunk_up.clear();
	hist_trunk_down.clear();
}

#endif
