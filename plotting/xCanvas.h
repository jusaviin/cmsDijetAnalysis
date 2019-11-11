
#ifndef xCanvas_H
#define xCanvas_H

class xCanvas : public TCanvas{
	public :
		xCanvas(TString cname, TString title, float w, float h): 
			name(cname), cw(w), ch(h),TCanvas(cname, title, w, h){
			}
		void divide(int, int);
		void CD(int);
	public :
		TString name;
		TPad **pad;
		float cw, ch;
};

void xCanvas::divide(int nx, int ny){
	float ml = this->GetLeftMargin();
	float mr = this->GetRightMargin();
	float mt = this->GetTopMargin();
	float mb = this->GetBottomMargin();
	this->SetMargin(0,0,0,0);
	float h = (1-mt-mb)/nx;
	float w = (1-ml-mr)/ny;
	float absl = cw*ml;
	float absb = ch*mb;
	pad = new TPad*[nx*ny];
	for(int i=0; i<nx; ++i){
		for(int j=0; j<ny; ++j){
			float sl=0, sb =0;
			if(j==0) sl = ml;
			if(i==nx-1) sb = mb;
			pad[i*ny+j] = new TPad(name, "", w*j+ml-sl, h*(nx-i-1)+mb-sb, ml+w*(j+1), h*(nx-i)+mb);
			if(j==0) sl = absl/(absl+w*cw);
			if(i==nx-1) sb = absb/(absb+h*ch);
			pad[i*ny+j]->SetMargin(sl,0,sb,0);
			pad[i*ny+j]->Draw();
		}
	}
}

void xCanvas::CD(int i){
	pad[i-1]->cd();
}

#endif

#ifndef auxi_canvas_H
#define auxi_canvas_H
class auxi_canvas: public TCanvas {
		public :
				auxi_canvas(TString cname, TString title, float w, float h): 
						name(cname), cw(w), ch(h),TCanvas(cname, title, w, h){
						};
				void divide(int, int);
  void divide2();
				void CD(int);
		public :
				TString name;
				TPad **pad;
				float cw, ch;
};

void auxi_canvas::divide(int nx, int ny){
  if(nx == 2 && ny == 2){
    divide2();
    return;
  }
		float ml = this->GetLeftMargin();
		float mr = this->GetRightMargin();
		float mt = this->GetTopMargin();
		float mb = this->GetBottomMargin();
		this->SetMargin(0,0,0,0);
		float h = (1-mt-mb)/nx;
		float w = (1-ml-mr)/ny;
		float absl = cw*ml;
		float absb = ch*mb;
		pad = new TPad*[nx*ny];
		for(int i=0; i<nx; ++i){
				for(int j=0; j<ny; ++j){
						float sl=0, sb =0;
						if(j==0) sl = ml;
						if(i==nx-1) sb = mb;
						pad[i*ny+j] = new TPad(name, "", w*j+ml-sl, h*(nx-i-1)+mb-sb, ml+w*(j+1), h*(nx-i)+mb);
						if(j==0) sl = absl/(absl+w*cw);
						if(i==nx-1) sb = absb/(absb+h*ch);
						pad[i*ny+j]->SetMargin(sl,0,sb,0);
						pad[i*ny+j]->Draw();
				}
		}
}

void auxi_canvas::divide2(){
  int nx = 2;
  int ny = 2;
  float ml = this->GetLeftMargin();
  float mr = this->GetRightMargin();
  float mt = this->GetTopMargin();
  float mb = this->GetBottomMargin();
  this->SetMargin(0,0,0,0);
  float h = (1-mt-mb)/nx;
  float w = (1-ml-mr)/ny;
  float absl = cw*ml;
  float absb = ch*mb;
  pad = new TPad*[nx*ny];
  for(int i=0; i<nx; ++i){
    for(int j=0; j<ny; ++j){
      float sl=0, sb =0;
      float extraMarginLeft = 0, extraMarginRight = 0, extraMarginBottom = 0, extraMarginTop = 0;
      if(j==0) sl = ml;
      if(i==nx-1) sb = mb;
      if(i == 0 && j == 0) extraMarginBottom = 0.1;
      if(i == 1 && j == 1) extraMarginLeft = 0.1;
      if(i == 1 && j == 0) {extraMarginTop = 0.1; extraMarginRight = 0.1;}
      pad[i*ny+j] = new TPad(name, "", w*j+ml-sl-extraMarginLeft, h*(nx-i-1)+mb-sb-extraMarginBottom, ml+w*(j+1)-extraMarginRight, h*(nx-i)+mb-extraMarginTop);
      if(j==0) sl = absl/(absl+w*cw);
      if(i==nx-1) sb = absb/(absb+h*ch);
      pad[i*ny+j]->SetMargin(sl+extraMarginLeft*1.83,0,sb+extraMarginBottom*1.96,0);
      pad[i*ny+j]->Draw();
    }
  }
}

void auxi_canvas::CD(int i){
		pad[i-1]->cd();
}
#endif
