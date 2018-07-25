#include "XiaoTrkCorr.h"

XiaoTrkCorr::XiaoTrkCorr(TString f){
	file = TFile::Open(f);
	for (int i=0; i<ncentbin; ++i)centbin[i]= i*10;
	centbin[ncentbin] = 200;
	for(int i=0; i<nptbin; ++i){
		for(int j=0;j<ncentbin; ++j){
			corrTable[i][j]=(TH2F*)file->Get(Form("corr_%d_%d",i,j));
		}
	}
}

XiaoTrkCorr::~XiaoTrkCorr()
{
  // Destructor
}

int XiaoTrkCorr::binarySearch(float key, float* arr, int i_max, int i_min ){
		if(key> arr[i_max] ) return -1;
		if(key< arr[i_min] ) return -1;
		int mid = floor(float(i_max +i_min)/2);
		//	cout<<mid<<endl;
		if(mid == i_min ) return mid;
		if( arr[mid]> key) return binarySearch(key, arr, mid, i_min);
		else if( arr[mid] < key) return binarySearch(key, arr, i_max, mid);
		else return mid;
}

double XiaoTrkCorr::getTrkCorr(float pt, float eta, float phi, float cent, float rmin, float jtpt, int correction){
	int jpt = binarySearch(pt, ptbin, nptbin,0);
	int jcent = binarySearch(cent, centbin, ncentbin,0);
	if(jpt <0 || jcent <0) {
//		std::cout<<"error!"<<std::endl;
		std::cout<<"jpt="<<jpt<<", jcent="<<jcent<<std::endl;
		std::cout<<"pt="<<pt<<", cent="<<cent<<std::endl;
		return 0;
	}
	return corrTable[jpt][jcent]->GetBinContent(corrTable[jpt][jcent]->FindBin(eta,phi));
}

