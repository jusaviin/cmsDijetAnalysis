
#ifndef ROOT_TH2F
#include  "TH2F.h"
#endif

#ifndef ROOT_TFile
#include "TFile.h"
#endif

#include <iostream>
#include "TrkCorrInterface.h"

class XiaoTrkCorr : public TrkCorrInterface{
	public: 
		XiaoTrkCorr(TString f);
        virtual ~XiaoTrkCorr();
		double getTrkCorr(float pt, float eta, float phi, float hibin, float rmin=99, float jtpt=0, int correction=0);
	private:
        int binarySearch(float key, float* arr, int i_max, int i_min );
		TFile* file;
		TH2F* corrTable[22][17];
		int nptbin = 22;
		int ncentbin = 17;
		float ptbin[23] = {0.7,0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2,
			2.5, 3, 3.5, 4,5, 6, 7, 8, 10, 12, 16, 20, 50, 999};
		float centbin[18];
};


