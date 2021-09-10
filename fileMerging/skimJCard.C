void skimJCard(TString fileName, int nEntries){
  TFile *skimmedFile = new TFile(fileName,"UPDATE");
  skimmedFile->cd("JCard");
  char vectorName[100];
  const int nCardEntries = 45;
  const char *cardEntries[nCardEntries] = {"DataType","McCorrelationType","MatchJets","ForestType","ReadMode","FilledHistograms","JetType","JetAxis","JetEtaCut","SearchEtaCut","MaxPtCut","MinPtCut","SubleadingPtCut","DeltaPhiCut","MinMaxTrackPtFraction","MaxMaxTrackPtFraction","JetUncertainty","TrackEtaCut","MinTrackPtCut","MaxTrackPtRelativeError","VertexMaxDistance","CalorimeterSignalLimitPt","HighPtEtFraction","Chi2QualityCut","MinimumTrackHits","SubeventCut","ZVertexCut","LowPtHatCut","HighPtHatCut","IncludeEventPlane","AsymmetryBinType","CentralityBinEdges","TrackPtBinEdges","AsymmetryBinEdges","PtHatBinEdges","DoEventMixing","MixWithPool","NMixedEventsPerDijet","VzTolerance","MixingVzBinWidth","MixingHiBinWidth","MixingPoolDepth","MixingFileIndex","DebugLevel","OnlyMix"};
  for(int iVector = 0; iVector < nCardEntries; iVector++){
    for(int j = 2; j <= nEntries; j++){
      sprintf(vectorName,"%s;%d",cardEntries[iVector],j);
      gDirectory->Delete(vectorName);
    }
  }
  skimmedFile->Close();
}
