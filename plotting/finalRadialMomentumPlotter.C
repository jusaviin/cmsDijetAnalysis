#include "paperFig2Plotter.C"

void finalRadialMomentumPlotter(){
  
  std::vector<TString> vlabel;
  vlabel.push_back("Leading jet radial momentum distribution");
  auto f = TFile::Open("hin-19-013_jetShapes_systematicUpdate.root");
  auto ld_js  = loadHist_cent_trk_asym("leadingJet_js_hist", f);
  auto sub_js = loadHist_cent_trk_asym("subleadingJet_js_hist", f);
  auto sub_sum= loadHistPair_cent_asym("subleadingJet_sum_hist", "subleadingJet_sum_systUncert", f);
  auto ld_sum = loadHistPair_cent_asym("leadingJet_sum_hist", "leadingJet_sum_systUncert", f);
  
  TString format = ".pdf";
  paperFig1Plotter("figures/radialMomentum_leadingAndSubleadingCombined_enlargedTextSize"+format, ld_js, ld_sum, sub_js, sub_sum, vlabel);
  paperFig2Plotter("figures/radialMomentumRatio_leadingAndSubleadingCombined_enlargedTextSize"+format, ld_sum, sub_sum);
}
