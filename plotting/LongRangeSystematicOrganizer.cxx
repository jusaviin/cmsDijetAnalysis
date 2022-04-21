/*
 * Implementation of LongRangeSystematicOrganizer
 */

#include "LongRangeSystematicOrganizer.h"

/*
 * Default constructor
 */
LongRangeSystematicOrganizer::LongRangeSystematicOrganizer() :
  fnGroupedUncertaintySources(knUncertaintySources),
  fGroupingStrategy(0)
{
  for(int iFlow = 0; iFlow < knMaxFlow; iFlow++){
    for(int iAsymmetry = 0; iAsymmetry <= knMaxXj; iAsymmetry++){
      for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
        fLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry] = NULL;
        fGroupedLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry] = NULL;
      } // Uncertainty loop
    } // Asymmetry loop
  } // Flow component loop
}

/*
 * Constructor
 */
LongRangeSystematicOrganizer::LongRangeSystematicOrganizer(TFile *inputFile) :
  LongRangeSystematicOrganizer()
{
  ReadInputFile(inputFile);
}

/*
 * Copy constructor
 */
LongRangeSystematicOrganizer::LongRangeSystematicOrganizer(const LongRangeSystematicOrganizer& in) :
  fnGroupedUncertaintySources(in.fnGroupedUncertaintySources),
  fGroupingStrategy(in.fGroupingStrategy)
{
  // Copy constructor
  
  for(int iFlow = 0; iFlow < knMaxFlow; iFlow++){
    for(int iAsymmetry = 0; iAsymmetry <= knMaxXj; iAsymmetry++){
      for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
        fLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry] = in.fLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry];
        fGroupedLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry] = in.fGroupedLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry] ;
      } // Uncertainty loop
    } // Asymmetry loop
  } // Flow component loop
}

/*
 * Destructor
 */
LongRangeSystematicOrganizer::~LongRangeSystematicOrganizer(){
  
}

// Read the input file containing the uncertainty histograms
void LongRangeSystematicOrganizer::ReadInputFile(TFile *inputFile){
  
  TString compactAsymmetryString[] = {"_A=0v0-0v6", "_A=0v6-0v8", "_A=0v8-1v0", ""};
  TString graphName;
  double errorY;
  
  const int nCentralityBins = 3;
  double summaryXaxis[nCentralityBins] = {1,2,3};
  double zeroArray[nCentralityBins] = {0,0,0};
  double smallError[nCentralityBins] = {0.1,0.1,0.1};
  
  // Read the uncertainties from the file
  for(int iFlow = 0; iFlow < knMaxFlow; iFlow++){
    for(int iAsymmetry = knMaxXj; iAsymmetry <= knMaxXj; iAsymmetry++){  // Asymmetry binning can be implemented easily from here
      for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
        
        graphName = Form("systematicUncertainty_v%d_%s%s", iFlow+1, fLongRangeUncertaintyName[iUncertainty].Data(), compactAsymmetryString[iAsymmetry].Data());
        fLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry] = (TGraphErrors*) inputFile->Get(graphName);
        
        // Set some width for the x-axis uncertainty
        if(fLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry] != NULL){
          for(int iCentrality = 0; iCentrality < fLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry]->GetN(); iCentrality++){
            errorY = fLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry]->GetErrorY(iCentrality);
            fLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry]->SetPointError(iCentrality, 0.1, errorY);
          }
        }  else {// If the graph does not exist, set the uncertainty to zero
          
          fLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry] = new TGraphErrors(nCentralityBins, summaryXaxis, zeroArray, smallError, zeroArray);
          
        }
        
      } // Uncertainty loop
    } // Asymmetry loop
  } // Flow component loop
  
  // Create groups from the individual uncertainty histograms
  GroupUncertaintyHistograms();
  
}

// Read the input file containing the uncertainty histograms
void LongRangeSystematicOrganizer::GroupUncertaintyHistograms(){
  
  // Define the grouping
  const int maxGroupSize = 5;
  int groupMap[knUncertaintySources][maxGroupSize+1] = {{0}};
  
  if(fGroupingStrategy == 0){
    // Group similar sources together. Good for producing tables for paper
    
    fnGroupedUncertaintySources = 10; // In this grouping scheme, there will be 10 different groups
    
    // Acceptance correction is read directly from deltaEtaSide
    fGroupedUncertaintyName[0] = "Acceptance correction";
    fGroupedUncertaintyAxisName[0] = "Acceptance";
    groupMap[0][0] = 1;
    groupMap[0][1] = kDeltaEtaSide;
    
    // For long range extraction, combine jet-hadron and dihadron
    fGroupedUncertaintyName[1] = "Long range extraction";
    fGroupedUncertaintyAxisName[1] = "Long range";
    groupMap[1][0] = 2;
    groupMap[1][1] = kDeltaEtaRegion;
    groupMap[1][2] = kDeltaEtaRegionDihadron;
    
    // For vz, do not combine anything
    fGroupedUncertaintyName[2] = "vz";
    fGroupedUncertaintyAxisName[2] = "v_{z}";
    groupMap[2][0] = 1;
    groupMap[2][1] = kVz;
    
    // For angle smearing, do not combine anything
    fGroupedUncertaintyName[3] = "Jet axis resolution";
    fGroupedUncertaintyAxisName[3] = "Angle smear";
    groupMap[3][0] = 1;
    groupMap[3][1] = kAngleSmear;
    
    // For jet reconstruction bias, combine all MC based correction related sources
    fGroupedUncertaintyName[4] = "Jet reconstruction bias";
    fGroupedUncertaintyAxisName[4] = "Jet reco bias";
    groupMap[4][0] = 5;
    groupMap[4][1] = kJetCollection;
    groupMap[4][2] = kMCTuning;
    groupMap[4][3] = kMCMethod;
    groupMap[4][4] = kMCFit;
    groupMap[4][5] = kQuarkGluonFraction;
    
    // For minimum bias, do not combine anything
    fGroupedUncertaintyName[5] = "Dijet bias for dihadron";
    fGroupedUncertaintyAxisName[5] = "MinBias";
    groupMap[5][0] = 1;
    groupMap[5][1] = kMinBias;
    
    // For tracking, do not combine anything
    fGroupedUncertaintyName[6] = "Tracking";
    fGroupedUncertaintyAxisName[6] = "Tracking";
    groupMap[6][0] = 1;
    groupMap[6][1] = kTracking;
    
    // For jet energy correction, do not combine anything
    fGroupedUncertaintyName[7] = "JEC";
    fGroupedUncertaintyAxisName[7] = "JEC";
    groupMap[7][0] = 1;
    groupMap[7][1] = kJEC;
    
    // For jet energy resolution, do not combine anything
    fGroupedUncertaintyName[8] = "JER";
    fGroupedUncertaintyAxisName[8] = "JER";
    groupMap[8][0] = 1;
    groupMap[8][1] = kJER;
    
    // For total uncertainty sum, do not combine anything
    fGroupedUncertaintyName[9] = "Total";
    fGroupedUncertaintyAxisName[9] = "Total";
    groupMap[9][0] = 1;
    groupMap[9][1] = kAll;
    
  } else {
    // Group things such that old obsolete sources are removed. Good for suppressing zero histograms if closer source inspection if needed
    
    fnGroupedUncertaintySources = 12; // In this grouping scheme, there will be 12 different groups
    
    // Acceptance correction is read directly from deltaEtaSide
    fGroupedUncertaintyName[0] = "deltaEtaSide";
    fGroupedUncertaintyAxisName[0] = "#Delta#eta side";
    groupMap[0][0] = 1;
    groupMap[0][1] = kDeltaEtaSide;
    
    // Combine obsolete vz to deltaEtaRegion to hide it
    fGroupedUncertaintyName[1] = "deltaEtaRegion";
    fGroupedUncertaintyAxisName[1] = "#Delta#eta region";
    groupMap[1][0] = 2;
    groupMap[1][1] = kDeltaEtaRegion;
    groupMap[1][2] = kVz;
    
    // Do not combine anything with hadron-hadron deltaEta region
    fGroupedUncertaintyName[2] = "deltaEtaRegionDihadron";
    fGroupedUncertaintyAxisName[2] = "#Delta#eta region hh";
    groupMap[2][0] = 1;
    groupMap[2][1] = kDeltaEtaRegionDihadron;
    
    // For angle smearing, do not combine anything
    fGroupedUncertaintyName[3] = "Jet axis resolution";
    fGroupedUncertaintyAxisName[3] = "Angle smear";
    groupMap[3][0] = 1;
    groupMap[3][1] = kAngleSmear;
    
    // Combine centrality shift and jet collection uncertainties together to MCMethod
    fGroupedUncertaintyName[4] = "mcMethod";
    fGroupedUncertaintyAxisName[4] = "MC method";
    groupMap[4][0] = 3;
    groupMap[4][1] = kJetCollection;
    groupMap[4][2] = kMCTuning;
    groupMap[4][3] = kMCMethod;
    
    // For MC fit, do not combine anything
    fGroupedUncertaintyName[5] = "mcFit";
    fGroupedUncertaintyAxisName[5] = "MC Fit";
    groupMap[5][0] = 1;
    groupMap[5][1] = kMCFit;
    
    // For quark/gluon fraction, do not combine anything
    fGroupedUncertaintyName[6] = "quarkGluonFraction";
    fGroupedUncertaintyAxisName[6] = "q/g fraction";
    groupMap[6][0] = 1;
    groupMap[6][1] = kTracking;
    
    // For minimum bias, do not combine anything
    fGroupedUncertaintyName[7] = "Dijet bias for dihadron";
    fGroupedUncertaintyAxisName[7] = "MinBias";
    groupMap[7][0] = 1;
    groupMap[7][1] = kMinBias;
    
    // For tracking, do not combine anything
    fGroupedUncertaintyName[8] = "Tracking";
    fGroupedUncertaintyAxisName[8] = "Tracking";
    groupMap[8][0] = 1;
    groupMap[8][1] = kTracking;
    
    // For jet energy correction, do not combine anything
    fGroupedUncertaintyName[9] = "JEC";
    fGroupedUncertaintyAxisName[9] = "JEC";
    groupMap[9][0] = 1;
    groupMap[9][1] = kJEC;
    
    // For jet energy resolution, do not combine anything
    fGroupedUncertaintyName[10] = "JER";
    fGroupedUncertaintyAxisName[10] = "JER";
    groupMap[10][0] = 1;
    groupMap[10][1] = kJER;
    
    // For total uncertainty sum, do not combine anything
    fGroupedUncertaintyName[11] = "Total";
    fGroupedUncertaintyAxisName[11] = "Total";
    groupMap[11][0] = 1;
    groupMap[11][1] = kAll;
    
  }
  
  // Find the number of centrality bins in the graphs
  const int nCentralityBins = fLongRangeUncertaintyGraph[kDeltaEtaSide][0][knMaxXj]->GetN();
  
  // Helper variable to add the uncertainties in quadrature
  double errorYSquareSum[nCentralityBins];
  double errorY;
  
  // Group the uncertainties together according to the grouping map
  for(int iFlow = 0; iFlow < knMaxFlow; iFlow++){
    for(int iAsymmetry = knMaxXj; iAsymmetry <= knMaxXj; iAsymmetry++){  // Asymmetry binning can be implemented easily from here
      for(int iUncertainty = 0; iUncertainty < fnGroupedUncertaintySources; iUncertainty++){
        
        // Reset the error sum
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          errorYSquareSum[iCentrality] = 0;
        }
        
        // First, clone the regular uncertainty graph as a starting point of the grouped graph
        fGroupedLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry] = (TGraphErrors*) fLongRangeUncertaintyGraph[groupMap[iUncertainty][1]][iFlow][iAsymmetry]->Clone(Form("groupedUncertaintyGraph%d%d%d", iFlow, iAsymmetry, iUncertainty));
        
        // Calculate the error sum from all the sources described in the grouping map
        for(int iSource = 1; iSource <= groupMap[iUncertainty][0]; iSource++){
          for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
            errorY = fLongRangeUncertaintyGraph[groupMap[iUncertainty][iSource]][iFlow][iAsymmetry]->GetErrorY(iCentrality);
            errorYSquareSum[iCentrality] += errorY*errorY;
          } // Centrlaity loop
        } // Grouped uncertainty sources loop
          
        // After all the sources are summed, set the y-error as the quadrature sum of the error sources
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          fGroupedLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry]->SetPointError(iCentrality, 0.1, TMath::Sqrt(errorYSquareSum[iCentrality]));
        }
        
      } // Uncertainty loop
    } // Asymmetry loop
  } // Flow component loop
  
}

// Getter for absolute systematic uncertainty for long range correlations
TGraphErrors* LongRangeSystematicOrganizer::GetLongRangeSystematicUncertainty(const int iFlow, const int iUncertainty, int iAsymmetry, const bool groupedUncertainty) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated uncertainty
  if(iAsymmetry < 0 || iAsymmetry >= knMaxXj) iAsymmetry = knMaxXj;
  
  // If grouped uncertainty is asked, return that
  if(groupedUncertainty) return fGroupedLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry];
  
  // Return the uncertainty in the selected bin
  return fLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry];
}

// Getter for absolute systematic uncertainty value for long range correlations
double LongRangeSystematicOrganizer::GetLongRangeSystematicUncertaintyValue(const int iFlow, const int iCentrality, const int iUncertainty, int iAsymmetry) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated uncertainty
  if(iAsymmetry < 0 || iAsymmetry >= knMaxXj) iAsymmetry = knMaxXj;
  
  // Return the uncertainty in the selected bin
  return fLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry]->GetErrorY(iCentrality);
}

// Getter for a name for the source of long range uncertainty
TString LongRangeSystematicOrganizer::GetLongRangeUncertaintyName(const int iUncertainty, const bool groupedUncertainty) const{
  
  // If grouped uncertainty name is asked, return that
  if(groupedUncertainty) return fGroupedUncertaintyName[iUncertainty];
  
  return fLongRangeUncertaintyName[iUncertainty];
}

// Getter for an axis name for the source of long range uncertainty
TString LongRangeSystematicOrganizer::GetUncertaintyAxisName(const int iUncertainty, const bool groupedUncertainty) const{
  
  // If grouped uncertainty axis name is asked, return that
  if(groupedUncertainty) return fGroupedUncertaintyAxisName[iUncertainty];
  
  return fUncertaintyAxisName[iUncertainty];
}

// Getter for an axis name for the source of long range uncertainty
int LongRangeSystematicOrganizer::GetNUncertaintySources(const bool groupedUncertainty) const{
  
  if(groupedUncertainty) return fnGroupedUncertaintySources;
  
  return knUncertaintySources;
}

// Define how to do histogram grouping and perform it
void LongRangeSystematicOrganizer::SetGroupingStrategy(const int groupingStrategy){
  
  fGroupingStrategy = groupingStrategy;
  GroupUncertaintyHistograms();
  
}

// Adjust the central points for uncertainty graphs according to an input graph
void LongRangeSystematicOrganizer::AdjustCentralPoints(TGraphErrors* resultGraph[knMaxFlow]){
  
  int nPoints;
  double resultX, resultY;
  
  for(int iFlow = 0; iFlow < knMaxFlow; iFlow++){
    if(resultGraph[iFlow] == NULL) continue;
    nPoints = resultGraph[iFlow]->GetN();
    for(int iAsymmetry = knMaxXj; iAsymmetry <= knMaxXj; iAsymmetry++){  // Asymmetry binning can be implemented easily from here
      for(int iPoint = 0; iPoint < nPoints; iPoint++){
        resultGraph[iFlow]->GetPoint(iPoint, resultX, resultY);
        for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
          fLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry]->SetPointY(iPoint, resultY);
        } // Uncertainty loop
        for(int iUncertainty = 0; iUncertainty < fnGroupedUncertaintySources; iUncertainty++){
          fGroupedLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry]->SetPointY(iPoint, resultY);
        } // Grouped uncertainty loop
      } // Point loop
    } // Asymmetry loop
  } // Flow component loop
  
}
