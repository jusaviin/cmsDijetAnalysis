/*
 * Implementation of LongRangeSystematicOrganizer
 */

#include "LongRangeSystematicOrganizer.h"

/*
 * Default constructor
 */
LongRangeSystematicOrganizer::LongRangeSystematicOrganizer()
{
  for(int iFlow = 0; iFlow < knMaxFlow; iFlow++){
    for(int iAsymmetry = 0; iAsymmetry <= knMaxXj; iAsymmetry++){
      for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
        fLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry] = NULL;
      } // Uncertainty loop
      for(int iUncertainty = 0; iUncertainty < knGroupedUncertaintySources; iUncertainty++){
        fGroupedLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry] = NULL;
      } // Grouped uncertainty loop
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
LongRangeSystematicOrganizer::LongRangeSystematicOrganizer(const LongRangeSystematicOrganizer& in)
{
  // Copy constructor
  
  for(int iFlow = 0; iFlow < knMaxFlow; iFlow++){
    for(int iAsymmetry = 0; iAsymmetry <= knMaxXj; iAsymmetry++){
      for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
        fLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry] = in.fLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry];
      } // Uncertainty loop
      for(int iUncertainty = 0; iUncertainty < knGroupedUncertaintySources; iUncertainty++){
        fGroupedLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry] = NULL;
      } // Grouped uncertainty loop
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
  int groupMap[knGroupedUncertaintySources][maxGroupSize+1] = {{0}};
  
  // Acceptance correction is read directly from deltaEtaSide
  groupMap[kAcceptanceCorrection][0] = 1;
  groupMap[kAcceptanceCorrection][1] = kDeltaEtaSide;
  
  // For long range extraction, combine jet-hadron and dihadron
  groupMap[kLongRangeExtraction][0] = 2;
  groupMap[kLongRangeExtraction][1] = kDeltaEtaRegion;
  groupMap[kLongRangeExtraction][2] = kDeltaEtaRegionDihadron;
  
  // For vz, do not combine anything
  groupMap[kVertexZ][0] = 1;
  groupMap[kVertexZ][1] = kVz;
  
  // For angle smearing, do not combine anything
  groupMap[kAngleSmearing][0] = 1;
  groupMap[kAngleSmearing][1] = kAngleSmear;
  
  // For jet reconstruction bias, combine all MC based correction related sources
  groupMap[kJetReconstructionBias][0] = 5;
  groupMap[kJetReconstructionBias][1] = kJetCollection;
  groupMap[kJetReconstructionBias][2] = kMCTuning;
  groupMap[kJetReconstructionBias][3] = kMCMethod;
  groupMap[kJetReconstructionBias][4] = kMCFit;
  groupMap[kJetReconstructionBias][5] = kQuarkGluonFraction;
  
  // For minimum bias, do not combine anything
  groupMap[kMinimumBias][0] = 1;
  groupMap[kMinimumBias][1] = kMinBias;
  
  // For tracking, do not combine anything
  groupMap[kTrackingGroups][0] = 1;
  groupMap[kTrackingGroups][1] = kTracking;
  
  // For jet energy correction, do not combine anything
  groupMap[kJetEnergyCorrection][0] = 1;
  groupMap[kJetEnergyCorrection][1] = kJEC;
  
  // For jet energy resolution, do not combine anything
  groupMap[kJetEnergyResolution][0] = 1;
  groupMap[kJetEnergyResolution][1] = kJER;
  
  // For total uncertainty sum, do not combine anything
  groupMap[kAllGroups][0] = 1;
  groupMap[kAllGroups][1] = kAll;
  
  // Find the number of centrality bins in the graphs
  const int nCentralityBins = fLongRangeUncertaintyGraph[kDeltaEtaSide][0][knMaxXj]->GetN();
  
  // Helper variable to add the uncertainties in quadrature
  double errorYSquareSum[nCentralityBins];
  double errorY;
  
  // Group the uncertainties together according to the grouping map
  for(int iFlow = 0; iFlow < knMaxFlow; iFlow++){
    for(int iAsymmetry = knMaxXj; iAsymmetry <= knMaxXj; iAsymmetry++){  // Asymmetry binning can be implemented easily from here
      for(int iUncertainty = 0; iUncertainty < knGroupedUncertaintySources; iUncertainty++){
        
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
  
  if(groupedUncertainty) return knGroupedUncertaintySources;
  
  return knUncertaintySources;
}
