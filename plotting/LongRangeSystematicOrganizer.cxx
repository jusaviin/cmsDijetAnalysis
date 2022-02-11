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
    } // Asymmetry loop
  } // Flow component loop
}

/*
 * Destructor
 */
LongRangeSystematicOrganizer::~LongRangeSystematicOrganizer(){
  
}

// Read the input file containing
void LongRangeSystematicOrganizer::ReadInputFile(TFile *inputFile){
  
  TString compactAsymmetryString[] = {"_A=0v0-0v6", "_A=0v6-0v8", "_A=0v8-1v0", ""};
  TString graphName;
  double errorY;
  
  const int nCentralityBins = 3;
  double summaryXaxis[nCentralityBins] = {1,2,3};
  double zeroArray[nCentralityBins] = {0,0,0};
  double smallError[nCentralityBins] = {0.1,0.1,0.1};
  
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
}

// Getter for absolute systematic uncertainty for long range correlations
TGraphErrors* LongRangeSystematicOrganizer::GetLongRangeSystematicUncertainty(const int iFlow, const int iUncertainty, int iAsymmetry) const{
  
  // If asymmetry bin is outside of the asymmetry bin range, return asymmetry integrated uncertainty
  if(iAsymmetry < 0 || iAsymmetry >= knMaxXj) iAsymmetry = knMaxXj;
  
  // Return the uncertainty in the selected bin
  return fLongRangeUncertaintyGraph[iUncertainty][iFlow][iAsymmetry];
}

// Getter for a name for the source of long range uncertainty
TString LongRangeSystematicOrganizer::GetLongRangeUncertaintyName(const int iUncertainty) const{
  return fLongRangeUncertaintyName[iUncertainty];
}

// Getter for an axis name for the source of long range uncertainty
TString LongRangeSystematicOrganizer::GetUncertaintyAxisName(const int iUncertainty) const{
  return fUncertaintyAxisName[iUncertainty];
}
