#
# Macro for creating HepData entry for dijet analysis
#
# Code is writted based on examples in https://github.com/HEPData/hepdata_lib/blob/master/examples/
#

# Use hepdata_lib package for easier data handling
import hepdata_lib

# First, create a HepData submission object
# Later all the tables from different figures are added to this object
from hepdata_lib import Submission
submission = Submission()

# Add basic paper info to the submission object
submission.read_abstract("abstract.txt")
submission.add_link("Webpage with all figures and tables", "https://cms-results.web.cern.ch/cms-results/public-results/publications/HIN-19-013/")
submission.add_link("arXiv", "http://arxiv.org/abs/arXiv:20XX.XXXXXX")
submission.add_record_id(1010101, "inspire")

# All the figures are inserted as table objects to HepData
from hepdata_lib import Table

# Read the root files with the histograms
from hepdata_lib import RootFileReader

##############################################################################
#      Read the root files containing final results from dijet analysis      #
##############################################################################

# Common naming convention in the final results files
centralityString = ["0-10", "10-30", "30-50", "50-90", "pp"]
xjString = ["0<xj<06", "06<xj<08", "08<xj<1", "allXj"]
asymmetryString = ["Imbalanced", "Balanced"]
trackPtString = ["","_07-1","_1-2","_2-3","_3-4","_4-8","_8-12","_12-300"]
centralityDescription = ["the 0-10 % centrality bin in PbPb", "the 10-30 % centrality bin in PbPb", "the 30-50 % centrality bin in PbPb", "the 50-90 % centrality bin in PbPb", "pp"]
trackPtDescription = ["", " for the charged particle $p_{T}$ bin $0.7 < p_{T}^{ch} <1$ GeV", " for the charged particle $p_{T}$ bin $1 < p_{T}^{ch} <2$ GeV", " for the charged particle $p_{T}$ bin $2 < p_{T}^{ch} <3$ GeV", " for the charged particle $p_{T}$ bin $3 < p_{T}^{ch} <4$ GeV", " for the charged particle $p_{T}$ bin $4 < p_{T}^{ch} <8$ GeV", " for the charged particle $p_{T}$ bin $8 < p_{T}^{ch} <12$ GeV", " for the charged particle $p_{T}$ bin $12 < p_{T}^{ch} <300$ GeV"]

# Numbers of bins
nCentrality = 5       # Number of centrality bins, including pp
nAsymmetry = 4        # Number of xj bins, including xj integrated
nTrackPt = 7          # Number of track pT bins for jet shapes
nTrackPtDeltaEta = 6  # Number of track pT bins for deltaEta yields

# File for particle yields as a function of deltaEta For paper figures 1 and 2.
deltaEtaReader = RootFileReader("/Users/jviinika/cms/development/dijet5TeV/hepdata/hepdata_deltaEta_update_hin-19-013.root")

# Read the deltaEta yield histograms from the root file
deltaEtaLeading = []
deltaEtaErrorLeading = []
deltaEtaSubleading = []
deltaEtaErrorSubleading = []

for iTrackPt in trackPtString:
    # Skip the highest bin for yields:
    if iTrackPt == "12-300":
        continue
    for iCentrality in centralityString:
        for iXj in xjString:
            deltaEtaLeading.append(deltaEtaReader.read_hist_1d("deltaEta_trackLeadingJet_" + iCentrality + iTrackPt + "_" + iXj))
            deltaEtaErrorLeading.append(deltaEtaReader.read_hist_1d("deltaEtaError_trackLeadingJet_" + iCentrality + iTrackPt + "_" + iXj))
            deltaEtaSubleading.append(deltaEtaReader.read_hist_1d("deltaEta_trackSubleadingJet_" + iCentrality + iTrackPt + "_" + iXj))
            deltaEtaErrorSubleading.append(deltaEtaReader.read_hist_1d("deltaEtaError_trackSubleadingJet_" + iCentrality + iTrackPt + "_" + iXj))

# File for jet radial momentum distributions. For paper figures 3 and 4.
jetMomentumReader = RootFileReader("/Users/jviinika/cms/development/dijet5TeV/hepdata/hepdata_jetRadialMomentum_update_hin-19-013.root")

# Read the jet radial momentum histograms from the root file
jetMomentumLeading = []
jetMomentumErrorLeading = []
jetMomentumRatioLeading = []
jetMomentumRatioErrorLeading = []
jetMomentumSubleading = []
jetMomentumErrorSubleading = []
jetMomentumRatioSubleading = []
jetMomentumRatioErrorSubleading = []

for iTrackPt in trackPtString:
    for iCentrality in centralityString:
        jetMomentumLeading.append(jetMomentumReader.read_hist_1d("jetRadialMomentum_trackLeadingJet_" + iCentrality + iTrackPt + "_allXj"))
        jetMomentumErrorLeading.append(jetMomentumReader.read_hist_1d("jetRadialMomentumError_trackLeadingJet_" + iCentrality + iTrackPt + "_allXj"))
        jetMomentumSubleading.append(jetMomentumReader.read_hist_1d("jetRadialMomentum_trackSubleadingJet_" + iCentrality + iTrackPt + "_allXj"))
        jetMomentumErrorSubleading.append(jetMomentumReader.read_hist_1d("jetRadialMomentumError_trackSubleadingJet_" + iCentrality + iTrackPt + "_allXj"))
        
        # As the ratios are between PbPb and pp, there is one less centrality bin compared to just distributions
        # Also the ratio is shown only for pT summed distributions
        if iCentrality != "pp" and iTrackPt == "":
            jetMomentumRatioLeading.append(jetMomentumReader.read_hist_1d("jetRadialMomentumRatio_trackLeadingJet_" + iCentrality + iTrackPt + "_allXj"))
            jetMomentumRatioErrorLeading.append(jetMomentumReader.read_hist_1d("jetRadialMomentumRatioError_trackLeadingJet_" + iCentrality + iTrackPt + "_allXj"))
            jetMomentumRatioSubleading.append(jetMomentumReader.read_hist_1d("jetRadialMomentumRatio_trackSubleadingJet_" + iCentrality + iTrackPt + "_allXj"))
            jetMomentumRatioErrorSubleading.append(jetMomentumReader.read_hist_1d("jetRadialMomentumRatioError_trackSubleadingJet_" + iCentrality + iTrackPt + "_allXj"))
        
# File for jet shapes. For paper figures 5, 6, 7 and 8.
jetShapeReader = RootFileReader("/Users/jviinika/cms/development/dijet5TeV/hepdata/hepdata_jetShapes_update_hin-19-013.root")

# Read the jet shape histograms from the root file
jetShapeLeading = []
jetShapeErrorLeading = []
jetShapeRatioLeading = []
jetShapeRatioErrorLeading = []
jetShapeSubleading = []
jetShapeErrorSubleading = []
jetShapeRatioSubleading = []
jetShapeRatioErrorSubleading = []

for iTrackPt in trackPtString:
    for iCentrality in centralityString:
        for iXj in xjString:
            jetShapeLeading.append(jetShapeReader.read_hist_1d("jetShape_trackLeadingJet_" + iCentrality + iTrackPt + "_" + iXj))
            jetShapeErrorLeading.append(jetShapeReader.read_hist_1d("jetShapeError_trackLeadingJet_" + iCentrality + iTrackPt + "_" + iXj))
            jetShapeSubleading.append(jetShapeReader.read_hist_1d("jetShape_trackSubleadingJet_" + iCentrality + iTrackPt + "_" + iXj))
            jetShapeErrorSubleading.append(jetShapeReader.read_hist_1d("jetShapeError_trackSubleadingJet_" + iCentrality + iTrackPt + "_" + iXj))
        
            # As the ratios are between PbPb and pp, there are one less of the ocmpared to just distributions
            if iCentrality != "pp" and iTrackPt == "":
                jetShapeRatioLeading.append(jetShapeReader.read_hist_1d("jetShapeRatio_trackLeadingJet_" + iCentrality + "_" + iXj))
                jetShapeRatioErrorLeading.append(jetShapeReader.read_hist_1d("jetShapeRatioError_trackLeadingJet_" + iCentrality + "_" + iXj))
                jetShapeRatioSubleading.append(jetShapeReader.read_hist_1d("jetShapeRatio_trackSubleadingJet_" + iCentrality + "_" + iXj))
                jetShapeRatioErrorSubleading.append(jetShapeReader.read_hist_1d("jetShapeRatioError_trackSubleadingJet_" + iCentrality + "_" + iXj))

# File for asymmetry ratios. For paper figures 9 and 10.
asymmetryRatioReader = RootFileReader("/Users/jviinika/cms/development/dijet5TeV/hepdata/hepdata_asymmetryRatio_update_hin-19-013.root")

# Read the asymmetry ratio histograms from the root file
asymmetryRatioLeading = []
asymmetryRatioErrorLeading = []
asymmetryRatioSubleading = []
asymmetryRatioErrorSubleading = []

for iCentrality in centralityString:
    for iAsymmetry in asymmetryString:
        asymmetryRatioLeading.append(asymmetryRatioReader.read_hist_1d("asymmetryRatio" + iAsymmetry + "_trackLeadingJet_" + iCentrality)) 
        asymmetryRatioErrorLeading.append(asymmetryRatioReader.read_hist_1d("asymmetryRatioError" + iAsymmetry + "_trackLeadingJet_" + iCentrality)) 
        asymmetryRatioSubleading.append(asymmetryRatioReader.read_hist_1d("asymmetryRatio" + iAsymmetry + "_trackSubleadingJet_" + iCentrality)) 
        asymmetryRatioErrorSubleading.append(asymmetryRatioReader.read_hist_1d("asymmetryRatioError" + iAsymmetry + "_trackSubleadingJet_" + iCentrality)) 

# File for xj matrices. For paper figures 11 and 12.
xjMatrixReader = RootFileReader("/Users/jviinika/cms/development/dijet5TeV/hepdata/xjMatrix_hepdata.root")

# Read the xj matrices from the root file
xjMatrix = []
xjMatrixReverse = []

for iCentrality in centralityString:
    xjMatrix.append(xjMatrixReader.read_hist_2d("xjMatrix_" + iCentrality))
    xjMatrixReverse.append(xjMatrixReader.read_hist_2d("xjMatrix_reverse_" + iCentrality))


# Read the variables from the histograms
from hepdata_lib import Variable, Uncertainty

# Function for finding x and y-axis variables from input histograms when the histograms are binned in deltaR
#
# Arguments:
#  valueHistogram: Histogram containing the actual values and statistical errors
#  errorHistogram: Histogram containing the systematic errors
#  yAxisName: Label given to the histogram y-axis
#  includeAsymmetry: True = Make data for all asymmetry bins. False = Only use one asymmetry bin
#  centralityBins: Tuple where first index gives the number of centrality bins and the second the selected one, if number is 1
#  deltaEtaMode: True = Assume deltaEta histogram. False = Assume deltaR histogram
#  trackPtBin: Given track pT bin. Negative value implies that pT integrated is used instead
#
#  return: Returns HepData variables for x and y-axes from the input histograms
#
def findVariables(valueHistogram, errorHistogram, yAxisName, includeAsymmetry, centralityBins, deltaEtaMode, trackPtBin = -1):

    # Define reaction and centrality labels.
    reactionLabel = ["PB PB --> CHARGED X", "PB PB --> CHARGED X", "PB PB --> CHARGED X", "PB PB --> CHARGED X", "P P --> CHARGED X"]
    centralityLabel = ["0-10%", "10-30%", "30-50%", "50-90%", "pp"]
    asymmetryLabel = ["$0 < x_{j} < 0.6$", "$0.6 < x_{j} < 0.8$", "$0.8 < x_{j} < 1.0$", "All $x_{j}$"]
    trackPtLabel = ["$0.7 < p_{T}^{ch} < 300$ GeV", "$0.7 < p_{T}^{ch} < 1$ GeV", "$1 < p_{T}^{ch} < 2$ GeV", "$2 < p_{T}^{ch} < 3$ GeV", "$3 < p_{T}^{ch} < 4$ GeV", "$4 < p_{T}^{ch} < 8$ GeV", "$8 < p_{T}^{ch} < 12$ GeV", "$12 < p_{T}^{ch} < 300$ GeV"]
    
    # Figure out the asymmetry bin range
    nAsymmetry = 1
    if includeAsymmetry:
        nAsymmetry = 4
        
    # Figure out the centrality bin range
    firstCentralityBin = centralityBins[1]
    lastCentralityBin = centralityBins[1]+1
    
    if centralityBins[0] > 1:
        firstCentralityBin = 0
        lastCentralityBin = centralityBins[0]
        
    if deltaEtaMode:
        trackPtLabel[0] = "$0.7 < p_{T}^{ch} < 12$ GeV"
    
    # x-axis: deltaR value
    variableName = "$\Delta r$"
    if deltaEtaMode:
        variableName = "$\Delta\eta$"
    
    xAxis = Variable(variableName, is_independent=True, is_binned=True, units="")
    xAxis.values = valueHistogram[0]["x_edges"]

    # Determine how many bins in the histogram is outside of the analysis range of [0.1]
    overRange = 0
    cutValue = 1.1
    if deltaEtaMode:
        cutValue = 1.6
      
    for (i,j) in xAxis.values:
        if j > cutValue:
            overRange = overRange + 1
                
    # For deltaEta histograms, we want to remove the negative points from the beginning
    underRange = 0
    if deltaEtaMode:
        for (i,j) in xAxis.values:
            if j < 0.01:
                underRange = underRange + 1
            
    # Remove the points that are outside of the range
    for i in range(overRange):
        xAxis.values.pop()
        
    for i in range(underRange):
        xAxis.values.pop(0)
        
    # Define a list for y-axis values
    yAxis = []

    # Read the values
    for i in range(firstCentralityBin,lastCentralityBin):
        for iXj in range(nAsymmetry):
            myVariable = Variable(yAxisName, is_independent=False, is_binned=False, units="")
            myVariable.values = valueHistogram[(i-firstCentralityBin)*nAsymmetry+iXj]["y"]
            myVariable.add_qualifier("reaction",reactionLabel[i])
            myVariable.add_qualifier("centrality",centralityLabel[i])
            myVariable.add_qualifier("$p_{T}^{ch}$",trackPtLabel[trackPtBin+1])
            if includeAsymmetry:
                myVariable.add_qualifier("$x_{j}$",asymmetryLabel[iXj])
        
            # If there are points outside of the analysis range, remove them from the values table:
            for j in range(overRange):
                myVariable.values.pop()
                
            for j in range(underRange):
                myVariable.values.pop(0)
        
            yAxis.append(myVariable)
        
    # Read the uncertainties
    for i in range(firstCentralityBin,lastCentralityBin):
        for iXj in range(nAsymmetry):
            statUncertainty = Uncertainty("stat", is_symmetric=True)
            statUncertainty.values = valueHistogram[(i-firstCentralityBin)*nAsymmetry+iXj]["dy"]

            systUncertainty = Uncertainty("sys", is_symmetric=True)
            systUncertainty.values = errorHistogram[(i-firstCentralityBin)*nAsymmetry+iXj]["dy"]
        
            # If there are points outside of the analysis range, remove their errors from the table:
            for j in range(overRange):
                statUncertainty.values.pop()
                systUncertainty.values.pop()
                
            for j in range(underRange):
                statUncertainty.values.pop(0)
                systUncertainty.values.pop(0)

            yAxis[(i-firstCentralityBin)*nAsymmetry+iXj].add_uncertainty(statUncertainty)
            yAxis[(i-firstCentralityBin)*nAsymmetry+iXj].add_uncertainty(systUncertainty)
        
    return xAxis, yAxis

# Function for finding x, y and z-axis variables from input histograms when the histograms are binned in deltaR
#
# Arguments:
#  valueHistogram: Histogram containing the actual values and statistical errors
#
#  return: Returns HepData variables for x and y-axes from the input histograms
#
def findVariablesXjMatrix(valueHistogram):

    # Define reaction and centrality labels.
    reactionLabel = ["PB PB --> CHARGED X", "PB PB --> CHARGED X", "PB PB --> CHARGED X", "PB PB --> CHARGED X", "P P --> CHARGED X"]
    centralityLabel = ["0-10%", "10-30%", "30-50%", "50-90%", "pp"]
    
    # x-axis:
    xAxis = Variable("Reconstructed $x_{j}$", is_independent=True, is_binned=True, units="")
    xAxis.values = valueHistogram[0]["x_edges"]
    
    # y-axis:
    yAxis = Variable("Generator level $x_{j}$", is_independent=True, is_binned=True, units="")
    yAxis.values = valueHistogram[0]["y_edges"]
    
    # Define a list for z-axis values
    zAxis = []
    
    # Read the values
    for iCentrality in range(0,len(centralityLabel)):
        myVariable = Variable("Probability", is_independent=False, is_binned=False, units="")
        myVariable.values = valueHistogram[iCentrality]["z"]
        myVariable.add_qualifier("reaction",reactionLabel[iCentrality])
        myVariable.add_qualifier("centrality",centralityLabel[iCentrality])
        zAxis.append(myVariable)
        
    # Read the uncertainties
    for iCentrality in range(0,len(centralityLabel)):
        statUncertainty = Uncertainty("stat", is_symmetric=True)
        statUncertainty.values = valueHistogram[iCentrality]["dz"]
        zAxis[iCentrality].add_uncertainty(statUncertainty)
        
    # Return extracted vatiables
    return xAxis, yAxis, zAxis

######################################################################################
#               Table for leading jet deltaEta yield. Paper figure 1.                #
######################################################################################

# Make separate tables from each centrality and track pT bin
for iTrackPt in range(-1,nTrackPtDeltaEta):
    for iCentrality in range(0,nCentrality):
        table1 = Table("Figure 1-{:d}".format(iCentrality+nCentrality*(iTrackPt+1)))
        table1.description = "The distribution of charged particle yields within $|\Delta\\varphi| < 1.0$ correlated with the leading jets as a function of $\Delta\eta$ in {:s} collisions. The results are shown in different dijet momentum balance bins{:s}.".format(centralityDescription[iCentrality],trackPtDescription[iTrackPt+1])
        table1.location = "Data from figure 1, located on page 9."
        table1.keywords["observables"] = ["Particle yield correlated to leading jets"]
        #table1.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

        # Extract x- and y-axis information from the histograms
        xDeltaEta, yDeltaEta = findVariables(deltaEtaLeading[nAsymmetry*iCentrality+(iTrackPt+1)*nAsymmetry*nCentrality : nAsymmetry+nAsymmetry*iCentrality+(iTrackPt+1)*nAsymmetry*nCentrality], deltaEtaErrorLeading[nAsymmetry*iCentrality+(iTrackPt+1)*nAsymmetry*nCentrality : nAsymmetry+nAsymmetry*iCentrality+(iTrackPt+1)*nAsymmetry*nCentrality], "$\\frac{1}{N_{\mathrm{dijet}}}\\frac{\mathrm{d}N}{\mathrm{d}\Delta\eta}$", True, (1,iCentrality), True, iTrackPt)

        # Add the variables to the table
        table1.add_variable(xDeltaEta)
        for variable in yDeltaEta:
            table1.add_variable(variable)

        # Add the table to submission object
        submission.add_table(table1)

######################################################################################
#             Table for subleading jet deltaEta yield. Paper figure 2.               #
######################################################################################

# Make separate tables from each centrality and track pT bin
for iTrackPt in range(-1,nTrackPtDeltaEta):
    for iCentrality in range(0,nCentrality):
        table2 = Table("Figure 2-{:d}".format(iCentrality+nCentrality*(iTrackPt+1)))
        table2.description = "The distribution of charged particle yields within $|\Delta\\varphi| < 1.0$ correlated with the subleading jets as a function of $\Delta\eta$ in {:s} collisions. The results are shown in different dijet momentum balance bins{:s}.".format(centralityDescription[iCentrality],trackPtDescription[iTrackPt+1])
        table2.location = "Data from figure 2, located on page 10."
        table2.keywords["observables"] = ["Particle yield correlated to subleading jets"]
        #table2.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

        # Extract x- and y-axis information from the histograms
        xDeltaEta, yDeltaEta = findVariables(deltaEtaSubleading[nAsymmetry*iCentrality+(iTrackPt+1)*nAsymmetry*nCentrality : nAsymmetry+nAsymmetry*iCentrality+(iTrackPt+1)*nAsymmetry*nCentrality], deltaEtaErrorSubleading[nAsymmetry*iCentrality+(iTrackPt+1)*nAsymmetry*nCentrality : nAsymmetry+nAsymmetry*iCentrality+(iTrackPt+1)*nAsymmetry*nCentrality], "$\\frac{1}{N_{\mathrm{dijet}}}\\frac{\mathrm{d}N}{\mathrm{d}\Delta\eta}$", True, (1,iCentrality), True, iTrackPt)

        # Add the variables to the table
        table2.add_variable(xDeltaEta)
        for variable in yDeltaEta:
            table2.add_variable(variable)

        # Add the table to submission object
        submission.add_table(table2)

##########################################################################################
#            Table for leading jet radial momuntum profile. Paper figure 3a.             #
##########################################################################################

# Make separate tables from each track pT bin
for iTrackPt in range(-1,nTrackPt):
    table3a = Table("Figure 3a-{:d}".format(iTrackPt+1))
    table3a.description = "The leading jet radial momentum profiles in pp and PbPb collisions and a function of $\Delta r${:s}. The PbPb results are shown for different centrality regions.".format(trackPtDescription[iTrackPt+1])
    table3a.location = "Data from the top row of figure 3, located on page 11."
    table3a.keywords["observables"] = ["Leading jet radial momentum profile"]
    #table3a.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

    # Extract x- and y-axis information from the histograms
    xDeltaR, yDeltaR = findVariables(jetMomentumLeading[(iTrackPt+1)*nCentrality : nCentrality + (iTrackPt+1)*nCentrality], jetMomentumErrorLeading[(iTrackPt+1)*nCentrality : nCentrality + (iTrackPt+1)*nCentrality], "$P(\Delta r)$", False, (5,1), False, iTrackPt)

    # Add the variables to the table
    table3a.add_variable(xDeltaR)
    for variable in yDeltaR:
        table3a.add_variable(variable)

    # Add the table to submission object
    submission.add_table(table3a)

##########################################################################################
#          Table for subleading jet radial momuntum profile. Paper figure 3b.            #
##########################################################################################

for iTrackPt in range(-1,nTrackPt):
    table3b = Table("Figure 3b-{:d}".format(iTrackPt+1))
    table3b.description = "The subleading jet radial momentum profiles in pp and PbPb collisions as a function of $\Delta r${:s}. The PbPb results are shown for different centrality regions.".format(trackPtDescription[iTrackPt+1])
    table3b.location = "Data from the bottom row of figure 3, located on page 11."
    table3b.keywords["observables"] = ["Subleading jet radial momentum profile"]
    #table3b.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

    # Extract x- and y-axis information from the histograms
    xDeltaR, yDeltaR = findVariables(jetMomentumSubleading[(iTrackPt+1)*nCentrality : nCentrality + (iTrackPt+1)*nCentrality], jetMomentumErrorSubleading[(iTrackPt+1)*nCentrality : nCentrality + (iTrackPt+1)*nCentrality], "$P(\Delta r)$", False, (5,1), False, iTrackPt)

    # Add the variables to the table
    table3b.add_variable(xDeltaR)
    for variable in yDeltaR:
        table3b.add_variable(variable)

    # Add the table to submission object
    submission.add_table(table3b)

##########################################################################################
#         Table for leading jet radial momuntum profile ratio. Paper figure 4a.          #
##########################################################################################

table4a = Table("Figure 4a")
table4a.description = "The ratio between leading jet radial momentum profiles in PbPb and pp collisions as a function of $\Delta r$."
table4a.location = "Data from the top row of figure 4, located on page 11."
table4a.keywords["observables"] = ["Leading jet radial momentum profile ratio"]
#table4a.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(jetMomentumRatioLeading, jetMomentumRatioErrorLeading, "$P(\Delta r)_{\mathrm{PbPb}} / P(\Delta r)_{\mathrm{pp}}$", False, (4,1), False)

# Add the variables to the table
table4a.add_variable(xDeltaR)
for variable in yDeltaR:
    table4a.add_variable(variable)

# Add the table to submission object
submission.add_table(table4a)

##########################################################################################
#        Table for subleading jet radial momuntum profile ratio. Paper figure 4b.        #
##########################################################################################

table4b = Table("Figure 4b")
table4b.description = "The ratio between subleading jet radial momentum profiles in PbPb and pp collisions as a function of $\Delta r$."
table4b.location = "Data from the bottom row of figure 4, located on page 11."
table4b.keywords["observables"] = ["Subleading jet radial momentum profile ratio"]
#table4b.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(jetMomentumRatioSubleading, jetMomentumRatioErrorSubleading, "$P(\Delta r)_{\mathrm{PbPb}} / P(\Delta r)_{\mathrm{pp}}$", False, (4,1), False)

# Add the variables to the table
table4b.add_variable(xDeltaR)
for variable in yDeltaR:
    table4b.add_variable(variable)

# Add the table to submission object
submission.add_table(table4b)

##################################################################################
#              Table for leading jet shapes. Paper figure 5.                     #
##################################################################################

# Make separate tables from each centrality and track pT bin
for iTrackPt in range(-1,nTrackPtDeltaEta):
    for iCentrality in range(0,nCentrality):
        table5 = Table("Figure 5-{:d}".format(iCentrality+nCentrality*(iTrackPt+1)))
        table5.description = "Jet shapes for leading jets in {:s} collisions. The results are shown in different dijet momentum balance bins{:s}.".format(centralityDescription[iCentrality],trackPtDescription[iTrackPt+1])
        table5.location = "Data from figure 5, located on page 13."
        table5.keywords["observables"] = ["Leading jet shape"]
        #table5.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

        # Extract x- and y-axis information from the histograms
        xDeltaR, yDeltaR = findVariables(jetShapeLeading[nAsymmetry*iCentrality+(iTrackPt+1)*nAsymmetry*nCentrality : nAsymmetry+nAsymmetry*iCentrality+(iTrackPt+1)*nAsymmetry*nCentrality], jetShapeErrorLeading[nAsymmetry*iCentrality+(iTrackPt+1)*nAsymmetry*nCentrality : nAsymmetry+nAsymmetry*iCentrality+(iTrackPt+1)*nAsymmetry*nCentrality], "$\\rho(\Delta r)$", True, (1,iCentrality), False, iTrackPt)

        # Add the variables to the table
        table5.add_variable(xDeltaR)
        for variable in yDeltaR:
            table5.add_variable(variable)

        # Add the table to submission object
        submission.add_table(table5)

##################################################################################
#           Table for leading jet shape ratios. Paper figure 6.                  #
##################################################################################

# Make separate tables from each centrality bin
for iCentrality in range(0,nCentrality-1):
    table6 = Table("Figure 6-{:d}".format(iCentrality))
    table6.description = "Ratios of leading jet shapes between PbPb and pp collisions. The results from {:s} % centrality bin in PbPb are compared to pp using several dijet momentum balance selections.".format(centralityString[iCentrality])
    table6.location = "Data from figure 6, located on page 14."
    table6.keywords["observables"] = ["Leading jet shape ratio"]
    #table6.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

    # Extract x- and y-axis information from the histograms
    xDeltaR, yDeltaR = findVariables(jetShapeRatioLeading[iCentrality*nAsymmetry : nAsymmetry + iCentrality*nAsymmetry], jetShapeRatioErrorLeading[iCentrality*nAsymmetry : nAsymmetry + iCentrality*nAsymmetry], "$\\rho(\Delta r)_{\mathrm{PbPb}} / \\rho(\Delta r)_{\mathrm{pp}}$", True, (1,iCentrality), False)

    # Add the variables to the table
    table6.add_variable(xDeltaR)
    for variable in yDeltaR:
        table6.add_variable(variable)

    # Add the table to submission object
    submission.add_table(table6)

##################################################################################
#             Table for subleading jet shapes. Paper figure 7.                   #
##################################################################################

# Make separate tables from each centrality and track pT bin
for iTrackPt in range(-1,nTrackPtDeltaEta):
    for iCentrality in range(0,nCentrality):
        table7 = Table("Figure 7-{:d}".format(iCentrality+nCentrality*(iTrackPt+1)))
        table7.description = "Jet shapes for subleading jets in {:s} collisions. The results are shown in different dijet momentum balance bins{:s}.".format(centralityDescription[iCentrality],trackPtDescription[iTrackPt+1])
        table7.location = "Data from figure 7, located on page 15."
        table7.keywords["observables"] = ["Subleading jet shape"]
        #table7.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

        # Extract x- and y-axis information from the histograms
        xDeltaR, yDeltaR = findVariables(jetShapeSubleading[nAsymmetry*iCentrality+(iTrackPt+1)*nAsymmetry*nCentrality : nAsymmetry+nAsymmetry*iCentrality+(iTrackPt+1)*nAsymmetry*nCentrality], jetShapeErrorSubleading[nAsymmetry*iCentrality+(iTrackPt+1)*nAsymmetry*nCentrality : nAsymmetry+nAsymmetry*iCentrality+(iTrackPt+1)*nAsymmetry*nCentrality], "$\\rho(\Delta r)$", True, (1,iCentrality), False, iTrackPt)

        # Add the variables to the table
        table7.add_variable(xDeltaR)
        for variable in yDeltaR:
            table7.add_variable(variable)

        # Add the table to submission object
        submission.add_table(table7)

##################################################################################
#           Table for leading jet shape ratios. Paper figure 8.                  #
##################################################################################

# Make separate tables from each centrality bin
for iCentrality in range(0,nCentrality-1):
    table8 = Table("Figure 8-{:d}".format(iCentrality))
    table8.description = "Ratios of subleading jet shapes between PbPb and pp collisions. The results from {:s} % centrality bin in PbPb are compared to pp using several dijet momentum balance selections.".format(centralityString[iCentrality])
    table8.location = "Data from figure 8, located on page 16."
    table8.keywords["observables"] = ["Subleading jet shape ratio"]
    #table8.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

    # Extract x- and y-axis information from the histograms
    xDeltaR, yDeltaR = findVariables(jetShapeRatioSubleading[iCentrality*nAsymmetry : nAsymmetry + iCentrality*nAsymmetry], jetShapeRatioErrorSubleading[iCentrality*nAsymmetry : nAsymmetry + iCentrality*nAsymmetry], "$\\rho(\Delta r)_{\mathrm{PbPb}} / \\rho(\Delta r)_{\mathrm{pp}}$", True, (1,iCentrality), False)

    # Add the variables to the table
    table8.add_variable(xDeltaR)
    for variable in yDeltaR:
        table8.add_variable(variable)

    # Add the table to submission object
    submission.add_table(table8)

##################################################################################
#        Table for unbalanced leading jet shape ratio. Paper figure 9a.          #
##################################################################################

table9a = Table("Figure 9a")
table9a.description = "Ratio between unbalanced selection of leading jet shapes to all leading jet shapes in pp and PbPb collisions. The PbPb results are shown for different centrality regions."
table9a.location = "Data from the top row of figure 9, located on page 17."
table9a.keywords["observables"] = ["Leading jet shape unbalanced ratio"]
#table9.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Find the unbalanced histograms from the list of all asymmetry histograms
unbalancedHistogram = []
unbalancedError = []
for i in range(5):
    unbalancedHistogram.append(asymmetryRatioLeading[2*i])
    unbalancedError.append(asymmetryRatioErrorLeading[2*i])
    
# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(unbalancedHistogram, unbalancedError, "$\\rho(\Delta r)_{x_{j} < 0.6} / \\rho(\Delta r)_{\mathrm{all}}$", False, (5,1), False)

# Table 9a: Add the variables to the table
table9a.add_variable(xDeltaR)
for variable in yDeltaR:
    table9a.add_variable(variable)

# Add the table to submission object
submission.add_table(table9a)

################################################################################
#        Table for balanced leading jet shape ratio. Paper figure 9b.          #
################################################################################

table9b = Table("Figure 9b")
table9b.description = "Ratio between balanced selection of leading jet shapes to all leading jet shapes in pp and PbPb collisions. The PbPb results are shown for different centrality regions."
table9b.location = "Data from the bottom row of figure 9, located on page 17."
table9b.keywords["observables"] = ["Leading jet shape balanced ratio"]
#table9b.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Find the balanced histograms from the list of all asymmetry histograms
balancedHistogram = []
balancedError = []
for i in range(5):
    balancedHistogram.append(asymmetryRatioLeading[2*i+1])
    balancedError.append(asymmetryRatioErrorLeading[2*i+1])
    
# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(balancedHistogram, balancedError, "$\\rho(\Delta r)_{x_{j} > 0.8} / \\rho(\Delta r)_{\mathrm{all}}$", False, (5,1), False)

# Table 9a: Add the variables to the table
table9b.add_variable(xDeltaR)
for variable in yDeltaR:
    table9b.add_variable(variable)

# Add the table to submission object
submission.add_table(table9b)

######################################################################################
#        Table for unbalanced subleading jet shape ratio. Paper figure 10a.          #
######################################################################################

table10a = Table("Figure 10a")
table10a.description = "Ratio between unbalanced selection of subleading jet shapes to all subleading jet shapes in pp and PbPb collisions. The PbPb results are shown for different centrality regions."
table10a.location = "Data from the top row of figure 10, located on page 22."
table10a.keywords["observables"] = ["Subleading jet shape unbalanced ratio"]
#table10a.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Find the unbalanced histograms from the list of all asymmetry histograms
unbalancedHistogram = []
unbalancedError = []
for i in range(5):
    unbalancedHistogram.append(asymmetryRatioSubleading[2*i])
    unbalancedError.append(asymmetryRatioErrorSubleading[2*i])
    
# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(unbalancedHistogram, unbalancedError, "$\\rho(\Delta r)_{x_{j} < 0.6} / \\rho(\Delta r)_{\mathrm{all}}$", False, (5,1), False)

# Table 10a: Add the variables to the table
table10a.add_variable(xDeltaR)
for variable in yDeltaR:
    table10a.add_variable(variable)

# Add the table to submission object
submission.add_table(table10a)

####################################################################################
#        Table for balanced subleading jet shape ratio. Paper figure 10b.          #
####################################################################################

table10b = Table("Figure 10b")
table10b.description = "Ratio between balanced selection of subleading jet shapes to all subleading jet shapes in pp and PbPb collisions. The PbPb results are shown for different centrality regions."
table10b.location = "Data from the bottom row of figure 10, located on page 22."
table10b.keywords["observables"] = ["Subleading jet shape balanced ratio"]
#table10b.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Find the balanced histograms from the list of all asymmetry histograms
balancedHistogram = []
balancedError = []
for i in range(5):
    balancedHistogram.append(asymmetryRatioSubleading[2*i+1])
    balancedError.append(asymmetryRatioErrorSubleading[2*i+1])
    
# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(balancedHistogram, balancedError, "$\\rho(\Delta r)_{x_{j} > 0.8} / \\rho(\Delta r)_{\mathrm{all}}$", False, (5,1), False)

# Table 10b: Add the variables to the table
table10b.add_variable(xDeltaR)
for variable in yDeltaR:
    table10b.add_variable(variable)

# Add the table to submission object
submission.add_table(table10b)

######################################################################################
#                Table for normalized xj matrixes. Paper figure 11.                  #
######################################################################################

table11 = Table("Figure 11")
table11.description = "Generator-level vs. reconstructed $x_{j}$ values in the analysis $x_{j}$ bins. The plots show the probability to find a generator level $x_{j}$ for a given reconstructed $x_{j}$."
table11.location = "Data from the figure 11, located on appendix A."
table11.keywords["observables"] = ["$x_{j}$ matrix"]
#table11.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Extract axis information from the histograms
xAxis, yAxis, zAxis = findVariablesXjMatrix(xjMatrix)

# Table 11: Add the variables to the table
table11.add_variable(xAxis)
table11.add_variable(yAxis)
for variable in zAxis:
    table11.add_variable(variable)

# Add the table to submission object
submission.add_table(table11)

######################################################################################
#            Table for reverse normalized xj matrixes. Paper figure 12.              #
######################################################################################

table12 = Table("Figure 12")
table12.description = "Generator-level vs. reconstructed $x_{j}$ values in the analysis $x_{j}$ bins. The plots show the probability to find a reconstructed $x_{j}$ for a given generator level $x_{j}$."
table12.location = "Data from the figure 12, located on appendix A."
table12.keywords["observables"] = ["$x_{j}$ matrix"]
#table12.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Extract axis information from the histograms
xAxis, yAxis, zAxis = findVariablesXjMatrix(xjMatrixReverse)

# Table 12: Add the variables to the table
table12.add_variable(xAxis)
table12.add_variable(yAxis)
for variable in zAxis:
    table12.add_variable(variable)

# Add the table to submission object
submission.add_table(table12)

###########################################################################
#                      Finalize the submission                            #
###########################################################################


# Add common keywords for all tables
for table in submission.tables:
    table.keywords["cmenergies"] = [5020]   # Center-of-mass energy

# Create the submission file for upload
outputDirectory = "dijetHepData"
submission.create_files(outputDirectory)
