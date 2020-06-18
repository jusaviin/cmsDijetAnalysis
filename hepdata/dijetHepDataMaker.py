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
submission.add_record_id(01010101, "inspire")

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

# File for particle yields as a function of deltaEta For paper figures 1 and 2.
deltaEtaReader = RootFileReader("/Users/jviinika/cms/development/dijet5TeV/hepdata/hepdata_deltaEta_hin-19-013.root")

# Read the deltaEta yield histograms from the root file
deltaEtaLeading = []
deltaEtaErrorLeading = []
deltaEtaSubleading = []
deltaEtaErrorSubleading = []

for iCentrality in centralityString:
    for iXj in xjString:
        deltaEtaLeading.append(deltaEtaReader.read_hist_1d("deltaEta_trackLeadingJet_" + iCentrality + "_" + iXj))
        deltaEtaErrorLeading.append(deltaEtaReader.read_hist_1d("deltaEtaError_trackLeadingJet_" + iCentrality + "_" + iXj))
        deltaEtaSubleading.append(deltaEtaReader.read_hist_1d("deltaEta_trackSubleadingJet_" + iCentrality + "_" + iXj))
        deltaEtaErrorSubleading.append(deltaEtaReader.read_hist_1d("deltaEtaError_trackSubleadingJet_" + iCentrality + "_" + iXj))

# File for jet radial momentum distributions. For paper figures 3 and 4.
jetMomentumReader = RootFileReader("/Users/jviinika/cms/development/dijet5TeV/hepdata/hepdata_jetRadialMomentum_hin-19-013.root")

# Read the jet radial momentum histograms from the root file
jetMomentumLeading = []
jetMomentumErrorLeading = []
jetMomentumRatioLeading = []
jetMomentumRatioErrorLeading = []
jetMomentumSubleading = []
jetMomentumErrorSubleading = []
jetMomentumRatioSubleading = []
jetMomentumRatioErrorSubleading = []

for iCentrality in centralityString:
    jetMomentumLeading.append(jetMomentumReader.read_hist_1d("jetRadialMomentum_trackLeadingJet_" + iCentrality + "_allXj"))
    jetMomentumErrorLeading.append(jetMomentumReader.read_hist_1d("jetRadialMomentumError_trackLeadingJet_" + iCentrality + "_allXj"))
    jetMomentumSubleading.append(jetMomentumReader.read_hist_1d("jetRadialMomentum_trackSubleadingJet_" + iCentrality + "_allXj"))
    jetMomentumErrorSubleading.append(jetMomentumReader.read_hist_1d("jetRadialMomentumError_trackSubleadingJet_" + iCentrality + "_allXj"))
        
    # As the ratios are between PbPb and pp, there are one less of the ocmpared to just distributions
    if iCentrality != "pp":
        jetMomentumRatioLeading.append(jetMomentumReader.read_hist_1d("jetRadialMomentumRatio_trackLeadingJet_" + iCentrality + "_allXj"))
        jetMomentumRatioErrorLeading.append(jetMomentumReader.read_hist_1d("jetRadialMomentumRatioError_trackLeadingJet_" + iCentrality + "_allXj"))
        jetMomentumRatioSubleading.append(jetMomentumReader.read_hist_1d("jetRadialMomentumRatio_trackSubleadingJet_" + iCentrality + "_allXj"))
        jetMomentumRatioErrorSubleading.append(jetMomentumReader.read_hist_1d("jetRadialMomentumRatioError_trackSubleadingJet_" + iCentrality + "_allXj"))

# File for jet shapes. For paper figures 5, 6, 7 and 8.
jetShapeReader = RootFileReader("/Users/jviinika/cms/development/dijet5TeV/hepdata/hepdata_jetShapes_hin-19-013.root")

# Read the jet shape histograms from the root file
jetShapeLeading = []
jetShapeErrorLeading = []
jetShapeRatioLeading = []
jetShapeRatioErrorLeading = []
jetShapeSubleading = []
jetShapeErrorSubleading = []
jetShapeRatioSubleading = []
jetShapeRatioErrorSubleading = []

for iCentrality in centralityString:
    for iXj in xjString:
        jetShapeLeading.append(jetShapeReader.read_hist_1d("jetShape_trackLeadingJet_" + iCentrality + "_" + iXj))
        jetShapeErrorLeading.append(jetShapeReader.read_hist_1d("jetShapeError_trackLeadingJet_" + iCentrality + "_" + iXj))
        jetShapeSubleading.append(jetShapeReader.read_hist_1d("jetShape_trackSubleadingJet_" + iCentrality + "_" + iXj))
        jetShapeErrorSubleading.append(jetShapeReader.read_hist_1d("jetShapeError_trackSubleadingJet_" + iCentrality + "_" + iXj))
        
        # As the ratios are between PbPb and pp, there are one less of the ocmpared to just distributions
        if iCentrality != "pp":
            jetShapeRatioLeading.append(jetShapeReader.read_hist_1d("jetShapeRatio_trackLeadingJet_" + iCentrality + "_" + iXj))
            jetShapeRatioErrorLeading.append(jetShapeReader.read_hist_1d("jetShapeRatioError_trackLeadingJet_" + iCentrality + "_" + iXj))
            jetShapeRatioSubleading.append(jetShapeReader.read_hist_1d("jetShapeRatio_trackSubleadingJet_" + iCentrality + "_" + iXj))
            jetShapeRatioErrorSubleading.append(jetShapeReader.read_hist_1d("jetShapeRatioError_trackSubleadingJet_" + iCentrality + "_" + iXj))

# File for asymmetry ratios. For paper figures 9 and 10.
asymmetryRatioReader = RootFileReader("/Users/jviinika/cms/development/dijet5TeV/hepdata/hepdata_asymmetryRatio_hin-19-013.root")

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



# Read the variables from the histograms
from hepdata_lib import Variable, Uncertainty

# Function for finding x and y-axis variables from input histograms when the histograms are binned in deltaR
#
# Arguments:
#  valueHistogram: Histogram containing the actual values and statistical errors
#  errorHistogram: Histogram containing the systematic errors
#  yAxisName: Label given to the histogram y-axis
#  includeAsymmetry: True = Make data for all asymmetry bins. False = Only use one asymmetry bin
#  skipPP: True = pp histograms included. False = pp histograms not included (for example for ratio plots)
#  deltaEtaMode: True = Assume deltaEta histogram. False = Assume deltaR histogram
#
#  return: Returns HepData variables for x and y-axes from the input histograms
#
def findVariables(valueHistogram, errorHistogram, yAxisName, includeAsymmetry, skipPP, deltaEtaMode):

    # Define reaction and centrality labels.
    reactionLabel = ["PB PB --> CHARGED X", "PB PB --> CHARGED X", "PB PB --> CHARGED X", "PB PB --> CHARGED X", "P P --> CHARGED X"]
    centralityLabel = ["0-10%", "10-30%", "30-50%", "50-90%", "pp"]
    asymmetryLabel = ["$0 < x_{j} < 0.6$", "$0.6 < x_{j} < 0.8$", "$0.8 < x_{j} < 1.0$", "All $x_{j}$"]
    
    nAsymmetry = 1
    if includeAsymmetry:
        nAsymmetry = 4
        
    nCentrality = 5
    if skipPP:
        nCentrality = 4
    
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
    for i in range(nCentrality):
        for iXj in range(nAsymmetry):
            myVariable = Variable(yAxisName, is_independent=False, is_binned=False, units="")
            myVariable.values = valueHistogram[i*nAsymmetry+iXj]["y"]
            myVariable.add_qualifier("reaction",reactionLabel[i])
            myVariable.add_qualifier("centrality",centralityLabel[i])
            if includeAsymmetry:
                myVariable.add_qualifier("$x_{j}$",asymmetryLabel[iXj])
        
            # If there are points outside of the analysis range, remove them from the values table:
            for j in range(overRange):
                myVariable.values.pop()
                
            for j in range(underRange):
                myVariable.values.pop(0)
        
            yAxis.append(myVariable)
        
    # Read the uncertainties
    for i in range(nCentrality):
        for iXj in range(nAsymmetry):
            statUncertainty = Uncertainty("stat", is_symmetric=True)
            statUncertainty.values = valueHistogram[i*nAsymmetry+iXj]["dy"]

            systUncertainty = Uncertainty("sys", is_symmetric=True)
            systUncertainty.values = errorHistogram[i*nAsymmetry+iXj]["dy"]
        
            # If there are points outside of the analysis range, remove their errors from the table:
            for j in range(overRange):
                statUncertainty.values.pop()
                systUncertainty.values.pop()
                
            for j in range(underRange):
                statUncertainty.values.pop(0)
                systUncertainty.values.pop(0)

            yAxis[i*nAsymmetry+iXj].add_uncertainty(statUncertainty)
            yAxis[i*nAsymmetry+iXj].add_uncertainty(systUncertainty)
        
    return xAxis, yAxis

######################################################################################
#               Table for leading jet deltaEta yield. Paper figure 1.                #
######################################################################################

table1 = Table("Figure 1")
table1.description = "This table tell whatever is shown in the figure 1 of our very best paper."
table1.location = "Data from Figure 1, located on page NN."
table1.keywords["observables"] = ["Particle yield correlated to leading jets"]
#table1.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Extract x- and y-axis information from the histograms
xDeltaEta, yDeltaEta = findVariables(deltaEtaLeading, deltaEtaErrorLeading, "$Y = \\frac{1}{N_{\mathrm{dijet}}}\\frac{\mathrm{d}N}{\mathrm{d}\Delta\eta}$", True, False, True)

# Add the variables to the table
table1.add_variable(xDeltaEta)
for variable in yDeltaEta:
    table1.add_variable(variable)

# Add the table to submission object
submission.add_table(table1)

######################################################################################
#             Table for subleading jet deltaEta yield. Paper figure 2.               #
######################################################################################

table2 = Table("Figure 2")
table2.description = "This table tell whatever is shown in the figure 2 of our very best paper."
table2.location = "Data from Figure 2, located on page NN."
table2.keywords["observables"] = ["Particle yield correlated to subleading jets"]
#table1.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Extract x- and y-axis information from the histograms
xDeltaEta, yDeltaEta = findVariables(deltaEtaSubleading, deltaEtaErrorSubleading, "$Y = \\frac{1}{N_{\mathrm{dijet}}}\\frac{\mathrm{d}N}{\mathrm{d}\Delta\eta}$", True, False, True)

# Add the variables to the table
table2.add_variable(xDeltaEta)
for variable in yDeltaEta:
    table2.add_variable(variable)

# Add the table to submission object
submission.add_table(table2)

##########################################################################################
#            Table for leading jet radial momuntum profile. Paper figure 3a.             #
##########################################################################################

table3a = Table("Figure 3a")
table3a.description = "This table tell whatever is shown in the figure 3a of our very best paper."
table3a.location = "Data from Figure 3a, located on page NN."
table3a.keywords["observables"] = ["Leading jet radial momentum profile"]
#table3a.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(jetMomentumLeading, jetMomentumErrorLeading, "$P(\Delta r)$", False, False, False)

# Add the variables to the table
table3a.add_variable(xDeltaR)
for variable in yDeltaR:
    table3a.add_variable(variable)

# Add the table to submission object
submission.add_table(table3a)

##########################################################################################
#          Table for subleading jet radial momuntum profile. Paper figure 3b.            #
##########################################################################################

table3b = Table("Figure 3b")
table3b.description = "This table tell whatever is shown in the figure 3b of our very best paper."
table3b.location = "Data from Figure 3b, located on page NN."
table3b.keywords["observables"] = ["Subleading jet radial momentum profile"]
#table3b.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(jetMomentumSubleading, jetMomentumErrorSubleading, "$P(\Delta r)$", False, False, False)

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
table4a.description = "This table tell whatever is shown in the figure 4a of our very best paper."
table4a.location = "Data from Figure 4a, located on page NN."
table4a.keywords["observables"] = ["Leading jet radial momentum profile ratio"]
#table4a.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(jetMomentumRatioLeading, jetMomentumRatioErrorLeading, "$P(\Delta r)_{\mathrm{PbPb}} / P(\Delta r)_{\mathrm{pp}}$", False, True, False)

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
table4b.description = "This table tell whatever is shown in the figure 4b of our very best paper."
table4b.location = "Data from Figure 4b, located on page NN."
table4b.keywords["observables"] = ["Subleading jet radial momentum profile ratio"]
#table4b.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(jetMomentumRatioSubleading, jetMomentumRatioErrorSubleading, "$P(\Delta r)_{\mathrm{PbPb}} / P(\Delta r)_{\mathrm{pp}}$", False, True, False)

# Add the variables to the table
table4b.add_variable(xDeltaR)
for variable in yDeltaR:
    table4b.add_variable(variable)

# Add the table to submission object
submission.add_table(table4b)

##################################################################################
#              Table for leading jet shapes. Paper figure 5.                     #
##################################################################################

table5 = Table("Figure 5")
table5.description = "This table tell whatever is shown in the figure 5 of our very best paper."
table5.location = "Data from Figure 5, located on page NN."
table5.keywords["observables"] = ["Leading jet shape"]
#table5.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(jetShapeLeading, jetShapeErrorLeading, "$\\rho(\Delta r)$", True, False, False)

# Add the variables to the table
table5.add_variable(xDeltaR)
for variable in yDeltaR:
    table5.add_variable(variable)

# Add the table to submission object
submission.add_table(table5)

##################################################################################
#           Table for leading jet shape ratios. Paper figure 6.                  #
##################################################################################

table6 = Table("Figure 6")
table6.description = "This table tell whatever is shown in the figure 6 of our very best paper."
table6.location = "Data from Figure 6, located on page NN."
table6.keywords["observables"] = ["Leading jet shape ratio"]
#table6.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(jetShapeRatioLeading, jetShapeRatioErrorLeading, "$\\rho(\Delta r)_{\mathrm{PbPb}}$ / \\rho(\Delta r)_{\mathrm{pp}}$", True, True, False)

# Add the variables to the table
table6.add_variable(xDeltaR)
for variable in yDeltaR:
    table6.add_variable(variable)

# Add the table to submission object
submission.add_table(table6)

##################################################################################
#             Table for subleading jet shapes. Paper figure 7.                   #
##################################################################################

table7 = Table("Figure 7")
table7.description = "This table tell whatever is shown in the figure 7 of our very best paper."
table7.location = "Data from Figure 7, located on page NN."
table7.keywords["observables"] = ["Subleading jet shape"]
#table7.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(jetShapeSubleading, jetShapeErrorSubleading, "$\\rho(\Delta r)$", True, False, False)

# Add the variables to the table
table7.add_variable(xDeltaR)
for variable in yDeltaR:
    table7.add_variable(variable)

# Add the table to submission object
submission.add_table(table7)

##################################################################################
#           Table for leading jet shape ratios. Paper figure 8.                  #
##################################################################################

table8 = Table("Figure 8")
table8.description = "This table tell whatever is shown in the figure 8 of our very best paper."
table8.location = "Data from Figure 8, located on page NN."
table8.keywords["observables"] = ["Subleading jet shape ratio"]
#table8.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(jetShapeRatioSubleading, jetShapeRatioErrorSubleading, "$\\rho(\Delta r)_{\mathrm{PbPb}}$ / \\rho(\Delta r)_{\mathrm{pp}}$", True, True, False)

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
table9a.description = "This table tell whatever is shown in the figure 9 of our very best paper."
table9a.location = "Data from Figure 9a, located on page NN."
table9a.keywords["observables"] = ["Leading jet shape unbalanced ratio"]
#table9.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Find the unbalanced histograms from the list of all asymmetry histograms
unbalancedHistogram = []
unbalancedError = []
for i in range(5):
    unbalancedHistogram.append(asymmetryRatioLeading[2*i])
    unbalancedError.append(asymmetryRatioErrorLeading[2*i])
    
# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(unbalancedHistogram, unbalancedError, "$\\rho(\Delta r)_{x_{j} < 0.6} / \\rho(\Delta r)_{\mathrm{all}}$", False, False, False)

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
table9b.description = "This table tell whatever is shown in the figure 9 of our very best paper."
table9b.location = "Data from Figure 9b, located on page NN."
table9b.keywords["observables"] = ["Leading jet shape balanced ratio"]
#table9b.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Find the balanced histograms from the list of all asymmetry histograms
balancedHistogram = []
balancedError = []
for i in range(5):
    balancedHistogram.append(asymmetryRatioLeading[2*i+1])
    balancedError.append(asymmetryRatioErrorLeading[2*i+1])
    
# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(balancedHistogram, balancedError, "$\\rho(\Delta r)_{x_{j} > 0.8} / \\rho(\Delta r)_{\mathrm{all}}$", False, False, False)

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
table10a.description = "This table tell whatever is shown in the figure 10 of our very best paper."
table10a.location = "Data from Figure 10a, located on page NN."
table10a.keywords["observables"] = ["Subleading jet shape unbalanced ratio"]
#table10a.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Find the unbalanced histograms from the list of all asymmetry histograms
unbalancedHistogram = []
unbalancedError = []
for i in range(5):
    unbalancedHistogram.append(asymmetryRatioSubleading[2*i])
    unbalancedError.append(asymmetryRatioErrorSubleading[2*i])
    
# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(unbalancedHistogram, unbalancedError, "$\\rho(\Delta r)_{x_{j} < 0.6} / \\rho(\Delta r)_{\mathrm{all}}$", False, False, False)

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
table10b.description = "This table tell whatever is shown in the figure 10 of our very best paper."
table10b.location = "Data from Figure 10b, located on page NN."
table10b.keywords["observables"] = ["Subleading jet shape balanced ratio"]
#table10b.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Find the balanced histograms from the list of all asymmetry histograms
balancedHistogram = []
balancedError = []
for i in range(5):
    balancedHistogram.append(asymmetryRatioSubleading[2*i+1])
    balancedError.append(asymmetryRatioErrorSubleading[2*i+1])
    
# Extract x- and y-axis information from the histograms
xDeltaR, yDeltaR = findVariables(balancedHistogram, balancedError, "$\\rho(\Delta r)_{x_{j} > 0.8} / \\rho(\Delta r)_{\mathrm{all}}$", False, False, False)

# Table 9a: Add the variables to the table
table10b.add_variable(xDeltaR)
for variable in yDeltaR:
    table10b.add_variable(variable)

# Add the table to submission object
submission.add_table(table10b)

###########################################################################
#                      Finalize the submission                            #
###########################################################################


# Add common keywords for all tables
for table in submission.tables:
    table.keywords["cmenergies"] = [5020]   # Center-of-mass energy

# Create the submission file for upload
outputDirectory = "dijetHepData"
submission.create_files(outputDirectory)
