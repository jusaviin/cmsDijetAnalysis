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
submission.read_abstract("vnAbstract.txt")
submission.add_link("Webpage with all figures and tables", "http://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/HIN-21-002/index.html")
#submission.add_link("arXiv", "http://arxiv.org/abs/arXiv:2101.04720")
#submission.add_record_id(1840683, "inspire")

# All the figures are inserted as table objects to HepData
from hepdata_lib import Table

# Read the root files with the histograms
from hepdata_lib import RootFileReader

################################################################################
#      Read the root files containing final results from dijet vnanalysis      #
################################################################################

# Common naming convention in the final results files
centralityString = ["0-10", "10-30", "30-50"]

# Numbers of bins
nCentrality = 3       # Number of centrality bins in results
nFlow = 3             # Number of vn components that are analyzed (2, 3 and 4)

# File for dijet vn results as a function of hadron pT. For paper figure 1.
vnVsPtReader = RootFileReader("hepdata_dijetVnPt_hin-21-002_test.root")

# Read the dijet vn graphs as a function of hadron pT from the root file
dijetVnHadronPt = []
dijetVnHadronPtError = []

# In the analysis we measure the vn components from 2 to 4
for centralityLabel in centralityString:
    for iFlow in range(2,5):
        dijetVnHadronPt.append(vnVsPtReader.read_graph("dijetV{0}VsPt_C{1}".format(iFlow, centralityLabel)))
        dijetVnHadronPtError.append(vnVsPtReader.read_graph("dijetV{0}VsPtError_C{1}".format(iFlow, centralityLabel)))

# File for dijet vn results as a function of centrality. For paper figure 2.
centralityReader = RootFileReader("hepdata_dijetVnCentrality_hin-21-002_test.root")

# Read the dijet vn graphs as a function of centrality from the root file
dijetVnCentrality = []
dijetVnCentralityError = []

# In the analysis we measure the vn components from 2 to 4
for iFlow in range(2,5):
    dijetVnCentrality.append(centralityReader.read_graph("dijetV{}Centrality".format(iFlow)))
    dijetVnCentralityError.append(centralityReader.read_graph("dijetV{}CentralityError".format(iFlow)))

# Reaction labels to be added as keywords
pbpbReactionLabel = "PB PB --> DIJET CHARGED X"


# Read the variables from the histograms
from hepdata_lib import Variable, Uncertainty

# Function for finding x and y-axis variables from input graphs
#
# Arguments:
#  valueGraph: Graph containing the actual values and statistical errors
#  errorGraph: Graph containing the systematic errors
#  centralityBin: -1 if binned in centrality, otherwise the centrality bin of the histogram
#
#  return: Returns HepData variables for x and y-axes from the input histograms
#
def findVariables(valueGraph, errorGraph, centralityBin = -1):

    # Define reaction and centrality labels.
    pbpbReactionLabel = "PB PB --> DIJET CHARGED X"
    pbpbEnergyLabel = "$\sqrt{s_{\mathrm{NN}}}$"
    centralityLabel = ["0-10%", "10-30%", "30-50%"]
    trackPtBins = [(0.7,1),(1,2),(2,3)]
    centralityBins = [(0,10),(10,30),(30,50)]
    
    # x-axis: Either centrality or track pT
    if centralityBin < 0:
        variableName = "Centrality"
        xAxis = Variable(variableName, is_independent=True, is_binned=True, units="%")
        xAxis.values = centralityBins
    else:
        variableName = "Hadron $p_{\mathrm{T}}$"
        xAxis = Variable(variableName, is_independent=True, is_binned=True, units="GeV")
        xAxis.values = trackPtBins
        
        # Commented out option to put the point to the mean pT within the bin
        #xAxis = Variable(variableName, is_independent=True, is_binned=False, units="GeV")
        #xAxis.values = valueGraph[0]["x"]

    # Define a list for y-axis values
    yAxis = []

    nFlow = 3 # Three flow components are analyzed: 2, 3 and 4
    # Read the values
    for iFlow in range(0,nFlow):
        myVariable = Variable("Dijet $v_{{{}}}$".format(iFlow+2), is_independent=False, is_binned=False, units="")
        myVariable.values = valueGraph[iFlow]["y"]
        myVariable.add_qualifier(pbpbEnergyLabel,"5.02 TeV")
        myVariable.add_qualifier("Reaction",pbpbReactionLabel)
        if centralityBin >= 0:
            myVariable.add_qualifier("Centrality",centralityLabel[centralityBin])
        myVariable.add_qualifier("Jet algorithm", "Anti-k$_{\mathrm{T}}$ R = 0.4")
        myVariable.add_qualifier("Leading jet $p_{\mathrm{T}}$", "> 120 GeV")
        myVariable.add_qualifier("Subleading jet $p_{\mathrm{T}}$", "> 50 GeV")
        myVariable.add_qualifier("$|\eta^{\mathrm{jet}}|$", "< 1.3")
        myVariable.add_qualifier("$\Delta\\varphi_{\mathrm{subleading}}^{\mathrm{leading}}$","$> 5\pi/6$")
        if centralityBin < 0:
            myVariable.add_qualifier("$p_{\mathrm{T}}^{\mathrm{ch}}$","$0.7 < p_{\mathrm{T}}^{\mathrm{ch}} < 3$ GeV")
        myVariable.add_qualifier("$|\eta^{\mathrm{ch}}|$", "< 2.4")
        
        yAxis.append(myVariable)
        
    # Read the uncertainties
    for iFlow in range(0,nFlow):
        statUncertainty = Uncertainty("stat", is_symmetric=True)
        statUncertainty.values = valueGraph[iFlow]["dy"]

        systUncertainty = Uncertainty("sys", is_symmetric=True)
                                
        # If we are looking vn as a function of pT, each bin has the same systematic uncertainty and there is only one value for it in the file
        if centralityBin >= 0:
            systematicUncertaintyValue = errorGraph[iFlow]["dy"][0]
            systematicUncertaintyList = []
            for iBin in range(0,len(valueGraph[iFlow]["dy"])):
                systematicUncertaintyList.append(systematicUncertaintyValue)
            systUncertainty.values = systematicUncertaintyList
        
        # For vn as a function of centrality, there is a different systematic uncertainty for each bin
        else:
            systUncertainty.values = errorGraph[iFlow]["dy"]

        yAxis[iFlow].add_uncertainty(statUncertainty)
        yAxis[iFlow].add_uncertainty(systUncertainty)
        
    return xAxis, yAxis



######################################################################################
#           Table for dijet vn as a function of hadron pT. Paper figure 1.           #
######################################################################################

# Make separate tables for each centrality bin, combine all vn:s within a centrality bin
for iCentrality in range(0,nCentrality):
    table1 = Table("Figure 1-{:d}".format(iCentrality))
    table1.description = "The dijet $v_{{n}}$ data points factorized using different associated hadron pT bins for {} % centrality bin. The data points are corrected for the jet reconstruction bias effects.".format(centralityString[iCentrality])
    table1.location = "Data for the {} % centrality bin from figure 1, located on page 7.".format(centralityString[iCentrality])
    table1.keywords["observables"] = ["Dijet $v_{2}$","Dijet $v_{3}$","Dijet $v_{4}$"]
    table1.keywords["reactions"] = [pbpbReactionLabel]
    #table1.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

    # Extract x- and y-axis information from the graphs
    xFlow, yFlow = findVariables(dijetVnHadronPt[nFlow*iCentrality : nFlow*iCentrality+nFlow], dijetVnHadronPtError[nFlow*iCentrality : nFlow*iCentrality+nFlow], iCentrality)

    # Add the variables to the table
    table1.add_variable(xFlow)
    for variable in yFlow:
        table1.add_variable(variable)

    # Add the table to submission object
    submission.add_table(table1)

######################################################################################
#          Table for dijet vn as a function of centrality. Paper figure 2.           #
######################################################################################

# Combine all vn values to one table
table2 = Table("Figure 2")
table2.description = "Final dijet $v_{2}$, $v_{3}$ and $v_{4}$ results presented as a function of centrality."
table2.location = "Data from figure 2, located on page 8."
table2.keywords["observables"] = ["Dijet $v_{2}$","Dijet $v_{3}$","Dijet $v_{4}$"]
table2.keywords["reactions"] = [pbpbReactionLabel]
#table2.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Extract x- and y-axis information from the graphs
xFlow, yFlow = findVariables(dijetVnCentrality, dijetVnCentralityError)

# Add the variables to the table
table2.add_variable(xFlow)
for variable in yFlow:
    table2.add_variable(variable)

# Add the table to submission object
submission.add_table(table2)


###########################################################################
#                      Finalize the submission                            #
###########################################################################


# Add common keywords for all tables
for table in submission.tables:
    table.keywords["cmenergies"] = [5020]   # Center-of-mass energy

# Create the submission file for upload
outputDirectory = "dijetVnHepData"
submission.create_files(outputDirectory)
