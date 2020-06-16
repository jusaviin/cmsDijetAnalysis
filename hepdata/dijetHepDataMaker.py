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

# Test file with old deltaEta histograms
reader = RootFileReader("/Users/jviinika/cms/development/dijet5TeV/data/publishedResults/officialHist_py_deta_16_020.root")

# Read the histograms for deltaEta yield and uncertainty from the root file
deltaEtaHistogram = []
deltaEtaError = []
for i in range(5):
    deltaEtaHistogram.append(reader.read_hist_1d("py_deta_all_"+str(i)))
    deltaEtaError.append(reader.read_hist_1d("py_deta_err_all_"+str(i)))

# File for asymmetry ratios. For paper figures 9 and 10.
asymmetryRatioReader = RootFileReader("/Users/jviinika/cms/development/dijet5TeV/hepdata/hepdata_asymmetryRatio_hin-19-013.root")

# Read the asymmetry ratio histograms from the root file
centralityString = ["0-10","10-30","30-50","50-90","pp"]
asymmetryString = ["Imbalanced", "Balanced"]

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
def findVariablesDeltaR(valueHistogram, errorHistogram, yAxisName):

    # Define reaction and centrality labels.
    reactionLabel = ["PB PB --> CHARGED X", "PB PB --> CHARGED X", "PB PB --> CHARGED X", "PB PB --> CHARGED X", "P P --> CHARGED X"]
    centralityLabel = ["0-10%", "10-30%", "30-50%", "50-90%", "pp"]
    
    # x-axis: deltaR value
    xDeltaR = Variable("$\Delta r$", is_independent=True, is_binned=True, units="")
    xDeltaR.values = valueHistogram[0]["x_edges"]

    # Determine how many bins in the histogram is outside of the analysis range of [0.1]
    overRange = 0
    for (i,j) in xDeltaR.values:
        if j > 1.1:
            overRange = overRange + 1
            
    # Remove the points that are outside of the range
    for i in range(overRange):
        xDeltaR.values.pop()
        
    # Define a list for y-axis values
    yDeltaR = []

    # Read the values
    for i in range(5):
        myVariable = Variable(yAxisName, is_independent=False, is_binned=False, units="")
        myVariable.values = valueHistogram[i]["y"]
        myVariable.add_qualifier("reaction",reactionLabel[i])
        myVariable.add_qualifier("centrality",centralityLabel[i])
        
        # If there are points outside of the analysis range, remove them from the values table:
        for j in range(overRange):
            myVariable.values.pop()
        
        yDeltaR.append(myVariable)
        
    # Read the uncertainties
    for i in range(5):
        statUncertainty = Uncertainty("stat", is_symmetric=True)
        statUncertainty.values = valueHistogram[i]["dy"]

        systUncertainty = Uncertainty("sys", is_symmetric=True)
        systUncertainty.values = errorHistogram[i]["dy"]
        
        # If there are points outside of the analysis range, remove their errors from the table:
        for j in range(overRange):
            statUncertainty.values.pop()
            systUncertainty.values.pop()

        yDeltaR[i].add_uncertainty(statUncertainty)
        yDeltaR[i].add_uncertainty(systUncertainty)
        
    return xDeltaR, yDeltaR

###########################################################################
# Test table: old data for testing purposes. Remove before submission!!!! #
###########################################################################

#table1 = Table("Figure 1")
#table1.description = "This table tell whatever is shown in the figure 1 of our very best paper."
#table1.location = "Data from Figure 1, located on page NN."
#table1.keywords["observables"] = ["Yield"]
##table1.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image
#
## Table1: x-axis: deltaEta value
#xDeltaEta = Variable("$\Delta\eta$", is_independent=True, is_binned=True, units="")
#xDeltaEta.values = deltaEtaHistogram[0]["x_edges"]
#
## Table1: y-axis Yield in each deltaEta bin
#yDeltaEta = []
#
#for i in range(5):
#   myVariable = Variable("$Y=\\frac{1}{N_{\\text{jet}}}\\frac{dN}{d\Delta\eta}$", is_independent=False, is_binned=False, units="")
#   myVariable.values = deltaEtaHistogram[i]["y"]
#   myVariable.add_qualifier("reaction",reactionLabel[i])
#   myVariable.add_qualifier("centrality",centralityLabel[i])
#   yDeltaEta.append(myVariable)
#
## Table1: uncertainties for histograms
#for i in range(5):
#    statUncertainty = Uncertainty("stat", is_symmetric=True)
#    statUncertainty.values = deltaEtaHistogram[i]["dy"]
#
#    systUncertainty = Uncertainty("sys", is_symmetric=True)
#    systUncertainty.values = deltaEtaError[i]["dy"]
#
#    yDeltaEta[i].add_uncertainty(statUncertainty)
#    yDeltaEta[i].add_uncertainty(systUncertainty)
#
## Table 1: Add the variables to the table
#table1.add_variable(xDeltaEta)
#for variable in yDeltaEta:
#    table1.add_variable(variable)
#
## Add the table to submission object
#submission.add_table(table1)

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
xDeltaR, yDeltaR = findVariablesDeltaR(unbalancedHistogram, unbalancedError, "$\\rho(\Delta r)_{x_{j} < 0.6} / \\rho(\Delta r)_{all}$")

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
xDeltaR, yDeltaR = findVariablesDeltaR(balancedHistogram, balancedError, "$\\rho(\Delta r)_{x_{j} > 0.8} / \\rho(\Delta r)_{all}$")

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
xDeltaR, yDeltaR = findVariablesDeltaR(unbalancedHistogram, unbalancedError, "$\\rho(\Delta r)_{x_{j} < 0.6} / \\rho(\Delta r)_{all}$")

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
xDeltaR, yDeltaR = findVariablesDeltaR(balancedHistogram, balancedError, "$\\rho(\Delta r)_{x_{j} > 0.8} / \\rho(\Delta r)_{all}$")

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
