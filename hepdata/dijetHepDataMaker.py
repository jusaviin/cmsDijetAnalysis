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
# Add the first table
from hepdata_lib import Table
table1 = Table("Figure 1")
table1.description = "This table tell whatever is shown in the figure 1 of our very best paper."
table1.location = "Data from Figure 1, located on page NN."
table1.keywords["observables"] = ["Yield"]
#table1.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

# Read the root file with the histograms
from hepdata_lib import RootFileReader

reader = RootFileReader("/Users/jviinika/cms/development/dijet5TeV/data/publishedResults/officialHist_py_deta_16_020.root")

# Read the histograms for deltaEta yield and uncertainty from the root file
deltaEtaHistogram = []
deltaEtaError = []
for i in range(5):
    deltaEtaHistogram.append(reader.read_hist_1d("py_deta_all_"+str(i)))
    deltaEtaError.append(reader.read_hist_1d("py_deta_err_all_"+str(i)))

# Read the variables from the histograms
from hepdata_lib import Variable, Uncertainty

# Table1: x-axis: deltaEta value
xDeltaEta = Variable("$\Delta\eta$", is_independent=True, is_binned=True, units="")
xDeltaEta.values = deltaEtaHistogram[0]["x_edges"]

# Table1: y-axis Yield in each deltaEta bin
yDeltaEta = []
reactionLabel = ["PB PB --> CHARGED X", "PB PB --> CHARGED X", "PB PB --> CHARGED X", "PB PB --> CHARGED X", "P P --> CHARGED X"]
centralityLabel = ["0-10%", "10-30%", "30-50%", "50-90%", "pp"]

for i in range(5):
   myVariable = Variable("$Y=\\frac{1}{N_{\\text{jet}}}\\frac{dN}{d\Delta\eta}$", is_independent=False, is_binned=False, units="")
   myVariable.values = deltaEtaHistogram[i]["y"]
   myVariable.add_qualifier("reaction",reactionLabel[i])
   myVariable.add_qualifier("centrality",centralityLabel[i])
   yDeltaEta.append(myVariable)
   
# Table1: uncertainties for histograms
for i in range(5):
    statUncertainty = Uncertainty("Statistical uncertainty", is_symmetric=True)
    statUncertainty.values = deltaEtaHistogram[i]["dy"]

    systUncertainty = Uncertainty("Systematic uncertainty", is_symmetric=True)
    systUncertainty.values = deltaEtaError[i]["dy"]

    yDeltaEta[i].add_uncertainty(statUncertainty)
    yDeltaEta[i].add_uncertainty(systUncertainty)

# Table 1: Add the variables to the table
table1.add_variable(xDeltaEta)
for variable in yDeltaEta:
    table1.add_variable(variable)

# Add the table to submission object
submission.add_table(table1)

# Add common keywords for all tables
for table in submission.tables:
    table.keywords["cmenergies"] = [5020]   # Center-of-mass energy

# Create the submission file for upload
outputDirectory = "dijetHepData"
submission.create_files(outputDirectory)
