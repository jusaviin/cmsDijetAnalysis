PROGRAM       = dijetAnalysis

version       = dijets
CXX           = g++
CXXFLAGS      = -g -Wall -Wno-deprecated -D$(version) 
LD            = g++
LDFLAGS       = -O2
SOFLAGS       = -shared
#############################################
# -bind_at_load helps to remove linker error
############################################
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS  += $(shell root-config --libs)

# Put to hdrsdict header files from all classes that inherit TObject
# HDRSDICT = 
        
# Use the following form if you have classes inherint TObject
# HDRS += $(HDRSDICT) src/Class.h ... nanoDict.h       
HDRS += src/ForestReader.h src/HighForestReader.h src/SkimForestReader.h src/GeneratorLevelForestReader.h src/GeneratorLevelSkimForestReader.h src/DijetHistograms.h src/DijetAnalyzer.h src/ConfigurationCard.h src/TrkCorr.h src/TrkSettings.h src/JffCorrection.h src/MixedEventLookoutTable.h src/XiaoTrkCorr.h src/TrkCorrInterface.h src/JetCorrector.h src/JetUncertainty.h src/trackingEfficiency2018PbPb.h src/trackingEfficiency2017pp.h src/TrackingEfficiencyInterface.h src/MixingForestReader.h src/GeneratorLevelMixingForestReader.h src/TrackPreCorrector.h

SRCS = $(HDRS:.h=.cxx)
OBJS = $(HDRS:.h=.o)

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS) $(PROGRAM).cxx
		@echo "Linking $(PROGRAM) ..."
		$(CXX) -lEG -lPhysics -L$(PWD) $(PROGRAM).cxx $(CXXFLAGS) $(OBJS) $(LDFLAGS) -o $(PROGRAM)
		@echo "done"

%.cxx:

%: %.cxx
#  commands to execute (built-in):
	$(LINK.cc) $^ $(CXXFLAGS) $(LOADLIBES) $(LDLIBS) -o $@

%.o: %.cxx %.h
#  commands to execute (built-in):
	$(COMPILE.cc) $(OUTPUT_OPTION) $<

# If dictionaries built, need to clean also them: *Dict*
clean:
		rm -rf $(OBJS) $(PROGRAM).o *.dSYM $(PROGRAM)

cl:  clean $(PROGRAM)

# Dictionary is needed for all classes inheriting TObject from root
# nanoDict.cc: $(HDRSDICT)
#		@echo "Generating dictionary ..."
#		@rm -f nanoDict.cc nanoDict.hh nanoDict.h
#		@rootcling nanoDict.cc -c -D$(version) $(HDRSDICT)
