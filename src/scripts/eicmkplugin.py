#!/usr/bin/env python3
from sys import argv
import subprocess
import sys
import os

# TODO: This needs to be enhanced to all other types of plugins like ones using just JEventProcessor instead
# TODO: of JEventProcessorSequentialRoot. Also, this should be able to makea factory plugin.

cmakelists_template = """
#
# EICrecon stand alone plugin
#
# This file was generated by the eicmkplugin.py script
# 
# Before using this, make sure your environment is set up
# for the EICrecon version you are using. The easiest way
# to do this is to source the bin/eicrecon-this.(c)sh
# script.
#

cmake_minimum_required(VERSION 3.16)
project({0}_project)

find_package(EICrecon REQUIRED)

set(CMAKE_CXX_STANDARD ${{EICrecon_CXX_STANDARD}})

# Automatically determine source file list.
file(GLOB mysourcefiles *.cpp *.cc *.c  *.hpp *.hh *.h)
set( {0}_PLUGIN_SOURCES ${{mysourcefiles}} )

# Create plugin
add_library({0}_plugin SHARED ${{{0}_PLUGIN_SOURCES}})
target_link_libraries({0}_plugin ${{JANA_LIB}} ${{ROOT_LIBRARIES}})
set_target_properties({0}_plugin PROPERTIES PREFIX "" OUTPUT_NAME "{0}" SUFFIX ".so")

# Install plugin USER_PLUGIN_OUTPUT_DIRECTORY is set depending on EICrecon_MY envar.
install(TARGETS {0}_plugin DESTINATION ${{USER_PLUGIN_OUTPUT_DIRECTORY}} )

"""

processor_sequentialroot_header_template = """
//
// Template for this file generated with eicmkplugin.py
//

#include <JANA/JEventProcessorSequentialRoot.h>
#include <TH2D.h>
#include <TFile.h>

// Include appropirate class headers. e.g.
// #include <edm4hep/SimCalorimeterHit.h>
// #include <detectors/BEMC/BEMCRawCalorimeterHit.h>


class {0}: public JEventProcessorSequentialRoot {{
private:

    // Data objects we will need from JANA e.g.
    // PrefetchT<edm4hep::SimCalorimeterHit> rawhits   = {{this, "EcalBarrelHits"}};
    // PrefetchT<BEMCRawCalorimeterHit>      digihits  = {{this}};

    // Declare histogram and tree pointers here. e.g.
    // TH1D* hEraw  = nullptr;
    // TH2D* hEdigi = nullptr ;

public:
    {0}() {{ SetTypeName(NAME_OF_THIS); }}
    
    void InitWithGlobalRootLock() override;
    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock() override;
}};
"""

processor_sequentialroot_implementation_template = """
//
// Template for this file generated with eicmkplugin.py
//

#include "{0}.h"
     
// The following just makes this a JANA plugin
extern "C" {{
    void InitPlugin(JApplication *app) {{
        InitJANAPlugin(app);
        app->Add(new {0});
    }}
}}

//-------------------------------------------
// InitWithGlobalRootLock
//-------------------------------------------
void {0}::InitWithGlobalRootLock(){{

    // Create histograms here. e.g.
    // hEraw  = new TH1D("Eraw",  "BEMC hit energy (raw)",  100, 0, 0.075);
    // hEdigi = new TH2D("Edigi", "BEMC hit energy (digi) vs. raw", 200, 0, 2000.0, 100, 0, 0.075);
}}

//-------------------------------------------
// ProcessSequential
//-------------------------------------------
void {0}::ProcessSequential(const std::shared_ptr<const JEvent>& event) {{

    // Fill histograms here. e.g.
    // for( auto hit : rawhits()  ) hEraw->Fill(  hit->getEnergy());
    // for( auto hit : digihits() ) hEdigi->Fill( hit->getAmplitude(), hit->getEnergy());
}}

//-------------------------------------------
// FinishWithGlobalRootLock
//-------------------------------------------
void {0}::FinishWithGlobalRootLock() {{

    // Do any final calculations here.
}}

"""

def print_usage():
    print("""
Usage:
       eicmkplugin.py name
    
    This will create a skeleton plugin for EICrecon with the given name. A directory with the
given name will be created and some files placed there that can be used to build your custom
plugin. The expected use case is to create plugins that create histograms or trees using the
eicrecon executable.

    The CMakeLists.txt file placed in the plugin directory will install the plugin either in
the plugins directory relative to your EICrecon_ROOT environment variable when cmake is run
or relative to the directory specified by your EICrecon_MY environment variable, if it is set.
The intent of the EICrecon_MY environment variable is to allow you to build a custom plugin
against a global build of EICrecon where you do not have write access.

Example usage:

 eicmkplugin.py MyCustomPlugin
 cmake -S MyCustomPlugin -B MyCustomPlugin/build
 cmake --build MyCustomPlugin/build --target install

 eicrecon -Pplugins=MyCustomPlugin file.root
   
    """)


# Print short usage statement if no arguments
if len(argv) < 2:
    print_usage()
    sys.exit()

pluginName = argv[1]
className = "{}Processor".format(pluginName)

os.mkdir(pluginName)
with open(pluginName+"/CMakeLists.txt", 'w') as f:
    print( "Writing "+pluginName+"/CMakeLists.txt ...")
    f.write( cmakelists_template.format( pluginName) )
    f.close()

with open(pluginName+"/{}.h".format(className), 'w') as f:
    print( "Writing "+pluginName+"/{}.h".format(className) +" ...")
    f.write( processor_sequentialroot_header_template.format(className) )
    f.close()

with open(pluginName+"/{}.cc".format(className), 'w') as f:
    print( "Writing "+pluginName+"/{}.cc".format(className) +" ...")
    f.write( processor_sequentialroot_implementation_template.format(className) )
    f.close()

print("""
Created plugin {0}.
Build with:

  cmake -S {0} -B {0}/build
  cmake --build {0}/build --target install
""".format(pluginName))