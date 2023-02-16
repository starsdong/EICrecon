#!/bin/bash

source /opt/detector/setup.sh

ddsim --steeringFile mysteer_etarange.py --compactFile $DETECTOR_PATH/epic.xml --outputFile output.edm4hep.root

source ../../../bin/eicrecon-this.sh

eicrecon \
-Pplugins=dump_flags,track_qa \
-Ppodio:output_file=eicrecon_out.root \
-Ppodio:output_include_collections=MCParticles,CentralTrackSeedingResults \
-Ptrack_qa:LogLevel=trace \
-Pjana:nevents=1000 \
-Pdd4hep:xml_files=epic.xml \
output.edm4hep.root | tee eicrecon_out.dat
