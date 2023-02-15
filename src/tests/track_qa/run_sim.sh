#!/bin/bash

source /opt/detector/setup.sh

ddsim --steeringFile mysteer.py --compactFile $DETECTOR_PATH/epic.xml --outputFile output.edm4hep.root

source /global/project/projectdirs/m3763/blianggi/DD4HEP/EICrecon/bin/eicrecon-this.sh

eicrecon \
-Pplugins=dump_flags,track_qa \
-Ppodio:output_file=eicrecon_out.root \
-Ptrack_qa:LogLevel=trace \
-Pjana:nevents=100 \
-Pdd4hep:xml_files=epic.xml \
output.edm4hep.root | tee eicrecon_out.dat
