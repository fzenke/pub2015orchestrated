#!/bin/sh

# Regenerates the monitor file "rf1.pat" which 
# stores the four different selective subpopulations.

# Adjust this variable to the output path of your network simulation
DIR="$HOME/tmp/sim"


# Leave unchanged in most cases
T_START=3000
T_STOP=3500

# Create trigger file for PSTH computation
./scripts/mktrigfile.awk $DIR/rf1.0.stimtimes > data/rf1.trig

# Create human-readable ras file from binary output files 
# (Assumes the auryn tool aube is accessible within your PATH)
aube --input $DIR/rf1.*.e.bras --from $T_START --to $T_STOP > $DIR/rf1.e.ras

# Creates PSTH using our trigger file
scripts/rassta.py -t data/rf1.trig -w 0.5 -f $DIR/rf1.e.ras -s $T_START -m $T_STOP -q 0.0 -o data/rf1.sta

# Generate pattern file through threshold computation
scripts/sta2pat.sh data/rf1.sta data/rf1.pat
