#!/bin/sh

OUTDIR="$HOME/data/sim"
RECFILE="./data/rf_discsR8.mtx"

DIR="."
BIN="sim_rc_p10c"
NP=4
BETA="0.04"
ETA="1e-3"
ALPHA=4
WEIGHTA="0.0"
XI="1.0"
BGRATE="5"
SCALE="25"
TAUD="0.15" # Was 0.2 (200ms) in original sim, but 0.15 (150ms) works better

