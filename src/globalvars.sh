#!/bin/sh

OUTDIR="$HOME/data/sim"
RECFILE="./data/rf_discsR8.mtx"

DIR="."
# BIN="sim_rc_p10c"
BIN="sim_rc_p11"
NP=4
BETA="0.05"
ETA="1e-3"
ALPHA=4
WEIGHTA="0.0"
XI="0.5"
BGRATE="5"
SCALE="25"
INTSPARSENESS=0.05
SIMTIME=1800
TAUF="0.6" # Was 0.6 (600ms) here we are trying something shorter
TAUD="0.15" # Was 0.2 (200ms) in original sim, but this seemingly works better

