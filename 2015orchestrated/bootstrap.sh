#!/bin/sh

DIR="$HOME/tmp/sim"

./scripts/mktrigfile.awk $DIR/rf1.0.stimtimes > data/rf1.trig

aube --input $DIR/rf1.*.ras --from 3000 > $DIR/rf1.ras

scripts/rassta.py -t data/rf1.trig -w 0.5 -f $DIR/rf1.ras -s 3000 -m 3500 -q 0.0 -o data/rf1.sta

scripts/sta2pat.sh data/rf1.sta data/rf1.pat
