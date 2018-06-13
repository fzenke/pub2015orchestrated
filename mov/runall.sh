#!/bin/bash

INPUTDIR="$HOME/data/sim"

# Auryn +v0.7
RASBIN="aube --to 400 --input "
RASEXT="spk"

# Auryn v0.5
# RASBIN="sort -g -m "
# RASEXT="ras"

# Create tmpdir
TMPDIR=`mktemp -d`
echo "tmpdir=$TMPDIR"
trap "{ rm -rf $TMPDIR; }" EXIT
# clean up mov dir

PREFIX=rf1
$RASBIN $INPUTDIR/$PREFIX.*.e.$RASEXT > $TMPDIR/$PREFIX.e.ras &
$RASBIN $INPUTDIR/$PREFIX.*.s.$RASEXT > $TMPDIR/$PREFIX.s.ras &
cp $INPUTDIR/$PREFIX.0.stimtimes $INPUTDIR/$PREFIX.0.pact $TMPDIR &
wait
./toolchain.sh $TMPDIR $PREFIX 0 1 shapes.img

sleep 5

PREFIX=rf2
$RASBIN $INPUTDIR/$PREFIX.*.e.$RASEXT > $TMPDIR/$PREFIX.e.ras &
$RASBIN $INPUTDIR/$PREFIX.*.s.$RASEXT > $TMPDIR/$PREFIX.s.ras &
cp $INPUTDIR/$PREFIX.0.stimtimes $INPUTDIR/$PREFIX.0.pact $TMPDIR &
wait
./toolchain.sh $TMPDIR $PREFIX 1200 1 shapes.img

sleep 5

PREFIX=rf3
$RASBIN $INPUTDIR/$PREFIX.*.e.$RASEXT > $TMPDIR/$PREFIX.e.ras &
$RASBIN $INPUTDIR/$PREFIX.*.s.$RASEXT > $TMPDIR/$PREFIX.s.ras &
cp $INPUTDIR/$PREFIX.0.stimtimes $INPUTDIR/$PREFIX.0.pact $TMPDIR &
wait
./toolchain.sh $TMPDIR $PREFIX 4800 4 cues.img

sleep 5

# PREFIX=rf4
# $RASBIN $INPUTDIR/$PREFIX.*.e.$RASEXT > $TMPDIR/$PREFIX.e.ras &
# $RASBIN $INPUTDIR/$PREFIX.*.s.$RASEXT > $TMPDIR/$PREFIX.s.ras &
# cp $INPUTDIR/$PREFIX.0.stimtimes $INPUTDIR/$PREFIX.0.pact $TMPDIR &
# wait
# ./toolchain.sh $TMPDIR $PREFIX 7200 4 cues.img


