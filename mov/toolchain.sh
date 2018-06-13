#!/bin/bash

SLICESIZE=10 # seconds
# PREFIX="rf2"
TMPDIR=$1
PREFIX=$2
TIMEOFFSET=$3
A=$4
IMGFILE=$5


echo "Slicing ras and pat files ..."
python tiser_slice.py $TMPDIR/$PREFIX.e.ras $TMPDIR $SLICESIZE &
python tiser_slice.py $TMPDIR/$PREFIX.0.pact $TMPDIR $SLICESIZE &
wait

echo "Creating act files ..."
# python activity_movie_fluo.py $TMPDIR/$PREFIX.e.ras $TMPDIR $SLICESIZE 
python activity_movie.py $TMPDIR/$PREFIX.s.ras $TMPDIR $SLICESIZE 
wait

echo "Making frames ..."
./mkframes_plain.sh $TMPDIR $PREFIX $TIMEOFFSET $A $IMGFILE

echo "Encoding video ..."
mencoder mf://$TMPDIR/*.png -mf w=1280:h=720:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o $PREFIX.avi

