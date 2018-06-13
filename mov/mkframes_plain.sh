#!/bin/bash

OUTPUTDIR=$1
PREFIX=$2
TIMEOFFSET=$3

if [ "$4" != "" ]; then
	A=$4
else
	A=4
fi

if [ "$5" != "" ]; then
	IMGFILE=$5
else
	IMGFILE="cues.img"
fi


for SLICE in {0..12}; do
	for FRAME in {0..249}; do
		OUTPUTFILE=`printf "%s/mov%09d.%03d.png" $OUTPUTDIR $SLICE $FRAME`
		echo $OUTPUTFILE
cat << __EOF | gnuplot &
slice = $SLICE
subframe = $FRAME
framerate = 25
slicesize = 10.0
color_break_point = 0.1
tmpdir = '$OUTPUTDIR'
prefix = '$PREFIX'

inputActivityFile = sprintf('%s/%s.s.%09d.%03d.act',tmpdir,prefix,slice,subframe)
outputActivityFile = sprintf('%s/%s.e.%09d.%03d.act',tmpdir,prefix,slice,subframe)
rasFile = sprintf('%s/%s.e.%09d.ras',tmpdir,prefix,slice)
stimFile = sprintf('%s/%s.0.stimtimes',tmpdir,prefix)
pactFile = sprintf('%s/%s.0.%09d.pact',tmpdir,prefix,slice)

# every
a = $A
imgfile = '$IMGFILE'


# frameselect
frame_select(x) = x <= slice*slicesize+1.0*subframe/framerate ? x : 1/0

# set hot palette
# set palette rgb 21,22,23
set palette defined (0 0 0 0, color_break_point 239./255 69./255 35./255, 1 254./255 228./255 102./255)


# Line styles
set style line 1 lc rgb '#E41A1C' # red
set style line 2 lc rgb '#377EB8' # blue
set style line 3 lc rgb '#4DAF4A' # green
set style line 4 lc rgb '#984EA3' # purple
set style line 5 lc rgb '#FF7F00' # orange

unset border
unset xtics
unset ytics
unset key
unset colorbox

set out '$OUTPUTFILE'
set term pngcairo size 1024,576 fontscale 1.5

set cbrange [0:40]
set cbtics 20

set multiplot

set tmargin screen 0.93
set bmargin screen 0.73
set lmargin screen 0.15
set rmargin screen 0.26

# do for [i=1:7] {
# 	set arrow i from screen 0.28, graph 0.5 to screen 0.38, screen 0.55+i*0.05 lc rgb "grey" lw 2
# }

set title 'Input units' offset 0, -0.5
set colorbox
set cblabel 'Firing rate [Hz]' offset 1,0

plot inputActivityFile with image

unset arrow
unset colorbox
unset cblabel 


set tmargin screen 0.93
set bmargin screen 0.73
set lmargin screen 0.45
set rmargin screen 0.56
# set cbrange [0:4]


set cbrange [0:*]
set palette defined ( 0 0 0 0, 1 1 1 1 )
# set label 1 at screen 0.45, graph 0.9 "Last stim.: " right
set title 'Last stimulus' offset 0, -0.5

# unset title 
cmd = sprintf("awk -v timepoint=%f -f laststim.awk < %s",slice*slicesize+1.0*subframe/framerate,stimFile)
stimno = system(cmd)+0
plot imgfile index stimno with image 

unset label 1


# Bottom plots with time series and RAS

unset title
set lmargin screen 0.15
set rmargin screen 0.98


set tmargin screen 0.62
set bmargin screen 0.60

set xrange [slice*slicesize:(slice+1)*slicesize]

set yrange [-1:1]

# color palette
set palette defined ( 0 '#E41A1C',\
	1 '#377EB8',\
	2 '#4DAF4A',\
	3 '#984EA3',\
	4 '#000000')
set cbrange [0:4]

set label 1 at screen 0.15, graph 0.5 "Stimulus: " right

plot stimFile every 1 using 1:(column(2)>0?0:1/0):(int(\$3/a)) with lines palette lw 4

unset label 1

set bmargin screen 0.30
set tmargin screen 0.59

set yrange [:500]

plot rasFile every 1 using (frame_select(\$1)):2 with dots lc -1


set bmargin screen 0.13
set tmargin screen 0.25

set border 2
set xtics out nomirror
set xtics 5
unset xtics

set ytics out nomirror
set ytics 20
set yrange [0:40]

# set xlabel 'Time [s]' offset 0,0.2
unset xlabel

set arrow 1 from graph 0.8, screen 0.1 to graph 0.9, screen 0.1 lw 5 nohead
set label 1 at graph 0.85, screen 0.07 sprintf('%ds',1.0*slicesize/10) center

timeoffsettext = sprintf("%imin",($TIMEOFFSET+slice*slicesize)/60) 
set label 2 at graph 0.0, graph -0.4 timeoffsettext left
# set label 2 at graph 0.0, screen 0.07 "$TIMEOFFSET" center

set ylabel 'Activity [Hz]'

plot for [i=1:4] pactFile using (frame_select(\$1)):(column(i+1)):(i-1) with lines ls i lw 2

unset multiplot
__EOF
	if ! ((FRAME % 4)); then
		wait
	fi
	# adds keys to all frames
	# exit 
done
wait
for FRAME in {0..249}; do
	OUTPUTFILE=`printf "%s/mov%09d.%03d.png" $OUTPUTDIR $SLICE $FRAME`
	composite -compose atop -geometry +80+300 nice_key.png $OUTPUTFILE $OUTPUTFILE &
	if ! ((FRAME % 4)); then
		wait
	fi
done
wait
done
