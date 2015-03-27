#!/bin/sh

. ./globalvars.sh

make -C $DIR -j8 $BIN && mpirun -n 4 $DIR/$BIN \
	--dir $OUTDIR \
	--prefix rf1 --size 4096 --save \
	--monf rf1.pat \
	--recfile $RECFILE --xi $XI \
	--stim ~/stim/shapes.pat \
	--wie 0.15 --wee 0.1 --wext 0.1 \
	--simtime 3600 --tauf 0.6 \
	--extsparse 0.05 \
	--off 2.0 --on 1.0 \
	--beta $BETA --eta $ETA --bgrate $BGRATE --scale $SCALE --weight_a $WEIGHTA --alpha $ALPHA --delta 0.02

