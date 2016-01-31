#!/bin/sh

. ./globalvars.sh

make -C $DIR -j8 $BIN && mpirun -n 4 $DIR/$BIN \
	--dir $OUTDIR \
	--load $OUTDIR/rf2 \
	--prefix rf3 --size 4096 --save \
	--monf ./data/rf1.pat \
	--stim ./data/shapes_cues.pat \
	--wie 0.15 --wee 0.1 --wext 0.1 \
	--simtime 3600 --tauf 0.6 \
	--extsparse 0.05 \
	--off 20.0 --on 0.5 \
	--beta $BETA --eta $ETA --bgrate $BGRATE --scale $SCALE --weight_a $WEIGHTA --alpha $ALPHA --delta 0.02



