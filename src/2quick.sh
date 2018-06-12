#!/bin/sh

. ./globalvars.sh

make -C $DIR -j8 $BIN && mpirun -n $NP $DIR/$BIN \
	--dir $OUTDIR \
	--load $OUTDIR/rf1 \
	--prefix rf2 --size 4096 --save \
	--monf ./data/rf1.pat \
	--stim ./data/shapes.pat \
	--wie 0.15 --wee 0.1 --wext 0.1 \
	--simtime 10 --tauf 0.6 --taud $TAUD \
	--extsparse 0.05 \
	--off 20.0 --on 1.0 \
	--beta $BETA --eta $ETA --bgrate $BGRATE --scale $SCALE --weight_a $WEIGHTA --alpha $ALPHA --delta 0.02

