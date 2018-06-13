#!/bin/sh

. ./globalvars.sh

make -C $DIR -j8 $BIN && mpirun -n $NP $DIR/$BIN \
	--dir $OUTDIR \
	--load $OUTDIR/rf1 \
	--prefix rf2 --size 4096 --save \
	--monf ./data/rf1.pat \
	--stim ./data/shapes.pat \
	--wie 0.2 --wee 0.1 --wext 0.1 \
	--simtime $SIMTIME --tauf $TAUF --taud $TAUD \
	--intsparse $INTSPARSENESS \
	--extsparse 0.10 \
	--off 5.0 --on 0.2 \
	--beta $BETA --eta $ETA --bgrate $BGRATE --scale $SCALE --weight_a $WEIGHTA --alpha $ALPHA --delta 0.02

