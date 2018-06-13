#!/bin/sh

. ./globalvars.sh

make -C $DIR -j8 $BIN && mpirun -n $NP $DIR/$BIN \
	--dir $OUTDIR \
	--load $OUTDIR/rf2 \
	--prefix rf3 --size 4096 --save \
	--monf ./data/rf1.pat \
	--stim ./data/shapes_cues.pat \
	--wie 0.15 --wee 0.1 --wext 0.1 \
	--simtime $SIMTIME --tauf $TAUF --taud $TAUD \
	--intsparse $INTSPARSENESS \
	--extsparse 0.10 \
	--off 5.0 --on 0.2 \
	--beta $BETA --eta $ETA --bgrate $BGRATE --scale 50 --weight_a $WEIGHTA --alpha $ALPHA --delta 0.02



