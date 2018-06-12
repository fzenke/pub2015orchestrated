#!/bin/sh

. ./globalvars.sh

make -C $DIR -j8 $BIN && mpirun -n $NP valgrind $DIR/$BIN \
	--dir $OUTDIR \
	--prefix tst1 --size 4096 --save \
	--monf ./data/rf1.pat \
	--stim ./data/shapes.pat \
	--wie 0.2 --wee 0.1 --wext 0.1 \
	--simtime 5 --tauf 0.6 --taud $TAUD \
	--extsparse 0.05 \
	--off 10.0 --on 1.0 \
	--beta $BETA --eta $ETA --bgrate $BGRATE --scale $SCALE --weight_a $WEIGHTA --alpha $ALPHA --delta 0.02

make -C $DIR -j8 $BIN && mpirun -n $NP valgrind $DIR/$BIN \
	--dir $OUTDIR \
	--load $OUTDIR/tst1 \
	--prefix tst2 --size 4096 --save \
	--monf ./data/rf1.pat \
	--stim ./data/shapes.pat \
	--wie 0.2 --wee 0.1 --wext 0.1 \
	--simtime 5 --tauf 0.6 --taud $TAUD \
	--extsparse 0.05 \
	--off 10.0 --on 1.0 \
	--beta $BETA --eta $ETA --bgrate $BGRATE --scale $SCALE --weight_a $WEIGHTA --alpha $ALPHA --delta 0.02

