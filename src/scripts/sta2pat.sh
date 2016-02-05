IN=$1
OUT=$2

# threshold is simply the mean
# THR=`awk '{sum+=$2} END {print 1.0*sum/NR}' $IN`
THR=10

echo "Threshold is $THR"

awk "{ if ( NF == 0 ) print ; else { if ( \$2>$THR ) print }}" $IN > $OUT
