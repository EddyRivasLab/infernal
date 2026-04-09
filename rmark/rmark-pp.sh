#/bin/bash

# get running time
ls $2/*.time | perl $4/rmark-time.pl > $2/$2.time

# Sort hit list once and reuse (was being recomputed 4x)
SORTED=$2/$2.em$3.sorted
cat $2/*out | perl $4/rmark-multiply-evalues.pl $3 | sort -g > $SORTED

# get MER
perl $4/rmark-mer.pl $1.pos $2/$2.time < $SORTED > $2/$2.em$3.mer

# .xy ROC plot skipped (never used)

# get mer from rmark-rocplot
$4/rmark-rocplot -N 10000 --mer --seed 181 $1 $SORTED > $2/$2.em$3.bmer
# get numbers of false negatives and false positives at E-threshold of 0.1 from rmark-rocplot (after E-value inflation)
$4/rmark-rocplot -N 10000 --Ethresh 0.1 --seed 181 $1 $SORTED > $2/$2.em$3.bEthresh

rm -f $SORTED

# summarize files to stdout
echo -n $2 | awk '{printf("%-70s  ", $0)}' > $2/$2.em$3.sum
echo -n $3 | awk '{printf("%4s  ", $0)}' >> $2/$2.em$3.sum
cat $2/$2.em$3.bmer     | awk '{printf("MER:  %5s  %5s  %5s  %5s   ", $3, $7, $11, $15)}' >> $2/$2.em$3.sum
cat $2/$2.em$3.bEthresh | awk '{printf("ETHRESH:0.1  %5s  %5s   ", $5, $9)}' >> $2/$2.em$3.sum
grep ummary $2/$2.em$3.mer | awk '{print $7}' >> $2/$2.em$3.sum
cp $2/$2.em$3.sum ./
cat $2.em$3.sum
