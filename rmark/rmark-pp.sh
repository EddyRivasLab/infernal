#/bin/bash
# get running time
ls $2/*.time   | perl rmark-time.pl > $2/$2.time
# get MER
cat $2/*out | sort -g | perl rmark-mer.pl $1.pos $2/$2.time > $2/$2.mer
# get ROC
cat $2/*out | sort -g | rmark-rocplot --seed 181 $1 - > $2/$2.xy

# copy files to cwd
cp $2/$2.mer ./
cp $2/$2.time ./
cp $2/$2.xy ./

# summarize files to stdout
echo -n $2
grep ummary $2/$2.mer
tail -n1 $2/$2.time
echo ' '
