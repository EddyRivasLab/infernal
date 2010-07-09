#/bin/bash
# get MER
cat $2/*out | sort -g | perl rmark-mer.pl $1.pos > $2/$2.mer
# get running time
grep "^# CPU" $2/*.search | cpu2hrs.pl > $2/$2.time
# get ROC
cat $2/*out | sort -g | rocplot-rmark $1 - > $2/$2.xy

# copy files to cwd
cp $2/$2.mer ./
cp $2/$2.time ./
cp $2/$2.xy ./

# summarize files to stdout
echo $2
grep ummary $2/$2.mer
cat $2/$2.time
echo ' '
