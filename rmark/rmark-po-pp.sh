#/bin/bash
# get MER
cat $2/*out | sort -g | perl mer.pl $1.ppos > $2/$2.mer
# get running time
grep "^# CPU" $2/*.search | cpu2hrs.pl > $2/$2.time

# copy files to cwd
cp $2/$2.mer ./
cp $2/$2.time ./

# summarize files to stdou
echo $2 "(positives only)"
grep ummary $2.mer
cat $2.time
echo ' '

