#/bin/bash
# get MER
cat $2/*out | sort -g | perl rmark-mer.pl $1.ppos > $2/$2.mer
# get hitlist
cat $2/*out | sort -g | perl rmark-hitlist.pl $1.ppos > $2/$2.hitlist
# get running time
ls $2/*.search | perl rmark-time.pl > $2/$2.time

# copy files to cwd
cp $2/$2.mer ./
cp $2/$2.hitlist ./
cp $2/$2.time ./

# summarize files to stdout
echo -n $2 "(positives only)"
grep ummary $2.mer
tail -n1 $2/$2.time
echo ' '

