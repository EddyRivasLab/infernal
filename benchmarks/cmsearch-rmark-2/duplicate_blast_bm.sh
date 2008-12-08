# duplicate_blast_bm.sh: Duplicate the blast benchmark reported in the
# Infernal 1.0 applications note.  See 00README for more information.  
#
# IMPORTANT: To run this benchmark, the command 'blastn' must execute
# wu-blastn version BLASTN 2.0MP-WashU, and the command 'xdformat'
# must execute XDFORMAT-WashU 1.0.  
#
# The blast benchmark differs from the infernal benchmarks in that we
# use E-values to rank scores from blast, while we use bit-scores to
# rank scores from Infernal (E-values were not yet implemented in
# Infernal 0.71). This is the -E option to rmark.pl and the 'E' 
# argument to rmark_process_glbf.pl shown below.
# 
# The current directory should have the following files:
# sre.pl, rmark.pl, blast.rmm and a file called blast_w7.rmk
# should be in a subdir called rmk_files/
#
# Also, the current dir should have a subdir called rmark-2/ containing 
# the following files: 
# rmark-2.idx, rmark-2.ebd, rmark-2.fa, and 51 sets of RFXXXXX.test, RFXXXXX.ali,
# RFXXXXX.idx, where XXXXX is a family specific identifier (see 00README). 
#
#
# Step 1: run the benchmark, outputting blast.glbf 
#         and blast.time
perl rmark.pl -E 1 blast.rmm rmk_files/blast_w7.rmk rmark-2/ rmark-2/rmark-2.idx rmark-2.fa blast
#
# Step 3: create blast.fam, blast.all, and blast.roc from the 
#         blast.glbf file.
perl rmark_process_glbf.pl E blast.rmm rmk_files/blast_w7.rmk rmark-2/ rmark-2/rmark-2.idx rmark-2 blast.glbf blast 

