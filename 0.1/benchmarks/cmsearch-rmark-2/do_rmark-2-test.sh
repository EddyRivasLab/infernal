# do_rmark-2-test.sh: Perform a trial run of the rmark-2 benchmark.
# See 00README for more information.
#
#
# The current directory should have the following files: sre.pl,
# rmark.pl, infernal-given-cm.rmm, infernal.pm, infernal2glbf.pl, and
# rmark_process_glbf.pl.
#
# The file inf-1p0.rmk should be in a subdir called rmk_files/ 
#
# Also, the current dir should have a subdir called rmark-2-test/
# containing the following files: rmark-2-test.idx, rmark-2-test.ebd,
# rmark-2-test.fa, RF00031.test, RF00031.ali, RF00031.cm,
# RF00031.idx, RF00177.test, RF00177.ali, RF00177.cm, and
# RF00177.idx.  
#
# Step 1: run the benchmark, outputting rmark-2-test_out.glbf 
#         and rmark-test_out.time
perl rmark.pl -M rmark-2-test infernal.rmm rmk_files/inf-1p0.rmk rmark-2-test/ rmark-2-test/rmark-2-test.idx rmark-2-test.fa rmark-2-test_out
#
# Step 2: create rmark-2-test_out.fam, rmark-2-test_out.all, and 
#         rmark-2-test_out.roc from the rmark-2-test_out.glbf file.
#         
perl rmark_process_glbf.pl E infernal.rmm rmk_files/inf-1p0.rmk rmark-2-test/ rmark-2-test/rmark-2-test.idx rmark-2-test rmark-2-test_out.glbf rmark-2-test_out

