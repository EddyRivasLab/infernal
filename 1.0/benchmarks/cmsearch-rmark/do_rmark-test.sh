# do_rmark-test.sh: Perform a trial run of the rmark benchmark.
# See 00README for more information.
#
# The current directory should have the following files:
# sre.pl, rmark.pl, infernal.rmm, infernal.pm, and infernal2glbf.pl.
# The file inf_qdb-71.rmk should be in a subdir called rmk_files/
#
# Also, the current dir should have a subdir called rmark-test/ containing the 
# following files: 
# rmark-test.idx, rmark-test.ebd, rmark-test.fa, RF00005.test, RF00005.ali,
# RF00005.idx, RF00031.test, RF00031.ali, RF00031.idx.
#
# Step 1: run the benchmark, outputting rmark-test_out.glbf 
#         and rmark-test_out.time
perl rmark.pl -B 8 infernal.rmm rmk_files/inf_qdb-71.rmk rmark-test/ rmark-test/rmark-test.idx rmark-test.fa rmark-test_out
#
# Step 2: create rmark-test_out.fam, rmark-test_out.all, and 
#         rmark-test_out.roc from the rmark-test_out.glbf file.
#         
perl rmark_process_glbf.pl B infernal.rmm rmk_files/inf_qdb-71.rmk rmark-test/ rmark-test/rmark-test.idx rmark-test rmark-test_out.glbf rmark-test_out

