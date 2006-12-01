# See 00README for details.
# These 6 rmark_clusterfy calls will create 6 subdirectories in the current directory, named:
#  inf-55_rmark-1_out_dir/
#  inf_p1_noent-71_rmark-1_out_dir/
#  inf_p1-71_rmark-1_out_dir/
#  inf_noent-71_rmark-1_out_dir/
#  inf-71_rmark-1_out_dir/
#  inf_qdb-71_rmark-1_out_dir/
# 
# Inside each dir there will be two scripts *.com and *pp.sh. To submit the jobs 
# to the cluster, run the *.com script. Then wait for all the jobs to finish and
# run the *pp.sh script.

perl rmark_clusterfy.pl -B 8 infernal_55_W.rmm rmk_files/inf-55.rmk          rmark-1/  rmark-1/rmark-1.idx rmark-1 inf-55
perl rmark_clusterfy.pl -B 8 infernal.rmm      rmk_files/inf_p1_noent-71.rmk rmark-1/  rmark-1/rmark-1.idx rmark-1 inf_p1_noent-71
perl rmark_clusterfy.pl -B 8 infernal.rmm      rmk_files/inf_p1-71.rmk       rmark-1/  rmark-1/rmark-1.idx rmark-1 inf_p1-71
perl rmark_clusterfy.pl -B 8 infernal.rmm      rmk_files/inf_noent-71.rmk    rmark-1/  rmark-1/rmark-1.idx rmark-1 inf_noent-71
perl rmark_clusterfy.pl -B 8 infernal.rmm      rmk_files/inf-71.rmk          rmark-1/  rmark-1/rmark-1.idx rmark-1 inf-71
perl rmark_clusterfy.pl -B 8 infernal.rmm      rmk_files/inf_qdb-71.rmk      rmark-1/  rmark-1/rmark-1.idx rmark-1 inf_qdb-71

