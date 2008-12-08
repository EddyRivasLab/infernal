# See 00README for details.
# These 3 rmark_clusterfy calls will create 3 subdirectories in the current directory, named:
#  inf-72_realmark_out_dir/
#  inf-1p0_realmark_out_dir/
#  inf-1p0_nofilter_realmark_out_dir/
# 
# Inside each dir there will be two scripts *.com and *pp.sh. To submit the jobs 
# to the cluster, run the *.com script. Then wait for all the jobs to finish and
# run the *pp.sh script.

perl rmark_clusterfy.pl            -B 8 infernal.rmm rmk_files/inf-72.rmk           rmark-2/  rmark-2/rmark-2.idx rmark-2 inf-72
perl rmark_clusterfy.pl -M cms-1p0 -E 1 infernal.rmm rmk_files/inf-1p0.rmk          rmark-2/  rmark-2/rmark-2.idx rmark-2 inf-1p0
perl rmark_clusterfy.pl -M cms-1p0 -E 1 infernal.rmm rmk_files/inf-1p0_nofilter.rmk rmark-2/  rmark-2/rmark-2.idx rmark-2 inf-1p0_nofilter
