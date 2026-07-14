#!/bin/bash
# brief 26_0610-093: global-mode (-g) base-vs-post regression on struct+aggr+bps0.
set -u
cd /net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-trcyk-dnc-localbug/local-mode-validation
SP=/tmp/claude-12754/-home-nawrocke-notebook-home-26-0610-inf-dp-mem-cost/bdae483c-3bda-4c59-97d8-561ab1d5b046/scratchpad
run() { # $1=binary $2=tag
  local B=$1 T=$2
  for cm in tRNA ar45 BjrC174 IRES_Picorna CsrB; do
    $B -g r090-struct/${cm}.cm r090-struct/${cm}_panel.fa > $SP/g_struct_${cm}_${T}.log 2>&1 &
  done
  for cm in ar45 BjrC174 IRES_Picorna CsrB tRNA; do
    $B -g r090-struct/${cm}.cm r090-struct/${cm}_aggr.fa > $SP/g_aggr_${cm}_${T}.log 2>&1 &
  done
  for c in 0 1 2 3 4 5 6 7; do
    $B -g rl2-data/matl300.cm /tmp/mp_chunk${c}.fa > $SP/g_bps0_chunk${c}_${T}.log 2>&1 &
  done
  wait
  cat $SP/g_struct_*_${T}.log > $SP/g_struct_${T}.log
  cat $SP/g_aggr_*_${T}.log   > $SP/g_aggr_${T}.log
  cat $SP/g_bps0_chunk*_${T}.log > $SP/g_bps0_${T}.log
}
run ./r093_argmax_base base
run ./r093_argmax      post
echo "GLOBAL DONE"
for panel in struct aggr bps0; do
  echo "$panel: base FAIL=$(grep -c '^FAIL' $SP/g_${panel}_base.log) post FAIL=$(grep -c '^FAIL' $SP/g_${panel}_post.log)  base CRASH=$(grep -c '^CRASH' $SP/g_${panel}_base.log) post CRASH=$(grep -c '^CRASH' $SP/g_${panel}_post.log)"
done
