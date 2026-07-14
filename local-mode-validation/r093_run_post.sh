#!/bin/bash
# brief 26_0610-093: run all post-fix panels with r093_argmax (fixed lib), chunked.
set -u
cd /net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-trcyk-dnc-localbug/local-mode-validation
SP=/tmp/claude-12754/-home-nawrocke-notebook-home-26-0610-inf-dp-mem-cost/bdae483c-3bda-4c59-97d8-561ab1d5b046/scratchpad
D=./r093_argmax
GF="${1:-}"; SUF=""; [ "$GF" = "-g" ] && SUF="_global"

# structured (240 each)
for cm in tRNA ar45 BjrC174 IRES_Picorna CsrB; do
  ( echo "=== $cm ==="; $D $GF r090-struct/${cm}.cm r090-struct/${cm}_panel.fa ) > $SP/post_struct_${cm}${SUF}.log 2>&1 &
done
# aggressive (640 each)
for cm in ar45 BjrC174 IRES_Picorna CsrB tRNA; do
  ( echo "=== $cm ==="; $D $GF r090-struct/${cm}.cm r090-struct/${cm}_aggr.fa ) > $SP/post_aggr_${cm}${SUF}.log 2>&1 &
done
# bps0 chunked (8)
if [ -z "$SUF" ]; then
  for c in 0 1 2 3 4 5 6 7; do
    $D rl2-data/matl300.cm /tmp/mp_chunk${c}.fa > $SP/post_bps0_chunk${c}.log 2>&1 &
  done
fi
wait
cat $SP/post_struct_tRNA${SUF}.log $SP/post_struct_ar45${SUF}.log $SP/post_struct_BjrC174${SUF}.log $SP/post_struct_IRES_Picorna${SUF}.log $SP/post_struct_CsrB${SUF}.log > $SP/post_struct${SUF}.log
cat $SP/post_aggr_ar45${SUF}.log $SP/post_aggr_BjrC174${SUF}.log $SP/post_aggr_IRES_Picorna${SUF}.log $SP/post_aggr_CsrB${SUF}.log $SP/post_aggr_tRNA${SUF}.log > $SP/post_aggr${SUF}.log
[ -z "$SUF" ] && cat $SP/post_bps0_chunk*.log > $SP/post_bps0_all.log
echo "POST${SUF} DONE"
echo "struct: PASS=$(grep -c '^PASS' $SP/post_struct${SUF}.log) FAIL=$(grep -c '^FAIL' $SP/post_struct${SUF}.log) CRASH=$(grep -c '^CRASH' $SP/post_struct${SUF}.log)"
echo "aggr:   PASS=$(grep -c '^PASS' $SP/post_aggr${SUF}.log) FAIL=$(grep -c '^FAIL' $SP/post_aggr${SUF}.log) CRASH=$(grep -c '^CRASH' $SP/post_aggr${SUF}.log)"
[ -z "$SUF" ] && echo "bps0:   PASS=$(grep -c '^PASS' $SP/post_bps0_all.log) FAIL=$(grep -c '^FAIL' $SP/post_bps0_all.log) CRASH=$(grep -c '^CRASH' $SP/post_bps0_all.log)"
