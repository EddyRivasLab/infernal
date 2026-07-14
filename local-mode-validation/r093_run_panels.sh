#!/bin/bash
# brief 26_0610-093 panel runner. Usage: r093_run_panels.sh <tag> [-g]
# Runs bps0 + structured + aggressive panels with r090_forked_drv, writing
# per-panel logs tagged <tag>. Runs CMs concurrently; waits for all.
set -u
TAG=$1; shift || true
GFLAG="${1:-}"
D=${DRV:-./r090_forked_drv}
SUF=""
[ "$GFLAG" = "-g" ] && SUF="_global"

# bps=0 (single CM, 1800 seqs) -- chunk into 8 for speed
$D $GFLAG rl2-data/matl300.cm r088-data/matl300_panel.fa > r093_bps0_${TAG}${SUF}.log 2>&1 &

# structured panel (240 each): order tRNA ar45 BjrC174 IRES_Picorna CsrB
for cm in tRNA ar45 BjrC174 IRES_Picorna CsrB; do
  ( echo "=== $cm ==="; $D $GFLAG r090-struct/${cm}.cm r090-struct/${cm}_panel.fa ) > r093_struct_${cm}_${TAG}${SUF}.log 2>&1 &
done

# aggressive panel (640 each): order ar45 BjrC174 IRES_Picorna CsrB tRNA
for cm in ar45 BjrC174 IRES_Picorna CsrB tRNA; do
  ( echo "=== $cm ==="; $D $GFLAG r090-struct/${cm}.cm r090-struct/${cm}_aggr.fa ) > r093_aggr_${cm}_${TAG}${SUF}.log 2>&1 &
done

wait
# assemble combined logs
cat r093_struct_tRNA_${TAG}${SUF}.log r093_struct_ar45_${TAG}${SUF}.log r093_struct_BjrC174_${TAG}${SUF}.log r093_struct_IRES_Picorna_${TAG}${SUF}.log r093_struct_CsrB_${TAG}${SUF}.log > r093_struct_${TAG}${SUF}.log
cat r093_aggr_ar45_${TAG}${SUF}.log r093_aggr_BjrC174_${TAG}${SUF}.log r093_aggr_IRES_Picorna_${TAG}${SUF}.log r093_aggr_CsrB_${TAG}${SUF}.log r093_aggr_tRNA_${TAG}${SUF}.log > r093_aggr_${TAG}${SUF}.log
echo "TAG=$TAG${SUF} DONE"
echo "bps0:   PASS=$(grep -c '^PASS' r093_bps0_${TAG}${SUF}.log) FAIL=$(grep -c '^FAIL' r093_bps0_${TAG}${SUF}.log) CRASH=$(grep -c '^CRASH' r093_bps0_${TAG}${SUF}.log)"
echo "struct: PASS=$(grep -c '^PASS' r093_struct_${TAG}${SUF}.log) FAIL=$(grep -c '^FAIL' r093_struct_${TAG}${SUF}.log) CRASH=$(grep -c '^CRASH' r093_struct_${TAG}${SUF}.log)"
echo "aggr:   PASS=$(grep -c '^PASS' r093_aggr_${TAG}${SUF}.log) FAIL=$(grep -c '^FAIL' r093_aggr_${TAG}${SUF}.log) CRASH=$(grep -c '^CRASH' r093_aggr_${TAG}${SUF}.log)"
