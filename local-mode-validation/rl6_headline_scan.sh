#!/bin/bash
#$ -N rl6_hl_scan
#$ -o /net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt/local-mode-validation/rl6-data/rl6_headline_scan.out
#$ -e /net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt/local-mode-validation/rl6-data/rl6_headline_scan.err
#$ -cwd
#$ -V
#$ -j n
#$ -l h_rt=86400
#$ -l m_mem_free=64G
#$ -l h_vmem=64G
#$ -m n
# Brief 26_0610-068 (R-L.6) HEADLINE search: find a genome-scale LOCAL case where the
# CM-DP cube dominates total RSS (so the sqrt(M) win shows end-to-end), by
# cross-aligning divergent caliciviruses to calici-NC_001959.cm (default local
# truncated).  Phase 1: --ckpt scan of cube sizes (cheap, sqrt(M)).  Phase 2:
# stock-vs-ckpt total-RSS on the largest-cube genome.

BIN=/net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt/src/cmalign
CM=/net/intdev/oblast01/infernal/notebook/26_0610_inf_dp_mem_cost/benchmark-runs/segment-census/cms/vadr/calici-NC_001959.cm
GEN=/net/intdev/oblast01/infernal/notebook/26_0610_inf_dp_mem_cost/benchmark-runs/checkpoint-proto/genomes
OUT=/net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt/local-mode-validation/rl6-data
C="--cpu 0 --fixedtau --tau 1e-7 --mxsize 60000"

echo "############ HOST: $(hostname)  DATE: $(date) ############"
echo "### Phase 1: --ckpt cube scan (calici-NC_001959.cm x each divergent calicivirus, LOCAL trunc)"
best=""; bestmb=0
printf "%-14s %-6s %-8s %-12s %-12s\n" genome mode cubeMB ckptPeakMB win
for f in $GEN/NC_*.fa; do
  g=$(basename $f .fa)
  line=$(INFERNAL_CKPT_VERBOSE=1 $BIN --ckpt $C -o /dev/null $CM $f 2>&1 | grep -m1 "engaged")
  mode=$(echo "$line"  | grep -oE "mode=[JLRT]" | cut -d= -f2)
  cube=$(echo "$line"  | grep -oE "full-cube\(2x\)=[0-9.]+" | cut -d= -f2)
  peak=$(echo "$line"  | grep -oE "peak=[0-9.]+" | head -1 | cut -d= -f2)
  win=$(echo "$line"   | grep -oE "win~[0-9.]+x")
  [ -z "$cube" ] && cube=0
  printf "%-14s %-6s %-8s %-12s %-12s\n" "$g" "${mode:-?}" "$cube" "$peak" "$win"
  # track max cube (bash float compare via awk)
  bigger=$(awk -v a=$cube -v b=$bestmb 'BEGIN{print (a>b)?1:0}')
  if [ "$bigger" = "1" ]; then bestmb=$cube; best=$g; bestf=$f; fi
done
echo
echo "### Largest-cube genome: $best  (full-cube(2x)=${bestmb} MB)"
echo
echo "### Phase 2: stock-vs-ckpt total RSS on $best (LOCAL, default trunc)"
INFERNAL_CKPT_VERBOSE=1 /usr/bin/time -v $BIN --ckpt --cpu 0 --fixedtau --tau 1e-7 --mxsize 60000 -o $OUT/hl_best_ckpt.sto $CM $bestf > $OUT/hl_best_ckpt.stdout 2> $OUT/hl_best_ckpt.time
echo "  ckpt exit=$?"; grep -E "Maximum resident|Elapsed" $OUT/hl_best_ckpt.time; grep -m1 engaged $OUT/hl_best_ckpt.time
/usr/bin/time -v $BIN --cpu 0 --fixedtau --tau 1e-7 --mxsize 60000 -o $OUT/hl_best_stock.sto $CM $bestf > $OUT/hl_best_stock.stdout 2> $OUT/hl_best_stock.time
echo "  stock exit=$?"; grep -E "Maximum resident|Elapsed" $OUT/hl_best_stock.time
CK=$(grep "Maximum resident" $OUT/hl_best_ckpt.time  | grep -oE "[0-9]+")
ST=$(grep "Maximum resident" $OUT/hl_best_stock.time | grep -oE "[0-9]+")
echo "  ckpt=${CK}kB stock=${ST}kB"
[ -n "$CK" ] && [ -n "$ST" ] && [ "$CK" -gt 0 ] && awk -v c=$CK -v s=$ST 'BEGIN{printf "  >>> LOCAL total-RSS win = %.2fx  (stock %.0f MB -> ckpt %.0f MB)\n", s/c, s/1024, c/1024}'
echo "############ DONE $(date) ############"
