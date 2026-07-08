#!/bin/bash
#$ -N rl6_headline
#$ -o /net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt/local-mode-validation/rl6-data/rl6_headline.out
#$ -e /net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt/local-mode-validation/rl6-data/rl6_headline.err
#$ -cwd
#$ -V
#$ -j n
#$ -l h_rt=86400
#$ -l m_mem_free=64G
#$ -m n
# Brief 26_0610-068 (R-L.6) HEADLINE: end-to-end total-RSS win through cmalign --ckpt in
# LOCAL (default) mode at wide-band scale, mirroring 056's global CH479288 measurement
# (global: stock 10555 MB -> ckpt 1406 MB = 7.5x).  Here: LOCAL config (no -g).

BIN=/net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt/src/cmalign
CM=/net/intdev/oblast01/infernal/notebook/26_0610_inf_dp_mem_cost/benchmark-runs/stage0_cyk_baseline/LSU_rRNA_eukarya.cm
SEQ=/net/intdev/oblast01/infernal/notebook/26_0610_inf_dp_mem_cost/subagent-work/040/seqs/CH479288.fa
OUT=/net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt/local-mode-validation/rl6-data
C="--cpu 0 --mxsize 60000"

echo "############ HOST: $(hostname)  DATE: $(date) ############"
echo "CM=$CM"; echo "SEQ=$SEQ (LSU-divergent CH479288, 2890 nt)"; echo

echo "==================== A: cmalign --ckpt (LOCAL, default trunc) ===================="
INFERNAL_CKPT_VERBOSE=1 /usr/bin/time -v $BIN --ckpt $C -o $OUT/headline_ckpt.sto $CM $SEQ \
   > $OUT/headline_ckpt.stdout 2> $OUT/headline_ckpt.time
echo "ckpt exit=$?"
grep -E "Maximum resident|Elapsed" $OUT/headline_ckpt.time
grep -E "engaged" $OUT/headline_ckpt.time | head
echo

echo "==================== B: cmalign STOCK (LOCAL, default trunc, no --ckpt) ===================="
/usr/bin/time -v $BIN $C -o $OUT/headline_stock.sto $CM $SEQ \
   > $OUT/headline_stock.stdout 2> $OUT/headline_stock.time
echo "stock exit=$?"
grep -E "Maximum resident|Elapsed" $OUT/headline_stock.time
echo

echo "==================== SUMMARY ===================="
CK=$(grep "Maximum resident" $OUT/headline_ckpt.time  | grep -oE "[0-9]+")
ST=$(grep "Maximum resident" $OUT/headline_stock.time | grep -oE "[0-9]+")
echo "ckpt  max RSS = ${CK} kB"
echo "stock max RSS = ${ST} kB"
if [ -n "$CK" ] && [ -n "$ST" ] && [ "$CK" -gt 0 ]; then
  awk -v c=$CK -v s=$ST 'BEGIN{printf "LOCAL total-RSS win = %.2fx  (stock %.0f MB -> ckpt %.0f MB)\n", s/c, s/1024, c/1024}'
fi
echo "############ DONE $(date) ############"
