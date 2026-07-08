#!/bin/bash
#$ -N rl6_hl_genome
#$ -o /net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt/local-mode-validation/rl6-data/rl6_headline_genome.out
#$ -e /net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt/local-mode-validation/rl6-data/rl6_headline_genome.err
#$ -cwd
#$ -V
#$ -j n
#$ -l h_rt=86400
#$ -l m_mem_free=64G
#$ -l h_vmem=64G
#$ -m n
# Brief 26_0610-068 (R-L.6) HEADLINE (genome scale): cmalign --ckpt vs stock in LOCAL
# (default) mode on a genome-scale bps=0 VADR viral CM (calici NC_001959,
# M=22966, clen=7654) aligning its own 7654 nt genome.  Per R-L.2b the sqrt(M)
# CM-DP cube win here is ~27x; the default-truncated cube is the large one
# (10-15 GB per the composition memory note) so total-RSS win shows in trunc.

BIN=/net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt/src/cmalign
CM=/net/intdev/oblast01/infernal/notebook/26_0610_inf_dp_mem_cost/benchmark-runs/segment-census/cms/vadr/calici-NC_001959.cm
SEQ=/net/intdev/oblast01/infernal/notebook/26_0610_inf_dp_mem_cost/benchmark-runs/checkpoint-proto/genomes/NC_001959.fa
OUT=/net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt/local-mode-validation/rl6-data
C="--cpu 0 --mxsize 60000"

echo "############ HOST: $(hostname)  DATE: $(date) ############"
echo "CM=$CM (calici NC_001959, bps=0, M=22966)"; echo "SEQ=$SEQ (7654 nt)"; echo

runpair() {  # modetag  modeflag
  local mt=$1 mf=$2
  echo "==================== $mt : --ckpt (LOCAL) ===================="
  INFERNAL_CKPT_VERBOSE=1 /usr/bin/time -v $BIN $mf --ckpt $C -o $OUT/hg_${mt}_ckpt.sto $CM $SEQ \
     > $OUT/hg_${mt}_ckpt.stdout 2> $OUT/hg_${mt}_ckpt.time
  echo "  ckpt exit=$?"; grep -E "Maximum resident|Elapsed" $OUT/hg_${mt}_ckpt.time; grep -m1 "engaged" $OUT/hg_${mt}_ckpt.time
  echo "==================== $mt : STOCK (LOCAL) ===================="
  /usr/bin/time -v $BIN $mf $C -o $OUT/hg_${mt}_stock.sto $CM $SEQ \
     > $OUT/hg_${mt}_stock.stdout 2> $OUT/hg_${mt}_stock.time
  echo "  stock exit=$?"; grep -E "Maximum resident|Elapsed" $OUT/hg_${mt}_stock.time
  local CK=$(grep "Maximum resident" $OUT/hg_${mt}_ckpt.time  | grep -oE "[0-9]+")
  local ST=$(grep "Maximum resident" $OUT/hg_${mt}_stock.time | grep -oE "[0-9]+")
  echo "  --> $mt  ckpt=${CK}kB  stock=${ST}kB"
  if [ -n "$CK" ] && [ -n "$ST" ] && [ "$CK" -gt 0 ] && [ "$ST" -gt 0 ]; then
    awk -v c=$CK -v s=$ST -v m=$mt 'BEGIN{printf "  >>> %s LOCAL total-RSS win = %.2fx  (stock %.0f MB -> ckpt %.0f MB)\n", m, s/c, s/1024, c/1024}'
  else
    echo "  >>> $mt : one side failed (likely stock OOM / mxsize exceeded) -> ckpt succeeds where stock does not"
  fi
  echo
}

runpair trunc    ""
runpair notrunc  "--notrunc"
echo "############ DONE $(date) ############"
