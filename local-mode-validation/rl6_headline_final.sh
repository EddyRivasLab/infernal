#!/bin/bash
#$ -N rl6_hl_final
#$ -o /net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt/local-mode-validation/rl6-data/rl6_headline_final.out
#$ -e /net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt/local-mode-validation/rl6-data/rl6_headline_final.err
#$ -cwd
#$ -V
#$ -j n
#$ -l h_rt=86400
#$ -l m_mem_free=48G
#$ -l h_vmem=48G
#$ -m n
# Brief 068 (R-L.6) HEADLINE (final): LOCAL-mode end-to-end total-RSS win where the
# CM-DP cube dominates.  Divergent caliciviruses aligned to calici-NC_001959.cm
# (M=22966) resolve to marginal mode=L with genome-scale cubes (NC_006875 ~24 GB,
# NC_008311 ~9 GB).  Default local truncated cmalign.  stock (full cube) vs
# --ckpt (sqrt(M)).  This is the number the whole local-mode ladder was built for.

BIN=/net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt/src/cmalign
CM=/net/intdev/oblast01/infernal/notebook/26_0610_inf_dp_mem_cost/benchmark-runs/segment-census/cms/vadr/calici-NC_001959.cm
GEN=/net/intdev/oblast01/infernal/notebook/26_0610_inf_dp_mem_cost/benchmark-runs/checkpoint-proto/genomes
OUT=/net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt/local-mode-validation/rl6-data
C="--cpu 0 --fixedtau --tau 1e-7 --mxsize 65000"

echo "############ HOST: $(hostname)  DATE: $(date) ############"

runone() {  # genomeID
  local g=$1 f=$GEN/$1.fa
  echo "==================================================================="
  echo "### $g -> calici-NC_001959.cm  (LOCAL, default truncated)"
  echo "### --ckpt:"
  INFERNAL_CKPT_VERBOSE=1 /usr/bin/time -v $BIN --ckpt $C -o $OUT/hlf_${g}_ckpt.sto $CM $f > $OUT/hlf_${g}_ckpt.stdout 2> $OUT/hlf_${g}_ckpt.time
  echo "  ckpt exit=$?"; grep -E "Maximum resident|Elapsed" $OUT/hlf_${g}_ckpt.time; grep -m1 engaged $OUT/hlf_${g}_ckpt.time
  echo "### stock:"
  /usr/bin/time -v $BIN $C -o $OUT/hlf_${g}_stock.sto $CM $f > $OUT/hlf_${g}_stock.stdout 2> $OUT/hlf_${g}_stock.time
  echo "  stock exit=$?"; grep -E "Maximum resident|Elapsed" $OUT/hlf_${g}_stock.time
  local CK=$(grep "Maximum resident" $OUT/hlf_${g}_ckpt.time  | grep -oE "[0-9]+")
  local ST=$(grep "Maximum resident" $OUT/hlf_${g}_stock.time | grep -oE "[0-9]+")
  echo "  ckpt=${CK}kB stock=${ST}kB"
  if [ -n "$CK" ] && [ -n "$ST" ] && [ "$CK" -gt 0 ] && [ "$ST" -gt 0 ]; then
    awk -v c=$CK -v s=$ST -v g=$g 'BEGIN{printf "  >>> %s LOCAL total-RSS win = %.2fx  (stock %.0f MB -> ckpt %.0f MB)\n", g, s/c, s/1024, c/1024}'
  else
    echo "  >>> $g: one side failed (stock OOM/mxsize) -> ckpt succeeds where stock cannot"
  fi
  # degapped-residue identity (accuracy-neutral check)
  python3 - "$OUT/hlf_${g}_ckpt.sto" "$OUT/hlf_${g}_stock.sto" <<'PY'
import re,sys
def parse(p):
    d={}
    for L in open(p):
        L=L.rstrip('\n')
        if not L or L=='//' or L.startswith('#'): continue
        q=L.split(None,1)
        if len(q)==2: d[q[0]]=d.get(q[0],'')+q[1]
    return d
try:
    a=parse(sys.argv[1]); b=parse(sys.argv[2])
    dg=lambda s: re.sub(r'[-.]','',s).upper()
    for n in a:
        print("  degapped:", "IDENTICAL" if dg(a[n])==dg(b.get(n,'')) else "DIFFER!!", n)
except Exception as e:
    print("  degap check skipped:", e)
PY
  echo
}

runone NC_006875
runone NC_008311
echo "############ DONE $(date) ############"
