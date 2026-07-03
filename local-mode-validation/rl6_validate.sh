#!/bin/bash
# Brief 068 (R-L.6) capstone validation: cmalign --ckpt in LOCAL (default) mode,
# actual binary (not a driver), all four bps-class x truncation-mode combos.
# Compares --ckpt vs stock (no --ckpt), same config.
set -u
cd "$(dirname "$0")"
BIN=../src/cmalign
BPS0=rl2-data/matl300.cm         # pure MATL chain, bps=0
STR=../testsuite/rnaseP-eubact.c.cm   # structured, bps>0
BPS0SEQ=rl2-data/el300.fa        # MATL seqs (incl EL inserts)
STRSEQ=rl4-data/rnasep_full.fa   # rnaseP seqs
OUT=rl6-data
mkdir -p $OUT
COMMON="--fixedtau --tau 1e-7 --mxsize 16384"

run() {   # tag  cm  seq  extra_opts
  local tag=$1 cm=$2 seq=$3 extra=$4
  echo "======================================================================"
  echo "### $tag :  $cm  x  $seq   [opts: $extra $COMMON]  (LOCAL, default config)"
  INFERNAL_CKPT_VERBOSE=1 $BIN $extra --ckpt $COMMON -o $OUT/${tag}_ck.sto $cm $seq 2> $OUT/${tag}_ck.verbose
  local rc_ck=$?
  $BIN $extra          $COMMON -o $OUT/${tag}_st.sto $cm $seq 2> $OUT/${tag}_st.verbose
  local rc_st=$?
  echo "  exit codes: ckpt=$rc_ck stock=$rc_st"
  echo "  --- engine engaged (verbose):"
  grep -iE "engaged|checkpt" $OUT/${tag}_ck.verbose | sed 's/^/      /'
  if diff -q $OUT/${tag}_ck.sto $OUT/${tag}_st.sto >/dev/null 2>&1; then
    echo "  RESULT: BYTE-IDENTICAL to stock"
  else
    echo "  RESULT: DIFFERS from stock (expected for structured/pin cases) — quantifying:"
    # count alignment lines that differ
    diff $OUT/${tag}_ck.sto $OUT/${tag}_st.sto | grep -c '^[<>]' | sed 's/^/      diff-lines: /'
  fi
}

run bps0_notrunc  $BPS0 $BPS0SEQ "--notrunc"
run bps0_trunc    $BPS0 $BPS0SEQ ""
run struct_notrunc $STR $STRSEQ  "--notrunc"
run struct_trunc   $STR $STRSEQ  ""

echo "======================================================================"
echo "DONE"
