#!/bin/bash
# Compare structured --ckpt-vs-stock divergence in LOCAL vs GLOBAL, and verify
# degapped-residue identity (the accuracy-neutral test: same residues, possibly
# different column assignment / near-tie flips; NOT sequence corruption).
set -u
cd "$(dirname "$0")"
BIN=../src/cmalign
STR=../testsuite/rnaseP-eubact.c.cm
STRSEQ=rl4-data/rnasep_full.fa
OUT=rl6-data
COMMON="--fixedtau --tau 1e-7 --mxsize 16384"

# esl-alimanip / esl-reformat to pull degapped seqs; fall back to awk if absent
degap () {   # stofile -> fasta of degapped aligned seqs (sorted by name)
  ../easel/miniapps/esl-reformat -u fasta "$1" 2>/dev/null | \
    awk '/^>/{name=$0; next}{seq[name]=seq[name]$0} END{for(n in seq){gsub(/[-.]/,"",seq[n]); print n"\t"toupper(seq[n])}}' | sort
}

for cfg in "local:" "global:-g"; do
  name=${cfg%%:*}; gflag=${cfg##*:}
  for mode in "notrunc:--notrunc" "trunc:"; do
    mname=${mode%%:*}; mflag=${mode##*:}
    tag=${name}_${mname}
    $BIN $gflag $mflag --ckpt $COMMON -o $OUT/cmp_${tag}_ck.sto $STR $STRSEQ >/dev/null 2>&1
    $BIN $gflag $mflag        $COMMON -o $OUT/cmp_${tag}_st.sto $STR $STRSEQ >/dev/null 2>&1
    dl=$(diff $OUT/cmp_${tag}_ck.sto $OUT/cmp_${tag}_st.sto | grep -c '^[<>]')
    # degapped-residue identity
    degap $OUT/cmp_${tag}_ck.sto > $OUT/cmp_${tag}_ck.degap
    degap $OUT/cmp_${tag}_st.sto > $OUT/cmp_${tag}_st.degap
    if diff -q $OUT/cmp_${tag}_ck.degap $OUT/cmp_${tag}_st.degap >/dev/null 2>&1; then
      degid="DEGAPPED-IDENTICAL"
    else
      degid="DEGAPPED-DIFFERS(!!)"
    fi
    printf "%-16s diff-lines=%-4s  %s\n" "$tag" "$dl" "$degid"
  done
done
