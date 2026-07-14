# brief 26_0610-090 — marginal (L/R) local-end EL in the truncated D&C: BUILT + validated

Builds the deferred marginal EL that brief 26_0610-088 root-caused (brief 26_0610-049
deferral): `tr_outside_hb` built the local-end (EL) outside deck only on the J plane
(`beta[cm->M]`); `betaL[cm->M]`/`betaR[cm->M]` were never built. Engine change is in
`src/cm_dpsmall.c` (commit c543f1ae). See the brief-090 summary for the full writeup.

## VERDICT
BUILT. 6 bps=0 EL-terminus cases fixed byte-exact; **0 regressions** across ~6200
structured/aggressive/bps=0 cases. sample9_5ptr is NOT fixed — it is a SEPARATE
pre-existing inside-begin inflation (own follow-on brief), not the marginal-EL gap.

## Reproduce
```
R=/net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-trcyk-dnc-localbug
cd $R/src && make libinfernal.a
cd $R/local-mode-validation
gcc -O3 -DHAVE_CONFIG_H -I$R/src -I$R/src/impl_sse -I$R/easel -I$R/hmmer/src \
    -L$R/src -L$R/hmmer/src -L$R/easel -o r090_forked_drv r090_forked_drv.c \
    -linfernal -lhmmer -leasel -lpthread -lm
# 2 known cases (sample5 now PASS byte-exact, sample9 still FAIL = separate bug):
./r090_forked_drv rl2-data/matl300.cm r088-data/two.fa
# full bps=0 panel (crash-tolerant; ~30 min):
./r090_forked_drv rl2-data/matl300.cm r088-data/matl300_panel.fa
```

## r090_forked_drv.c
Crash-tolerant byte-exact panel driver: same gate as r086_argmax_drv (argmax over
root-valid modes vs the monolithic oracle cm_TrAlignHB; PASS iff mode match AND
|dev|<=0.01 AND ParsetreeCompare==0), but FORKS per sequence so the pre-existing
stock segfault (067/068, short truncated-local seqs) on one seq does not abort the
panel. Prints PASS/FAIL/CRASH per seq. Built pre-fix (via `git stash push
src/cm_dpsmall.c` + rebuild) and post-fix to diff fail-sets.

## Validation logs (in this dir, uncommitted)
- `r090_{pre,post}fix_forked.log`   -- full bps=0 1800-panel (diff: 0 new, 6 fixed)
- `r090_struct_{pre,post}fix.log`   -- structured MP panel (diff: 48==48 identical)
- `r090_aggr_{pre,post}fix.log`     -- aggressive helix-crossing panel (113==113 identical)
- `r090_valgrind_{bps0,struct}.err` -- valgrind (0 new errors vs baseline)
- structured CMs + panels under `r090-struct/`

## Key measurement for the separate sample9 bug (inside-begin inflation)
oracle forced-L Ldp[606][180][180] = -15.175 (true begin@606 = -30.74, rejected);
oracle picks begin@594 (Ldp[594]=-13.437 -> -29.011). D&C forced-L returns -29.002
via a phantom begin@606, computing 606's inside as -13.436 (inflated 1.738 bits =
emission mass of states 594..603 above 606). Not fed by the outside-EL decks ->
unchanged by this brief's fix.
