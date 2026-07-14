# brief 26_0610-088 — L-mode undershoot in the truncated D&C, ROOT-CAUSED (not fixed)

Investigates the 2/1800 `matl300` (bps=0) L-mode undershoots brief 086 flagged
(`sample5_3ptr`, `sample9_5ptr`): the D&C forced-L result collapses to a 4-node
parse instead of the oracle's ~150-node one.

## VERDICT
Root cause **found and proven** (byte-exact); **not fixed** (a deferred feature,
own follow-on brief). Engine (`src/cm_dpsmall.c`) left pristine except a
doc-comment breadcrumb.

**Root cause:** the truncated D&C `tr_outside_hb()` builds the local-end (EL)
outside deck ONLY on the J plane (`beta[cm->M]`). The marginal L/R EL decks
(`betaL[cm->M]`/`betaR[cm->M]`) are **never built** — a gap the original 049
author documented (`cm_dpsmall.c:7166-7167`: "marginal local-end ... NOT yet
built; valid only for CMH_LOCAL_END==off (global)" and `:7594-7596`). A truncated
L(/R) parse whose optimum terminates via EL is therefore unrepresentable in
`tr_wedge_splitter_hb`/`tr_generic_splitter_hb`; they fall back to the J-plane EL
deck (or a mode-converting split) and undershoot.

**Proof (sample5_3ptr):** inside is byte-exact (`Ldp[471]=Ldp[3]=9.403411`, ==
the D&C's begin-at-471 sub-wedge score). The oracle picks begin@3
(`trpenalty(3)+9.403 = -6.184`); the D&C's top-level inside-begin only sees the
lower model half, so it offers begin@471 (`trpenalty(471)+9.403 = -6.229`), and
the begin@3 parse must come via the std-split `betaL` candidate, which maxes at
only -6.237 (0.05 short) because the parse ends at state 420->EL entirely ABOVE
the midnode split (state 450) and needs the missing L-EL candidate. The
hypothetical L-EL feed `betaL[420][j][d+1] + endsc[420] + el_selfsc*d + esc`
= **-6.18375**, byte-exact vs the oracle **-6.183784**; the J-plane EL cand the
code actually uses maxes at **-43.68**. (sample9_5ptr is the same family via
nested recursion: its EL is deep, reached after a mode-converting L->J split that
then falls to the J-plane EL.)

## Reproduce
```
R=/net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-trcyk-dnc-localbug
# 1. regenerate the exact 086 panel (deterministic):
$R/src/cmemit --seed 86 -N 600 rl2-data/matl300.cm > r088-data/matl300_full.fa
python3 -c "
seqs=[];name=None;buf=[]
for line in open('r088-data/matl300_full.fa'):
    if line.startswith('>'):
        if name: seqs.append((name,''.join(buf)))
        name=line[1:].split()[0];buf=[]
    else: buf.append(line.strip())
if name: seqs.append((name,''.join(buf)))
o=open('r088-data/matl300_panel.fa','w')
for nm,s in seqs:
    L=len(s); o.write('>%s_full\n%s\n'%(nm,s))
    o.write('>%s_5ptr\n%s\n'%(nm,s[int(L*0.45):])); o.write('>%s_3ptr\n%s\n'%(nm,s[:int(L*0.55)]))
"
# 2. build + run the argmax driver (the two failures):
gcc -O3 -DHAVE_CONFIG_H -I$R/src -I$R/src/impl_sse -I$R/easel -I$R/hmmer/src \
    -L$R/src -L$R/hmmer/src -L$R/easel -o r086_argmax_drv r086_argmax_drv.c \
    -linfernal -lhmmer -leasel -lpthread -lm
./r086_argmax_drv rl2-data/matl300.cm r088-data/two.fa   # 2 FAILs, 4-node parses
```
`r088-data/two.fa` (the two failing fragments) is committed for convenience.

## Investigation tooling (committed)
- `r088_probe_drv.c` — forces a mode (J/L/R), dumps BOTH the D&C
  `TrCYKDivideAndConquerHB` parse and the oracle `cm_TrAlignHB` parse
  (+ `ParsetreeScore`) for eyeball comparison of the structural collapse.
- `r088_oracleL_dump.c` — dumps the oracle forced-L inside-HB `Ldp[v][L][L]` at
  chosen states (used to show the inside is byte-exact: `Ldp[471]=9.403411`).

## Fix spec (for the follow-on brief)
1. `tr_outside_hb`: allocate + init + seed + propagate `betaL[cm->M]` (and
   `betaR[cm->M]`) mirroring the J EL feed (`cm_dpsmall.c:7253-7294` seed,
   `:7518-7570` propagate) but with per-mode marginal emission rules (left-only
   for L, right-only for R; cf. `sdl`/`sdr` and `Lelbeta`/`Relbeta` in
   `cm_dpalign_trunc.c`). The ML/IL L-feed formula is already proven correct
   (== the J feed with `beta`->`betaL`).
2. `tr_wedge_splitter_hb` + `tr_generic_splitter_hb`: add L/R EL candidates
   reading `betaL[cm->M]`/`betaR[cm->M]`, setting `p_mode=L/R`, `c_mode=T`.
3. Traceback: route the `best_v==-1` EL case through the mode-aware
   `tr_v_splitter_hb` (pass `z_allow_L`/`z_allow_R`), which delegates to
   `tr_vinsideT_hb` (already EL + mode aware).
4. Fix the deck-M free logic (`:7595-7604`; the note there says deck M is not
   allocated for L/R — it now would be).
5. **Validate byte-exact on BOTH bps=0 (this 1800-panel: 2->0 fails, zero new)
   AND structured J/L/R (tRNA/5S/LSU, the 079/086 drivers)** — the bps=0 panel
   does NOT exercise MP states, so the per-mode MP EL rule needs the structured
   panel to catch regressions.
