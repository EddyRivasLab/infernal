# R-L.2 — bps=0 non-trunc LOCAL: EL (local-end) in cm_CheckptAlignHB (brief 061)

First checkpointed engine to get EL — the **template** for R-L.3/4/5. All new EL
code is gated on `CMH_LOCAL_END`; the global (`-g`) path is byte-identical to the
pre-R-L.2 engine.

## What was built (cm_dpalign.c, gated on CMH_LOCAL_END)
Per the brief-060 recipe:
- **Forward** (`ckpt_inside_deck`): on-the-fly `el_scA[d-sd]+endsc` base-case
  re-init. **No deck** (the fixed ramp; mirrors cm_InsideAlignHB:3576-3583).
- **Outside** (`ckpt_outside_deck`): v->EL accumulation into a **banded** EL deck
  `elbeta[j][0..eldmax[j]]` (upper-d-edge banded; cells above `eldmax[j]` are
  provably IMPOSSIBLE). Mirrors cm_OutsideAlignHB:6356-6402.
- **Posterior** (`cm_CheckptAlignHB` Step B tail): EL->EL self-transition
  (band-edge: top cell keeps its value, absent `d+1` neighbour reads IMPOSSIBLE)
  + fold to the 1-D `l_pp[cm->M]`; EL added to `sum[]` normalization LAST and the
  EL row normalized (mirror cm_EmitterPosteriorHB:6961-6999/6832-6843).
- **OptAcc** (`ckpt_optacc_deck`): `elalpha` prefix-sum built from the normalized
  `l_pp[cm->M]` (**SUBSTITUTE, don't omit** — the 51/053 segfault); per-state EL
  re-init `alpha[v]=elalpha[j-sdr][d-sd]`; d==0 EL routing via precomputed
  `el_esc`/`el_endsc` (the have_el branch of the DZero shadow init); the two
  `USED_EL->IMPOSSIBLE` resets gated on `!have_el` so local EL cells keep value.
- **Traceback**: `USED_EL -> cm->M` was already pre-wired; lit up by the OA forward.

**Banding is byte-exact**: `FLogsum(x, IMPOSSIBLE) == x` (IMPOSSIBLE = -1e36, not
-inf), so dropping out-of-band EL cells matches stock cell-for-cell.

## Driver: rl2_ckptalign_drv.c
Calls the **library** `cm_CheckptAlignHB()` directly (unlike the 023-era
checkpoint-proto/ckptoa_drv.c, which reimplemented the engine inline) and compares
vs stock `cm_AlignHB(do_optacc=TRUE)`. Tests the actual R-L.2 deliverable.

**Config note:** `cm_CheckptAlignHB` supports local **ENDS** (EL) but not local
**BEGINS** (its qualifier rejects `CMH_LOCAL_BEGIN`; local begins are a separate
rung). The default (non `-g`) driver config is therefore **local-ends-only**: it
un-does the local-begin half of `cm_localize()` (restores ROOT transitions, clears
the begin flag), keeping local ends + the EL-configured cp9 bands. EL is
independent of begins, so this fully exercises the new EL code. `-g` = pure global
(regression).

## Build
```
R=/net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt
gcc -O3 -DHAVE_CONFIG_H -I$R/src -I$R/src/impl_sse -I$R/easel -I$R/hmmer/src \
    -L$R/src -L$R/hmmer/src -L$R/easel \
    -o rl2_ckptalign_drv rl2_ckptalign_drv.c -linfernal -lhmmer -leasel -lpthread -lm
```

## Gate (per case)
Inside score `|dZ|` (byte-exact), `ParsetreeCompare==0`, per-residue PP string
`strcmp==0`, EL-node count (`state==cm->M`) match between ckpt and stock.

## Test data (rl2-data/)
- `matl60.{sto,cm}` — synthetic bps=0 MATL chain, clen 60, M 184.
- `matl300.{sto,cm}` — synthetic bps=0 MATL chain, clen 300, M 904.
- `eltest.fa` — 6 seqs for matl60 (prefix-match + divergent tail → favors EL).
- `el300.fa` — 2 seqs for matl300.
Genome-scale viral CMs (bps=0): `../../segment-census/cms/vadr/{calici,flavi}-*.cm`
with cmemit'd genomes (`../rl2-el/{calici,flavi}_emit.fa`).

## Result: R-L.2 GREEN
- **matl60** (6 seqs): 6/6 byte-exact PASS in EL-only local AND global; **3 cases
  take 1 EL** (max for bps=0: a local end terminates the single chain),
  byte-identical to stock. Robust across tau 1e-3…1e-10 (band-edge stress).
- **matl300** (M=904): byte-exact PASS, 1 EL fires; √M CM-DP win ~6× WITH EL on;
  banded EL deck = 300–1785 cells vs unbanded O(L²/2)=45150 → **25–151× smaller**.
- valgrind clean (no leaks / invalid reads) on EL-firing case.
- Genome-scale (calici/flavi) — see `rl2_viral.out`.

See logs `rl2_matl60_local.log`, `rl2_matl60_global.log`, `rl2_matl300_local.log`,
`rl2_tau_sweep.log`, `rl2_viral.out`.
