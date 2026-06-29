# R-L.1 — D&C CYK local-mode validation (brief 059)

Validation-only (NO engine edits). Confirms the D&C CYK **pin source**
(`CYKDivideAndConquerHB` + its `*_hb` kernels) produces byte-exact parses in
**local** config (HMM-local + EL), the foundation for the local-mode `--ckpt` ladder.

## What this is
- `hbdnc_drv.c` — the existing isolated correctness driver, here extended with
  **EL local-end coverage** instrumentation (R-L.1 change):
  - `count_el(tr, cm->M)` — counts USED_EL nodes (state == `cm->M`, recorded by
    `InsertTraceNode` at `cm_dpsmall.c:7705/8497`); the unambiguous "EL taken" signal.
  - `entry_state(tr)` — reports the root-child entry state (local-begin diagnostic;
    see caveat below).
  - per-case `elA=/elB=/entryB=` columns + a panel `EL coverage:` summary line;
    `n_el_mismatch` is folded into the PASS/FAIL gate.
- The driver already defaulted to **local** config with a `-g` toggle
  (`hbdnc_drv.c` header + lines 87-93) — no "requires -g" assert to remove
  (unlike the ckpt drivers, `ckptoa_drv.c:516` / `ckpttr_drv.c:741`).

## Build (against this worktree's libinfernal)
```
R=/net/intdev/oblast01/infernal/git/EddyRivasLab/infernal-local-ckpt
gcc -O3 -DHAVE_CONFIG_H -I$R/src -I$R/src/impl_sse -I$R/easel -I$R/hmmer/src \
    -L$R/src -L$R/hmmer/src -L$R/easel \
    -o hbdnc_drv_local hbdnc_drv.c -linfernal -lhmmer -leasel -lpthread -lm
```

## Gate
Per case: `ParsetreeCompare`-equivalent (`trees_equal`: state/emitl/emitr byte-identical)
== oracle `cm_AlignHB` parse AND CYK score match (≤0.01; float-sum-order rounding only).
PLUS EL structure: `elA == elB` on non-clip cases, and ≥1 panel case must actually take EL.

## Result (rl1_*.log): R-L.1 GREEN
All panel CMs PASS byte-exact in BOTH local and global config; EL local-ends are
exercised in local (and absent in global), byte-identical between D&C(B) and stock(A):

| CM     | M     | local PASS | global PASS | EL nodes (local, A==B) |
|--------|-------|-----------|-------------|------------------------|
| tRNA   | 227   | yes       | yes         | 1 (all 8 seqs)         |
| 5S     | 369   | yes       | yes         | 1 (sample3); +clip tau=0.2 ok |
| RNaseP | 960   | yes       | yes         | 1-2                    |
| LSU    | 10559 | yes (71 bif) | yes      | 3 (bif × EL × local)   |
| dengue | 32284 | global ok; local genome-scale slow (see brief) | yes | 0 (bps=0) |

## Caveat
`entryB` reads 3 in BOTH modes — state 3 is the canonical ROOT child, not a
local-begin signal. Local-begin *equivalence* between A and B is nonetheless
guaranteed by `trees_equal` (byte-identity covers node 1). The panel's full-length
seqs enter normally; the load-bearing local feature exercised here is EL (local-end).
