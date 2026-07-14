# Infernal 1.2 release notes (in progress)

### Infernal 1.2 is under development; this file collects notable changes
### as they land on the development branch, ahead of an eventual release.

## Notable changes from 1.1.5:

 * The `cmsearch`/`cmscan` per-hit alignment display's base-pair-quality
   trailer line, previously labeled `NC`, is now labeled `PS` ("Pair
   Status") on every structured-CM hit, not just pseudoknot-containing
   ones. The line's existing negative-only nested-pair markup (`v`
   broken, `?` truncated half) is unchanged; pseudoknot pairs additionally
   get positive+negative pair-status marks (`=` maintained, `$` covarying,
   `x` broken). **Any user script that greps the human-readable alignment
   block for the ` NC` trailer will silently stop matching** — this is a
   plain-text/alidisplay-only change; machine-readable formats (`--tblout`,
   Stockholm output from `cmalign`) are unaffected.

 * **`cmbuild` now calibrates E-value statistical parameters
   automatically**, using a fast analytical predictor
   (`cm_FastCalibrate`) in place of `cmcalibrate`'s Monte Carlo
   simulation. A default-construction CM out of `cmbuild` is
   immediately ready for `cmsearch`/`cmscan` — `cmcalibrate` is no
   longer a required step in the standard build-then-search workflow.
   `cmcalibrate` remains available as a fallback: it is still required
   for CMs built with non-default `--pbegin`/`--pebegin`/`--pend`/
   `--pfend`/`--null`, for anyone who wants true simulated E-value
   statistics rather than the analytical prediction, and for accurate
   `--nonull3` E-values on large or glocal-mode models (see the known
   limitations below). Models built with `--ere`/`--enone` entropy
   weighting are calibrated correctly by a dedicated entropy-aware
   variant of the predictor.

 * **CMs now store both a null3-on and a null3-off E-value parameter
   set ("store-both").** `cmbuild` writes both by default, so
   `cmsearch --nonull3`/`cmscan --nonull3` work directly off of a
   normally-built CM with no extra calibration step. A CM lacking the
   null3-off set (e.g. an older or externally-produced CM) falls back
   to the null3-on parameters for a `--nonull3` search, with a
   one-time warning.

 * **New CM file format `INFERNAL1/d`**, which adds the null3-off
   parameter block above and consensus pseudoknot annotation (`PKNOT`)
   to the model header and node lines. Older Infernal binaries cleanly
   reject an `INFERNAL1/d` binary CM file rather than misreading it;
   the ASCII format identifier and CM files without either addition
   are unaffected.

 * **Known limitations of the new analytically-predicted E-values:**
   - Predicted `--nonull3` E-values are accurate for local-mode search
     of most models, but can be over-confident (too low, by up to
     several orders of magnitude for large, AT-rich models such as the
     large ribosomal RNAs) for glocal-mode search and for larger models
     (consensus length greater than a few hundred). Run
     `cmcalibrate --nonull3` for accurate `--nonull3` E-values on such
     models.
   - Glocal-mode E-values for large and huge models (consensus length
     greater than ~1500) are less reliable than local-mode E-values for
     both the predictor and `cmcalibrate`, because the underlying score
     distributions are sparsely sampled at that scale.
   - Default (null3-on) local-Inside E-values for strongly AT-rich or
     low-complexity models carry a small residual over-confidence
     (roughly two- to three-fold near the significance threshold),
     within the intrinsic seed-to-seed variation of `cmcalibrate`
     itself and not expected to materially affect hit ranking.
   - Default (null3-on), local-mode E-values are unaffected by any of
     the above and were validated to be effectively identical to prior
     `cmcalibrate`-derived values, including on the `rmark4h` benchmark.

Infernal 1.2 is not yet released; this notes file will be finalized at
release time.
