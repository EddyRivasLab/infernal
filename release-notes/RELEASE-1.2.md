# Infernal 1.2 release notes (in progress)

### Infernal 1.2 is under development; this file collects notable changes
### as they land on the development branch, ahead of an eventual release.

## Notable changes from 1.1.5:

 * **Infernal 1.2 reads only Infernal 1.2 model files: convert each
   Infernal 1.1 CM file once, with `cmconvert --upgrade`, before using
   it.** Every program except `cmconvert` refuses a CM file written by
   Infernal 1.1 (format `INFERNAL1/a`) and says to convert it, as
   Infernal 1.1 did for Infernal 1.0 files. This includes the CM files
   of current Rfam releases.
   - `cmconvert --upgrade old.cm > new.cm` writes the current format,
     `INFERNAL1/e` (see the next item). For a model built with
     `cmbuild --prior <f>`, add `--prior <f>`, or `--noprior` to
     convert it without that prior file. Converting takes seconds for
     most models and about a minute for the largest; all of Rfam takes
     a few CPU-minutes.
   - Converting leaves the CM's parameters and its `cmcalibrate`
     E-value parameters unchanged. It replaces only the filter HMM
     stored in the file with the profile HMMs Infernal 1.2 derives
     from the CM, so `cmsearch` and `cmscan` results can differ from
     Infernal 1.1's results on the original file. For most models,
     `cmalign` alignments are unchanged by the conversion. For some
     models, notably unstructured models built with `cmbuild --noh3pri`
     such as VADR's influenza models, bit scores and posterior
     probabilities can change slightly, and sequences that match the
     model poorly can align differently. This is because `cmalign`
     computes its alignment bands from the filter HMM, which
     `cmconvert --upgrade` recomputes from the model.
   - Press `cmscan` databases again: a database pressed by Infernal
     1.1 is refused too. Convert its CM file, then run `cmpress -F` on
     the converted file.
   - To write a model that Infernal 1.1 can read, use
     `cmconvert --outfmt 1/a`. Plain `cmconvert` of an `INFERNAL1/a`
     file writes it back in `INFERNAL1/a`.

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

 * **HMM-only search changes in `cmsearch` and `cmscan`.** These
   apply to `--hmmonly` and to models with no basepairs, which are
   searched in HMM-only mode by default.
   - The default HMM-only filter P-value thresholds are looser:
     `--hmmF1` 0.05 (was 0.02), `--hmmF2` 3e-3 (was 1e-3) and
     `--hmmF3` 1e-4 (was 1e-5), tuned for the HMM that 1.2 derives
     from the CM. HMM-only searches can report more hits, and take
     somewhat longer, than in 1.1.
   - `-g` combined with `--hmmonly`, which 1.1 rejected, now runs a
     glocal HMM-only search. With `-g` alone, models with no basepairs
     are now searched with a glocal HMM; in 1.1 `-g` switched them to
     a glocal CM search, which needed `cmcalibrate` E-value parameters.
     Use `--nohmmonly` for the 1.1 behavior.
   - New option `--hmmwindows`: with `--hmmonly` and `--noali`, stop
     after the third filter stage (Forward) and report the target
     windows that survive it instead of hits.

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

 * **New CM file format `INFERNAL1/e`, written by `cmbuild`.** CM
   files no longer store a separately built filter HMM after each CM.
   The filter HMM is now the CM's maximum likelihood p7 HMM; it and
   the HMM used for HMM-only search (`--hmmonly`) are derived from the
   CM whenever the file is read, and the file stores only their
   E-value parameters. The format also adds the null3-off E-value
   parameters above and consensus pseudoknot annotation (`PKNOT`).
   Infernal 1.1 cannot read `INFERNAL1/e` files: it stops with a file
   format error rather than misreading them.

   Only `cmconvert` reads Infernal 1.1 (`INFERNAL1/a`) CM files (see
   the previous item); `cmconvert --upgrade` converts them to
   `INFERNAL1/e`, deriving and calibrating the HMMs as `cmbuild`
   would. Infernal v1.0 CM files must still be converted with
   `cmconvert`, which now writes them in `INFERNAL1/e`. The `cmbuild` options
   `--p7ere`, `--p7eent`, `--p7prior`, `--p7hprior` and `--p7hemit`, which
   controlled how the separate filter HMM was built, have been
   removed; `--p7ml` is still accepted and has no effect.

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

 * **`cmalign` can align whole viral genomes with default settings.**
   In 1.1.5, aligning a SARS-CoV-2 genome (about 30 kb) to a model
   built from one needed about 28 Gb for the HMM band calculation
   alone, and was then refused at the default `--mxsize` after almost
   an hour; longer genomes crashed (github issue #74). In 1.2 a
   SARS-CoV-2 genome aligns to such a model in about 1 second using
   0.2 Gb (with the default HMM alignment for a model with no
   basepairs, below), or in 2 seconds using 0.4 Gb when aligned with
   the CM (`--nohmm`). The changes that make this possible:
   - For models with consensus length (CLEN) greater than 4000, the
     bands that constrain the alignment are derived from exact k-mer
     matches between each sequence and the model's consensus sequence,
     chained into a colinear set of anchors (`--p7band --p7kmerchain`).
     For the genome above, CM alignment with these bands was about 100
     times faster than with the plan 9 HMM bands used for smaller
     models, which `--hbanded` still selects.
   - `--mxsize` now also limits the memory of the plan 9 HMM band
     calculation, which is done in a checkpointed, lower-memory form
     when the full matrices do not fit. When the alignment matrices
     for a sequence do not fit, `cmalign` switches automatically, for
     that sequence only, to a checkpointed version of the optimal
     accuracy algorithm, and then to the divide-and-conquer CYK
     algorithm, before it stops with an error.
   - HMM alignment computes its Forward and Backward scores in double
     precision, and CM alignment rescales its scores during the
     calculation, so that single-precision rounding error no longer
     accumulates along very long sequences.

   The alignment method and bands used are reported in the `cmalign`
   output header.

 * **`cmalign` aligns models with no basepairs with their profile
   HMM** instead of the CM, which is faster and needs less memory. The
   new `--nohmm` option keeps such models on the CM. The new `--hmm`
   option selects HMM alignment for any model, and `--hmmvit` and
   `--hmmnoband` select its Viterbi alignment or its unbanded optimal
   accuracy alignment. HMM alignment scores are HMM bit scores, and
   are not comparable to CM bit scores.

 * **Models with no basepairs are now built the same way as models
   with basepairs**, with the same prior and entropy weighting target,
   instead of as if they were profile HMMs, so adding or removing a
   basepair no longer changes the rest of the model. Rebuilding a
   model that has no basepairs with 1.2 therefore changes its
   parameters and bit scores, and score thresholds set with an
   earlier version for such a model (such as Rfam-style gathering
   thresholds) may need to be re-checked. CM files built with earlier
   versions keep their parameters. `cmbuild --noh3pri`, which used to
   give this behavior, is still accepted and has no effect.

 * **New `cmalign` band options.** `--p7band` derives the bands for CM
   alignment from the model's profile HMM instead of its plan 9 HMM;
   `--p7ibv` (bands from Forward/Backward, with `--p7ibv-mem` for a
   low-memory version) and `--p7kmerchain` (k-mer chain bands) choose
   how. They also apply to HMM alignment (`--hmm`). See the `cmalign`
   man page.

 * **`cmalign --sub` and `cmbuild --refine --sub` have been removed.**
   `--sub` aligned each sequence to a sub-CM covering the consensus
   columns that an HMM predicted for its first and last residues. It
   was an older method for partial sequences than truncated alignment,
   which has been the default since 1.1 and handles partial sequences
   as accurately, much faster and in far less memory (`--sub` could
   need about 28 GB for one 30 kb viral genome). For the glocal
   equivalent of `-g --sub --notrunc`, use `-g` without `--notrunc`.
   Command lines that pass `--sub` now fail with an unknown-option
   error.

## Bug fixes:

 * fixes `cmalign` crashing (segmentation fault, no error message)
   when (sequence length + 1) x (model CLEN + 1) exceeds 2^31, for
   example when aligning a genome to a genome-length model of more
   than about 46 kb (github issue #74).

 * fixes truncated alignment silently failing for models with a CLEN
   of 46,341 or more, because of an integer overflow (github issue
   #59).

 * fixes `cmalign`'s HMM-banded Outside algorithm writing local-end
   (EL) matrix cells outside their valid range, which changed some
   local alignments (github issue #64).

 * fixes `cmalign --nonbanded` crashing on a NaN score in truncated
   alignment (github issue #52).

 * fixes `cmalign --mapali` corrupting rows of the included alignment
   that begin with `~` (github issue #61).

 * fixes `cmbuild -O` intermittently writing one row of the resaved
   alignment truncated and filled with stray bytes (github issue #66).

 * fixes `cmbuild --eset 0` writing models that could not be read,
   for alignments with no basepairs and for `--p7ml` builds (github
   issue #57).

 * fixes `cmbuild` aborting with "xl,xr do not bracket a root" during
   entropy weighting for some alignments (github issue #71).

 * fixes the `bias` column of HMM-only `cmsearch` and `cmscan` hits
   being reported in nats instead of bits (github issue #55).

 * fixes `cmsearch` and `cmscan` `--incT` working only when `-T` was
   also used (github issue #40).

 * fixes a `cmbuild --refine` bug related to fragments (github issue
   #42).

 * fixes `cmalign --mpi` workers quitting early when they were given
   no sequence, for example when there were fewer sequences than
   workers. Open MPI then printed "unexpected tag-receive descriptor
   ... not matched" warnings to standard output, and in rare cases the
   run could hang. Alignments were not affected.

 * consensus structure lines now annotate single-stranded positions
   between top-level stems as `:` (full WUSS notation) instead of `-`
   or `,` (github issue #69).

Infernal 1.2 is not yet released; this notes file will be finalized at
release time.
