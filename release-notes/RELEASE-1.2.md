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

Infernal 1.2 is not yet released; this notes file will be finalized at
release time.
