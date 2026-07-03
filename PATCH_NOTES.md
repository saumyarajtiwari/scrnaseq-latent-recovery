# Patch Notes

## SymSim — BiocGenerics S4 rank() conflict

**Affected file:** `vendor/SymSim_patched/R/simulation_functions.R`
**Upstream repo:** https://github.com/YosefLab/SymSim

### Problem

When BiocGenerics is loaded (via scran/scater), it registers an S4 generic
for `rank()`. SymSim's `Get_params()` calls `alply()` (from plyr) followed by
`do.call(c, temp)`, which produces an S4 list-like object. BiocGenerics' S4
method dispatch then validates this object before calling the rank method and
fails with:

    "error in evaluating the argument 'x' in selecting a method for function
    'rank': 'x' must be an array of at least two dimensions"

This error occurs on every `SimulateTrueCounts()` call regardless of
parameters because `Get_params()` is called unconditionally.

### Fix

Two changes to `R/simulation_functions.R`:

1. **`Get_params()` (~line 96):** Replace `alply(X, 1, ...)` + `do.call(c, temp)`
   with `as.numeric(t(as.matrix(X)))` — mathematically identical (row-major
   unrolling of matrix X) but produces a plain numeric vector that bypasses S4
   dispatch. Changed `rank(values)` to `base::rank(values)`.

2. **`SimulateTrueCounts()` (~line 666):** Changed
   `rank(rowSums(params[[3]][chosen_hge,]))` to
   `base::rank(as.numeric(rowSums(as.matrix(params[[3]][chosen_hge, , drop=FALSE]))))`.

### Reproduction

```r
Rscript R/01_simulation/install_symsim_patched.R
```

This script clones SymSim, applies both patches, and installs from source.
The `vendor/SymSim_patched/` directory is gitignored (231MB) but regenerated
by this script.
