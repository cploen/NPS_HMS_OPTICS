# Autonomous Clustering Optics Testbed

## Current branch

`testbed/autonomous-clustering-optics-20260604`

## Scope

This branch is for autonomous/semi-autonomous ytar, FP/PFP clustering, sieve row/column assignment, and event quality-score development for HMS/NPS optics.

## Purpose of this cleanup

This branch organizes exploratory optics scripts, diagnostics, and notes without deleting useful work. The goal is to make the active workflow easier to find while preserving older one-off tests in an attic directory.

## Organized script locations

### FP/PFP cluster and band-selection scripts

```text
scripts/optics/fp_pfp/
```

This directory contains the current exploratory scripts for focal-plane / focal-plane-prime selection, including:

- DBSCAN tests
- angle-scan band selection
- strong/weak band selection
- ridge-band and v-band methods
- PC1 / local-PC1 estimates
- early event-quality scoring tests

### Ytar / foil-selection notes

```text
scripts/optics/ytar/
```

This directory contains notes and summaries related to ytar foil/ridge selection, including multifoil ytar cut work.

Current guiding rule:

- Initial foil/ridge selection should use the full useful delta range, roughly delta = ±10%.
- Do not loop through delta slices for the initial foil-selection step.
- Delta-slice-dependent logic belongs later, for sieve-hole ownership and quality scoring.

### Experimental one-off tests

```text
attic/optics_exploration/macros/
```

This directory contains preserved one-off tests, debugging macros, and exploratory scripts that are not currently part of the active workflow.

These are kept for reference, not treated as active production scripts.

## Generated files

The following are treated as generated outputs and should normally remain untracked:

- PDF diagnostics
- ROOT files
- plots
- histogram outputs
- temporary scan outputs
- CINT / ACLiC build products

The repository `.gitignore` already ignores common generated outputs such as:

- `plots`
- `cuts`
- `hist`
- `ROOTfiles`
- `*.pdf`
- `*.png`
- `*.root`
- `*.pcm`
- `*_C.so`
- `*_C.d`

Force-add generated files only when they are intentionally being preserved as reference artifacts.

## Current cleanup commits

Recent cleanup commits include:

- pre-cleanup file inventory
- organization of exploratory optics scripts
- generated output inventory
- this workflow document

## Next real workflow task

Pick one active FP/PFP clustering method and promote it from exploratory status to current workflow.

Candidate methods currently preserved include:

- `assign_yfp_ypfp_angleScanBands_strongWeak_multifoil_ready.C`
- `assign_yfp_ypfp_angleScanBands_strongWeak_with_pages.C`
- `assign_yfp_ypfp_angleScanBands_interPeakOcc_with_pages.C`
- `assign_yfp_ypfp_angleScanBands_qpStripe_with_pages.C`
- `assign_xfp_xpfp_angleScanBands_strongWeak_with_pages.C`
- `dbscan_fp_pfp_islands.C`
- `event-quality-yfp_ypfp_onefoil.C`

The next decision should be based on which script most cleanly supports:

1. multiple foils,
2. delta slices,
3. ytar-cut input,
4. clear diagnostic output,
5. eventual quality scoring before SVD fitting.

