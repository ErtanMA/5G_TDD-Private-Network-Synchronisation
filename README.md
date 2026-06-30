# Passive Relative Synchronisation Fault Detection in 5G TDD Networks

MATLAB implementation of a passive, IQ-based method for detecting **relative**
synchronisation faults between visible 5G NR TDD base stations (gNBs) in a
private network, using USRP B210 IQ captures.

The receiver is driven by its internal clock only (no GPSDO / GNSS reference),
so the method checks whether visible cells are aligned **relative to one
another** within the measured survey geometry. It does **not** claim absolute
UTC or regulatory timing compliance.

## What it does

1. **Single-capture receiver chain** — blind cell detection from raw IQ:
   SSB centre-frequency search over the n78 SS raster, handmade PSS
   matched-filter search, CFO estimation/correction, SSS/PCI recovery, and
   per-PCI frame-arrival timing.
2. **Multi-position survey** — with no known gNB locations, jointly fits the
   unknown effective cell positions and relative transmit offsets to the
   reference-subtracted arrivals (weighted nonlinear least squares), separating
   transmit-timing offset from propagation delay.
3. **Verdict** — judges each cell on its static relative offset and its relative
   timing instability against a configurable threshold, and reports
   `PASS` / `SUSPECT` / `FAIL` / `NOT_ASSESSABLE`. Survey conditioning (Jacobian
   condition number and offset uncertainty) decides whether a verdict can be
   trusted at all.

## Repository layout

```
src/        Core implementation (receiver chain, survey estimator, verdict logic)
tests/      Phase-by-phase regression tests of the pipeline
scripts/    Standalone diagnostic runners
examples/   Capture-manifest and metadata templates
docs/        Implementation notes
```

## Requirements

- MATLAB R2023a or newer
- 5G Toolbox (`nrPSS`, `nrSSS`, `nrPBCHDMRS`, `nrOFDMModulate`/`nrOFDMDemodulate`)
- Signal Processing Toolbox

## Quick start

```matlab
addpath('src');

% Synthetic validation (no hardware needed): generate a scenario, run the
% receiver chain, and inspect the recovered PCI and timing.
result = runThesisSyntheticResults();

% Process a real multi-position survey from saved captures described by a
% capture manifest (see examples/capture_manifest_template.csv).
survey = noLocationSurveyCheck(measurementTable, captureInfo);
```

Run the test suite from the repository root:

```matlab
addpath('src'); addpath('tests');
results = runtests('tests');
```
