# Passive Relative Synchronisation Fault Detection in 5G TDD Networks

MATLAB implementation of a passive, IQ-based method for detecting **relative**
synchronisation faults between visible 5G NR TDD base stations (gNBs) in a
private network, using USRP B210 IQ captures.

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
