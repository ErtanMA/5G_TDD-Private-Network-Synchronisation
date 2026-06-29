# Implementation Notes

This file tracks the working code plan and phase decisions while the thesis tool
is implemented.

## Current Constraint

Known gNB/operator/station locations will not be available. Therefore the thesis
must not include a "known gNB location" mode.

The supported workflow is:

1. Run the single-capture receiver chain on each B210 IQ capture.
2. Store one timing row per detected PCI per capture.
3. Record the receiver position for each capture point.
4. Subtract a reference PCI inside each capture to remove the B210 internal
   clock's arbitrary capture phase.
5. Fit unknown gNB positions and unknown relative transmit timing offsets.
6. Median-center the fitted transmit offsets and classify them with uncertainty.
7. Optionally replace the default static verdict thresholds with a project or
   operator-defined `SyncThresholdNs`. RMS and peak-to-peak timing-instability
   limits remain separate because those metrics are not physically equivalent
   to static offset.

## Implemented Public Entry Points

- `analyzeCapture`: analyze one IQ capture and produce detected-cell timing evidence.
- `analyzeCaptureAcrossSSBRaster`: search all complete legal n78 SSB centres in one wideband IQ capture and retain every confirmed `GSCN + PCI` cell.
- `enumerateN78SSBRaster`: enumerate legal 1.44 MHz-spaced n78 GSCN centres that fit inside the capture.
- `timingMeasurementsFromCaptureResult`: convert one capture result into survey rows.
- `noLocationSurveyCheck`: estimate relative synchronization from several capture positions without gNB locations.
- `generateSurveyReport`: write survey reports.
- `runSyntheticStressTest`: stress-test the survey-only workflow.
- `detectCellsWithSIC`: iterative near-far/overlap detection using PSS+SSS cancellation.
- `runMultiCaptureSurvey`: load real IQ/metadata files from a manifest and run the full survey workflow.
- `estimateRssPriorsAndPlanStops`: optional RSS-assisted rough gNB prior estimation and receiver-stop planning.
- `generateThesisArtifacts`: synthetic thesis figure/report package generator.
- `runThesisSyntheticResults`: deterministic thesis Results-section synthetic survey and receiver overview package.
- `run_identifiability_diagnostics`: processed-data scaled-Jacobian and multi-start offset-consensus analysis.
- `run_jackknife_from_processed`: leave-one-contributing-position-out sensitivity analysis without reprocessing IQ.

## Removed as a Supported Path

Known-location propagation correction has been removed from the public workflow.
The codebase should not ask for, require, or document private gNB coordinates.

## Phase 9 Status

Implemented modules:

- `analyzeCapture`: single-capture receiver chain.
- `detectCellsWithSIC`: near-far/overlap detection with local cancellation helpers.
- `estimatePSSCandidates`, `detectSSSAndPCI`, `estimateCellTiming`, `checkTDDPattern`: receiver-stage modules with local defaults.
- `noLocationSurveyCheck` and `estimateSurveyTimingNoLocations`: no-location survey fitting.
- `generateSyntheticScenario` and `buildSyntheticSurveyDataset`: synthetic validation data.
- `runMultiCaptureSurvey`: real-data manifest orchestration.
- `estimateRssPriorsAndPlanStops`: Phase 13 SS-RSRP-like power, RSS rough priors, and adaptive stop ranking.
- `generateThesisArtifacts`: Phase 14 final wording/report packaging and synthetic figures.
- Phase 9-14 regression tests.
- Phase 16 configurable verdict-threshold tests.
- Phase 17 thesis synthetic result package tests.
- Phase 18 multi-GSCN one-capture raster and timing-preservation tests.

Validated synthetic cases:

- aligned 5-gNB survey
- one 3 us timing fault
- one 250 ns warning-scale timing offset
- all gNBs shifted together by 3 us
- insufficient receiver positions
- receiver position input as latitude/longitude

Important limitation:

This phase removes the need for private gNB locations, but it does not remove
the need for known receiver locations. If neither gNB nor receiver positions are
known, static timing offset and propagation delay remain inseparable.

## Real Campaign Outcome

The completed P01-P15 campaign processed 45 accepted captures from 15
positions. Only 10 positions contributed accepted reference-target geometry
rows. With reference PCI 446, target support was 5 positions for PCI 559, 5 for
PCI 9, and 4 for PCI 489. The fit used 14 relative measurements for 11
parameters.

The scaled weighted Jacobian was full rank, but multi-start and
leave-one-position-out diagnostics showed that the static offsets were not
stable at the 1 us target. The correct result is therefore:

- receiver-chain and stationary-repeatability evidence: usable;
- public-campaign static no-location offsets: inconclusive;
- no public cell is declared synchronized or desynchronized from this fit.

## Next Priorities

1. If one short private-network test is possible, collect repeated stationary
   captures of a known nominal network. This strengthens false-alarm and
   stability evidence, but does not create an absolute timing reference.
2. For a future static-offset campaign, process each stop before moving and keep
   only routes that preserve common cells while improving scaled-Jacobian
   conditioning.
3. Consider a graph-based changing-reference estimator so positions connected
   through different overlapping PCIs are not discarded. This is useful only
   when the overlap graph remains connected.
4. Treat stronger SIC/PBCH reconstruction as lower priority. Better cell
   recovery can increase overlap, but it cannot repair weak survey geometry by
   itself.

## Ten-Issue Priority State

1. Absolute UTC/common-mode timing: not solvable with the current internal-clock
   B210 unless external timing hardware is added.
2. No gNB locations: implemented with multi-position survey fitting. Receiver
   positions are still required.
3. Overlapping SSB: started with iterative PSS+SSS+PBCH-DMRS cancellation.
4. Near-far masking: started with the same SIC implementation and synthetic
   same-NID2 near-far test.
5. B210 calibration: still to document and test when the physical device is
   available.
6. PBCH/MIB anchoring: PBCH DM-RS is implemented for cancellation; PBCH payload/MIB decode remains later.
7. More realistic RF impairments: implemented in synthetic generator and Phase 11 tests.
8. Real-data orchestration: implemented with `runMultiCaptureSurvey`.
9. Final threshold/report wording: implemented through
   `docs/final_verdict_wording.md` and Markdown survey summaries.
10. Thesis result packaging and figures: implemented through
   `generateThesisArtifacts` and `main_generate_thesis_artifacts.m`.
