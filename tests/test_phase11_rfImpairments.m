% test_phase11_rfImpairments
% Phase 11 tests for more realistic RF impairment stress cases.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot, "src")));

fprintf("Running Phase 11 RF impairment tests...\n");

testModerateRFImpairmentDetection();
testImpairmentControlsChangeWaveform();
testHarshRFCaseIsCharacterized();

fprintf("Phase 11 RF impairment tests passed.\n");

function testModerateRFImpairmentDetection()
[iq, meta, truth] = generateSyntheticScenario("rf_impairment_moderate");
result = analyzeCapture(iq, meta, struct("SIC", struct("MaxIterations", 5)));

detectedPCIs = unique(result.SSSDetections.PCI(result.SSSDetections.IsUsable));
assert(numel(detectedPCIs) >= 3, ...
    "Moderate RF impairment case should still recover several PCIs.");
assert(any(ismember(truth.PCI, detectedPCIs)), ...
    "Expected at least one true PCI under moderate RF impairments.");

timing = result.CellTiming;
assert(height(timing) >= 3, ...
    "Moderate RF impairment case should produce multiple timing rows.");

faultPCI = truth.PCI(truth.FrameOffsetNs == 250);
if any(timing.PCI == faultPCI)
    row = timing(timing.PCI == faultPCI, :);
    assert(abs(row.RelativeArrivalOffsetNs) < 800, ...
        "250 ns impaired case should remain sub-microsecond at arrival-timing level.");
end
end

function testImpairmentControlsChangeWaveform()
baseOpts = struct("RandomSeed", 909, "NumGNBs", 1, "PCIs", 11, ...
    "FrameOffsetsNs", 0, "CFOHz", 0, "GNBPowerdB", 0, ...
    "LocationUncertaintyM", 30, "SiteDistanceM", 100, "SNRdB", 50);
impairedOpts = baseOpts;
impairedOpts.MultipathDelaysSamples = [0 3.25];
impairedOpts.MultipathGainsdB = [0 -6];
impairedOpts.MultipathPhasesRad = [0 1.1];
impairedOpts.IQGainImbalanceDB = 1;
impairedOpts.IQPhaseImbalanceDeg = 5;
impairedOpts.DCOffset = 0.03 + 0.02i;
impairedOpts.PhaseNoiseStepStdDeg = 0.03;
impairedOpts.ClippingLevelRMS = 3;

[iqBase, ~, ~] = generateSyntheticScenario("aligned_5gnb", baseOpts);
[iqImpaired, ~, ~] = generateSyntheticScenario("aligned_5gnb", impairedOpts);
delta = norm(iqBase - iqImpaired) / norm(iqBase);
assert(delta > 0.05, "RF impairment controls should materially change the waveform.");
end

function testHarshRFCaseIsCharacterized()
[iq, meta, ~] = generateSyntheticScenario("rf_impairment_harsh");
result = analyzeCapture(iq, meta, struct("SIC", struct("MaxIterations", 5)));

assert(isfield(result, "SICDebug"));
assert(isfield(result.SICDebug, "NumIterations"));
assert(result.SICDebug.NumIterations >= 1);
% The harsh case is not a guaranteed pass; the test ensures it runs and
% produces bounded diagnostic tables instead of crashing.
assert(istable(result.SSSDetections));
assert(istable(result.CellTiming));
end
