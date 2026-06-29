% test_phase16_configurableThresholds
% Focused tests for user-configurable synchronization decision thresholds.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot, "src")));

fprintf("Running Phase 16 configurable threshold tests...\n");

testStaticSyncThresholdCanBeRelaxed();
testStaticSyncThresholdCanBeTightened();
testNearThresholdBecomesSuspect();
testSeparateJitterThresholdAliases();

fprintf("Phase 16 configurable threshold tests passed.\n");

function testStaticSyncThresholdCanBeRelaxed()
[meas, caps, truth] = offsetSurveyLocal(3000, 1601, 5);
faultPCI = truth.GNB.PCI(4);

opts = baseSurveyOptsLocal();
opts.SyncThresholdNs = 5000;

result = noLocationSurveyCheck(meas, caps, opts);
row = result.TimingEstimates(result.TimingEstimates.PCI == faultPCI, :);

assert(row.StaticTimingStatus == "PASS", ...
    "A 3 us relative offset should pass when the configured static limit is 5 us.");
assert(isinf(result.FitInfo.HardFailOffsetNs), ...
    "Supplying SyncThresholdNs should disable the legacy hard-fail threshold unless explicitly configured.");
end

function testStaticSyncThresholdCanBeTightened()
[meas, caps, truth] = offsetSurveyLocal(3000, 1602, 5);
faultPCI = truth.GNB.PCI(4);

opts = baseSurveyOptsLocal();
opts.SyncThresholdNs = 2000;

result = noLocationSurveyCheck(meas, caps, opts);
row = result.TimingEstimates(result.TimingEstimates.PCI == faultPCI, :);

assert(row.StaticTimingStatus == "FAIL", ...
    "A 3 us relative offset should fail when the configured static limit is 2 us.");
end

function testNearThresholdBecomesSuspect()
[meas, caps, truth] = offsetSurveyLocal(990, 1603, 15);
faultPCI = truth.GNB.PCI(4);

opts = baseSurveyOptsLocal();
opts.SyncThresholdNs = 1000;
opts.MaxFitRMSNsForPass = 250;

result = noLocationSurveyCheck(meas, caps, opts);
row = result.TimingEstimates(result.TimingEstimates.PCI == faultPCI, :);

assert(row.StaticTimingStatus == "SUSPECT", ...
    "An offset whose uncertainty interval touches the configured limit should be SUSPECT.");
end

function testSeparateJitterThresholdAliases()
[meas, caps, truth] = buildSyntheticSurveyDataset("aligned_5gnb", struct( ...
    "MeasurementNoiseNs", 5, "RandomSeed", 1604, "NumCaptures", 10));
faultPCI = truth.GNB.PCI(4);
jitterNs = [650; -520; 430; -610; 360; -560; 720; -380; 540; -640];
meas = addPerCaptureJitterLocal(meas, faultPCI, jitterNs);

opts = baseSurveyOptsLocal();
opts.JitterRmsThresholdNs = 1000;
opts.JitterPeakToPeakThresholdNs = 2500;

result = noLocationSurveyCheck(meas, caps, opts);
row = result.TimingEstimates(result.TimingEstimates.PCI == faultPCI, :);

assert(row.TimingInstabilityStatus ~= "FAIL", ...
    "Large relative jitter should not fail when explicit jitter limits are relaxed.");
assert(result.FitInfo.FailInstabilityRmsNs == 1000);
assert(result.FitInfo.FailInstabilityPeakToPeakNs == 2500);
end

function [meas, caps, truth] = offsetSurveyLocal(offsetNs, seed, noiseNs)
tx = [0 0 0 offsetNs 0];
[meas, caps, truth] = buildSyntheticSurveyDataset("aligned_5gnb", struct( ...
    "TxOffsetsNs", tx, "MeasurementNoiseNs", noiseNs, ...
    "RandomSeed", seed, "NumCaptures", 10));
end

function opts = baseSurveyOptsLocal()
opts = struct();
opts.ReferencePCI = 11;
opts.NumRandomStarts = 12;
opts.MinCaptures = 5;
opts.MaxFitRMSNsForPass = 200;
opts.MaxBiasUncertaintyForPassNs = 180;
opts.RequireInstabilityAssessableForPass = false;
end

function meas = addPerCaptureJitterLocal(meas, pci, jitterNs)
captureIDs = unique(meas.CaptureID, "stable");
assert(numel(jitterNs) == numel(captureIDs), ...
    "Jitter vector must match number of captures.");
frameNs = 10e6;
for k = 1:numel(captureIDs)
    idx = meas.CaptureID == captureIDs(k) & meas.PCI == pci;
    meas.FramePhaseNs(idx) = mod(meas.FramePhaseNs(idx) + jitterNs(k), frameNs);
end
end
