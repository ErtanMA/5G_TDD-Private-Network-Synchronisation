addpath(genpath(fullfile(fileparts(mfilename("fullpath")), "..", "src")));

fprintf("Running Phase 15 relative timing instability tests...\n");

surveyOpts = struct();
surveyOpts.NumRandomStarts = 12;
surveyOpts.MinCaptures = 5;
surveyOpts.MaxFitRMSNsForPass = 120;
surveyOpts.MaxBiasUncertaintyForPassNs = 180;
surveyOpts.ReferencePCI = 11;
surveyOpts.WarningInstabilityRmsNs = 100;
surveyOpts.FailInstabilityRmsNs = 500;
surveyOpts.WarningInstabilityPeakToPeakNs = 250;
surveyOpts.FailInstabilityPeakToPeakNs = 1000;

testAlignedSurveyHasNoInstability(surveyOpts);
testRelativeJitterFaultFails(surveyOpts);
testModerateRelativeJitterIsSuspect(surveyOpts);
testCommonModeMovementIsNotRelativeJitter(surveyOpts);

fprintf("Phase 15 relative timing instability tests passed.\n");

function testAlignedSurveyHasNoInstability(surveyOpts)
[meas, caps, ~] = buildSyntheticSurveyDataset("aligned_5gnb", struct( ...
    "MeasurementNoiseNs", 5, "RandomSeed", 1501, "NumCaptures", 10));
result = noLocationSurveyCheck(meas, caps, surveyOpts);

targets = result.TimingEstimates(result.TimingEstimates.SurveyRole == "TARGET", :);
assert(~any(targets.TimingInstabilityStatus == "FAIL"), ...
    "Aligned survey should not produce timing-instability failures.");
assert(all(targets.ExcessTimingJitterRmsNs < 80 | isnan(targets.ExcessTimingJitterRmsNs)), ...
    "Aligned survey produced excessive residual jitter.");
end

function testRelativeJitterFaultFails(surveyOpts)
[meas, caps, truth] = buildSyntheticSurveyDataset("aligned_5gnb", struct( ...
    "MeasurementNoiseNs", 5, "RandomSeed", 1502, "NumCaptures", 10));
faultPCI = truth.GNB.PCI(4);
jitterNs = [650; -520; 430; -610; 360; -560; 720; -380; 540; -640];
meas = addPerCaptureJitterLocal(meas, faultPCI, jitterNs);

result = noLocationSurveyCheck(meas, caps, surveyOpts);
row = result.TimingEstimates(result.TimingEstimates.PCI == faultPCI, :);

assert(row.TimingInstabilityStatus == "FAIL", ...
    "Large target-relative timing jitter should fail instability check.");
assert(row.TimingStatus == "FAIL", ...
    "Combined timing verdict should fail when instability fails.");
assert(row.ExcessTimingJitterRmsNs > 350, ...
    "Estimated excess jitter RMS should be large for injected jitter fault.");
end

function testModerateRelativeJitterIsSuspect(surveyOpts)
[meas, caps, truth] = buildSyntheticSurveyDataset("aligned_5gnb", struct( ...
    "MeasurementNoiseNs", 5, "RandomSeed", 1503, "NumCaptures", 10));
faultPCI = truth.GNB.PCI(4);
jitterNs = [180; -140; 90; -180; 150; -120; 160; -130; 120; -110];
meas = addPerCaptureJitterLocal(meas, faultPCI, jitterNs);

result = noLocationSurveyCheck(meas, caps, surveyOpts);
row = result.TimingEstimates(result.TimingEstimates.PCI == faultPCI, :);

assert(row.TimingInstabilityStatus == "SUSPECT", ...
    "Moderate target-relative timing jitter should be suspect, not clean pass.");
assert(row.TimingStatus == "SUSPECT", ...
    "Combined timing verdict should be suspect when instability is suspect.");
end

function testCommonModeMovementIsNotRelativeJitter(surveyOpts)
[meas, caps, ~] = buildSyntheticSurveyDataset("aligned_5gnb", struct( ...
    "MeasurementNoiseNs", 5, "RandomSeed", 1504, "NumCaptures", 10));
commonNs = [700; -550; 420; -630; 340; -510; 680; -360; 500; -600];
meas = addCommonModeJitterLocal(meas, commonNs);

result = noLocationSurveyCheck(meas, caps, surveyOpts);
targets = result.TimingEstimates(result.TimingEstimates.SurveyRole == "TARGET", :);

assert(~any(targets.TimingInstabilityStatus == "FAIL"), ...
    "Common-mode capture timing movement should cancel in relative measurements.");
assert(all(targets.ExcessTimingJitterRmsNs < 80 | isnan(targets.ExcessTimingJitterRmsNs)), ...
    "Common-mode timing movement should not appear as relative excess jitter.");
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

function meas = addCommonModeJitterLocal(meas, jitterNs)
captureIDs = unique(meas.CaptureID, "stable");
assert(numel(jitterNs) == numel(captureIDs), ...
    "Jitter vector must match number of captures.");
frameNs = 10e6;
for k = 1:numel(captureIDs)
    idx = meas.CaptureID == captureIDs(k);
    meas.FramePhaseNs(idx) = mod(meas.FramePhaseNs(idx) + jitterNs(k), frameNs);
end
end
