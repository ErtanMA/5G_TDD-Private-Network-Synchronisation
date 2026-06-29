% test_phase10_sicNearFar
% Phase 10 tests for overlapping/near-far SSB robustness.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot, "src")));

fprintf("Running Phase 10 SIC near-far tests...\n");

testNearFarSameNID2Recovery();
testAnalyzeCaptureUsesSIC();
testSICCanBeDisabled();
testPBCHDMRSReferenceIsUsed();

fprintf("Phase 10 SIC near-far tests passed.\n");

function testNearFarSameNID2Recovery()
[iq, meta, truth] = nearFarScenario();

pssOneShot = estimatePSSCandidates(iq, meta);
sssOneShot = detectSSSAndPCI(iq, pssOneShot, meta);
oneShotPCIs = unique(sssOneShot.PCI(sssOneShot.IsUsable));

[sicDetections, sicDebug] = detectCellsWithSIC(iq, meta, struct("MaxIterations", 4));
sicPCIs = unique(sicDetections.PCI(sicDetections.IsUsable));

assert(ismember(truth.PCI(1), oneShotPCIs), ...
    "Expected one-shot detector to recover the strong PCI.");
assert(~ismember(truth.PCI(2), oneShotPCIs), ...
    "This synthetic near-far case should hide the weak PCI in one-shot detection.");
assert(all(ismember(truth.PCI, sicPCIs)), ...
    "SIC should recover both strong and weak PCIs.");
assert(numel(unique(sicDebug.CancelledPCIs)) >= 2, ...
    "Expected at least two cancellation iterations.");

weakRows = sicDetections(sicDetections.PCI == truth.PCI(2), :);
assert(~isempty(weakRows), "Expected weak PCI rows after cancellation.");
expectedWeakStart = truth.SSBStartOffsetSamples(2) + ...
    truth.FrameOffsetNs(2) * 1e-9 * meta.SampleRateHz;
minErr = min(abs(weakRows.StartSample0 - expectedWeakStart));
assert(minErr < 4, "Expected weak PCI start timing near the injected 3 us offset.");
end

function testAnalyzeCaptureUsesSIC()
[iq, meta, truth] = nearFarScenario();
result = analyzeCapture(iq, meta);
pcis = unique(result.SSSDetections.PCI(result.SSSDetections.IsUsable));
assert(all(ismember(truth.PCI, pcis)), ...
    "analyzeCapture should use SIC by default and recover both PCIs.");
assert(isfield(result, "SICDebug"));
assert(numel(unique(result.SICDebug.CancelledPCIs)) >= 2);
end

function testSICCanBeDisabled()
[iq, meta, truth] = nearFarScenario();
result = analyzeCapture(iq, meta, struct("SIC", struct("Enable", false)));
pcis = unique(result.SSSDetections.PCI(result.SSSDetections.IsUsable));
assert(ismember(truth.PCI(1), pcis));
assert(~ismember(truth.PCI(2), pcis), ...
    "Disabling SIC should reproduce the one-shot near-far limitation.");
end

function testPBCHDMRSReferenceIsUsed()
[iq, meta, ~] = nearFarScenario();
[sicDetections, sicDebug] = detectCellsWithSIC(iq, meta, struct("MaxIterations", 2));
assert(~isempty(sicDetections), "Expected SIC detections for DM-RS usage test.");
assert(~isempty(sicDebug.Cancellations), "Expected cancellation rows.");
assert(all(sicDebug.Cancellations.IncludedPBCHDMRS), ...
    "Default SIC cancellation should include PBCH DM-RS.");
end

function [iq, meta, truth] = nearFarScenario()
custom = struct();
custom.NumGNBs = 2;
custom.PCIs = [11 104];          % Same NID2, different NID1.
custom.FrameOffsetsNs = [0 3000];
custom.CFOHz = [0 0];
custom.GNBPowerdB = [0 -15];
custom.LocationUncertaintyM = [30 30];
custom.SiteDistanceM = [100 100];
custom.SNRdB = 40;

[iq, meta, truth] = generateSyntheticScenario("aligned_5gnb", custom);
end
