% test_phase4_pssSearch
% Phase 4 smoke tests for handmade PSS matched-filter search.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot, "src")));

% Single-gNB case: PSS peaks should repeat every 10 ms.
[iqSingle, metaSingle, truthSingle] = generateSyntheticScenario("single_gnb_not_assessable");
[candSingle, dbgSingle] = estimatePSSCandidates(iqSingle, metaSingle);

assert(~isempty(candSingle), "Expected PSS candidates for single-gNB scenario.");
expectedNID2 = mod(truthSingle.PCI(1), 3);
assert(any(candSingle.NID2 == expectedNID2), "Expected correct NID2 in candidates.");
assert(dbgSingle.MetricMax(expectedNID2 + 1) > dbgSingle.MetricMedian(expectedNID2 + 1));

frameSamples = round(10e-3 * metaSingle.SampleRateHz);
expectedStarts = truthSingle.SSBStartOffsetSamples(1) + (0:3).' * frameSamples;
tolSamples = 3;
for k = 1:numel(expectedStarts)
    sameNID2 = candSingle(candSingle.NID2 == expectedNID2, :);
    minErr = min(abs(sameNID2.StartSample0 - expectedStarts(k)));
    assert(minErr <= tolSamples, ...
        "Expected PSS peak near frame start %.0f samples; min error was %.3f.", ...
        expectedStarts(k), minErr);
end

% Custom two-gNB timing case: one NID2 at 0 ns, one NID2 at +3 us.
custom = struct();
custom.NumGNBs = 2;
custom.PCIs = [11 12];
custom.FrameOffsetsNs = [0 3000];
custom.CFOHz = [0 0];
custom.GNBPowerdB = [0 0];
custom.LocationUncertaintyM = [30 30];
custom.SiteDistanceM = [100 100];
custom.SNRdB = 30;

[iqTwo, metaTwo, truthTwo] = generateSyntheticScenario("aligned_5gnb", custom);
[candTwo, ~] = estimatePSSCandidates(iqTwo, metaTwo);

for g = 1:height(truthTwo)
    expectedNID2 = mod(truthTwo.PCI(g), 3);
    expectedStart = truthTwo.SSBStartOffsetSamples(g) + ...
        truthTwo.FrameOffsetNs(g) * 1e-9 * metaTwo.SampleRateHz;
    sameNID2 = candTwo(candTwo.NID2 == expectedNID2, :);
    assert(~isempty(sameNID2), "Expected candidates for NID2 %d.", expectedNID2);
    minErr = min(abs(sameNID2.StartSample0 - expectedStart));
    assert(minErr <= tolSamples, ...
        "Expected NID2 %d PSS peak near %.3f samples; min error %.3f.", ...
        expectedNID2, expectedStart, minErr);
end

% Low-SNR case should still produce at least one plausible candidate.
[iqLow, metaLow, ~] = generateSyntheticScenario("low_snr_250ns");
candLow = estimatePSSCandidates(iqLow, metaLow, struct("MinPeakToMedianRatio", 8));
assert(~isempty(candLow), "Expected at least one PSS candidate in low-SNR scenario.");

disp("Phase 4 PSS search tests passed.");
