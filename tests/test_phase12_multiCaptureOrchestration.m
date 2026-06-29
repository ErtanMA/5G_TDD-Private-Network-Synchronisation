% test_phase12_multiCaptureOrchestration
% Phase 12 tests for real-data-style multi-capture orchestration.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot, "src")));

fprintf("Running Phase 12 multi-capture orchestration tests...\n");

tmpRoot = tempname;
mkdir(tmpRoot);
cleanup = onCleanup(@() rmdir(tmpRoot, "s"));

[manifest, truth] = writeSurveyIQFixture(tmpRoot);

opts = struct();
opts.Analyze = struct("SIC", struct("MaxIterations", 3));
opts.Survey = struct("NumRandomStarts", 6, "MinCaptures", 5, ...
    "MaxFitRMSNsForPass", 200);

result = runMultiCaptureSurvey(manifest, opts);

assert(result.Mode == "multi_capture_survey_from_iq");
assert(height(result.LoadSummary) == height(manifest));
assert(height(result.MeasurementTable) >= height(manifest) * 2, ...
    "Expected at least two detected PCI timing rows per capture.");
assert(isfield(result, "SurveyResult"));
assert(height(result.SurveyResult.TimingEstimates) >= 2);

detectedPCIs = unique(result.MeasurementTable.PCI);
assert(all(ismember(truth.PCIs(:), detectedPCIs(:))), ...
    "Expected orchestration to recover all synthetic fixture PCIs.");

fprintf("Phase 12 multi-capture orchestration tests passed.\n");

function [manifest, truth] = writeSurveyIQFixture(tmpRoot)
pcis = [11; 12; 13];
gnbXY = [-350 -150; 300 -250; 120 420];
txOffsetsNs = [0; 3000; 0];
rxXY = [ ...
    -900 -700; ...
     900 -700; ...
     950  650; ...
    -850  750; ...
       0 -950; ...
       0  900];
c = 299792458;

rows = cell(size(rxXY,1), 1);
for m = 1:size(rxXY,1)
    commonPhaseNs = 700 * m;
    distances = sqrt(sum((gnbXY - rxXY(m,:)).^2, 2));
    arrivalNs = commonPhaseNs + txOffsetsNs + distances / c * 1e9;

    custom = struct();
    custom.NumGNBs = numel(pcis);
    custom.PCIs = pcis.';
    custom.FrameOffsetsNs = arrivalNs.';
    custom.CFOHz = [0 0 0];
    custom.GNBPowerdB = [0 -1 -2];
    custom.LocationUncertaintyM = 30 * ones(1, numel(pcis));
    custom.SiteDistanceM = 100 * ones(1, numel(pcis));
    custom.SNRdB = 38;
    custom.DurationMs = 20;
    custom.RandomSeed = 7000 + m;

    [iq, meta, capTruth] = generateSyntheticScenario("aligned_5gnb", custom);
    paths = writeSyntheticCapture(tmpRoot, "cap_" + string(m), iq, meta, capTruth);

    rows{m} = table( ...
        "cap_" + string(m), string(paths.IQ), string(paths.Metadata), ...
        rxXY(m,1), rxXY(m,2), 0, 5, ...
        'VariableNames', ["CaptureID","IQPath","MetadataPath", ...
        "RxX_m","RxY_m","RxZ_m","RxPositionUncertaintyM"]);
end

manifest = vertcat(rows{:});
truth = struct();
truth.PCIs = pcis;
truth.GNBXY = gnbXY;
truth.TxOffsetsNs = txOffsetsNs;
truth.RxXY = rxXY;
end
