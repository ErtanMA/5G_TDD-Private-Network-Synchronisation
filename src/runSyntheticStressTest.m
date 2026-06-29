function [summary, caseResults] = runSyntheticStressTest(opts)
%runSyntheticStressTest Run survey-only synthetic stress tests.
%
%   [summary, caseResults] = runSyntheticStressTest(opts) validates the
%   supported thesis workflow: unknown gNB locations, known receiver survey
%   positions, relative transmit-timing estimation, and receiver PCI/timing
%   extraction before survey assembly.

arguments
    opts struct = struct()
end

opts = mergeStructsLocal(defaultStressTestOptions(), opts);
rng(opts.RandomSeed, "twister");

[surveyResults, surveySummary] = runSurveyStressLocal(opts);
if opts.RunReceiverHybridCases
    [receiverResults, receiverSummary] = runReceiverHybridStressLocal(opts);
else
    receiverResults = table();
    receiverSummary = struct("NumCases", 0);
end

caseResults = struct();
caseResults.SurveyLayer = surveyResults;
caseResults.ReceiverHybrid = receiverResults;

summary = struct();
summary.Mode = "survey_only_synthetic_stress";
summary.SurveyLayer = surveySummary;
summary.ReceiverHybrid = receiverSummary;
summary.Statement = "Stress tests use no gNB-location input. Survey-layer cases estimate unknown gNB positions and timing; receiver-hybrid cases test PCI/timing extraction before survey fitting.";

if opts.WriteArtifacts
    if strlength(string(opts.OutputDir)) == 0
        error("runSyntheticStressTest:MissingOutputDir", ...
            "OutputDir must be set when WriteArtifacts is true.");
    end
    writeStressArtifactsLocal(opts.OutputDir, summary, caseResults);
end

end

function [results, summary] = runSurveyStressLocal(opts)
rows = cell(opts.NumSurveyCases, 1);
surveyOpts = defaultSurveyOptions();
surveyOpts.NumRandomStarts = opts.SurveyRandomStarts;
surveyOpts.MinCaptures = min(opts.SurveyCaptureCountPool);

for caseIdx = 1:opts.NumSurveyCases
    scenarioName = opts.SurveyScenarios(mod(caseIdx-1, numel(opts.SurveyScenarios)) + 1);
    noiseNs = opts.SurveyNoisePoolNs(randi(numel(opts.SurveyNoisePoolNs)));
    numCaptures = opts.SurveyCaptureCountPool(randi(numel(opts.SurveyCaptureCountPool)));

    synthOpts = struct();
    synthOpts.MeasurementNoiseNs = noiseNs;
    synthOpts.NumCaptures = numCaptures;
    synthOpts.RandomSeed = opts.RandomSeed + caseIdx;

    [meas, caps, truth] = buildSyntheticSurveyDataset(scenarioName, synthOpts);
    result = noLocationSurveyCheck(meas, caps, surveyOpts);

    trueOffsets = medianCenteredTruthLocal(truth.GNB);
    est = joinTruthEstimateLocal(trueOffsets, result.TimingEstimates);
    maxAbsTruth = max(abs(est.TrueCenteredOffsetNs), [], "omitnan");
    maxAbsError = max(abs(est.EstimatedRelativeTxOffsetNs - est.TrueCenteredOffsetNs), [], "omitnan");
    expectedHardFault = maxAbsTruth >= surveyOpts.HardFailOffsetNs;
    actualFail = any(est.TimingStatus == "FAIL");
    actualNotAssessable = any(est.TimingStatus == "NOT_ASSESSABLE") || result.FitInfo.Status == "NOT_ASSESSABLE";
    commonMode = string(scenarioName) == "common_mode_3us";

    rows{caseIdx} = table( ...
        caseIdx, string(scenarioName), noiseNs, numCaptures, ...
        result.FitInfo.Status, result.FitInfo.FitRMSNs, ...
        result.FitInfo.ConditionNumber, maxAbsTruth, maxAbsError, ...
        expectedHardFault, actualFail, actualNotAssessable, commonMode, ...
        verdictListLocal(result.TimingEstimates), ...
        'VariableNames', ["CaseID","Scenario","NoiseNs","NumCaptures", ...
        "FitStatus","FitRMSNs","ConditionNumber", ...
        "MaxTrueCenteredOffsetNs","MaxAbsoluteEstimateErrorNs", ...
        "ExpectedHardFault","ActualAnyFail","ActualAnyNotAssessable", ...
        "IsCommonMode","Verdicts"]);
end

results = vertcat(rows{:});

hardCases = results.ExpectedHardFault & ~results.ActualAnyNotAssessable;
commonModeCases = results.IsCommonMode & ~results.ActualAnyNotAssessable;
summary = struct();
summary.NumCases = height(results);
summary.NumHardFaultCases = nnz(hardCases);
summary.HardFaultDetectionRate = safeRateLocal(nnz(results.ActualAnyFail & hardCases), nnz(hardCases));
summary.NumCommonModeCases = nnz(commonModeCases);
summary.CommonModeFalseFailRate = safeRateLocal(nnz(results.ActualAnyFail & commonModeCases), nnz(commonModeCases));
summary.NumNotAssessableCases = nnz(results.ActualAnyNotAssessable);
summary.MedianAbsEstimateErrorNs = median(results.MaxAbsoluteEstimateErrorNs, "omitnan");
end

function [results, summary] = runReceiverHybridStressLocal(opts)
rows = {};
rowCount = 0;

for caseName = opts.ReceiverHybridCases
    [pcis, offsets] = hybridCaseDefinitionLocal(caseName);
    for snrDb = opts.ReceiverHybridSNRdB
        for cfoHz = opts.ReceiverHybridCFOHz
            [sssCombined, truthCombined, meta] = decodeHybridIsolatedGNBsLocal(pcis, offsets, snrDb, cfoHz);
            [cellTiming, ~] = estimateCellTiming(sssCombined, meta);

            recoveredPCIs = unique(sssCombined.PCI(sssCombined.IsUsable));
            pciRecoveryOK = all(ismember(pcis(:), recoveredPCIs(:)));

            timingRows = timingRowsFromCellTimingLocal(cellTiming);
            expected = table(truthCombined.PCI(:), truthCombined.ExpectedRelativeOffsetNs(:), ...
                'VariableNames', ["PCI","ExpectedRelativeOffsetNs"]);
            joined = outerjoin(expected, timingRows, "Keys", "PCI", ...
                "MergeKeys", true, "Type", "left");
            timingErrorNs = joined.RelativeArrivalOffsetNs - joined.ExpectedRelativeOffsetNs;
            maxTimingErrorNs = max(abs(timingErrorNs), [], "omitnan");
            expectedHard = max(abs(joined.ExpectedRelativeOffsetNs), [], "omitnan") >= 3000;

            rowCount = rowCount + 1;
            rows{rowCount,1} = table( ...
                string(caseName), snrDb, cfoHz, numel(pcis), pciRecoveryOK, ...
                expectedHard, maxTimingErrorNs, ...
                strjoin(string(recoveredPCIs), ";"), ...
                'VariableNames', ["CaseName","SNRdB","CFOHz","NumGNBs", ...
                "PCIRecoveryOK","ExpectedAnyHardFault","MaxTimingErrorNs", ...
                "RecoveredPCIs"]);
        end
    end
end

if isempty(rows)
    results = table();
else
    results = vertcat(rows{:});
end

summary = struct();
summary.NumCases = height(results);
if height(results) == 0
    summary.PCIRecoveryRate = NaN;
    summary.MedianMaxTimingErrorNs = NaN;
else
    summary.PCIRecoveryRate = mean(results.PCIRecoveryOK);
    summary.MedianMaxTimingErrorNs = median(results.MaxTimingErrorNs, "omitnan");
end
end

function [sssCombined, truthCombined, meta] = decodeHybridIsolatedGNBsLocal(pcis, offsetsNs, snrDb, cfoHz)
sssRows = {};
truthRows = {};
meta = struct();
for k = 1:numel(pcis)
    custom = struct();
    custom.NumGNBs = 1;
    custom.PCIs = pcis(k);
    custom.FrameOffsetsNs = offsetsNs(k);
    custom.CFOHz = cfoHz;
    custom.GNBPowerdB = 0;
    custom.LocationUncertaintyM = 30;
    custom.SiteDistanceM = 100;
    custom.SNRdB = snrDb;
    custom.RandomSeed = 20260603 + k;

    [iq, meta, truth] = generateSyntheticScenario("aligned_5gnb", custom);
    pss = estimatePSSCandidates(iq, meta);
    sss = detectSSSAndPCI(iq, pss, meta);
    sss = sss(sss.IsUsable & sss.PCI == pcis(k), :);
    sssRows{end+1,1} = sss; %#ok<AGROW>
    truthRows{end+1,1} = truth; %#ok<AGROW>
end

sssCombined = vertcat(sssRows{:});
truthCombined = vertcat(truthRows{:});
end

function [pcis, offsets] = hybridCaseDefinitionLocal(caseName)
switch string(caseName)
    case "hybrid_aligned_3gnb"
        pcis = [11; 104; 257];
        offsets = [0; 0; 0];
    case "hybrid_250ns_one"
        pcis = [11; 104; 257];
        offsets = [0; 250; 0];
    case "hybrid_3us_one"
        pcis = [11; 104; 257];
        offsets = [0; 3000; 0];
    case "hybrid_common_mode_3us"
        pcis = [11; 104; 257];
        offsets = [3000; 3000; 3000];
    otherwise
        error("runSyntheticStressTest:UnknownHybridCase", ...
            "Unknown receiver-hybrid case: %s", caseName);
end
end

function truthCentered = medianCenteredTruthLocal(gnbTruth)
truthCentered = gnbTruth(:, ["PCI","RelativeTxOffsetNs"]);
center = median(truthCentered.RelativeTxOffsetNs, "omitnan");
truthCentered.TrueCenteredOffsetNs = truthCentered.RelativeTxOffsetNs - center;
truthCentered.RelativeTxOffsetNs = [];
end

function joined = joinTruthEstimateLocal(trueOffsets, estimates)
keep = estimates(:, ["PCI","EstimatedRelativeTxOffsetNs","TimingStatus"]);
joined = outerjoin(trueOffsets, keep, "Keys", "PCI", "MergeKeys", true, "Type", "left");
end

function rows = timingRowsFromCellTimingLocal(cellTiming)
if isempty(cellTiming) || height(cellTiming) == 0
    rows = table(zeros(0,1), zeros(0,1), ...
        'VariableNames', ["PCI","RelativeArrivalOffsetNs"]);
    return;
end
rows = cellTiming(:, ["PCI","RelativeArrivalOffsetNs"]);
end

function s = verdictListLocal(timingEstimates)
if isempty(timingEstimates) || height(timingEstimates) == 0
    s = "";
else
    s = strjoin(string(timingEstimates.PCI) + ":" + string(timingEstimates.TimingStatus), ";");
end
end

function r = safeRateLocal(num, den)
if den == 0
    r = NaN;
else
    r = num / den;
end
end

function writeStressArtifactsLocal(outputDir, summary, caseResults)
if ~isfolder(outputDir)
    mkdir(outputDir);
end

save(fullfile(outputDir, "synthetic_stress_result.mat"), "summary", "caseResults", "-v7.3");
writetable(caseResults.SurveyLayer, fullfile(outputDir, "synthetic_stress_survey_layer.csv"));
if istable(caseResults.ReceiverHybrid) && ~isempty(caseResults.ReceiverHybrid)
    writetable(caseResults.ReceiverHybrid, fullfile(outputDir, "synthetic_stress_receiver_hybrid.csv"));
end

fid = fopen(fullfile(outputDir, "synthetic_stress_summary.json"), "w");
if fid < 0
    error("runSyntheticStressTest:OpenFailed", "Could not write stress summary JSON.");
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s", jsonencode(summary, "PrettyPrint", true));
clear cleanup;
end

function opts = defaultStressTestOptions()
% Defaults for survey-only stress testing. Kept local so the stress harness is
% self-contained instead of scattering small option files across src.
opts = struct();
opts.RandomSeed = 20260603;
opts.NumSurveyCases = 80;
opts.SurveyScenarios = ["aligned_5gnb", "offset_250ns_one", ...
    "offset_1us_one", "offset_3us_one", "offset_3us_two", ...
    "common_mode_3us", "noisy_250ns_one", "weak_geometry_3us"];
opts.SurveyNoisePoolNs = [5 10 20 50 80];
opts.SurveyCaptureCountPool = [5 6 8 10];
opts.SurveyRandomStarts = 10;
opts.RunReceiverHybridCases = true;
opts.ReceiverHybridSNRdB = [35 20 10];
opts.ReceiverHybridCFOHz = [0 250 -250];
opts.ReceiverHybridCases = ["hybrid_aligned_3gnb", ...
    "hybrid_250ns_one", "hybrid_3us_one", "hybrid_common_mode_3us"];
opts.OutputDir = "";
opts.WriteArtifacts = false;
end

function opts = defaultSurveyOptions()
% Survey thresholds used by synthetic stress cases.
opts = struct();
opts.SpeedOfLightMps = 299792458;
opts.FramePeriodMs = 10;
opts.ReferencePCI = NaN;
opts.MinVisibleGNBs = 2;
opts.MinCaptures = 5;
opts.MinMeasurementsPerPCI = 4;
opts.RequireReferenceInEveryCapture = true;
opts.TimingStdNsDefault = 30;
opts.TimingStdFloorNs = 10;
opts.ReceiverPositionUncertaintyM = 5;
opts.PositionSearchRadiusM = 1500;
opts.NumRandomStarts = 24;
opts.RandomSeed = 20260603;
opts.UseLsqnonlinIfAvailable = true;
opts.MaxIterations = 2500;
opts.FunctionTolerance = 1e-10;
opts.MaxFitRMSNsForPass = 150;
opts.MaxBiasUncertaintyForPassNs = 150;
opts.MaxConditionNumber = 1e10;
opts.WarningOffsetNs = 250;
opts.FailOffsetNs = 1000;
opts.HardFailOffsetNs = 3000;
opts.UncertaintySigmaMultiplier = 3;
end

function out = mergeStructsLocal(base, override)
out = base;
if isempty(fieldnames(override))
    return;
end
names = fieldnames(override);
for k = 1:numel(names)
    out.(names{k}) = override.(names{k});
end
end
