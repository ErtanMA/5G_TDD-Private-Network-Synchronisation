% test_phase17_thesisSyntheticResults
% Checks thesis-ready synthetic result generation and combined-envelope verdicts.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot, "src")));

fprintf("Running Phase 17 thesis synthetic result tests...\n");

tmpRoot = tempname;
mkdir(tmpRoot);
cleanup = onCleanup(@() rmdir(tmpRoot, "s"));

opts = struct();
opts.Survey = struct("NumRandomStarts", 6, "SyncThresholdNs", 1000, ...
    "MaxFitRMSNsForPass", 250, "RequireInstabilityAssessableForPass", false);

artifacts = runThesisSyntheticResults(tmpRoot, opts);

assert(isfile(artifacts.SummaryJSON));
assert(isfile(artifacts.SurveyCasesCSV));
assert(isfile(artifacts.SurveyEstimatesCSV));
assert(isfile(artifacts.ReceiverOverviewCSV));
assert(isfile(artifacts.Figures.ThresholdProfile));
assert(isfile(artifacts.Figures.StaticSweep));
assert(isfile(artifacts.Figures.WorstCaseEnvelope));
assert(isfile(artifacts.Figures.JitterResiduals));
assert(isfile(artifacts.Figures.CommonMode));
assert(isfile(artifacts.Figures.Assessability));

cases = artifacts.SurveyCases;
est = artifacts.SurveyEstimates;

assert(any(cases.CaseName == "static_plus_jitter_fail"));
combo = cases(cases.CaseName == "static_plus_jitter_fail", :);
assert(combo.ObservedOverallStatus == "FAIL", ...
    "Combined static-plus-jitter envelope should fail the 1 us synthetic threshold.");
assert(combo.TargetWorstCaseRelativeTimingNs > combo.SyncThresholdNs, ...
    "Combined envelope should exceed configured threshold in static-plus-jitter case.");

near = cases(cases.CaseName == "near_threshold_suspect", :);
assert(near.ObservedOverallStatus == "SUSPECT" || near.TargetTimingStatus == "SUSPECT", ...
    "Near-threshold case should be suspect due to uncertainty overlap.");

commonMode = cases(cases.CaseName == "common_mode_3us", :);
assert(commonMode.ObservedOverallStatus ~= "FAIL", ...
    "Common-mode 3 us timing shift must not fail relative survey mode.");

assert(all(ismember(["TimingInstabilityMaxAbsNs","WorstCaseRelativeTimingNs"], ...
    string(est.Properties.VariableNames))), ...
    "Per-PCI estimates must expose max residual and combined envelope fields.");

assert(artifacts.Summary.SyncThresholdNs == 1000);
assert(artifacts.Summary.NumSurveyCases >= 12);
assert(artifacts.Summary.ReceiverOverviewPassRate > 0.5);

fprintf("Phase 17 thesis synthetic result tests passed.\n");
