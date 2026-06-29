% test_syntheticStressHarness
% Smoke test for the survey-only synthetic stress-test harness.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot, "src")));

opts = struct();
opts.NumSurveyCases = 16;
opts.SurveyRandomStarts = 6;
opts.SurveyNoisePoolNs = [5 20];
opts.SurveyCaptureCountPool = [5 8];
opts.ReceiverHybridSNRdB = [35 20];
opts.ReceiverHybridCFOHz = [0 250];
opts.ReceiverHybridCases = ["hybrid_aligned_3gnb", "hybrid_3us_one", "hybrid_common_mode_3us"];

[summary, caseResults] = runSyntheticStressTest(opts);

assert(summary.SurveyLayer.NumCases == opts.NumSurveyCases);
assert(summary.SurveyLayer.NumHardFaultCases > 0);
assert(summary.SurveyLayer.HardFaultDetectionRate > 0);
assert(summary.SurveyLayer.NumCommonModeCases > 0);
assert(summary.SurveyLayer.CommonModeFalseFailRate == 0 | isnan(summary.SurveyLayer.CommonModeFalseFailRate));

assert(height(caseResults.ReceiverHybrid) == ...
    numel(opts.ReceiverHybridCases) * numel(opts.ReceiverHybridSNRdB) * numel(opts.ReceiverHybridCFOHz));
assert(summary.ReceiverHybrid.PCIRecoveryRate > 0.8);

commonMode = caseResults.SurveyLayer.IsCommonMode;
assert(~any(caseResults.SurveyLayer.ActualAnyFail(commonMode)), ...
    "Common-mode timing shifts should not fail relative survey mode.");

disp("Synthetic stress harness tests passed.");
