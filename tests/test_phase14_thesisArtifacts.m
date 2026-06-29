% test_phase14_thesisArtifacts
% Checks synthetic thesis artifact/report generation.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot, "src")));

fprintf("Running Phase 14 thesis artifact tests...\n");

tmpRoot = tempname;
mkdir(tmpRoot);
cleanup = onCleanup(@() rmdir(tmpRoot, "s"));

opts = struct();
opts.Survey = struct("NumRandomStarts", 4, "MinCaptures", 5);
opts.RSS = struct("GridSpacingM", 150, "SearchMarginM", 500, ...
    "MinMeasurementsPerPCI", 4, "MaxSuggestedStops", 4);

artifacts = generateThesisArtifacts(tmpRoot, opts);

assert(isfile(artifacts.IndexMarkdown));
assert(isfile(artifacts.SystemLimitsReport));
assert(isfile(artifacts.SurveyReportPaths.Markdown));
assert(isfile(artifacts.Figures.TimingThresholds));
assert(isfile(artifacts.Figures.GNBRequirements));
assert(isfile(artifacts.Figures.SyntheticSurvey));
assert(isfile(artifacts.Figures.TimingOffsets));
assert(isfile(artifacts.Figures.RSSPlanning));
assert(isfield(artifacts, "ThesisSyntheticResults"));
assert(isfile(artifacts.ThesisSyntheticResults.SummaryJSON));
assert(isfile(artifacts.ThesisSyntheticResults.SurveyCasesCSV));
assert(isfile(artifacts.ThesisSyntheticResults.ReceiverOverviewCSV));
assert(isfile(artifacts.ThesisSyntheticResults.Figures.WorstCaseEnvelope));
assert(artifacts.SurveyResult.FitInfo.Status ~= "NOT_ASSESSABLE");
assert(height(artifacts.RSSPlanning.RssPriors) >= 1);

testCustomThresholdReport(tmpRoot);

fprintf("Phase 14 thesis artifact tests passed.\n");

function testCustomThresholdReport(tmpRoot)
[meas, caps, ~] = buildSyntheticSurveyDataset("aligned_5gnb", struct( ...
    "MeasurementNoiseNs", 5, "RandomSeed", 1414));
surveyOpts = struct("NumRandomStarts", 4, "WarningOffsetNs", 123, ...
    "FailOffsetNs", 456, "HardFailOffsetNs", 789, ...
    "UncertaintySigmaMultiplier", 2);
result = noLocationSurveyCheck(meas, caps, surveyOpts);
result.IsSynthetic = true;
paths = generateSurveyReport(result, tmpRoot, "custom_threshold_smoke");
txt = string(fileread(paths.Markdown));
jsonTxt = string(fileread(paths.JSON));
assert(contains(txt, "Synthetic validation artifact"));
assert(contains(txt, "Static warning threshold: 123 ns"));
assert(contains(txt, "Static relative timing desynchronization threshold: 456 ns"));
assert(contains(txt, "Static hard-fail threshold: 789 ns"));
assert(contains(txt, "Uncertainty margin: 2 sigma"));
assert(contains(jsonTxt, "SyntheticValidationArtifact"));
assert(contains(jsonTxt, "123"));
end
