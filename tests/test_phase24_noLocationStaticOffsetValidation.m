% test_phase24_noLocationStaticOffsetValidation
% Checks model-level no-location static-offset validation artifacts.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot, "src")));

fprintf("Running Phase 24 no-location static-offset validation tests...\n");

tmpRoot = tempname;
mkdir(tmpRoot);
cleanup = onCleanup(@() rmdir(tmpRoot, "s"));

opts = struct();
opts.Survey = struct("NumRandomStarts", 8, "SyncThresholdNs", 1000, ...
    "MaxFitRMSNsForPass", 250, "RequireInstabilityAssessableForPass", false);

artifacts = runNoLocationStaticOffsetValidation(tmpRoot, opts);

assert(isfile(artifacts.CasesCSV));
assert(isfile(artifacts.EstimatesCSV));
assert(isfile(artifacts.GeometryCSV));
assert(isfile(artifacts.SummaryJSON));
assert(isfile(artifacts.Figures.InjectedVsEstimated));
assert(isfile(artifacts.Figures.AbsoluteError));
assert(isfile(artifacts.Figures.GeometryComparison));
assert(isfile(artifacts.Figures.VerdictSummary));

cases = artifacts.Cases;
est = artifacts.Estimates;

good3us = cases(cases.CaseName == "good_static_3000ns", :);
assert(good3us.ObservedOverallStatus == "FAIL", ...
    "Good-geometry 3 us static offset should fail.");
assert(good3us.TargetAbsoluteErrorNs < 400, ...
    "Good-geometry 3 us offset should be recovered close to injected truth.");

common = cases(cases.CaseName == "common_mode_3000ns", :);
assert(common.ObservedOverallStatus ~= "FAIL", ...
    "Common-mode 3 us shift must not fail relative no-location mode.");

tooFew = cases(cases.CaseName == "too_few_positions_1000ns", :);
assert(tooFew.ObservedOverallStatus == "NOT_ASSESSABLE", ...
    "Too few positions should be not assessable.");

weak = cases(cases.CaseName == "weak_collinear_1000ns", :);
assert(ismember(weak.ObservedOverallStatus, ["SUSPECT","NOT_ASSESSABLE"]), ...
    "Weak geometry should not produce a clean pass/fail overclaim.");

assert(all(ismember(["TrueCenteredOffsetNs","EstimatedRelativeTxOffsetNs", ...
    "AbsoluteErrorNs"], string(est.Properties.VariableNames))), ...
    "Estimate table must compare hidden injected truth against estimates.");

assert(artifacts.Summary.SyncThresholdNs == 1000);
assert(artifacts.Summary.ExpectationMatchRate > 0.70);

fprintf("Phase 24 no-location static-offset validation tests passed.\n");
