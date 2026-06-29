function artifacts = runThesisSyntheticResults(outputDir, opts)
%runThesisSyntheticResults Generate thesis-ready synthetic validation results.
%
%   artifacts = runThesisSyntheticResults(outputDir) writes deterministic
%   synthetic survey verdict results, a compact receiver-chain overview, and
%   figures intended for the thesis Results section. These are synthetic
%   validation artifacts, not real B210 measurement results.

arguments
    outputDir (1,1) string = "reports/thesis_artifacts"
    opts struct = struct()
end

opts = mergeStructsLocal(defaultThesisSyntheticOptionsLocal(), opts);
if ~isfolder(outputDir)
    mkdir(outputDir);
end
figureDir = fullfile(outputDir, "figures");
if ~isfolder(figureDir)
    mkdir(figureDir);
end

[surveyCases, surveyEstimates, residualRows] = runSurveyResultCasesLocal(opts);
receiverOverview = runReceiverOverviewLocal(opts);

surveyCsv = fullfile(outputDir, "synthetic_results_survey_cases.csv");
surveyEstCsv = fullfile(outputDir, "synthetic_results_survey_estimates.csv");
receiverCsv = fullfile(outputDir, "synthetic_results_receiver_overview.csv");
writetable(surveyCases, surveyCsv);
writetable(surveyEstimates, surveyEstCsv);
writetable(receiverOverview, receiverCsv);

figures = struct();
figures.ThresholdProfile = fullfile(figureDir, "synthetic_threshold_profile.png");
figures.StaticSweep = fullfile(figureDir, "synthetic_static_offset_sweep.png");
figures.WorstCaseEnvelope = fullfile(figureDir, "synthetic_worst_case_envelope.png");
figures.JitterResiduals = fullfile(figureDir, "synthetic_jitter_residuals.png");
figures.CommonMode = fullfile(figureDir, "synthetic_common_mode_invisibility.png");
figures.Assessability = fullfile(figureDir, "synthetic_assessability_cases.png");

plotThresholdProfileLocal(opts.Survey, figures.ThresholdProfile);
plotStaticSweepLocal(surveyCases, opts.Survey.SyncThresholdNs, figures.StaticSweep);
plotEnvelopeByPciLocal(surveyEstimates, "static_plus_jitter_fail", ...
    opts.Survey.SyncThresholdNs, figures.WorstCaseEnvelope);
plotResidualTraceLocal(residualRows, "jitter_only_fail", figures.JitterResiduals);
plotCommonModeLocal(surveyEstimates, "common_mode_3us", ...
    opts.Survey.SyncThresholdNs, figures.CommonMode);
plotAssessabilityLocal(surveyCases, figures.Assessability);

summary = buildSummaryLocal(opts, surveyCases, receiverOverview);
summaryPath = fullfile(outputDir, "synthetic_results_summary.json");
fid = fopen(summaryPath, "w");
if fid < 0
    error("runThesisSyntheticResults:OpenFailed", ...
        "Could not write synthetic summary JSON: %s", summaryPath);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s", jsonencode(summary, "PrettyPrint", true));
clear cleanup;

artifacts = struct();
artifacts.Mode = "thesis_synthetic_results";
artifacts.Statement = "Synthetic validation results for relative B210 survey mode. Not real measurement data.";
artifacts.OutputDir = outputDir;
artifacts.FigureDir = figureDir;
artifacts.SummaryJSON = summaryPath;
artifacts.SurveyCasesCSV = surveyCsv;
artifacts.SurveyEstimatesCSV = surveyEstCsv;
artifacts.ReceiverOverviewCSV = receiverCsv;
artifacts.Figures = figures;
artifacts.Summary = summary;
artifacts.SurveyCases = surveyCases;
artifacts.SurveyEstimates = surveyEstimates;
artifacts.ReceiverOverview = receiverOverview;

end

function [caseTable, estimateTable, residualTable] = runSurveyResultCasesLocal(opts)
caseNames = [
    "aligned_pass"
    "static_100ns"
    "static_250ns"
    "static_500ns"
    "static_1000ns"
    "static_1500ns"
    "static_3000ns"
    "near_threshold_suspect"
    "jitter_only_fail"
    "static_plus_jitter_fail"
    "common_mode_3us"
    "insufficient_captures"
    "one_visible_pci"
    "missing_reference_pci"
    "weak_geometry"
    ];

caseRows = cell(numel(caseNames), 1);
estimateRows = {};
residualRows = {};

for k = 1:numel(caseNames)
    [meas, caps, truth, meta] = buildSurveyCaseLocal(caseNames(k), opts, k);
    result = noLocationSurveyCheck(meas, caps, opts.Survey);
    result.IsSynthetic = true;
    result.ArtifactType = "THESIS_SYNTHETIC_RESULT";

    [caseRows{k}, targetPCI] = summarizeSurveyCaseLocal(caseNames(k), result, truth, meta, opts);

    est = result.TimingEstimates;
    if istable(est) && height(est) > 0
        est.CaseName = repmat(caseNames(k), height(est), 1);
        est = movevars(est, "CaseName", "Before", 1);
        estimateRows{end+1,1} = est; %#ok<AGROW>
    end

    res = residualRowsFromResultLocal(caseNames(k), result, targetPCI);
    if istable(res) && height(res) > 0
        residualRows{end+1,1} = res; %#ok<AGROW>
    end
end

caseTable = vertcat(caseRows{:});
if isempty(estimateRows)
    estimateTable = table();
else
    estimateTable = vertcat(estimateRows{:});
end
if isempty(residualRows)
    residualTable = table();
else
    residualTable = vertcat(residualRows{:});
end
end

function [meas, caps, truth, meta] = buildSurveyCaseLocal(caseName, opts, caseIdx)
base = struct();
base.NumCaptures = 10;
base.MeasurementNoiseNs = 5;
base.RandomSeed = opts.RandomSeed + caseIdx;
base.TxOffsetsNs = [0 0 0 0 0];
meta = struct();
meta.TargetPCI = 503;
meta.InjectedStaticOffsetNs = 0;
meta.InjectedMaxAbsJitterNs = 0;
meta.ExpectedOutcome = "PASS";
meta.Description = "";

switch string(caseName)
    case "aligned_pass"
        meta.Description = "Aligned five-cell survey.";
    case "static_100ns"
        base.TxOffsetsNs(4) = 100;
        meta.InjectedStaticOffsetNs = 100;
        meta.ExpectedOutcome = "PASS";
    case "static_250ns"
        base.TxOffsetsNs(4) = 250;
        meta.InjectedStaticOffsetNs = 250;
        meta.ExpectedOutcome = "PASS";
    case "static_500ns"
        base.TxOffsetsNs(4) = 500;
        meta.InjectedStaticOffsetNs = 500;
        meta.ExpectedOutcome = "PASS";
    case "static_1000ns"
        base.TxOffsetsNs(4) = 1000;
        base.MeasurementNoiseNs = 15;
        meta.InjectedStaticOffsetNs = 1000;
        meta.ExpectedOutcome = "SUSPECT_OR_FAIL";
    case "static_1500ns"
        base.TxOffsetsNs(4) = 1500;
        meta.InjectedStaticOffsetNs = 1500;
        meta.ExpectedOutcome = "FAIL";
    case "static_3000ns"
        base.TxOffsetsNs(4) = 3000;
        meta.InjectedStaticOffsetNs = 3000;
        meta.ExpectedOutcome = "FAIL";
    case "near_threshold_suspect"
        base.TxOffsetsNs(4) = 990;
        base.MeasurementNoiseNs = 20;
        meta.InjectedStaticOffsetNs = 990;
        meta.ExpectedOutcome = "SUSPECT";
    case "jitter_only_fail"
        meta.ExpectedOutcome = "FAIL";
        meta.Description = "Target PCI has capture-to-capture relative jitter.";
    case "static_plus_jitter_fail"
        base.TxOffsetsNs(4) = 800;
        meta.InjectedStaticOffsetNs = 800;
        meta.ExpectedOutcome = "FAIL";
    case "common_mode_3us"
        base.TxOffsetsNs(:) = 3000;
        meta.InjectedStaticOffsetNs = 0;
        meta.ExpectedOutcome = "NOT_FAIL";
    case "insufficient_captures"
        base.NumCaptures = 3;
        meta.ExpectedOutcome = "NOT_ASSESSABLE";
    case "one_visible_pci"
        meta.ExpectedOutcome = "NOT_ASSESSABLE";
    case "missing_reference_pci"
        meta.ExpectedOutcome = "NOT_ASSESSABLE";
    case "weak_geometry"
        base.RxPositionsM = [-900 -20; -450 10; 0 -10; 450 20; 900 -15];
        meta.ExpectedOutcome = "SUSPECT_OR_NOT_ASSESSABLE";
    otherwise
        error("runThesisSyntheticResults:UnknownSurveyCase", ...
            "Unknown survey case: %s", caseName);
end

[meas, caps, truth] = buildSyntheticSurveyDataset("aligned_5gnb", base);

switch string(caseName)
    case "jitter_only_fail"
        jitterNs = [650; -520; 430; -610; 360; -560; 720; -380; 540; -640];
        meas = addPerCaptureJitterLocal(meas, meta.TargetPCI, jitterNs);
        meta.InjectedMaxAbsJitterNs = max(abs(jitterNs));
    case "static_plus_jitter_fail"
        jitterNs = [720; -650; 610; -690; 750; -620; 700; -740; 630; -660];
        meas = addPerCaptureJitterLocal(meas, meta.TargetPCI, jitterNs);
        meta.InjectedMaxAbsJitterNs = max(abs(jitterNs));
    case "one_visible_pci"
        meas = meas(meas.PCI == truth.GNB.PCI(1), :);
    case "missing_reference_pci"
        refPCI = opts.Survey.ReferencePCI;
        dropCap = meas.CaptureID(1);
        meas = meas(~(meas.CaptureID == dropCap & meas.PCI == refPCI), :);
end
end

function [caseRow, targetPCI] = summarizeSurveyCaseLocal(caseName, result, truth, meta, opts)
targetPCI = meta.TargetPCI;
est = result.TimingEstimates;
target = table();
if istable(est) && height(est) > 0 && any(est.PCI == targetPCI)
    target = est(est.PCI == targetPCI, :);
end

overall = overallStatusLocal(result);
if isempty(target)
    targetStatus = "NOT_ASSESSABLE";
    targetStatic = NaN;
    targetJitter = NaN;
    targetEnvelope = NaN;
    targetUnc = NaN;
else
    targetStatus = string(target.TimingStatus(1));
    targetStatic = target.StaticOffsetNs(1);
    targetJitter = target.TimingInstabilityMaxAbsNs(1);
    targetEnvelope = target.WorstCaseRelativeTimingNs(1);
    targetUnc = target.WorstCaseRelativeTimingUncertaintyNs(1);
end

caseRow = table( ...
    string(caseName), targetPCI, opts.Survey.SyncThresholdNs, ...
    meta.InjectedStaticOffsetNs, meta.InjectedMaxAbsJitterNs, ...
    meta.ExpectedOutcome, overall, targetStatus, ...
    string(result.FitInfo.Status), string(result.FitInfo.Reason), ...
    targetStatic, targetJitter, targetEnvelope, targetUnc, ...
    outcomeMatchesLocal(meta.ExpectedOutcome, overall), ...
    truth.Statement, ...
    'VariableNames', ["CaseName","TargetPCI","SyncThresholdNs", ...
    "InjectedStaticOffsetNs","InjectedMaxAbsJitterNs", ...
    "ExpectedOutcome","ObservedOverallStatus","TargetTimingStatus", ...
    "FitStatus","FitReason","TargetStaticOffsetNs", ...
    "TargetMaxAbsInstabilityNs","TargetWorstCaseRelativeTimingNs", ...
    "TargetWorstCaseUncertaintyNs","OutcomeMatchesExpectation", ...
    "SyntheticTruthStatement"]);
end

function status = overallStatusLocal(result)
if ~isfield(result, "TimingEstimates") || isempty(result.TimingEstimates) || height(result.TimingEstimates) == 0
    status = "NOT_ASSESSABLE";
    return;
end
statuses = string(result.TimingEstimates.TimingStatus);
if any(statuses == "FAIL")
    status = "FAIL";
elseif any(statuses == "SUSPECT")
    status = "SUSPECT";
elseif all(statuses == "NOT_ASSESSABLE")
    status = "NOT_ASSESSABLE";
elseif any(statuses == "PASS")
    status = "PASS";
else
    status = statuses(1);
end
end

function tf = outcomeMatchesLocal(expected, observed)
switch string(expected)
    case "PASS"
        tf = observed == "PASS";
    case "FAIL"
        tf = observed == "FAIL";
    case "SUSPECT"
        tf = observed == "SUSPECT";
    case "NOT_ASSESSABLE"
        tf = observed == "NOT_ASSESSABLE";
    case "NOT_FAIL"
        tf = observed ~= "FAIL";
    case "SUSPECT_OR_FAIL"
        tf = observed == "SUSPECT" || observed == "FAIL";
    case "SUSPECT_OR_NOT_ASSESSABLE"
        tf = observed == "SUSPECT" || observed == "NOT_ASSESSABLE";
    otherwise
        tf = false;
end
end

function rows = residualRowsFromResultLocal(caseName, result, targetPCI)
if ~isfield(result, "RelativeMeasurements") || isempty(result.RelativeMeasurements) || ...
        height(result.RelativeMeasurements) == 0 || isempty(result.TimingEstimates)
    rows = table();
    return;
end

rel = result.RelativeMeasurements;
est = result.TimingEstimates;
refPCI = result.FitInfo.ReferencePCI;
ref = est(est.PCI == refPCI, :);
if isempty(ref)
    rows = table();
    return;
end

outRows = {};
rowCount = 0;
frameNs = 10e6;
for k = 1:height(rel)
    target = est(est.PCI == rel.PCI(k), :);
    if isempty(target)
        continue;
    end
    rx = [rel.RxX_m(k), rel.RxY_m(k)];
    targetXY = [target.EstimatedX_m(1), target.EstimatedY_m(1)];
    refXY = [ref.EstimatedX_m(1), ref.EstimatedY_m(1)];
    bRel = target.EstimatedRelativeTxOffsetNs(1) - ref.EstimatedRelativeTxOffsetNs(1);
    geomNs = (norm(targetXY-rx) - norm(refXY-rx)) / 299792458 * 1e9;
    residualNs = wrapNsLocal(rel.RelativeArrivalNs(k) - (bRel + geomNs), frameNs);
    rowCount = rowCount + 1;
    outRows{rowCount,1} = table( ...
        string(caseName), rel.CaptureID(k), rel.PCI(k), targetPCI, ...
        rel.RelativeArrivalNs(k), residualNs, ...
        'VariableNames', ["CaseName","CaptureID","PCI","TargetPCI", ...
        "RelativeArrivalNs","ResidualNs"]);
end

if isempty(outRows)
    rows = table();
    return;
end
rows = vertcat(outRows{:});
pcis = unique(rows.PCI, "stable");
rows.CenteredResidualNs = NaN(height(rows), 1);
for p = 1:numel(pcis)
    idx = rows.PCI == pcis(p);
    rows.CenteredResidualNs(idx) = rows.ResidualNs(idx) - mean(rows.ResidualNs(idx), "omitnan");
end
end

function overview = runReceiverOverviewLocal(opts)
rows = {
    receiverPciRecoveryCaseLocal("clean_pci_recovery", "aligned_5gnb", ...
        singleCellCustomLocal(35, 0), "Single-cell PCI recovery in clean synthetic IQ.")
    receiverPciRecoveryCaseLocal("moderate_snr_cfo_recovery", "aligned_5gnb", ...
        singleCellCustomLocal(20, 250), "Single-cell PCI recovery with moderate SNR and CFO.")
    receiverSicCaseLocal()
    receiverTddCaseLocal("tdd_dddsu_pass", "aligned_5gnb", "Known DDDSU x4 pattern passes.")
    receiverTddCaseLocal("wrong_tdd_pattern_fail", "wrong_tdd_pattern", "Expected uplink slots contain downlink-like energy.")
    receiverTddCaseLocal("wrong_special_slot_fail", "wrong_special_slot", "Special slot tail contains downlink-like energy.")
    };
overview = vertcat(rows{:});
end

function custom = singleCellCustomLocal(snrDb, cfoHz)
custom = struct();
custom.NumGNBs = 1;
custom.PCIs = 503;
custom.FrameOffsetsNs = 0;
custom.CFOHz = cfoHz;
custom.GNBPowerdB = 0;
custom.LocationUncertaintyM = 30;
custom.SiteDistanceM = 100;
custom.SNRdB = snrDb;
end

function row = receiverPciRecoveryCaseLocal(caseName, scenarioName, custom, purpose)
[iq, meta, truth] = generateSyntheticScenario(scenarioName, custom);
result = analyzeCapture(iq, meta, struct("SIC", struct("MaxIterations", 6)));
truthPCIs = unique(truth.PCI);
detected = unique(result.SSSDetections.PCI(result.SSSDetections.IsUsable));
passed = all(ismember(truthPCIs, detected));
row = receiverRowLocal(caseName, purpose, "PCI recovery", ...
    strjoin(string(truthPCIs), ";"), strjoin(string(detected), ";"), passed, ...
    "Synthetic IQ receiver overview.");
end

function row = receiverSicCaseLocal()
custom = struct();
custom.NumGNBs = 2;
custom.PCIs = [11 104];
custom.FrameOffsetsNs = [0 3000];
custom.CFOHz = [0 0];
custom.GNBPowerdB = [0 -15];
custom.LocationUncertaintyM = [30 30];
custom.SiteDistanceM = [100 100];
custom.SNRdB = 40;
[iq, meta, truth] = generateSyntheticScenario("aligned_5gnb", custom);
pssOneShot = estimatePSSCandidates(iq, meta);
sssOneShot = detectSSSAndPCI(iq, pssOneShot, meta);
oneShot = unique(sssOneShot.PCI(sssOneShot.IsUsable));
[sicDetections, ~] = detectCellsWithSIC(iq, meta, struct("MaxIterations", 4));
sicPCIs = unique(sicDetections.PCI(sicDetections.IsUsable));
passed = ~all(ismember(truth.PCI, oneShot)) && all(ismember(truth.PCI, sicPCIs));
observed = "one-shot=" + strjoin(string(oneShot), ";") + "; SIC=" + strjoin(string(sicPCIs), ";");
row = receiverRowLocal("sic_near_far_recovery", ...
    "Near-far same-NID2 weak-cell recovery with SIC.", ...
    "SIC recovery", strjoin(string(truth.PCI), ";"), observed, passed, ...
    "Timing refinement uses PCI-specific PSS+SSS only; PBCH DM-RS is not used for timing refinement.");
end

function row = receiverTddCaseLocal(caseName, scenarioName, purpose)
custom = struct();
custom.NumGNBs = 1;
custom.PCIs = 11;
custom.FrameOffsetsNs = 0;
custom.CFOHz = 0;
custom.GNBPowerdB = 0;
custom.LocationUncertaintyM = 30;
custom.SiteDistanceM = 100;
custom.SNRdB = 35;
[iq, meta, ~] = generateSyntheticScenario(scenarioName, custom);
pss = estimatePSSCandidates(iq, meta);
sss = detectSSSAndPCI(iq, pss, meta);
timing = estimateCellTiming(sss, meta);
[tdd, ~, special, ~] = checkTDDPattern(iq, meta, timing);
observed = "TDD=" + join(unique(tdd.TDDPatternStatus), ";") + ...
    "; Special=" + join(unique(special.SpecialSlotStatus), ";");
switch string(scenarioName)
    case "aligned_5gnb"
        passed = all(tdd.TDDPatternStatus == "PASS") && all(special.SpecialSlotStatus == "PASS");
        expected = "TDD PASS; Special PASS";
    case "wrong_tdd_pattern"
        passed = all(tdd.TDDPatternStatus == "FAIL");
        expected = "TDD FAIL";
    case "wrong_special_slot"
        passed = all(special.SpecialSlotStatus == "FAIL");
        expected = "Special FAIL";
end
row = receiverRowLocal(caseName, purpose, "TDD envelope check", expected, observed, passed, ...
    "Known-pattern energy check, not blind TDD-pattern discovery.");
end

function row = receiverRowLocal(caseName, purpose, metric, expected, observed, passed, note)
row = table(string(caseName), string(purpose), string(metric), ...
    string(expected), string(observed), logical(passed), string(note), ...
    'VariableNames', ["CaseName","Purpose","Metric","Expected", ...
    "Observed","Passed","Note"]);
end

function summary = buildSummaryLocal(opts, surveyCases, receiverOverview)
summary = struct();
summary.Mode = "thesis_synthetic_results";
summary.Statement = "Synthetic validation only. The 1 us default is an engineering/private-network threshold, not regulatory or absolute UTC compliance.";
summary.SyncThresholdNs = opts.Survey.SyncThresholdNs;
summary.NumSurveyCases = height(surveyCases);
summary.NumSurveyFailures = nnz(surveyCases.ObservedOverallStatus == "FAIL");
summary.NumSurveySuspect = nnz(surveyCases.ObservedOverallStatus == "SUSPECT");
summary.NumSurveyNotAssessable = nnz(surveyCases.ObservedOverallStatus == "NOT_ASSESSABLE");
summary.SurveyExpectationMatchRate = mean(surveyCases.OutcomeMatchesExpectation);
summary.NumReceiverOverviewCases = height(receiverOverview);
summary.ReceiverOverviewPassRate = mean(receiverOverview.Passed);
summary.MainVerdictMetric = "WorstCaseRelativeTimingNs = abs(static offset) + max absolute post-fit residual";
end

function plotThresholdProfileLocal(surveyOpts, outputPath)
sigma = [0; 10; 20; 50; 100; 150; 250; 400];
k = surveyOpts.UncertaintySigmaMultiplier;
sync = surveyOpts.SyncThresholdNs;
passUpTo = max(0, sync - k*sigma);
failFrom = sync + k*sigma;
fig = figure("Visible", "off", "Color", "w");
ax = axes(fig);
plot(ax, sigma, passUpTo, "-o", "LineWidth", 1.8, "DisplayName", "Clean PASS boundary");
hold(ax, "on");
plot(ax, sigma, failFrom, "-s", "LineWidth", 1.8, "DisplayName", "Definite FAIL boundary");
yline(ax, sync, "--", "Configured threshold", "LineWidth", 1.2, ...
    "HandleVisibility", "off");
grid(ax, "on");
xlabel(ax, "Timing uncertainty sigma (ns)");
ylabel(ax, "Worst-case relative timing envelope (ns)");
title(ax, "Configurable synchronization threshold with uncertainty overlap");
legend(ax, "Location", "northwest");
exportgraphics(fig, outputPath, "Resolution", 180);
close(fig);
end

function plotStaticSweepLocal(cases, syncThresholdNs, outputPath)
idx = startsWith(cases.CaseName, "static_") & cases.InjectedMaxAbsJitterNs == 0;
data = cases(idx, :);
data = sortrows(data, "InjectedStaticOffsetNs");
fig = figure("Visible", "off", "Color", "w");
ax = axes(fig);
plot(ax, data.InjectedStaticOffsetNs, data.TargetWorstCaseRelativeTimingNs, ...
    "-o", "LineWidth", 1.8, "MarkerSize", 6);
hold(ax, "on");
yline(ax, syncThresholdNs, "--", "Threshold", "LineWidth", 1.2);
grid(ax, "on");
xlabel(ax, "Injected static relative offset (ns)");
ylabel(ax, "Estimated worst-case relative timing (ns)");
title(ax, "Synthetic static-offset sweep");
exportgraphics(fig, outputPath, "Resolution", 180);
close(fig);
end

function plotEnvelopeByPciLocal(estimates, caseName, syncThresholdNs, outputPath)
rows = estimates(estimates.CaseName == string(caseName), :);
fig = figure("Visible", "off", "Color", "w");
ax = axes(fig);
if isempty(rows)
    text(ax, 0.5, 0.5, "No envelope rows", "HorizontalAlignment", "center");
else
    bar(ax, categorical(string(rows.PCI)), rows.WorstCaseRelativeTimingNs, ...
        "FaceColor", [0.3 0.55 0.85]);
    hold(ax, "on");
    yline(ax, syncThresholdNs, "--", "Threshold", "LineWidth", 1.2);
    grid(ax, "on");
    xlabel(ax, "PCI");
    ylabel(ax, "Worst-case relative timing (ns)");
    title(ax, "Worst-case timing envelope per PCI");
end
exportgraphics(fig, outputPath, "Resolution", 180);
close(fig);
end

function plotResidualTraceLocal(residuals, caseName, outputPath)
rows = residuals(residuals.CaseName == string(caseName) & residuals.PCI == residuals.TargetPCI, :);
fig = figure("Visible", "off", "Color", "w");
ax = axes(fig);
if isempty(rows)
    text(ax, 0.5, 0.5, "No residual rows", "HorizontalAlignment", "center");
else
    plot(ax, 1:height(rows), rows.CenteredResidualNs, "-o", "LineWidth", 1.7);
    grid(ax, "on");
    xlabel(ax, "Capture index");
    ylabel(ax, "Centered residual timing (ns)");
    title(ax, "Injected relative timing instability after geometry fit");
end
exportgraphics(fig, outputPath, "Resolution", 180);
close(fig);
end

function plotCommonModeLocal(estimates, caseName, syncThresholdNs, outputPath)
rows = estimates(estimates.CaseName == string(caseName), :);
fig = figure("Visible", "off", "Color", "w");
ax = axes(fig);
if isempty(rows)
    text(ax, 0.5, 0.5, "No common-mode rows", "HorizontalAlignment", "center");
else
    bar(ax, categorical(string(rows.PCI)), abs(rows.EstimatedRelativeTxOffsetNs), ...
        "FaceColor", [0.25 0.65 0.45]);
    hold(ax, "on");
    yline(ax, syncThresholdNs, "--", "Threshold", "LineWidth", 1.2);
    grid(ax, "on");
    xlabel(ax, "PCI");
    ylabel(ax, "Estimated relative offset magnitude (ns)");
    title(ax, "Common-mode 3 us shift is invisible to relative mode");
end
exportgraphics(fig, outputPath, "Resolution", 180);
close(fig);
end

function plotAssessabilityLocal(cases, outputPath)
names = ["aligned_pass"; "insufficient_captures"; "one_visible_pci"; ...
    "missing_reference_pci"; "weak_geometry"];
idx = ismember(cases.CaseName, names);
data = cases(idx, :);
score = zeros(height(data), 1);
for k = 1:height(data)
    switch string(data.ObservedOverallStatus(k))
        case "PASS"
            score(k) = 1;
        case "SUSPECT"
            score(k) = 0.5;
        case "FAIL"
            score(k) = -0.5;
        otherwise
            score(k) = 0;
    end
end
fig = figure("Visible", "off", "Color", "w");
ax = axes(fig);
bar(ax, categorical(data.CaseName), score, "FaceColor", [0.45 0.45 0.75]);
ylim(ax, [-0.75 1.25]);
yticks(ax, [-0.5 0 0.5 1]);
yticklabels(ax, ["FAIL","NOT ASSESSABLE","SUSPECT","PASS"]);
grid(ax, "on");
ylabel(ax, "Observed survey outcome");
title(ax, "Assessability and geometry edge cases");
exportgraphics(fig, outputPath, "Resolution", 180);
close(fig);
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

function y = wrapNsLocal(x, periodNs)
y = mod(x + periodNs/2, periodNs) - periodNs/2;
end

function opts = defaultThesisSyntheticOptionsLocal()
opts = struct();
opts.RandomSeed = 20260612;
opts.Survey = struct();
opts.Survey.ReferencePCI = 11;
opts.Survey.SyncThresholdNs = 1000;
opts.Survey.NumRandomStarts = 16;
opts.Survey.MinCaptures = 5;
opts.Survey.MaxFitRMSNsForPass = 200;
opts.Survey.MaxBiasUncertaintyForPassNs = 180;
opts.Survey.RequireInstabilityAssessableForPass = false;
opts.Survey.UncertaintySigmaMultiplier = 3;
end

function out = mergeStructsLocal(base, override)
out = base;
if isempty(fieldnames(override))
    return;
end
names = fieldnames(override);
for k = 1:numel(names)
    if isstruct(override.(names{k})) && isfield(out, names{k}) && isstruct(out.(names{k}))
        out.(names{k}) = mergeStructsLocal(out.(names{k}), override.(names{k}));
    else
        out.(names{k}) = override.(names{k});
    end
end
end
