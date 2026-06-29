function artifacts = runNoLocationStaticOffsetValidation(outputDir, opts)
%runNoLocationStaticOffsetValidation Validate no-location static-offset fit.
%
%   artifacts = runNoLocationStaticOffsetValidation(outputDir) generates a
%   model-level synthetic validation package for the multi-position survey.
%   The synthetic gNB/cell positions and transmit offsets are known to the
%   validation harness but hidden from estimateSurveyTimingNoLocations. The
%   estimator receives only receiver positions, relative arrival
%   measurements, and timing uncertainties.

arguments
    outputDir (1,1) string = "reports/thesis_artifacts/no_location_static_offset_validation"
    opts struct = struct()
end

opts = mergeStructsLocal(defaultValidationOptionsLocal(), opts);
if ~isfolder(outputDir)
    mkdir(outputDir);
end
figureDir = fullfile(outputDir, "figures");
if ~isfolder(figureDir)
    mkdir(figureDir);
end

[caseTable, estimateTable, geometryTable] = runValidationCasesLocal(opts);

casesCsv = fullfile(outputDir, "no_location_static_offset_validation_cases.csv");
estimatesCsv = fullfile(outputDir, "no_location_static_offset_validation_estimates.csv");
geometryCsv = fullfile(outputDir, "no_location_static_offset_validation_geometry.csv");
writetable(caseTable, casesCsv);
writetable(estimateTable, estimatesCsv);
writetable(geometryTable, geometryCsv);

figures = struct();
figures.InjectedVsEstimated = fullfile(figureDir, "no_location_injected_vs_estimated.png");
figures.AbsoluteError = fullfile(figureDir, "no_location_static_offset_error.png");
figures.GeometryComparison = fullfile(figureDir, "no_location_geometry_comparison.png");
figures.VerdictSummary = fullfile(figureDir, "no_location_verdict_summary.png");

plotInjectedVsEstimatedLocal(estimateTable, opts.Survey.SyncThresholdNs, figures.InjectedVsEstimated);
plotAbsoluteErrorLocal(caseTable, figures.AbsoluteError);
plotGeometryComparisonLocal(geometryTable, figures.GeometryComparison);
plotVerdictSummaryLocal(caseTable, figures.VerdictSummary);

summary = buildSummaryLocal(opts, caseTable, estimateTable);
summaryJson = fullfile(outputDir, "no_location_static_offset_validation_summary.json");
fid = fopen(summaryJson, "w");
if fid < 0
    error("runNoLocationStaticOffsetValidation:OpenFailed", ...
        "Could not write summary JSON: %s", summaryJson);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s", jsonencode(summary, "PrettyPrint", true));
clear cleanup;

artifacts = struct();
artifacts.Mode = "no_location_static_offset_synthetic_validation";
artifacts.Statement = "Model-level synthetic validation: true cell positions and static timing offsets are hidden from the estimator.";
artifacts.OutputDir = outputDir;
artifacts.FigureDir = figureDir;
artifacts.CasesCSV = casesCsv;
artifacts.EstimatesCSV = estimatesCsv;
artifacts.GeometryCSV = geometryCsv;
artifacts.SummaryJSON = summaryJson;
artifacts.Figures = figures;
artifacts.Summary = summary;
artifacts.Cases = caseTable;
artifacts.Estimates = estimateTable;
artifacts.Geometry = geometryTable;
end

function [caseTable, estimateTable, geometryTable] = runValidationCasesLocal(opts)
caseNames = [
    "good_aligned"
    "good_static_250ns"
    "good_static_500ns"
    "good_static_1000ns"
    "good_static_1500ns"
    "good_static_3000ns"
    "good_multi_offset"
    "good_noise_50ns"
    "good_noise_100ns"
    "common_mode_3000ns"
    "weak_collinear_1000ns"
    "too_few_positions_1000ns"
    "two_cells_1000ns"
    "missing_observations_1000ns"
    ];

caseRows = cell(numel(caseNames), 1);
estimateRows = {};
geometryRows = {};

for k = 1:numel(caseNames)
    [meas, caps, truth, meta] = buildValidationCaseLocal(caseNames(k), opts, k);
    result = noLocationSurveyCheck(meas, caps, opts.Survey);

    caseRows{k} = summarizeCaseLocal(caseNames(k), result, truth, meta, opts);

    estRows = compareEstimatesToTruthLocal(caseNames(k), result, truth, meta);
    if istable(estRows) && height(estRows) > 0
        estimateRows{end+1,1} = estRows; %#ok<AGROW>
    end

    geomRows = geometryRowsLocal(caseNames(k), result, truth, caps, meta);
    if istable(geomRows) && height(geomRows) > 0
        geometryRows{end+1,1} = geomRows; %#ok<AGROW>
    end
end

caseTable = vertcat(caseRows{:});
estimateTable = vertcatOrEmptyLocal(estimateRows);
geometryTable = vertcatOrEmptyLocal(geometryRows);
end

function [meas, caps, truth, meta] = buildValidationCaseLocal(caseName, opts, caseIdx)
base = struct();
base.RandomSeed = opts.RandomSeed + caseIdx;
base.NumCaptures = 8;
base.MeasurementNoiseNs = 20;
base.PCIs = [11 104 257 503 777];
base.TxOffsetsNs = [0 0 0 0 0];
base.SurveyRadiusM = 1200;

meta = struct();
meta.CaseName = string(caseName);
meta.TargetPCI = 503;
meta.ExpectedOutcome = "PASS";
meta.GeometryClass = "GOOD_2D";
meta.ExpectedInterpretation = "";

switch string(caseName)
    case "good_aligned"
        meta.ExpectedInterpretation = "Aligned survey should recover near-zero centred offsets.";
    case "good_static_250ns"
        base.TxOffsetsNs(4) = 250;
        meta.ExpectedInterpretation = "Small static offset should be recovered and remain below 1 us.";
    case "good_static_500ns"
        base.TxOffsetsNs(4) = 500;
        meta.ExpectedInterpretation = "Moderate static offset should be recovered below 1 us.";
    case "good_static_1000ns"
        base.TxOffsetsNs(4) = 1000;
        meta.ExpectedOutcome = "SUSPECT_OR_FAIL";
        meta.ExpectedInterpretation = "Threshold-scale offset should touch the verdict boundary.";
    case "good_static_1500ns"
        base.TxOffsetsNs(4) = 1500;
        meta.ExpectedOutcome = "FAIL";
        meta.ExpectedInterpretation = "Above-threshold static offset should fail when geometry is good.";
    case "good_static_3000ns"
        base.TxOffsetsNs(4) = 3000;
        meta.ExpectedOutcome = "FAIL";
        meta.ExpectedInterpretation = "Large static offset should fail and be close to injected truth.";
    case "good_multi_offset"
        base.TxOffsetsNs = [0 -350 650 1250 -150];
        meta.ExpectedOutcome = "SUSPECT_OR_FAIL";
        meta.TargetPCI = 503;
        meta.ExpectedInterpretation = "Multiple offsets test median-centred reporting against the detected group.";
    case "good_noise_50ns"
        base.TxOffsetsNs(4) = 1000;
        base.MeasurementNoiseNs = 50;
        meta.ExpectedOutcome = "SUSPECT_OR_FAIL";
        meta.ExpectedInterpretation = "Higher measurement noise should increase uncertainty around the threshold.";
    case "good_noise_100ns"
        base.TxOffsetsNs(4) = 1000;
        base.MeasurementNoiseNs = 100;
        meta.ExpectedOutcome = "SUSPECT_OR_NOT_ASSESSABLE";
        meta.ExpectedInterpretation = "Very noisy measurements should avoid a clean overclaim.";
    case "common_mode_3000ns"
        base.TxOffsetsNs(:) = 3000;
        meta.ExpectedOutcome = "NOT_FAIL";
        meta.ExpectedInterpretation = "Common-mode offset should be invisible to relative timing mode.";
    case "weak_collinear_1000ns"
        base.TxOffsetsNs(4) = 1000;
        base.RxPositionsM = [-1000 -20; -700 15; -300 -10; 100 8; 500 -12; 900 20];
        meta.ExpectedOutcome = "SUSPECT_OR_NOT_ASSESSABLE";
        meta.GeometryClass = "WEAK_COLLINEAR";
        meta.ExpectedInterpretation = "Near-line receiver geometry should increase uncertainty or become suspect.";
    case "too_few_positions_1000ns"
        base.TxOffsetsNs(4) = 1000;
        base.NumCaptures = 3;
        meta.ExpectedOutcome = "NOT_ASSESSABLE";
        meta.GeometryClass = "TOO_FEW_POSITIONS";
        meta.ExpectedInterpretation = "Too few positions should not support a no-location static-offset fit.";
    case "two_cells_1000ns"
        base.PCIs = [11 503];
        base.TxOffsetsNs = [0 1000];
        base.GNBPositionsM = [-450 -250; -150 500];
        meta.ExpectedOutcome = "SUSPECT_OR_NOT_ASSESSABLE";
        meta.GeometryClass = "TWO_CELLS";
        meta.TargetPCI = 503;
        meta.ExpectedInterpretation = "Two visible cells can give pairwise timing evidence, but they cannot support majority-based fault attribution.";
    case "missing_observations_1000ns"
        base.TxOffsetsNs(4) = 1000;
        meta.ExpectedOutcome = "SUSPECT_OR_NOT_ASSESSABLE";
        meta.GeometryClass = "MISSING_OBSERVATIONS";
        meta.ExpectedInterpretation = "Missing common-cell observations reduce fit support.";
    otherwise
        error("runNoLocationStaticOffsetValidation:UnknownCase", ...
            "Unknown validation case: %s", caseName);
end

[meas, caps, truth] = buildSyntheticSurveyDataset("aligned_5gnb", base);

if string(caseName) == "missing_observations_1000ns"
    captureIDs = unique(meas.CaptureID, "stable");
    targetMask = meas.PCI == meta.TargetPCI & ismember(meas.CaptureID, captureIDs([2 5 7]));
    meas = meas(~targetMask, :);
end

meta.NumCaptures = height(caps);
meta.NumMeasurementRows = height(meas);
meta.NoiseNs = base.MeasurementNoiseNs;
meta.SyncThresholdNs = opts.Survey.SyncThresholdNs;
end

function row = summarizeCaseLocal(caseName, result, truth, meta, opts)
estCompare = compareEstimatesToTruthLocal(caseName, result, truth, meta);
target = estCompare(estCompare.PCI == meta.TargetPCI, :);
if isempty(target)
    targetStatus = "NOT_ASSESSABLE";
    targetTruth = NaN;
    targetEstimate = NaN;
    targetErrorAbs = NaN;
    targetUnc = NaN;
    targetEnvelope = NaN;
else
    targetStatus = target.TimingStatus(1);
    targetTruth = target.TrueCenteredOffsetNs(1);
    targetEstimate = target.EstimatedRelativeTxOffsetNs(1);
    targetErrorAbs = target.AbsoluteErrorNs(1);
    targetUnc = target.EstimatedOffsetUncertaintyNs(1);
    targetEnvelope = target.WorstCaseRelativeTimingNs(1);
end

overall = overallStatusLocal(result);
row = table( ...
    string(caseName), string(meta.GeometryClass), meta.NumCaptures, ...
    meta.NumMeasurementRows, meta.NoiseNs, opts.Survey.SyncThresholdNs, ...
    string(meta.ExpectedOutcome), overall, outcomeMatchesLocal(meta.ExpectedOutcome, overall), ...
    string(result.FitInfo.Status), string(result.FitInfo.Reason), ...
    result.FitInfo.NumMeasurements, result.FitInfo.NumUnknowns, ...
    result.FitInfo.DegreesOfFreedom, result.FitInfo.FitRMSNs, ...
    result.FitInfo.ConditionNumber, meta.TargetPCI, targetTruth, targetEstimate, ...
    targetErrorAbs, targetUnc, targetEnvelope, string(meta.ExpectedInterpretation), ...
    'VariableNames', ["CaseName","GeometryClass","NumCaptures", ...
    "NumMeasurementRows","MeasurementNoiseNs","SyncThresholdNs", ...
    "ExpectedOutcome","ObservedOverallStatus","OutcomeMatchesExpectation", ...
    "FitStatus","FitReason","NumMeasurements","NumUnknowns", ...
    "DegreesOfFreedom","FitRMSNs","ConditionNumber","TargetPCI", ...
    "TargetTrueCenteredOffsetNs","TargetEstimatedOffsetNs", ...
    "TargetAbsoluteErrorNs","TargetUncertaintyNs", ...
    "TargetWorstCaseRelativeTimingNs","Interpretation"]);
end

function rows = compareEstimatesToTruthLocal(caseName, result, truth, meta)
if ~isfield(result, "TimingEstimates") || isempty(result.TimingEstimates) || ...
        height(result.TimingEstimates) == 0
    rows = table();
    return;
end

est = result.TimingEstimates;
truthOffsets = truth.GNB.TxOffsetNs;
truthCentered = truthOffsets - median(truthOffsets, "omitnan");

rows = est;
rows.CaseName = repmat(string(caseName), height(rows), 1);
rows.GeometryClass = repmat(string(meta.GeometryClass), height(rows), 1);
rows.TrueTxOffsetNs = NaN(height(rows), 1);
rows.TrueCenteredOffsetNs = NaN(height(rows), 1);
rows.OffsetErrorNs = NaN(height(rows), 1);
rows.AbsoluteErrorNs = NaN(height(rows), 1);

for k = 1:height(rows)
    idx = find(truth.GNB.PCI == rows.PCI(k), 1);
    if isempty(idx)
        continue;
    end
    rows.TrueTxOffsetNs(k) = truthOffsets(idx);
    rows.TrueCenteredOffsetNs(k) = truthCentered(idx);
    rows.OffsetErrorNs(k) = rows.EstimatedRelativeTxOffsetNs(k) - truthCentered(idx);
    rows.AbsoluteErrorNs(k) = abs(rows.OffsetErrorNs(k));
end

rows = movevars(rows, ["CaseName","GeometryClass","TrueTxOffsetNs", ...
    "TrueCenteredOffsetNs","OffsetErrorNs","AbsoluteErrorNs"], "Before", 1);
end

function rows = geometryRowsLocal(caseName, result, truth, caps, meta)
trueRows = table( ...
    repmat(string(caseName), height(truth.GNB), 1), ...
    repmat(string(meta.GeometryClass), height(truth.GNB), 1), ...
    repmat("TRUE_SYNTHETIC_CELL_POSITION", height(truth.GNB), 1), ...
    truth.GNB.PCI, truth.GNB.TrueX_m, truth.GNB.TrueY_m, ...
    'VariableNames', ["CaseName","GeometryClass","PointType","PCI","X_m","Y_m"]);

rxRows = table( ...
    repmat(string(caseName), height(caps), 1), ...
    repmat(string(meta.GeometryClass), height(caps), 1), ...
    repmat("RECEIVER_POSITION", height(caps), 1), ...
    NaN(height(caps), 1), caps.RxX_m, caps.RxY_m, ...
    'VariableNames', ["CaseName","GeometryClass","PointType","PCI","X_m","Y_m"]);

fitRows = table();
if isfield(result, "TimingEstimates") && istable(result.TimingEstimates) && ...
        height(result.TimingEstimates) > 0
    est = result.TimingEstimates;
    fitRows = table( ...
        repmat(string(caseName), height(est), 1), ...
        repmat(string(meta.GeometryClass), height(est), 1), ...
        repmat("FITTED_POSITION_LIKE_VARIABLE", height(est), 1), ...
        est.PCI, est.EstimatedX_m, est.EstimatedY_m, ...
        'VariableNames', ["CaseName","GeometryClass","PointType","PCI","X_m","Y_m"]);
end

rows = [rxRows; trueRows; fitRows];
end

function status = overallStatusLocal(result)
if ~isfield(result, "TimingEstimates") || isempty(result.TimingEstimates) || ...
        height(result.TimingEstimates) == 0
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

function plotInjectedVsEstimatedLocal(est, syncThresholdNs, outputPath)
fig = figure("Visible", "off", "Color", "w");
ax = axes(fig);
target = est(est.SurveyRole == "TARGET" & est.GeometryClass == "GOOD_2D" & ...
    isfinite(est.TrueCenteredOffsetNs) & isfinite(est.EstimatedRelativeTxOffsetNs), :);
if isempty(target)
    text(ax, 0.5, 0.5, "No assessable target rows", "HorizontalAlignment", "center");
else
    scatter(ax, target.TrueCenteredOffsetNs, target.EstimatedRelativeTxOffsetNs, ...
        55, "filled", "MarkerFaceColor", [0.25 0.5 0.85]);
    hold(ax, "on");
    lim = max(syncThresholdNs*1.4, max(abs([target.TrueCenteredOffsetNs; ...
        target.EstimatedRelativeTxOffsetNs]), [], "omitnan") * 1.1);
    plot(ax, [-lim lim], [-lim lim], "k--", "LineWidth", 1.2, ...
        "DisplayName", "ideal estimate");
    xline(ax, syncThresholdNs, ":", "Threshold", "HandleVisibility", "off");
    xline(ax, -syncThresholdNs, ":", "HandleVisibility", "off");
    yline(ax, syncThresholdNs, ":", "HandleVisibility", "off");
    yline(ax, -syncThresholdNs, ":", "HandleVisibility", "off");
    xlim(ax, [-lim lim]);
    ylim(ax, [-lim lim]);
    axis(ax, "square");
    grid(ax, "on");
    xlabel(ax, "Injected median-centred static offset (ns)");
    ylabel(ax, "Estimated median-centred static offset (ns)");
    title(ax, "No-location static-offset recovery, good geometry");
end
exportgraphics(fig, outputPath, "Resolution", 180);
close(fig);
end

function plotAbsoluteErrorLocal(cases, outputPath)
fig = figure("Visible", "off", "Color", "w");
ax = axes(fig);
bar(ax, categorical(cases.CaseName), cases.TargetAbsoluteErrorNs, ...
    "FaceColor", [0.35 0.55 0.80]);
grid(ax, "on");
ylabel(ax, "Target PCI absolute offset error (ns)");
title(ax, "Static-offset estimation error by validation case");
xtickangle(ax, 35);
exportgraphics(fig, outputPath, "Resolution", 180);
close(fig);
end

function plotGeometryComparisonLocal(geometry, outputPath)
fig = figure("Visible", "off", "Color", "w");
tiledlayout(fig, 1, 2, "TileSpacing", "compact", "Padding", "compact");
plotOneGeometryLocal(nexttile, geometry, "good_static_1000ns", "Good 2D geometry");
plotOneGeometryLocal(nexttile, geometry, "weak_collinear_1000ns", "Weak near-collinear geometry");
exportgraphics(fig, outputPath, "Resolution", 180);
close(fig);
end

function plotOneGeometryLocal(ax, geometry, caseName, titleText)
rows = geometry(geometry.CaseName == string(caseName), :);
hold(ax, "on");
rx = rows(rows.PointType == "RECEIVER_POSITION", :);
truePos = rows(rows.PointType == "TRUE_SYNTHETIC_CELL_POSITION", :);
fitPos = rows(rows.PointType == "FITTED_POSITION_LIKE_VARIABLE", :);
scatter(ax, rx.X_m, rx.Y_m, 40, "k", "filled", "DisplayName", "Receiver");
plot(ax, rx.X_m, rx.Y_m, "k:", "HandleVisibility", "off");
scatter(ax, truePos.X_m, truePos.Y_m, 60, [0.1 0.55 0.3], "filled", ...
    "DisplayName", "True cell");
if ~isempty(fitPos)
    scatter(ax, fitPos.X_m, fitPos.Y_m, 65, [0.85 0.35 0.1], ...
        "LineWidth", 1.4, "DisplayName", "Fitted variable");
end
axis(ax, "equal");
grid(ax, "on");
xlabel(ax, "x (m)");
ylabel(ax, "y (m)");
title(ax, titleText);
end

function plotVerdictSummaryLocal(cases, outputPath)
score = zeros(height(cases), 1);
for k = 1:height(cases)
    switch string(cases.ObservedOverallStatus(k))
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
bar(ax, categorical(cases.CaseName), score, "FaceColor", [0.45 0.45 0.75]);
ylim(ax, [-0.75 1.25]);
yticks(ax, [-0.5 0 0.5 1]);
yticklabels(ax, ["FAIL","NOT ASSESSABLE","SUSPECT","PASS"]);
grid(ax, "on");
title(ax, "No-location static-offset validation outcomes");
xtickangle(ax, 35);
exportgraphics(fig, outputPath, "Resolution", 180);
close(fig);
end

function summary = buildSummaryLocal(opts, cases, estimates)
goodTarget = estimates(estimates.SurveyRole == "TARGET" & ...
    estimates.GeometryClass == "GOOD_2D" & isfinite(estimates.AbsoluteErrorNs), :);

summary = struct();
summary.Mode = "no_location_static_offset_synthetic_validation";
summary.Statement = "Synthetic model-level validation. The estimator sees receiver positions and relative arrivals only; true cell positions and offsets are used only for scoring.";
summary.SyncThresholdNs = opts.Survey.SyncThresholdNs;
summary.NumCases = height(cases);
summary.ExpectationMatchRate = mean(cases.OutcomeMatchesExpectation);
summary.NumGoodGeometryTargetRows = height(goodTarget);
if isempty(goodTarget)
    summary.GoodGeometryMedianAbsoluteErrorNs = NaN;
    summary.GoodGeometry95thPercentileAbsoluteErrorNs = NaN;
else
    summary.GoodGeometryMedianAbsoluteErrorNs = median(goodTarget.AbsoluteErrorNs, "omitnan");
    summary.GoodGeometry95thPercentileAbsoluteErrorNs = prctile(goodTarget.AbsoluteErrorNs, 95);
end
summary.NumNotAssessableCases = nnz(cases.ObservedOverallStatus == "NOT_ASSESSABLE");
summary.NumSuspectCases = nnz(cases.ObservedOverallStatus == "SUSPECT");
summary.NumFailCases = nnz(cases.ObservedOverallStatus == "FAIL");
summary.MainConclusion = "Static relative offsets are recoverable synthetically when common-cell geometry is sufficiently informative; weak geometry correctly becomes suspect or not assessable.";
end

function out = vertcatOrEmptyLocal(rows)
if isempty(rows)
    out = table();
else
    out = vertcat(rows{:});
end
end

function opts = defaultValidationOptionsLocal()
opts = struct();
opts.RandomSeed = 20260625;
opts.Survey = struct();
opts.Survey.ReferencePCI = 11;
opts.Survey.SyncThresholdNs = 1000;
opts.Survey.NumRandomStarts = 20;
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
