function artifacts = generateThesisArtifacts(outputDir, opts)
%generateThesisArtifacts Create thesis-ready synthetic reports and figures.
%
%   artifacts = generateThesisArtifacts(outputDir) writes deterministic
%   pre-real-data artifacts: threshold tables, synthetic survey reports, and
%   figures that explain the method and its limits. These figures are
%   synthetic-validation artifacts, not real measurement results.

arguments
    outputDir (1,1) string = "reports/thesis_artifacts"
    opts struct = struct()
end

opts = mergeStructsLocal(defaultArtifactOptionsLocal(), opts);
if ~isfolder(outputDir)
    mkdir(outputDir);
end
figureDir = fullfile(outputDir, "figures");
if ~isfolder(figureDir)
    mkdir(figureDir);
end

limits = analyzeSystemLimits(fullfile(outputDir, "system_limits_report.md"));
thresholdFig = fullfile(figureDir, "timing_thresholds.png");
gnbFig = fullfile(figureDir, "receiver_positions_vs_visible_gnbs.png");
plotTimingThresholdsLocal(limits.DecisionThresholds, thresholdFig);
plotGNBRequirementsLocal(limits.GNBObservability, gnbFig);

[measurementTable, captureInfo, truth] = buildSyntheticSurveyDataset( ...
    opts.SyntheticScenario, opts.SyntheticSurvey);
surveyResult = noLocationSurveyCheck(measurementTable, captureInfo, opts.Survey);
surveyResult.IsSynthetic = true;
surveyResult.ArtifactType = "SYNTHETIC_VALIDATION";
surveyReportPaths = generateSurveyReport(surveyResult, outputDir, ...
    "synthetic_" + opts.SyntheticScenario + "_survey");

surveyFig = fullfile(figureDir, "synthetic_survey_geometry_and_offsets.png");
offsetFig = fullfile(figureDir, "synthetic_relative_timing_offsets.png");
plotSyntheticSurveyLocal(surveyResult, truth, surveyFig);
plotTimingOffsetsLocal(surveyResult, offsetFig);

powerTable = syntheticRssPowerTableLocal(captureInfo, truth, opts.RSSSynthetic);
rssPlanning = estimateRssPriorsAndPlanStops("all", powerTable, captureInfo, table(), opts.RSS);
rssFig = fullfile(figureDir, "rss_assisted_receiver_stop_planning.png");
plotRssPlanningLocal(powerTable, truth, rssPlanning, rssFig);
thesisSyntheticResults = runThesisSyntheticResults(outputDir, opts.ThesisSynthetic);
noLocationStaticOffsetValidation = runNoLocationStaticOffsetValidation( ...
    fullfile(outputDir, "no_location_static_offset_validation"), ...
    opts.NoLocationStaticOffsetValidation);

indexPath = fullfile(outputDir, "thesis_artifact_index.md");
writeArtifactIndexLocal(indexPath, opts, surveyResult, rssPlanning, limits, ...
    thesisSyntheticResults, noLocationStaticOffsetValidation);

artifacts = struct();
artifacts.Mode = "synthetic_thesis_artifact_package";
artifacts.Statement = "Artifacts are synthetic and support thesis explanation before real B210 data is available.";
artifacts.OutputDir = outputDir;
artifacts.FigureDir = figureDir;
artifacts.IndexMarkdown = indexPath;
artifacts.SystemLimitsReport = fullfile(outputDir, "system_limits_report.md");
artifacts.SurveyReportPaths = surveyReportPaths;
artifacts.Figures = struct( ...
    "TimingThresholds", thresholdFig, ...
    "GNBRequirements", gnbFig, ...
    "SyntheticSurvey", surveyFig, ...
    "TimingOffsets", offsetFig, ...
    "RSSPlanning", rssFig);
artifacts.SurveyResult = surveyResult;
artifacts.RSSPlanning = rssPlanning;
artifacts.ThesisSyntheticResults = thesisSyntheticResults;
artifacts.NoLocationStaticOffsetValidation = noLocationStaticOffsetValidation;

end

function opts = defaultArtifactOptionsLocal()
opts = struct();
opts.SyntheticScenario = "offset_3us_one";
opts.SyntheticSurvey = struct("RandomSeed", 20260604, ...
    "NumCaptures", 8, "MeasurementNoiseNs", 20);
opts.Survey = struct("NumRandomStarts", 16, "MinCaptures", 5);
opts.RSS = struct("GridSpacingM", 75, "SearchMarginM", 900, ...
    "MinMeasurementsPerPCI", 4, "VisibilityPowerDbFS", -95, ...
    "NoiseFloorDbFS", -110, "MaxSuggestedStops", 8);
opts.RSSSynthetic = struct("PathLossExponent", 2.4, ...
    "ReferencePowerDbFS", -8, "NoiseFloorDbFS", -110, ...
    "ShadowingDb", 2.5);
opts.ThesisSynthetic = struct();
opts.NoLocationStaticOffsetValidation = struct();
end

function plotTimingThresholdsLocal(thresholds, outputPath)
fig = figure("Visible", "off", "Color", "w");
ax = axes(fig);
hold(ax, "on");
plot(ax, thresholds.OffsetUncertaintySigmaNs, thresholds.CleanPassUpToNs, ...
    "-o", "LineWidth", 1.7, "DisplayName", "Clean PASS up to");
plot(ax, thresholds.OffsetUncertaintySigmaNs, thresholds.FailFromNs, ...
    "-s", "LineWidth", 1.7, "DisplayName", "FAIL from");
yline(ax, thresholds.HardFailFromNs(1), "--", "Hard FAIL", ...
    "LineWidth", 1.3, "LabelHorizontalAlignment", "left");
grid(ax, "on");
xlabel(ax, "Offset uncertainty sigma (ns)");
ylabel(ax, "Relative timing offset (ns)");
title(ax, "Decision thresholds vs timing uncertainty");
legend(ax, "Location", "northwest");
exportgraphics(fig, outputPath, "Resolution", 180);
close(fig);
end

function plotGNBRequirementsLocal(tbl, outputPath)
fig = figure("Visible", "off", "Color", "w");
ax = axes(fig);
valid = isfinite(tbl.MinimumReceiverPositions);
hold(ax, "on");
bar(ax, tbl.VisibleGNBs(valid), tbl.MinimumReceiverPositions(valid), ...
    "FaceColor", [0.25 0.45 0.85], "DisplayName", "Minimum");
plot(ax, tbl.VisibleGNBs(valid), tbl.RecommendedReceiverPositions(valid), ...
    "-o", "Color", [0.8 0.25 0.2], "LineWidth", 1.7, ...
    "DisplayName", "Recommended");
grid(ax, "on");
xlabel(ax, "Visible gNBs");
ylabel(ax, "Receiver positions");
title(ax, "Receiver positions needed for no-location survey fitting");
legend(ax, "Location", "northeast");
exportgraphics(fig, outputPath, "Resolution", 180);
close(fig);
end

function plotSyntheticSurveyLocal(surveyResult, truth, outputPath)
fig = figure("Visible", "off", "Color", "w");
ax = axes(fig);
hold(ax, "on");
cap = surveyResult.CaptureInfo;
scatter(ax, cap.RxX_m, cap.RxY_m, 70, "k", "filled", ...
    "DisplayName", "Receiver stops");
plot(ax, cap.RxX_m, cap.RxY_m, "k:", "HandleVisibility", "off");
scatter(ax, truth.GNB.TrueX_m, truth.GNB.TrueY_m, 90, ...
    [0.15 0.55 0.35], "filled", "DisplayName", "True synthetic gNB");
est = surveyResult.TimingEstimates;
if istable(est) && height(est) > 0 && ismember("EstimatedX_m", string(est.Properties.VariableNames))
    scatter(ax, est.EstimatedX_m, est.EstimatedY_m, 90, ...
        [0.85 0.35 0.1], "LineWidth", 1.5, ...
        "DisplayName", "Fitted nuisance position");
    for k = 1:height(est)
        text(ax, est.EstimatedX_m(k), est.EstimatedY_m(k), ...
            " PCI " + string(est.PCI(k)), "FontSize", 8);
    end
end
axis(ax, "equal");
grid(ax, "on");
xlabel(ax, "x position (m)");
ylabel(ax, "y position (m)");
title(ax, "Synthetic survey geometry: fitted nuisance positions, not verified gNB sites");
legend(ax, "Location", "bestoutside");
exportgraphics(fig, outputPath, "Resolution", 180);
close(fig);
end

function plotTimingOffsetsLocal(surveyResult, outputPath)
fig = figure("Visible", "off", "Color", "w");
ax = axes(fig);
est = surveyResult.TimingEstimates;
if isempty(est) || height(est) == 0
    text(ax, 0.5, 0.5, "No timing estimates", "HorizontalAlignment", "center");
else
    offsets = est.EstimatedRelativeTxOffsetNs;
    unc = est.EstimatedOffsetUncertaintyNs;
    bar(ax, categorical(string(est.PCI)), offsets, "FaceColor", [0.3 0.55 0.85]);
    hold(ax, "on");
    errorbar(ax, categorical(string(est.PCI)), offsets, 3*unc, ...
        "k.", "LineWidth", 1.2, "DisplayName", "3 sigma");
    yline(ax, 250, "--", "Warning +250 ns", "LineWidth", 1.1);
    yline(ax, -250, "--", "Warning -250 ns", "LineWidth", 1.1);
    yline(ax, 1000, "-.", "Fail +1000 ns", "LineWidth", 1.1);
    yline(ax, -1000, "-.", "Fail -1000 ns", "LineWidth", 1.1);
    yline(ax, 3000, ":", "Hard +3000 ns", "LineWidth", 1.1);
    yline(ax, -3000, ":", "Hard -3000 ns", "LineWidth", 1.1);
    grid(ax, "on");
    xlabel(ax, "PCI");
    ylabel(ax, "Estimated relative offset (ns)");
    title(ax, "Synthetic relative timing offsets with uncertainty");
end
exportgraphics(fig, outputPath, "Resolution", 180);
close(fig);
end

function powerTable = syntheticRssPowerTableLocal(captureInfo, truth, opts)
rng(20260604, "twister");
rows = cell(height(captureInfo) * height(truth.GNB), 1);
rowCount = 0;
for m = 1:height(captureInfo)
    rx = [captureInfo.RxX_m(m), captureInfo.RxY_m(m)];
    for g = 1:height(truth.GNB)
        tx = [truth.GNB.TrueX_m(g), truth.GNB.TrueY_m(g)];
        d = max(norm(rx - tx), 10);
        powerDb = opts.ReferencePowerDbFS - 10*opts.PathLossExponent*log10(d) + ...
            opts.ShadowingDb * randn();
        rowCount = rowCount + 1;
        rows{rowCount} = table( ...
            captureInfo.CaptureID(m), truth.GNB.PCI(g), ...
            captureInfo.RxX_m(m), captureInfo.RxY_m(m), captureInfo.RxZ_m(m), ...
            powerDb, opts.NoiseFloorDbFS, powerDb - opts.NoiseFloorDbFS, 1, "OK", ...
            'VariableNames', ["CaptureID","PCI","RxX_m","RxY_m","RxZ_m", ...
            "PowerDbFS","NoiseDbFS","SNRdB","NumDetections","QualityFlag"]);
    end
end
powerTable = vertcat(rows{:});
end

function plotRssPlanningLocal(powerTable, truth, rssPlanning, outputPath)
fig = figure("Visible", "off", "Color", "w");
ax = axes(fig);
hold(ax, "on");
rx = unique(powerTable(:, ["RxX_m","RxY_m"]), "rows", "stable");
scatter(ax, rx.RxX_m, rx.RxY_m, 65, "k", "filled", ...
    "DisplayName", "Existing receiver stops");
scatter(ax, truth.GNB.TrueX_m, truth.GNB.TrueY_m, 90, ...
    [0.15 0.55 0.35], "filled", "DisplayName", "True synthetic gNB");

pri = rssPlanning.RssPriors;
if istable(pri) && height(pri) > 0
    scatter(ax, pri.EstimatedX_m, pri.EstimatedY_m, 90, ...
        [0.65 0.2 0.75], "LineWidth", 1.5, "DisplayName", "RSS rough prior");
end

stops = rssPlanning.SuggestedReceiverStops;
if istable(stops) && height(stops) > 0
    top = stops(1:min(5,height(stops)), :);
    scatter(ax, top.X_m, top.Y_m, 85, [0.95 0.55 0.05], "filled", ...
        "DisplayName", "Suggested next stops");
end
axis(ax, "equal");
grid(ax, "on");
xlabel(ax, "x position (m)");
ylabel(ax, "y position (m)");
title(ax, "RSS-assisted rough priors and suggested receiver stops");
legend(ax, "Location", "bestoutside");
exportgraphics(fig, outputPath, "Resolution", 180);
close(fig);
end

function writeArtifactIndexLocal(indexPath, opts, surveyResult, rssPlanning, limits, thesisSyntheticResults, noLocationStaticOffsetValidation)
fid = fopen(indexPath, "w");
if fid < 0
    error("generateThesisArtifacts:OpenFailed", ...
        "Could not write artifact index: %s", indexPath);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, "# Thesis Artifact Index\n\n");
fprintf(fid, "These artifacts are synthetic and intended for thesis explanation before real B210 data is available.\n\n");
fprintf(fid, "## Scenario\n\n");
fprintf(fid, "- Synthetic scenario: `%s`\n", opts.SyntheticScenario);
fprintf(fid, "- Survey fit status: `%s`\n", string(surveyResult.FitInfo.Status));
fprintf(fid, "- Survey fit reason: %s\n", string(surveyResult.FitInfo.Reason));
fprintf(fid, "- RSS priors produced: %d\n", height(rssPlanning.RssPriors));
fprintf(fid, "- Suggested receiver stops produced: %d\n\n", height(rssPlanning.SuggestedReceiverStops));

fprintf(fid, "## Included Files\n\n");
fprintf(fid, "- `system_limits_report.md`: threshold and observability limits.\n");
fprintf(fid, "- `synthetic_%s_survey_human_summary.md`: human-readable synthetic survey summary.\n", opts.SyntheticScenario);
fprintf(fid, "- `figures/timing_thresholds.png`: PASS/SUSPECT/FAIL threshold behavior.\n");
fprintf(fid, "- `figures/receiver_positions_vs_visible_gnbs.png`: receiver-position requirements.\n");
fprintf(fid, "- `figures/synthetic_survey_geometry_and_offsets.png`: synthetic survey geometry and fitted nuisance positions, not verified gNB sites.\n");
fprintf(fid, "- `figures/synthetic_relative_timing_offsets.png`: synthetic timing offsets with uncertainty.\n");
fprintf(fid, "- `figures/rss_assisted_receiver_stop_planning.png`: RSS rough priors and suggested receiver stops.\n\n");
fprintf(fid, "- `synthetic_results_summary.json`: thesis-ready synthetic result summary.\n");
fprintf(fid, "- `synthetic_results_survey_cases.csv`: deterministic survey-verdict cases.\n");
fprintf(fid, "- `synthetic_results_survey_estimates.csv`: per-PCI synthetic survey estimates.\n");
fprintf(fid, "- `synthetic_results_receiver_overview.csv`: compact receiver-chain overview.\n");
fprintf(fid, "- `figures/synthetic_threshold_profile.png`: configurable 1 us threshold profile.\n");
fprintf(fid, "- `figures/synthetic_static_offset_sweep.png`: static-offset sweep result.\n");
fprintf(fid, "- `figures/synthetic_worst_case_envelope.png`: combined static-plus-instability envelope.\n");
fprintf(fid, "- `figures/synthetic_jitter_residuals.png`: residual instability case.\n");
fprintf(fid, "- `figures/synthetic_common_mode_invisibility.png`: common-mode invisibility case.\n");
fprintf(fid, "- `figures/synthetic_assessability_cases.png`: not-assessable and geometry edge cases.\n\n");
fprintf(fid, "- `no_location_static_offset_validation/no_location_static_offset_validation_summary.json`: model-level no-location static-offset validation summary.\n");
fprintf(fid, "- `no_location_static_offset_validation/no_location_static_offset_validation_cases.csv`: static-offset recovery and assessability cases.\n");
fprintf(fid, "- `no_location_static_offset_validation/no_location_static_offset_validation_estimates.csv`: per-PCI comparison between hidden injected truth and fitted offset.\n");
fprintf(fid, "- `no_location_static_offset_validation/figures/no_location_injected_vs_estimated.png`: injected-vs-estimated static offset under good geometry.\n");
fprintf(fid, "- `no_location_static_offset_validation/figures/no_location_static_offset_error.png`: target offset error by case.\n");
fprintf(fid, "- `no_location_static_offset_validation/figures/no_location_geometry_comparison.png`: good versus weak receiver geometry.\n");
fprintf(fid, "- `no_location_static_offset_validation/figures/no_location_verdict_summary.png`: static-offset validation verdict outcomes.\n\n");

fprintf(fid, "## Thesis Synthetic Result Summary\n\n");
fprintf(fid, "- Synthetic sync threshold: %.0f ns\n", thesisSyntheticResults.Summary.SyncThresholdNs);
fprintf(fid, "- Survey cases: %d\n", thesisSyntheticResults.Summary.NumSurveyCases);
fprintf(fid, "- Survey expectation match rate: %.3f\n", thesisSyntheticResults.Summary.SurveyExpectationMatchRate);
fprintf(fid, "- Receiver overview pass rate: %.3f\n\n", thesisSyntheticResults.Summary.ReceiverOverviewPassRate);
fprintf(fid, "## No-Location Static-Offset Validation Summary\n\n");
fprintf(fid, "- Cases: %d\n", noLocationStaticOffsetValidation.Summary.NumCases);
fprintf(fid, "- Expectation match rate: %.3f\n", noLocationStaticOffsetValidation.Summary.ExpectationMatchRate);
fprintf(fid, "- Good-geometry median absolute error: %.3f ns\n", noLocationStaticOffsetValidation.Summary.GoodGeometryMedianAbsoluteErrorNs);
fprintf(fid, "- Good-geometry 95th percentile absolute error: %.3f ns\n\n", noLocationStaticOffsetValidation.Summary.GoodGeometry95thPercentileAbsoluteErrorNs);

fprintf(fid, "## Correct Use\n\n");
fprintf(fid, "Use these figures to explain and validate the method. Do not present them as real measurement results. The default synthetic threshold is 1 us for private-network engineering assessment; it is configurable and is not regulatory or UTC compliance.\n\n");

fprintf(fid, "## Key Limit Reminder\n\n");
for k = 1:numel(limits.TheoreticalLimits)
    fprintf(fid, "- %s\n", limits.TheoreticalLimits(k));
end
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
