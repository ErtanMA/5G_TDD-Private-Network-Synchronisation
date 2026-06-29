function limits = analyzeSystemLimits(outputPath, opts)
%analyzeSystemLimits Summarize timing and survey observability limits.
%
%   limits = analyzeSystemLimits() returns tables that explain when the tool
%   passes, suspects, fails, or cannot assess a result.
%
%   limits = analyzeSystemLimits("reports/system_limits_report.md") also
%   writes a Markdown report. This is a deterministic theory/reporting helper:
%   it does not rerun the full IQ receiver.

arguments
    outputPath (1,1) string = ""
    opts.SigmaNs (:,1) double = [10; 20; 50; 100; 150; 300]
    opts.OffsetGridNs (:,1) double = [0; 50; 100; 150; 200; 250; 300; 500; 750; 1000; 1500; 2000; 3000]
    opts.GNBCounts (:,1) double = (1:10).'
    opts.SyncThresholdNs (1,1) double = NaN
    opts.WarningOffsetNs (1,1) double = NaN
    opts.FailOffsetNs (1,1) double = 1000
    opts.HardFailOffsetNs (1,1) double = Inf
    opts.UncertaintySigmaMultiplier (1,1) double = 3
    opts.BiasUncertaintyPassFraction (1,1) double = 0.20
    opts.MaxBiasUncertaintyForPassNs (1,1) double = NaN
    opts.DefaultMinCaptures (1,1) double = 5
    opts.MinMeasurementsPerPCI (1,1) double = 4
end

opts = normalizeLimitOptionsLocal(opts);
thresholdTable = buildThresholdTableLocal(opts);
statusTable = buildStatusGridLocal(opts);
gnbTable = buildGNBObservabilityTableLocal(opts);

limits = struct();
limits.Mode = "relative_unknown_gnb_location_survey";
limits.Statement = "Limits are for internal-clock B210 relative survey mode, not absolute UTC compliance.";
limits.DecisionThresholds = thresholdTable;
limits.OffsetStatusGrid = statusTable;
limits.GNBObservability = gnbTable;
limits.TheoreticalLimits = theoreticalLimitsLocal();

if strlength(outputPath) > 0
    writeLimitReportLocal(outputPath, limits, opts);
    limits.ReportPath = outputPath;
end

end

function tbl = buildThresholdTableLocal(opts)
sigma = opts.SigmaNs(:);
k = opts.UncertaintySigmaMultiplier;

if isfinite(opts.WarningOffsetNs)
    cleanPassUpTo = max(0, opts.WarningOffsetNs - k*sigma);
    definiteWarningFrom = opts.WarningOffsetNs + k*sigma;
else
    cleanPassUpTo = max(0, opts.FailOffsetNs - k*sigma);
    definiteWarningFrom = max(0, opts.FailOffsetNs - k*sigma);
end
failFrom = opts.FailOffsetNs + k*sigma;
if isfinite(opts.HardFailOffsetNs)
    failFrom = min(failFrom, opts.HardFailOffsetNs * ones(size(sigma)));
end
hardFailFrom = opts.HardFailOffsetNs * ones(size(sigma));

if isfinite(opts.WarningOffsetNs)
    passBoundary = opts.WarningOffsetNs;
else
    passBoundary = opts.FailOffsetNs;
end
zeroOffsetCanPass = (k*sigma <= passBoundary) & ...
    (sigma <= opts.MaxBiasUncertaintyForPassNs);

tbl = table(sigma, cleanPassUpTo, definiteWarningFrom, failFrom, ...
    hardFailFrom, zeroOffsetCanPass, ...
    'VariableNames', ["OffsetUncertaintySigmaNs","CleanPassUpToNs", ...
    "DefiniteWarningFromNs","FailFromNs","HardFailFromNs", ...
    "ZeroOffsetCanPass"]);
end

function tbl = buildStatusGridLocal(opts)
rows = {};
rowCount = 0;
for s = 1:numel(opts.SigmaNs)
    sigma = opts.SigmaNs(s);
    for o = 1:numel(opts.OffsetGridNs)
        offset = opts.OffsetGridNs(o);
        [status, reason] = classifyOffsetLocal(offset, sigma, opts);
        rowCount = rowCount + 1;
        rows{rowCount,1} = table(sigma, offset, status, reason, ...
            'VariableNames', ["OffsetUncertaintySigmaNs","OffsetNs", ...
            "Status","Reason"]);
    end
end
tbl = vertcat(rows{:});
end

function tbl = buildGNBObservabilityTableLocal(opts)
rows = cell(numel(opts.GNBCounts), 1);
for i = 1:numel(opts.GNBCounts)
    n = opts.GNBCounts(i);
    if n < 2
        unknowns = NaN;
        minCaptures = Inf;
        status = "IMPOSSIBLE";
        reason = "At least two visible gNBs are required for relative synchronization.";
    else
        unknowns = 3*n - 1;             % 2D position per gNB + N-1 timing biases
        minCaptures = floor(unknowns / (n-1)) + 1; % require observations > unknowns
        minCaptures = max(minCaptures, opts.MinMeasurementsPerPCI);
        if minCaptures <= opts.DefaultMinCaptures
            status = "SUPPORTED_BY_DEFAULT";
            reason = "Default five-position survey is mathematically overdetermined, assuming all gNBs are visible with good geometry.";
        else
            status = "NEEDS_MORE_CAPTURES";
            reason = "More receiver positions are needed before the nonlinear survey fit is overdetermined.";
        end
    end

    recommended = minCaptures;
    if isfinite(recommended)
        recommended = max(recommended + 2, opts.DefaultMinCaptures);
    end

    rows{i} = table(n, unknowns, minCaptures, recommended, status, reason, ...
        'VariableNames', ["VisibleGNBs","UnknownParameters2D", ...
        "MinimumReceiverPositions","RecommendedReceiverPositions", ...
        "Assessability","Reason"]);
end
tbl = vertcat(rows{:});
end

function [status, reason] = classifyOffsetLocal(offsetNs, sigmaNs, opts)
margin = opts.UncertaintySigmaMultiplier * sigmaNs;

if isfinite(opts.HardFailOffsetNs) && offsetNs >= opts.HardFailOffsetNs
    status = "FAIL";
    reason = "Offset exceeds hard-fail threshold regardless of uncertainty.";
elseif offsetNs - margin >= opts.FailOffsetNs
    status = "FAIL";
    reason = "Offset exceeds synchronization threshold after uncertainty margin.";
elseif offsetNs + margin >= opts.FailOffsetNs
    status = "SUSPECT";
    reason = "Offset is near the synchronization threshold after uncertainty handling.";
elseif isfinite(opts.WarningOffsetNs) && offsetNs - margin >= opts.WarningOffsetNs
    status = "SUSPECT";
    reason = "Offset exceeds warning threshold after uncertainty margin.";
elseif isfinite(opts.WarningOffsetNs) && offsetNs + margin >= opts.WarningOffsetNs
    status = "SUSPECT";
    reason = "Offset is near the warning threshold after uncertainty handling.";
elseif sigmaNs <= opts.MaxBiasUncertaintyForPassNs
    status = "PASS";
    reason = "Offset is below the synchronization threshold and uncertainty is acceptable.";
else
    status = "SUSPECT";
    reason = "Offset is near a threshold or uncertainty is too high for a clean pass.";
end
end

function limits = theoreticalLimitsLocal()
limits = [
    "Absolute/common-mode timing cannot be detected without external timing such as GPSDO or 10 MHz + 1 PPS."
    "A single visible gNB cannot be checked for relative synchronization."
    "If receiver positions are unknown, static transmit timing and propagation delay are inseparable."
    "Too few receiver positions, collinear survey geometry, or missing common reference PCI can make the survey fit not assessable."
    "NLOS multipath can bias arrival timing; the tool reports confidence but cannot know the true path without more information."
    "A gNB outside the captured bandwidth or below detection SNR is invisible to the method."
    "PCI reuse or unresolved overlapping SSBs can confuse association; SIC reduces but does not eliminate this risk."
    "The method assumes gNBs are stationary and timing offsets are constant during the survey."
    ];
end

function writeLimitReportLocal(outputPath, limits, opts)
folder = fileparts(outputPath);
if strlength(folder) > 0 && ~isfolder(folder)
    mkdir(folder);
end

fid = fopen(outputPath, "w");
if fid < 0
    error("analyzeSystemLimits:OpenFailed", ...
        "Could not write report: %s", outputPath);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, "# System Limits Report\n\n");
fprintf(fid, "Mode: `%s`\n\n", limits.Mode);
fprintf(fid, "%s\n\n", limits.Statement);

fprintf(fid, "## Decision Thresholds\n\n");
fprintf(fid, "The implemented timing decision uses these defaults:\n\n");
fprintf(fid, "- Static relative timing desynchronization threshold: %s ns\n", formatThresholdLocal(opts.FailOffsetNs));
if isfinite(opts.WarningOffsetNs)
    fprintf(fid, "- Static warning threshold: %s ns\n", formatThresholdLocal(opts.WarningOffsetNs));
else
    fprintf(fid, "- Static warning mode: uncertainty interval overlaps the desynchronization threshold\n");
end
fprintf(fid, "- Static hard-fail threshold: %s ns\n", formatThresholdLocal(opts.HardFailOffsetNs));
fprintf(fid, "- Uncertainty margin: %.0f sigma\n\n", opts.UncertaintySigmaMultiplier);
fprintf(fid, "%s\n\n", tableToMarkdownLocal(limits.DecisionThresholds));

fprintf(fid, "Interpretation: with a good survey uncertainty of about 20 ns, a clean PASS is possible up to about %.0f ns and a definite FAIL starts around %.0f ns.\n", ...
    limits.DecisionThresholds.CleanPassUpToNs(limits.DecisionThresholds.OffsetUncertaintySigmaNs == 20), ...
    limits.DecisionThresholds.FailFromNs(limits.DecisionThresholds.OffsetUncertaintySigmaNs == 20));
if isfinite(opts.HardFailOffsetNs)
    fprintf(fid, "Any %.0f ns offset hard-fails regardless of uncertainty.\n\n", opts.HardFailOffsetNs);
else
    fprintf(fid, "No separate hard-fail threshold is enabled in this threshold profile.\n\n");
end

fprintf(fid, "## Offset Status Grid\n\n");
fprintf(fid, "%s\n\n", tableToMarkdownLocal(limits.OffsetStatusGrid));

fprintf(fid, "## Number of Base Stations and Receiver Positions\n\n");
fprintf(fid, "For `N` visible gNBs, the 2D no-location survey estimates `2N` position parameters plus `N-1` relative timing biases. Each receiver position contributes `N-1` relative timing observations. The fit should be overdetermined, so `M(N-1) > 3N-1`.\n\n");
fprintf(fid, "%s\n\n", tableToMarkdownLocal(limits.GNBObservability));

fprintf(fid, "## When It Theoretically Does Not Work\n\n");
for k = 1:numel(limits.TheoreticalLimits)
    fprintf(fid, "- %s\n", limits.TheoreticalLimits(k));
end
fprintf(fid, "\n");

fprintf(fid, "## Practical Recommendation\n\n");
fprintf(fid, "Use at least three visible gNBs when possible, and collect at least six to eight receiver positions around the site rather than along one straight line. Two visible gNBs can be mathematically assessable only with at least six good receiver positions, but it is much less robust than three or more.\n");
end

function opts = normalizeLimitOptionsLocal(opts)
if isfinite(opts.SyncThresholdNs)
    opts.FailOffsetNs = opts.SyncThresholdNs;
    opts.WarningOffsetNs = NaN;
    opts.HardFailOffsetNs = Inf;
end
if ~isfinite(opts.MaxBiasUncertaintyForPassNs)
    opts.MaxBiasUncertaintyForPassNs = opts.BiasUncertaintyPassFraction * opts.FailOffsetNs;
end
end

function md = tableToMarkdownLocal(tbl)
names = string(tbl.Properties.VariableNames);
lines = strings(0,1);
lines(end+1,1) = "| " + strjoin(names, " | ") + " |";
lines(end+1,1) = "| " + strjoin(repmat("---", 1, numel(names)), " | ") + " |";
for r = 1:height(tbl)
    vals = strings(1, numel(names));
    for c = 1:numel(names)
        value = tbl.(names(c))(r);
        if isstring(value)
            vals(c) = value;
        elseif iscell(value)
            vals(c) = string(value{1});
        elseif islogical(value)
            vals(c) = string(value);
        elseif isnumeric(value)
            if isinf(value)
                vals(c) = "Inf";
            elseif isnan(value)
                vals(c) = "NaN";
            else
                vals(c) = string(round(value, 3));
            end
        else
            vals(c) = string(value);
        end
    end
    lines(end+1,1) = "| " + strjoin(vals, " | ") + " |"; %#ok<AGROW>
end
md = strjoin(lines, newline);
end

function txt = formatThresholdLocal(value)
if isnan(value) || isinf(value)
    txt = "disabled";
else
    txt = string(round(value, 3));
end
end
