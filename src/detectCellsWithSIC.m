function [detections, debug] = detectCellsWithSIC(iq, meta, opts)
%detectCellsWithSIC Iterative PSS/SSS detection with PSS+SSS cancellation.
%
%   [detections, debug] = detectCellsWithSIC(iq, meta, opts) repeatedly runs
%   PSS/SSS detection on the current residual, selects the strongest not-yet
%   cancelled PCI, subtracts reconstructed PSS+SSS bursts for that PCI, and
%   then searches the residual again.

arguments
    iq (:,1) {mustBeNumeric}
    meta struct
    opts struct = struct()
end

[meta, ~] = validateMetadata(meta);
base = defaultSICOptions(meta);
opts = mergeStructsLocal(base, opts);

if ~opts.Enable
    pss = estimatePSSCandidates(iq, meta, opts.PSS);
    [detections, sssDebug] = detectSSSAndPCI(iq, pss, meta, opts.SSS);
    detections = addSICColumnsLocal(detections, 1, false);
    debug = struct("Options", opts, "Stages", {stageCellLocal(pss, detections, sssDebug, 1)}, ...
        "Cancellations", emptyCancellationTableForDebugLocal(), "Residual", iq(:));
    return;
end

residual = iq(:);
cancelledPCIs = zeros(0,1);
selectedRows = {};
stageRows = {};
cancellationRows = {};

for iter = 1:opts.MaxIterations
    pss = estimatePSSCandidates(residual, meta, opts.PSS);
    [sss, sssDebug] = detectSSSAndPCI(residual, pss, meta, opts.SSS);
    sss = addSICColumnsLocal(sss, iter, false);
    stageRows{end+1,1} = struct( ...
        "Iteration", iter, ...
        "PSSCandidates", pss, ...
        "SSSDetections", sss, ...
        "SSSDebug", sssDebug); %#ok<AGROW>

    usable = sss(sss.IsUsable & sss.SSSMetric >= opts.MinSelectedSSSMetric, :);
    if isempty(usable) || height(usable) < opts.MinNewUsableDetections
        break;
    end

    if opts.StopWhenNoNewPCI
        usable = usable(~ismember(usable.PCI, cancelledPCIs), :);
    end
    if isempty(usable) || height(usable) == 0
        break;
    end

    usable = sortrows(usable, "SSSMetric", "descend");
    selectedPCI = usable.PCI(1);
    if opts.CancelAllDetectionsForSelectedPCI
        toCancel = usable(usable.PCI == selectedPCI, :);
    else
        toCancel = usable(1, :);
    end
    toCancel = clusterCancellationRowsLocal(toCancel, ...
        opts.CancellationMinDistanceSamples, opts.MaxCancellationRowsPerIteration);
    toCancel.WasCancelled = true(height(toCancel), 1);
    selectedRows{end+1,1} = toCancel; %#ok<AGROW>

    [residual, cancellationTable] = cancelDetectedSSBs(residual, toCancel, meta, opts);
    if ~isempty(cancellationTable) && height(cancellationTable) > 0
        cancellationTable.Iteration = repmat(iter, height(cancellationTable), 1);
        cancellationRows{end+1,1} = cancellationTable; %#ok<AGROW>
    end

    cancelledPCIs(end+1,1) = selectedPCI; %#ok<AGROW>
end

if isempty(selectedRows)
    detections = emptySSSTableLikeLocal();
else
    detections = vertcat(selectedRows{:});
    detections.GlobalDetectionID = (1:height(detections)).';
    detections = movevars(detections, "GlobalDetectionID", "Before", 1);
    detections = sortrows(detections, ["SICIteration","StartSample0"], ["ascend","ascend"]);
end

if isempty(cancellationRows)
    cancellations = emptyCancellationTableForDebugLocal();
else
    cancellations = vertcat(cancellationRows{:});
    cancellations = movevars(cancellations, "Iteration", "Before", 1);
end

debug = struct();
debug.Options = opts;
debug.NumIterations = numel(stageRows);
debug.CancelledPCIs = cancelledPCIs;
debug.Stages = stageRows;
debug.Cancellations = cancellations;
debug.Residual = residual;

end

function rows = clusterCancellationRowsLocal(rows, minDistanceSamples, maxRows)
if isempty(rows) || height(rows) == 0
    return;
end

rows = sortrows(rows, "SSSMetric", "descend");
selected = false(height(rows), 1);
selectedStarts = zeros(0,1);

for k = 1:height(rows)
    start = rows.StartSample0(k);
    if isempty(selectedStarts) || all(abs(start - selectedStarts) >= minDistanceSamples)
        selected(k) = true;
        selectedStarts(end+1,1) = start; %#ok<AGROW>
        if nnz(selected) >= maxRows
            break;
        end
    end
end

rows = rows(selected, :);
rows = sortrows(rows, "StartSample0", "ascend");
end

function detections = addSICColumnsLocal(detections, iteration, wasCancelled)
if isempty(detections) || height(detections) == 0
    detections = emptySSSTableLikeLocal();
    return;
end

detections.SICIteration = repmat(iteration, height(detections), 1);
detections.WasCancelled = repmat(wasCancelled, height(detections), 1);
detections = movevars(detections, ["SICIteration","WasCancelled"], "After", "CandidateID");
end

function stages = stageCellLocal(pss, detections, sssDebug, iteration)
stages = {struct("Iteration", iteration, "PSSCandidates", pss, ...
    "SSSDetections", detections, "SSSDebug", sssDebug)};
end

function tbl = emptySSSTableLikeLocal()
tbl = table( ...
    zeros(0,1), zeros(0,1), false(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), false(0,1), false(0,1), ...
    'VariableNames', ["CandidateID","SICIteration","WasCancelled","NID2","NID1","PCI", ...
    "StartSample0","StartSample1","StartTimeUs","PSSPeakMetric", ...
    "PSSPeakToMedian","SegmentStartIndex","CFOHz","CFOQuality", ...
    "SSSMetric","SSSSecondMetric","SSSPeakToSecond","SSBIndex", ...
    "PBCHDMRSMetric","PBCHDMRSSecondMetric","PBCHDMRSPeakToSecond", ...
    "IsSSBIndexUsable","IsUsable"]);
end

function tbl = emptyCancellationTableForDebugLocal()
tbl = table( ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), complex(zeros(0,1)), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), false(0,1), zeros(0,1), ...
    'VariableNames', ["Iteration","PCI","StartSample0","CFOHz","ComplexGain", ...
    "GainMagnitude","SegmentPowerBefore","SegmentPowerAfter","SSSMetric", ...
    "IncludedPBCHDMRS","SSBIndex"]);
end

function [residual, cancellationTable] = cancelDetectedSSBs(iq, detections, meta, opts)
% Subtract the known SSB symbols of selected detections from the residual IQ.
% The complex gain absorbs unknown channel amplitude/phase over the SSB.
[meta, ~] = validateMetadata(meta);
base = defaultSICOptions(meta);
opts = mergeStructsLocal(base, opts);

residual = iq(:);
if isempty(detections) || height(detections) == 0
    cancellationTable = emptyCancellationTableLocal();
    return;
end

required = ["PCI","StartSample0","CFOHz","SSSMetric"];
missing = required(~ismember(required, string(detections.Properties.VariableNames)));
if ~isempty(missing)
    error("detectCellsWithSIC:MissingCancellationColumns", ...
        "detections table is missing columns: %s", strjoin(missing, ", "));
end

rows = cell(height(detections), 1);
rowCount = 0;

for k = 1:height(detections)
    det = detections(k, :);
    refOpts = mergeStructsLocal(opts.SSS, opts.Reference);
    actualSSBIndex = opts.Reference.SSBIndex;
    if ismember("SSBIndex", string(det.Properties.VariableNames)) && ...
            ismember("IsSSBIndexUsable", string(det.Properties.VariableNames)) && ...
            det.IsSSBIndexUsable
        actualSSBIndex = det.SSBIndex;
    end
    refOpts.SSBIndex = actualSSBIndex;
    ref = buildSSBSyncReference(det.PCI, meta, refOpts);
    ref = applyFractionalDelayLocal(ref, det.StartSample0 - floor(det.StartSample0));

    startIdx = floor(det.StartSample0) + 1;
    endIdx = startIdx + numel(ref) - 1;
    if endIdx < 1 || startIdx > numel(residual)
        continue;
    end

    srcStart = max(1, 2 - startIdx);
    dstStart = max(1, startIdx);
    dstEnd = min(numel(residual), endIdx);
    srcEnd = srcStart + (dstEnd - dstStart);

    refSeg = ref(srcStart:srcEnd);
    n = (0:numel(refSeg)-1).' + (srcStart-1);
    cfoHz = det.CFOHz;
    if ~isfinite(cfoHz)
        cfoHz = 0;
    end
    refSeg = refSeg .* exp(1i*2*pi*cfoHz*n/meta.SampleRateHz);
    rxSeg = residual(dstStart:dstEnd);

    alpha = (refSeg' * rxSeg) / (refSeg' * refSeg + eps);
    subtractSeg = opts.SubtractionScale * alpha * refSeg;
    beforePower = mean(abs(rxSeg).^2);
    residual(dstStart:dstEnd) = residual(dstStart:dstEnd) - subtractSeg;
    afterPower = mean(abs(residual(dstStart:dstEnd)).^2);

    rowCount = rowCount + 1;
    rows{rowCount} = table( ...
        det.PCI, det.StartSample0, cfoHz, alpha, abs(alpha), ...
        beforePower, afterPower, det.SSSMetric, ...
        opts.Reference.IncludePBCHDMRS, actualSSBIndex, ...
        'VariableNames', ["PCI","StartSample0","CFOHz","ComplexGain", ...
        "GainMagnitude","SegmentPowerBefore","SegmentPowerAfter", ...
        "SSSMetric","IncludedPBCHDMRS","SSBIndex"]);
end

if rowCount == 0
    cancellationTable = emptyCancellationTableLocal();
else
    cancellationTable = vertcat(rows{1:rowCount});
end
end

function ref = buildSSBSyncReference(pci, meta, opts)
% Build a four-symbol SS/PBCH-shaped reference from public, standards-known
% symbols. PBCH payload is not reconstructed; PBCH DM-RS is included by
% default because it is generated from PCI and SSB index.
[meta, ~] = validateMetadata(meta);
if ~isfield(opts, "SCSkHz"); opts.SCSkHz = 30; end
if ~isfield(opts, "NSizeGrid"); opts.NSizeGrid = 20; end
if ~isfield(opts, "RemoveMean"); opts.RemoveMean = true; end
if ~isfield(opts, "NormalizeReferencePower"); opts.NormalizeReferencePower = true; end
if ~isfield(opts, "IncludePSS"); opts.IncludePSS = true; end
if ~isfield(opts, "IncludeSSS"); opts.IncludeSSS = true; end
if ~isfield(opts, "IncludePBCHDMRS"); opts.IncludePBCHDMRS = true; end
if ~isfield(opts, "SSBIndex"); opts.SSBIndex = 0; end

ref = buildSSBSyncWaveform(pci, meta, opts);
end

function opts = defaultSICOptions(meta)
% Defaults for iterative near-far/overlap mitigation.
[meta, ~] = validateMetadata(meta);

opts = struct();
opts.Enable = true;
opts.MaxIterations = 6;
opts.MinNewUsableDetections = 1;
opts.MinSelectedSSSMetric = 0.08;
opts.CancellationMinDistanceSamples = round(0.25e-3 * meta.SampleRateHz);
opts.CancelAllDetectionsForSelectedPCI = true;
opts.MaxCancellationRowsPerIteration = 12;
opts.SubtractionScale = 1.0;
opts.StopWhenNoNewPCI = true;
opts.PSS = struct();
opts.SSS = struct();
opts.Reference = struct();
opts.Reference.IncludePSS = true;
opts.Reference.IncludeSSS = true;
opts.Reference.IncludePBCHDMRS = true;
opts.Reference.SSBIndex = 0;
opts.Reference.SSBStartSymbol = 2;
end

function y = applyFractionalDelayLocal(x, delaySamples)
if abs(delaySamples) < 1e-12
    y = x;
    return;
end

halfLen = 32;
n = (-halfLen:halfLen).';
h = sinc(n - delaySamples);
window = 0.5 - 0.5*cos(2*pi*(0:numel(n)-1).'/(numel(n)-1));
h = h .* window;
h = h ./ sum(h);
y = conv(x, h, "same");
end

function tbl = emptyCancellationTableLocal()
tbl = table( ...
    zeros(0,1), zeros(0,1), zeros(0,1), complex(zeros(0,1)), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), false(0,1), zeros(0,1), ...
    'VariableNames', ["PCI","StartSample0","CFOHz","ComplexGain", ...
    "GainMagnitude","SegmentPowerBefore","SegmentPowerAfter","SSSMetric", ...
    "IncludedPBCHDMRS","SSBIndex"]);
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
