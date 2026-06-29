function [cellTiming, detectionTiming] = estimateCellTiming(sssDetections, meta, opts, iq)
%estimateCellTiming Combine SSB detections into per-PCI frame timing.
%
%   [cellTiming, detectionTiming] = estimateCellTiming(sssDetections, meta)
%   optionally refines PCI-confirmed SSB timing, folds repeated detections
%   onto a 10 ms frame period, applies the configured SSB-to-frame anchor
%   offset, and combines detections by PCI.
%
%   This function estimates arrival timing at the B210. Propagation delay is
%   not corrected here.

arguments
    sssDetections table
    meta struct
    opts struct = struct()
    iq {mustBeNumeric} = []
end

[meta, ~] = validateMetadata(meta);
base = defaultTimingOptions(meta);
opts = mergeStructsLocal(base, opts);

if isempty(sssDetections) || height(sssDetections) == 0
    cellTiming = emptyCellTimingTable();
    detectionTiming = emptyDetectionTimingTable();
    return;
end

required = ["PCI","StartSample0","SSSMetric","IsUsable"];
missing = required(~ismember(required, string(sssDetections.Properties.VariableNames)));
if ~isempty(missing)
    error("estimateCellTiming:MissingColumns", ...
        "SSS detections table is missing required columns: %s", strjoin(missing, ", "));
end

det = sssDetections;
if opts.RequireUsableSSS
    det = det(det.IsUsable, :);
end
if opts.RequireUsableSSBIndex
    names = string(det.Properties.VariableNames);
    if ~ismember("IsSSBIndexUsable", names)
        cellTiming = emptyCellTimingTable();
        detectionTiming = emptyDetectionTimingTable();
        return;
    end
    det = det(det.IsSSBIndexUsable, :);
end

if isempty(det) || height(det) == 0
    cellTiming = emptyCellTimingTable();
    detectionTiming = emptyDetectionTimingTable();
    return;
end

fs = meta.SampleRateHz;
frameSamples = opts.FramePeriodMs * 1e-3 * fs;
fixedSSBOffsetSamples = opts.SSBStartOffsetUs * 1e-6 * fs;

[det, refinementSummary] = refineDetectionsLocal(det, iq, meta, opts);
[detectedSSBOffsetSamples, anchorMode] = ...
    detectedSSBOffsetSamplesLocal(det, meta, opts);

rawFrameStart = det.TimingStartSample0 - ...
    fixedSSBOffsetSamples - detectedSSBOffsetSamples;
framePhase = canonicalFramePhase(mod(rawFrameStart, frameSamples), frameSamples);
frameIndex = floor(rawFrameStart ./ frameSamples);
framePhaseNs = framePhase ./ fs * 1e9;

detectionTiming = det;
detectionTiming.FrameIndexEstimate = frameIndex;
detectionTiming.FramePhaseSamples = framePhase;
detectionTiming.FramePhaseNs = framePhaseNs;
detectionTiming.SSBIndexOffsetSamples = detectedSSBOffsetSamples;
detectionTiming.AnchorMode = repmat(anchorMode, height(det), 1);

pcis = unique(det.PCI, "stable");
rows = cell(numel(pcis), 1);
rowCount = 0;

for k = 1:numel(pcis)
    pci = pcis(k);
    idx = find(det.PCI == pci);
    if numel(idx) < opts.MinDetectionsPerCell
        continue;
    end

    if isfinite(opts.MaxDetectionsPerPCI) && numel(idx) > opts.MaxDetectionsPerPCI
        [~, order] = sort(det.SSSMetric(idx), "descend");
        idx = idx(order(1:opts.MaxDetectionsPerPCI));
    end

    phases = framePhase(idx);
    if opts.UseMetricWeights
        weights = det.SSSMetric(idx);
    else
        weights = ones(numel(idx), 1);
    end
    weights = weights(:);
    if all(weights <= 0) || any(~isfinite(weights))
        weights = ones(numel(idx), 1);
    end

    phaseEstimate = canonicalFramePhase(circularMeanSamples(phases, weights, frameSamples), frameSamples);
    residualSamples = wrapSamples(phases - phaseEstimate, frameSamples);

    if numel(residualSamples) > 1
        timingStdSamples = sqrt(sum(weights .* residualSamples.^2) / sum(weights));
    else
        timingStdSamples = 0;
    end

    timingStdNs = max(timingStdSamples / fs * 1e9, opts.TimingUncertaintyFloorNs);
    meanCFO = weightedMeanIfPresent(det, idx, "CFOHz", weights);
    meanSSSMetric = weightedMeanIfPresent(det, idx, "SSSMetric", weights);
    meanPSSMetric = weightedMeanIfPresent(det, idx, "PSSPeakMetric", weights);
    meanRefineMetric = weightedMeanIfPresent(det, idx, "TimingRefinementMetric", weights);
    meanRefineOffset = weightedMeanIfPresent(det, idx, "TimingRefinementOffsetSamples", weights);
    numRefined = nnz(det.TimingRefined(idx));

    rowCount = rowCount + 1;
    rows{rowCount} = table( ...
        pci, ...
        numel(idx), ...
        phaseEstimate, ...
        phaseEstimate / fs * 1e9, ...
        timingStdNs, ...
        meanCFO, ...
        meanPSSMetric, ...
        meanSSSMetric, ...
        min(det.StartSample0(idx)), ...
        max(det.StartSample0(idx)), ...
        min(det.TimingStartSample0(idx)), ...
        max(det.TimingStartSample0(idx)), ...
        meanRefineMetric, ...
        meanRefineOffset, ...
        numRefined, ...
        representativeStringLocal(det.TimingRefinementReference(idx)), ...
        anchorMode, ...
        'VariableNames', ["PCI","NumDetections","FramePhaseSamples", ...
        "FramePhaseNs","TimingStdNs","MeanCFOHz","MeanPSSMetric", ...
        "MeanSSSMetric","FirstDetectionStartSample0","LastDetectionStartSample0", ...
        "FirstTimingStartSample0","LastTimingStartSample0", ...
        "MeanTimingRefinementMetric","MeanTimingRefinementOffsetSamples", ...
        "NumTimingRefinedDetections","TimingRefinementReference","AnchorMode"]);
end

if rowCount == 0
    cellTiming = emptyCellTimingTable();
    return;
end

cellTiming = vertcat(rows{1:rowCount});
cellTiming = addRelativeArrivalOffsets(cellTiming, frameSamples, fs, opts);
cellTiming.TimingRefinementSummary = repmat(refinementSummary, height(cellTiming), 1);

end

function [offsetSamples, anchorMode] = detectedSSBOffsetSamplesLocal(det, meta, opts)
offsetSamples = zeros(height(det),1);
anchorMode = string(opts.AnchorMode);
if ~opts.UseDetectedSSBIndex
    return;
end

names = string(det.Properties.VariableNames);
if ~ismember("SSBIndex", names)
    return;
end
if string(opts.SSBBlockPattern) ~= "Case C"
    error("estimateCellTiming:UnsupportedSSBPattern", ...
        "Detected SSB-index timing currently supports Case C only.");
end

caseCStartSymbols = [2 8 16 22 30 36 44 50];
carrier = nrCarrierConfig;
carrier.NSizeGrid = opts.NSizeGrid;
carrier.SubcarrierSpacing = opts.SCSkHz;
ofdmInfo = nrOFDMInfo(carrier, "SampleRate", meta.SampleRateHz);
slotSymbolLengths = double(ofdmInfo.SymbolLengths(:));
symbolsPerSlot = numel(slotSymbolLengths);
slotLength = sum(slotSymbolLengths);

for k = 1:height(det)
    ssbIndex = det.SSBIndex(k);
    if ~isfinite(ssbIndex) || ssbIndex < 0 || ssbIndex > 7
        continue;
    end
    absoluteSymbol = caseCStartSymbols(ssbIndex+1);
    slotIndex = floor(absoluteSymbol/symbolsPerSlot);
    symbolInSlot = mod(absoluteSymbol,symbolsPerSlot);
    offsetSamples(k) = slotIndex*slotLength + ...
        sum(slotSymbolLengths(1:symbolInSlot));
end
anchorMode = "pbch_dmrs_ssb_index_case_c_half_frame";
end

function [det, summary] = refineDetectionsLocal(det, iq, meta, opts)
det.OriginalStartSample0 = det.StartSample0;
det.TimingStartSample0 = det.StartSample0;
det.TimingRefined = false(height(det), 1);
det.TimingRefinementOffsetSamples = zeros(height(det), 1);
det.TimingRefinementMetric = NaN(height(det), 1);
det.TimingRefinementReference = repmat("none", height(det), 1);
det.TimingRefinementReason = repmat("Timing refinement was not run.", height(det), 1);
det.TimingAnchorMode = repmat(string(opts.AnchorMode), height(det), 1);

if isempty(iq) || ~opts.EnableTimingRefinement
    if isempty(iq)
        summary = "No IQ waveform supplied; using PCI-confirmed PSS candidate start samples.";
    else
        summary = "Timing refinement disabled; using PCI-confirmed PSS candidate start samples.";
    end
    det.TimingRefinementReason(:) = summary;
    return;
end

iq = iq(:);
for r = 1:height(det)
    [refinedStart, metric, didRefine, reason, refName] = refineOneDetectionLocal( ...
        iq, det(r,:), meta, opts);
    det.TimingStartSample0(r) = refinedStart;
    det.TimingRefined(r) = didRefine;
    det.TimingRefinementOffsetSamples(r) = refinedStart - det.StartSample0(r);
    det.TimingRefinementMetric(r) = metric;
    det.TimingRefinementReference(r) = refName;
    det.TimingRefinementReason(r) = reason;
    if didRefine
        det.TimingAnchorMode(r) = "pci_confirmed_ssb_correlation";
    end
end

summary = "PCI-confirmed local SSB timing refinement used where valid.";
end

function [refinedStart0, bestMetric, didRefine, reason, refName] = refineOneDetectionLocal(iq, row, meta, opts)
refinedStart0 = row.StartSample0;
bestMetric = NaN;
didRefine = false;

[ref, refName] = buildTimingReferenceLocal(row.PCI, meta, opts);
refLen = numel(ref);
halfWindow = max(0, round(opts.TimingRefinementSearchHalfWindowUs * 1e-6 * meta.SampleRateHz));

searchStart0 = floor(row.StartSample0) - halfWindow;
searchEnd0 = ceil(row.StartSample0) + halfWindow;
if searchStart0 < 0 || searchEnd0 + refLen > numel(iq)
    reason = "Timing refinement skipped because the local search window extends outside the capture.";
    return;
end

windowStart0 = searchStart0;
windowEnd0 = searchEnd0 + refLen - 1;
segment = iq(windowStart0+1:windowEnd0+1);

cfoHz = 0;
if ismember("CFOHz", string(row.Properties.VariableNames)) && isfinite(row.CFOHz)
    cfoHz = row.CFOHz;
end
if cfoHz ~= 0
    n = (windowStart0:windowEnd0).';
    segment = segment .* exp(-1i*2*pi*cfoHz*n/meta.SampleRateHz);
end

[metric, ~] = normalizedMatchedFilterLocal(segment, ref);
if isempty(metric) || all(~isfinite(metric))
    reason = "Timing refinement skipped because the local correlation metric was invalid.";
    return;
end

[bestMetric, bestIdx] = max(metric);
frac = parabolicFractionLocal(metric, bestIdx);
starts0 = (searchStart0:searchEnd0).';
refinedStart0 = starts0(bestIdx) + frac;
didRefine = true;
reason = "Timing refined by local PCI-specific SSB correlation.";
end

function [ref, refName] = buildTimingReferenceLocal(pci, meta, opts)
refOpts = struct( ...
    "SCSkHz", opts.SCSkHz, ...
    "NSizeGrid", opts.NSizeGrid, ...
    "SSBStartSymbol", opts.SSBStartSymbol, ...
    "IncludePSS", opts.TimingRefinementIncludePSS, ...
    "IncludeSSS", opts.TimingRefinementIncludeSSS, ...
    "IncludePBCHDMRS", false, ...
    "RemoveMean", opts.TimingRefinementRemoveMean, ...
    "NormalizeReferencePower", opts.TimingRefinementNormalizeReferencePower);
ref = buildSSBSyncWaveform(pci, meta, refOpts);

parts = strings(0,1);
if opts.TimingRefinementIncludePSS; parts(end+1,1) = "PSS"; end
if opts.TimingRefinementIncludeSSS; parts(end+1,1) = "SSS"; end
refName = strjoin(parts, "+");
end

function [metric, corrValid] = normalizedMatchedFilterLocal(x, ref)
ref = ref(:);
refEnergy = sum(abs(ref).^2);
refLen = numel(ref);

corrFull = filter(conj(flipud(ref)), 1, x(:));
corrValid = corrFull(refLen:end);

windowEnergy = movsum(abs(x(:)).^2, [refLen-1 0]);
windowEnergy = windowEnergy(refLen:end);

metric = abs(corrValid).^2 ./ (refEnergy .* windowEnergy + eps);
metric(~isfinite(metric)) = 0;
end

function frac = parabolicFractionLocal(metric, idx)
if idx <= 1 || idx >= numel(metric)
    frac = 0;
    return;
end

y1 = metric(idx-1);
y2 = metric(idx);
y3 = metric(idx+1);
den = y1 - 2*y2 + y3;
if abs(den) < eps
    frac = 0;
else
    frac = 0.5 * (y1 - y3) / den;
    frac = max(min(frac, 0.5), -0.5);
end
end

function cellTiming = addRelativeArrivalOffsets(cellTiming, frameSamples, fs, opts)
phases = cellTiming.FramePhaseSamples;

switch string(opts.RelativeReference)
    case "median"
        groupCenter = circularMeanSamples(phases, ones(size(phases)), frameSamples);
        centered = wrapSamples(phases - groupCenter, frameSamples);
        referenceCentered = median(centered);
        relativeSamples = centered - referenceCentered;
    case "first"
        relativeSamples = wrapSamples(phases - phases(1), frameSamples);
    otherwise
        error("estimateCellTiming:UnknownReference", ...
            "Unknown RelativeReference value: %s", opts.RelativeReference);
end

cellTiming.RelativeArrivalOffsetSamples = relativeSamples;
cellTiming.RelativeArrivalOffsetNs = relativeSamples / fs * 1e9;
cellTiming.ReferenceMode = repmat("relative_arrival_median", height(cellTiming), 1);
end

function phase = circularMeanSamples(samples, weights, periodSamples)
samples = samples(:);
weights = weights(:);
weights = weights ./ sum(weights);
angles = 2*pi*samples/periodSamples;
z = sum(weights .* exp(1i*angles));
phase = mod(angle(z) / (2*pi) * periodSamples, periodSamples);
end

function phase = canonicalFramePhase(phase, periodSamples)
% Treat numerical wraparound at the frame boundary as zero. The 1-sample
% tolerance removes OFDM/correlation roundoff without masking real
% sub-microsecond early arrivals such as -3 us.
phase = mod(phase, periodSamples);
phase(phase >= periodSamples - 1) = 0;
end

function wrapped = wrapSamples(samples, periodSamples)
wrapped = mod(samples + periodSamples/2, periodSamples) - periodSamples/2;
end

function y = weightedMeanIfPresent(tbl, idx, varName, weights)
if ismember(varName, string(tbl.Properties.VariableNames))
    x = tbl.(varName)(idx);
    y = sum(weights(:) .* x(:)) / sum(weights);
else
    y = NaN;
end
end

function value = representativeStringLocal(values)
values = string(values);
if isempty(values)
    value = "";
    return;
end
nonempty = values(strlength(values) > 0);
if isempty(nonempty)
    value = values(1);
else
    value = nonempty(1);
end
end

function tbl = emptyCellTimingTable()
tbl = table( ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    strings(0,1), strings(0,1), zeros(0,1), zeros(0,1), strings(0,1), strings(0,1), ...
    'VariableNames', ["PCI","NumDetections","FramePhaseSamples", ...
    "FramePhaseNs","TimingStdNs","MeanCFOHz","MeanPSSMetric", ...
    "MeanSSSMetric","FirstDetectionStartSample0","LastDetectionStartSample0", ...
    "FirstTimingStartSample0","LastTimingStartSample0", ...
    "MeanTimingRefinementMetric","MeanTimingRefinementOffsetSamples", ...
    "NumTimingRefinedDetections","TimingRefinementReference","AnchorMode", ...
    "RelativeArrivalOffsetSamples","RelativeArrivalOffsetNs", ...
    "ReferenceMode","TimingRefinementSummary"]);
end

function tbl = emptyDetectionTimingTable()
tbl = table();
end

function opts = defaultTimingOptions(meta)
% Defaults for converting usable SSB detections into per-PCI frame phase.
% Propagation is deliberately not corrected here; this function estimates
% arrival timing at the receiver.
[meta, ~] = validateMetadata(meta);

opts = struct();
opts.SampleRateHz = meta.SampleRateHz;
opts.SCSkHz = 30;
opts.NSizeGrid = 20;
opts.SSBStartSymbol = 2;
opts.FramePeriodMs = 10;
opts.SSBStartOffsetUs = 0;
opts.AnchorMode = "assumed_ssb_offset";
opts.RequireUsableSSS = true;
opts.UseDetectedSSBIndex = true;
opts.RequireUsableSSBIndex = true;
opts.SSBBlockPattern = "Case C";
opts.MinDetectionsPerCell = 2;
opts.MaxDetectionsPerPCI = Inf;
opts.UseMetricWeights = true;
% Numerical reporting floor; not a calibrated B210 timing-accuracy claim.
opts.TimingUncertaintyFloorNs = 20;
opts.RelativeReference = "median";
% The interpolated PSS peak is the operational timing estimate. The
% optional PSS+SSS waveform refinement remains available for experiments,
% but real-capture validation showed that unknown PBCH content can make its
% local time-domain metric less stable than the original PSS estimate.
opts.EnableTimingRefinement = false;
opts.TimingRefinementSearchHalfWindowUs = 2;
opts.TimingRefinementIncludePSS = true;
opts.TimingRefinementIncludeSSS = true;
opts.TimingRefinementRemoveMean = true;
opts.TimingRefinementNormalizeReferencePower = true;
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
