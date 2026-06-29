function [candidates, debug] = estimatePSSCandidates(iq, meta, opts)
%estimatePSSCandidates Handmade matched-filter NR PSS search.
%
%   [candidates, debug] = estimatePSSCandidates(iq, meta, opts) searches a
%   baseband IQ capture for the three legal NR PSS hypotheses. It returns a
%   table of candidate SSB starts. PCI is not recovered here; Phase 5 uses
%   SSS to recover NID1 and PCI from these PSS candidates.

arguments
    iq (:,1) {mustBeNumeric}
    meta struct
    opts struct = struct()
end

[meta, ~] = validateMetadata(meta);
base = defaultPSSSearchOptions(meta);
opts = mergeStructsLocal(base, opts);

x = iq(:);
if opts.RemoveMean
    x = x - mean(x);
end

refs = buildPSSReferences(meta, opts);

allRows = cell(numel(refs.NID2Values), 1);
debug = struct();
debug.Options = opts;
debug.ReferenceInfo = rmfield(refs, "Waveforms");
debug.MetricMedian = zeros(numel(refs.NID2Values), 1);
debug.MetricMax = zeros(numel(refs.NID2Values), 1);

for k = 1:numel(refs.NID2Values)
    nid2 = refs.NID2Values(k);
    ref = refs.Waveforms{k};
    [metric, corrValid] = normalizedMatchedFilter(x, ref);

    metricMedian = median(metric);
    adaptiveThreshold = max(opts.MinPeakMetric, ...
        metricMedian * opts.MinPeakToMedianRatio);

    [peakIdx, peakMetric] = pickClusteredPeaks(metric, ...
        adaptiveThreshold, opts.MinPeakDistanceSamples, opts.MaxCandidatesPerNID2);

    numPeaks = numel(peakIdx);
    rows = table();
    if numPeaks > 0
        frac = zeros(numPeaks, 1);
        corrMag = zeros(numPeaks, 1);
        for p = 1:numPeaks
            frac(p) = parabolicFraction(metric, peakIdx(p));
            corrMag(p) = abs(corrValid(peakIdx(p)));
        end

        startSample0 = double(peakIdx(:) - 1) + frac;
        rows = table( ...
            repmat(nid2, numPeaks, 1), ...
            startSample0, ...
            startSample0 + 1, ...
            startSample0 ./ meta.SampleRateHz * 1e6, ...
            peakMetric(:), ...
            peakMetric(:) ./ (metricMedian + eps), ...
            corrMag(:), ...
            repmat(refs.ReferenceLength, numPeaks, 1), ...
            frac(:), ...
            'VariableNames', ["NID2","StartSample0","StartSample1", ...
            "StartTimeUs","PeakMetric","PeakToMedian","CorrelationMagnitude", ...
            "ReferenceLength","FractionalPeakSamples"]);
    end

    allRows{k} = rows;
    debug.MetricMedian(k) = metricMedian;
    debug.MetricMax(k) = max(metric);
end

candidates = vertcat(allRows{:});
if isempty(candidates)
    candidates = emptyCandidateTable();
    return;
end

switch string(opts.SortBy)
    case "PeakMetric"
        candidates = sortrows(candidates, "PeakMetric", "descend");
    case "StartSample0"
        candidates = sortrows(candidates, "StartSample0", "ascend");
    otherwise
        error("estimatePSSCandidates:UnknownSort", ...
            "Unknown SortBy value: %s", opts.SortBy);
end

if height(candidates) > opts.MaxTotalCandidates
    candidates = candidates(1:opts.MaxTotalCandidates, :);
end

candidates.CandidateID = (1:height(candidates)).';
candidates = movevars(candidates, "CandidateID", "Before", 1);

end

function [metric, corrValid] = normalizedMatchedFilter(x, ref)
ref = ref(:);
refEnergy = sum(abs(ref).^2);
refLen = numel(ref);

b = conj(flipud(ref));
if exist("fftfilt", "file") == 2 && numel(x) > 10*numel(b)
    corrFull = fftfilt(double(b), double(x));
else
    corrFull = filter(b, 1, x);
end
corrValid = corrFull(refLen:end);

windowEnergy = movsum(abs(x).^2, [refLen-1 0]);
windowEnergy = windowEnergy(refLen:end);

metric = abs(corrValid).^2 ./ (refEnergy .* windowEnergy + eps);
metric(~isfinite(metric)) = 0;
end

function [selectedIdx, selectedMetric] = pickClusteredPeaks(metric, threshold, minDistance, maxPeaks)
metric = metric(:);
if numel(metric) < 3
    selectedIdx = [];
    selectedMetric = [];
    return;
end

isPeak = false(size(metric));
isPeak(1) = metric(1) >= metric(2) && metric(1) >= threshold;
isPeak(2:end-1) = metric(2:end-1) > metric(1:end-2) & ...
    metric(2:end-1) >= metric(3:end) & ...
    metric(2:end-1) >= threshold;
isPeak(end) = metric(end) > metric(end-1) && metric(end) >= threshold;

candidateIdx = find(isPeak);
if isempty(candidateIdx)
    selectedIdx = [];
    selectedMetric = [];
    return;
end

[~, order] = sort(metric(candidateIdx), "descend");
candidateIdx = candidateIdx(order);

selectedIdx = zeros(0,1);
for k = 1:numel(candidateIdx)
    idx = candidateIdx(k);
    if isempty(selectedIdx) || all(abs(idx - selectedIdx) >= minDistance)
        selectedIdx(end+1,1) = idx; 
        if numel(selectedIdx) >= maxPeaks
            break;
        end
    end
end

selectedMetric = metric(selectedIdx);
[selectedIdx, order] = sort(selectedIdx, "ascend");
selectedMetric = selectedMetric(order);
end

function frac = parabolicFraction(metric, idx)
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

function tbl = emptyCandidateTable()
tbl = table( ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', ["CandidateID","NID2","StartSample0","StartSample1", ...
    "StartTimeUs","PeakMetric","PeakToMedian","CorrelationMagnitude", ...
    "ReferenceLength","FractionalPeakSamples"]);
end

function refs = buildPSSReferences(meta, opts)
% Build one time-domain SSB-shaped PSS reference per NID2 hypothesis.
% The waveform starts at the SSB boundary, so matched-filter peak locations
% are directly interpretable as candidate SSB start samples.
[meta, ~] = validateMetadata(meta);
base = defaultPSSSearchOptions(meta);
opts = mergeStructsLocal(base, opts);

nid2Values = opts.NID2Values(:).';
waveforms = cell(numel(nid2Values), 1);
energies = zeros(numel(nid2Values), 1);

for k = 1:numel(nid2Values)
    nid2 = nid2Values(k);
    refOpts = struct( ...
        "SCSkHz", opts.SCSkHz, ...
        "NSizeGrid", opts.NSizeGrid, ...
        "SSBStartSymbol", opts.SSBStartSymbol, ...
        "NumOutputSymbols", 1, ...
        "IncludePSS", true, ...
        "IncludeSSS", false, ...
        "IncludePBCHDMRS", false, ...
        "RemoveMean", opts.RemoveMean, ...
        "NormalizeReferencePower", opts.NormalizeReferencePower);
    waveform = buildSSBSyncWaveform(nid2, meta, refOpts);

    waveforms{k} = waveform;
    energies(k) = sum(abs(waveform).^2);
end

refs = struct();
refs.NID2Values = nid2Values;
refs.Waveforms = waveforms;
refs.Energies = energies;
refs.ReferenceLength = numel(waveforms{1});
refs.SampleRateHz = meta.SampleRateHz;
refs.SCSkHz = opts.SCSkHz;
refs.NSizeGrid = opts.NSizeGrid;
end

function opts = defaultPSSSearchOptions(meta)
% Defaults for blind NR PSS candidate search.
[meta, ~] = validateMetadata(meta);

opts = struct();
opts.SampleRateHz = meta.SampleRateHz;
opts.SCSkHz = 30;
opts.NSizeGrid = 20;
opts.SSBStartSymbol = 2;
opts.NID2Values = 0:2;
opts.MaxCandidatesPerNID2 = 24;
opts.MaxTotalCandidates = 72;
% Empirical acceptance settings. These are calibrated implementation
% parameters, not values prescribed by 3GPP.
opts.MinPeakMetric = 0.005;
opts.MinPeakToMedianRatio = 12;
opts.MinPeakDistanceSamples = round(0.25e-3 * meta.SampleRateHz);
opts.RemoveMean = true;
opts.NormalizeReferencePower = true;
opts.SortBy = "PeakMetric";
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
