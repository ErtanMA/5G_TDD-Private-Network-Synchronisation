function [detections, debug] = detectSSSAndPCI(iq, pssCandidates, meta, opts)
%detectSSSAndPCI Recover NID1 and PCI from PSS candidates using SSS search.
%
%   [detections, debug] = detectSSSAndPCI(iq, pssCandidates, meta, opts)
%   processes PSS candidates from estimatePSSCandidates. For each candidate,
%   it extracts the SS/PBCH-block-sized time segment, estimates and corrects
%   CFO, OFDM-demodulates the candidate, then correlates the received SSS
%   resource elements against all 336 NID1 hypotheses for the detected NID2.

arguments
    iq (:,1) {mustBeNumeric}
    pssCandidates table
    meta struct
    opts struct = struct()
end

[meta, ~] = validateMetadata(meta);
base = defaultSSSDetectionOptions(meta);
opts = mergeStructsLocal(base, opts);
opts.MaxAbsCFOHz = effectiveCFOClampHzLocal(opts);

if isempty(pssCandidates) || height(pssCandidates) == 0
    detections = emptySSSTable();
    debug = struct("Options", opts, "NumProcessed", 0);
    return;
end

refs = buildPSSReferences(meta, opts);
numProcess = min(height(pssCandidates), opts.MaxCandidatesToProcess);
rows = cell(numProcess, 1);
processed = 0;

for k = 1:numProcess
    cand = pssCandidates(k, :);
    [segment, startIdx, isValid] = extractCandidateSegment( ...
        iq, cand.StartSample0, refs.SSBSegmentLength);
    if ~isValid
        continue;
    end

    refIdx = find(refs.NID2Values == cand.NID2, 1);
    if isempty(refIdx)
        continue;
    end

    ref = refs.Waveforms{refIdx};
    if opts.EstimateCFO && string(opts.CFOEstimationMode) == "sss_grid"
        try
            [cfoHz, cfoQuality, nid1, pci, sssMetric, secondMetric, rxGrid] = ...
                recoverSSSWithCFOGridLocal(segment, ref, cand.NID2, meta, opts);
        catch
            continue;
        end
    else
        if opts.EstimateCFO
            [cfoHz, cfoQuality] = estimateCFOFromReference(segment, ref, meta.SampleRateHz, opts);
            cfoHz = max(min(cfoHz, opts.MaxAbsCFOHz), -opts.MaxAbsCFOHz);
        else
            cfoHz = 0;
            cfoQuality = NaN;
        end
        n = (0:numel(segment)-1).';
        segmentCorrected = segment .* exp(-1i*2*pi*cfoHz*n/meta.SampleRateHz);
        try
            rxGrid = demodulateSSBWaveform(segmentCorrected, meta, demodOptionsLocal(opts));
        catch
            continue;
        end
        sssRx = rxGrid(nrSSSIndices);
        [nid1, pci, sssMetric, secondMetric] = searchSSSHypotheses(sssRx, cand.NID2);
    end

    peakToSecond = sssMetric / (secondMetric + eps);
    isUsable = sssMetric >= opts.MinSSSMetric && peakToSecond >= opts.MinSSSPeakToSecond;
    [ssbIndex, dmrsMetric, dmrsSecondMetric] = ...
        searchPBCHDMRSSSBIndexLocal(rxGrid, pci);
    dmrsPeakToSecond = dmrsMetric / (dmrsSecondMetric + eps);
    isSSBIndexUsable = dmrsMetric >= opts.MinPBCHDMRSMetric && ...
        dmrsPeakToSecond >= opts.MinPBCHDMRSPeakToSecond;

    processed = processed + 1;
    rows{processed} = table( ...
        cand.CandidateID, ...
        cand.NID2, ...
        nid1, ...
        pci, ...
        cand.StartSample0, ...
        cand.StartSample1, ...
        cand.StartTimeUs, ...
        cand.PeakMetric, ...
        cand.PeakToMedian, ...
        startIdx, ...
        cfoHz, ...
        cfoQuality, ...
        sssMetric, ...
        secondMetric, ...
        peakToSecond, ...
        ssbIndex, ...
        dmrsMetric, ...
        dmrsSecondMetric, ...
        dmrsPeakToSecond, ...
        isSSBIndexUsable, ...
        isUsable, ...
        'VariableNames', ["CandidateID","NID2","NID1","PCI", ...
        "StartSample0","StartSample1","StartTimeUs","PSSPeakMetric", ...
        "PSSPeakToMedian","SegmentStartIndex","CFOHz","CFOQuality", ...
        "SSSMetric","SSSSecondMetric","SSSPeakToSecond","SSBIndex", ...
        "PBCHDMRSMetric","PBCHDMRSSecondMetric","PBCHDMRSPeakToSecond", ...
        "IsSSBIndexUsable","IsUsable"]);
end

function [cfoHz, quality, nid1, pci, bestMetric, secondMetric, rxGrid] = ...
    recoverSSSWithCFOGridLocal(segment, pssRef, nid2, meta, opts)
if ~isfinite(opts.MaxAbsCFOHz)
    error("detectSSSAndPCI:FiniteCFOGridRequired", ...
        "SSS-grid CFO estimation requires a finite MaxAbsCFOHz.");
end

cfoGrid = (-opts.MaxAbsCFOHz:opts.CFOGridStepHz:opts.MaxAbsCFOHz).';
if isempty(cfoGrid) || ~any(abs(cfoGrid) < eps)
    cfoGrid(end+1,1) = 0;
end
cfoGrid = unique(cfoGrid);

[gridMetrics, ~, ~, ~] = ...
    evaluateCFOGridLocal(segment, nid2, meta, opts, cfoGrid);
[~, coarseIdx] = max(gridMetrics);
coarseBestHz = cfoGrid(coarseIdx);

fineMinHz = max(-opts.MaxAbsCFOHz, coarseBestHz - opts.CFOGridStepHz);
fineMaxHz = min(opts.MaxAbsCFOHz, coarseBestHz + opts.CFOGridStepHz);
fineGrid = (fineMinHz:opts.CFOFineStepHz:fineMaxHz).';
if isempty(fineGrid) || ~any(abs(fineGrid-coarseBestHz) < eps)
    fineGrid(end+1,1) = coarseBestHz;
end
fineGrid = unique(fineGrid);
[fineMetrics, fineNID1, finePCI, fineSecond] = ...
    evaluateCFOGridLocal(segment, nid2, meta, opts, fineGrid);

[bestMetric, idx] = max(fineMetrics);
cfoHz = fineGrid(idx);
nid1 = fineNID1(idx);
pci = finePCI(idx);
secondMetric = fineSecond(idx);
sortedGridMetrics = sort(fineMetrics, "descend");
if numel(sortedGridMetrics) > 1
    quality = bestMetric / (sortedGridMetrics(2) + eps);
else
    quality = NaN;
end

[phaseCFOHz, phaseQuality] = estimateCFOFromReference( ...
    segment, pssRef, meta.SampleRateHz, opts);
phaseCFOHz = max(min(phaseCFOHz, opts.MaxAbsCFOHz), -opts.MaxAbsCFOHz);
if phaseQuality >= opts.MinCFOPhaseQuality && ...
        abs(phaseCFOHz - cfoHz) <= opts.CFOPhaseGridAgreementHz
    [phaseMetric, phaseNID1, phasePCI, phaseSecond] = ...
        evaluateCFOGridLocal(segment, nid2, meta, opts, phaseCFOHz);
    cfoHz = phaseCFOHz;
    nid1 = phaseNID1;
    pci = phasePCI;
    bestMetric = phaseMetric;
    secondMetric = phaseSecond;
    quality = phaseQuality;
end

n = (0:numel(segment)-1).';
corrected = segment .* exp(-1i*2*pi*cfoHz*n/meta.SampleRateHz);
rxGrid = demodulateSSBWaveform(corrected, meta, demodOptionsLocal(opts));
end

function [metrics, nid1Values, pciValues, secondValues] = ...
    evaluateCFOGridLocal(segment, nid2, meta, opts, cfoGrid)
metrics = -Inf(size(cfoGrid));
nid1Values = zeros(size(cfoGrid));
pciValues = zeros(size(cfoGrid));
secondValues = zeros(size(cfoGrid));
n = (0:numel(segment)-1).';
for g = 1:numel(cfoGrid)
    corrected = segment .* exp(-1i*2*pi*cfoGrid(g)*n/meta.SampleRateHz);
    rxGrid = demodulateSSBWaveform(corrected, meta, demodOptionsLocal(opts));
    sssRx = rxGrid(nrSSSIndices);
    [nid1Values(g), pciValues(g), metrics(g), secondValues(g)] = ...
        searchSSSHypotheses(sssRx, nid2);
end
end

function demodOpts = demodOptionsLocal(opts)
demodOpts = struct( ...
    "SCSkHz", opts.SCSkHz, ...
    "NSizeGrid", opts.NSizeGrid, ...
    "SSBStartSymbol", opts.SSBStartSymbol, ...
    "NCellID", 0);
end

if processed == 0
    detections = emptySSSTable();
else
    detections = vertcat(rows{1:processed});
    switch string(opts.SortBy)
        case "SSSMetric"
            detections = sortrows(detections, "SSSMetric", "descend");
        case "StartSample0"
            detections = sortrows(detections, "StartSample0", "ascend");
        case "PSSPeakMetric"
            detections = sortrows(detections, "PSSPeakMetric", "descend");
        otherwise
            error("detectSSSAndPCI:UnknownSort", ...
                "Unknown SortBy value: %s", opts.SortBy);
    end
end

debug = struct();
debug.Options = opts;
debug.NumProcessed = processed;
debug.PSSReferenceLength = refs.PSSReferenceLength;
debug.SSBSegmentLength = refs.SSBSegmentLength;

end

function [segment, startIdx, isValid] = extractCandidateSegment(iq, startSample0, segmentLength)
startIdx = round(startSample0) + 1;
endIdx = startIdx + segmentLength - 1;
isValid = startIdx >= 1 && endIdx <= numel(iq);
if isValid
    segment = iq(startIdx:endIdx);
else
    segment = complex(zeros(segmentLength, 1));
end
end

function [cfoHz, quality] = estimateCFOFromReference(segment, ref, fs, opts)
segment = segment(:);
ref = ref(:);
len = min(numel(segment), numel(ref));
segment = segment(1:len);
ref = ref(1:len);

weights = abs(ref).^2;
mask = weights >= opts.CFOReferencePowerThreshold * max(weights);
if nnz(mask) < 8
    cfoHz = 0;
    quality = 0;
    return;
end

n = (0:len-1).';
z = segment(mask) .* conj(ref(mask));
t = n(mask) / fs;
w = weights(mask);
w = w ./ sum(w);

phase = unwrap(angle(z));
t0 = sum(w .* t);
p0 = sum(w .* phase);
den = sum(w .* (t - t0).^2);
if den <= eps
    cfoHz = 0;
else
    slope = sum(w .* (t - t0) .* (phase - p0)) / den;
    cfoHz = slope / (2*pi);
end

quality = abs(sum(z)) / sqrt(sum(abs(segment(mask)).^2) * sum(abs(ref(mask)).^2) + eps);
end

function [bestNID1, bestPCI, bestMetric, secondMetric] = searchSSSHypotheses(sssRx, nid2)
sssRx = sssRx(:);
rxEnergy = sum(abs(sssRx).^2);
metrics = zeros(336, 1);

for nid1 = 0:335
    pci = 3*nid1 + nid2;
    sssRef = nrSSS(pci);
    sssRef = sssRef(:);
    metrics(nid1 + 1) = abs(sssRef' * sssRx)^2 / ...
        (sum(abs(sssRef).^2) * rxEnergy + eps);
end

[sortedMetrics, order] = sort(metrics, "descend");
bestNID1 = order(1) - 1;
bestPCI = 3*bestNID1 + nid2;
bestMetric = sortedMetrics(1);
if numel(sortedMetrics) >= 2
    secondMetric = sortedMetrics(2);
else
    secondMetric = 0;
end
end

function [bestIndex, bestMetric, secondMetric] = ...
    searchPBCHDMRSSSBIndexLocal(rxGrid, pci)
dmrsIndices = nrPBCHDMRSIndices(pci);
[~, dmrsSymbols] = ind2sub(size(rxGrid), dmrsIndices);
dmrsRx = rxGrid(dmrsIndices);
activeSymbols = unique(dmrsSymbols, "stable");
metrics = zeros(8,1);
for ssbIndex = 0:7
    ref = nrPBCHDMRS(pci, ssbIndex);
    symbolMetrics = zeros(numel(activeSymbols),1);
    symbolWeights = zeros(numel(activeSymbols),1);
    for s = 1:numel(activeSymbols)
        mask = dmrsSymbols == activeSymbols(s);
        rxPart = dmrsRx(mask);
        refPart = ref(mask);
        symbolMetrics(s) = abs(refPart' * rxPart)^2 / ...
            (sum(abs(refPart).^2)*sum(abs(rxPart).^2) + eps);
        symbolWeights(s) = nnz(mask);
    end
    % Residual CFO and common phase error can rotate separate OFDM symbols
    % differently. Combining normalized correlation powers preserves the
    % SSB-index evidence without requiring those symbol phases to agree.
    metrics(ssbIndex+1) = sum(symbolWeights.*symbolMetrics) / ...
        sum(symbolWeights);
end
[sortedMetrics, order] = sort(metrics, "descend");
bestIndex = order(1)-1;
bestMetric = sortedMetrics(1);
secondMetric = sortedMetrics(2);
end

function tbl = emptySSSTable()
tbl = table( ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), false(0,1), false(0,1), ...
    'VariableNames', ["CandidateID","NID2","NID1","PCI", ...
    "StartSample0","StartSample1","StartTimeUs","PSSPeakMetric", ...
        "PSSPeakToMedian","SegmentStartIndex","CFOHz","CFOQuality", ...
        "SSSMetric","SSSSecondMetric","SSSPeakToSecond","SSBIndex", ...
        "PBCHDMRSMetric","PBCHDMRSSecondMetric","PBCHDMRSPeakToSecond", ...
        "IsSSBIndexUsable","IsUsable"]);
end

function refs = buildPSSReferences(meta, opts)
% Build the same PSS references used by the coarse search so CFO estimation
% and candidate extraction stay aligned with the detector.
[meta, ~] = validateMetadata(meta);
base = defaultSSSDetectionOptions(meta);
opts = mergeStructsLocal(base, opts);

nid2Values = 0:2;
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
        "RemoveMean", true, ...
        "NormalizeReferencePower", true);
    waveform = buildSSBSyncWaveform(nid2, meta, refOpts);

    waveforms{k} = waveform;
    energies(k) = sum(abs(waveform).^2);
end

refs = struct();
refs.NID2Values = nid2Values;
refs.Waveforms = waveforms;
refs.Energies = energies;
refs.PSSReferenceLength = numel(waveforms{1});
segmentOpts = struct( ...
    "SCSkHz", opts.SCSkHz, ...
    "NSizeGrid", opts.NSizeGrid, ...
    "SSBStartSymbol", opts.SSBStartSymbol, ...
    "NumOutputSymbols", 4, ...
    "IncludePSS", true, ...
    "IncludeSSS", false, ...
    "IncludePBCHDMRS", false);
segmentWaveform = buildSSBSyncWaveform(0, meta, segmentOpts);
refs.SSBSegmentLength = numel(segmentWaveform);
refs.SampleRateHz = meta.SampleRateHz;
refs.SCSkHz = opts.SCSkHz;
refs.NSizeGrid = opts.NSizeGrid;
end

function opts = defaultSSSDetectionOptions(meta)
% Defaults for candidate CFO correction and SSS/PCI recovery.
[meta, ~] = validateMetadata(meta);

opts = struct();
opts.SampleRateHz = meta.SampleRateHz;
opts.SCSkHz = 30;
opts.NSizeGrid = 20;
opts.SSBStartSymbol = 2;
opts.MaxCandidatesToProcess = 24;
opts.EstimateCFO = true;
opts.MaxAbsCFOHz = NaN;
opts.CFOEstimationMode = "sss_grid";
% Search resolution and agreement gates are implementation settings.
opts.CFOGridStepHz = 1500;
opts.CFOFineStepHz = 50;
opts.CFOPhaseGridAgreementHz = 750;
opts.MinCFOPhaseQuality = 0.1;
opts.CFOReferencePowerThreshold = 0.05;
% Empirical normalized-correlation acceptance gates, not 3GPP limits.
opts.MinSSSMetric = 0.06;
opts.MinSSSPeakToSecond = 1.15;
opts.MinPBCHDMRSMetric = 0.06;
opts.MinPBCHDMRSPeakToSecond = 1.15;
opts.SortBy = "SSSMetric";
end

function limitHz = effectiveCFOClampHzLocal(opts)
% NaN means derive the clamp from the SSB numerology. A half-subcarrier
% default is less arbitrary than a fixed Hz value and scales with SCS.
limitHz = opts.MaxAbsCFOHz;
if isempty(limitHz) || isnan(limitHz)
    limitHz = 0.5 * opts.SCSkHz * 1e3;
end
if limitHz < 0
    error("detectSSSAndPCI:InvalidCFOClamp", ...
        "MaxAbsCFOHz must be non-negative, Inf, or NaN for SCS-derived default.");
end
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
