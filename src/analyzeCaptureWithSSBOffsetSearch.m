function result = analyzeCaptureWithSSBOffsetSearch(iq, meta, opts)
%analyzeCaptureWithSSBOffsetSearch Analyze IQ with SSB offset hypotheses.
%
%   Public-network captures may contain the SSB away from the tuned center
%   frequency. This wrapper digitally shifts the saved IQ by several
%   hypothesized SSB offsets, runs the normal receiver chain for each
%   hypothesis, and returns the strongest usable result. It does not change
%   timing sample indices because a frequency shift does not move samples in
%   time.

arguments
    iq (:,1) {mustBeNumeric}
    meta struct
    opts struct = struct()
end

[meta, ~] = validateMetadata(meta);
searchOpts = defaultOffsetSearchOptionsLocal(meta);
if isfield(opts, "SSBOffsetSearch")
    searchOpts = mergeStructsLocal(searchOpts, opts.SSBOffsetSearch);
end

if ~searchOpts.Enable
    result = analyzeCapture(iq, meta, removeOffsetSearchOptionsLocal(opts));
    result.SSBOffsetSearch = disabledSearchSummaryLocal(meta);
    return;
end

coarseOffsetsHz = getCoarseOffsetHypothesesLocal(searchOpts);
if isempty(coarseOffsetsHz)
    error("analyzeCaptureWithSSBOffsetSearch:NoOffsets", ...
        "No SSB frequency-offset hypotheses were generated.");
end

analysisOpts = removeOffsetSearchOptionsLocal(opts);
results = {};
rows = {};
testedOffsetsHz = zeros(0,1);

[results, rows, testedOffsetsHz] = analyzeOffsetsLocal( ...
    iq, meta, analysisOpts, searchOpts, coarseOffsetsHz, ...
    results, rows, testedOffsetsHz, "coarse");

if searchOpts.EnableFineSearch && ~isempty(rows)
    coarseSummary = vertcat(rows{:});
    coarseSummary = sortrows(coarseSummary, ...
        ["CoarseScore","MaxPSSMetric","BestAnySSSMetric"], ...
        ["descend","descend","descend"]);
    numSeeds = min(searchOpts.NumFineSeeds, height(coarseSummary));
    seedOffsets = coarseSummary.DigitalShiftHz(1:numSeeds);
    fineOffsetsHz = getFineOffsetHypothesesLocal(seedOffsets, searchOpts);
    fineOffsetsHz = setdiff(fineOffsetsHz, testedOffsetsHz, "stable");
    [results, rows, testedOffsetsHz] = analyzeOffsetsLocal( ...
        iq, meta, analysisOpts, searchOpts, fineOffsetsHz, ...
        results, rows, testedOffsetsHz, "fine");
end

summary = vertcat(rows{:});
summary = sortrows(summary, ["IsConfirmed","Score","NumCellTimingRows", ...
    "MaxUsableDetectionsPerPCI","NumUsableDetections","BestUsableSSSMetric"], ...
    ["descend","descend","descend","descend","descend","descend"]);
bestIndex = summary.ResultIndex(1);
result = results{bestIndex};
result.SSBOffsetSearch.Enable = true;
result.SSBOffsetSearch.SearchTable = summary;
result.SSBOffsetSearch.SelectedOffsetHz = summary.DigitalShiftHz(1);
result.SSBOffsetSearch.SelectedSSBCenterHz = summary.HypothesizedSSBCenterHz(1);
result.SSBOffsetSearch.SelectedCaptureCenterHz = meta.CenterFrequencyHz;
result.SSBOffsetSearch.NumOffsetsTested = numel(testedOffsetsHz);
result.SSBOffsetSearch.IsConfirmed = summary.IsConfirmed(1);
result.SSBOffsetSearch.ConfirmationReason = summary.ConfirmationReason(1);
result.SSBOffsetSearch.Statement = "Selected the SSB frequency-offset hypothesis with the strongest usable PCI/timing result.";
if ~result.SSBOffsetSearch.IsConfirmed
    result.SSBOffsetSearch.Statement = "Offset search found only an unconfirmed PCI hypothesis. Treat this as diagnostics, not a validated cell detection.";
end

if searchOpts.RequireUsable && summary.NumUsableDetections(1) == 0
    result.SSBOffsetSearch.Statement = "No offset produced a usable PCI detection; returned the strongest non-usable hypothesis for diagnostics.";
end
end

function [results, rows, testedOffsetsHz] = analyzeOffsetsLocal( ...
    iq, meta, analysisOpts, searchOpts, offsetsHz, results, rows, testedOffsetsHz, stageName)
for k = 1:numel(offsetsHz)
    offsetHz = offsetsHz(k);
    if searchOpts.Verbose
        fprintf("SSB offset search (%s): %d/%d, shift %.0f Hz\n", ...
            stageName, k, numel(offsetsHz), offsetHz);
    end
    shifted = applyDigitalShiftLocal(iq, meta.SampleRateHz, offsetHz);
    candidate = analyzeCapture(shifted, meta, analysisOpts);
    candidate.SSBOffsetSearch = struct( ...
        "Enable", true, ...
        "SelectedOffsetHz", offsetHz, ...
        "SelectedSSBCenterHz", meta.CenterFrequencyHz + offsetHz, ...
        "Statement", "IQ was digitally shifted before PSS/SSS search; reported timing sample indices remain in the original capture time base.");

    resultIndex = numel(results) + 1;
    results{resultIndex,1} = candidate;
    rows{resultIndex,1} = summarizeCandidateLocal( ...
        candidate, meta, offsetHz, resultIndex, searchOpts);
    testedOffsetsHz(end+1,1) = offsetHz; %#ok<AGROW>
end
end

function shifted = applyDigitalShiftLocal(iq, fs, offsetHz)
% Positive offsetHz means the SSB is assumed above the tuned center; multiply
% by a negative complex exponential to move that SSB toward baseband.
n = (0:numel(iq)-1).';
rot = exp(-1i * 2*pi * offsetHz/fs .* n);
shifted = iq(:) .* rot;
end

function row = summarizeCandidateLocal(result, meta, offsetHz, resultIndex, searchOpts)
usable = table();
if isfield(result, "SSSDetections") && ~isempty(result.SSSDetections)
    usable = result.SSSDetections(result.SSSDetections.IsUsable, :);
end

allSSS = table();
if isfield(result, "SSSDetections") && ~isempty(result.SSSDetections)
    allSSS = result.SSSDetections;
end

numUsable = height(usable);
numPSSCandidates = 0;
maxPSSMetric = NaN;
if isfield(result, "PSSCandidates") && ~isempty(result.PSSCandidates)
    numPSSCandidates = height(result.PSSCandidates);
    if ismember("PeakMetric", string(result.PSSCandidates.Properties.VariableNames))
        maxPSSMetric = max(result.PSSCandidates.PeakMetric, [], "omitnan");
    end
end

numTiming = 0;
if isfield(result, "CellTiming") && ~isempty(result.CellTiming)
    numTiming = height(result.CellTiming);
end

bestPCI = NaN;
bestUsableSSSMetric = NaN;
bestUsablePeakToSecond = NaN;
medianAbsCFOHz = NaN;
usablePCIs = "";
dominantUsablePCI = NaN;
maxUsableDetectionsPerPCI = 0;
if numUsable > 0
    [~, idx] = max(usable.SSSMetric);
    bestPCI = usable.PCI(idx);
    bestUsableSSSMetric = usable.SSSMetric(idx);
    bestUsablePeakToSecond = usable.SSSPeakToSecond(idx);
    medianAbsCFOHz = median(abs(usable.CFOHz), "omitnan");
    uniqueUsablePCIs = unique(usable.PCI, "stable");
    countsPerPCI = arrayfun(@(pci) nnz(usable.PCI == pci), uniqueUsablePCIs);
    [maxUsableDetectionsPerPCI, dominantIdx] = max(countsPerPCI);
    dominantUsablePCI = uniqueUsablePCIs(dominantIdx);
    usablePCIs = strjoin(string(uniqueUsablePCIs), "|");
end

bestAnyPCI = NaN;
bestAnySSSMetric = NaN;
bestAnyPeakToSecond = NaN;
if height(allSSS) > 0
    [~, idx] = max(allSSS.SSSMetric);
    bestAnyPCI = allSSS.PCI(idx);
    bestAnySSSMetric = allSSS.SSSMetric(idx);
    bestAnyPeakToSecond = allSSS.SSSPeakToSecond(idx);
end

timingSpreadNs = NaN;
if numTiming > 0
    timingSpreadNs = median(result.CellTiming.TimingStdNs, "omitnan");
end

distanceToSearchEdgeHz = min(abs(offsetHz - searchOpts.MinOffsetHz), ...
    abs(offsetHz - searchOpts.MaxOffsetHz));
isBoundaryOffset = distanceToSearchEdgeHz <= searchOpts.BoundaryToleranceHz;
maxAbsCFOHz = getMaxAbsCFOLocal(result);
isCFOAtClamp = isfinite(maxAbsCFOHz) && isfinite(medianAbsCFOHz) && ...
    medianAbsCFOHz >= searchOpts.CFOClampWarningFraction * maxAbsCFOHz;
isConfirmed = maxUsableDetectionsPerPCI >= searchOpts.MinConfirmedUsableDetections && ...
    numTiming > 0 && (~searchOpts.RejectBoundaryForConfirmation || ~isBoundaryOffset) && ...
    (~searchOpts.RejectCFOClampForConfirmation || ~isCFOAtClamp);
confirmationReason = getConfirmationReasonLocal( ...
    isConfirmed, maxUsableDetectionsPerPCI, searchOpts, ...
    isBoundaryOffset, isCFOAtClamp, numTiming);

% Usable PCI/timing rows dominate the score. SSS metrics are only tie-breaks.
score = 1e9*double(isConfirmed) + 1e6*numTiming + 1e4*numUsable + ...
    1e5*maxUsableDetectionsPerPCI + ...
    100*nanToZeroLocal(bestUsableSSSMetric) + ...
    nanToZeroLocal(bestUsablePeakToSecond) + ...
    0.01*nanToZeroLocal(bestAnySSSMetric);
coarseScore = 1000*nanToZeroLocal(maxPSSMetric) + ...
    10*nanToZeroLocal(bestAnySSSMetric) + ...
    nanToZeroLocal(bestAnyPeakToSecond);

row = table( ...
    resultIndex, offsetHz, meta.CenterFrequencyHz + offsetHz, ...
    numPSSCandidates, maxPSSMetric, numUsable, usablePCIs, numTiming, bestPCI, ...
    dominantUsablePCI, maxUsableDetectionsPerPCI, ...
    bestUsableSSSMetric, bestUsablePeakToSecond, ...
    medianAbsCFOHz, timingSpreadNs, bestAnyPCI, ...
    bestAnySSSMetric, bestAnyPeakToSecond, isBoundaryOffset, ...
    isCFOAtClamp, isConfirmed, confirmationReason, coarseScore, score, ...
    'VariableNames', ["ResultIndex","DigitalShiftHz", ...
    "HypothesizedSSBCenterHz","NumPSSCandidates","MaxPSSMetric", ...
    "NumUsableDetections","UsablePCIs", ...
    "NumCellTimingRows","BestUsablePCI","DominantUsablePCI", ...
    "MaxUsableDetectionsPerPCI","BestUsableSSSMetric", ...
    "BestUsableSSSPeakToSecond","MedianAbsCFOHz", ...
    "MedianTimingSpreadNs","BestAnyPCI","BestAnySSSMetric", ...
    "BestAnySSSPeakToSecond","IsBoundaryOffset","IsCFOAtClamp", ...
    "IsConfirmed","ConfirmationReason","CoarseScore","Score"]);
end

function maxAbsCFOHz = getMaxAbsCFOLocal(result)
maxAbsCFOHz = NaN;
if isfield(result, "SSSDebug") && isfield(result.SSSDebug, "Options") && ...
        isfield(result.SSSDebug.Options, "MaxAbsCFOHz")
    maxAbsCFOHz = result.SSSDebug.Options.MaxAbsCFOHz;
end
end

function reason = getConfirmationReasonLocal( ...
    isConfirmed, maxUsableDetectionsPerPCI, searchOpts, ...
    isBoundaryOffset, isCFOAtClamp, numTiming)
if isConfirmed
    reason = "CONFIRMED";
elseif maxUsableDetectionsPerPCI < searchOpts.MinConfirmedUsableDetections
    reason = "TOO_FEW_REPEATED_DETECTIONS_FOR_SAME_PCI";
elseif numTiming == 0
    reason = "NO_CELL_TIMING_ROW";
elseif searchOpts.RejectBoundaryForConfirmation && isBoundaryOffset
    reason = "SELECTED_OFFSET_AT_SEARCH_BOUNDARY";
elseif searchOpts.RejectCFOClampForConfirmation && isCFOAtClamp
    reason = "CFO_ESTIMATE_AT_CLAMP";
else
    reason = "UNCONFIRMED";
end
end

function y = nanToZeroLocal(x)
if isempty(x) || ~isfinite(x)
    y = 0;
else
    y = x;
end
end

function offsetsHz = getCoarseOffsetHypothesesLocal(opts)
if isfield(opts, "OffsetsHz") && ~isempty(opts.OffsetsHz)
    offsetsHz = double(opts.OffsetsHz(:));
else
    offsetsHz = (opts.MinOffsetHz:opts.StepHz:opts.MaxOffsetHz).';
end

if opts.IncludeZero && ~any(abs(offsetsHz) < eps)
    offsetsHz(end+1,1) = 0;
end
offsetsHz = unique(round(offsetsHz(:)), "stable");
if string(opts.SearchOrder) == "center_first"
    [~, order] = sort(abs(offsetsHz), "ascend");
    offsetsHz = offsetsHz(order);
end
end

function offsetsHz = getFineOffsetHypothesesLocal(seedOffsets, opts)
offsetsHz = zeros(0,1);
for k = 1:numel(seedOffsets)
    seed = seedOffsets(k);
    local = (seed-opts.FineSpanHz:opts.FineStepHz:seed+opts.FineSpanHz).';
    local = local(local >= opts.MinOffsetHz & local <= opts.MaxOffsetHz);
    offsetsHz = [offsetsHz; local]; %#ok<AGROW>
end
if opts.IncludeZero && ~any(abs(offsetsHz) < eps)
    offsetsHz(end+1,1) = 0;
end
offsetsHz = unique(round(offsetsHz(:)), "stable");
end

function opts = defaultOffsetSearchOptionsLocal(meta)
fs = meta.SampleRateHz;
spanHz = min(3e6, 0.40*fs);
opts = struct();
opts.Enable = true;
opts.MinOffsetHz = -spanHz;
opts.MaxOffsetHz = spanHz;
opts.StepHz = 250e3;
opts.OffsetsHz = [];
opts.IncludeZero = true;
opts.SearchOrder = "ascending";
opts.RequireUsable = false;
opts.Verbose = false;
opts.EnableFineSearch = true;
opts.NumFineSeeds = 4;
opts.FineSpanHz = 150e3;
opts.FineStepHz = 15e3;
opts.MinConfirmedUsableDetections = 2;
opts.RejectBoundaryForConfirmation = true;
opts.BoundaryToleranceHz = 1;
opts.RejectCFOClampForConfirmation = true;
opts.CFOClampWarningFraction = 0.98;
end

function opts = removeOffsetSearchOptionsLocal(opts)
if isfield(opts, "SSBOffsetSearch")
    opts = rmfield(opts, "SSBOffsetSearch");
end
end

function summary = disabledSearchSummaryLocal(meta)
summary = struct();
summary.Enable = false;
summary.SelectedOffsetHz = 0;
summary.SelectedSSBCenterHz = meta.CenterFrequencyHz;
summary.NumOffsetsTested = 1;
summary.Statement = "SSB offset search disabled; normal centered receiver chain used.";
end

function out = mergeStructsLocal(base, override)
out = base;
if isempty(override) || isempty(fieldnames(override))
    return;
end
names = fieldnames(override);
for k = 1:numel(names)
    out.(names{k}) = override.(names{k});
end
end
