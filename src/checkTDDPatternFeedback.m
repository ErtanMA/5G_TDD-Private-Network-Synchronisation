function [tddResults, slotMeasurements] = checkTDDPatternFeedback(iq, meta, cellTiming, opts)
%checkTDDPatternFeedback Optional non-verdict TDD energy feedback.
%
%   This function checks whether the observed IQ power envelope is broadly
%   consistent with one configured TDD pattern. It is auxiliary feedback only
%   and must not be used as a synchronization PASS/FAIL verdict.

arguments
    iq (:,1) {mustBeNumeric}
    meta struct
    cellTiming table
    opts struct = struct()
end

[meta, ~] = validateMetadata(meta);
base = defaultTDDPatternFeedbackOptions(meta);
opts = mergeStructsLocal(base, opts);
opts = normalizeTDDPatternFeedbackOptions(opts);

if isempty(cellTiming) || height(cellTiming) == 0
    tddResults = emptyTDDResults();
    slotMeasurements = emptySlotMeasurements();
    return;
end

required = ["PCI","FramePhaseSamples"];
missing = required(~ismember(required, string(cellTiming.Properties.VariableNames)));
if ~isempty(missing)
    error("checkTDDPatternFeedback:MissingColumns", ...
        "cellTiming table is missing required columns: %s", strjoin(missing, ", "));
end

fs = meta.SampleRateHz;
slotSamples = round(opts.SlotPeriodMs * 1e-3 * fs);
frameSamples = slotSamples * numel(opts.TDDPattern);
expectedFrameSamples = round(opts.FramePeriodMs * 1e-3 * fs);
if abs(frameSamples - expectedFrameSamples) > 1
    error("checkTDDPatternFeedback:InconsistentFrame", ...
        "TDD pattern length and frame period do not match.");
end

power = abs(iq(:)).^2;
smoothSamples = max(1, round(opts.EnvelopeSmoothUs * 1e-6 * fs));
if smoothSamples > 1
    power = movmean(power, smoothSamples);
end

noiseFloor = prctile(power, opts.NoisePercentile);

resultRows = cell(height(cellTiming), 1);
slotRows = {};
resultCount = 0;

for c = 1:height(cellTiming)
    pci = cellTiming.PCI(c);
    framePhase = cellTiming.FramePhaseSamples(c);
    frameStarts0 = validFrameStartsLocal(framePhase, frameSamples, numel(power));

    if numel(frameStarts0) < opts.MinFrames
        resultCount = resultCount + 1;
        resultRows{resultCount} = tddResultRowLocal(pci, "NOT_ASSESSABLE", ...
            "Not enough complete frames for optional TDD-pattern feedback.", ...
            numel(frameStarts0), NaN, NaN, NaN, NaN, noiseFloor);
        continue;
    end

    localSlotRows = cell(numel(frameStarts0) * numel(opts.TDDPattern), 1);
    localCount = 0;
    slotPowers = zeros(numel(frameStarts0), numel(opts.TDDPattern));

    for f = 1:numel(frameStarts0)
        frameStart0 = round(frameStarts0(f));
        for s = 1:numel(opts.TDDPattern)
            slotStart = frameStart0 + (s-1)*slotSamples + 1;
            slotEnd = slotStart + slotSamples - 1;
            slotPower = robustPowerLocal(power(slotStart:slotEnd), opts.PowerStatistic);
            slotPowers(f, s) = slotPower;

            localCount = localCount + 1;
            localSlotRows{localCount} = table( ...
                pci, f, s, opts.TDDPattern(s), slotStart, slotEnd, ...
                slotPower, pow2dbLocal(slotPower), ...
                'VariableNames', ["PCI","FrameNumber","SlotIndex", ...
                "ExpectedSlotType","StartSample1","EndSample1", ...
                "SlotPower","SlotPowerDb"]);
        end
    end

    slotRows = [slotRows; localSlotRows(1:localCount)]; %#ok<AGROW>

    dPowers = slotPowers(:, opts.TDDPattern == "D");
    uPowers = slotPowers(:, opts.TDDPattern == "U");
    sPowers = slotPowers(:, opts.TDDPattern == "S");

    medianD = median(dPowers(:));
    medianU = median(uPowers(:));
    medianS = median(sPowers(:));
    dToURatioDb = pow2dbLocal(medianD / (medianU + eps));
    dToNoiseDb = pow2dbLocal(medianD / (noiseFloor + eps));

    if dToNoiseDb < opts.MinDToNoiseDb
        status = "NOT_ASSESSABLE";
        reason = "Downlink energy is too close to the estimated noise floor.";
    elseif dToURatioDb >= opts.MinDToURatioDb
        status = "PASS";
        reason = "Optional feedback: expected D slots are stronger than expected U slots.";
    elseif dToURatioDb <= opts.FailDToURatioDb
        status = "FAIL";
        reason = "Optional feedback: expected U slots contain downlink-like energy comparable to D slots.";
    else
        status = "SUSPECT";
        reason = "Optional feedback: D/U energy separation is weak.";
    end

    resultCount = resultCount + 1;
    resultRows{resultCount} = tddResultRowLocal(pci, status, reason, ...
        numel(frameStarts0), medianD, medianS, medianU, dToURatioDb, noiseFloor);
    resultRows{resultCount}.DToNoiseDb = dToNoiseDb;
end

tddResults = vertcat(resultRows{1:resultCount});
if isempty(slotRows)
    slotMeasurements = emptySlotMeasurements();
else
    slotMeasurements = vertcat(slotRows{:});
end
end

function frameStarts0 = validFrameStartsLocal(framePhase, frameSamples, numSamples)
firstFrame = ceil((-framePhase) / frameSamples);
lastFrame = floor((numSamples - frameSamples - framePhase) / frameSamples);
if lastFrame < firstFrame
    frameStarts0 = zeros(0,1);
else
    frameStarts0 = framePhase + (firstFrame:lastFrame).' * frameSamples;
end
end

function p = robustPowerLocal(x, statistic)
switch string(statistic)
    case "mean"
        p = mean(x);
    case "median"
        p = median(x);
    otherwise
        error("checkTDDPatternFeedback:UnknownStatistic", ...
            "Unknown PowerStatistic: %s", statistic);
end
end

function row = tddResultRowLocal(pci, status, reason, numFrames, medianD, medianS, medianU, dToURatioDb, noiseFloor)
row = table( ...
    pci, string(status), string(reason), numFrames, ...
    medianD, medianS, medianU, dToURatioDb, NaN, noiseFloor, ...
    false, string("Auxiliary feedback only; not used for synchronization verdict."), ...
    'VariableNames', ["PCI","TDDPatternStatus","Reason","NumFrames", ...
    "MedianDPower","MedianSPower","MedianUPower","DToURatioDb", ...
    "DToNoiseDb","NoiseFloorPower","UsedForTimingVerdict", ...
    "VerdictUseStatement"]);
end

function tbl = emptyTDDResults()
tbl = table( ...
    zeros(0,1), strings(0,1), strings(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    false(0,1), strings(0,1), ...
    'VariableNames', ["PCI","TDDPatternStatus","Reason","NumFrames", ...
    "MedianDPower","MedianSPower","MedianUPower","DToURatioDb", ...
    "DToNoiseDb","NoiseFloorPower","UsedForTimingVerdict", ...
    "VerdictUseStatement"]);
end

function tbl = emptySlotMeasurements()
tbl = table( ...
    zeros(0,1), zeros(0,1), zeros(0,1), strings(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', ["PCI","FrameNumber","SlotIndex", ...
    "ExpectedSlotType","StartSample1","EndSample1", ...
    "SlotPower","SlotPowerDb"]);
end

function opts = defaultTDDPatternFeedbackOptions(meta)
[meta, ~] = validateMetadata(meta);

opts = struct();
opts.SampleRateHz = meta.SampleRateHz;
opts.SCSkHz = 30;
opts.FramePeriodMs = 10;
opts.SlotPeriodMs = NaN;
opts.TDDPattern = ["D","D","D","S","U","D","D","D","S","U", ...
                   "D","D","D","S","U","D","D","D","S","U"];
opts.EnvelopeSmoothUs = 2;
opts.PowerStatistic = "mean";
opts.NoisePercentile = 10;
opts.MinFrames = 1;
% Auxiliary empirical energy-separation gates, not 3GPP conformance limits.
opts.MinDToURatioDb = 3;
opts.FailDToURatioDb = 1;
opts.MinDToNoiseDb = 1;
end

function opts = normalizeTDDPatternFeedbackOptions(opts)
if isempty(opts.SlotPeriodMs) || isnan(opts.SlotPeriodMs)
    opts.SlotPeriodMs = 15 / opts.SCSkHz;
end
if opts.SlotPeriodMs <= 0 || ~isfinite(opts.SlotPeriodMs)
    error("checkTDDPatternFeedback:InvalidSlotPeriod", ...
        "SlotPeriodMs must be positive finite, or NaN for SCS-derived default.");
end
end

function y = pow2dbLocal(x)
y = 10*log10(x + eps);
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
