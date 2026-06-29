function [tddResults, slotMeasurements, specialResults, specialMeasurements] = checkTDDPattern(iq, meta, cellTiming, opts)
%checkTDDPattern Check expected TDD slot energy pattern from IQ envelope.
%
%   The function uses per-PCI frame timing from estimateCellTiming to align
%   the power envelope to the expected 10 ms frame and 0.5 ms slot grid.
%   It reports whether the observed energy follows the configured TDD pattern
%   and separately returns special-slot analysis.

arguments
    iq (:,1) {mustBeNumeric}
    meta struct
    cellTiming table
    opts struct = struct()
end

[meta, ~] = validateMetadata(meta);
base = defaultTDDPatternOptions(meta);
opts = mergeStructsLocal(base, opts);
opts = normalizeTDDPatternOptionsLocal(opts);

if isempty(cellTiming) || height(cellTiming) == 0
    tddResults = emptyTDDResults();
    slotMeasurements = emptySlotMeasurements();
    [specialResults, specialMeasurements] = checkSpecialSlot(iq, meta, cellTiming, opts);
    return;
end

required = ["PCI","FramePhaseSamples"];
missing = required(~ismember(required, string(cellTiming.Properties.VariableNames)));
if ~isempty(missing)
    error("checkTDDPattern:MissingColumns", ...
        "cellTiming table is missing required columns: %s", strjoin(missing, ", "));
end

fs = meta.SampleRateHz;
slotSamples = round(opts.SlotPeriodMs * 1e-3 * fs);
frameSamples = slotSamples * numel(opts.TDDPattern);
expectedFrameSamples = round(opts.FramePeriodMs * 1e-3 * fs);
if abs(frameSamples - expectedFrameSamples) > 1
    error("checkTDDPattern:InconsistentFrame", ...
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
    frameStarts0 = validFrameStarts(framePhase, frameSamples, numel(power));

    if numel(frameStarts0) < opts.MinFrames
        resultCount = resultCount + 1;
        resultRows{resultCount} = tddResultRow(pci, "NOT_ASSESSABLE", ...
            "Not enough complete frames for TDD-pattern analysis.", ...
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
            slotPower = robustPower(power(slotStart:slotEnd), opts.PowerStatistic);
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
        reason = "Expected D slots are clearly stronger than expected U slots.";
    elseif dToURatioDb <= opts.FailDToURatioDb
        status = "FAIL";
        reason = "Expected U slots contain downlink-like energy comparable to D slots.";
    else
        status = "SUSPECT";
        reason = "D/U energy separation is weak.";
    end

    resultCount = resultCount + 1;
    resultRows{resultCount} = tddResultRow(pci, status, reason, ...
        numel(frameStarts0), medianD, medianS, medianU, dToURatioDb, noiseFloor);
    resultRows{resultCount}.DToNoiseDb = dToNoiseDb;
end

tddResults = vertcat(resultRows{1:resultCount});
if isempty(slotRows)
    slotMeasurements = emptySlotMeasurements();
else
    slotMeasurements = vertcat(slotRows{:});
end

[specialResults, specialMeasurements] = checkSpecialSlot(iq, meta, cellTiming, opts);

end

function frameStarts0 = validFrameStarts(framePhase, frameSamples, numSamples)
firstFrame = ceil((-framePhase) / frameSamples);
lastFrame = floor((numSamples - frameSamples - framePhase) / frameSamples);
if lastFrame < firstFrame
    frameStarts0 = zeros(0,1);
else
    frameStarts0 = framePhase + (firstFrame:lastFrame).' * frameSamples;
end
end

function p = robustPower(x, statistic)
switch string(statistic)
    case "mean"
        p = mean(x);
    case "median"
        p = median(x);
    otherwise
        error("checkTDDPattern:UnknownStatistic", ...
            "Unknown PowerStatistic: %s", statistic);
end
end

function row = tddResultRow(pci, status, reason, numFrames, medianD, medianS, medianU, dToURatioDb, noiseFloor)
row = table( ...
    pci, string(status), string(reason), numFrames, ...
    medianD, medianS, medianU, dToURatioDb, NaN, noiseFloor, ...
    'VariableNames', ["PCI","TDDPatternStatus","Reason","NumFrames", ...
    "MedianDPower","MedianSPower","MedianUPower","DToURatioDb", ...
    "DToNoiseDb","NoiseFloorPower"]);
end

function tbl = emptyTDDResults()
tbl = table( ...
    zeros(0,1), strings(0,1), strings(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', ["PCI","TDDPatternStatus","Reason","NumFrames", ...
    "MedianDPower","MedianSPower","MedianUPower","DToURatioDb", ...
    "DToNoiseDb","NoiseFloorPower"]);
end

function tbl = emptySlotMeasurements()
tbl = table( ...
    zeros(0,1), zeros(0,1), zeros(0,1), strings(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', ["PCI","FrameNumber","SlotIndex", ...
    "ExpectedSlotType","StartSample1","EndSample1", ...
    "SlotPower","SlotPowerDb"]);
end

function [specialResults, specialMeasurements] = checkSpecialSlot(iq, meta, cellTiming, opts)
% Check the energy split inside expected special slots. This is kept local to
% TDD checking because it uses the same frame/slot alignment and envelope.
[meta, ~] = validateMetadata(meta);
base = defaultTDDPatternOptions(meta);
opts = mergeStructsLocal(base, opts);
opts = normalizeTDDPatternOptionsLocal(opts);

if isempty(cellTiming) || height(cellTiming) == 0
    specialResults = emptySpecialResults();
    specialMeasurements = emptySpecialMeasurements();
    return;
end

fs = meta.SampleRateHz;
slotSamples = round(opts.SlotPeriodMs * 1e-3 * fs);
frameSamples = slotSamples * numel(opts.TDDPattern);
specialSlotIdx = find(opts.TDDPattern == "S");

power = abs(iq(:)).^2;
smoothSamples = max(1, round(opts.EnvelopeSmoothUs * 1e-6 * fs));
if smoothSamples > 1
    power = movmean(power, smoothSamples);
end

resultRows = cell(height(cellTiming), 1);
measurementRows = {};
resultCount = 0;

for c = 1:height(cellTiming)
    pci = cellTiming.PCI(c);
    framePhase = cellTiming.FramePhaseSamples(c);
    frameStarts0 = validFrameStarts(framePhase, frameSamples, numel(power));

    if isempty(specialSlotIdx) || isempty(frameStarts0)
        resultCount = resultCount + 1;
        resultRows{resultCount} = specialResultRow(pci, "NOT_ASSESSABLE", ...
            "No complete special-slot windows available.", 0, NaN, NaN, NaN);
        continue;
    end

    dlPowers = [];
    tailPowers = [];
    localRows = {};
    localCount = 0;

    for f = 1:numel(frameStarts0)
        frameStart0 = round(frameStarts0(f));
        for s = specialSlotIdx(:).'
            slotStart = frameStart0 + (s-1)*slotSamples + 1;
            slotEnd = slotStart + slotSamples - 1;

            [dlRange, guardRange, ulRange, tailRange] = specialRanges( ...
                slotStart, slotEnd, opts.SpecialSlotSymbols, meta, opts);
            dlPower = mean(power(dlRange));
            guardPower = mean(power(guardRange));
            ulPower = mean(power(ulRange));
            tailPower = mean(power(tailRange));

            dlPowers(end+1,1) = dlPower; %#ok<AGROW>
            tailPowers(end+1,1) = tailPower; %#ok<AGROW>

            localCount = localCount + 1;
            localRows{localCount,1} = table( ...
                pci, f, s, dlRange(1), dlRange(end), guardRange(1), guardRange(end), ...
                ulRange(1), ulRange(end), dlPower, guardPower, ulPower, tailPower, ...
                pow2dbLocal(dlPower / (tailPower + eps)), ...
                'VariableNames', ["PCI","FrameNumber","SlotIndex", ...
                "DLStartSample1","DLEndSample1","GuardStartSample1","GuardEndSample1", ...
                "ULStartSample1","ULEndSample1","DLPower","GuardPower", ...
                "ULPower","TailPower","DLToTailRatioDb"]);
        end
    end

    measurementRows = [measurementRows; localRows]; %#ok<AGROW>

    medianDL = median(dlPowers);
    medianTail = median(tailPowers);
    dlToTailRatioDb = pow2dbLocal(medianDL / (medianTail + eps));

    if dlToTailRatioDb >= opts.MinSpecialDLToTailRatioDb
        status = "PASS";
        reason = "Special-slot DL portion is clearly stronger than guard/uplink tail.";
    elseif dlToTailRatioDb <= opts.FailSpecialDLToTailRatioDb
        status = "FAIL";
        reason = "Special-slot tail contains downlink-like energy comparable to the DL portion.";
    else
        status = "SUSPECT";
        reason = "Special-slot DL/tail energy separation is weak.";
    end

    resultCount = resultCount + 1;
    resultRows{resultCount} = specialResultRow(pci, status, reason, ...
        numel(dlPowers), medianDL, medianTail, dlToTailRatioDb);
end

specialResults = vertcat(resultRows{1:resultCount});
if isempty(measurementRows)
    specialMeasurements = emptySpecialMeasurements();
else
    specialMeasurements = vertcat(measurementRows{:});
end
end

function [dlRange, guardRange, ulRange, tailRange] = ...
        specialRanges(slotStart, slotEnd, symbols, meta, opts)
numDLSym = symbols(1);
numGuardSym = symbols(2);
numULSym = symbols(3);

if numDLSym + numGuardSym + numULSym ~= 14
    error("checkTDDPattern:InvalidSpecialSlotSymbols", ...
        "SpecialSlotSymbols must sum to 14 normal-CP OFDM symbols.");
end

carrier = nrCarrierConfig;
carrier.NSizeGrid = opts.NSizeGrid;
carrier.SubcarrierSpacing = opts.SCSkHz;
info = nrOFDMInfo(carrier, "SampleRate", meta.SampleRateHz);
symbolsPerSlot = double(info.SymbolsPerSlot);
symbolLengths = double(info.SymbolLengths(1:symbolsPerSlot));
symbolLengths = symbolLengths(:);
if symbolsPerSlot ~= 14
    error("checkTDDPattern:UnsupportedCyclicPrefix", ...
        "Special-slot checking requires 14-symbol normal-CP slots.");
end
slotSamples = slotEnd-slotStart+1;
relativeEdges = [0; cumsum(symbolLengths)] / sum(symbolLengths);
edges = slotStart + round(relativeEdges*slotSamples);
edges(end) = slotEnd+1;
edges = edges(:).';

dlRange = edges(1):edges(numDLSym+1)-1;
guardStart = numDLSym + 1;
guardEnd = numDLSym + numGuardSym + 1;
ulStart = numDLSym + numGuardSym + 1;
ulEnd = numDLSym + numGuardSym + numULSym + 1;

guardRange = edges(guardStart):edges(guardEnd)-1;
ulRange = edges(ulStart):edges(ulEnd)-1;
tailRange = edges(guardStart):edges(ulEnd)-1;
end

function row = specialResultRow(pci, status, reason, numSpecialSlots, medianDL, medianTail, dlToTailRatioDb)
row = table( ...
    pci, string(status), string(reason), numSpecialSlots, ...
    medianDL, medianTail, dlToTailRatioDb, ...
    'VariableNames', ["PCI","SpecialSlotStatus","Reason","NumSpecialSlots", ...
    "MedianDLPower","MedianTailPower","DLToTailRatioDb"]);
end

function tbl = emptySpecialResults()
tbl = table( ...
    zeros(0,1), strings(0,1), strings(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', ["PCI","SpecialSlotStatus","Reason","NumSpecialSlots", ...
    "MedianDLPower","MedianTailPower","DLToTailRatioDb"]);
end

function tbl = emptySpecialMeasurements()
tbl = table( ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', ["PCI","FrameNumber","SlotIndex", ...
    "DLStartSample1","DLEndSample1","GuardStartSample1","GuardEndSample1", ...
    "ULStartSample1","ULEndSample1","DLPower","GuardPower", ...
    "ULPower","TailPower","DLToTailRatioDb"]);
end

function opts = defaultTDDPatternOptions(meta)
% Defaults for envelope-based DDDSU and special-slot checks.
[meta, ~] = validateMetadata(meta);

opts = struct();
opts.SampleRateHz = meta.SampleRateHz;
opts.SCSkHz = 30;
opts.NSizeGrid = 20;
opts.FramePeriodMs = 10;
opts.SlotPeriodMs = NaN;
opts.TDDPattern = ["D","D","D","S","U","D","D","D","S","U", ...
                   "D","D","D","S","U","D","D","D","S","U"];
opts.SpecialSlotSymbols = [10 2 2];
opts.EnvelopeSmoothUs = 2;
opts.PowerStatistic = "mean";
opts.NoisePercentile = 10;
opts.MinFrames = 1;
opts.MinDToURatioDb = 3;
opts.FailDToURatioDb = 1;
opts.MinDToNoiseDb = 1;
opts.MinSpecialDLToTailRatioDb = 3;
opts.FailSpecialDLToTailRatioDb = 1;
end

function opts = normalizeTDDPatternOptionsLocal(opts)
% NaN SlotPeriodMs means derive NR slot duration from numerology:
% 15 kHz -> 1 ms, 30 kHz -> 0.5 ms, 60 kHz -> 0.25 ms.
if isempty(opts.SlotPeriodMs) || isnan(opts.SlotPeriodMs)
    opts.SlotPeriodMs = 15 / opts.SCSkHz;
end
if opts.SlotPeriodMs <= 0 || ~isfinite(opts.SlotPeriodMs)
    error("checkTDDPattern:InvalidSlotPeriod", ...
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
