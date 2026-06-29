function [specialResults, specialMeasurements] = checkSpecialSlotFeedback(iq, meta, cellTiming, opts)
%checkSpecialSlotFeedback Optional non-verdict special-slot energy feedback.
%
%   This function checks whether configured special-slot DL/tail energy
%   behavior is broadly visible in the IQ envelope. It is auxiliary feedback
%   only and must not be used as a synchronization PASS/FAIL verdict.

arguments
    iq (:,1) {mustBeNumeric}
    meta struct
    cellTiming table
    opts struct = struct()
end

[meta, ~] = validateMetadata(meta);
base = defaultSpecialSlotFeedbackOptions(meta);
opts = mergeStructsLocal(base, opts);
opts = normalizeSpecialSlotFeedbackOptions(opts);

if isempty(cellTiming) || height(cellTiming) == 0
    specialResults = emptySpecialResults();
    specialMeasurements = emptySpecialMeasurements();
    return;
end

required = ["PCI","FramePhaseSamples"];
missing = required(~ismember(required, string(cellTiming.Properties.VariableNames)));
if ~isempty(missing)
    error("checkSpecialSlotFeedback:MissingColumns", ...
        "cellTiming table is missing required columns: %s", strjoin(missing, ", "));
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
    frameStarts0 = validFrameStartsLocal(framePhase, frameSamples, numel(power));

    if isempty(specialSlotIdx) || isempty(frameStarts0)
        resultCount = resultCount + 1;
        resultRows{resultCount} = specialResultRowLocal(pci, "NOT_ASSESSABLE", ...
            "No complete special-slot windows available for optional feedback.", ...
            0, NaN, NaN, NaN);
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

            [dlRange, guardRange, ulRange, tailRange] = specialRangesLocal( ...
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
        reason = "Optional feedback: special-slot DL portion is stronger than the guard/uplink tail.";
    elseif dlToTailRatioDb <= opts.FailSpecialDLToTailRatioDb
        status = "FAIL";
        reason = "Optional feedback: special-slot tail contains downlink-like energy comparable to the DL portion.";
    else
        status = "SUSPECT";
        reason = "Optional feedback: special-slot DL/tail energy separation is weak.";
    end

    resultCount = resultCount + 1;
    resultRows{resultCount} = specialResultRowLocal(pci, status, reason, ...
        numel(dlPowers), medianDL, medianTail, dlToTailRatioDb);
end

specialResults = vertcat(resultRows{1:resultCount});
if isempty(measurementRows)
    specialMeasurements = emptySpecialMeasurements();
else
    specialMeasurements = vertcat(measurementRows{:});
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

function [dlRange, guardRange, ulRange, tailRange] = ...
        specialRangesLocal(slotStart, slotEnd, symbols, meta, opts)
numDLSym = symbols(1);
numGuardSym = symbols(2);
numULSym = symbols(3);
if numDLSym + numGuardSym + numULSym ~= 14
    error("checkSpecialSlotFeedback:InvalidSpecialSlotSymbols", ...
        "SpecialSlotSymbols must sum to 14 normal-CP OFDM symbols.");
end

edges = symbolEdgesLocal(slotStart, slotEnd, meta, opts);
dlRange = edges(1):edges(numDLSym+1)-1;
guardStart = numDLSym + 1;
guardEnd = numDLSym + numGuardSym + 1;
ulStart = numDLSym + numGuardSym + 1;
ulEnd = numDLSym + numGuardSym + numULSym + 1;

guardRange = edges(guardStart):edges(guardEnd)-1;
ulRange = edges(ulStart):edges(ulEnd)-1;
tailRange = edges(guardStart):edges(ulEnd)-1;
end

function edges = symbolEdgesLocal(slotStart, slotEnd, meta, opts)
carrier = nrCarrierConfig;
carrier.NSizeGrid = opts.NSizeGrid;
carrier.SubcarrierSpacing = opts.SCSkHz;
info = nrOFDMInfo(carrier, "SampleRate", meta.SampleRateHz);
symbolsPerSlot = double(info.SymbolsPerSlot);
symbolLengths = double(info.SymbolLengths(1:symbolsPerSlot));
symbolLengths = symbolLengths(:);
if symbolsPerSlot ~= 14
    error("checkSpecialSlotFeedback:UnsupportedCyclicPrefix", ...
        "Special-slot feedback requires 14-symbol normal-CP slots.");
end

% Preserve the real CP-dependent symbol-length proportions while forcing the
% final edge to the already established integer-sample slot boundary.
slotSamples = slotEnd-slotStart+1;
relativeEdges = [0; cumsum(symbolLengths)] / sum(symbolLengths);
edges = slotStart + round(relativeEdges*slotSamples);
edges(end) = slotEnd+1;
edges = edges(:).';
end

function row = specialResultRowLocal(pci, status, reason, numSpecialSlots, medianDL, medianTail, dlToTailRatioDb)
row = table( ...
    pci, string(status), string(reason), numSpecialSlots, ...
    medianDL, medianTail, dlToTailRatioDb, ...
    false, string("Auxiliary feedback only; not used for synchronization verdict."), ...
    'VariableNames', ["PCI","SpecialSlotStatus","Reason","NumSpecialSlots", ...
    "MedianDLPower","MedianTailPower","DLToTailRatioDb", ...
    "UsedForTimingVerdict","VerdictUseStatement"]);
end

function tbl = emptySpecialResults()
tbl = table( ...
    zeros(0,1), strings(0,1), strings(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), false(0,1), strings(0,1), ...
    'VariableNames', ["PCI","SpecialSlotStatus","Reason","NumSpecialSlots", ...
    "MedianDLPower","MedianTailPower","DLToTailRatioDb", ...
    "UsedForTimingVerdict","VerdictUseStatement"]);
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

function opts = defaultSpecialSlotFeedbackOptions(meta)
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
% Auxiliary empirical gates; special-slot feedback never enters the verdict.
opts.MinSpecialDLToTailRatioDb = 3;
opts.FailSpecialDLToTailRatioDb = 1;
end

function opts = normalizeSpecialSlotFeedbackOptions(opts)
if isempty(opts.SlotPeriodMs) || isnan(opts.SlotPeriodMs)
    opts.SlotPeriodMs = 15 / opts.SCSkHz;
end
if opts.SlotPeriodMs <= 0 || ~isfinite(opts.SlotPeriodMs)
    error("checkSpecialSlotFeedback:InvalidSlotPeriod", ...
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
