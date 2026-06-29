function result = analyzeCapture(iq, meta, opts)
%analyzeCapture Analyze one B210 IQ capture without gNB-location input.
%
%   result = analyzeCapture(iq, meta, opts) runs the receiver chain for a
%   single capture: PSS search, SSS/PCI recovery, and per-PCI timing.
%   Optional TDD/special-slot feedback can be enabled, but it is auxiliary
%   only and is not used for the synchronization verdict. This function does
%   not produce a strict synchronization verdict because static transmit
%   timing and propagation delay cannot be separated from one internally
%   clocked capture without a survey geometry.

arguments
    iq (:,1) {mustBeNumeric}
    meta struct
    opts struct = struct()
end

[meta, metaReport] = validateMetadata(meta);
[analysisIQ, analysisTrim] = applyAnalysisTrimLocal(iq, meta, opts);

pssOpts = getSubOptionsLocal(opts, "PSS");
sssOpts = getSubOptionsLocal(opts, "SSS");
sicOpts = getSubOptionsLocal(opts, "SIC");
timingOpts = getSubOptionsLocal(opts, "Timing");
tddOpts = getSubOptionsLocal(opts, "TDD");
powerOpts = getSubOptionsLocal(opts, "Power");

if ~isfield(sicOpts, "PSS")
    sicOpts.PSS = pssOpts;
end
if ~isfield(sicOpts, "SSS")
    sicOpts.SSS = sssOpts;
end

[sssDetections, sicDebug] = detectCellsWithSIC(analysisIQ, meta, sicOpts);
[sssDetections, sicDebug] = restoreTrimmedSampleOriginLocal( ...
    sssDetections, sicDebug, analysisTrim, meta.SampleRateHz);
[pssCandidates, pssDebug, sssDebug] = firstStageDebugLocal(sicDebug);
[cellTiming, detectionTiming] = estimateCellTiming(sssDetections, meta, timingOpts, iq);
[auxFeedback, tddResults, slotMeasurements, specialResults, specialMeasurements] = ...
    optionalFeedbackLocal(iq, meta, cellTiming, tddOpts, opts);

result = struct();
result.Mode = "single_capture_analysis_no_gnb_locations";
result.Statement = "Single-capture analysis only. Strict synchronization verdicts require multi-position survey fitting because gNB locations are unavailable. TDD and special-slot checks, when enabled, are auxiliary feedback only.";
result.Metadata = meta;
result.MetadataReport = metaReport;
result.AnalysisTrim = analysisTrim;
result.PSSCandidates = pssCandidates;
result.PSSDebug = pssDebug;
result.SSSDetections = sssDetections;
result.SSSDebug = sssDebug;
result.SICDebug = sicDebug;
result.CellTiming = cellTiming;
result.DetectionTiming = detectionTiming;
result.AuxiliaryFeedback = auxFeedback;
result.TDDResults = tddResults;
result.SlotMeasurements = slotMeasurements;
result.SpecialResults = specialResults;
result.SpecialMeasurements = specialMeasurements;

if isfield(opts, "EnablePowerMeasurements") && opts.EnablePowerMeasurements
    if isfield(opts, "CaptureID")
        powerOpts.CaptureID = opts.CaptureID;
    end
    if isfield(opts, "RxX_m")
        powerOpts.RxX_m = opts.RxX_m;
    end
    if isfield(opts, "RxY_m")
        powerOpts.RxY_m = opts.RxY_m;
    end
    if isfield(opts, "RxZ_m")
        powerOpts.RxZ_m = opts.RxZ_m;
    end
    result.PowerMeasurements = estimateRssPriorsAndPlanStops("power", iq, meta, result, powerOpts);
end

end

function [analysisIQ, trim] = applyAnalysisTrimLocal(iq, meta, opts)
trimSeconds = getNumericOptionLocal(opts, ...
    ["TrimStartSeconds","AnalysisTrimStartSeconds"], 0);
if trimSeconds == 0
    trimSeconds = getNumericOptionLocal(opts, ...
        ["TrimStartMs","AnalysisTrimStartMs"], 0) * 1e-3;
end

if ~isfinite(trimSeconds) || trimSeconds < 0
    error("analyzeCapture:InvalidTrimStart", ...
        "TrimStartSeconds must be a non-negative finite scalar.");
end

trimSamples = round(trimSeconds * meta.SampleRateHz);
if trimSamples >= numel(iq)
    error("analyzeCapture:TrimTooLong", ...
        "TrimStartSeconds removes the entire capture.");
end

durationSeconds = getNumericOptionLocal(opts, ...
    ["AnalysisDurationSeconds","MaxAnalysisSeconds"], Inf);
if ~isfinite(durationSeconds)
    analysisEndSample = numel(iq);
elseif durationSeconds <= 0
    error("analyzeCapture:InvalidAnalysisDuration", ...
        "AnalysisDurationSeconds must be positive when specified.");
else
    analysisEndSample = min(numel(iq), ...
        trimSamples + round(durationSeconds*meta.SampleRateHz));
end
analysisIQ = iq(trimSamples+1:analysisEndSample);
trim = struct();
trim.TrimStartSeconds = trimSamples / meta.SampleRateHz;
trim.TrimStartMs = trim.TrimStartSeconds * 1e3;
trim.TrimStartSamples = trimSamples;
trim.OriginalNumSamples = numel(iq);
trim.AnalyzedNumSamples = numel(analysisIQ);
trim.AnalysisEndSample1 = analysisEndSample;
trim.AnalysisDurationSeconds = numel(analysisIQ)/meta.SampleRateHz;
trim.SampleOriginRestored = true;
if trimSamples > 0 && analysisEndSample < numel(iq)
    trim.Statement = "A bounded IQ interval was analyzed after startup trimming; reported sample indices were restored to the original capture origin.";
elseif trimSamples > 0
    trim.Statement = "Initial IQ samples were ignored during detection to avoid receiver startup/DC transients; reported sample indices were restored to the original capture origin.";
elseif analysisEndSample < numel(iq)
    trim.Statement = "A bounded IQ interval was analyzed; reported sample indices remain in the original capture origin.";
else
    trim.Statement = "No analysis trim was applied.";
end
end

function [detections, sicDebug] = restoreTrimmedSampleOriginLocal(detections, sicDebug, trim, fs)
offset = trim.TrimStartSamples;
if offset == 0
    return;
end

detections = addSampleOffsetToTableLocal(detections, offset, fs);
if isfield(sicDebug, "Stages")
    for k = 1:numel(sicDebug.Stages)
        if isfield(sicDebug.Stages{k}, "PSSCandidates")
            sicDebug.Stages{k}.PSSCandidates = addSampleOffsetToTableLocal( ...
                sicDebug.Stages{k}.PSSCandidates, offset, fs);
        end
        if isfield(sicDebug.Stages{k}, "SSSDetections")
            sicDebug.Stages{k}.SSSDetections = addSampleOffsetToTableLocal( ...
                sicDebug.Stages{k}.SSSDetections, offset, fs);
        end
    end
end
if isfield(sicDebug, "Cancellations")
    sicDebug.Cancellations = addSampleOffsetToTableLocal( ...
        sicDebug.Cancellations, offset, fs);
end
sicDebug.AnalysisTrim = trim;
end

function tbl = addSampleOffsetToTableLocal(tbl, offset, fs)
if ~istable(tbl) || isempty(tbl) || height(tbl) == 0
    return;
end

sampleFields = ["StartSample0","StartSample1","SegmentStartIndex"];
names = string(tbl.Properties.VariableNames);
for k = 1:numel(sampleFields)
    field = sampleFields(k);
    if ismember(field, names)
        tbl.(field) = tbl.(field) + offset;
    end
end
if ismember("StartTimeUs", names) && ismember("StartSample0", names)
    tbl.StartTimeUs = tbl.StartSample0 ./ fs * 1e6;
end
end

function value = getNumericOptionLocal(opts, names, defaultValue)
value = defaultValue;
for k = 1:numel(names)
    name = char(names(k));
    if isfield(opts, name)
        value = double(opts.(name));
        return;
    end
end
end

function [feedback, tddResults, slotMeasurements, specialResults, specialMeasurements] = optionalFeedbackLocal(iq, meta, cellTiming, tddOpts, opts)
enableAll = getLogicalOptionLocal(opts, ["EnableAuxiliaryFeedback","EnableTDDFeedback"], false);
enableTDD = enableAll || getLogicalOptionLocal(opts, "EnableTDDPatternFeedback", false);
enableSpecial = enableAll || getLogicalOptionLocal(opts, "EnableSpecialSlotFeedback", false);

feedback = struct();
feedback.Statement = "Auxiliary feedback only. TDD-pattern and special-slot energy checks are not used for synchronization PASS/FAIL decisions.";
feedback.UsedForSynchronizationVerdict = false;
feedback.TDDPatternFeedbackEnabled = enableTDD;
feedback.SpecialSlotFeedbackEnabled = enableSpecial;

if enableTDD
    [tddResults, slotMeasurements] = checkTDDPatternFeedback(iq, meta, cellTiming, tddOpts);
else
    tddResults = emptyTDDFeedbackResultsLocal();
    slotMeasurements = emptySlotFeedbackMeasurementsLocal();
end

if enableSpecial
    [specialResults, specialMeasurements] = checkSpecialSlotFeedback(iq, meta, cellTiming, tddOpts);
else
    specialResults = emptySpecialFeedbackResultsLocal();
    specialMeasurements = emptySpecialFeedbackMeasurementsLocal();
end
end

function tf = getLogicalOptionLocal(opts, names, defaultValue)
tf = defaultValue;
for k = 1:numel(names)
    name = char(names(k));
    if isfield(opts, name)
        tf = logical(opts.(name));
        return;
    end
end
end

function tbl = emptyTDDFeedbackResultsLocal()
tbl = table( ...
    zeros(0,1), strings(0,1), strings(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    false(0,1), strings(0,1), ...
    'VariableNames', ["PCI","TDDPatternStatus","Reason","NumFrames", ...
    "MedianDPower","MedianSPower","MedianUPower","DToURatioDb", ...
    "DToNoiseDb","NoiseFloorPower","UsedForTimingVerdict", ...
    "VerdictUseStatement"]);
end

function tbl = emptySlotFeedbackMeasurementsLocal()
tbl = table( ...
    zeros(0,1), zeros(0,1), zeros(0,1), strings(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', ["PCI","FrameNumber","SlotIndex", ...
    "ExpectedSlotType","StartSample1","EndSample1", ...
    "SlotPower","SlotPowerDb"]);
end

function tbl = emptySpecialFeedbackResultsLocal()
tbl = table( ...
    zeros(0,1), strings(0,1), strings(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), false(0,1), strings(0,1), ...
    'VariableNames', ["PCI","SpecialSlotStatus","Reason","NumSpecialSlots", ...
    "MedianDLPower","MedianTailPower","DLToTailRatioDb", ...
    "UsedForTimingVerdict","VerdictUseStatement"]);
end

function tbl = emptySpecialFeedbackMeasurementsLocal()
tbl = table( ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', ["PCI","FrameNumber","SlotIndex", ...
    "DLStartSample1","DLEndSample1","GuardStartSample1","GuardEndSample1", ...
    "ULStartSample1","ULEndSample1","DLPower","GuardPower", ...
    "ULPower","TailPower","DLToTailRatioDb"]);
end

function sub = getSubOptionsLocal(opts, name)
if isfield(opts, name)
    sub = opts.(name);
else
    sub = struct();
end
end

function [pssCandidates, pssDebug, sssDebug] = firstStageDebugLocal(sicDebug)
if isfield(sicDebug, "Stages") && ~isempty(sicDebug.Stages)
    first = sicDebug.Stages{1};
    pssCandidates = first.PSSCandidates;
    sssDebug = first.SSSDebug;
else
    pssCandidates = table();
    sssDebug = struct();
end

pssDebug = struct();
pssDebug.Statement = "PSS candidates shown for first residual-search stage. Full SIC stage data is in result.SICDebug.";
if isfield(sicDebug, "Stages")
    pssDebug.NumStages = numel(sicDebug.Stages);
else
    pssDebug.NumStages = 0;
end
end
