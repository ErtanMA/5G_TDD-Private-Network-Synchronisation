function measurementTable = timingMeasurementsFromCaptureResult(captureResult, captureID)
%timingMeasurementsFromCaptureResult Convert one capture result to survey rows.

arguments
    captureResult struct
    captureID (1,1) string
end

if ~isfield(captureResult, "CellTiming") || isempty(captureResult.CellTiming) || ...
        height(captureResult.CellTiming) == 0
    measurementTable = emptyMeasurementTableLocal();
    return;
end

cellTiming = captureResult.CellTiming;
num = height(cellTiming);

if isfield(captureResult, "Metadata") && isfield(captureResult.Metadata, "CenterFrequencyHz")
    centerFrequencyHz = captureResult.Metadata.CenterFrequencyHz;
else
    centerFrequencyHz = NaN;
end

[ssbOffsetHz, ssbCenterHz] = getSSBOffsetSearchInfoLocal(captureResult, centerFrequencyHz);
[gscn, cellKey, rowSSBCenterHz, rowOffsetHz] = ...
    getFrequencyIdentityLocal(cellTiming, ssbCenterHz, ssbOffsetHz);

measurementTable = table( ...
    repmat(captureID, num, 1), ...
    cellKey, ...
    gscn, ...
    cellTiming.PCI, ...
    cellTiming.FramePhaseNs, ...
    cellTiming.TimingStdNs, ...
    repmat(centerFrequencyHz, num, 1), ...
    rowOffsetHz, ...
    rowSSBCenterHz, ...
    cellTiming.MeanPSSMetric, ...
    cellTiming.MeanSSSMetric, ...
    cellTiming.MeanCFOHz, ...
    'VariableNames', ["CaptureID","CellKey","GSCN","PCI", ...
    "FramePhaseNs","TimingStdNs", ...
    "CenterFrequencyHz","SelectedSSBOffsetHz","HypothesizedSSBCenterHz", ...
    "MeanPSSMetric","MeanSSSMetric","MeanCFOHz"]);

end

function tbl = emptyMeasurementTableLocal()
tbl = table( ...
    strings(0,1), strings(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', ["CaptureID","CellKey","GSCN","PCI", ...
    "FramePhaseNs","TimingStdNs", ...
    "CenterFrequencyHz","SelectedSSBOffsetHz","HypothesizedSSBCenterHz", ...
    "MeanPSSMetric","MeanSSSMetric","MeanCFOHz"]);
end

function [gscn,cellKey,ssbCenterHz,offsetHz] = ...
    getFrequencyIdentityLocal(cellTiming,defaultCenterHz,defaultOffsetHz)
num = height(cellTiming);
names = string(cellTiming.Properties.VariableNames);
if ismember("GSCN",names)
    gscn = cellTiming.GSCN;
else
    gscn = NaN(num,1);
end
if ismember("CellKey",names)
    cellKey = string(cellTiming.CellKey);
else
    cellKey = "PCI"+string(cellTiming.PCI);
end
if ismember("SSBCenterHz",names)
    ssbCenterHz = cellTiming.SSBCenterHz;
else
    ssbCenterHz = repmat(defaultCenterHz,num,1);
end
if ismember("DigitalShiftHz",names)
    offsetHz = cellTiming.DigitalShiftHz;
else
    offsetHz = repmat(defaultOffsetHz,num,1);
end
end

function [ssbOffsetHz, ssbCenterHz] = getSSBOffsetSearchInfoLocal(captureResult, centerFrequencyHz)
ssbOffsetHz = 0;
ssbCenterHz = centerFrequencyHz;
if isfield(captureResult, "SSBOffsetSearch")
    search = captureResult.SSBOffsetSearch;
    if isfield(search, "SelectedOffsetHz")
        ssbOffsetHz = search.SelectedOffsetHz;
    end
    if isfield(search, "SelectedSSBCenterHz")
        ssbCenterHz = search.SelectedSSBCenterHz;
    else
        ssbCenterHz = centerFrequencyHz + ssbOffsetHz;
    end
end
end
