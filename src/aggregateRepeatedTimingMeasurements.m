function result = aggregateRepeatedTimingMeasurements(windowMeasurements, repeatManifest, opts)
%aggregateRepeatedTimingMeasurements Collapse stationary repeats by position.
%
% The input contains one timing row per detected cell and analysis window.
% Relative timing is formed inside each window, so the arbitrary B210 phase
% cancels before observations from different captures are combined.

arguments
    windowMeasurements table
    repeatManifest table
    opts struct = struct()
end

opts = defaultsLocal(opts);
validateInputsLocal(windowMeasurements,repeatManifest);
positionInfo = aggregatePositionInfoLocal(repeatManifest,opts);

measurements = windowMeasurements;
measurements.CaptureID = string(measurements.CaptureID);
measurements.PositionID = string(measurements.PositionID);
measurements.ObservationID = string(measurements.ObservationID);
measurements.CellKey = string(measurements.CellKey);

[measurements,rejectedMeasurements,reliableCellSupport] = ...
    qualityFilterLocal(measurements,opts);
if isempty(measurements) || height(measurements)==0
    result = emptyQualityResultLocal( ...
        windowMeasurements,rejectedMeasurements,reliableCellSupport, ...
        positionInfo);
    return;
end

[referenceCellKey,referencePCI,referenceGSCN,referenceRanking] = ...
    chooseReferenceLocal(measurements,reliableCellSupport,opts);
relative = formRelativeObservationsLocal( ...
    measurements,referenceCellKey,referencePCI,opts);
if isempty(relative)
    result = emptyRelativeResultLocal( ...
        measurements,positionInfo,referenceCellKey,referencePCI, ...
        referenceGSCN,rejectedMeasurements,reliableCellSupport, ...
        referenceRanking);
    return;
end
[perPosition,relative] = summarizePerPositionLocal(relative,opts);
stationary = summarizePerCellLocal(relative,opts);
aggregatedMeasurements = buildGeometryMeasurementsLocal( ...
    perPosition,referenceCellKey,referencePCI,referenceGSCN,opts);

result = struct();
result.Mode = "stationary_repeats_aggregated_by_position";
result.Statement = ...
    "Repeated captures are reference-subtracted within each stationary window, summarized as temporal relative timing instability, then collapsed to one geometry observation per receiver position.";
result.ReferenceCellKey = referenceCellKey;
result.ReferencePCI = referencePCI;
result.ReferenceGSCN = referenceGSCN;
result.RawWindowMeasurements = windowMeasurements;
result.WindowMeasurements = measurements;
result.RejectedWindowMeasurements = rejectedMeasurements;
result.ReliableCellSupport = reliableCellSupport;
result.ReferenceCandidateRanking = referenceRanking;
result.RelativeObservations = relative;
result.PerPositionTiming = perPosition;
result.StationaryInstability = stationary;
result.AggregatedMeasurementTable = aggregatedMeasurements;
result.PositionInfo = positionInfo;
end

function opts = defaultsLocal(opts)
defaults = struct();
defaults.FramePeriodNs = 10e6;
defaults.ReferencePCI = NaN;
defaults.ReferenceGSCN = NaN;
defaults.MinRepeatsPerPosition = 3;
defaults.MinObservationsPerPositionCell = 4;
defaults.MinRepeatsPerPositionCell = 2;
defaults.TimingUncertaintyFloorNs = 10;
defaults.DefaultPositionUncertaintyM = 5;
defaults.MaxTimingStdNs = 1000;
defaults.MinMeanSSSMetric = 0.06;
defaults.MinWindowsPerCaptureCell = 2;
defaults.MinCapturesPerPositionCell = 2;
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts,names{k})
        opts.(names{k}) = defaults.(names{k});
    end
end
end

function validateInputsLocal(measurements,manifest)
requiredMeasurements = ["CaptureID","PositionID","RepeatIndex", ...
    "ObservationID","CellKey","GSCN","PCI","FramePhaseNs", ...
    "TimingStdNs","MeanPSSMetric","MeanSSSMetric","MeanCFOHz"];
missing = requiredMeasurements(~ismember( ...
    requiredMeasurements,string(measurements.Properties.VariableNames)));
if ~isempty(missing)
    error("aggregateRepeatedTimingMeasurements:MissingMeasurementColumns", ...
        "windowMeasurements is missing: %s",strjoin(missing,", "));
end

requiredManifest = ["CaptureID","PositionID","RepeatIndex", ...
    "RxLat","RxLon","RxAltM","RxPositionUncertaintyM"];
missing = requiredManifest(~ismember( ...
    requiredManifest,string(manifest.Properties.VariableNames)));
if ~isempty(missing)
    error("aggregateRepeatedTimingMeasurements:MissingManifestColumns", ...
        "repeatManifest is missing: %s",strjoin(missing,", "));
end
if isempty(measurements) || height(measurements)==0
    error("aggregateRepeatedTimingMeasurements:NoMeasurements", ...
        "No stationary-window timing measurements were supplied.");
end
end

function info = aggregatePositionInfoLocal(manifest,opts)
manifest.CaptureID = string(manifest.CaptureID);
manifest.PositionID = string(manifest.PositionID);
positions = unique(manifest.PositionID,"stable");
parts = cell(numel(positions),1);
for k = 1:numel(positions)
    positionID = positions(k);
    rows = manifest(manifest.PositionID==positionID,:);
    repeats = unique(rows.RepeatIndex);
    if numel(repeats)<opts.MinRepeatsPerPosition
        error("aggregateRepeatedTimingMeasurements:TooFewRepeats", ...
            "Position %s has %d repeat(s); at least %d are required.", ...
            positionID,numel(repeats),opts.MinRepeatsPerPosition);
    end
    if numel(repeats)~=height(rows)
        error("aggregateRepeatedTimingMeasurements:DuplicateRepeatIndex", ...
            "Position %s contains duplicate RepeatIndex values.",positionID);
    end

    lat = median(double(rows.RxLat),"omitnan");
    lon = median(double(rows.RxLon),"omitnan");
    alt = median(double(rows.RxAltM),"omitnan");
    [dx,dy] = localOffsetsLocal(double(rows.RxLat),double(rows.RxLon),lat,lon);
    gpsSpreadM = max(hypot(dx,dy),[],"omitnan");
    reportedUncertainty = median( ...
        double(rows.RxPositionUncertaintyM),"omitnan");
    uncertaintyM = max([opts.DefaultPositionUncertaintyM, ...
        reportedUncertainty,gpsSpreadM]);

    parts{k} = table(positionID,numel(repeats),lat,lon,alt, ...
        uncertaintyM,gpsSpreadM, ...
        'VariableNames',["CaptureID","NumRepeats","RxLat","RxLon", ...
        "RxAltM","RxPositionUncertaintyM","RepeatGpsSpreadM"]);
end
info = vertcat(parts{:});
end

function [accepted,rejected,support] = qualityFilterLocal(measurements,opts)
reason = strings(height(measurements),1);
invalidNumeric = ~isfinite(measurements.FramePhaseNs) | ...
    ~isfinite(measurements.TimingStdNs) | ...
    ~isfinite(measurements.MeanSSSMetric);
reason(invalidNumeric) = "NONFINITE_TIMING_OR_QUALITY";

tooUncertain = reason=="" & ...
    measurements.TimingStdNs>opts.MaxTimingStdNs;
reason(tooUncertain) = "TIMING_SPREAD_ABOVE_LIMIT";

weakSSS = reason=="" & ...
    measurements.MeanSSSMetric<opts.MinMeanSSSMetric;
reason(weakSSS) = "SSS_METRIC_BELOW_LIMIT";

provisional = measurements(reason=="",:);
captureKeys = unique(provisional(:,["CaptureID","CellKey"]),"rows","stable");
for k = 1:height(captureKeys)
    idx = measurements.CaptureID==captureKeys.CaptureID(k) & ...
        measurements.CellKey==captureKeys.CellKey(k) & reason=="";
    if numel(unique(measurements.ObservationID(idx))) < ...
            opts.MinWindowsPerCaptureCell
        reason(idx) = "TOO_FEW_WINDOWS_IN_CAPTURE";
    end
end

provisional = measurements(reason=="",:);
positionKeys = unique(provisional(:,["PositionID","CellKey"]),"rows","stable");
for k = 1:height(positionKeys)
    idx = measurements.PositionID==positionKeys.PositionID(k) & ...
        measurements.CellKey==positionKeys.CellKey(k) & reason=="";
    if numel(unique(measurements.CaptureID(idx))) < ...
            opts.MinCapturesPerPositionCell
        reason(idx) = "TOO_FEW_CAPTURES_AT_POSITION";
    end
end

accepted = measurements(reason=="",:);
rejected = measurements(reason~="",:);
rejected.RejectionReason = reason(reason~="");
support = reliableCellSupportLocal(accepted);
end

function support = reliableCellSupportLocal(measurements)
if isempty(measurements) || height(measurements)==0
    support = emptyReliableCellSupportLocal();
    return;
end
cells = unique(measurements(:,["CellKey","GSCN","PCI"]),"rows","stable");
parts = cell(height(cells),1);
for k = 1:height(cells)
    rows = measurements(measurements.CellKey==cells.CellKey(k),:);
    positionCapture = unique(rows(:,["PositionID","CaptureID"]),"rows");
    parts{k} = table( ...
        cells.CellKey(k),cells.GSCN(k),cells.PCI(k), ...
        numel(unique(rows.PositionID)),height(positionCapture), ...
        numel(unique(rows.ObservationID)),height(rows), ...
        mean(rows.MeanSSSMetric,"omitnan"), ...
        median(rows.TimingStdNs,"omitnan"), ...
        'VariableNames',["CellKey","GSCN","PCI", ...
        "ReliablePositions","ReliableCaptures","ReliableObservations", ...
        "NumTimingRows","MeanSSSMetric","MedianTimingStdNs"]);
end
support = vertcat(parts{:});
end

function [cellKey,pci,gscn,ranking] = chooseReferenceLocal( ...
        measurements,support,opts)
ranking = sortrows(support, ...
    ["ReliablePositions","ReliableCaptures","ReliableObservations", ...
    "MeanSSSMetric","MedianTimingStdNs"], ...
    ["descend","descend","descend","descend","ascend"]);
if isfinite(opts.ReferencePCI)
    candidates = measurements(measurements.PCI==opts.ReferencePCI,:);
    if isfinite(opts.ReferenceGSCN)
        candidates = candidates(candidates.GSCN==opts.ReferenceGSCN,:);
    end
    if isempty(candidates)
        error("aggregateRepeatedTimingMeasurements:ReferenceMissing", ...
            "Configured reference PCI/GSCN was not detected.");
    end
    keys = unique(candidates.CellKey,"stable");
    cellKey = keys(1);
else
    cellKey = ranking.CellKey(1);
end
row = measurements(find(measurements.CellKey==cellKey,1),:);
pci = row.PCI;
gscn = row.GSCN;
end

function relative = formRelativeObservationsLocal(measurements,referenceCellKey,referencePCI,opts)
observationIDs = unique(measurements.ObservationID,"stable");
parts = {};
count = 0;
for o = 1:numel(observationIDs)
    obsID = observationIDs(o);
    rows = measurements(measurements.ObservationID==obsID,:);
    refRows = rows(rows.CellKey==referenceCellKey,:);
    if isempty(refRows)
        continue;
    end
    ref = bestRowLocal(refRows);
    targets = rows(rows.CellKey~=referenceCellKey,:);
    for k = 1:height(targets)
        target = targets(k,:);
        count = count+1;
        relativeNs = wrapLocal( ...
            target.FramePhaseNs-ref.FramePhaseNs,opts.FramePeriodNs);
        sigmaNs = hypot(target.TimingStdNs,ref.TimingStdNs);
        parts{count,1} = table( ...
            target.PositionID,target.CaptureID,target.RepeatIndex, ...
            target.ObservationID,target.CellKey,target.GSCN,target.PCI, ...
            referenceCellKey,referencePCI,relativeNs,sigmaNs, ...
            target.MeanPSSMetric,target.MeanSSSMetric,target.MeanCFOHz, ...
            NaN, ...
            'VariableNames',["PositionID","CaptureID","RepeatIndex", ...
            "ObservationID","CellKey","GSCN","PCI","ReferenceCellKey", ...
            "ReferencePCI","RelativeArrivalNs","SigmaNs", ...
            "MeanPSSMetric","MeanSSSMetric","MeanCFOHz", ...
            "CenteredResidualNs"]); %#ok<AGROW>
    end
end
if count==0
    relative = emptyRelativeObservationsLocal();
    return;
end
relative = vertcat(parts{:});
end

function [summary,relative] = summarizePerPositionLocal(relative,opts)
keys = unique(relative(:,["PositionID","CellKey"]),"rows","stable");
parts = cell(height(keys),1);
keep = false(height(keys),1);
for k = 1:height(keys)
    idx = relative.PositionID==keys.PositionID(k) & ...
        relative.CellKey==keys.CellKey(k);
    rows = relative(idx,:);
    numRepeats = numel(unique(rows.RepeatIndex));
    if height(rows)<opts.MinObservationsPerPositionCell || ...
            numRepeats<opts.MinRepeatsPerPositionCell
        continue;
    end
    [center,residual] = robustCircularCenterLocal( ...
        rows.RelativeArrivalNs,opts.FramePeriodNs);
    relative.CenteredResidualNs(idx) = residual;
    rmsNs = sqrt(mean(residual.^2,"omitnan"));
    ppNs = max(residual)-min(residual);
    maxAbsNs = max(abs(residual));
    expectedNoiseNs = sqrt(mean(rows.SigmaNs.^2,"omitnan"));
    excessRmsNs = sqrt(max(0,rmsNs^2-expectedNoiseNs^2));
    robustStdNs = 1.4826*median(abs(residual-median(residual)),"omitnan");
    effectiveIndependentSamples = max(1,numRepeats);
    centerUncertaintyNs = max([opts.TimingUncertaintyFloorNs, ...
        expectedNoiseNs/sqrt(effectiveIndependentSamples), ...
        robustStdNs/sqrt(effectiveIndependentSamples)]);
    meanPSSMetric = mean(rows.MeanPSSMetric,"omitnan");
    meanSSSMetric = mean(rows.MeanSSSMetric,"omitnan");
    meanCFOHz = mean(rows.MeanCFOHz,"omitnan");
    keep(k) = true;
    parts{k} = table( ...
        rows.PositionID(1),rows.CellKey(1),rows.GSCN(1),rows.PCI(1), ...
        height(rows),numRepeats,center,centerUncertaintyNs,rmsNs,ppNs, ...
        maxAbsNs,expectedNoiseNs,excessRmsNs, ...
        meanPSSMetric,meanSSSMetric,meanCFOHz, ...
        'VariableNames',["PositionID","CellKey","GSCN","PCI", ...
        "NumObservations","NumRepeats","RelativeArrivalCenterNs", ...
        "CenterUncertaintyNs","StationaryRmsNs", ...
        "StationaryPeakToPeakNs","StationaryMaxAbsNs", ...
        "ExpectedMeasurementNoiseRmsNs","ExcessStationaryRmsNs", ...
        "MeanPSSMetric","MeanSSSMetric","MeanCFOHz"]);
end
if ~any(keep)
    summary = emptyPerPositionTimingLocal();
    return;
end
summary = vertcat(parts{keep});
end

function stationary = summarizePerCellLocal(relative,~)
cells = unique(relative(:,["CellKey","GSCN","PCI"]),"rows","stable");
parts = cell(height(cells),1);
for k = 1:height(cells)
    rows = relative(relative.CellKey==cells.CellKey(k) & ...
        isfinite(relative.CenteredResidualNs),:);
    if isempty(rows)
        rmsNs = NaN; ppNs = NaN; maxAbsNs = NaN;
        expectedNoiseNs = NaN; excessRmsNs = NaN;
        status = "NOT_ASSESSABLE";
        reason = "Too few repeated stationary timing observations.";
    else
        residual = rows.CenteredResidualNs;
        rmsNs = sqrt(mean(residual.^2,"omitnan"));
        ppNs = max(residual)-min(residual);
        maxAbsNs = max(abs(residual));
        expectedNoiseNs = sqrt(mean(rows.SigmaNs.^2,"omitnan"));
        excessRmsNs = sqrt(max(0,rmsNs^2-expectedNoiseNs^2));
        status = "ASSESSABLE";
        reason = "Reference-subtracted timing variation measured while the receiver position was fixed.";
    end
    parts{k} = table( ...
        cells.CellKey(k),cells.GSCN(k),cells.PCI(k),height(rows), ...
        numel(unique(rows.PositionID)),numel(unique(rows.CaptureID)), ...
        rmsNs,ppNs,maxAbsNs,expectedNoiseNs,excessRmsNs,status,reason, ...
        'VariableNames',["CellKey","GSCN","PCI","NumObservations", ...
        "NumPositions","NumCaptures","StationaryTimingInstabilityRmsNs", ...
        "StationaryTimingInstabilityPeakToPeakNs", ...
        "StationaryTimingInstabilityMaxAbsNs", ...
        "ExpectedMeasurementNoiseRmsNs","ExcessStationaryRmsNs", ...
        "StationaryAssessmentStatus","StationaryAssessmentReason"]);
end
stationary = vertcat(parts{:});
end

function measurements = buildGeometryMeasurementsLocal( ...
        perPosition,referenceCellKey,referencePCI,referenceGSCN,opts)
positions = unique(perPosition.PositionID,"stable");
parts = {};
count = 0;
for p = 1:numel(positions)
    positionID = positions(p);
    rows = perPosition(perPosition.PositionID==positionID,:);
    if isempty(rows)
        continue;
    end
    relativeSigma = rows.CenterUncertaintyNs;
    refSigma = max(opts.TimingUncertaintyFloorNs, ...
        min(relativeSigma)/sqrt(2));
    count = count+1;
    parts{count,1} = geometryRowLocal( ...
        positionID,referenceCellKey,referenceGSCN,referencePCI,0, ...
        refSigma,NaN,NaN,NaN); %#ok<AGROW>
    for k = 1:height(rows)
        targetSigma = sqrt(max( ...
            rows.CenterUncertaintyNs(k)^2-refSigma^2, ...
            opts.TimingUncertaintyFloorNs^2));
        count = count+1;
        parts{count,1} = geometryRowLocal( ...
            positionID,rows.CellKey(k),rows.GSCN(k),rows.PCI(k), ...
            mod(rows.RelativeArrivalCenterNs(k),opts.FramePeriodNs), ...
            targetSigma,rows.MeanPSSMetric(k),rows.MeanSSSMetric(k), ...
            rows.MeanCFOHz(k));
    end
end
if count==0
    measurements = emptyGeometryMeasurementsLocal();
else
    measurements = vertcat(parts{:});
end
end

function row = geometryRowLocal(captureID,cellKey,gscn,pci,phaseNs, ...
        timingStdNs,pssMetric,sssMetric,cfoHz)
row = table(captureID,cellKey,gscn,pci,phaseNs,timingStdNs, ...
    NaN,NaN,NaN,pssMetric,sssMetric,cfoHz, ...
    'VariableNames',["CaptureID","CellKey","GSCN","PCI", ...
    "FramePhaseNs","TimingStdNs","CenterFrequencyHz", ...
    "SelectedSSBOffsetHz","HypothesizedSSBCenterHz", ...
    "MeanPSSMetric","MeanSSSMetric","MeanCFOHz"]);
end

function row = bestRowLocal(rows)
[~,idx] = max(rows.MeanSSSMetric);
row = rows(idx,:);
end

function [center,residual] = robustCircularCenterLocal(values,period)
values = values(:);
seed = mod(angle(sum(exp(1i*2*pi*values/period)))*period/(2*pi),period);
residual = wrapLocal(values-seed,period);
center = mod(seed+median(residual,"omitnan"),period);
residual = wrapLocal(values-center,period);
end

function y = wrapLocal(x,period)
y = mod(x+period/2,period)-period/2;
end

function [x,y] = localOffsetsLocal(lat,lon,lat0,lon0)
earthRadiusM = 6371000;
x = deg2rad(lon-lon0).*earthRadiusM.*cos(deg2rad(lat0));
y = deg2rad(lat-lat0).*earthRadiusM;
end

function result = emptyRelativeResultLocal(measurements,positionInfo, ...
        referenceCellKey,referencePCI,referenceGSCN,rejected,support,ranking)
cells = unique(measurements(:,["CellKey","GSCN","PCI"]),"rows","stable");
num = height(cells);
stationary = table( ...
    cells.CellKey,cells.GSCN,cells.PCI,zeros(num,1),zeros(num,1), ...
    zeros(num,1),NaN(num,1),NaN(num,1),NaN(num,1),NaN(num,1), ...
    NaN(num,1),repmat("NOT_ASSESSABLE",num,1), ...
    repmat("No analysis window contained both the reference cell and this cell.",num,1), ...
    'VariableNames',["CellKey","GSCN","PCI","NumObservations", ...
    "NumPositions","NumCaptures","StationaryTimingInstabilityRmsNs", ...
    "StationaryTimingInstabilityPeakToPeakNs", ...
    "StationaryTimingInstabilityMaxAbsNs", ...
    "ExpectedMeasurementNoiseRmsNs","ExcessStationaryRmsNs", ...
    "StationaryAssessmentStatus","StationaryAssessmentReason"]);

result = struct();
result.Mode = "stationary_repeats_aggregated_by_position";
result.Statement = ...
    "Stationary repeats were loaded, but no target-relative timing observations could be formed.";
result.ReferenceCellKey = referenceCellKey;
result.ReferencePCI = referencePCI;
result.ReferenceGSCN = referenceGSCN;
result.RawWindowMeasurements = [measurements; ...
    rejected(:,measurements.Properties.VariableNames)];
result.WindowMeasurements = measurements;
result.RejectedWindowMeasurements = rejected;
result.ReliableCellSupport = support;
result.ReferenceCandidateRanking = ranking;
result.RelativeObservations = emptyRelativeObservationsLocal();
result.PerPositionTiming = emptyPerPositionTimingLocal();
result.StationaryInstability = stationary;
result.AggregatedMeasurementTable = emptyGeometryMeasurementsLocal();
result.PositionInfo = positionInfo;
end

function result = emptyQualityResultLocal(raw,rejected,support,positionInfo)
result = struct();
result.Mode = "stationary_repeats_aggregated_by_position";
result.Statement = ...
    "No timing rows passed the repeated-detection and timing-quality gates.";
result.ReferenceCellKey = "";
result.ReferencePCI = NaN;
result.ReferenceGSCN = NaN;
result.RawWindowMeasurements = raw;
result.WindowMeasurements = raw([],:);
result.RejectedWindowMeasurements = rejected;
result.ReliableCellSupport = support;
result.ReferenceCandidateRanking = support;
result.RelativeObservations = emptyRelativeObservationsLocal();
result.PerPositionTiming = emptyPerPositionTimingLocal();
result.StationaryInstability = emptyStationaryInstabilityLocal();
result.AggregatedMeasurementTable = emptyGeometryMeasurementsLocal();
result.PositionInfo = positionInfo;
end

function tbl = emptyRelativeObservationsLocal()
tbl = table( ...
    strings(0,1),strings(0,1),zeros(0,1),strings(0,1), ...
    strings(0,1),zeros(0,1),zeros(0,1),strings(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',["PositionID","CaptureID","RepeatIndex", ...
    "ObservationID","CellKey","GSCN","PCI","ReferenceCellKey", ...
    "ReferencePCI","RelativeArrivalNs","SigmaNs","MeanPSSMetric", ...
    "MeanSSSMetric","MeanCFOHz","CenteredResidualNs"]);
end

function tbl = emptyPerPositionTimingLocal()
tbl = table( ...
    strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',["PositionID","CellKey","GSCN","PCI", ...
    "NumObservations","NumRepeats","RelativeArrivalCenterNs", ...
    "CenterUncertaintyNs","StationaryRmsNs", ...
    "StationaryPeakToPeakNs","StationaryMaxAbsNs", ...
    "ExpectedMeasurementNoiseRmsNs","ExcessStationaryRmsNs", ...
    "MeanPSSMetric","MeanSSSMetric","MeanCFOHz"]);
end

function tbl = emptyGeometryMeasurementsLocal()
tbl = table( ...
    strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1), ...
    'VariableNames',["CaptureID","CellKey","GSCN","PCI", ...
    "FramePhaseNs","TimingStdNs","CenterFrequencyHz", ...
    "SelectedSSBOffsetHz","HypothesizedSSBCenterHz", ...
    "MeanPSSMetric","MeanSSSMetric","MeanCFOHz"]);
end

function tbl = emptyReliableCellSupportLocal()
tbl = table( ...
    strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',["CellKey","GSCN","PCI","ReliablePositions", ...
    "ReliableCaptures","ReliableObservations","NumTimingRows", ...
    "MeanSSSMetric","MedianTimingStdNs"]);
end

function tbl = emptyStationaryInstabilityLocal()
tbl = table( ...
    strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),strings(0,1),strings(0,1), ...
    'VariableNames',["CellKey","GSCN","PCI","NumObservations", ...
    "NumPositions","NumCaptures","StationaryTimingInstabilityRmsNs", ...
    "StationaryTimingInstabilityPeakToPeakNs", ...
    "StationaryTimingInstabilityMaxAbsNs", ...
    "ExpectedMeasurementNoiseRmsNs","ExcessStationaryRmsNs", ...
    "StationaryAssessmentStatus","StationaryAssessmentReason"]);
end
