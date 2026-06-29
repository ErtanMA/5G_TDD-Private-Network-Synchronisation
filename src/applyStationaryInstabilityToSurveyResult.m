function result = applyStationaryInstabilityToSurveyResult(result,stationary,opts)
%applyStationaryInstabilityToSurveyResult Use stationary instability in verdict.
%
% Survey residuals remain available as fit-consistency metrics. The final
% worst-case envelope is rebuilt from fitted static offset plus the maximum
% stationary reference-subtracted timing deviation.

arguments
    result struct
    stationary table
    opts struct = struct()
end

if ~isfield(result,"TimingEstimates") || isempty(result.TimingEstimates)
    result.StationaryInstability = stationary;
    return;
end

timing = result.TimingEstimates;
num = height(timing);
timing.SurveyModelResidualSamples = timing.TimingInstabilitySamples;
timing.SurveyModelResidualRmsNs = timing.TimingInstabilityRmsNs;
timing.SurveyModelResidualPeakToPeakNs = ...
    timing.TimingInstabilityPeakToPeakNs;
timing.SurveyModelResidualMaxAbsNs = timing.TimingInstabilityMaxAbsNs;
timing.SurveyModelExpectedTimingNoiseRmsNs = ...
    timing.ExpectedTimingNoiseRmsNs;
timing.StationaryTimingInstabilityRmsNs = NaN(num,1);
timing.StationaryTimingInstabilityPeakToPeakNs = NaN(num,1);
timing.StationaryTimingInstabilityMaxAbsNs = NaN(num,1);
timing.StationaryExpectedMeasurementNoiseRmsNs = NaN(num,1);
timing.StationaryExcessTimingInstabilityRmsNs = NaN(num,1);
timing.StationaryNumObservations = zeros(num,1);
timing.StationaryInstabilityStatus = repmat("NOT_ASSESSABLE",num,1);
timing.StationaryInstabilityReason = repmat( ...
    "No stationary repeated-capture evidence was available.",num,1);

syncThresholdNs = optionLocal(opts,"SyncThresholdNs", ...
    fieldLocal(result.FitInfo,"SyncThresholdNs",1000));
hardFailNs = optionLocal(opts,"HardFailOffsetNs", ...
    fieldLocal(result.FitInfo,"HardFailOffsetNs",Inf));
sigmaMultiplier = optionLocal(opts,"UncertaintySigmaMultiplier", ...
    fieldLocal(result.FitInfo,"UncertaintySigmaMultiplier",3));
rmsThresholdNs = firstFiniteLocal(opts, ...
    ["JitterRmsThresholdNs","InstabilityRmsThresholdNs"],NaN);
ppThresholdNs = firstFiniteLocal(opts, ...
    ["JitterPeakToPeakThresholdNs","InstabilityPeakToPeakThresholdNs"],NaN);

for k = 1:num
    role = string(timing.SurveyRole(k));
    if isReferenceRoleLocal(role)
        timing.StationaryTimingInstabilityRmsNs(k) = 0;
        timing.StationaryTimingInstabilityPeakToPeakNs(k) = 0;
        timing.StationaryTimingInstabilityMaxAbsNs(k) = 0;
        timing.StationaryExpectedMeasurementNoiseRmsNs(k) = 0;
        timing.StationaryExcessTimingInstabilityRmsNs(k) = 0;
        timing.StationaryInstabilityStatus(k) = "REFERENCE";
        timing.StationaryInstabilityReason(k) = ...
            "Gauge reference; target-relative instability is reported on other cells.";
        continue;
    end

    row = stationary(stationary.PCI==timing.PCI(k),:);
    if isempty(row)
        continue;
    end
    row = row(1,:);
    timing.StationaryTimingInstabilityRmsNs(k) = ...
        row.StationaryTimingInstabilityRmsNs;
    timing.StationaryTimingInstabilityPeakToPeakNs(k) = ...
        row.StationaryTimingInstabilityPeakToPeakNs;
    timing.StationaryTimingInstabilityMaxAbsNs(k) = ...
        row.StationaryTimingInstabilityMaxAbsNs;
    timing.StationaryExpectedMeasurementNoiseRmsNs(k) = ...
        row.ExpectedMeasurementNoiseRmsNs;
    timing.StationaryExcessTimingInstabilityRmsNs(k) = ...
        row.ExcessStationaryRmsNs;
    timing.StationaryNumObservations(k) = row.NumObservations;
    [status,reason] = stationaryStatusLocal(row,rmsThresholdNs, ...
        ppThresholdNs,sigmaMultiplier);
    timing.StationaryInstabilityStatus(k) = status;
    timing.StationaryInstabilityReason(k) = reason;
end

timing.TimingInstabilityRmsNs = ...
    timing.StationaryTimingInstabilityRmsNs;
timing.TimingInstabilityPeakToPeakNs = ...
    timing.StationaryTimingInstabilityPeakToPeakNs;
timing.TimingInstabilityMaxAbsNs = ...
    timing.StationaryTimingInstabilityMaxAbsNs;
timing.TimingInstabilitySamples = timing.StationaryNumObservations;
timing.ExpectedTimingNoiseRmsNs = ...
    timing.StationaryExpectedMeasurementNoiseRmsNs;
timing.ExcessTimingJitterRmsNs = ...
    timing.StationaryExcessTimingInstabilityRmsNs;
timing.TimingInstabilityStatus = timing.StationaryInstabilityStatus;
timing.TimingInstabilityReason = timing.StationaryInstabilityReason;

for k = 1:num
    if isReferenceRoleLocal(timing.SurveyRole(k))
        timing.WorstCaseRelativeTimingNs(k) = abs(timing.StaticOffsetNs(k));
        timing.WorstCaseRelativeTimingUncertaintyNs(k) = ...
            timing.EstimatedOffsetUncertaintyNs(k);
        timing.EnvelopeTimingStatus(k) = timing.StaticTimingStatus(k);
        timing.EnvelopeTimingReason(k) = ...
            "Gauge reference; no target-relative stationary envelope is assigned.";
        continue;
    end

    if string(timing.StaticTimingStatus(k))=="NOT_ASSESSABLE"
        timing.WorstCaseRelativeTimingNs(k) = NaN;
        timing.WorstCaseRelativeTimingUncertaintyNs(k) = NaN;
        timing.EnvelopeTimingStatus(k) = "NOT_ASSESSABLE";
        timing.EnvelopeTimingReason(k) = string(timing.StaticTimingReason(k));
        timing.TimingStatus(k) = "NOT_ASSESSABLE";
        timing.TimingReason(k) = string(timing.StaticTimingReason(k));
        continue;
    end

    maxAbsNs = timing.StationaryTimingInstabilityMaxAbsNs(k);
    if ~isfinite(maxAbsNs)
        timing.WorstCaseRelativeTimingNs(k) = NaN;
        timing.WorstCaseRelativeTimingUncertaintyNs(k) = NaN;
        timing.EnvelopeTimingStatus(k) = "NOT_ASSESSABLE";
        timing.EnvelopeTimingReason(k) = ...
            "Stationary timing instability was not assessable.";
    else
        timing.WorstCaseRelativeTimingNs(k) = ...
            abs(timing.StaticOffsetNs(k))+maxAbsNs;
        timing.WorstCaseRelativeTimingUncertaintyNs(k) = hypot( ...
            timing.EstimatedOffsetUncertaintyNs(k), ...
            timing.StationaryExpectedMeasurementNoiseRmsNs(k));
        [timing.EnvelopeTimingStatus(k),timing.EnvelopeTimingReason(k)] = ...
            envelopeStatusLocal(timing(k,:),result.FitInfo, ...
            syncThresholdNs,hardFailNs,sigmaMultiplier);
    end
    [timing.TimingStatus(k),timing.TimingReason(k)] = ...
        finalStatusLocal(timing(k,:),result.FitInfo);
end

result.TimingEstimates = timing;
result.StationaryInstability = stationary;
result.VerdictBasis = ...
    "Timing-only: fitted static relative offset plus stationary repeated-capture timing instability. Survey-model residuals are fit-quality evidence, not transmitter jitter.";
result.DecisionTable = decisionTableLocal(timing);
result.FitInfo.MaxSurveyModelResidualRmsNs = max( ...
    timing.SurveyModelResidualRmsNs,[],"omitnan");
result.FitInfo.MaxSurveyModelResidualPeakToPeakNs = max( ...
    timing.SurveyModelResidualPeakToPeakNs,[],"omitnan");
result.FitInfo.MaxSurveyModelResidualMaxAbsNs = max( ...
    timing.SurveyModelResidualMaxAbsNs,[],"omitnan");
result.FitInfo.MaxStationaryTimingInstabilityRmsNs = max( ...
    timing.StationaryTimingInstabilityRmsNs,[],"omitnan");
result.FitInfo.MaxStationaryTimingInstabilityPeakToPeakNs = max( ...
    timing.StationaryTimingInstabilityPeakToPeakNs,[],"omitnan");
result.FitInfo.MaxStationaryTimingInstabilityMaxAbsNs = max( ...
    timing.StationaryTimingInstabilityMaxAbsNs,[],"omitnan");
result.FitInfo.MaxTimingInstabilityRmsNs = ...
    result.FitInfo.MaxStationaryTimingInstabilityRmsNs;
result.FitInfo.MaxTimingInstabilityPeakToPeakNs = ...
    result.FitInfo.MaxStationaryTimingInstabilityPeakToPeakNs;
result.FitInfo.MaxTimingInstabilityMaxAbsNs = ...
    result.FitInfo.MaxStationaryTimingInstabilityMaxAbsNs;
result.FitInfo.MaxExcessTimingJitterRmsNs = max( ...
    timing.StationaryExcessTimingInstabilityRmsNs,[],"omitnan");
result.FitInfo.MaxWorstCaseRelativeTimingNs = max( ...
    timing.WorstCaseRelativeTimingNs,[],"omitnan");
end

function [status,reason] = stationaryStatusLocal(row,rmsLimit,ppLimit,kappa)
if string(row.StationaryAssessmentStatus)~="ASSESSABLE"
    status = "NOT_ASSESSABLE";
    reason = string(row.StationaryAssessmentReason);
    return;
end
if ~isfinite(rmsLimit) && ~isfinite(ppLimit)
    status = "PASS";
    reason = ...
        "Stationary timing instability is assessable and is used in the combined timing envelope; no separate instability threshold was configured.";
    return;
end

num = max(1,row.NumObservations);
rmsMargin = kappa*row.ExpectedMeasurementNoiseRmsNs/sqrt(num);
ppMargin = kappa*row.ExpectedMeasurementNoiseRmsNs*sqrt(2);
fail = (isfinite(rmsLimit) && ...
    row.StationaryTimingInstabilityRmsNs-rmsMargin>=rmsLimit) || ...
    (isfinite(ppLimit) && ...
    row.StationaryTimingInstabilityPeakToPeakNs-ppMargin>=ppLimit);
near = (isfinite(rmsLimit) && ...
    row.StationaryTimingInstabilityRmsNs+rmsMargin>=rmsLimit) || ...
    (isfinite(ppLimit) && ...
    row.StationaryTimingInstabilityPeakToPeakNs+ppMargin>=ppLimit);
if fail
    status = "FAIL";
    reason = "Stationary relative timing instability exceeds its configured threshold with uncertainty margin.";
elseif near
    status = "SUSPECT";
    reason = "Stationary relative timing instability overlaps its configured threshold after uncertainty handling.";
else
    status = "PASS";
    reason = "Stationary relative timing instability is below its configured threshold.";
end
end

function [status,reason] = envelopeStatusLocal(row,fitInfo,threshold,hardFail,kappa)
envelope = row.WorstCaseRelativeTimingNs;
margin = kappa*row.WorstCaseRelativeTimingUncertaintyNs;
if isfinite(hardFail) && envelope>=hardFail
    status = "FAIL";
    reason = "Static offset plus stationary maximum timing deviation exceeds the hard-fail threshold.";
elseif envelope-margin>=threshold
    status = "FAIL";
    reason = "Static offset plus stationary maximum timing deviation exceeds the synchronization threshold with uncertainty margin.";
elseif envelope+margin>=threshold
    status = "SUSPECT";
    reason = "The stationary worst-case timing envelope overlaps the synchronization threshold.";
elseif string(fitInfo.Status)~="OK"
    status = "SUSPECT";
    reason = "The stationary timing envelope is below threshold, but survey-model fit quality is not clean.";
else
    status = "PASS";
    reason = "The stationary worst-case timing envelope is below the synchronization threshold.";
end
end

function [status,reason] = finalStatusLocal(row,fitInfo)
if string(row.StaticTimingStatus)=="NOT_ASSESSABLE"
    status = "NOT_ASSESSABLE";
    reason = string(row.StaticTimingReason);
    return;
end
if string(fitInfo.Status)=="NOT_ASSESSABLE"
    status = "NOT_ASSESSABLE";
    reason = string(fitInfo.Reason);
    return;
end
statuses = [string(row.StaticTimingStatus), ...
    string(row.StationaryInstabilityStatus), ...
    string(row.EnvelopeTimingStatus)];
if any(statuses=="FAIL")
    status = "FAIL";
    reason = "At least one static, stationary-instability, or combined-envelope check failed.";
elseif any(statuses=="NOT_ASSESSABLE")
    status = "SUSPECT";
    reason = "The survey fit exists, but stationary timing instability was not fully assessable.";
elseif any(statuses=="SUSPECT") || string(fitInfo.Status)~="OK"
    status = "SUSPECT";
    reason = "At least one timing component or survey-fit check is uncertain.";
else
    status = "PASS";
    reason = "Static offset and stationary worst-case timing envelope are below threshold with acceptable uncertainty.";
end
end

function tbl = decisionTableLocal(timing)
names = ["PCI","SurveyRole","ReferencePCI","StaticOffsetNs", ...
    "EstimatedOffsetUncertaintyNs","StaticTimingStatus", ...
    "StationaryTimingInstabilityRmsNs", ...
    "StationaryTimingInstabilityPeakToPeakNs", ...
    "StationaryTimingInstabilityMaxAbsNs", ...
    "StationaryInstabilityStatus","SurveyModelResidualRmsNs", ...
    "SurveyModelResidualMaxAbsNs", ...
    "WorstCaseRelativeTimingNs","WorstCaseRelativeTimingUncertaintyNs", ...
    "EnvelopeTimingStatus","TimingStatus","TimingReason"];
tbl = timing(:,names);
tbl.Properties.VariableNames("TimingStatus") = "Verdict";
tbl.Properties.VariableNames("TimingReason") = "Reason";
end

function value = optionLocal(opts,name,defaultValue)
if isfield(opts,name) && isnumeric(opts.(name)) && ...
        isscalar(opts.(name)) && isfinite(opts.(name))
    value = double(opts.(name));
else
    value = defaultValue;
end
end

function value = firstFiniteLocal(opts,names,defaultValue)
value = defaultValue;
for k = 1:numel(names)
    name = char(names(k));
    if isfield(opts,name) && isnumeric(opts.(name)) && ...
            isscalar(opts.(name)) && isfinite(opts.(name))
        value = double(opts.(name));
        return;
    end
end
end

function value = fieldLocal(s,name,defaultValue)
if isfield(s,name)
    value = s.(name);
else
    value = defaultValue;
end
end

function tf = isReferenceRoleLocal(role)
role = string(role);
tf = role=="GAUGE_REFERENCE" || role=="REFERENCE";
end
