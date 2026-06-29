function result = noLocationSurveyCheck(measurementTable, captureInfo, opts)
%noLocationSurveyCheck Multi-position relative sync check without gNB locations.
%
%   result = noLocationSurveyCheck(measurementTable, captureInfo, opts)
%   estimates relative gNB transmit timing from multiple captures at known
%   receiver positions. gNB locations are not input. The method remains
%   relative: a common timing error shared by all gNBs is not visible with an
%   internal-clock B210.

arguments
    measurementTable table
    captureInfo table
    opts struct = struct()
end

[timingEstimates, fitInfo, relativeMeasurements] = estimateSurveyTimingNoLocations( ...
    measurementTable, captureInfo, opts);

result = struct();
result.Mode = "multi_position_unknown_gnb_locations";
result.Statement = "Relative survey mode. gNB locations are estimated, not provided. Absolute UTC phase compliance is not assessed with internal-clock B210 captures.";
result.VerdictBasis = "Synchronization verdict uses timing-only evidence: relative static offset, relative timing instability, and combined worst-case timing envelope.";
result.ExcludedFromVerdict = ["TDD pattern energy feedback","Special-slot energy feedback","RSS/RSSI planning feedback","CFO diagnostic reporting"];
result.Measurements = measurementTable;
result.CaptureInfo = captureInfo;
result.RelativeMeasurements = relativeMeasurements;
result.TimingEstimates = timingEstimates;
result.FitInfo = fitInfo;
result.DecisionTable = surveyDecisionTableLocal(timingEstimates);

end

function tbl = surveyDecisionTableLocal(timingEstimates)
if isempty(timingEstimates) || height(timingEstimates) == 0
    tbl = table();
    return;
end

tbl = timingEstimates(:, ["PCI","SurveyRole","ReferencePCI", ...
    "EstimatedRelativeTxOffsetNs","EstimatedOffsetUncertaintyNs", ...
    "StaticTimingStatus","StaticTimingReason", ...
    "StaticOffsetNs","TimingInstabilityRmsNs","TimingInstabilityPeakToPeakNs", ...
    "TimingInstabilityMaxAbsNs","WorstCaseRelativeTimingNs", ...
    "WorstCaseRelativeTimingUncertaintyNs","EnvelopeTimingStatus", ...
    "EnvelopeTimingReason", ...
    "ExcessTimingJitterRmsNs","TimingInstabilityStatus", ...
    "TimingInstabilityReason", ...
    "TimingStatus","TimingReason","FitRMSNs","ConditionNumber"]);
tbl.Properties.VariableNames("TimingStatus") = "Verdict";
tbl.Properties.VariableNames("TimingReason") = "Reason";
tbl.VerdictBasis = repmat("Timing-only: static offset + residual instability + combined envelope.", height(tbl), 1);
tbl.TDDPatternUsedForVerdict = false(height(tbl), 1);
tbl.SpecialSlotUsedForVerdict = false(height(tbl), 1);
end
