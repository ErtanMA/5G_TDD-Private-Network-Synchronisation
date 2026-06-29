function stability = jackknifeSurveyStability(measurementTable, captureInfo, surveyOpts, jkOpts)
%jackknifeSurveyStability Leave-one-position-out stability of survey offsets.
%
%   stability = jackknifeSurveyStability(measurementTable, captureInfo, ...
%       surveyOpts, jkOpts) re-fits the no-location survey once per receiver
%   position, each time dropping that position, and reports how much each
%   cell's estimated relative transmit offset moves.
%
%   This directly tests the offset-vs-geometry identifiability risk: the
%   linearized covariance (EstimatedOffsetUncertaintyNs) is a local, optimistic
%   measure under weak survey geometry. If a cell's offset stays put as
%   positions are removed, the geometry genuinely supports it. If it swings by
%   microseconds, that cell is geometry-limited regardless of what the
%   covariance says, and its PASS/FAIL offset verdict should not be trusted.
%
%   The reference PCI is held fixed across every fold and all offsets are
%   re-anchored to it, so the leave-one-out values are directly comparable
%   (the per-fold median-centering gauge is removed).
%
%   measurementTable/captureInfo are the same inputs as noLocationSurveyCheck.
%   A runRepeatedPositionSurvey or noLocationSurveyCheck result struct can be
%   passed as the first argument instead, in which case the aggregated geometry
%   measurements, position info, and reference PCI are extracted from it.
%
%   surveyOpts is forwarded to noLocationSurveyCheck (e.g. SyncThresholdNs).
%   jkOpts controls the stability classification:
%       SyncThresholdNs        - scale for STABLE/MODERATE/UNSTABLE bands
%                                (default: taken from the survey fit, else 1000)
%       StableFraction         - range <= StableFraction*threshold -> STABLE
%                                (default 0.10)
%       ModerateFraction       - range <= ModerateFraction*threshold -> MODERATE
%                                (default 0.50)
%       MinFoldsForAssessment  - minimum assessable folds to classify (default 3)

if nargin < 3 || isempty(surveyOpts)
    surveyOpts = struct();
end
if nargin < 4 || isempty(jkOpts)
    jkOpts = struct();
end

if isstruct(measurementTable)
    [measurementTable, captureInfo, surveyOpts] = ...
        extractFromResultLocal(measurementTable, surveyOpts);
end

if ~istable(measurementTable) || ~istable(captureInfo)
    error("jackknifeSurveyStability:InvalidInput", ...
        "measurementTable and captureInfo must be tables (or pass a survey result struct).");
end

jk = defaultsLocal(jkOpts);

allPositionIDs = unique(string(captureInfo.CaptureID),"stable");
measurementPositionIDs = unique(string(measurementTable.CaptureID),"stable");
if jk.OnlyMeasuredPositions
    usedPositionIDs = allPositionIDs(ismember(allPositionIDs,measurementPositionIDs));
else
    usedPositionIDs = allPositionIDs;
end
unusedPositionIDs = allPositionIDs(~ismember(allPositionIDs,usedPositionIDs));

% Baseline full fit. Fix the chosen reference PCI for every fold so the
% gauge cell does not change between folds.
baseResult = noLocationSurveyCheck(measurementTable, captureInfo, surveyOpts);
refPCI = baseResult.FitInfo.ReferencePCI;
surveyOpts.ReferencePCI = refPCI;

threshold = thresholdLocal(jk, baseResult);
baseAnchored = reAnchoredOffsetsLocal(baseResult.TimingEstimates, refPCI);

captureIDs = usedPositionIDs;
numFolds = numel(captureIDs);

measIDs = string(measurementTable.CaptureID);
capIDs = string(captureInfo.CaptureID);

foldParts = cell(numFolds, 1);
for f = 1:numFolds
    dropped = captureIDs(f);
    capSub = captureInfo(capIDs ~= dropped, :);
    measSub = measurementTable(measIDs ~= dropped, :);
    foldResult = noLocationSurveyCheck(measSub, capSub, surveyOpts);
    foldParts{f} = foldRowsLocal(foldResult, refPCI, dropped);
end
folds = vertcat(foldParts{:});

perCell = aggregatePerCellLocal(folds, baseAnchored, threshold, jk, refPCI);

stability = struct();
stability.ReferencePCI = refPCI;
stability.SyncThresholdNs = threshold;
stability.NumFolds = numFolds;
stability.NumPositions = numFolds;
stability.NumInputPositions = numel(allPositionIDs);
stability.UsedPositionIDs = usedPositionIDs;
stability.UnusedPositionIDs = unusedPositionIDs;
stability.Baseline = baseAnchored;
stability.BaselineFitInfo = baseResult.FitInfo;
stability.Folds = folds;
stability.PerCell = perCell;
stability.JackknifeOptions = jk;
stability.Statement = ...
    "Leave-one-contributing-position-out jackknife: each fold drops one receiver position that supplied an accepted geometry measurement and re-fits the no-location survey. Wide swings indicate sensitivity to the retained position set; they are evidence of weak geometry, not an independent proof of transmitter timing.";

end

function jk = defaultsLocal(jkOpts)
jk = struct();
jk.SyncThresholdNs = NaN;
jk.StableFraction = 0.10;
jk.ModerateFraction = 0.50;
jk.MinFoldsForAssessment = 3;
jk.OnlyMeasuredPositions = true;
names = fieldnames(jkOpts);
for k = 1:numel(names)
    jk.(names{k}) = jkOpts.(names{k});
end
end

function threshold = thresholdLocal(jk, baseResult)
if isfield(jk, "SyncThresholdNs") && isfinite(jk.SyncThresholdNs)
    threshold = double(jk.SyncThresholdNs);
    return;
end
threshold = 1000;
if isfield(baseResult, "FitInfo") && isfield(baseResult.FitInfo, "SyncThresholdNs") && ...
        isfinite(baseResult.FitInfo.SyncThresholdNs)
    threshold = double(baseResult.FitInfo.SyncThresholdNs);
end
end

function tbl = reAnchoredOffsetsLocal(timing, refPCI)
if isempty(timing) || height(timing) == 0
    tbl = table(zeros(0,1), zeros(0,1), zeros(0,1), strings(0,1), strings(0,1), ...
        'VariableNames', ["PCI","ReAnchoredOffsetNs","OffsetUncertaintyNs", ...
        "StaticTimingStatus","SurveyRole"]);
    return;
end
pci = double(timing.PCI);
off = double(timing.EstimatedRelativeTxOffsetNs);
unc = double(timing.EstimatedOffsetUncertaintyNs);
status = string(timing.StaticTimingStatus);
if ismember("SurveyRole", string(timing.Properties.VariableNames))
    role = string(timing.SurveyRole);
else
    role = repmat("TARGET", height(timing), 1);
end

anchorIdx = find(pci == refPCI, 1);
if isempty(anchorIdx) || ~isfinite(off(anchorIdx))
    anchor = 0;
else
    anchor = off(anchorIdx);
end
reanch = off - anchor;

tbl = table(pci, reanch, unc, status, role, ...
    'VariableNames', ["PCI","ReAnchoredOffsetNs","OffsetUncertaintyNs", ...
    "StaticTimingStatus","SurveyRole"]);
end

function rows = foldRowsLocal(foldResult, refPCI, droppedID)
timing = foldResult.TimingEstimates;
anch = reAnchoredOffsetsLocal(timing, refPCI);
n = height(anch);
if n == 0
    rows = emptyFoldTableLocal();
    return;
end

fitRMS = NaN;
condNum = NaN;
fitStatus = "UNKNOWN";
if isfield(foldResult, "FitInfo")
    if isfield(foldResult.FitInfo, "FitRMSNs"); fitRMS = foldResult.FitInfo.FitRMSNs; end
    if isfield(foldResult.FitInfo, "ConditionNumber"); condNum = foldResult.FitInfo.ConditionNumber; end
    if isfield(foldResult.FitInfo, "Status"); fitStatus = string(foldResult.FitInfo.Status); end
end

rows = table( ...
    repmat(string(droppedID), n, 1), anch.PCI, anch.ReAnchoredOffsetNs, ...
    anch.OffsetUncertaintyNs, anch.StaticTimingStatus, ...
    repmat(fitRMS, n, 1), repmat(condNum, n, 1), repmat(fitStatus, n, 1), ...
    'VariableNames', ["DroppedCaptureID","PCI","ReAnchoredOffsetNs", ...
    "OffsetUncertaintyNs","StaticTimingStatus","FoldFitRMSNs", ...
    "FoldConditionNumber","FoldFitStatus"]);
end

function tbl = emptyFoldTableLocal()
tbl = table(strings(0,1), zeros(0,1), zeros(0,1), zeros(0,1), strings(0,1), ...
    zeros(0,1), zeros(0,1), strings(0,1), ...
    'VariableNames', ["DroppedCaptureID","PCI","ReAnchoredOffsetNs", ...
    "OffsetUncertaintyNs","StaticTimingStatus","FoldFitRMSNs", ...
    "FoldConditionNumber","FoldFitStatus"]);
end

function perCell = aggregatePerCellLocal(folds, baseAnchored, threshold, jk, refPCI)
pcis = baseAnchored.PCI;
if ~isempty(folds) && height(folds) > 0
    pcis = unique([pcis; folds.PCI], "stable");
end

parts = cell(numel(pcis), 1);
for k = 1:numel(pcis)
    pci = pcis(k);

    baseIdx = find(baseAnchored.PCI == pci, 1);
    if isempty(baseIdx)
        baseOffset = NaN;
        baseUnc = NaN;
        role = "TARGET";
    else
        baseOffset = baseAnchored.ReAnchoredOffsetNs(baseIdx);
        baseUnc = baseAnchored.OffsetUncertaintyNs(baseIdx);
        role = baseAnchored.SurveyRole(baseIdx);
    end

    if pci == refPCI
        role = "GAUGE_REFERENCE";
    end

    foldVals = zeros(0,1);
    foldDropped = strings(0,1);
    if ~isempty(folds) && height(folds) > 0
        sel = folds.PCI == pci & isfinite(folds.ReAnchoredOffsetNs) & ...
            folds.StaticTimingStatus ~= "NOT_ASSESSABLE";
        foldVals = folds.ReAnchoredOffsetNs(sel);
        foldDropped = folds.DroppedCaptureID(sel);
    end
    numAssessable = numel(foldVals);

    if numAssessable >= 1
        jkMedian = median(foldVals);
        jkMin = min(foldVals);
        jkMax = max(foldVals);
        jkRange = jkMax - jkMin;
        if numAssessable >= 2
            jkStd = std(foldVals);
        else
            jkStd = 0;
        end
        if isfinite(baseOffset)
            devFromBase = abs(foldVals - baseOffset);
        else
            devFromBase = abs(foldVals - jkMedian);
        end
        [maxDev, worstIdx] = max(devFromBase);
        mostInfluential = foldDropped(worstIdx);
    else
        jkMedian = NaN; jkMin = NaN; jkMax = NaN; jkRange = NaN;
        jkStd = NaN; maxDev = NaN; mostInfluential = "";
    end

    [status, reason] = classifyLocal(pci, refPCI, numAssessable, jkRange, ...
        maxDev, threshold, jk);

    parts{k} = table(pci, string(role), baseOffset, baseUnc, ...
        jkMedian, jkStd, jkRange, maxDev, numAssessable, ...
        string(mostInfluential), string(status), string(reason), ...
        'VariableNames', ["PCI","SurveyRole","BaselineOffsetNs", ...
        "LinearUncertaintyNs","JackknifeMedianNs","JackknifeStdNs", ...
        "JackknifeRangeNs","MaxAbsDevFromBaselineNs","NumFoldsAssessable", ...
        "MostInfluentialDroppedID","StabilityStatus","StabilityReason"]);
end
perCell = vertcat(parts{:});
end

function [status, reason] = classifyLocal(pci, refPCI, numAssessable, jkRange, maxDev, threshold, jk)
if pci == refPCI
    status = "REFERENCE";
    reason = "Gauge reference cell; offset is fixed to zero by construction and is not jackknifed.";
    return;
end
if numAssessable < jk.MinFoldsForAssessment || ~isfinite(jkRange)
    status = "NOT_ASSESSABLE";
    reason = sprintf("Only %d leave-one-out fold(s) produced an assessable offset; need %d to judge stability.", ...
        numAssessable, jk.MinFoldsForAssessment);
    return;
end

stableLimit = jk.StableFraction * threshold;
moderateLimit = jk.ModerateFraction * threshold;
if jkRange <= stableLimit
    status = "STABLE";
    reason = sprintf("Offset moves at most %.0f ns across leave-one-out fits (<= %.0f ns); geometry supports this offset.", ...
        jkRange, stableLimit);
elseif jkRange <= moderateLimit
    status = "MODERATE";
    reason = sprintf("Offset moves %.0f ns across leave-one-out fits (between %.0f and %.0f ns); treat the offset as indicative, not strict.", ...
        jkRange, stableLimit, moderateLimit);
else
    status = "UNSTABLE";
    reason = sprintf("Offset swings %.0f ns across leave-one-out fits (> %.0f ns); not separable from gNB geometry with this survey. Worst single position changes it by %.0f ns.", ...
        jkRange, moderateLimit, maxDev);
end
end

function [measurementTable, captureInfo, surveyOpts] = extractFromResultLocal(res, surveyOpts)
if isfield(res, "AggregatedMeasurementTable") && isfield(res, "PositionInfo")
    measurementTable = res.AggregatedMeasurementTable;
    captureInfo = res.PositionInfo;
elseif isfield(res, "Measurements") && isfield(res, "CaptureInfo")
    measurementTable = res.Measurements;
    captureInfo = res.CaptureInfo;
elseif isfield(res, "SurveyResult")
    [measurementTable, captureInfo, surveyOpts] = ...
        extractFromResultLocal(res.SurveyResult, surveyOpts);
    if isfield(res, "AggregatedMeasurementTable") && isfield(res, "PositionInfo")
        measurementTable = res.AggregatedMeasurementTable;
        captureInfo = res.PositionInfo;
    end
else
    error("jackknifeSurveyStability:UnsupportedResult", ...
        "Result struct must contain AggregatedMeasurementTable/PositionInfo or Measurements/CaptureInfo.");
end

if ~isfield(surveyOpts, "ReferencePCI") || ~isfinite(surveyOpts.ReferencePCI)
    if isfield(res, "ReferencePCI") && isfinite(res.ReferencePCI)
        surveyOpts.ReferencePCI = double(res.ReferencePCI);
    end
end
end
