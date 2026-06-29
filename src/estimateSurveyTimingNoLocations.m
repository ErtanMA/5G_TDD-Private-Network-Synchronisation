function [timingEstimates, fitInfo, relativeMeasurements] = estimateSurveyTimingNoLocations(measurementTable, captureInfo, opts)
%estimateSurveyTimingNoLocations Estimate relative timing without gNB locations.
%
%   This is the no-private-gNB-location phase. It takes per-capture timing
%   measurements from the existing receiver chain, subtracts a reference cell
%   inside each capture to remove the B210's arbitrary capture clock phase,
%   and fits unknown gNB positions plus relative transmit timing offsets.

arguments
    measurementTable table
    captureInfo table
    opts struct = struct()
end

userOpts = opts;
opts = mergeStructsLocal(defaultSurveyOptions(), opts);
opts = normalizeDecisionOptionsLocal(opts, userOpts);
if height(measurementTable) == 0
    relativeMeasurements = emptyRelativeMeasurementsLocal();
    surveyInfo = emptySurveyInfoLocal(captureInfo, "No detected PCIs.");
    timingEstimates = emptyTimingEstimateTableLocal();
    fitInfo = baseFitInfoLocal(opts, surveyInfo, "NOT_ASSESSABLE", ...
        "No detected PCIs.");
    return;
end

[relativeMeasurements, surveyInfo] = prepareSurveyRelativeMeasurements(measurementTable, captureInfo, opts);

visiblePCIs = unique(double(measurementTable.PCI), "stable");
visiblePCIs = orderedPCIsLocal(visiblePCIs, surveyInfo.ReferencePCI);

if isfield(surveyInfo, "InputValidationStatus") && surveyInfo.InputValidationStatus == "NOT_ASSESSABLE"
    timingEstimates = notAssessableRowsLocal(visiblePCIs, surveyInfo.ReferencePCI, ...
        surveyInfo.InputValidationReason);
    fitInfo = baseFitInfoLocal(opts, surveyInfo, "NOT_ASSESSABLE", ...
        surveyInfo.InputValidationReason);
    fitInfo.NumGNBs = numel(visiblePCIs);
    return;
end

if opts.RequireReferenceInEveryCapture && ~isempty(surveyInfo.MissingReferenceCaptureIDs)
    reason = "Reference PCI is missing in at least one capture.";
    timingEstimates = notAssessableRowsLocal(visiblePCIs, surveyInfo.ReferencePCI, reason);
    fitInfo = baseFitInfoLocal(opts, surveyInfo, "NOT_ASSESSABLE", reason);
    fitInfo.NumGNBs = numel(visiblePCIs);
    fitInfo.NumCapturesUsed = surveyInfo.NumCapturesWithReference;
    return;
end

if isempty(relativeMeasurements) || height(relativeMeasurements) == 0
    if numel(visiblePCIs) < opts.MinVisibleGNBs
        reason = "Fewer than two gNBs are visible.";
    else
        reason = "No relative measurements could be formed.";
    end
    timingEstimates = notAssessableRowsLocal(visiblePCIs, surveyInfo.ReferencePCI, reason);
    fitInfo = baseFitInfoLocal(opts, surveyInfo, "NOT_ASSESSABLE", ...
        reason);
    fitInfo.NumGNBs = numel(visiblePCIs);
    return;
end

referencePCI = surveyInfo.ReferencePCI;
[fitMeasurements, targetSupport, excludedPCIs] = ...
    selectEstimableTargetsLocal(relativeMeasurements, opts);

if isempty(fitMeasurements) || height(fitMeasurements) == 0
    reason = "No target PCI has enough distinct receiver positions for the no-location timing model.";
    timingEstimates = notAssessableRowsLocal(visiblePCIs, referencePCI, reason);
    fitInfo = baseFitInfoLocal(opts, surveyInfo, "NOT_ASSESSABLE", reason);
    fitInfo.NumVisibleGNBs = numel(visiblePCIs);
    fitInfo.TargetSupport = targetSupport;
    fitInfo.ExcludedPCIs = excludedPCIs;
    return;
end

pcisTarget = unique(fitMeasurements.PCI, "stable");
pcisAll = [referencePCI; pcisTarget(:)];
numG = numel(pcisAll);
numTargets = numel(pcisTarget);
numCaptures = numel(unique(fitMeasurements.CaptureID));
numMeasurements = height(fitMeasurements);
numUnknowns = 2*numG + numTargets;
degreesOfFreedom = numMeasurements - numUnknowns;

fitInfo = baseFitInfoLocal(opts, surveyInfo, "OK", "");
fitInfo.NumGNBs = numG;
fitInfo.NumVisibleGNBs = numel(visiblePCIs);
fitInfo.NumTargetGNBs = numTargets;
fitInfo.NumCapturesUsed = numCaptures;
fitInfo.NumMeasurements = numMeasurements;
fitInfo.NumUnknowns = numUnknowns;
fitInfo.DegreesOfFreedom = degreesOfFreedom;
fitInfo.TargetSupport = targetSupport;
fitInfo.ExcludedPCIs = excludedPCIs;

if numG < opts.MinVisibleGNBs
    timingEstimates = notAssessableRowsLocal(pcisAll, referencePCI, "Fewer than two gNBs are visible.");
    fitInfo.Status = "NOT_ASSESSABLE";
    fitInfo.Reason = "Fewer than two gNBs are visible.";
    return;
end

if opts.MinCaptures > 0 && numCaptures < opts.MinCaptures
    timingEstimates = notAssessableRowsLocal(pcisAll, referencePCI, "Too few receiver positions for no-location geometry fitting.");
    fitInfo.Status = "NOT_ASSESSABLE";
    fitInfo.Reason = "Too few receiver positions for no-location geometry fitting.";
    return;
end

if numMeasurements < numUnknowns
    timingEstimates = notAssessableRowsLocal(visiblePCIs, referencePCI, "Survey is underdetermined; more shared receiver positions are needed.");
    fitInfo.Status = "NOT_ASSESSABLE";
    fitInfo.Reason = "Survey is underdetermined; more shared receiver positions are needed.";
    return;
end

problem = buildProblemLocal(fitMeasurements, pcisAll, pcisTarget, opts);
[bestX, bestResidual, bestCost, optimizerName, startDiagnostics] = ...
    solveProblemLocal(problem, opts);

[biasStdNs, positionStdM, conditionNumber, numericalDiagnostics] = ...
    estimateUncertaintyLocal(bestX, problem, bestResidual);
fitRMSNs = sqrt(mean(bestResidual.^2));
weightedRMS = sqrt(mean((bestResidual ./ problem.sigmaNs).^2));
timingInstability = timingInstabilityFromResidualsLocal(bestResidual, problem, opts);

fitInfo.Optimizer = optimizerName;
fitInfo.Cost = bestCost;
fitInfo.FitRMSNs = fitRMSNs;
fitInfo.WeightedRMS = weightedRMS;
fitInfo.ConditionNumber = conditionNumber;
fitInfo.ConditionNumberBasis = "scaled weighted Jacobian J";
fitInfo.NumericalDiagnostics = numericalDiagnostics;
fitInfo.StartDiagnostics = startDiagnostics;
fitInfo.MaxNearOptimalOffsetSpreadNs = ...
    startDiagnostics.MaxNearOptimalOffsetSpreadNs;
fitInfo.MaxTimingInstabilityRmsNs = max(timingInstability.TimingInstabilityRmsNs, [], "omitnan");
fitInfo.MaxExcessTimingJitterRmsNs = max(timingInstability.ExcessTimingJitterRmsNs, [], "omitnan");
fitInfo.MaxTimingInstabilityPeakToPeakNs = max(timingInstability.TimingInstabilityPeakToPeakNs, [], "omitnan");
fitInfo.MaxTimingInstabilityMaxAbsNs = max(timingInstability.TimingInstabilityMaxAbsNs, [], "omitnan");

if fitInfo.ConditionNumber > opts.MaxConditionNumber || ~isfinite(fitInfo.ConditionNumber)
    fitInfo.Status = "SUSPECT";
    fitInfo.Reason = "Survey geometry is ill-conditioned; timing estimates may not be unique.";
elseif degreesOfFreedom == 0
    fitInfo.Status = "SUSPECT";
    fitInfo.Reason = "Survey fit is exactly determined. Offsets were estimated, but no redundant geometry measurement remains to validate the fit residual.";
elseif fitInfo.FitRMSNs > opts.MaxFitRMSNsForPass
    fitInfo.Status = "SUSPECT";
    fitInfo.Reason = "Model residual is high; multipath, bad detections, or weak survey geometry may be present.";
else
    fitInfo.Status = "OK";
    fitInfo.Reason = "Survey fit completed.";
end

timingEstimates = estimatesFromStateLocal(bestX, problem, referencePCI, ...
    biasStdNs, positionStdM, fitInfo, opts, timingInstability);
if ~isempty(excludedPCIs)
    excludedReason = "Excluded from the geometry fit because fewer than the structurally required distinct receiver positions were available for this PCI.";
    timingEstimates = [timingEstimates; ...
        notAssessableRowsLocal(excludedPCIs, referencePCI, excludedReason)];
    timingEstimates = orderTimingRowsLocal(timingEstimates, visiblePCIs);
end
fitInfo.MaxWorstCaseRelativeTimingNs = max(timingEstimates.WorstCaseRelativeTimingNs, [], "omitnan");

end

function problem = buildProblemLocal(rel, pcisAll, pcisTarget, opts)
numG = numel(pcisAll);
targetIndex = zeros(height(rel), 1);
for k = 1:height(rel)
    targetIndex(k) = find(pcisAll == rel.PCI(k), 1);
end

sigmaRxNs = sqrt(2) * rel.RxPositionUncertaintyM ./ opts.SpeedOfLightMps * 1e9;
sigmaNs = hypot(rel.SigmaNs, sigmaRxNs);
sigmaNs = max(sigmaNs, opts.TimingStdFloorNs);

problem = struct();
problem.pcisAll = pcisAll(:);
problem.pcisTarget = pcisTarget(:);
problem.numG = numG;
problem.numTargets = numel(pcisTarget);
problem.targetIndex = targetIndex;
problem.rxXY = [rel.RxX_m(:), rel.RxY_m(:)];
problem.yNs = rel.RelativeArrivalNs(:);
problem.sigmaNs = sigmaNs(:);
problem.c = opts.SpeedOfLightMps;
problem.frameNs = opts.FramePeriodMs * 1e6;
problem.rxCenter = mean(problem.rxXY, 1);
end

function [bestX, bestResidual, bestCost, optimizerName, diagnostics] = ...
        solveProblemLocal(problem, opts)
rng(opts.RandomSeed, "twister");
starts = initialStatesLocal(problem, opts);
bestCost = Inf;
bestX = starts(:,1);
bestResidual = residualUnweightedLocal(bestX, problem);
numStarts = size(starts,2);
allCosts = Inf(numStarts,1);
allOffsetsNs = NaN(numStarts,problem.numTargets);

useLsq = opts.UseLsqnonlinIfAvailable && exist("lsqnonlin", "file") == 2;
if useLsq
    optimizerName = "lsqnonlin";
    lsqOpts = optimoptions("lsqnonlin", ...
        "Display", "off", ...
        "MaxIterations", opts.MaxIterations, ...
        "FunctionTolerance", opts.FunctionTolerance);
else
    optimizerName = "fminsearch";
    fmOpts = optimset("Display", "off", ...
        "MaxIter", opts.MaxIterations, ...
        "MaxFunEvals", opts.MaxIterations * max(10, numel(bestX)), ...
        "TolFun", opts.FunctionTolerance);
end

for s = 1:size(starts, 2)
    x0 = starts(:,s);
    if useLsq
        x = lsqnonlin(@(x) residualWeightedLocal(x, problem), x0, [], [], lsqOpts);
    else
        x = fminsearch(@(x) sum(residualWeightedLocal(x, problem).^2), x0, fmOpts);
    end

    r = residualUnweightedLocal(x, problem);
    cost = sum((r ./ problem.sigmaNs).^2);
    [~,biasNs] = unpackStateLocal(x,problem);
    allCosts(s) = cost;
    allOffsetsNs(s,:) = (biasNs(2:end)-biasNs(1)).';
    if cost < bestCost
        bestCost = cost;
        bestX = x;
        bestResidual = r;
    end
end

valid = isfinite(allCosts) & all(isfinite(allOffsetsNs),2);
costTolerance = max(opts.ConsensusAbsoluteCostTolerance, ...
    opts.ConsensusRelativeCostTolerance*max(1,bestCost));
nearOptimal = valid & allCosts <= bestCost+costTolerance;
if ~any(nearOptimal) && any(valid)
    [~,idx] = min(allCosts);
    nearOptimal(idx) = true;
end

nearOffsets = allOffsetsNs(nearOptimal,:);
if isempty(nearOffsets)
    spreadNs = NaN(problem.numTargets,1);
    stdNs = NaN(problem.numTargets,1);
    medianNs = NaN(problem.numTargets,1);
else
    spreadNs = (max(nearOffsets,[],1)-min(nearOffsets,[],1)).';
    stdNs = std(nearOffsets,0,1).';
    medianNs = median(nearOffsets,1).';
end

diagnostics = struct();
diagnostics.Statement = ...
    "Reference-anchored offset agreement across independent optimizer starts. Large spread among near-equal-cost starts indicates that the survey does not uniquely determine the offsets.";
diagnostics.TargetPCIs = problem.pcisTarget;
diagnostics.NumStarts = numStarts;
diagnostics.NumValidStarts = nnz(valid);
diagnostics.NumNearOptimalStarts = nnz(nearOptimal);
diagnostics.BestCost = bestCost;
diagnostics.CostTolerance = costTolerance;
diagnostics.AllCosts = allCosts;
diagnostics.NearOptimalMask = nearOptimal;
diagnostics.AllReferenceAnchoredOffsetsNs = allOffsetsNs;
diagnostics.PerTarget = table(problem.pcisTarget,medianNs,stdNs,spreadNs, ...
    'VariableNames',["PCI","NearOptimalMedianOffsetNs", ...
    "NearOptimalOffsetStdNs","NearOptimalOffsetSpreadNs"]);
diagnostics.MaxNearOptimalOffsetSpreadNs = max(spreadNs,[],"omitnan");
end

function starts = initialStatesLocal(problem, opts)
numG = problem.numG;
numParams = 2*numG;
numStarts = max(1, opts.NumRandomStarts);
starts = zeros(numParams, numStarts);

baseAngles = linspace(0, 2*pi, numG+1).';
baseAngles(end) = [];
for s = 1:numStarts
    if s == 1
        radius = 0.35 * opts.PositionSearchRadiusM;
        angles = baseAngles;
    elseif s == 2
        radius = 0.75 * opts.PositionSearchRadiusM;
        angles = baseAngles + pi/numG;
    else
        radius = opts.PositionSearchRadiusM * (0.15 + 0.85*rand(numG,1));
        angles = 2*pi*rand(numG,1);
    end

    if isscalar(radius)
        radius = radius * ones(numG,1);
    end
    xy = problem.rxCenter + [radius(:).*cos(angles(:)), radius(:).*sin(angles(:))];

    x0 = zeros(numParams, 1);
    x0(1:numG) = xy(:,1);
    x0(numG+1:2*numG) = xy(:,2);
    starts(:,s) = x0;
end
end

function r = residualWeightedLocal(x, problem)
r = residualUnweightedLocal(x, problem) ./ problem.sigmaNs;
end

function r = residualUnweightedLocal(x, problem)
[xy, biasNs] = unpackStateLocal(x, problem);
refXY = xy(1,:);
targetXY = xy(problem.targetIndex,:);
rxXY = problem.rxXY;

distTarget = sqrt(sum((targetXY - rxXY).^2, 2));
distRef = sqrt(sum((refXY - rxXY).^2, 2));

targetBias = biasNs(problem.targetIndex);
predNs = targetBias + (distTarget - distRef) ./ problem.c * 1e9;
r = wrapNsLocal(problem.yNs - predNs, problem.frameNs);
end

function [xy, biasNs] = unpackStateLocal(x, problem)
numG = problem.numG;
xy = [x(1:numG), x(numG+1:2*numG)];
biasNs = zeros(numG, 1);
if numel(x) > 2*numG
    biasNs(2:end) = x(2*numG+1:end);
else
    biasNs = computeBestBiasesLocal(xy, problem);
end
end

function biasNs = computeBestBiasesLocal(xy, problem)
numG = problem.numG;
biasNs = zeros(numG, 1);
refXY = xy(1,:);
for g = 2:numG
    idx = problem.targetIndex == g;
    if ~any(idx)
        biasNs(g) = NaN;
        continue;
    end

    rxXY = problem.rxXY(idx,:);
    distTarget = sqrt(sum((xy(g,:) - rxXY).^2, 2));
    distRef = sqrt(sum((refXY - rxXY).^2, 2));
    geomNs = (distTarget - distRef) ./ problem.c * 1e9;
    samples = wrapNsLocal(problem.yNs(idx) - geomNs, problem.frameNs);
    weights = 1 ./ (problem.sigmaNs(idx).^2);
    biasNs(g) = sum(weights .* samples) / sum(weights);
end
end

function [biasStdNs, positionStdM, conditionNumber, diagnostics] = ...
        estimateUncertaintyLocal(x, problem, residualNs)
[xy, biasNs] = unpackStateLocal(x, problem);
fullX = [xy(:,1); xy(:,2); biasNs(2:end)];
numParams = numel(fullX);
numMeas = numel(residualNs);
diagnostics = emptyNumericalDiagnosticsLocal(problem,numParams,numMeas);
if numMeas < numParams
    biasStdNs = NaN(problem.numTargets, 1);
    positionStdM = NaN(problem.numG, 1);
    conditionNumber = Inf;
    diagnostics.ConditionNumber = Inf;
    diagnostics.Statement = ...
        "The survey has fewer measurements than fitted parameters.";
    return;
end

% Positions are parameterized in propagation-delay-equivalent nanoseconds:
% q = p*(1e9/c). This makes position and timing-offset columns comparable
% and avoids a condition number that changes merely because metres are
% replaced by kilometres.
Jscaled = analyticScaledJacobianLocal(fullX,problem);
Jw = Jscaled ./ problem.sigmaNs;
[~,S,V] = svd(Jw,0);
singularValues = diag(S);
if isempty(singularValues)
    conditionNumber = Inf;
    rankTolerance = NaN;
    numericalRank = 0;
else
    rankTolerance = max(size(Jw))*eps(max(singularValues));
    numericalRank = nnz(singularValues>rankTolerance);
    if singularValues(end)<=0
        conditionNumber = Inf;
    else
        conditionNumber = singularValues(1)/singularValues(end);
    end
end
degreesOfFreedom = numMeas - numParams;
if degreesOfFreedom > 0
    sigma2 = max(1, sum((residualNs ./ problem.sigmaNs).^2) / degreesOfFreedom);
else
    sigma2 = 1;
end

normalMatrix = Jw.'*Jw;
covScaled = pinv(normalMatrix) * sigma2;
stdScaled = sqrt(max(diag(covScaled), 0));

numG = problem.numG;
nsToM = problem.c/1e9;
positionStdM = hypot(stdScaled(1:numG), ...
    stdScaled(numG+1:2*numG))*nsToM;
biasStdNs = stdScaled(2*numG+1:end);

diagnostics.Statement = ...
    "Analytic weighted Jacobian using delay-equivalent nanosecond scaling for position and timing parameters.";
diagnostics.ParameterNames = parameterNamesLocal(problem);
diagnostics.SingularValues = singularValues;
diagnostics.NumericalRank = numericalRank;
diagnostics.RankTolerance = rankTolerance;
diagnostics.FullColumnRank = numericalRank==numParams;
diagnostics.ConditionNumber = conditionNumber;
diagnostics.NormalMatrixRcond = rcond(normalMatrix);
if ~isempty(V)
    diagnostics.WeakestRightSingularVector = V(:,end);
else
    diagnostics.WeakestRightSingularVector = NaN(numParams,1);
end
end

function [fitMeasurements, support, excludedPCIs] = selectEstimableTargetsLocal(rel, opts)
pcis = unique(rel.PCI, "stable");
numMeasurements = zeros(numel(pcis), 1);
numPositions = zeros(numel(pcis), 1);
for k = 1:numel(pcis)
    rows = rel(rel.PCI == pcis(k), :);
    numMeasurements(k) = height(rows);
    numPositions(k) = height(unique(rows(:,["RxX_m","RxY_m"]),"rows"));
end
minimumPositions = max(3, opts.MinDistinctReceiverPositionsPerTarget);
minimumMeasurements = max(minimumPositions, opts.MinMeasurementsPerPCI);
isEligible = numPositions >= minimumPositions & ...
    numMeasurements >= minimumMeasurements;
support = table(pcis(:), numMeasurements, numPositions, isEligible, ...
    'VariableNames', ["PCI","NumMeasurements","NumDistinctPositions", ...
    "IncludedInGeometryFit"]);
includedPCIs = support.PCI(support.IncludedInGeometryFit);
excludedPCIs = support.PCI(~support.IncludedInGeometryFit);
fitMeasurements = rel(ismember(rel.PCI, includedPCIs), :);
end

function J = analyticScaledJacobianLocal(fullX,problem)
[xy,~] = unpackStateLocal(fullX,problem);
numRows = numel(problem.yNs);
numG = problem.numG;
J = zeros(numRows,2*numG+problem.numTargets);
refXY = xy(1,:);
for row = 1:numRows
    g = problem.targetIndex(row);
    rx = problem.rxXY(row,:);
    refDelta = refXY-rx;
    targetDelta = xy(g,:)-rx;
    refDistance = max(norm(refDelta),eps);
    targetDistance = max(norm(targetDelta),eps);
    refUnit = refDelta/refDistance;
    targetUnit = targetDelta/targetDistance;

    % Residual = measured - predicted. In delay-equivalent position units,
    % the propagation scale 1e9/c cancels.
    J(row,1) = refUnit(1);
    J(row,numG+1) = refUnit(2);
    J(row,g) = J(row,g)-targetUnit(1);
    J(row,numG+g) = J(row,numG+g)-targetUnit(2);
    J(row,2*numG+g-1) = -1;
end
end

function names = parameterNamesLocal(problem)
names = strings(2*problem.numG+problem.numTargets,1);
for g = 1:problem.numG
    names(g) = "PCI"+string(problem.pcisAll(g))+"_XdelayNs";
    names(problem.numG+g) = "PCI"+string(problem.pcisAll(g))+"_YdelayNs";
end
for g = 1:problem.numTargets
    names(2*problem.numG+g) = ...
        "PCI"+string(problem.pcisTarget(g))+"_RelativeOffsetNs";
end
end

function diagnostics = emptyNumericalDiagnosticsLocal(problem,numParams,numMeas)
diagnostics = struct();
diagnostics.Statement = "";
diagnostics.ParameterNames = parameterNamesLocal(problem);
diagnostics.NumMeasurements = numMeas;
diagnostics.NumParameters = numParams;
diagnostics.SingularValues = NaN(min(numMeas,numParams),1);
diagnostics.NumericalRank = 0;
diagnostics.RankTolerance = NaN;
diagnostics.FullColumnRank = false;
diagnostics.ConditionNumber = NaN;
diagnostics.NormalMatrixRcond = NaN;
diagnostics.WeakestRightSingularVector = NaN(numParams,1);
end

function tbl = estimatesFromStateLocal(x, problem, referencePCI, biasStdNs, positionStdM, fitInfo, opts, timingInstability)
[xy, biasNs] = unpackStateLocal(x, problem);
pcisAll = problem.pcisAll;
centerNs = median(biasNs, "omitnan");
centeredBiasNs = biasNs - centerNs;

biasStdAll = NaN(numel(pcisAll), 1);
biasStdAll(2:end) = biasStdNs;
if all(~isfinite(biasStdNs))
    biasStdAll(1) = NaN;
else
    biasStdAll(1) = median(biasStdNs, "omitnan");
end

rows = cell(numel(pcisAll), 1);
for k = 1:numel(pcisAll)
    pci = pcisAll(k);
    if pci == referencePCI
        role = "GAUGE_REFERENCE";
    else
        role = "TARGET";
    end

    relOffset = centeredBiasNs(k);
    biasStd = biasStdAll(k);
    [staticStatus, staticReason] = classifySurveyTimingLocal(relOffset, biasStd, fitInfo, opts);
    instability = timingInstability(timingInstability.PCI == pci, :);

    rows{k} = table( ...
        pci, role, referencePCI, ...
        relOffset, biasStd, xy(k,1), xy(k,2), positionStdM(k), ...
        staticStatus, staticReason, ...
        instability.TimingInstabilitySamples, ...
        instability.TimingInstabilityRmsNs, ...
        instability.TimingInstabilityPeakToPeakNs, ...
        instability.TimingInstabilityMaxAbsNs, ...
        instability.ExpectedTimingNoiseRmsNs, ...
        instability.ExcessTimingJitterRmsNs, ...
        instability.TimingInstabilityStatus, ...
        instability.TimingInstabilityReason, ...
        abs(relOffset), ...
        abs(relOffset) + instability.TimingInstabilityMaxAbsNs, ...
        combinedEnvelopeUncertaintyLocal(biasStd, instability.ExpectedTimingNoiseRmsNs), ...
        "", ...
        "", ...
        "", "", ...
        fitInfo.FitRMSNs, fitInfo.ConditionNumber, ...
        'VariableNames', ["PCI","SurveyRole","ReferencePCI", ...
        "EstimatedRelativeTxOffsetNs","EstimatedOffsetUncertaintyNs", ...
        "EstimatedX_m","EstimatedY_m","EstimatedPositionUncertaintyM", ...
        "StaticTimingStatus","StaticTimingReason", ...
        "TimingInstabilitySamples","TimingInstabilityRmsNs", ...
        "TimingInstabilityPeakToPeakNs","TimingInstabilityMaxAbsNs", ...
        "ExpectedTimingNoiseRmsNs", ...
        "ExcessTimingJitterRmsNs","TimingInstabilityStatus", ...
        "TimingInstabilityReason","StaticOffsetNs", ...
        "WorstCaseRelativeTimingNs","WorstCaseRelativeTimingUncertaintyNs", ...
        "EnvelopeTimingStatus","EnvelopeTimingReason", ...
        "TimingStatus","TimingReason","FitRMSNs","ConditionNumber"]);
    [rows{k}.EnvelopeTimingStatus, rows{k}.EnvelopeTimingReason] = ...
        classifyCombinedEnvelopeLocal(rows{k}, fitInfo, opts);
    [rows{k}.TimingStatus, rows{k}.TimingReason] = combineTimingStatusLocal( ...
        rows{k}.StaticTimingStatus, rows{k}.StaticTimingReason, ...
        rows{k}.TimingInstabilityStatus, rows{k}.TimingInstabilityReason, ...
        rows{k}.EnvelopeTimingStatus, rows{k}.EnvelopeTimingReason, role, opts);
end
tbl = vertcat(rows{:});
end

function tbl = timingInstabilityFromResidualsLocal(residualNs, problem, opts)
% Estimate relative timing instability after static timing and geometry fit.
% Residuals are reference-subtracted time errors; their per-PCI scatter is
% the observable relative instability. A common-mode receiver or network
% timing movement remains invisible in this relative mode.
pcisAll = problem.pcisAll(:);
num = numel(pcisAll);
samples = NaN(num, 1);
rmsNs = NaN(num, 1);
peakToPeakNs = NaN(num, 1);
maxAbsNs = NaN(num, 1);
expectedNoiseNs = NaN(num, 1);
excessJitterNs = NaN(num, 1);
status = repmat("NOT_ASSESSABLE", num, 1);
reason = repmat("No target-relative residuals are available.", num, 1);

if num > 0
    reason(1) = "Reference PCI is the subtraction gauge; target rows carry relative instability evidence.";
end

for g = 2:num
    idx = problem.targetIndex == g;
    res = residualNs(idx);
    sigma = problem.sigmaNs(idx);
    samples(g) = numel(res);

    if numel(res) < opts.MinInstabilitySamplesPerPCI
        status(g) = "NOT_ASSESSABLE";
        reason(g) = "Too few residual samples to assess timing instability.";
        continue;
    end

    weights = 1 ./ max(sigma(:), opts.TimingStdFloorNs).^2;
    meanResidual = sum(weights .* res(:)) ./ sum(weights);
    centered = res(:) - meanResidual;

    rmsNs(g) = sqrt(mean(centered.^2));
    peakToPeakNs(g) = max(centered) - min(centered);
    maxAbsNs(g) = max(abs(centered));
    expectedNoiseNs(g) = sqrt(mean(sigma(:).^2));
    excessJitterNs(g) = sqrt(max(0, rmsNs(g)^2 - expectedNoiseNs(g)^2));

    [status(g), reason(g)] = classifyTimingInstabilityLocal( ...
        samples(g), excessJitterNs(g), peakToPeakNs(g), expectedNoiseNs(g), opts);
end

tbl = table(pcisAll, samples, rmsNs, peakToPeakNs, maxAbsNs, expectedNoiseNs, ...
    excessJitterNs, status, reason, ...
    'VariableNames', ["PCI","TimingInstabilitySamples", ...
    "TimingInstabilityRmsNs","TimingInstabilityPeakToPeakNs", ...
    "TimingInstabilityMaxAbsNs", ...
    "ExpectedTimingNoiseRmsNs","ExcessTimingJitterRmsNs", ...
    "TimingInstabilityStatus","TimingInstabilityReason"]);
end

function [status, reason] = classifyTimingInstabilityLocal(numSamples, excessRmsNs, peakToPeakNs, expectedNoiseRmsNs, opts)
if numSamples < opts.MinInstabilitySamplesPerPCI
    status = "NOT_ASSESSABLE";
    reason = "Too few residual samples to assess timing instability.";
    return;
end

rmsMargin = opts.UncertaintySigmaMultiplier * expectedNoiseRmsNs / sqrt(max(1, numSamples));
peakMargin = opts.UncertaintySigmaMultiplier * expectedNoiseRmsNs * sqrt(2);

if isfinite(opts.HardFailInstabilityPeakToPeakNs) && ...
        peakToPeakNs >= opts.HardFailInstabilityPeakToPeakNs
    status = "FAIL";
    reason = "Relative timing instability exceeds the hard-fail peak-to-peak threshold.";
elseif isfinite(opts.FailInstabilityRmsNs) && ...
        excessRmsNs - rmsMargin >= opts.FailInstabilityRmsNs
    status = "FAIL";
    reason = "Relative timing instability RMS exceeds the fail threshold with uncertainty margin.";
elseif isfinite(opts.FailInstabilityPeakToPeakNs) && ...
        peakToPeakNs - peakMargin >= opts.FailInstabilityPeakToPeakNs
    status = "FAIL";
    reason = "Relative timing instability peak-to-peak exceeds the fail threshold with uncertainty margin.";
elseif isfinite(opts.FailInstabilityRmsNs) && ...
        excessRmsNs + rmsMargin >= opts.FailInstabilityRmsNs
    status = "SUSPECT";
    reason = "Relative timing instability RMS is near the fail threshold after uncertainty handling.";
elseif isfinite(opts.FailInstabilityPeakToPeakNs) && ...
        peakToPeakNs + peakMargin >= opts.FailInstabilityPeakToPeakNs
    status = "SUSPECT";
    reason = "Relative timing instability peak-to-peak is near the fail threshold after uncertainty handling.";
elseif isfinite(opts.WarningInstabilityRmsNs) && ...
        excessRmsNs - rmsMargin >= opts.WarningInstabilityRmsNs
    status = "SUSPECT";
    reason = "Relative timing instability RMS exceeds the warning threshold with uncertainty margin.";
elseif isfinite(opts.WarningInstabilityPeakToPeakNs) && ...
        peakToPeakNs - peakMargin >= opts.WarningInstabilityPeakToPeakNs
    status = "SUSPECT";
    reason = "Relative timing instability peak-to-peak exceeds the warning threshold with uncertainty margin.";
elseif isfinite(opts.WarningInstabilityRmsNs) && ...
        excessRmsNs + rmsMargin >= opts.WarningInstabilityRmsNs
    status = "SUSPECT";
    reason = "Relative timing instability RMS is near the warning threshold after uncertainty handling.";
elseif isfinite(opts.WarningInstabilityPeakToPeakNs) && ...
        peakToPeakNs + peakMargin >= opts.WarningInstabilityPeakToPeakNs
    status = "SUSPECT";
    reason = "Relative timing instability peak-to-peak is near the warning threshold after uncertainty handling.";
elseif thresholdsAreSatisfiedLocal(excessRmsNs + rmsMargin, opts.FailInstabilityRmsNs) && ...
        thresholdsAreSatisfiedLocal(peakToPeakNs + peakMargin, opts.FailInstabilityPeakToPeakNs)
    status = "PASS";
    reason = "No relative timing instability reaches the configured desynchronization thresholds.";
else
    status = "SUSPECT";
    reason = "Relative timing instability is near threshold or uncertainty is too high for a strict pass.";
end
end

function [status, reason] = combineTimingStatusLocal(staticStatus, staticReason, instabilityStatus, instabilityReason, envelopeStatus, envelopeReason, role, opts)
if role == "GAUGE_REFERENCE"
    status = staticStatus;
    reason = staticReason;
    return;
end

if envelopeStatus == "FAIL"
    status = "FAIL";
    reason = envelopeReason;
elseif instabilityStatus == "FAIL"
    status = "FAIL";
    reason = "Relative timing instability fault detected. " + instabilityReason;
elseif staticStatus == "FAIL"
    status = "FAIL";
    reason = staticReason;
elseif envelopeStatus == "SUSPECT"
    status = "SUSPECT";
    reason = envelopeReason;
elseif instabilityStatus == "SUSPECT"
    status = "SUSPECT";
    reason = "Relative timing instability is suspect. " + instabilityReason;
elseif staticStatus == "SUSPECT"
    status = "SUSPECT";
    reason = staticReason;
elseif instabilityStatus == "NOT_ASSESSABLE" && opts.RequireInstabilityAssessableForPass
    status = "SUSPECT";
    reason = "Static timing is acceptable, but relative timing instability was not assessable. " + instabilityReason;
else
    status = envelopeStatus;
    reason = envelopeReason;
end
end

function [status, reason] = classifyCombinedEnvelopeLocal(row, fitInfo, opts)
if row.SurveyRole == "GAUGE_REFERENCE"
    status = row.StaticTimingStatus;
    reason = "Gauge reference only; combined target envelope is not interpreted as a target verdict.";
    return;
end

if fitInfo.Status == "NOT_ASSESSABLE"
    status = "NOT_ASSESSABLE";
    reason = fitInfo.Reason;
    return;
end

if row.TimingInstabilityStatus == "NOT_ASSESSABLE" && opts.RequireInstabilityAssessableForPass
    status = "SUSPECT";
    reason = "Worst-case relative timing envelope cannot be cleanly assessed because residual instability was not assessable. " + row.TimingInstabilityReason;
    return;
end

if ~isfinite(row.WorstCaseRelativeTimingNs)
    status = "SUSPECT";
    reason = "Worst-case relative timing envelope is unavailable.";
    return;
end

envelope = row.WorstCaseRelativeTimingNs;
margin = opts.UncertaintySigmaMultiplier * row.WorstCaseRelativeTimingUncertaintyNs;
if isfinite(opts.HardFailOffsetNs) && envelope >= opts.HardFailOffsetNs
    status = "FAIL";
    reason = "Worst-case relative timing envelope exceeds the hard-fail threshold.";
elseif envelope - margin >= opts.FailOffsetNs
    status = "FAIL";
    reason = "Worst-case relative timing envelope exceeds the configured synchronization threshold with uncertainty margin.";
elseif envelope + margin >= opts.FailOffsetNs
    status = "SUSPECT";
    reason = "Worst-case relative timing envelope is near the configured synchronization threshold after uncertainty handling.";
elseif isfield(fitInfo, "NumGNBs") && fitInfo.NumGNBs < opts.MinCellsForCleanGroupPass
    status = "SUSPECT";
    reason = "Only two cells are available. Pairwise relative timing can be estimated, but the group-centred pass/fail attribution is not defensible without a third cell.";
elseif fitInfo.FitRMSNs > opts.MaxFitRMSNsForPass || ...
        fitInfo.ConditionNumber > opts.MaxConditionNumber || ...
        ~isfinite(fitInfo.ConditionNumber)
    status = "SUSPECT";
    reason = "Worst-case relative timing envelope is below threshold, but survey fit quality is weak.";
elseif row.EstimatedOffsetUncertaintyNs <= opts.MaxBiasUncertaintyForPassNs
    status = "PASS";
    reason = "Worst-case relative timing envelope is below the configured synchronization threshold with acceptable uncertainty.";
else
    status = "SUSPECT";
    reason = "Worst-case relative timing envelope is below threshold, but timing uncertainty is too high for a clean pass.";
end
end

function sigmaNs = combinedEnvelopeUncertaintyLocal(offsetUncertaintyNs, expectedNoiseNs)
sigmaNs = hypot(offsetUncertaintyNs, expectedNoiseNs);
if ~isfinite(sigmaNs)
    sigmaNs = offsetUncertaintyNs;
end
if ~isfinite(sigmaNs)
    sigmaNs = expectedNoiseNs;
end
end

function [status, reason] = classifySurveyTimingLocal(offsetNs, uncertaintyNs, fitInfo, opts)
if fitInfo.Status == "NOT_ASSESSABLE"
    status = "NOT_ASSESSABLE";
    reason = fitInfo.Reason;
    return;
end

if fitInfo.FitRMSNs > opts.MaxFitRMSNsForPass || ...
        fitInfo.ConditionNumber > opts.MaxConditionNumber || ...
        ~isfinite(fitInfo.ConditionNumber)
    status = "SUSPECT";
    reason = "Survey fit quality is weak, so timing is not strict enough for PASS/FAIL.";
    return;
end

offset = abs(offsetNs);
margin = opts.UncertaintySigmaMultiplier * uncertaintyNs;
if isfinite(opts.HardFailOffsetNs) && offset >= opts.HardFailOffsetNs
    status = "FAIL";
    reason = "Estimated relative transmit timing exceeds hard-fail threshold.";
elseif offset - margin >= opts.FailOffsetNs
    status = "FAIL";
    reason = "Estimated relative transmit timing exceeds the configured synchronization threshold with uncertainty margin.";
elseif offset + margin >= opts.FailOffsetNs
    status = "SUSPECT";
    reason = "Estimated relative transmit timing is near the configured synchronization threshold after uncertainty handling.";
elseif isfinite(opts.WarningOffsetNs) && offset - margin >= opts.WarningOffsetNs
    status = "SUSPECT";
    reason = "Estimated relative transmit timing exceeds the configured warning threshold with uncertainty margin.";
elseif isfinite(opts.WarningOffsetNs) && offset + margin >= opts.WarningOffsetNs
    status = "SUSPECT";
    reason = "Estimated relative transmit timing is near the configured warning threshold after uncertainty handling.";
elseif uncertaintyNs <= opts.MaxBiasUncertaintyForPassNs
    status = "PASS";
    reason = "Estimated relative transmit timing is below the configured synchronization threshold and survey uncertainty is acceptable.";
else
    status = "SUSPECT";
    reason = "Estimated timing is near threshold or survey uncertainty is too high for a strict pass.";
end
end

function fitInfo = baseFitInfoLocal(opts, surveyInfo, status, reason)
fitInfo = struct();
fitInfo.Mode = "multi_position_unknown_gnb_locations";
fitInfo.Statement = "No gNB locations are required, but several captures from known receiver positions are required. Absolute UTC phase is still not assessed.";
fitInfo.ReferencePCI = surveyInfo.ReferencePCI;
fitInfo.NumCaptures = surveyInfo.NumCaptures;
fitInfo.NumCapturesWithReference = surveyInfo.NumCapturesWithReference;
fitInfo.NumRelativeMeasurements = surveyInfo.NumRelativeMeasurements;
fitInfo.MissingReferenceCaptureIDs = surveyInfo.MissingReferenceCaptureIDs;
fitInfo.MinCaptures = opts.MinCaptures;
fitInfo.WarningOffsetNs = opts.WarningOffsetNs;
fitInfo.FailOffsetNs = opts.FailOffsetNs;
fitInfo.HardFailOffsetNs = opts.HardFailOffsetNs;
fitInfo.SyncThresholdNs = opts.FailOffsetNs;
fitInfo.UncertaintySigmaMultiplier = opts.UncertaintySigmaMultiplier;
fitInfo.WarningInstabilityRmsNs = opts.WarningInstabilityRmsNs;
fitInfo.FailInstabilityRmsNs = opts.FailInstabilityRmsNs;
fitInfo.WarningInstabilityPeakToPeakNs = opts.WarningInstabilityPeakToPeakNs;
fitInfo.FailInstabilityPeakToPeakNs = opts.FailInstabilityPeakToPeakNs;
fitInfo.HardFailInstabilityPeakToPeakNs = opts.HardFailInstabilityPeakToPeakNs;
fitInfo.MaxFitRMSNsForPass = opts.MaxFitRMSNsForPass;
fitInfo.MaxBiasUncertaintyForPassNs = opts.MaxBiasUncertaintyForPassNs;
fitInfo.FitRMSPassFraction = opts.FitRMSPassFraction;
fitInfo.BiasUncertaintyPassFraction = opts.BiasUncertaintyPassFraction;
fitInfo.Status = status;
fitInfo.Reason = reason;
fitInfo.NumGNBs = 0;
fitInfo.NumVisibleGNBs = 0;
fitInfo.NumTargetGNBs = 0;
fitInfo.NumCapturesUsed = 0;
fitInfo.NumMeasurements = 0;
fitInfo.NumUnknowns = 0;
fitInfo.DegreesOfFreedom = NaN;
fitInfo.TargetSupport = table();
fitInfo.ExcludedPCIs = zeros(0,1);
fitInfo.Optimizer = "";
fitInfo.Cost = NaN;
fitInfo.FitRMSNs = NaN;
fitInfo.WeightedRMS = NaN;
fitInfo.ConditionNumber = NaN;
fitInfo.ConditionNumberBasis = "";
fitInfo.NumericalDiagnostics = struct();
fitInfo.StartDiagnostics = struct();
fitInfo.MaxNearOptimalOffsetSpreadNs = NaN;
fitInfo.MaxTimingInstabilityRmsNs = NaN;
fitInfo.MaxExcessTimingJitterRmsNs = NaN;
fitInfo.MaxTimingInstabilityPeakToPeakNs = NaN;
fitInfo.MaxTimingInstabilityMaxAbsNs = NaN;
fitInfo.MaxWorstCaseRelativeTimingNs = NaN;
end

function tbl = notAssessableRowsLocal(pcisAll, referencePCI, reason)
num = numel(pcisAll);
role = repmat("TARGET", num, 1);
role(pcisAll == referencePCI) = "REFERENCE";
tbl = table( ...
    pcisAll(:), role, repmat(referencePCI, num, 1), ...
    NaN(num,1), NaN(num,1), NaN(num,1), NaN(num,1), NaN(num,1), ...
    repmat("NOT_ASSESSABLE", num, 1), repmat(string(reason), num, 1), ...
    NaN(num,1), NaN(num,1), NaN(num,1), NaN(num,1), NaN(num,1), NaN(num,1), ...
    repmat("NOT_ASSESSABLE", num, 1), repmat(string(reason), num, 1), ...
    NaN(num,1), NaN(num,1), NaN(num,1), ...
    repmat("NOT_ASSESSABLE", num, 1), repmat(string(reason), num, 1), ...
    repmat("NOT_ASSESSABLE", num, 1), repmat(string(reason), num, 1), ...
    NaN(num,1), NaN(num,1), ...
    'VariableNames', ["PCI","SurveyRole","ReferencePCI", ...
    "EstimatedRelativeTxOffsetNs","EstimatedOffsetUncertaintyNs", ...
    "EstimatedX_m","EstimatedY_m","EstimatedPositionUncertaintyM", ...
    "StaticTimingStatus","StaticTimingReason", ...
    "TimingInstabilitySamples","TimingInstabilityRmsNs", ...
    "TimingInstabilityPeakToPeakNs","TimingInstabilityMaxAbsNs", ...
    "ExpectedTimingNoiseRmsNs", ...
    "ExcessTimingJitterRmsNs","TimingInstabilityStatus", ...
    "TimingInstabilityReason","StaticOffsetNs", ...
    "WorstCaseRelativeTimingNs","WorstCaseRelativeTimingUncertaintyNs", ...
    "EnvelopeTimingStatus","EnvelopeTimingReason", ...
    "TimingStatus","TimingReason","FitRMSNs","ConditionNumber"]);
end

function tbl = orderTimingRowsLocal(tbl, pciOrder)
order = zeros(height(tbl),1);
for k = 1:height(tbl)
    idx = find(pciOrder == tbl.PCI(k), 1);
    if isempty(idx)
        idx = numel(pciOrder) + k;
    end
    order(k) = idx;
end
[~,idx] = sort(order);
tbl = tbl(idx,:);
end

function tbl = emptyTimingEstimateTableLocal()
tbl = table( ...
    zeros(0,1), strings(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), strings(0,1), strings(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    strings(0,1), strings(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    strings(0,1), strings(0,1), strings(0,1), strings(0,1), ...
    zeros(0,1), zeros(0,1), ...
    'VariableNames', ["PCI","SurveyRole","ReferencePCI", ...
    "EstimatedRelativeTxOffsetNs","EstimatedOffsetUncertaintyNs", ...
    "EstimatedX_m","EstimatedY_m","EstimatedPositionUncertaintyM", ...
    "StaticTimingStatus","StaticTimingReason", ...
    "TimingInstabilitySamples","TimingInstabilityRmsNs", ...
    "TimingInstabilityPeakToPeakNs","TimingInstabilityMaxAbsNs", ...
    "ExpectedTimingNoiseRmsNs", ...
    "ExcessTimingJitterRmsNs","TimingInstabilityStatus", ...
    "TimingInstabilityReason","StaticOffsetNs", ...
    "WorstCaseRelativeTimingNs","WorstCaseRelativeTimingUncertaintyNs", ...
    "EnvelopeTimingStatus","EnvelopeTimingReason", ...
    "TimingStatus","TimingReason","FitRMSNs","ConditionNumber"]);
end

function y = wrapNsLocal(x, periodNs)
y = mod(x + periodNs/2, periodNs) - periodNs/2;
end

function [relativeMeasurements, surveyInfo] = prepareSurveyRelativeMeasurements(measurementTable, captureInfo, opts)
% Convert per-capture PCI timing into reference-subtracted TDOA rows. The
% subtraction removes each capture's arbitrary B210 internal-clock phase.
opts = mergeStructsLocal(defaultSurveyOptions(), opts);
validateSurveyInputsLocal(measurementTable, captureInfo);

captureXY = receiverXYLocal(captureInfo, opts);
meas = measurementTable;
meas.CaptureID = string(meas.CaptureID);
meas.PCI = double(meas.PCI);

[isValidInput, validationReason, badCaptureIDs] = validateReceiverCoverageLocal(meas, captureXY);

if ~ismember("TimingStdNs", string(meas.Properties.VariableNames))
    meas.TimingStdNs = opts.TimingStdNsDefault * ones(height(meas), 1);
end
meas.TimingStdNs = max(meas.TimingStdNs, opts.TimingStdFloorNs);

referencePCI = chooseReferencePCILocal(meas, opts);
captureIDs = unique(meas.CaptureID, "stable");
frameNs = opts.FramePeriodMs * 1e6;

if ~isValidInput
    relativeMeasurements = emptyRelativeMeasurementsLocal();
    surveyInfo = baseSurveyInfoFromInputsLocal(referencePCI, captureIDs, strings(0,1), ...
        0, "Input receiver-position validation failed before forming relative measurements.");
    surveyInfo.InputValidationStatus = "NOT_ASSESSABLE";
    surveyInfo.InputValidationReason = validationReason;
    surveyInfo.BadCaptureIDs = badCaptureIDs;
    return;
end

rows = {};
rowCount = 0;
missingReference = strings(0,1);

for m = 1:numel(captureIDs)
    capID = captureIDs(m);
    capRows = meas(meas.CaptureID == capID, :);
    refRows = capRows(capRows.PCI == referencePCI, :);
    if isempty(refRows) || height(refRows) == 0
        missingReference(end+1,1) = capID; %#ok<AGROW>
        continue;
    end

    ref = bestMeasurementRowLocal(refRows);
    capPos = captureXY(captureXY.CaptureID == capID, :);
    if isempty(capPos)
        continue;
    end

    for k = 1:height(capRows)
        if capRows.PCI(k) == referencePCI
            continue;
        end

        relNs = wrapNsLocal(capRows.FramePhaseNs(k) - ref.FramePhaseNs, frameNs);
        sigmaNs = hypot(capRows.TimingStdNs(k), ref.TimingStdNs);

        rowCount = rowCount + 1;
        rows{rowCount,1} = table( ...
            capID, capRows.PCI(k), referencePCI, ...
            relNs, sigmaNs, ...
            capPos.RxX_m, capPos.RxY_m, capPos.RxZ_m, ...
            capPos.RxPositionUncertaintyM, ...
            capRows.FramePhaseNs(k), ref.FramePhaseNs, ...
            'VariableNames', ["CaptureID","PCI","ReferencePCI", ...
            "RelativeArrivalNs","SigmaNs","RxX_m","RxY_m","RxZ_m", ...
            "RxPositionUncertaintyM","FramePhaseNs","ReferenceFramePhaseNs"]);
    end
end

if isempty(rows)
    relativeMeasurements = emptyRelativeMeasurementsLocal();
else
    relativeMeasurements = vertcat(rows{:});
end

surveyInfo = struct();
surveyInfo.ReferencePCI = referencePCI;
surveyInfo.NumCaptures = numel(captureIDs);
surveyInfo.NumCapturesWithReference = numel(captureIDs) - numel(missingReference);
surveyInfo.MissingReferenceCaptureIDs = missingReference;
surveyInfo.NumRelativeMeasurements = height(relativeMeasurements);
surveyInfo.Statement = "Relative measurements subtract one detected PCI inside each capture to remove the B210 capture-specific clock phase.";
surveyInfo.InputValidationStatus = "OK";
surveyInfo.InputValidationReason = "";
surveyInfo.BadCaptureIDs = strings(0,1);
end

function validateSurveyInputsLocal(measurementTable, captureInfo)
requiredMeas = ["CaptureID","PCI","FramePhaseNs"];
missingMeas = requiredMeas(~ismember(requiredMeas, string(measurementTable.Properties.VariableNames)));
if ~isempty(missingMeas)
    error("estimateSurveyTimingNoLocations:MissingMeasurementColumns", ...
        "measurementTable is missing columns: %s", strjoin(missingMeas, ", "));
end

if ~ismember("CaptureID", string(captureInfo.Properties.VariableNames))
    error("estimateSurveyTimingNoLocations:MissingCaptureID", ...
        "captureInfo must contain CaptureID.");
end
end

function [isValid, reason, badCaptureIDs] = validateReceiverCoverageLocal(meas, captureXY)
isValid = true;
reason = "";
badCaptureIDs = strings(0,1);
captureIDs = unique(string(meas.CaptureID), "stable");
for k = 1:numel(captureIDs)
    capID = captureIDs(k);
    idx = find(captureXY.CaptureID == capID);
    if numel(idx) ~= 1
        isValid = false;
        badCaptureIDs(end+1,1) = capID; %#ok<AGROW>
        if isempty(idx)
            reason = "Missing receiver-position row for at least one measurement CaptureID.";
        else
            reason = "Duplicate receiver-position rows for at least one measurement CaptureID.";
        end
        return;
    end

    row = captureXY(idx, :);
    if ~isfinite(row.RxX_m) || ~isfinite(row.RxY_m) || ~isfinite(row.RxZ_m) || ...
            ~isfinite(row.RxPositionUncertaintyM)
        isValid = false;
        badCaptureIDs(end+1,1) = capID; %#ok<AGROW>
        reason = "Receiver-position row contains non-finite coordinates or uncertainty.";
        return;
    end
end
end

function referencePCI = chooseReferencePCILocal(meas, opts)
if isfield(opts, "ReferencePCI") && isfinite(opts.ReferencePCI)
    referencePCI = double(opts.ReferencePCI);
    return;
end

pcis = unique(meas.PCI, "stable");
numCaps = zeros(numel(pcis), 1);
meanMetric = zeros(numel(pcis), 1);
for k = 1:numel(pcis)
    rows = meas(meas.PCI == pcis(k), :);
    numCaps(k) = numel(unique(rows.CaptureID));
    if ismember("SNRdB", string(rows.Properties.VariableNames))
        meanMetric(k) = mean(rows.SNRdB, "omitnan");
        if ~isfinite(meanMetric(k))
            meanMetric(k) = -mean(rows.TimingStdNs, "omitnan");
        end
    else
        meanMetric(k) = -mean(rows.TimingStdNs, "omitnan");
    end
end

[~, order] = sortrows([-numCaps(:), -meanMetric(:)]);
referencePCI = double(pcis(order(1)));
end

function row = bestMeasurementRowLocal(rows)
if ismember("SNRdB", string(rows.Properties.VariableNames))
    [~, idx] = max(rows.SNRdB);
elseif ismember("TimingStdNs", string(rows.Properties.VariableNames))
    [~, idx] = min(rows.TimingStdNs);
else
    idx = 1;
end
row = rows(idx, :);
end

function pcis = orderedPCIsLocal(pcis, referencePCI)
pcis = pcis(:);
if isempty(pcis) || ~isfinite(referencePCI)
    return;
end
pcis = [referencePCI; pcis(pcis ~= referencePCI)];
end

function surveyInfo = emptySurveyInfoLocal(captureInfo, statement)
if istable(captureInfo) && ismember("CaptureID", string(captureInfo.Properties.VariableNames))
    numCaptures = height(captureInfo);
else
    numCaptures = 0;
end
surveyInfo = baseSurveyInfoFromInputsLocal(NaN, strings(0,1), strings(0,1), 0, statement);
surveyInfo.NumCaptures = numCaptures;
surveyInfo.InputValidationStatus = "OK";
surveyInfo.InputValidationReason = "";
surveyInfo.BadCaptureIDs = strings(0,1);
end

function surveyInfo = baseSurveyInfoFromInputsLocal(referencePCI, captureIDs, missingReference, numRelativeMeasurements, statement)
surveyInfo = struct();
surveyInfo.ReferencePCI = referencePCI;
surveyInfo.NumCaptures = numel(captureIDs);
surveyInfo.NumCapturesWithReference = numel(captureIDs) - numel(missingReference);
surveyInfo.MissingReferenceCaptureIDs = missingReference;
surveyInfo.NumRelativeMeasurements = numRelativeMeasurements;
surveyInfo.Statement = statement;
surveyInfo.InputValidationStatus = "OK";
surveyInfo.InputValidationReason = "";
surveyInfo.BadCaptureIDs = strings(0,1);
end

function captureXY = receiverXYLocal(captureInfo, opts)
cap = captureInfo;
cap.CaptureID = string(cap.CaptureID);
names = string(cap.Properties.VariableNames);

if all(ismember(["RxX_m","RxY_m"], names))
    x = cap.RxX_m;
    y = cap.RxY_m;
    if ismember("RxZ_m", names)
        z = cap.RxZ_m;
    else
        z = zeros(height(cap), 1);
    end
elseif all(ismember(["RxLat","RxLon"], names))
    lat0 = cap.RxLat(1);
    lon0 = cap.RxLon(1);
    earthRadiusM = 6371000;
    x = deg2rad(cap.RxLon - lon0) .* earthRadiusM .* cos(deg2rad(lat0));
    y = deg2rad(cap.RxLat - lat0) .* earthRadiusM;
    if ismember("RxAltM", names)
        z = cap.RxAltM - cap.RxAltM(1);
    else
        z = zeros(height(cap), 1);
    end
else
    error("estimateSurveyTimingNoLocations:MissingReceiverPosition", ...
        "captureInfo must contain either RxX_m/RxY_m or RxLat/RxLon.");
end

if ismember("RxPositionUncertaintyM", names)
    uncertaintyM = cap.RxPositionUncertaintyM;
else
    uncertaintyM = opts.ReceiverPositionUncertaintyM * ones(height(cap), 1);
end

captureXY = table(cap.CaptureID, x(:), y(:), z(:), uncertaintyM(:), ...
    'VariableNames', ["CaptureID","RxX_m","RxY_m","RxZ_m", ...
    "RxPositionUncertaintyM"]);
end

function tbl = emptyRelativeMeasurementsLocal()
tbl = table( ...
    strings(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', ["CaptureID","PCI","ReferencePCI","RelativeArrivalNs", ...
    "SigmaNs","RxX_m","RxY_m","RxZ_m","RxPositionUncertaintyM", ...
    "FramePhaseNs","ReferenceFramePhaseNs"]);
end

function opts = defaultSurveyOptions()
% Defaults for no-gNB-location survey fitting.
opts = struct();
opts.SpeedOfLightMps = 299792458;
opts.FramePeriodMs = 10;
opts.ReferencePCI = NaN;
opts.MinVisibleGNBs = 2;
opts.MinCellsForCleanGroupPass = 3;
opts.MinCaptures = 0;
opts.MinMeasurementsPerPCI = 0;
opts.MinDistinctReceiverPositionsPerTarget = 3;
opts.RequireReferenceInEveryCapture = true;
opts.TimingStdNsDefault = 30;
opts.TimingStdFloorNs = 10;
opts.ReceiverPositionUncertaintyM = 5;
opts.PositionSearchRadiusM = 1500;
opts.NumRandomStarts = 24;
opts.RandomSeed = 20260603;
opts.ConsensusRelativeCostTolerance = 0.05;
opts.ConsensusAbsoluteCostTolerance = 1e-6;
opts.UseLsqnonlinIfAvailable = true;
opts.MaxIterations = 2500;
opts.FunctionTolerance = 1e-10;
% These fractions and verdict thresholds are thesis engineering settings.
% They are configurable and are not regulatory or 3GPP limits.
opts.FitRMSPassFraction = 0.20;
opts.BiasUncertaintyPassFraction = 0.20;
opts.MaxFitRMSNsForPass = NaN;
opts.MaxBiasUncertaintyForPassNs = NaN;
opts.MaxConditionNumber = 1e10;
opts.SyncThresholdNs = NaN;
opts.StaticSyncThresholdNs = NaN;
opts.StaticOffsetLimitNs = NaN;
opts.StaticDesyncThresholdNs = NaN;
opts.WarningOffsetNs = NaN;
opts.FailOffsetNs = 1000;
opts.HardFailOffsetNs = Inf;
opts.UncertaintySigmaMultiplier = 3;
opts.MinInstabilitySamplesPerPCI = 4;
opts.JitterRmsThresholdNs = NaN;
opts.InstabilityRmsThresholdNs = NaN;
opts.TimingInstabilityRmsLimitNs = NaN;
opts.JitterPeakToPeakThresholdNs = NaN;
opts.InstabilityPeakToPeakThresholdNs = NaN;
opts.TimingInstabilityPeakToPeakLimitNs = NaN;
opts.WarningInstabilityRmsNs = NaN;
opts.FailInstabilityRmsNs = NaN;
opts.WarningInstabilityPeakToPeakNs = NaN;
opts.FailInstabilityPeakToPeakNs = NaN;
opts.HardFailInstabilityPeakToPeakNs = Inf;
opts.RequireInstabilityAssessableForPass = true;
end

function opts = normalizeDecisionOptionsLocal(opts, userOpts)
% User-facing aliases:
% - SyncThresholdNs / Static*Limit*: static relative timing limit.
% - Jitter*/Instability*Threshold*: relative instability limits.
% When a single synchronization limit is supplied, the warning band is not a
% fixed percentage. It is the uncertainty-overlap region around the limit.
[hasStaticLimit, staticLimit] = firstFiniteOptionLocal(opts, ...
    ["SyncThresholdNs","StaticSyncThresholdNs","StaticOffsetLimitNs","StaticDesyncThresholdNs"]);
if hasStaticLimit
    opts.FailOffsetNs = staticLimit;
    opts.SyncThresholdNs = staticLimit;
    if ~optionProvidedLocal(userOpts, "WarningOffsetNs")
        opts.WarningOffsetNs = NaN;
    end
    if ~optionProvidedLocal(userOpts, "HardFailOffsetNs")
        opts.HardFailOffsetNs = Inf;
    end
end

if ~isfinite(opts.MaxFitRMSNsForPass)
    opts.MaxFitRMSNsForPass = opts.FitRMSPassFraction * opts.FailOffsetNs;
end
if ~isfinite(opts.MaxBiasUncertaintyForPassNs)
    opts.MaxBiasUncertaintyForPassNs = opts.BiasUncertaintyPassFraction * opts.FailOffsetNs;
end

[hasRmsLimit, rmsLimit] = firstFiniteOptionLocal(opts, ...
    ["JitterRmsThresholdNs","InstabilityRmsThresholdNs","TimingInstabilityRmsLimitNs"]);
if hasRmsLimit
    opts.FailInstabilityRmsNs = rmsLimit;
    if ~optionProvidedLocal(userOpts, "WarningInstabilityRmsNs")
        opts.WarningInstabilityRmsNs = NaN;
    end
end

[hasPpLimit, ppLimit] = firstFiniteOptionLocal(opts, ...
    ["JitterPeakToPeakThresholdNs","InstabilityPeakToPeakThresholdNs", ...
     "TimingInstabilityPeakToPeakLimitNs"]);
if hasPpLimit
    opts.FailInstabilityPeakToPeakNs = ppLimit;
    if ~optionProvidedLocal(userOpts, "WarningInstabilityPeakToPeakNs")
        opts.WarningInstabilityPeakToPeakNs = NaN;
    end
    if ~optionProvidedLocal(userOpts, "HardFailInstabilityPeakToPeakNs")
        opts.HardFailInstabilityPeakToPeakNs = Inf;
    end
end
end

function [hasValue, value] = firstFiniteOptionLocal(opts, names)
hasValue = false;
value = NaN;
for k = 1:numel(names)
    name = char(names(k));
    if isfield(opts, name) && isnumeric(opts.(name)) && isscalar(opts.(name)) && isfinite(opts.(name))
        hasValue = true;
        value = double(opts.(name));
        return;
    end
end
end

function tf = optionProvidedLocal(opts, name)
tf = isstruct(opts) && isfield(opts, char(name));
end

function tf = thresholdsAreSatisfiedLocal(valueWithMargin, threshold)
tf = ~isfinite(threshold) || valueWithMargin < threshold;
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
