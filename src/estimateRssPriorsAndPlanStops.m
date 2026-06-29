function varargout = estimateRssPriorsAndPlanStops(mode, varargin)
%estimateRssPriorsAndPlanStops RSS-assisted gNB priors and survey-stop planning.
%
%   This Phase 13 extension is deliberately advisory. It uses SS-RSRP-like
%   power measurements to estimate rough gNB position priors and rank future
%   B210 receiver stops. It does not feed RSS-derived positions into strict
%   synchronization verdicts by default.
%
%   Public modes:
%     "power"  -> powerTable = estimateRssPriorsAndPlanStops("power", iq, meta, captureResult, opts)
%     "priors" -> rssPriors = estimateRssPriorsAndPlanStops("priors", powerTable, captureInfo, opts)
%     "plan"   -> stopTable = estimateRssPriorsAndPlanStops("plan", powerTable, captureInfo, rssPriors, candidateStops, opts)
%     "all"    -> result = estimateRssPriorsAndPlanStops("all", powerTable, captureInfo, candidateStops, opts)

arguments
    mode (1,1) string
end
arguments (Repeating)
    varargin
end

switch lower(mode)
    case "power"
        varargout{1} = estimateSSBPowerMeasurements(varargin{:});
    case "priors"
        varargout{1} = estimateRssGnbPriors(varargin{:});
    case "plan"
        varargout{1} = planAdaptiveReceiverStops(varargin{:});
    case "all"
        result = runAllLocal(varargin{:});
        varargout{1} = result;
    otherwise
        error("estimateRssPriorsAndPlanStops:UnknownMode", ...
            "Unknown mode '%s'. Use power, priors, plan, or all.", mode);
end

end

function result = runAllLocal(powerTable, captureInfo, candidateStops, opts)
if nargin < 4
    opts = struct();
end
if getOptLocal(opts, "UseForTimingVerdict", false)
    error("estimateRssPriorsAndPlanStops:TimingUseUnsupported", ...
        "UseForTimingVerdict=true is not supported. RSS priors are planning aids only and are not validated as timing-grade propagation corrections.");
end
rssPriors = estimateRssGnbPriors(powerTable, captureInfo, opts);
suggestedStops = planAdaptiveReceiverStops(powerTable, captureInfo, rssPriors, candidateStops, opts);

result = struct();
result.Mode = "rss_assisted_survey_planning";
result.Statement = "RSS/SS-RSRP-like power is used only for rough gNB priors and receiver-stop planning; timing verdicts still come from the timing survey.";
result.PowerTable = powerTable;
result.RssPriors = rssPriors;
result.SuggestedReceiverStops = suggestedStops;
result.UseForTimingVerdict = false;
end

function powerTable = estimateSSBPowerMeasurements(iq, meta, captureResult, opts)
% Estimate one relative SS-RSRP-like power row per detected PCI.
if nargin < 4
    opts = struct();
end
opts = mergeStructsLocal(defaultRssOptionsLocal(), opts);
[meta, ~] = validateMetadata(meta);

if ~isfield(captureResult, "SSSDetections") || isempty(captureResult.SSSDetections) || ...
        height(captureResult.SSSDetections) == 0
    powerTable = emptyPowerTableLocal();
    return;
end

det = captureResult.SSSDetections;
if ismember("IsUsable", string(det.Properties.VariableNames))
    det = det(det.IsUsable, :);
end
if isempty(det) || height(det) == 0
    powerTable = emptyPowerTableLocal();
    return;
end

refs = buildPowerReferenceLocal(meta, opts);
rows = cell(height(det), 1);
rowCount = 0;
for k = 1:height(det)
    [powerDb, noiseDb, snrDb, quality] = powerFromDetectionLocal(iq, det(k,:), meta, refs, opts);
    if ~isfinite(powerDb)
        continue;
    end
    rowCount = rowCount + 1;
    rows{rowCount} = table( ...
        string(getOptLocal(opts, "CaptureID", "")), ...
        det.PCI(k), ...
        getOptLocal(opts, "RxX_m", NaN), ...
        getOptLocal(opts, "RxY_m", NaN), ...
        getOptLocal(opts, "RxZ_m", NaN), ...
        powerDb, noiseDb, snrDb, 1, quality, ...
        'VariableNames', ["CaptureID","PCI","RxX_m","RxY_m","RxZ_m", ...
        "PowerDbFS","NoiseDbFS","SNRdB","NumDetections","QualityFlag"]);
end

if rowCount == 0
    powerTable = emptyPowerTableLocal();
    return;
end

detPower = vertcat(rows{1:rowCount});
pcis = unique(detPower.PCI, "stable");
groupRows = cell(numel(pcis), 1);
for p = 1:numel(pcis)
    pciRows = detPower(detPower.PCI == pcis(p), :);
    powerLin = mean(db2powLocal(pciRows.PowerDbFS), "omitnan");
    noiseLin = mean(db2powLocal(pciRows.NoiseDbFS), "omitnan");
    powerDb = pow2dbLocal(powerLin);
    noiseDb = pow2dbLocal(noiseLin);
    snrDb = powerDb - noiseDb;
    if height(pciRows) < opts.MinPowerDetectionsPerPCI || snrDb < opts.MinPowerSNRdB
        quality = "SUSPECT";
    else
        quality = "OK";
    end

    groupRows{p} = table( ...
        pciRows.CaptureID(1), pcis(p), ...
        pciRows.RxX_m(1), pciRows.RxY_m(1), pciRows.RxZ_m(1), ...
        powerDb, noiseDb, snrDb, height(pciRows), quality, ...
        'VariableNames', ["CaptureID","PCI","RxX_m","RxY_m","RxZ_m", ...
        "PowerDbFS","NoiseDbFS","SNRdB","NumDetections","QualityFlag"]);
end
powerTable = vertcat(groupRows{:});
end

function rssPriors = estimateRssGnbPriors(powerTable, captureInfo, opts)
% Estimate rough gNB positions from relative power over known receiver stops.
if nargin < 3
    opts = struct();
end
opts = mergeStructsLocal(defaultRssOptionsLocal(), opts);

if isempty(powerTable) || height(powerTable) == 0
    rssPriors = emptyPriorTableLocal();
    return;
end

powerTable = attachReceiverXYLocal(powerTable, captureInfo, opts);
pcis = unique(powerTable.PCI, "stable");
candidateXY = makeSearchGridLocal(powerTable, opts);
rows = cell(numel(pcis), 1);

for k = 1:numel(pcis)
    rows{k} = estimateOnePriorLocal(powerTable(powerTable.PCI == pcis(k), :), ...
        pcis(k), candidateXY, opts);
end

rssPriors = vertcat(rows{:});
end

function stopTable = planAdaptiveReceiverStops(powerTable, captureInfo, rssPriors, candidateStops, opts)
% Rank candidate B210 positions using RSS priors and timing-geometry health.
if nargin < 5
    opts = struct();
end
opts = mergeStructsLocal(defaultRssOptionsLocal(), opts);

if isempty(rssPriors) || height(rssPriors) == 0
    stopTable = emptyStopTableLocal();
    return;
end

powerTable = attachReceiverXYLocal(powerTable, captureInfo, opts);
captureXY = unique(powerTable(:, ["RxX_m","RxY_m"]), "rows", "stable");
candidateStops = canonicalCandidateStopsLocal(candidateStops, captureXY, opts);

usablePriors = rssPriors(rssPriors.PriorStatus ~= "NOT_ASSESSABLE", :);
if isempty(usablePriors) || height(usablePriors) == 0 || isempty(candidateStops)
    stopTable = emptyStopTableLocal();
    return;
end

existingXY = [captureXY.RxX_m, captureXY.RxY_m];
rows = cell(height(candidateStops), 1);
for k = 1:height(candidateStops)
    cand = [candidateStops.X_m(k), candidateStops.Y_m(k)];
    [score, visible, minSnr, condAfter, reason] = scoreCandidateLocal(cand, ...
        existingXY, usablePriors, opts);
    rows{k} = table( ...
        candidateStops.CandidateID(k), cand(1), cand(2), score, visible, ...
        minSnr, condAfter, reason, ...
        'VariableNames', ["CandidateID","X_m","Y_m","Score", ...
        "ExpectedVisiblePCIs","PredictedMinSNRdB","ExpectedConditionAfter","Reason"]);
end

stopTable = vertcat(rows{:});
stopTable = sortrows(stopTable, "Score", "descend");
if height(stopTable) > opts.MaxSuggestedStops
    stopTable = stopTable(1:opts.MaxSuggestedStops, :);
end
end

function [powerDb, noiseDb, snrDb, quality] = powerFromDetectionLocal(iq, row, meta, refs, opts)
powerDb = NaN;
noiseDb = NaN;
snrDb = NaN;
quality = "NOT_ASSESSABLE";

startIdx = round(row.StartSample0) + 1;
endIdx = startIdx + refs.ReferenceLength - 1;
if startIdx < 1 || endIdx > numel(iq)
    return;
end

segment = iq(startIdx:endIdx);
if ismember("CFOHz", string(row.Properties.VariableNames)) && isfinite(row.CFOHz)
    n = (0:numel(segment)-1).';
    segment = segment .* exp(-1i*2*pi*row.CFOHz*n/meta.SampleRateHz);
end

try
    demodOpts = struct( ...
        "SCSkHz", opts.SCSkHz, ...
        "NSizeGrid", opts.NSizeGrid, ...
        "SSBStartSymbol", opts.SSBStartSymbol, ...
        "NCellID", row.PCI);
    rxGrid = demodulateSSBWaveform(segment, meta, demodOpts);
    measIdx = nrSSSIndices;
    if opts.UsePBCHDMRS
        dmrsIdx = nrPBCHDMRSIndices(row.PCI);
        dmrsIdx = dmrsIdx(dmrsIdx <= numel(rxGrid));
        measIdx = unique([measIdx(:); dmrsIdx(:)]);
    end
    measIdx = measIdx(measIdx <= numel(rxGrid));
    if numel(measIdx) < 32
        error("estimateRssPriorsAndPlanStops:TooFewRE", "Too few measurement REs.");
    end

    rePower = abs(rxGrid(measIdx)).^2;
    allPower = abs(rxGrid(:)).^2;
    mask = true(numel(allPower), 1);
    mask(measIdx) = false;
    noiseSamples = allPower(mask);
    noiseSamples = noiseSamples(isfinite(noiseSamples) & noiseSamples > 0);
    if isempty(noiseSamples)
        noisePower = percentileLocal(allPower, opts.NoisePercentile);
    else
        noisePower = percentileLocal(noiseSamples, opts.NoisePercentile);
    end
    signalPower = mean(rePower, "omitnan");
    quality = "SSS_RE";
catch
    [signalPower, noisePower] = fallbackWindowPowerLocal(iq, startIdx, endIdx, opts);
    quality = "WINDOW_FALLBACK";
end

powerDb = pow2dbLocal(signalPower);
noiseDb = pow2dbLocal(noisePower);
snrDb = powerDb - noiseDb;
if snrDb < opts.MinPowerSNRdB
    quality = "SUSPECT";
end
end

function [signalPower, noisePower] = fallbackWindowPowerLocal(iq, startIdx, endIdx, opts)
sig = iq(startIdx:endIdx);
signalPower = mean(abs(sig).^2, "omitnan");
guard = endIdx - startIdx + 1;
left = max(1, startIdx - 2*guard):max(1, startIdx - guard);
right = min(numel(iq), endIdx + guard):min(numel(iq), endIdx + 2*guard);
noiseSamples = [iq(left); iq(right)];
if isempty(noiseSamples)
    noisePower = percentileLocal(abs(iq).^2, opts.NoisePercentile);
else
    noisePower = percentileLocal(abs(noiseSamples).^2, opts.NoisePercentile);
end
end

function refs = buildPowerReferenceLocal(meta, opts)
refOpts = struct( ...
    "SCSkHz", opts.SCSkHz, ...
    "NSizeGrid", opts.NSizeGrid, ...
    "SSBStartSymbol", opts.SSBStartSymbol, ...
    "IncludePSS", true, ...
    "IncludeSSS", false, ...
    "IncludePBCHDMRS", false);
waveform = buildSSBSyncWaveform(0, meta, refOpts);
refs.ReferenceLength = numel(waveform);
end

function row = estimateOnePriorLocal(rowsIn, pci, candidateXY, opts)
rowsIn = rowsIn(isfinite(rowsIn.PowerDbFS) & isfinite(rowsIn.RxX_m) & isfinite(rowsIn.RxY_m), :);
numMeas = height(rowsIn);
if numMeas < opts.MinMeasurementsPerPCI
    row = priorRowLocal(pci, NaN, NaN, NaN, NaN, NaN, NaN, numMeas, NaN, ...
        "NOT_ASSESSABLE", "Too few RSS measurements for this PCI.");
    return;
end

rxXY = [rowsIn.RxX_m, rowsIn.RxY_m];
powerDb = rowsIn.PowerDbFS(:);
[isPoorGeom, geomReason] = poorGeometryLocal(rxXY, opts);

best = struct("cost", Inf, "xy", [NaN NaN], "ple", NaN, "a", NaN, ...
    "residual", [], "sigma", Inf);
for n = opts.PathLossExponentGrid(:).'
    logD = log10(max(distanceMatrixLocal(candidateXY, rxXY), opts.MinDistanceM));
    basis = -10 * n * logD;
    for c = 1:size(candidateXY, 1)
        predShape = basis(c,:).';
        aHat = robustLocationLocal(powerDb - predShape, opts);
        residual = powerDb - (aHat + predShape);
        cost = robustCostLocal(residual, opts.ShadowingSigmaDb);
        if cost < best.cost
            best.cost = cost;
            best.xy = candidateXY(c,:);
            best.ple = n;
            best.a = aHat;
            best.residual = residual;
            best.sigma = robustSigmaLocal(residual);
        end
    end
end

likelihood = priorLikelihoodLocal(candidateXY, rxXY, powerDb, best.ple, best.cost, opts);
[r50, r90] = radiusFromLikelihoodLocal(candidateXY, best.xy, likelihood);
status = "ROUGH_ONLY";
reason = "RSS prior only; do not use as timing-grade propagation correction.";
if isPoorGeom
    status = "SUSPECT";
    reason = "RSS geometry is weak: " + geomReason;
elseif best.sigma > opts.MaxResidualDbForRough
    status = "SUSPECT";
    reason = "RSS residual is high, likely due to NLOS, beamforming, or shadowing.";
elseif max(abs(best.residual)) > opts.MaxOutlierResidualDb
    status = "SUSPECT";
    reason = "At least one RSS measurement is an outlier, likely due to NLOS or beam mismatch.";
elseif r90 > opts.MaxRssRadius90ForRough
    status = "SUSPECT";
    reason = "RSS likelihood is broad, so the position prior is too uncertain for confident planning.";
end

row = priorRowLocal(pci, best.xy(1), best.xy(2), r50, r90, best.ple, ...
    best.sigma, numMeas, best.a, status, reason);
end

function likelihood = priorLikelihoodLocal(candidateXY, rxXY, powerDb, ple, bestCost, opts)
numC = size(candidateXY, 1);
costs = zeros(numC, 1);
logD = log10(max(distanceMatrixLocal(candidateXY, rxXY), opts.MinDistanceM));
for c = 1:numC
    predShape = -10 * ple * logD(c,:).';
    aHat = robustLocationLocal(powerDb - predShape, opts);
    residual = powerDb - (aHat + predShape);
    costs(c) = robustCostLocal(residual, opts.ShadowingSigmaDb);
end
delta = max(0, costs - bestCost);
likelihood = exp(-0.5 * delta);
likelihood = likelihood ./ (sum(likelihood) + eps);
end

function candidateXY = makeSearchGridLocal(powerTable, opts)
rxXY = unique([powerTable.RxX_m, powerTable.RxY_m], "rows", "stable");
minXY = min(rxXY, [], 1) - opts.SearchMarginM;
maxXY = max(rxXY, [], 1) + opts.SearchMarginM;
x = minXY(1):opts.GridSpacingM:maxXY(1);
y = minXY(2):opts.GridSpacingM:maxXY(2);
[xx, yy] = meshgrid(x, y);
candidateXY = [xx(:), yy(:)];
end

function tbl = attachReceiverXYLocal(tbl, captureInfo, opts)
names = string(tbl.Properties.VariableNames);
hasXY = all(ismember(["RxX_m","RxY_m"], names));
if hasXY && all(isfinite(tbl.RxX_m)) && all(isfinite(tbl.RxY_m))
    return;
end

capXY = receiverXYLocal(captureInfo, opts);
if ~ismember("CaptureID", names)
    error("estimateRssPriorsAndPlanStops:MissingCaptureID", ...
        "powerTable must include CaptureID when RxX_m/RxY_m are missing or incomplete.");
end

tbl.CaptureID = string(tbl.CaptureID);
for k = 1:height(tbl)
    idx = find(capXY.CaptureID == tbl.CaptureID(k), 1);
    if isempty(idx)
        continue;
    end
    tbl.RxX_m(k) = capXY.RxX_m(idx);
    tbl.RxY_m(k) = capXY.RxY_m(idx);
    if ismember("RxZ_m", string(capXY.Properties.VariableNames))
        tbl.RxZ_m(k) = capXY.RxZ_m(idx);
    end
end
end

function candidateStops = canonicalCandidateStopsLocal(candidateStops, captureXY, opts)
if istable(candidateStops) && height(candidateStops) > 0
    names = string(candidateStops.Properties.VariableNames);
    if ~all(ismember(["X_m","Y_m"], names))
        error("estimateRssPriorsAndPlanStops:BadCandidateStops", ...
            "candidateStops must include X_m and Y_m columns.");
    end
    if ~ismember("CandidateID", names)
        candidateStops.CandidateID = "cand_" + string((1:height(candidateStops)).');
    else
        candidateStops.CandidateID = string(candidateStops.CandidateID);
    end
    return;
end

rxXY = [captureXY.RxX_m, captureXY.RxY_m];
center = mean(rxXY, 1);
spread = max(vecnorm(rxXY - center, 2, 2));
radius = max(opts.AutoCandidateRadiusM, spread + opts.AutoCandidateMarginM);
angles = linspace(0, 2*pi, opts.NumAutoCandidates + 1).';
angles(end) = [];
xy = center + radius * [cos(angles), sin(angles)];
candidateStops = table("auto_" + string((1:size(xy,1)).'), xy(:,1), xy(:,2), ...
    'VariableNames', ["CandidateID","X_m","Y_m"]);
end

function [score, visible, minSnr, condAfter, reason] = scoreCandidateLocal(cand, existingXY, priors, opts)
distExisting = vecnorm(existingXY - cand, 2, 2);
if any(distExisting < opts.MinCandidateSpacingM)
    score = -Inf;
    visible = 0;
    minSnr = -Inf;
    condAfter = Inf;
    reason = "Rejected: too close to an existing receiver stop.";
    return;
end

numP = height(priors);
predPower = NaN(numP, 1);
for k = 1:numP
    d = max(norm(cand - [priors.EstimatedX_m(k), priors.EstimatedY_m(k)]), opts.MinDistanceM);
    predPower(k) = priors.ReferencePowerDb(k) - 10 * priors.PathLossExponent(k) * log10(d);
end
visibleMask = predPower >= opts.VisibilityPowerDbFS;
visible = nnz(visibleMask);
if visible > 0
    minSnr = min(predPower(visibleMask) - opts.NoiseFloorDbFS);
else
    minSnr = -Inf;
end

newXY = [existingXY; cand];
condAfter = geometryConditionLocal(newXY);
spreadScore = minDistanceToRouteScoreLocal(cand, existingXY, opts);
rssInfo = rssInformationScoreLocal(cand, priors, opts);
visibilityScore = visible / max(1, numP);
conditionScore = 1 / max(log10(condAfter + 10), 1);
score = 0.35*rssInfo + 0.25*visibilityScore + 0.25*conditionScore + 0.15*spreadScore;

if visible < opts.MinExpectedVisiblePCIs
    score = score - 0.5;
    reason = "Low expected common-PCI visibility; useful only if access constraints force this stop.";
else
    reason = "Good candidate: improves RSS geometry while preserving expected PCI visibility.";
end
end

function info = rssInformationScoreLocal(cand, priors, opts)
values = zeros(height(priors), 1);
for k = 1:height(priors)
    p = [priors.EstimatedX_m(k), priors.EstimatedY_m(k)];
    r = max(norm(cand - p), opts.MinDistanceM);
    values(k) = 1 / (r / opts.AutoCandidateRadiusM + 0.1)^2;
end
info = min(1, mean(values, "omitnan"));
end

function s = minDistanceToRouteScoreLocal(cand, existingXY, opts)
d = min(vecnorm(existingXY - cand, 2, 2));
s = min(1, d / max(opts.AutoCandidateRadiusM, opts.MinCandidateSpacingM));
end

function condVal = geometryConditionLocal(xy)
if size(xy,1) < 3
    condVal = Inf;
    return;
end
centered = xy - mean(xy, 1);
covXY = centered.' * centered / max(1, size(xy,1)-1);
condVal = cond(covXY + 1e-6*eye(2));
end

function [tf, reason] = poorGeometryLocal(rxXY, opts)
tf = false;
reason = "";
if size(rxXY,1) < opts.MinMeasurementsPerPCI
    tf = true;
    reason = "too few receiver stops";
    return;
end
condVal = geometryConditionLocal(rxXY);
if condVal > opts.MaxRssGeometryCondition || ~isfinite(condVal)
    tf = true;
    reason = "receiver stops are close to collinear";
end
end

function d = distanceMatrixLocal(candidateXY, rxXY)
numC = size(candidateXY, 1);
numR = size(rxXY, 1);
d = zeros(numC, numR);
for c = 1:numC
    d(c,:) = vecnorm(rxXY - candidateXY(c,:), 2, 2).';
end
end

function [r50, r90] = radiusFromLikelihoodLocal(candidateXY, bestXY, likelihood)
dist = vecnorm(candidateXY - bestXY, 2, 2);
[distSorted, order] = sort(dist);
cdf = cumsum(likelihood(order));
r50 = percentileRadiusLocal(distSorted, cdf, 0.50);
r90 = percentileRadiusLocal(distSorted, cdf, 0.90);
end

function r = percentileRadiusLocal(distSorted, cdf, p)
idx = find(cdf >= p, 1);
if isempty(idx)
    r = max(distSorted);
else
    r = distSorted(idx);
end
end

function y = robustLocationLocal(x, opts)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
    return;
end
y = median(x);
for iter = 1:opts.RobustIterations
    r = x - y;
    scale = max(robustSigmaLocal(r), 1);
    w = min(1, opts.HuberK * scale ./ (abs(r) + eps));
    y = sum(w .* x) / sum(w);
end
end

function sigma = robustSigmaLocal(r)
r = r(isfinite(r));
if isempty(r)
    sigma = Inf;
else
    sigma = 1.4826 * median(abs(r - median(r))) + eps;
end
end

function cost = robustCostLocal(r, sigmaDb)
r = r(isfinite(r));
if isempty(r)
    cost = Inf;
    return;
end
u = r / max(sigmaDb, eps);
cost = sum(min(u.^2, 4));
end

function capXY = receiverXYLocal(captureInfo, opts)
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
    error("estimateRssPriorsAndPlanStops:MissingReceiverPosition", ...
        "captureInfo must contain either RxX_m/RxY_m or RxLat/RxLon.");
end

capXY = table(cap.CaptureID, x(:), y(:), z(:), ...
    'VariableNames', ["CaptureID","RxX_m","RxY_m","RxZ_m"]);
end

function row = priorRowLocal(pci, x, y, r50, r90, ple, residualDb, numMeas, referencePowerDb, status, reason)
row = table(pci, x, y, r50, r90, ple, residualDb, numMeas, ...
    referencePowerDb, status, reason, ...
    'VariableNames', ["PCI","EstimatedX_m","EstimatedY_m","Radius50M", ...
    "Radius90M","PathLossExponent","ResidualDb","NumMeasurements", ...
    "ReferencePowerDb","PriorStatus","Reason"]);
end

function tbl = emptyPowerTableLocal()
tbl = table(strings(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), strings(0,1), ...
    'VariableNames', ["CaptureID","PCI","RxX_m","RxY_m","RxZ_m", ...
    "PowerDbFS","NoiseDbFS","SNRdB","NumDetections","QualityFlag"]);
end

function tbl = emptyPriorTableLocal()
tbl = table(zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), strings(0,1), strings(0,1), ...
    'VariableNames', ["PCI","EstimatedX_m","EstimatedY_m","Radius50M", ...
    "Radius90M","PathLossExponent","ResidualDb","NumMeasurements", ...
    "ReferencePowerDb","PriorStatus","Reason"]);
end

function tbl = emptyStopTableLocal()
tbl = table(strings(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), strings(0,1), ...
    'VariableNames', ["CandidateID","X_m","Y_m","Score", ...
    "ExpectedVisiblePCIs","PredictedMinSNRdB","ExpectedConditionAfter","Reason"]);
end

function opts = defaultRssOptionsLocal()
opts = struct();
opts.SCSkHz = 30;
opts.NSizeGrid = 20;
opts.SSBStartSymbol = 2;
opts.UsePBCHDMRS = true;
opts.NoisePercentile = 20;
opts.MinPowerSNRdB = 3;
opts.MinPowerDetectionsPerPCI = 1;
opts.MinMeasurementsPerPCI = 4;
opts.PathLossExponentGrid = 1.6:0.1:4.5;
opts.ShadowingSigmaDb = 8;
opts.GridSpacingM = 50;
opts.SearchMarginM = 1500;
opts.MinCandidateSpacingM = 100;
opts.MinExpectedVisiblePCIs = 2;
opts.UseForTimingVerdict = false;
opts.MinDistanceM = 10;
opts.MaxResidualDbForRough = 12;
opts.MaxOutlierResidualDb = 20;
opts.MaxRssRadius90ForRough = 800;
opts.MaxRssGeometryCondition = 50;
opts.RobustIterations = 4;
opts.HuberK = 1.5;
opts.VisibilityPowerDbFS = -55;
opts.NoiseFloorDbFS = -75;
opts.AutoCandidateRadiusM = 1200;
opts.AutoCandidateMarginM = 500;
opts.NumAutoCandidates = 24;
opts.MaxSuggestedStops = 10;
end

function value = getOptLocal(opts, name, defaultValue)
if isfield(opts, name)
    value = opts.(name);
else
    value = defaultValue;
end
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

function y = db2powLocal(x)
y = 10.^(x/10);
end

function y = pow2dbLocal(x)
y = 10*log10(x + eps);
end

function y = percentileLocal(x, p)
x = sort(x(isfinite(x)));
if isempty(x)
    y = NaN;
    return;
end
p = min(max(p, 0), 100);
idx = 1 + (numel(x)-1) * p/100;
lo = floor(idx);
hi = ceil(idx);
if lo == hi
    y = x(lo);
else
    y = x(lo) + (x(hi) - x(lo)) * (idx - lo);
end
end
