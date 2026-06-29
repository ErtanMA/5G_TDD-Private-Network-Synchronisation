% test_phase13_rssPlanning
% Phase 13 tests for RSS-assisted rough gNB priors and receiver-stop planning.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot, "src")));

fprintf("Running Phase 13 RSS-assisted planning tests...\n");

rng(20260603, "twister");

[powerTable, captureInfo, truth] = makeRssFixture(0, false, false);
opts = phase13Opts();
priors = estimateRssPriorsAndPlanStops("priors", powerTable, captureInfo, opts);

assert(height(priors) == numel(truth.PCIs));
assert(all(priors.PriorStatus ~= "NOT_ASSESSABLE"));
assert(maxLocationError(priors, truth) < 300, ...
    "Clean RSS fixture should produce rough locations within a few hundred meters.");

% Unknown EIRP is handled by fitting one effective power constant per PCI.
[eirpTable, eirpInfo, eirpTruth] = makeRssFixture(0, true, false);
eirpPriors = estimateRssPriorsAndPlanStops("priors", eirpTable, eirpInfo, opts);
assert(maxLocationError(eirpPriors, eirpTruth) < 350, ...
    "Unknown per-gNB power should not break rough localization.");

% Searching PLE should fit at least as well as forcing a fixed free-space-like PLE.
fixedOpts = opts;
fixedOpts.PathLossExponentGrid = 2.0;
fullResidual = mean(priors.ResidualDb, "omitnan");
fixedPriors = estimateRssPriorsAndPlanStops("priors", powerTable, captureInfo, fixedOpts);
fixedResidual = mean(fixedPriors.ResidualDb, "omitnan");
assert(fullResidual <= fixedResidual + 1e-6);

% Shadowing stress: the estimator should degrade gracefully and keep uncertainty explicit.
for sigmaDb = [4 8 12]
    [shadowTable, shadowInfo, ~] = makeRssFixture(sigmaDb, true, false);
    shadowPriors = estimateRssPriorsAndPlanStops("priors", shadowTable, shadowInfo, opts);
    assert(height(shadowPriors) == numel(truth.PCIs));
    assert(all(ismember(shadowPriors.PriorStatus, ["ROUGH_ONLY","SUSPECT"])));
end

% One bad NLOS-like row should not collapse every prior.
[outlierTable, outlierInfo, outlierTruth] = makeRssFixture(4, true, true);
outlierPriors = estimateRssPriorsAndPlanStops("priors", outlierTable, outlierInfo, opts);
assert(all(isfinite(outlierPriors.EstimatedX_m)));
assert(all(outlierPriors.PriorStatus ~= "NOT_ASSESSABLE"));
assert(any(outlierPriors.PriorStatus == "SUSPECT") || ...
    maxLocationError(outlierPriors, outlierTruth) < 600, ...
    "NLOS outlier should either remain roughly correct or lower confidence.");

% Poor receiver geometry is flagged.
[lineTable, lineInfo, ~] = makeCollinearFixture();
linePriors = estimateRssPriorsAndPlanStops("priors", lineTable, lineInfo, opts);
assert(all(linePriors.PriorStatus == "SUSPECT"));

% Too few RSS measurements for one PCI is not assessable.
fewRows = powerTable(powerTable.PCI == truth.PCIs(1), :);
fewRows = fewRows(1:3, :);
fewPrior = estimateRssPriorsAndPlanStops("priors", fewRows, captureInfo, opts);
assert(fewPrior.PriorStatus == "NOT_ASSESSABLE");
assert(fewPrior.NumMeasurements == 3);

% Direct use should work when receiver coordinates are already in the power
% table, even if CaptureID is not present.
directRows = powerTable(powerTable.PCI == truth.PCIs(1), ...
    ["PCI","RxX_m","RxY_m","PowerDbFS"]);
directPrior = estimateRssPriorsAndPlanStops("priors", directRows, table(), opts);
assert(height(directPrior) == 1);
assert(directPrior.PriorStatus ~= "NOT_ASSESSABLE");

% Planner sanity: top candidate should be usable, visible, and not too close.
candidateStops = table( ...
    ["near"; "north"; "east"; "southwest"], ...
    [truth.RxXY(1,1)+10; 0; 1700; -1500], ...
    [truth.RxXY(1,2)+10; 1700; 0; -1300], ...
    'VariableNames', ["CandidateID","X_m","Y_m"]);
stops = estimateRssPriorsAndPlanStops("plan", powerTable, captureInfo, priors, candidateStops, opts);
assert(height(stops) >= 1);
assert(stops.Score(1) > -Inf);
assert(stops.ExpectedVisiblePCIs(1) >= opts.MinExpectedVisiblePCIs);
assert(stops.CandidateID(1) ~= "near", ...
    "Candidate too close to an existing stop should not be ranked first.");

allResult = estimateRssPriorsAndPlanStops("all", powerTable, captureInfo, candidateStops, opts);
assert(isfield(allResult, "RssPriors"));
assert(isfield(allResult, "SuggestedReceiverStops"));
assert(~allResult.UseForTimingVerdict);

badOpts = opts;
badOpts.UseForTimingVerdict = true;
didError = false;
try
    estimateRssPriorsAndPlanStops("all", powerTable, captureInfo, candidateStops, badOpts);
catch ME
    didError = ME.identifier == "estimateRssPriorsAndPlanStops:TimingUseUnsupported";
end
assert(didError, "UseForTimingVerdict=true must be rejected until validated.");

% Integration test: GPS manifests must be converted to nonzero local XY
% coordinates before RSS planning.
tmpRoot = tempname;
mkdir(tmpRoot);
cleanup = onCleanup(@() rmdir(tmpRoot, "s"));
gpsXY = [0 0;180 20;40 190;-160 80];
lat0 = 52.245567;
lon0 = 6.853050;
earthRadiusM = 6371000;
manifestParts = cell(size(gpsXY,1),1);
for m = 1:size(gpsXY,1)
    custom = struct();
    custom.NumGNBs = 3;
    custom.PCIs = [11 104 257];
    custom.CFOHz = [0 0 0];
    custom.FrameOffsetsNs = [0 0 0];
    custom.GNBPowerdB = [0 -2 -4] + [0.4 -0.3 0.2]*m;
    custom.LocationUncertaintyM = [30 30 30];
    custom.SiteDistanceM = [100 100 100];
    custom.SNRdB = 38;
    custom.DurationMs = 20;
    custom.RandomSeed = 9130+m;
    [iq, meta, capTruth] = generateSyntheticScenario("aligned_5gnb", custom);
    capID = "rss_gps_"+string(m);
    paths = writeSyntheticCapture(tmpRoot,capID,iq,meta,capTruth);
    lat = lat0 + rad2deg(gpsXY(m,2)/earthRadiusM);
    lon = lon0 + rad2deg(gpsXY(m,1)/(earthRadiusM*cos(deg2rad(lat0))));
    manifestParts{m} = table(capID,string(paths.IQ),string(paths.Metadata), ...
        lat,lon,0,5, ...
        'VariableNames',["CaptureID","IQPath","MetadataPath", ...
        "RxLat","RxLon","RxAltM","RxPositionUncertaintyM"]);
end
manifest = vertcat(manifestParts{:});
runOpts = struct();
runOpts.EnableRSSPlanning = true;
runOpts.Analyze = struct("SIC", struct("MaxIterations", 3));
runOpts.Survey = struct("MinCaptures", 1);
runOpts.RSS = opts;
surveyRun = runMultiCaptureSurvey(manifest, runOpts);
assert(isfield(surveyRun, "PowerTable"));
assert(height(surveyRun.PowerTable) >= 1);
assert(isfield(surveyRun.SurveyResult, "RssPriors"));
assert(any(hypot(surveyRun.PowerTable.RxX_m, ...
    surveyRun.PowerTable.RxY_m)>100), ...
    "GPS receiver positions must not collapse to (0,0) in RSS rows.");
assert(height(surveyRun.SuggestedReceiverStops)>=1, ...
    "Four usable GPS positions should produce ranked receiver-stop options.");

fprintf("Phase 13 RSS-assisted planning tests passed.\n");

function opts = phase13Opts()
opts = struct();
opts.MinMeasurementsPerPCI = 4;
opts.PathLossExponentGrid = 1.6:0.1:4.2;
opts.ShadowingSigmaDb = 8;
opts.GridSpacingM = 50;
opts.SearchMarginM = 900;
opts.MaxResidualDbForRough = 12;
opts.VisibilityPowerDbFS = -95;
opts.NoiseFloorDbFS = -110;
opts.MinExpectedVisiblePCIs = 2;
opts.MinCandidateSpacingM = 100;
opts.AutoCandidateRadiusM = 1400;
opts.MaxSuggestedStops = 5;
end

function [powerTable, captureInfo, truth] = makeRssFixture(shadowSigmaDb, unknownEirp, addOutlier)
pcis = [11; 104; 257];
gnbXY = [-450 -350; 560 -240; 120 620];
rxXY = [ ...
    -1000 -900; ...
      900 -850; ...
     1050  850; ...
     -950  900; ...
        0 -1150; ...
        0  1150; ...
    -1250    0; ...
     1250    0];
ple = [2.0; 2.7; 2.3];
eirp = [-10; -8; -12];
if unknownEirp
    eirp = [-2; -18; -7];
end

captureIDs = "cap_" + string((1:size(rxXY,1)).');
captureInfo = table(captureIDs, rxXY(:,1), rxXY(:,2), zeros(size(rxXY,1),1), ...
    'VariableNames', ["CaptureID","RxX_m","RxY_m","RxZ_m"]);

rows = cell(numel(pcis)*size(rxXY,1), 1);
rowCount = 0;
for m = 1:size(rxXY,1)
    for g = 1:numel(pcis)
        d = max(norm(rxXY(m,:) - gnbXY(g,:)), 10);
        pDb = eirp(g) - 10 * ple(g) * log10(d);
        if shadowSigmaDb > 0
            pDb = pDb + shadowSigmaDb * randn();
        end
        if addOutlier && m == 2 && g == 2
            pDb = pDb - 25;
        end
        rowCount = rowCount + 1;
        rows{rowCount} = table(captureIDs(m), pcis(g), rxXY(m,1), rxXY(m,2), 0, ...
            pDb, -110, pDb + 110, 1, "OK", ...
            'VariableNames', ["CaptureID","PCI","RxX_m","RxY_m","RxZ_m", ...
            "PowerDbFS","NoiseDbFS","SNRdB","NumDetections","QualityFlag"]);
    end
end
powerTable = vertcat(rows{:});

truth = struct();
truth.PCIs = pcis;
truth.GNBXY = gnbXY;
truth.RxXY = rxXY;
end

function [powerTable, captureInfo, truth] = makeCollinearFixture()
[powerTable, captureInfo, truth] = makeRssFixture(0, false, false);
lineX = linspace(-1200, 1200, height(captureInfo)).';
captureInfo.RxX_m = lineX;
captureInfo.RxY_m = zeros(height(captureInfo), 1);
for k = 1:height(powerTable)
    capIdx = find(captureInfo.CaptureID == powerTable.CaptureID(k), 1);
    powerTable.RxX_m(k) = captureInfo.RxX_m(capIdx);
    powerTable.RxY_m(k) = captureInfo.RxY_m(capIdx);
end
end

function err = maxLocationError(priors, truth)
err = 0;
for k = 1:numel(truth.PCIs)
    idx = find(priors.PCI == truth.PCIs(k), 1);
    if isempty(idx)
        err = Inf;
        return;
    end
    err = max(err, norm([priors.EstimatedX_m(idx), priors.EstimatedY_m(idx)] - truth.GNBXY(k,:)));
end
end
