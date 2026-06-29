addpath(genpath(fullfile(fileparts(mfilename("fullpath")), "..", "src")));

fprintf("Running Phase 9 unknown-location survey tests...\n");

surveyOpts = struct();
surveyOpts.NumRandomStarts = 10;
surveyOpts.MinCaptures = 5;
surveyOpts.MaxFitRMSNsForPass = 120;
surveyOpts.MaxBiasUncertaintyForPassNs = 180;
surveyOpts.HardFailOffsetNs = 3000;

testAlignedSurvey(surveyOpts);
testThreeMicrosecondFault(surveyOpts);
testWarningScaleFault(surveyOpts);
testCommonModeInvisible(surveyOpts);
testInsufficientCaptures(surveyOpts);
testLatLonReceiverInput(surveyOpts);
testEmptyDetections(surveyOpts);
testOneVisibleGNB(surveyOpts);
testMissingReferenceEnforced(surveyOpts);
testMissingCaptureInfoNotAssessable(surveyOpts);
testDuplicateCaptureInfoNotAssessable(surveyOpts);

fprintf("Phase 9 unknown-location survey tests passed.\n");

function testAlignedSurvey(surveyOpts)
[meas, caps, ~] = buildSyntheticSurveyDataset("aligned_5gnb", struct( ...
    "MeasurementNoiseNs", 5, "RandomSeed", 101));
result = noLocationSurveyCheck(meas, caps, surveyOpts);
assert(result.FitInfo.Status == "OK" || result.FitInfo.Status == "SUSPECT");

targets = result.TimingEstimates(result.TimingEstimates.SurveyRole == "TARGET", :);
assert(all(abs(targets.EstimatedRelativeTxOffsetNs) < 180), ...
    "Aligned survey estimated a large relative timing error.");
assert(~any(targets.TimingStatus == "FAIL"), ...
    "Aligned survey should not fail any target.");
assert(ismember(result.FitInfo.ReferencePCI, unique(meas.PCI)), ...
    "Auto-selected reference PCI should be one of the detected PCIs.");
end

function testThreeMicrosecondFault(surveyOpts)
[meas, caps, truth] = buildSyntheticSurveyDataset("offset_3us_one", struct( ...
    "MeasurementNoiseNs", 5, "RandomSeed", 202));
result = noLocationSurveyCheck(meas, caps, surveyOpts);

faultPCI = truth.GNB.PCI(truth.GNB.RelativeTxOffsetNs == 3000);
row = result.TimingEstimates(result.TimingEstimates.PCI == faultPCI, :);
assert(~isempty(row), "Fault PCI missing from survey result.");
assert(abs(row.EstimatedRelativeTxOffsetNs - 3000) < 250, ...
    "3 us fault estimate is not close to truth.");
assert(row.TimingStatus == "FAIL" || abs(row.EstimatedRelativeTxOffsetNs) >= surveyOpts.HardFailOffsetNs, ...
    "3 us fault should be failed or exceed hard-fail estimate.");
end

function testWarningScaleFault(surveyOpts)
[meas, caps, truth] = buildSyntheticSurveyDataset("offset_250ns_one", struct( ...
    "MeasurementNoiseNs", 5, "RandomSeed", 303));
result = noLocationSurveyCheck(meas, caps, surveyOpts);

faultPCI = truth.GNB.PCI(truth.GNB.RelativeTxOffsetNs == 250);
row = result.TimingEstimates(result.TimingEstimates.PCI == faultPCI, :);
assert(~isempty(row), "250 ns PCI missing from survey result.");
assert(abs(row.EstimatedRelativeTxOffsetNs - 250) < 180, ...
    "250 ns warning-scale offset estimate is not close enough to truth.");
assert(ismember(row.TimingStatus, ["PASS","SUSPECT"]), ...
    "250 ns case should be pass/suspect depending uncertainty, not fail.");
end

function testCommonModeInvisible(surveyOpts)
[meas, caps, ~] = buildSyntheticSurveyDataset("common_mode_3us", struct( ...
    "MeasurementNoiseNs", 5, "RandomSeed", 404));
result = noLocationSurveyCheck(meas, caps, surveyOpts);

targets = result.TimingEstimates(result.TimingEstimates.SurveyRole == "TARGET", :);
assert(all(abs(targets.EstimatedRelativeTxOffsetNs) < 180), ...
    "Common-mode timing shift should be invisible in relative survey mode.");
assert(~any(targets.TimingStatus == "FAIL"), ...
    "Common-mode timing shift should not create a relative failure.");
end

function testInsufficientCaptures(surveyOpts)
[meas, caps, ~] = buildSyntheticSurveyDataset("insufficient_captures", struct( ...
    "MeasurementNoiseNs", 5, "RandomSeed", 505));
result = noLocationSurveyCheck(meas, caps, surveyOpts);
assert(result.FitInfo.Status == "NOT_ASSESSABLE", ...
    "Too few receiver positions must be not assessable.");
assert(all(result.TimingEstimates.TimingStatus == "NOT_ASSESSABLE"), ...
    "All rows should be not assessable for insufficient survey geometry.");
end

function testLatLonReceiverInput(surveyOpts)
[meas, capsXY, ~] = buildSyntheticSurveyDataset("aligned_5gnb", struct( ...
    "MeasurementNoiseNs", 5, "RandomSeed", 606));

lat0 = 52.0;
lon0 = 5.0;
earthRadiusM = 6371000;
capsLL = table();
capsLL.CaptureID = capsXY.CaptureID;
capsLL.RxLat = lat0 + rad2deg(capsXY.RxY_m / earthRadiusM);
capsLL.RxLon = lon0 + rad2deg(capsXY.RxX_m / (earthRadiusM * cos(deg2rad(lat0))));
capsLL.RxAltM = zeros(height(capsXY), 1);
capsLL.RxPositionUncertaintyM = capsXY.RxPositionUncertaintyM;

result = noLocationSurveyCheck(meas, capsLL, surveyOpts);
assert(height(result.RelativeMeasurements) > 0, ...
    "Lat/lon receiver positions should be converted into local survey coordinates.");
end

function testEmptyDetections(surveyOpts)
meas = table(strings(0,1), zeros(0,1), zeros(0,1), ...
    'VariableNames', ["CaptureID","PCI","FramePhaseNs"]);
caps = table("cap_1", 0, 0, 0, ...
    'VariableNames', ["CaptureID","RxX_m","RxY_m","RxZ_m"]);
result = noLocationSurveyCheck(meas, caps, surveyOpts);
assert(result.FitInfo.Status == "NOT_ASSESSABLE");
assert(result.FitInfo.Reason == "No detected PCIs.");
assert(height(result.TimingEstimates) == 0);
end

function testOneVisibleGNB(surveyOpts)
[meas, caps, truth] = buildSyntheticSurveyDataset("aligned_5gnb", struct( ...
    "MeasurementNoiseNs", 5, "RandomSeed", 707));
onePCI = truth.GNB.PCI(1);
meas = meas(meas.PCI == onePCI, :);
result = noLocationSurveyCheck(meas, caps, surveyOpts);
assert(result.FitInfo.Status == "NOT_ASSESSABLE");
assert(height(result.TimingEstimates) == 1);
assert(result.TimingEstimates.PCI == onePCI);
assert(result.DecisionTable.Verdict == "NOT_ASSESSABLE");
end

function testMissingReferenceEnforced(surveyOpts)
[meas, caps, truth] = buildSyntheticSurveyDataset("aligned_5gnb", struct( ...
    "MeasurementNoiseNs", 5, "RandomSeed", 808));
surveyOpts.ReferencePCI = truth.ReferencePCI;
dropID = meas.CaptureID(1);
meas = meas(~(meas.CaptureID == dropID & meas.PCI == truth.ReferencePCI), :);
result = noLocationSurveyCheck(meas, caps, surveyOpts);
assert(result.FitInfo.Status == "NOT_ASSESSABLE");
assert(numel(result.FitInfo.MissingReferenceCaptureIDs) == 1);
assert(all(result.TimingEstimates.TimingStatus == "NOT_ASSESSABLE"));
end

function testMissingCaptureInfoNotAssessable(surveyOpts)
[meas, caps, ~] = buildSyntheticSurveyDataset("aligned_5gnb", struct( ...
    "MeasurementNoiseNs", 5, "RandomSeed", 909));
caps = caps(2:end, :);
result = noLocationSurveyCheck(meas, caps, surveyOpts);
assert(result.FitInfo.Status == "NOT_ASSESSABLE");
assert(contains(result.FitInfo.Reason, "Missing receiver-position"));
end

function testDuplicateCaptureInfoNotAssessable(surveyOpts)
[meas, caps, ~] = buildSyntheticSurveyDataset("aligned_5gnb", struct( ...
    "MeasurementNoiseNs", 5, "RandomSeed", 1001));
caps = [caps; caps(1,:)];
result = noLocationSurveyCheck(meas, caps, surveyOpts);
assert(result.FitInfo.Status == "NOT_ASSESSABLE");
assert(contains(result.FitInfo.Reason, "Duplicate receiver-position"));
end
