% test_phase19_realSurveyWorkflow
% Regression tests for fixed-profile B210 survey capture validation.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot,"src")));

fprintf("Running Phase 19 real-survey workflow tests...\n");

cfg = validCampaignConfigLocal();
report = validateB210CampaignCapture(cfg,30720000,0);
assert(report.IsValid);
assert(report.ExpectedSampleCount==30720000);

assertErrorLocal(@() validateB210CampaignCapture(cfg,30719999,0), ...
    "B210Capture:RejectedSampleCount");
assertErrorLocal(@() validateB210CampaignCapture(cfg,30720000,1), ...
    "B210Capture:RejectedOverrun");

badRate = cfg;
badRate.DecimationFactor = 2;
badRate.SampleRateHz = 15.36e6;
assertErrorLocal(@() validateB210CampaignCapture(badRate), ...
    "B210Capture:CampaignProfileMismatch");

badGPS = cfg;
badGPS.RxLat = NaN;
assertErrorLocal(@() validateB210CampaignCapture(badGPS), ...
    "B210Capture:MissingGPS");

duplicateManifest = table( ...
    ["survey_20260621_p01";"survey_20260621_p02"], ...
    ["missing_1.mat";"missing_2.mat"], ...
    ["missing_1.json";"missing_2.json"], ...
    [52.245567;52.245567],[6.853050;6.853050],[33;33],[5;5], ...
    'VariableNames',["CaptureID","IQPath","MetadataPath", ...
    "RxLat","RxLon","RxAltM","RxPositionUncertaintyM"]);
assertErrorLocal(@() runMultiCaptureSurvey(duplicateManifest), ...
    "runMultiCaptureSurvey:DuplicateReceiverPosition");

tmpRoot = tempname;
mkdir(tmpRoot);
cleanup = onCleanup(@() rmdir(tmpRoot,"s"));
custom = struct("NumGNBs",2,"PCIs",[11 12],"CFOHz",[0 0], ...
    "FrameOffsetsNs",[0 0],"GNBPowerdB",[0 -2], ...
    "LocationUncertaintyM",[30 30],"SiteDistanceM",[100 100], ...
    "SNRdB",35,"DurationMs",20,"RandomSeed",1919);
[iq,meta,truth] = generateSyntheticScenario("aligned_5gnb",custom);
paths = writeSyntheticCapture(tmpRoot,"invalid_capture",iq,meta,truth);
raw = jsondecode(fileread(paths.Metadata));
raw.capture_validated = false;
raw.overrun_count = 1;
fid = fopen(paths.Metadata,"w");
assert(fid>=0);
fileCleanup = onCleanup(@() fclose(fid));
fprintf(fid,"%s\n",jsonencode(raw,"PrettyPrint",true));
clear fileCleanup;
invalidManifest = table("invalid_capture",string(paths.IQ), ...
    string(paths.Metadata),0,0,0,5, ...
    'VariableNames',["CaptureID","IQPath","MetadataPath", ...
    "RxX_m","RxY_m","RxZ_m","RxPositionUncertaintyM"]);
assertErrorLocal(@() runMultiCaptureSurvey(invalidManifest), ...
    "runMultiCaptureSurvey:UnvalidatedCapture");

fprintf("Phase 19 real-survey workflow tests passed.\n");

function cfg = validCampaignConfigLocal()
cfg = struct();
cfg.CampaignID = "survey_20260621";
cfg.PositionID = "p01";
cfg.RepeatIndex = 1;
cfg.CaptureID = "survey_20260621_p01_01";
cfg.MasterClockRateHz = 30.72e6;
cfg.DecimationFactor = 1;
cfg.SampleRateHz = 30.72e6;
cfg.ActualSampleRateHz = 30.72e6;
cfg.BandwidthHz = 30e6;
cfg.DurationSeconds = 1;
cfg.SamplesPerFrame = 153600;
cfg.CenterFrequencyHz = 3.70272e9;
cfg.Gain = 45;
cfg.RxLat = 52.245567;
cfg.RxLon = 6.853050;
cfg.RxAltM = 33;
cfg.RxPositionUncertaintyM = 5;
end

function assertErrorLocal(fn,expectedID)
didError = false;
try
    fn();
catch ME
    didError = string(ME.identifier)==string(expectedID);
end
assert(didError,"Expected error identifier %s.",expectedID);
end
