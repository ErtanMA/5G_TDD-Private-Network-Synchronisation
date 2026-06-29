% test_phase22_dataDrivenAssessability
% Assessability must follow model support, not a fixed capture-count rule.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot,"src")));

fprintf("Running Phase 22 data-driven assessability tests...\n");

[measurements,captures] = exactFourPositionFixtureLocal();
opts = struct("ReferencePCI",446,"NumRandomStarts",6, ...
    "SyncThresholdNs",1000,"MaxConditionNumber",Inf);

result = noLocationSurveyCheck(measurements,captures,opts);

assert(result.FitInfo.NumMeasurements==result.FitInfo.NumUnknowns);
assert(result.FitInfo.DegreesOfFreedom==0);
assert(result.FitInfo.Status=="SUSPECT", ...
    "An exactly determined four-position fit should be computed and marked SUSPECT.");
assert(all(ismember([446;559;9],result.TimingEstimates.PCI)));
assert(any(isfinite(result.TimingEstimates.EstimatedRelativeTxOffsetNs( ...
    ismember(result.TimingEstimates.PCI,[446;559;9])))));

row489 = result.TimingEstimates(result.TimingEstimates.PCI==489,:);
assert(height(row489)==1);
assert(row489.TimingStatus=="NOT_ASSESSABLE");
assert(ismember(489,result.FitInfo.ExcludedPCIs));

stationary = table(489,50,170,90,60,0,10,"ASSESSABLE","fixture", ...
    'VariableNames',["PCI","StationaryTimingInstabilityRmsNs", ...
    "StationaryTimingInstabilityPeakToPeakNs", ...
    "StationaryTimingInstabilityMaxAbsNs", ...
    "ExpectedMeasurementNoiseRmsNs","ExcessStationaryRmsNs", ...
    "NumObservations","StationaryAssessmentStatus", ...
    "StationaryAssessmentReason"]);
merged = applyStationaryInstabilityToSurveyResult(result,stationary,opts);
merged489 = merged.TimingEstimates(merged.TimingEstimates.PCI==489,:);
assert(merged489.TimingStatus=="NOT_ASSESSABLE", ...
    "Stationary evidence must not overwrite a geometry-excluded PCI verdict.");

fprintf("Phase 22 data-driven assessability tests passed.\n");

function [measurements,captures] = exactFourPositionFixtureLocal()
rx = [0 0; 140 20; 30 170; -120 80];
pcis = [446 559 9];
cellXY = [500 700; -650 300; 250 -800];
biasNs = [0 220 -180];
c = 299792458;
parts = cell(numel(pcis)*size(rx,1)+1,1);
count = 0;
for m = 1:size(rx,1)
    captureID = "p"+compose("%02d",m);
    for g = 1:numel(pcis)
        distanceNs = norm(cellXY(g,:)-rx(m,:))/c*1e9;
        phaseNs = mod(2e6 + distanceNs + biasNs(g),10e6);
        count = count+1;
        parts{count} = table(captureID,pcis(g),phaseNs,20, ...
            'VariableNames',["CaptureID","PCI","FramePhaseNs","TimingStdNs"]);
    end
end

% A one-position cell must be reported but must not invalidate the fit.
count = count+1;
parts{count} = table("p01",489,2.003e6,20, ...
    'VariableNames',["CaptureID","PCI","FramePhaseNs","TimingStdNs"]);
measurements = vertcat(parts{1:count});

captures = table("p"+compose("%02d",(1:size(rx,1)).'), ...
    rx(:,1),rx(:,2),zeros(size(rx,1),1),5*ones(size(rx,1),1), ...
    'VariableNames',["CaptureID","RxX_m","RxY_m","RxZ_m", ...
    "RxPositionUncertaintyM"]);
end
