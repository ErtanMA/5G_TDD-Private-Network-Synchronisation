% test_phase20_stationaryRepeatedSurvey
% Regression tests for fixed-position repeats and stationary instability.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot,"src")));

fprintf("Running Phase 20 stationary-repeat survey tests...\n");

[windowMeasurements,repeatManifest] = repeatedFixtureLocal();
aggregate = aggregateRepeatedTimingMeasurements( ...
    windowMeasurements,repeatManifest,struct( ...
    "ReferencePCI",11,"ReferenceGSCN",7987, ...
    "MinRepeatsPerPosition",3));

assert(height(aggregate.PositionInfo)==6);
assert(numel(unique(aggregate.PositionInfo.CaptureID))==6);
assert(height(aggregate.AggregatedMeasurementTable)==12);
assert(isempty(aggregate.RejectedWindowMeasurements));
assert(aggregate.ReferenceCandidateRanking.PCI(1)==11);
assert(numel(unique(aggregate.AggregatedMeasurementTable.CaptureID))==6, ...
    "Repeats must collapse to one geometry capture ID per physical position.");

targetStationary = aggregate.StationaryInstability( ...
    aggregate.StationaryInstability.PCI==12,:);
assert(height(targetStationary)==1);
assert(targetStationary.NumCaptures==18);
assert(targetStationary.NumPositions==6);
assert(targetStationary.StationaryTimingInstabilityMaxAbsNs>=40);
assert(targetStationary.StationaryTimingInstabilityMaxAbsNs<=50);

survey = noLocationSurveyCheck( ...
    aggregate.AggregatedMeasurementTable,aggregate.PositionInfo,struct( ...
    "ReferencePCI",11,"MinCaptures",5,"NumRandomStarts",3, ...
    "SyncThresholdNs",1000,"MaxFitRMSNsForPass",5000));
targetIndex = find(survey.TimingEstimates.PCI==12,1);
assert(~isempty(targetIndex));
survey.TimingEstimates.TimingInstabilityMaxAbsNs(targetIndex) = 9999;
survey.TimingEstimates.WorstCaseRelativeTimingNs(targetIndex) = ...
    abs(survey.TimingEstimates.StaticOffsetNs(targetIndex))+9999;

updated = applyStationaryInstabilityToSurveyResult( ...
    survey,aggregate.StationaryInstability,struct("SyncThresholdNs",1000));
target = updated.TimingEstimates(updated.TimingEstimates.PCI==12,:);
expectedEnvelope = abs(target.StaticOffsetNs)+ ...
    target.StationaryTimingInstabilityMaxAbsNs;
assert(abs(target.WorstCaseRelativeTimingNs-expectedEnvelope)<1e-9);
assert(target.WorstCaseRelativeTimingNs< ...
    abs(target.StaticOffsetNs)+1000, ...
    "The old survey residual must not be reused as stationary jitter.");
assert(ismember("SurveyModelResidualMaxAbsNs", ...
    string(updated.TimingEstimates.Properties.VariableNames)));
assert(contains(updated.VerdictBasis,"stationary repeated-capture"));

fprintf("Phase 20 stationary-repeat survey tests passed.\n");

function [measurements,manifest] = repeatedFixtureLocal()
numPositions = 6;
numRepeats = 3;
numWindows = 4;
frameNs = 10e6;
relativeCentersNs = [120 210 350 500 690 910];
jitterNs = [-40 -10 20 30];
measurementParts = cell(numPositions*numRepeats*numWindows*2,1);
manifestParts = cell(numPositions*numRepeats,1);
measurementCount = 0;
manifestCount = 0;

for p = 1:numPositions
    positionID = "p"+compose("%02d",p);
    lat = 52.242940+0.0007*(p-1);
    lon = 6.852839+0.0005*mod(p-1,3);
    for r = 1:numRepeats
        captureID = "survey_20260621_"+positionID+"_"+compose("%02d",r);
        manifestCount = manifestCount+1;
        manifestParts{manifestCount} = table( ...
            captureID,positionID,r,"C:\fixture\"+captureID+"_iq.mat", ...
            "C:\fixture\"+captureID+"_meta.json",lat,lon,30,5, ...
            'VariableNames',["CaptureID","PositionID","RepeatIndex", ...
            "IQPath","MetadataPath","RxLat","RxLon","RxAltM", ...
            "RxPositionUncertaintyM"]);

        for w = 1:numWindows
            observationID = captureID+"_w"+compose("%02d",w)+"_g7987";
            commonPhase = mod(1.1e6*p+2.3e5*r+7.1e4*w,frameNs);
            targetPhase = mod(commonPhase+relativeCentersNs(p)+ ...
                jitterNs(w),frameNs);
            measurementCount = measurementCount+1;
            measurementParts{measurementCount} = measurementRowLocal( ...
                captureID,positionID,r,observationID,11,commonPhase);
            measurementCount = measurementCount+1;
            measurementParts{measurementCount} = measurementRowLocal( ...
                captureID,positionID,r,observationID,12,targetPhase);
        end
    end
end
measurements = vertcat(measurementParts{1:measurementCount});
manifest = vertcat(manifestParts{1:manifestCount});
end

function row = measurementRowLocal(captureID,positionID,repeatIndex, ...
        observationID,pci,phaseNs)
cellKey = "GSCN7987_PCI"+string(pci);
row = table( ...
    captureID,positionID,repeatIndex,observationID,cellKey,7987,pci, ...
    phaseNs,10,3.70272e9,0,3.70272e9,0.20,0.40,100, ...
    'VariableNames',["CaptureID","PositionID","RepeatIndex", ...
    "ObservationID","CellKey","GSCN","PCI","FramePhaseNs", ...
    "TimingStdNs","CenterFrequencyHz","SelectedSSBOffsetHz", ...
    "HypothesizedSSBCenterHz","MeanPSSMetric","MeanSSSMetric", ...
    "MeanCFOHz"]);
end
