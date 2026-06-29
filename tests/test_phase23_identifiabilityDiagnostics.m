% test_phase23_identifiabilityDiagnostics

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot,"src")));

fprintf("Running Phase 23 identifiability diagnostics tests...\n");

[measurements,positions] = fixtureLocal();
opts = struct("ReferencePCI",446,"NumRandomStarts",5, ...
    "SyncThresholdNs",1000);
result = noLocationSurveyCheck(measurements,positions,opts);

assert(isfield(result.FitInfo,"NumericalDiagnostics"));
numeric = result.FitInfo.NumericalDiagnostics;
assert(numeric.NumParameters==result.FitInfo.NumUnknowns);
assert(numel(numeric.ParameterNames)==numeric.NumParameters);
assert(numel(numeric.SingularValues)==numeric.NumParameters);
assert(isfinite(numeric.ConditionNumber));

starts = result.FitInfo.StartDiagnostics;
assert(starts.NumStarts==5);
assert(height(starts.PerTarget)==2);
assert(all(ismember(starts.PerTarget.PCI,[559;9])));
assert(starts.NumNearOptimalStarts>=1);

extra = positions(1,:);
extra.CaptureID = "unused_position";
withUnused = [positions;extra];
stability = jackknifeSurveyStability(measurements,withUnused,opts, ...
    struct("OnlyMeasuredPositions",true));
assert(stability.NumInputPositions==height(withUnused));
assert(stability.NumFolds==height(positions));
assert(stability.UnusedPositionIDs=="unused_position");

fprintf("Phase 23 identifiability diagnostics tests passed.\n");

function [measurements,positions] = fixtureLocal()
rx = [0 0;120 30;-40 150;-130 -60;80 -140;210 100];
cellXY = [500 700;-600 250;200 -800];
pcis = [446 559 9];
biasNs = [0 250 -180];
c = 299792458;
parts = cell(size(rx,1)*numel(pcis),1);
count = 0;
for m = 1:size(rx,1)
    captureID = "p"+compose("%02d",m);
    for g = 1:numel(pcis)
        phaseNs = mod(2e6+norm(cellXY(g,:)-rx(m,:))/c*1e9+biasNs(g),10e6);
        count = count+1;
        parts{count} = table(captureID,pcis(g),phaseNs,20, ...
            'VariableNames',["CaptureID","PCI","FramePhaseNs","TimingStdNs"]);
    end
end
measurements = vertcat(parts{:});
positions = table("p"+compose("%02d",(1:size(rx,1)).'), ...
    rx(:,1),rx(:,2),zeros(size(rx,1),1),5*ones(size(rx,1),1), ...
    'VariableNames',["CaptureID","RxX_m","RxY_m","RxZ_m", ...
    "RxPositionUncertaintyM"]);
end
