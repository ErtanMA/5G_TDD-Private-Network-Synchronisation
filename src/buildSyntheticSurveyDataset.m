function [measurementTable, captureInfo, truth] = buildSyntheticSurveyDataset(scenarioName, opts)
%buildSyntheticSurveyDataset Create synthetic multi-position survey data.
%
%   [measurementTable, captureInfo, truth] =
%   buildSyntheticSurveyDataset(scenarioName, opts) generates timing
%   measurements that look like the output of several independent B210
%   captures after PSS/SSS/timing recovery. The gNB positions are hidden from
%   the estimator and kept only in truth for validation.

arguments
    scenarioName (1,1) string
    opts struct = struct()
end

base = defaultSyntheticSurveyOptionsLocal();
opts = mergeStructsLocal(base, opts);
opts = configureScenarioLocal(lower(scenarioName), opts);

rng(opts.RandomSeed, "twister");

pcis = opts.PCIs(:);
numG = numel(pcis);
numCaptures = opts.NumCaptures;
c = opts.SpeedOfLightMps;
frameNs = opts.FramePeriodMs * 1e6;

gnbXY = opts.GNBPositionsM;
rxXY = opts.RxPositionsM;
if isempty(rxXY)
    rxXY = receiverPathLocal(numCaptures, opts);
end
numCaptures = size(rxXY, 1);

rxClockNs = rand(numCaptures, 1) * frameNs;
rows = cell(numCaptures * numG, 1);
rowCount = 0;

for m = 1:numCaptures
    for g = 1:numG
        distanceM = norm(gnbXY(g,:) - rxXY(m,:));
        noiseNs = opts.MeasurementNoiseNs * randn();
        arrivalNs = rxClockNs(m) + opts.TxOffsetsNs(g) + distanceM / c * 1e9 + noiseNs;
        framePhaseNs = mod(arrivalNs, frameNs);

        rowCount = rowCount + 1;
        rows{rowCount} = table( ...
            "cap_" + string(m), ...
            pcis(g), ...
            framePhaseNs, ...
            opts.MeasurementNoiseNs, ...
            opts.CenterFrequencyHz, ...
            20 + 3*randn(), ...
            'VariableNames', ["CaptureID","PCI","FramePhaseNs", ...
            "TimingStdNs","CenterFrequencyHz","SNRdB"]);
    end
end

measurementTable = vertcat(rows{1:rowCount});
captureInfo = table( ...
    "cap_" + string((1:numCaptures).'), ...
    rxXY(:,1), rxXY(:,2), zeros(numCaptures, 1), ...
    repmat(opts.ReceiverPositionUncertaintyM, numCaptures, 1), ...
    'VariableNames', ["CaptureID","RxX_m","RxY_m","RxZ_m", ...
    "RxPositionUncertaintyM"]);

truth = struct();
truth.ScenarioName = scenarioName;
truth.ReferencePCI = pcis(1);
truth.GNB = table( ...
    pcis, gnbXY(:,1), gnbXY(:,2), opts.TxOffsetsNs(:), ...
    opts.TxOffsetsNs(:) - opts.TxOffsetsNs(1), ...
    'VariableNames', ["PCI","TrueX_m","TrueY_m","TxOffsetNs", ...
    "RelativeTxOffsetNs"]);
truth.Captures = captureInfo;
truth.RxClockPhaseNs = rxClockNs;
truth.Statement = "Synthetic survey data: gNB positions are generated for truth only and are not passed to the estimator.";

end

function opts = defaultSyntheticSurveyOptionsLocal()
opts = struct();
opts.RandomSeed = 20260603;
opts.SpeedOfLightMps = 299792458;
opts.FramePeriodMs = 10;
opts.CenterFrequencyHz = 3.77e9;
opts.NumCaptures = 8;
opts.PCIs = [11 104 257 503 777];
opts.TxOffsetsNs = [0 0 0 0 0];
opts.MeasurementNoiseNs = 20;
opts.ReceiverPositionUncertaintyM = 5;
opts.AreaRadiusM = 900;
opts.SurveyRadiusM = 1200;
opts.RxPositionsM = [];
opts.GNBPositionsM = [ ...
    -450 -250; ...
     350 -350; ...
     520  250; ...
    -150  500; ...
     -50   50];
end

function opts = configureScenarioLocal(name, opts)
switch name
    case "aligned_5gnb"
        % Defaults.
    case "offset_250ns_one"
        opts.TxOffsetsNs(4) = 250;
    case "offset_1us_one"
        opts.TxOffsetsNs(4) = 1000;
    case "offset_3us_one"
        opts.TxOffsetsNs(4) = 3000;
    case "offset_3us_two"
        opts.TxOffsetsNs(3) = 3000;
        opts.TxOffsetsNs(4) = -3000;
    case "common_mode_3us"
        opts.TxOffsetsNs(:) = 3000;
    case "insufficient_captures"
        opts.NumCaptures = 3;
    case "noisy_250ns_one"
        opts.TxOffsetsNs(4) = 250;
        opts.MeasurementNoiseNs = 80;
    case "weak_geometry_3us"
        opts.TxOffsetsNs(4) = 3000;
        opts.RxPositionsM = [-900 -20; -450 10; 0 -10; 450 20; 900 -15];
    otherwise
        error("buildSyntheticSurveyDataset:UnknownScenario", ...
            "Unknown synthetic survey scenario: %s", name);
end

opts.PCIs = opts.PCIs(:);
opts.TxOffsetsNs = opts.TxOffsetsNs(:);
if size(opts.GNBPositionsM, 1) ~= numel(opts.PCIs)
    opts.GNBPositionsM = opts.GNBPositionsM(1:numel(opts.PCIs), :);
end
end

function rxXY = receiverPathLocal(numCaptures, opts)
angles = linspace(0, 2*pi, numCaptures+1).';
angles(end) = [];
rxXY = opts.SurveyRadiusM * [cos(angles), sin(angles)];

if numCaptures >= 6
    rxXY(2:2:end,:) = 0.65 * rxXY(2:2:end,:);
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
