function [iq, meta, truth] = generateSyntheticScenario(scenarioName, opts)
%generateSyntheticScenario Generate controlled synthetic multi-gNB IQ.
%
%   [iq, meta, truth] = generateSyntheticScenario(scenarioName) creates a
%   complex baseband IQ capture and ground-truth table for the requested
%   scenario. The generator is intentionally transparent: it uses
%   standards-correct PSS/SSS sequences from 5G Toolbox, but the timing,
%   TDD-pattern, CFO, and multi-gNB composition are implemented here.

arguments
    scenarioName (1,1) string
    opts struct = struct()
end

base = defaultSyntheticOptions();
opts = mergeStructs(base, opts);
opts = configureScenario(scenarioName, opts);

rng(opts.RandomSeed, "twister");

fs = opts.SampleRateHz;
numSamples = round(opts.DurationMs * 1e-3 * fs);
numGNBs = opts.NumGNBs;

iq = complex(zeros(numSamples, 1));
componentPowers = zeros(numGNBs, 1);

for g = 1:numGNBs
    component = synthesizeGNBComponent(g, opts, numSamples);
    componentPowers(g) = mean(abs(component).^2);
    iq = iq + component;
end

iq = addToneInterference(iq, opts);

signalPower = mean(abs(iq).^2);
noisePower = signalPower / db2powLocal(opts.SNRdB);
noise = sqrt(noisePower/2) * (randn(numSamples,1) + 1i*randn(numSamples,1));
iq = iq + noise;
iq = applyReceiverImpairments(iq, opts);

meta = struct();
meta.SampleRateHz = opts.SampleRateHz;
meta.CenterFrequencyHz = opts.CenterFrequencyHz;
meta.BandwidthHz = opts.BandwidthHz;
meta.RxLat = 52.0;
meta.RxLon = 5.0;
meta.RxAltM = 0;
meta.ClockSource = "internal";
meta.GnssLocked = false;
meta.ScenarioName = scenarioName;
meta.DurationMs = opts.DurationMs;

truth = buildTruthTable(scenarioName, opts, componentPowers, noisePower);

end

function opts = configureScenario(name, opts)
name = lower(name);

switch name
    case "aligned_5gnb"
        % Defaults already describe this case.

    case "offset_100ns_one"
        opts.FrameOffsetsNs(4) = 100;

    case "offset_250ns_one"
        opts.FrameOffsetsNs(4) = 250;

    case "offset_500ns_one"
        opts.FrameOffsetsNs(4) = 500;

    case "offset_1us_one"
        opts.FrameOffsetsNs(4) = 1000;

    case "offset_1p5us_one"
        opts.FrameOffsetsNs(4) = 1500;

    case "offset_3us_two"
        opts.FrameOffsetsNs(3) = 3000;
        opts.FrameOffsetsNs(4) = -3000;

    case "offset_negative3us_one"
        opts.FrameOffsetsNs(4) = -3000;

    case "mixed_offsets"
        opts.FrameOffsetsNs = [-250 0 500 1500 3000];

    case "common_mode_3us"
        opts.FrameOffsetsNs = 3000 * ones(1, opts.NumGNBs);

    case "low_snr_250ns"
        opts.FrameOffsetsNs(4) = 250;
        opts.SNRdB = 5;

    case "cfo_anomaly_one"
        opts.CFOHz(5) = 250;

    case "wrong_tdd_pattern"
        opts.AddWrongUplinkDLEnergy = true;

    case "wrong_special_slot"
        opts.AddWrongSpecialSlotEnergy = true;

    case "single_gnb_not_assessable"
        opts.NumGNBs = 1;
        opts.PCIs = opts.PCIs(1);
        opts.CFOHz = 0;
        opts.FrameOffsetsNs = 0;
        opts.GNBPowerdB = 0;
        opts.LocationUncertaintyM = 30;
        opts.SiteDistanceM = 100;

    case "large_geometry_uncertainty"
        opts.FrameOffsetsNs(4) = 250;
        opts.LocationUncertaintyM = 500 * ones(1, opts.NumGNBs);

    case "rf_impairment_moderate"
        opts.NumGNBs = 3;
        opts.PCIs = [11 12 13];
        opts.SNRdB = 26;
        opts.FrameOffsetsNs = [0 250 1000];
        opts.CFOHz = [0 80 -60];
        opts.GNBPowerdB = [0 -1.5 -3];
        opts.LocationUncertaintyM = 30 * ones(1, opts.NumGNBs);
        opts.SiteDistanceM = 100 * ones(1, opts.NumGNBs);
        opts.MultipathDelaysSamples = [0 1.5 5.25];
        opts.MultipathGainsdB = [0 -12 -20];
        opts.MultipathPhasesRad = [0 0.8 -0.4];
        opts.IQGainImbalanceDB = 0.3;
        opts.IQPhaseImbalanceDeg = 1.5;
        opts.DCOffset = 0.01 + 0.005i;
        opts.PhaseNoiseStepStdDeg = 0.01;
        opts.ToneInterferenceOffsetHz = 2.8e6;
        opts.ToneInterferencePowerdB = -32;
        opts.ClippingLevelRMS = 5;

    case "rf_impairment_harsh"
        opts.NumGNBs = 3;
        opts.PCIs = [11 12 13];
        opts.SNRdB = 10;
        opts.FrameOffsetsNs = [0 1000 3000];
        opts.CFOHz = [0 250 -250];
        opts.GNBPowerdB = [0 -6 -10];
        opts.LocationUncertaintyM = 30 * ones(1, opts.NumGNBs);
        opts.SiteDistanceM = 100 * ones(1, opts.NumGNBs);
        opts.MultipathDelaysSamples = [0 4.5 14.75 31];
        opts.MultipathGainsdB = [0 -5 -10 -18];
        opts.MultipathPhasesRad = [0 -0.7 2.1 0.4];
        opts.IQGainImbalanceDB = 1.2;
        opts.IQPhaseImbalanceDeg = 6;
        opts.DCOffset = 0.04 - 0.03i;
        opts.PhaseNoiseStepStdDeg = 0.05;
        opts.ToneInterferenceOffsetHz = -3.1e6;
        opts.ToneInterferencePowerdB = -18;
        opts.ClippingLevelRMS = 2.8;

    otherwise
        error("generateSyntheticScenario:UnknownScenario", ...
            "Unknown synthetic scenario '%s'. Use listSyntheticScenarios().", name);
end
end

function component = synthesizeGNBComponent(g, opts, numSamples)
fs = opts.SampleRateHz;
component = complex(zeros(numSamples, 1));

powerScale = sqrt(db2powLocal(opts.GNBPowerdB(g)));
dl = synthesizeTDDEnergy(opts, numSamples) * powerScale;
component = component + dl;

ssb = buildSyntheticSSB(opts.PCIs(g), opts);
frameSamples = round(opts.FramePeriodMs * 1e-3 * fs);
offsetSamples = opts.FrameOffsetsNs(g) * 1e-9 * fs;
ssbBaseSamples = opts.SSBStartOffsetUs * 1e-6 * fs + ...
    caseCSSBIndexOffsetSamplesLocal(opts);

for start0 = 0:frameSamples:(numSamples-1)
    startSample0 = start0 + ssbBaseSamples + offsetSamples;
    component = addBurstAtSample(component, ssb * opts.SSBPower * powerScale, startSample0);
end

if opts.CFOHz(g) ~= 0
    n = (0:numSamples-1).';
    component = component .* exp(1i*2*pi*opts.CFOHz(g)*n/fs);
end

component = applyMultipathChannel(component, opts);
end

function offsetSamples = caseCSSBIndexOffsetSamplesLocal(opts)
caseCStartSymbols = [2 8 16 22 30 36 44 50];
carrier = nrCarrierConfig;
carrier.NSizeGrid = opts.NSizeGrid;
carrier.SubcarrierSpacing = opts.SCSkHz;
info = nrOFDMInfo(carrier, "SampleRate", opts.SampleRateHz);
symbolLengths = double(info.SymbolLengths(:));
symbolsPerSlot = numel(symbolLengths);
absoluteSymbol = caseCStartSymbols(opts.SSBIndex+1);
slotIndex = floor(absoluteSymbol/symbolsPerSlot);
symbolInSlot = mod(absoluteSymbol,symbolsPerSlot);
offsetSamples = slotIndex*sum(symbolLengths) + ...
    sum(symbolLengths(1:symbolInSlot));
end

function dl = synthesizeTDDEnergy(opts, numSamples)
fs = opts.SampleRateHz;
slotSamples = round(opts.SlotPeriodMs * 1e-3 * fs);
slotsPerFrame = numel(opts.TDDPattern);
frameSamples = slotSamples * slotsPerFrame;
dl = complex(zeros(numSamples,1));

numFrames = ceil(numSamples / frameSamples);
for frame = 0:numFrames-1
    frameStart = frame * frameSamples;
    for s = 1:slotsPerFrame
        slotStart = frameStart + (s-1)*slotSamples + 1;
        slotEnd = min(slotStart + slotSamples - 1, numSamples);
        if slotStart > numSamples
            break;
        end

        slotType = opts.TDDPattern(s);
        activeRange = [];
        if slotType == "D"
            activeRange = slotStart:slotEnd;
        elseif slotType == "S"
            numDLSym = opts.SpecialSlotSymbols(1);
            activeEnd = min(slotStart + floor(slotSamples * numDLSym/14) - 1, slotEnd);
            activeRange = slotStart:activeEnd;
        elseif opts.AddWrongUplinkDLEnergy && slotType == "U"
            activeRange = slotStart:slotEnd;
        end

        if opts.AddWrongSpecialSlotEnergy && slotType == "S"
            activeRange = slotStart:slotEnd;
        end

        if opts.SuppressOneDownlinkSlot && frame == 0 && s == 2
            activeRange = [];
        end

        if ~isempty(activeRange)
            dl(activeRange) = dl(activeRange) + sqrt(opts.DownlinkNoisePower/2) * ...
                (randn(numel(activeRange),1) + 1i*randn(numel(activeRange),1));
        end
    end
end
end

function ssb = buildSyntheticSSB(pci, opts)
meta = struct( ...
    "SampleRateHz", opts.SampleRateHz, ...
    "CenterFrequencyHz", opts.CenterFrequencyHz, ...
    "BandwidthHz", opts.BandwidthHz);
refOpts = struct( ...
    "SCSkHz", opts.SCSkHz, ...
    "NSizeGrid", opts.NSizeGrid, ...
    "SSBStartSymbol", opts.SSBStartSymbol, ...
    "IncludePSS", true, ...
    "IncludeSSS", true, ...
    "IncludePBCHDMRS", opts.IncludePBCHDMRS, ...
    "SSBIndex", opts.SSBIndex, ...
    "RemoveMean", true, ...
    "NormalizeReferencePower", true);
ssb = buildSSBSyncWaveform(pci, meta, refOpts);
end

function y = addBurstAtSample(y, burst, startSample0)
% startSample0 is zero-based and can be fractional.
integerStart0 = floor(startSample0);
frac = startSample0 - integerStart0;
burst = applyFractionalDelay(burst, frac);

startIdx = integerStart0 + 1;
endIdx = startIdx + numel(burst) - 1;

if endIdx < 1 || startIdx > numel(y)
    return;
end

srcStart = max(1, 2 - startIdx);
dstStart = max(1, startIdx);
dstEnd = min(numel(y), endIdx);
srcEnd = srcStart + (dstEnd - dstStart);

y(dstStart:dstEnd) = y(dstStart:dstEnd) + burst(srcStart:srcEnd);
end

function y = applyFractionalDelay(x, delaySamples)
if abs(delaySamples) < 1e-12
    y = x;
    return;
end

halfLen = 32;
n = (-halfLen:halfLen).';
h = sinc(n - delaySamples);
window = 0.5 - 0.5*cos(2*pi*(0:numel(n)-1).'/(numel(n)-1));
h = h .* window;
h = h ./ sum(h);
y = conv(x, h, "same");
end

function y = applyMultipathChannel(x, opts)
delays = opts.MultipathDelaysSamples(:);
gainsDb = opts.MultipathGainsdB(:);
phases = opts.MultipathPhasesRad(:);
numPaths = max([numel(delays), numel(gainsDb), numel(phases)]);
delays = expandVector(delays, numPaths, 0);
gainsDb = expandVector(gainsDb, numPaths, -Inf);
phases = expandVector(phases, numPaths, 0);

if numPaths == 1 && abs(delays(1)) < 1e-12 && gainsDb(1) == 0 && abs(phases(1)) < 1e-12
    y = x;
    return;
end

y = complex(zeros(size(x)));
for p = 1:numPaths
    if ~isfinite(gainsDb(p))
        continue;
    end
    delayed = applyWholeSignalFractionalDelay(x, delays(p));
    gain = db2powLocal(gainsDb(p)/2) * exp(1i*phases(p));
    y = y + gain * delayed;
end
end

function y = addToneInterference(x, opts)
y = x;
if ~isfinite(opts.ToneInterferencePowerdB)
    return;
end

fs = opts.SampleRateHz;
n = (0:numel(x)-1).';
signalPower = mean(abs(x).^2) + eps;
tonePower = signalPower * db2powLocal(opts.ToneInterferencePowerdB);
tone = sqrt(tonePower) * exp(1i*2*pi*opts.ToneInterferenceOffsetHz*n/fs);
y = y + tone;
end

function y = applyReceiverImpairments(x, opts)
y = x;
fs = opts.SampleRateHz; %#ok<NASGU>

if opts.PhaseNoiseStepStdDeg > 0
    phase = cumsum(deg2rad(opts.PhaseNoiseStepStdDeg) * randn(numel(y), 1));
    y = y .* exp(1i*phase);
end

if opts.IQGainImbalanceDB ~= 0 || opts.IQPhaseImbalanceDeg ~= 0
    gainI = db2powLocal(opts.IQGainImbalanceDB/2);
    gainQ = db2powLocal(-opts.IQGainImbalanceDB/2);
    phaseErr = deg2rad(opts.IQPhaseImbalanceDeg);
    iPart = real(y) * gainI;
    qPart = imag(y) * gainQ;
    y = iPart + 1i * qPart .* exp(1i*phaseErr);
end

if opts.DCOffset ~= 0
    y = y + opts.DCOffset;
end

if isfinite(opts.ClippingLevelRMS)
    rmsVal = sqrt(mean(abs(y).^2) + eps);
    clipLevel = opts.ClippingLevelRMS * rmsVal;
    mag = abs(y);
    idx = mag > clipLevel;
    y(idx) = clipLevel * y(idx) ./ (mag(idx) + eps);
end
end

function y = applyWholeSignalFractionalDelay(x, delaySamples)
integerDelay = floor(delaySamples);
fracDelay = delaySamples - integerDelay;
shifted = complex(zeros(size(x)));
srcStart = max(1, 1 - integerDelay);
dstStart = max(1, 1 + integerDelay);
num = min(numel(x) - srcStart + 1, numel(x) - dstStart + 1);
if num > 0
    shifted(dstStart:dstStart+num-1) = x(srcStart:srcStart+num-1);
end
y = applyFractionalDelay(shifted, fracDelay);
end

function v = expandVector(v, n, fill)
if isempty(v)
    v = fill * ones(n, 1);
elseif numel(v) < n
    v = [v(:); repmat(v(end), n - numel(v), 1)];
else
    v = v(1:n);
end
end

function truth = buildTruthTable(scenarioName, opts, componentPowers, noisePower)
numGNBs = opts.NumGNBs;
medianOffset = median(opts.FrameOffsetsNs);
relativeOffsets = opts.FrameOffsetsNs(:) - medianOffset;
uncertaintyNs = opts.LocationUncertaintyM(:) / 299792458 * 1e9;
snrDb = pow2dbLocal(componentPowers ./ noisePower);
ssbIndexOffsetSamples = caseCSSBIndexOffsetSamplesLocal(opts);

expectedTDD = repmat("PASS", numGNBs, 1);
expectedSpecial = repmat("PASS", numGNBs, 1);
if opts.AddWrongUplinkDLEnergy
    expectedTDD(:) = "FAIL";
end
if opts.AddWrongSpecialSlotEnergy
    expectedSpecial(:) = "FAIL";
end

truth = table( ...
    (1:numGNBs).', ...
    opts.PCIs(:), ...
    repmat(opts.SSBIndex, numGNBs, 1), ...
    repmat(ssbIndexOffsetSamples, numGNBs, 1), ...
    opts.FrameOffsetsNs(:), ...
    relativeOffsets, ...
    opts.CFOHz(:), ...
    opts.SiteDistanceM(:), ...
    opts.LocationUncertaintyM(:), ...
    uncertaintyNs, ...
    snrDb(:), ...
    expectedTDD, ...
    expectedSpecial, ...
    repmat(string(scenarioName), numGNBs, 1), ...
    'VariableNames', ["GNBIndex","PCI","SSBIndex","SSBStartOffsetSamples", ...
    "FrameOffsetNs","ExpectedRelativeOffsetNs", ...
    "CFOHz","SiteDistanceM","LocationUncertaintyM","TimingUncertaintyNs", ...
    "ApproxSNRdB","ExpectedTDDStatus","ExpectedSpecialSlotStatus","ScenarioName"]);
end

function opts = defaultSyntheticOptions()
% Defaults for synthetic IQ generation. These define the clean baseline plus
% optional RF impairments used by validation and stress tests.
opts = struct();
opts.SampleRateHz = 30.72e6;
opts.CenterFrequencyHz = 3.77e9;
opts.BandwidthHz = 15e6;
opts.DurationMs = 40;
opts.FramePeriodMs = 10;
opts.SlotPeriodMs = 0.5;
opts.SCSkHz = 30;
opts.NSizeGrid = 20;
opts.SSBStartSymbol = 2;
opts.NumGNBs = 5;
opts.PCIs = [11 104 257 503 777];
opts.SSBStartOffsetUs = 0;
opts.SSBPerFrame = 1;
opts.TDDPattern = ["D","D","D","S","U","D","D","D","S","U", ...
                   "D","D","D","S","U","D","D","D","S","U"];
opts.SpecialSlotSymbols = [10 2 2];
opts.DownlinkNoisePower = 1;
opts.SSBPower = 4;
opts.IncludePBCHDMRS = true;
opts.SSBIndex = 0;
opts.CFOHz = zeros(1, opts.NumGNBs);
opts.FrameOffsetsNs = zeros(1, opts.NumGNBs);
opts.GNBPowerdB = zeros(1, opts.NumGNBs);
opts.SNRdB = 25;
opts.RandomSeed = 42;
opts.LocationUncertaintyM = 30 * ones(1, opts.NumGNBs);
opts.SiteDistanceM = 100 * ones(1, opts.NumGNBs);
opts.AddWrongUplinkDLEnergy = false;
opts.AddWrongSpecialSlotEnergy = false;
opts.SuppressOneDownlinkSlot = false;
opts.MultipathDelaysSamples = 0;
opts.MultipathGainsdB = 0;
opts.MultipathPhasesRad = 0;
opts.IQGainImbalanceDB = 0;
opts.IQPhaseImbalanceDeg = 0;
opts.DCOffset = 0;
opts.PhaseNoiseStepStdDeg = 0;
opts.ToneInterferenceOffsetHz = 0;
opts.ToneInterferencePowerdB = -Inf;
opts.ClippingLevelRMS = Inf;
end

function out = mergeStructs(base, override)
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
