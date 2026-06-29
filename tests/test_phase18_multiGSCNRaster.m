% test_phase18_multiGSCNRaster
% One wideband IQ capture must recover cells on multiple legal n78 GSCNs
% without changing their original sample-time coordinates.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot,"src")));

fprintf("Running Phase 18 multi-GSCN raster tests...\n");

captureCenterHz = 3.70272e9; % GSCN 7987
fs = 30.72e6;
metaRaster = struct( ...
    "SampleRateHz",fs, ...
    "CenterFrequencyHz",captureCenterHz, ...
    "BandwidthHz",30.72e6);
raster = enumerateN78SSBRaster(metaRaster);
assert(height(raster)==17, ...
    "Expected 17 complete 30 kHz-SCS SSB centres in a 30.72 MHz capture.");
assert(any(raster.GSCN==7987 & abs(raster.DigitalShiftHz)<1), ...
    "Expected the capture-centre GSCN in the raster list.");
assert(all(raster.LowerEdgeMarginHz>=-1));
assert(all(raster.UpperEdgeMarginHz>=-1));

% Create two independent cells at different legal GSCN centres, then combine
% them into one wideband IQ capture centred at GSCN 7987.
[iqLow,metaLow,truthLow] = makeCellLocal(257,3.69984e9,0,9101);
[iqHigh,metaHigh,truthHigh] = makeCellLocal(503,3.70560e9,1000,9102);
n = (0:numel(iqLow)-1).';
wideIQ = iqLow.*exp(1i*2*pi*(metaLow.CenterFrequencyHz-captureCenterHz)/fs.*n) + ...
    0.9*iqHigh.*exp(1i*2*pi*(metaHigh.CenterFrequencyHz-captureCenterHz)/fs.*n);

metaWide = metaLow;
metaWide.CenterFrequencyHz = captureCenterHz;
metaWide.BandwidthHz = 15.36e6;

opts = struct();
opts.SSBRasterSearch = struct( ...
    "ScreenStartSeconds",0, ...
    "ScreenDurationSeconds",0.04, ...
    "FullAnalysisStartSeconds",0, ...
    "FullAnalysisDurationSeconds",0.04, ...
    "MinFullAnalysisCenters",2, ...
    "MaxFullAnalysisCenters",6);
opts.PSS = struct("MaxCandidatesPerNID2",16,"MaxTotalCandidates",48);
opts.SSS = struct("MaxCandidatesToProcess",48,"MaxAbsCFOHz",15e3);
opts.SIC = struct("Enable",false);
opts.Timing = struct("EnableTimingRefinement",false);

result = analyzeCaptureAcrossSSBRaster(wideIQ,metaWide,opts);
assert(result.SSBRasterSearch.SingleCaptureCommonTimeBase);
assert(result.SSBRasterSearch.NumRasterCenters==5);

expectedKeys = ["GSCN7985_PCI257";"GSCN7989_PCI503"];
assert(all(ismember(expectedKeys,result.CellTiming.CellKey)), ...
    "Expected both frequency-separated cells to be retained.");

lowRow = result.CellTiming(result.CellTiming.CellKey==expectedKeys(1),:);
highRow = result.CellTiming(result.CellTiming.CellKey==expectedKeys(2),:);
assert(~isempty(lowRow) && ~isempty(highRow));

directLow = analyzeCapture(iqLow,metaLow,optsWithoutRasterLocal(opts));
directHigh = analyzeCapture(iqHigh,metaHigh,optsWithoutRasterLocal(opts));
assert(abs(wrapNsLocal(lowRow.FramePhaseNs-directLow.CellTiming.FramePhaseNs,10e6))<80, ...
    "Digital GSCN search must preserve the low-frequency cell timing.");
assert(abs(wrapNsLocal(highRow.FramePhaseNs-directHigh.CellTiming.FramePhaseNs,10e6))<80, ...
    "Digital GSCN search must preserve the high-frequency cell timing.");

rows = timingMeasurementsFromCaptureResult(result,"wide_capture");
assert(all(ismember(expectedKeys,rows.CellKey)));
assert(all(ismember([7985;7989],rows.GSCN)));
assert(all(ismember([truthLow.PCI(1);truthHigh.PCI(1)],rows.PCI)));

% The manifest orchestrator must accept raster-search results and refuse a
% survey verdict when only one receiver position is available.
tmpRoot = tempname;
mkdir(tmpRoot);
cleanup = onCleanup(@() rmdir(tmpRoot,"s"));
paths = writeSyntheticCapture(tmpRoot,"wide_capture",wideIQ,metaWide, ...
    [truthLow;truthHigh]);
manifest = table("wide_capture",string(paths.IQ),string(paths.Metadata), ...
    0,0,0,5, ...
    'VariableNames',["CaptureID","IQPath","MetadataPath", ...
    "RxX_m","RxY_m","RxZ_m","RxPositionUncertaintyM"]);
surveyOpts = struct();
surveyOpts.Analyze = opts;
surveyOpts.Survey = struct("MinCaptures",3);
orchestrated = runMultiCaptureSurvey(manifest,surveyOpts);
assert(all(ismember(expectedKeys,orchestrated.MeasurementTable.CellKey)));
assert(all(orchestrated.SurveyResult.TimingEstimates.TimingStatus== ...
    "NOT_ASSESSABLE"));

fprintf("Phase 18 multi-GSCN raster tests passed.\n");

function [iq,meta,truth] = makeCellLocal(pci,centerHz,offsetNs,seed)
custom = struct();
custom.NumGNBs = 1;
custom.PCIs = pci;
custom.FrameOffsetsNs = offsetNs;
custom.CFOHz = 0;
custom.GNBPowerdB = 0;
custom.LocationUncertaintyM = 30;
custom.SiteDistanceM = 100;
custom.SNRdB = 38;
custom.DurationMs = 40;
custom.RandomSeed = seed;
custom.CenterFrequencyHz = centerHz;
custom.BandwidthHz = 15.36e6;
[iq,meta,truth] = generateSyntheticScenario("aligned_5gnb",custom);
end

function opts = optsWithoutRasterLocal(opts)
opts = rmfield(opts,"SSBRasterSearch");
end

function value = wrapNsLocal(value,period)
value = mod(value+period/2,period)-period/2;
end
