function result = runRepeatedPositionSurvey(manifest,opts)
%runRepeatedPositionSurvey Analyze stationary repeats, then fit survey timing.
%
% Each physical receiver position contains repeated B210 captures. Several
% separated windows from each capture are analyzed on the configured SSB
% raster centre. Relative cell timing is formed inside each window, temporal
% variation is measured while the receiver is stationary, and the repeats
% are collapsed to one geometry observation per physical position.

arguments
    manifest
    opts struct = struct()
end

opts = defaultsLocal(opts);
[manifestTable,manifestBaseDir] = readManifestLocal(manifest);
manifestTable = canonicalizeManifestLocal(manifestTable,manifestBaseDir);
validateManifestLocal(manifestTable);

positionInfo = positionInfoFromManifestLocal(manifestTable,opts.Repeats);
positionInfoLocal = receiverXYLocal(positionInfo);
numCaptures = height(manifestTable);
loadRows = cell(numCaptures,1);
windowParts = {};
powerParts = {};
observationParts = {};
windowCount = 0;
powerCount = 0;
observationCount = 0;

bootstrap = struct();
if strlength(string(opts.Cache.BootstrapResultPath))>0
    bootstrap = loadBootstrapResultLocal(opts.Cache.BootstrapResultPath);
    if opts.Verbose
        fprintf("Reusing saved receiver analysis from:\n  %s\n", ...
            opts.Cache.BootstrapResultPath);
        fprintf("Matching observations will be cached; only missing observations will be analyzed.\n");
    end
end

totalObservations = numCaptures*numel(opts.Repeats.WindowStartSeconds)* ...
    numel(opts.Repeats.TargetGSCNs);
processedObservations = 0;
analysisStart = tic;

for captureIndex = 1:numCaptures
    capID = manifestTable.CaptureID(captureIndex);
    positionID = manifestTable.PositionID(captureIndex);
    repeatIndex = manifestTable.RepeatIndex(captureIndex);
    [bootstrapCaptureHit,captureTimingRows,capturePowerRows, ...
        captureObservationRows,captureLoadRow] = ...
        bootstrapCaptureLocal(bootstrap,capID,positionID,repeatIndex, ...
        opts.Repeats);
    if bootstrapCaptureHit
        if ~isempty(captureTimingRows)
            windowCount = windowCount+1;
            windowParts{windowCount,1} = captureTimingRows; %#ok<AGROW>
        end
        if ~isempty(capturePowerRows)
            powerCount = powerCount+1;
            powerParts{powerCount,1} = capturePowerRows; %#ok<AGROW>
        end
        observationCount = observationCount+1;
        observationParts{observationCount,1} = ...
            captureObservationRows; %#ok<AGROW>
        loadRows{captureIndex} = captureLoadRow;
        processedObservations = processedObservations+ ...
            height(captureObservationRows);
        if opts.Verbose
            fprintf("[%d/%d] %s: saved-result hit for all %d observation(s).\n", ...
                processedObservations,totalObservations,capID, ...
                height(captureObservationRows));
        end
        continue;
    end
    meta = readCaptureMetadata(manifestTable.MetadataPath(captureIndex));
    [iq,loadInfo] = loadIQWithOptionsLocal( ...
        manifestTable.IQPath(captureIndex),opts.Loader);
    validateSavedCaptureLocal(meta,loadInfo,manifestTable(captureIndex,:));

    raster = enumerateN78SSBRaster(meta,struct( ...
        "SCSkHz",opts.Repeats.SCSkHz, ...
        "CaptureBandwidthHz",meta.BandwidthHz));
    selectedRaster = raster(ismember(raster.GSCN,opts.Repeats.TargetGSCNs),:);
    missingGSCN = setdiff(opts.Repeats.TargetGSCNs,selectedRaster.GSCN);
    if ~isempty(missingGSCN)
        error("runRepeatedPositionSurvey:TargetGSCNOutsideCapture", ...
            "Capture %s does not fully contain target GSCN(s): %s.", ...
            capID,strjoin(string(missingGSCN.'),", "));
    end

    captureMeasurementCount = 0;
    for windowIndex = 1:numel(opts.Repeats.WindowStartSeconds)
        startSample0 = round( ...
            opts.Repeats.WindowStartSeconds(windowIndex)*meta.SampleRateHz);
        numWindowSamples = round( ...
            opts.Repeats.WindowDurationSeconds*meta.SampleRateHz);
        endSample0 = startSample0+numWindowSamples-1;
        if startSample0<0 || endSample0>=numel(iq)
            error("runRepeatedPositionSurvey:WindowOutsideCapture", ...
                "Window %d for capture %s lies outside the saved IQ record.", ...
                windowIndex,capID);
        end
        segment = iq(startSample0+(1:numWindowSamples));

        for rasterIndex = 1:height(selectedRaster)
            rasterRow = selectedRaster(rasterIndex,:);
            observationID = capID+"_w"+compose("%02d",windowIndex)+ ...
                "_g"+string(rasterRow.GSCN);
            processedObservations = processedObservations+1;
            if opts.Verbose
                fprintf("[%d/%d] %s, window %d, GSCN %d: ", ...
                    processedObservations,totalObservations,capID, ...
                    windowIndex,rasterRow.GSCN);
                drawnow;
            end

            signature = cacheSignatureLocal( ...
                manifestTable(captureIndex,:),meta,observationID, ...
                windowIndex,rasterRow,opts);
            cachePath = cachePathLocal(opts.Cache.Directory,observationID);
            [cacheHit,timingRows,powerRows,observationRow] = ...
                loadObservationCacheLocal(cachePath,signature,opts.Cache);
            bootstrapHit = false;
            if ~cacheHit && ~isempty(fieldnames(bootstrap))
                [bootstrapHit,timingRows,powerRows,observationRow] = ...
                    bootstrapObservationLocal(bootstrap,observationID);
                if bootstrapHit
                    saveObservationCacheLocal( ...
                        cachePath,signature,timingRows,powerRows, ...
                        observationRow,opts.Cache);
                end
            end

            if cacheHit
                if opts.Verbose
                    fprintf("cache hit, %d timing row(s).\n",height(timingRows));
                end
            elseif bootstrapHit
                if opts.Verbose
                    fprintf("saved-result hit, %d timing row(s).\n", ...
                        height(timingRows));
                end
            else
                shifted = digitalShiftLocal(segment,meta.SampleRateHz, ...
                    rasterRow.DigitalShiftHz,startSample0);
                analyzeOpts = analysisOptionsLocal( ...
                    opts,observationID,positionInfoLocal,positionID);
                captureResult = analyzeCapture(shifted,meta,analyzeOpts);

                timingRows = timingMeasurementsFromCaptureResult( ...
                    captureResult,observationID);
                if ~isempty(timingRows) && height(timingRows)>0
                    timingRows = addObservationIdentityLocal( ...
                        timingRows,capID,positionID,repeatIndex, ...
                        observationID,rasterRow,startSample0, ...
                        meta.SampleRateHz,opts.Repeats.FramePeriodNs);
                end
                powerRows = powerRowsLocal(captureResult,capID,positionID, ...
                    repeatIndex,observationID,rasterRow);
                observationRow = table( ...
                    capID,positionID,repeatIndex,observationID,windowIndex, ...
                    startSample0,startSample0/meta.SampleRateHz, ...
                    rasterRow.GSCN,rasterRow.SSBCenterHz, ...
                    height(captureResult.PSSCandidates), ...
                    usableDetectionCountLocal(captureResult), ...
                    height(captureResult.CellTiming), ...
                    'VariableNames',["CaptureID","PositionID","RepeatIndex", ...
                    "ObservationID","WindowIndex","WindowStartSample0", ...
                    "WindowStartSeconds","GSCN","SSBCenterHz", ...
                    "NumPSSCandidates","NumUsableSSSDetections", ...
                    "NumCellTimingRows"]);
                saveObservationCacheLocal( ...
                    cachePath,signature,timingRows,powerRows, ...
                    observationRow,opts.Cache);
                if opts.Verbose
                    fprintf("analyzed in %.1f s, %d timing row(s).\n", ...
                        toc(analysisStart),height(timingRows));
                end
            end

            if ~isempty(timingRows) && height(timingRows)>0
                windowCount = windowCount+1;
                windowParts{windowCount,1} = timingRows; %#ok<AGROW>
                captureMeasurementCount = ...
                    captureMeasurementCount+height(timingRows);
            end
            if ~isempty(powerRows) && height(powerRows)>0
                powerCount = powerCount+1;
                powerParts{powerCount,1} = powerRows; %#ok<AGROW>
            end
            observationCount = observationCount+1;
            observationParts{observationCount,1} = observationRow; %#ok<AGROW>
        end
    end

    loadRows{captureIndex} = table( ...
        capID,positionID,repeatIndex,string(loadInfo.CapturePath), ...
        loadInfo.NumSamples,captureMeasurementCount, ...
        'VariableNames',["CaptureID","PositionID","RepeatIndex", ...
        "IQPath","NumSamples","NumWindowMeasurementRows"]);
end

windowMeasurements = vertcatOrEmptyLocal( ...
    windowParts,@emptyWindowMeasurementsLocal);
powerObservations = vertcatOrEmptyLocal( ...
    powerParts,@emptyPowerObservationsLocal);
observationSummary = vertcatOrEmptyLocal( ...
    observationParts,@emptyObservationSummaryLocal);

[windowMeasurements,powerObservations,ambiguousPCIs] = ...
    excludeAmbiguousPCIReuseLocal(windowMeasurements,powerObservations);

if isempty(windowMeasurements) || height(windowMeasurements)==0
    aggregate = emptyAggregateLocal(positionInfo);
else
    aggregate = aggregateRepeatedTimingMeasurements( ...
        windowMeasurements,manifestTable,opts.Repeats);
end
acceptedWindowMeasurements = aggregate.WindowMeasurements;
rejectedWindowMeasurements = aggregate.RejectedWindowMeasurements;
reliableCellSupport = aggregate.ReliableCellSupport;
referenceCandidateRanking = aggregate.ReferenceCandidateRanking;

surveyOpts = opts.Survey;
if isfinite(aggregate.ReferencePCI)
    surveyOpts.ReferencePCI = aggregate.ReferencePCI;
end
surveyResult = noLocationSurveyCheck( ...
    aggregate.AggregatedMeasurementTable,aggregate.PositionInfo,surveyOpts);
surveyResult = applyStationaryInstabilityToSurveyResult( ...
    surveyResult,aggregate.StationaryInstability,surveyOpts);

[filteredPowerObservations,rssRejectedPowerObservations] = ...
    filterPowerToReliableTimingCellsLocal( ...
    powerObservations,acceptedWindowMeasurements,reliableCellSupport,opts.RSS);
powerTable = aggregatePowerByPositionLocal( ...
    filteredPowerObservations,positionInfoLocal);
rssPlan = struct();
if opts.EnableRSSPlanning && ~isempty(powerTable) && height(powerTable)>0
    candidateStops = table();
    if isfield(opts.RSS,"CandidateStops")
        candidateStops = opts.RSS.CandidateStops;
    end
    rssPlan = estimateRssPriorsAndPlanStops( ...
        "all",powerTable,aggregate.PositionInfo,candidateStops,opts.RSS);
    surveyResult.PowerTable = powerTable;
    surveyResult.RssPriors = rssPlan.RssPriors;
    surveyResult.SuggestedReceiverStops = rssPlan.SuggestedReceiverStops;
    surveyResult.RssPlanningStatement = rssPlan.Statement;
end

result = struct();
result.Mode = "stationary_repeat_multi_position_survey";
result.Statement = ...
    "Repeated fixed-position captures measured temporal relative timing instability; repeat groups were then collapsed to one timing observation per receiver position for the no-location survey fit.";
result.RepeatManifest = manifestTable;
result.LoadSummary = vertcat(loadRows{:});
result.ObservationSummary = observationSummary;
result.RawWindowMeasurements = windowMeasurements;
result.WindowMeasurements = acceptedWindowMeasurements;
result.RejectedWindowMeasurements = rejectedWindowMeasurements;
result.ReliableCellSupport = reliableCellSupport;
result.ReferenceCandidateRanking = referenceCandidateRanking;
result.RelativeObservations = aggregate.RelativeObservations;
result.PerPositionTiming = aggregate.PerPositionTiming;
result.StationaryInstability = aggregate.StationaryInstability;
result.AggregatedMeasurementTable = aggregate.AggregatedMeasurementTable;
result.PositionInfo = aggregate.PositionInfo;
result.PositionInfoLocal = positionInfoLocal;
result.RawPowerObservations = powerObservations;
result.PowerObservations = filteredPowerObservations;
result.RejectedPowerObservations = rssRejectedPowerObservations;
result.PowerTable = powerTable;
result.ExcludedAmbiguousRasterPCIs = ambiguousPCIs;
result.ReferenceCellKey = aggregate.ReferenceCellKey;
result.ReferencePCI = aggregate.ReferencePCI;
result.ReferenceGSCN = aggregate.ReferenceGSCN;
result.SurveyResult = surveyResult;
if isfield(rssPlan,"RssPriors")
    result.RssPriors = rssPlan.RssPriors;
    result.SuggestedReceiverStops = rssPlan.SuggestedReceiverStops;
end

if opts.WriteReport
    result.ReportPaths = writeReportsLocal(result,opts);
end
end

function opts = defaultsLocal(opts)
if ~isfield(opts,"Analyze"); opts.Analyze = struct(); end
if ~isfield(opts,"Survey"); opts.Survey = struct(); end
if ~isfield(opts,"Loader"); opts.Loader = struct(); end
if ~isfield(opts,"RSS"); opts.RSS = struct(); end
if ~isfield(opts,"Repeats"); opts.Repeats = struct(); end
if ~isfield(opts,"Cache"); opts.Cache = struct(); end
if ~isfield(opts,"EnableRSSPlanning"); opts.EnableRSSPlanning = false; end
if ~isfield(opts,"WriteReport"); opts.WriteReport = false; end
if ~isfield(opts,"OutputDir"); opts.OutputDir = ""; end
if ~isfield(opts,"Verbose"); opts.Verbose = true; end
if ~isfield(opts,"ReportBaseName")
    opts.ReportBaseName = "stationary_repeat_survey";
end

repeatDefaults = struct( ...
    "TargetGSCNs",7987, ...
    "SCSkHz",30, ...
    "WindowStartSeconds",[0.05 0.31 0.57 0.83], ...
    "WindowDurationSeconds",0.06, ...
    "FramePeriodNs",10e6, ...
    "ReferencePCI",NaN, ...
    "ReferenceGSCN",7987, ...
    "MinRepeatsPerPosition",3, ...
    "MinObservationsPerPositionCell",4, ...
    "MinRepeatsPerPositionCell",2, ...
    "TimingUncertaintyFloorNs",10, ...
    "DefaultPositionUncertaintyM",5, ...
    "MaxTimingStdNs",1000, ...
    "MinMeanSSSMetric",0.06, ...
    "MinWindowsPerCaptureCell",2, ...
    "MinCapturesPerPositionCell",2);
opts.Repeats = mergeStructsLocal(repeatDefaults,opts.Repeats);

analyzeDefaults = struct( ...
    "PSS",struct("MaxCandidatesPerNID2",16,"MaxTotalCandidates",48), ...
    "SSS",struct("MaxCandidatesToProcess",48,"MaxAbsCFOHz",15e3), ...
    "SIC",struct("Enable",false), ...
    "Timing",struct("EnableTimingRefinement",false));
opts.Analyze = mergeStructsLocal(analyzeDefaults,opts.Analyze);

rssDefaults = struct( ...
    "MinReliablePositions",3, ...
    "MinReliableCaptures",6);
opts.RSS = mergeStructsLocal(rssDefaults,opts.RSS);

if strlength(string(opts.OutputDir))>0
    defaultCacheDir = fullfile(string(opts.OutputDir),"analysis_cache");
else
    defaultCacheDir = fullfile(string(pwd),"reports","analysis_cache");
end
cacheDefaults = struct( ...
    "Enable",true, ...
    "Rebuild",false, ...
    "Directory",defaultCacheDir, ...
    "BootstrapResultPath","");
opts.Cache = mergeStructsLocal(cacheDefaults,opts.Cache);
end

function analyzeOpts = analysisOptionsLocal(opts,observationID, ...
        positionInfoLocal,positionID)
analyzeOpts = opts.Analyze;
for name = ["SSBRasterSearch","SSBOffsetSearch","TrimStartSeconds", ...
        "AnalysisTrimStartSeconds","AnalysisDurationSeconds", ...
        "MaxAnalysisSeconds"]
    if isfield(analyzeOpts,name)
        analyzeOpts = rmfield(analyzeOpts,name);
    end
end
analyzeOpts.EnablePowerMeasurements = opts.EnableRSSPlanning;
analyzeOpts.CaptureID = observationID;
idx = find(positionInfoLocal.CaptureID==positionID,1);
if ~isempty(idx)
    analyzeOpts.RxX_m = positionInfoLocal.RxX_m(idx);
    analyzeOpts.RxY_m = positionInfoLocal.RxY_m(idx);
    analyzeOpts.RxZ_m = positionInfoLocal.RxZ_m(idx);
end
if isfield(analyzeOpts,"Power")
    analyzeOpts.Power = mergeStructsLocal(opts.RSS,analyzeOpts.Power);
else
    analyzeOpts.Power = opts.RSS;
end
end

function rows = addObservationIdentityLocal(rows,capID,positionID, ...
        repeatIndex,observationID,rasterRow,startSample0,fs,framePeriodNs)
num = height(rows);
originNs = startSample0/fs*1e9;
rows.FramePhaseNs = mod(rows.FramePhaseNs+originNs,framePeriodNs);
rows.CaptureID = repmat(capID,num,1);
rows.PositionID = repmat(positionID,num,1);
rows.RepeatIndex = repmat(repeatIndex,num,1);
rows.ObservationID = repmat(observationID,num,1);
rows.GSCN = repmat(rasterRow.GSCN,num,1);
rows.CellKey = "GSCN"+string(rows.GSCN)+"_PCI"+string(rows.PCI);
rows.SelectedSSBOffsetHz = repmat(rasterRow.DigitalShiftHz,num,1);
rows.HypothesizedSSBCenterHz = repmat(rasterRow.SSBCenterHz,num,1);
front = ["CaptureID","PositionID","RepeatIndex","ObservationID", ...
    "CellKey","GSCN","PCI"];
rows = movevars(rows,front,"Before",1);
end

function rows = powerRowsLocal(captureResult,capID,positionID, ...
        repeatIndex,observationID,rasterRow)
if ~isfield(captureResult,"PowerMeasurements") || ...
        isempty(captureResult.PowerMeasurements) || ...
        height(captureResult.PowerMeasurements)==0
    rows = emptyPowerObservationsLocal();
    return;
end
rows = captureResult.PowerMeasurements;
num = height(rows);
rows.CaptureID = repmat(capID,num,1);
rows.PositionID = repmat(positionID,num,1);
rows.RepeatIndex = repmat(repeatIndex,num,1);
rows.ObservationID = repmat(observationID,num,1);
rows.GSCN = repmat(rasterRow.GSCN,num,1);
rows.SSBCenterHz = repmat(rasterRow.SSBCenterHz,num,1);
rows.DigitalShiftHz = repmat(rasterRow.DigitalShiftHz,num,1);
rows.CellKey = "GSCN"+string(rows.GSCN)+"_PCI"+string(rows.PCI);
rows = rows(:,["CaptureID","PositionID","RepeatIndex","ObservationID", ...
    "CellKey","GSCN","SSBCenterHz","DigitalShiftHz","PCI", ...
    "PowerDbFS","NoiseDbFS","SNRdB","NumDetections","QualityFlag"]);
end

function count = usableDetectionCountLocal(captureResult)
count = 0;
if ~isfield(captureResult,"SSSDetections") || ...
        isempty(captureResult.SSSDetections)
    return;
end
det = captureResult.SSSDetections;
if ismember("IsUsable",string(det.Properties.VariableNames))
    count = nnz(det.IsUsable);
else
    count = height(det);
end
end

function shifted = digitalShiftLocal(iq,fs,offsetHz,originSample0)
n = originSample0+(0:numel(iq)-1).';
rotation = exp(-1i*2*pi*offsetHz/fs.*n);
if isa(iq,"single")
    rotation = single(rotation);
end
shifted = iq(:).*rotation;
end

function power = aggregatePowerByPositionLocal(observations,positionInfoLocal)
if isempty(observations) || height(observations)==0
    power = emptyPowerTableLocal();
    return;
end

keys = unique(observations(:,["PositionID","CellKey","GSCN", ...
    "SSBCenterHz","PCI"]),"rows","stable");
parts = cell(height(keys),1);
for k = 1:height(keys)
    idx = observations.PositionID==keys.PositionID(k) & ...
        observations.CellKey==keys.CellKey(k);
    rows = observations(idx,:);
    powerDb = linearMeanDbLocal(rows.PowerDbFS);
    noiseDb = linearMeanDbLocal(rows.NoiseDbFS);
    snrDb = powerDb-noiseDb;
    positionIndex = find( ...
        positionInfoLocal.CaptureID==keys.PositionID(k),1);
    if isempty(positionIndex)
        xyz = [NaN NaN NaN];
    else
        xyz = [positionInfoLocal.RxX_m(positionIndex), ...
            positionInfoLocal.RxY_m(positionIndex), ...
            positionInfoLocal.RxZ_m(positionIndex)];
    end
    if any(rows.QualityFlag=="SUSPECT")
        quality = "SUSPECT";
    else
        quality = "OK";
    end
    parts{k} = table( ...
        keys.PositionID(k),keys.CellKey(k),keys.GSCN(k), ...
        keys.SSBCenterHz(k),0,keys.PCI(k),xyz(1),xyz(2),xyz(3), ...
        powerDb,noiseDb,snrDb,sum(rows.NumDetections),quality, ...
        'VariableNames',["CaptureID","CellKey","GSCN","SSBCenterHz", ...
        "DigitalShiftHz","PCI","RxX_m","RxY_m","RxZ_m", ...
        "PowerDbFS","NoiseDbFS","SNRdB","NumDetections","QualityFlag"]);
end
power = vertcat(parts{:});
end

function [accepted,rejected] = filterPowerToReliableTimingCellsLocal( ...
        power,timing,support,rssOpts)
if isempty(power) || height(power)==0
    accepted = power;
    rejected = addPowerRejectionReasonLocal(power,strings(0,1));
    return;
end

reason = strings(height(power),1);
if isempty(timing) || height(timing)==0
    reason(:) = "NO_QUALITY_ACCEPTED_TIMING";
else
    timingKeys = unique(timing.ObservationID+"|"+timing.CellKey);
    powerKeys = power.ObservationID+"|"+power.CellKey;
    reason(~ismember(powerKeys,timingKeys)) = ...
        "NO_TIMING_CONFIRMATION_IN_OBSERVATION";
end

if isempty(support) || height(support)==0
    reason(reason=="") = "INSUFFICIENT_REPEAT_SUPPORT";
else
    reliable = support.CellKey( ...
        support.ReliablePositions>=rssOpts.MinReliablePositions & ...
        support.ReliableCaptures>=rssOpts.MinReliableCaptures);
    reason(reason=="" & ~ismember(power.CellKey,reliable)) = ...
        "INSUFFICIENT_REPEAT_SUPPORT";
end

accepted = power(reason=="",:);
rejected = addPowerRejectionReasonLocal( ...
    power(reason~="",:),reason(reason~=""));
end

function rows = addPowerRejectionReasonLocal(rows,reason)
rows.RejectionReason = reason;
end

function value = linearMeanDbLocal(values)
value = 10*log10(mean(10.^(double(values)/10),"omitnan")+eps);
end

function [measurements,power,ambiguousPCIs] = ...
        excludeAmbiguousPCIReuseLocal(measurements,power)
ambiguousPCIs = zeros(0,1);
if isempty(measurements) || height(measurements)==0
    return;
end
pcis = unique(measurements.PCI);
for k = 1:numel(pcis)
    if numel(unique(measurements.CellKey(measurements.PCI==pcis(k))))>1
        ambiguousPCIs(end+1,1) = pcis(k); %#ok<AGROW>
    end
end
if isempty(ambiguousPCIs)
    return;
end
warning("runRepeatedPositionSurvey:PCIReuseAcrossRasterCenters", ...
    "Excluded PCI(s) %s because the geometry estimator is PCI-indexed and cannot merge distinct GSCN/PCI cell identities.", ...
    strjoin(string(ambiguousPCIs.'),", "));
measurements = measurements(~ismember(measurements.PCI,ambiguousPCIs),:);
if ~isempty(power) && height(power)>0
    power = power(~ismember(power.PCI,ambiguousPCIs),:);
end
end

function paths = writeReportsLocal(result,opts)
if strlength(string(opts.OutputDir))==0
    error("runRepeatedPositionSurvey:MissingOutputDir", ...
        "OutputDir must be set when WriteReport is true.");
end
outputDir = string(opts.OutputDir);
if ~isfolder(outputDir); mkdir(outputDir); end
base = string(opts.ReportBaseName);
paths = generateSurveyReport(result.SurveyResult,outputDir,base);

extra = struct();
extra.WindowMeasurementsCSV = fullfile(outputDir,base+"_window_timing.csv");
extra.RawWindowMeasurementsCSV = fullfile(outputDir,base+"_raw_window_timing.csv");
extra.RejectedWindowMeasurementsCSV = fullfile(outputDir,base+"_rejected_window_timing.csv");
extra.ReliableCellSupportCSV = fullfile(outputDir,base+"_reliable_cell_support.csv");
extra.ReferenceCandidateRankingCSV = fullfile(outputDir,base+"_reference_ranking.csv");
extra.RelativeObservationsCSV = fullfile(outputDir,base+"_stationary_relative_observations.csv");
extra.PerPositionTimingCSV = fullfile(outputDir,base+"_per_position_timing.csv");
extra.StationaryInstabilityCSV = fullfile(outputDir,base+"_stationary_instability.csv");
extra.AggregatedGeometryCSV = fullfile(outputDir,base+"_geometry_measurements.csv");
extra.PositionInfoCSV = fullfile(outputDir,base+"_positions.csv");
extra.ObservationSummaryCSV = fullfile(outputDir,base+"_observation_summary.csv");
extra.RawPowerObservationsCSV = fullfile(outputDir,base+"_raw_power_observations.csv");
extra.FilteredPowerObservationsCSV = fullfile(outputDir,base+"_filtered_power_observations.csv");
extra.RejectedPowerObservationsCSV = fullfile(outputDir,base+"_rejected_power_observations.csv");
extra.RepeatedSurveyMAT = fullfile(outputDir,base+"_repeated_survey.mat");
writetable(result.WindowMeasurements,extra.WindowMeasurementsCSV);
writetable(result.RawWindowMeasurements,extra.RawWindowMeasurementsCSV);
writetable(result.RejectedWindowMeasurements,extra.RejectedWindowMeasurementsCSV);
writetable(result.ReliableCellSupport,extra.ReliableCellSupportCSV);
writetable(result.ReferenceCandidateRanking,extra.ReferenceCandidateRankingCSV);
writetable(result.RelativeObservations,extra.RelativeObservationsCSV);
writetable(result.PerPositionTiming,extra.PerPositionTimingCSV);
writetable(result.StationaryInstability,extra.StationaryInstabilityCSV);
writetable(result.AggregatedMeasurementTable,extra.AggregatedGeometryCSV);
writetable(result.PositionInfo,extra.PositionInfoCSV);
writetable(result.ObservationSummary,extra.ObservationSummaryCSV);
writetable(result.RawPowerObservations,extra.RawPowerObservationsCSV);
writetable(result.PowerObservations,extra.FilteredPowerObservationsCSV);
writetable(result.RejectedPowerObservations,extra.RejectedPowerObservationsCSV);
save(extra.RepeatedSurveyMAT,"result","-v7.3");
names = fieldnames(extra);
for k = 1:numel(names)
    paths.(names{k}) = extra.(names{k});
end
end

function [tbl,baseDir] = readManifestLocal(manifest)
if istable(manifest)
    tbl = manifest;
    baseDir = string(pwd);
elseif ischar(manifest) || isstring(manifest)
    manifestPath = string(manifest);
    if ~isfile(manifestPath)
        error("runRepeatedPositionSurvey:ManifestNotFound", ...
            "Repeat manifest file not found: %s",manifestPath);
    end
    tbl = readtable(manifestPath,"TextType","string", ...
        "VariableNamingRule","preserve","Delimiter",",", ...
        "ReadVariableNames",true);
    baseDir = string(fileparts(manifestPath));
else
    error("runRepeatedPositionSurvey:InvalidManifest", ...
        "Manifest must be a table or CSV path.");
end
end

function tbl = canonicalizeManifestLocal(tbl,baseDir)
required = ["CaptureID","PositionID","RepeatIndex","IQPath", ...
    "MetadataPath","RxLat","RxLon","RxAltM","RxPositionUncertaintyM"];
missing = required(~ismember(required,string(tbl.Properties.VariableNames)));
if ~isempty(missing)
    error("runRepeatedPositionSurvey:MissingManifestColumns", ...
        "Repeat manifest is missing: %s",strjoin(missing,", "));
end
tbl = tbl(:,required);
tbl.CaptureID = string(tbl.CaptureID);
tbl.PositionID = string(tbl.PositionID);
tbl.RepeatIndex = double(tbl.RepeatIndex);
tbl.IQPath = string(tbl.IQPath);
tbl.MetadataPath = string(tbl.MetadataPath);
for k = 1:height(tbl)
    tbl.IQPath(k) = resolvePathLocal(tbl.IQPath(k),baseDir);
    tbl.MetadataPath(k) = resolvePathLocal(tbl.MetadataPath(k),baseDir);
end
end

function validateManifestLocal(tbl)
if isempty(tbl) || height(tbl)==0
    error("runRepeatedPositionSurvey:EmptyManifest", ...
        "The repeated-position manifest has no accepted captures.");
end
if numel(unique(tbl.CaptureID))~=height(tbl)
    error("runRepeatedPositionSurvey:DuplicateCaptureID", ...
        "CaptureID values must be unique.");
end
keys = tbl.PositionID+"_"+string(tbl.RepeatIndex);
if numel(unique(keys))~=height(tbl)
    error("runRepeatedPositionSurvey:DuplicateRepeatIndex", ...
        "Each PositionID/RepeatIndex pair must be unique.");
end
if any(~isfinite(tbl.RxLat)) || any(~isfinite(tbl.RxLon)) || ...
        any(~isfinite(tbl.RxAltM)) || ...
        any(~isfinite(tbl.RxPositionUncertaintyM))
    error("runRepeatedPositionSurvey:InvalidReceiverPosition", ...
        "Every repeat manifest row requires finite GPS position fields.");
end
end

function info = positionInfoFromManifestLocal(manifest,repeatOpts)
positions = unique(manifest.PositionID,"stable");
parts = cell(numel(positions),1);
for k = 1:numel(positions)
    rows = manifest(manifest.PositionID==positions(k),:);
    repeats = unique(rows.RepeatIndex);
    if numel(repeats)<repeatOpts.MinRepeatsPerPosition
        error("runRepeatedPositionSurvey:TooFewRepeats", ...
            "Position %s has %d repeat(s); at least %d are required.", ...
            positions(k),numel(repeats),repeatOpts.MinRepeatsPerPosition);
    end
    lat = median(double(rows.RxLat),"omitnan");
    lon = median(double(rows.RxLon),"omitnan");
    alt = median(double(rows.RxAltM),"omitnan");
    [dx,dy] = localOffsetsLocal( ...
        double(rows.RxLat),double(rows.RxLon),lat,lon);
    spread = max(hypot(dx,dy),[],"omitnan");
    uncertainty = max([repeatOpts.DefaultPositionUncertaintyM, ...
        median(double(rows.RxPositionUncertaintyM),"omitnan"),spread]);
    parts{k} = table(positions(k),numel(repeats),lat,lon,alt, ...
        uncertainty,spread, ...
        'VariableNames',["CaptureID","NumRepeats","RxLat","RxLon", ...
        "RxAltM","RxPositionUncertaintyM","RepeatGpsSpreadM"]);
end
info = vertcat(parts{:});
end

function xy = receiverXYLocal(info)
lat = double(info.RxLat);
lon = double(info.RxLon);
earthRadiusM = 6371000;
x = deg2rad(lon-lon(1))*earthRadiusM*cos(deg2rad(lat(1)));
y = deg2rad(lat-lat(1))*earthRadiusM;
z = double(info.RxAltM)-double(info.RxAltM(1));
xy = table(string(info.CaptureID),x(:),y(:),z(:), ...
    double(info.RxPositionUncertaintyM), ...
    'VariableNames',["CaptureID","RxX_m","RxY_m","RxZ_m", ...
    "RxPositionUncertaintyM"]);
end

function validateSavedCaptureLocal(meta,loadInfo,manifestRow)
capID = manifestRow.CaptureID;
if ~isfield(meta,"Raw") || ~isstruct(meta.Raw)
    error("runRepeatedPositionSurvey:MissingValidationMetadata", ...
        "Capture %s lacks the validation metadata required for the real survey.", ...
        capID);
end
raw = meta.Raw;
required = ["capture_validated","overrun_count","actual_sample_count", ...
    "expected_sample_count","capture_id","position_id","repeat_index"];
missing = required(~isfield(raw,cellstr(required)));
if ~isempty(missing)
    error("runRepeatedPositionSurvey:MissingValidationMetadata", ...
        "Capture %s metadata is missing: %s.",capID,strjoin(missing,", "));
end
if ~logical(raw.capture_validated) || double(raw.overrun_count)~=0
    error("runRepeatedPositionSurvey:InvalidCapture", ...
        "Capture %s is not a validated zero-overrun record.",capID);
end
if double(raw.actual_sample_count)~=double(loadInfo.NumSamples) || ...
        double(raw.expected_sample_count)~=double(loadInfo.NumSamples)
    error("runRepeatedPositionSurvey:SampleCountMismatch", ...
        "Capture %s sample count does not match its validated metadata.",capID);
end
if string(raw.capture_id)~=capID || ...
        string(raw.position_id)~=manifestRow.PositionID || ...
        double(raw.repeat_index)~=manifestRow.RepeatIndex
    error("runRepeatedPositionSurvey:MetadataIdentityMismatch", ...
        "Capture %s metadata does not match the repeat manifest identity.",capID);
end
end

function [iq,info] = loadIQWithOptionsLocal(path,loaderOpts)
if isempty(fieldnames(loaderOpts))
    [iq,info] = loadIQCapture(path);
    return;
end
names = fieldnames(loaderOpts);
args = cell(1,2*numel(names));
for k = 1:numel(names)
    args{2*k-1} = names{k};
    args{2*k} = loaderOpts.(names{k});
end
[iq,info] = loadIQCapture(path,args{:});
end

function result = loadBootstrapResultLocal(path)
path = string(path);
if ~isfile(path)
    error("runRepeatedPositionSurvey:BootstrapResultNotFound", ...
        "Saved receiver-analysis result not found: %s",path);
end
loaded = load(path,"result");
if ~isfield(loaded,"result")
    error("runRepeatedPositionSurvey:InvalidBootstrapResult", ...
        "Saved MAT file does not contain a result structure.");
end
result = loaded.result;
if isfield(result,"RawWindowMeasurements")
    result.WindowMeasurements = result.RawWindowMeasurements;
end
if isfield(result,"RawPowerObservations")
    result.PowerObservations = result.RawPowerObservations;
end
required = ["WindowMeasurements","PowerObservations", ...
    "ObservationSummary"];
missing = required(~isfield(result,cellstr(required)));
if ~isempty(missing)
    error("runRepeatedPositionSurvey:InvalidBootstrapResult", ...
        "Saved result is missing: %s",strjoin(missing,", "));
end
end

function [hit,timingRows,powerRows,observationRow] = ...
        bootstrapObservationLocal(result,observationID)
timingRows = result.WindowMeasurements( ...
    result.WindowMeasurements.ObservationID==observationID,:);
powerRows = result.PowerObservations( ...
    result.PowerObservations.ObservationID==observationID,:);
if ~isempty(result.ObservationSummary) && ...
        ismember("ObservationID", ...
        string(result.ObservationSummary.Properties.VariableNames))
    observationRow = result.ObservationSummary( ...
        result.ObservationSummary.ObservationID==observationID,:);
else
    observationRow = emptyObservationSummaryLocal();
end
hit = ~isempty(timingRows) || ~isempty(powerRows) || ...
    ~isempty(observationRow);
if hit && isempty(observationRow)
    source = timingRows;
    if isempty(source)
        source = powerRows;
    end
    tokens = regexp(observationID,"_w(\d+)_g(\d+)$","tokens","once");
    if isempty(tokens)
        windowIndex = 0;
        gscn = source.GSCN(1);
    else
        windowIndex = str2double(tokens{1});
        gscn = str2double(tokens{2});
    end
    observationRow = table( ...
        source.CaptureID(1),source.PositionID(1),source.RepeatIndex(1), ...
        observationID,windowIndex,NaN,NaN,gscn,NaN,NaN,NaN, ...
        height(timingRows), ...
        'VariableNames',["CaptureID","PositionID","RepeatIndex", ...
        "ObservationID","WindowIndex","WindowStartSample0", ...
        "WindowStartSeconds","GSCN","SSBCenterHz","NumPSSCandidates", ...
        "NumUsableSSSDetections","NumCellTimingRows"]);
end
end

function [hit,timingRows,powerRows,observationRows,loadRow] = ...
        bootstrapCaptureLocal(result,capID,positionID,repeatIndex,repeatOpts)
hit = false;
timingRows = table();
powerRows = table();
observationRows = table();
loadRow = table();
if isempty(fieldnames(result))
    return;
end

expected = strings(0,1);
for windowIndex = 1:numel(repeatOpts.WindowStartSeconds)
    for gscn = reshape(repeatOpts.TargetGSCNs,1,[])
        expected(end+1,1) = capID+"_w"+compose("%02d",windowIndex)+ ...
            "_g"+string(gscn); %#ok<AGROW>
    end
end
available = unique([ ...
    string(result.WindowMeasurements.ObservationID); ...
    string(result.PowerObservations.ObservationID)]);
if ~all(ismember(expected,available))
    return;
end

timingRows = result.WindowMeasurements( ...
    result.WindowMeasurements.CaptureID==capID,:);
powerRows = result.PowerObservations( ...
    result.PowerObservations.CaptureID==capID,:);
observationParts = cell(numel(expected),1);
for k = 1:numel(expected)
    [~,~,~,observationParts{k}] = ...
        bootstrapObservationLocal(result,expected(k));
end
observationRows = vertcat(observationParts{:});

if isfield(result,"LoadSummary") && ~isempty(result.LoadSummary) && ...
        ismember("CaptureID", ...
        string(result.LoadSummary.Properties.VariableNames))
    loadRow = result.LoadSummary(result.LoadSummary.CaptureID==capID,:);
end
if isempty(loadRow)
    loadRow = table( ...
        capID,positionID,repeatIndex,"",NaN,height(timingRows), ...
        'VariableNames',["CaptureID","PositionID","RepeatIndex", ...
        "IQPath","NumSamples","NumWindowMeasurementRows"]);
end
hit = true;
end

function signature = cacheSignatureLocal(manifestRow,meta,observationID, ...
        windowIndex,rasterRow,opts)
fileInfo = dir(manifestRow.IQPath);
if isempty(fileInfo)
    fileBytes = NaN;
    fileDate = NaN;
else
    fileBytes = fileInfo.bytes;
    fileDate = fileInfo.datenum;
end
payload = struct();
payload.Version = "repeat_survey_analysis_cache_v1";
payload.CaptureID = manifestRow.CaptureID;
payload.ObservationID = observationID;
payload.IQBytes = fileBytes;
payload.IQDateNum = fileDate;
payload.SampleRateHz = meta.SampleRateHz;
payload.WindowIndex = windowIndex;
payload.WindowStartSeconds = opts.Repeats.WindowStartSeconds(windowIndex);
payload.WindowDurationSeconds = opts.Repeats.WindowDurationSeconds;
payload.GSCN = rasterRow.GSCN;
payload.SSBCenterHz = rasterRow.SSBCenterHz;
payload.DigitalShiftHz = rasterRow.DigitalShiftHz;
payload.Analyze = opts.Analyze;
payload.EnablePowerMeasurements = opts.EnableRSSPlanning;
signature = string(jsonencode(payload));
end

function path = cachePathLocal(cacheDir,observationID)
path = fullfile(string(cacheDir),observationID+"_analysis.mat");
end

function [hit,timingRows,powerRows,observationRow] = ...
        loadObservationCacheLocal(path,signature,cacheOpts)
hit = false;
timingRows = table();
powerRows = table();
observationRow = table();
if ~cacheOpts.Enable || cacheOpts.Rebuild || ~isfile(path)
    return;
end
cached = load(path,"signature","timingRows","powerRows","observationRow");
required = ["signature","timingRows","powerRows","observationRow"];
if any(~isfield(cached,cellstr(required))) || ...
        string(cached.signature)~=signature
    return;
end
timingRows = cached.timingRows;
powerRows = cached.powerRows;
observationRow = cached.observationRow;
hit = true;
end

function saveObservationCacheLocal(path,signature,timingRows,powerRows, ...
        observationRow,cacheOpts)
if ~cacheOpts.Enable
    return;
end
cacheDir = fileparts(path);
if ~isfolder(cacheDir)
    mkdir(cacheDir);
end
save(path,"signature","timingRows","powerRows","observationRow");
end

function value = vertcatOrEmptyLocal(parts,emptyFactory)
if isempty(parts)
    value = emptyFactory();
else
    value = vertcat(parts{:});
end
end

function aggregate = emptyAggregateLocal(positionInfo)
aggregate = struct();
aggregate.ReferenceCellKey = "";
aggregate.ReferencePCI = NaN;
aggregate.ReferenceGSCN = NaN;
aggregate.RawWindowMeasurements = emptyWindowMeasurementsLocal();
aggregate.WindowMeasurements = emptyWindowMeasurementsLocal();
aggregate.RejectedWindowMeasurements = addTimingRejectionReasonLocal( ...
    emptyWindowMeasurementsLocal(),strings(0,1));
aggregate.ReliableCellSupport = emptyReliableCellSupportLocal();
aggregate.ReferenceCandidateRanking = emptyReliableCellSupportLocal();
aggregate.RelativeObservations = emptyRelativeObservationsLocal();
aggregate.PerPositionTiming = emptyPerPositionTimingLocal();
aggregate.StationaryInstability = emptyStationaryInstabilityLocal();
aggregate.AggregatedMeasurementTable = emptyGeometryMeasurementsLocal();
aggregate.PositionInfo = positionInfo;
end

function rows = addTimingRejectionReasonLocal(rows,reason)
rows.RejectionReason = reason;
end

function tbl = emptyWindowMeasurementsLocal()
tbl = table( ...
    strings(0,1),strings(0,1),zeros(0,1),strings(0,1),strings(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',["CaptureID","PositionID","RepeatIndex", ...
    "ObservationID","CellKey","GSCN","PCI","FramePhaseNs", ...
    "TimingStdNs","CenterFrequencyHz","SelectedSSBOffsetHz", ...
    "HypothesizedSSBCenterHz","MeanPSSMetric","MeanSSSMetric", ...
    "MeanCFOHz","Unused"]);
tbl.Unused = [];
end

function tbl = emptyPowerObservationsLocal()
tbl = table( ...
    strings(0,1),strings(0,1),zeros(0,1),strings(0,1),strings(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),strings(0,1), ...
    'VariableNames',["CaptureID","PositionID","RepeatIndex", ...
    "ObservationID","CellKey","GSCN","SSBCenterHz","DigitalShiftHz", ...
    "PCI","PowerDbFS","NoiseDbFS","SNRdB","NumDetections","QualityFlag"]);
end

function tbl = emptyObservationSummaryLocal()
tbl = table( ...
    strings(0,1),strings(0,1),zeros(0,1),strings(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1), ...
    'VariableNames',["CaptureID","PositionID","RepeatIndex", ...
    "ObservationID","WindowIndex","WindowStartSample0", ...
    "WindowStartSeconds","GSCN","SSBCenterHz","NumPSSCandidates", ...
    "NumUsableSSSDetections","NumCellTimingRows"]);
end

function tbl = emptyPowerTableLocal()
tbl = table( ...
    strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),strings(0,1), ...
    'VariableNames',["CaptureID","CellKey","GSCN","SSBCenterHz", ...
    "DigitalShiftHz","PCI","RxX_m","RxY_m","RxZ_m", ...
    "PowerDbFS","NoiseDbFS","SNRdB","NumDetections","QualityFlag"]);
end

function tbl = emptyReliableCellSupportLocal()
tbl = table( ...
    strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',["CellKey","GSCN","PCI","ReliablePositions", ...
    "ReliableCaptures","ReliableObservations","NumTimingRows", ...
    "MeanSSSMetric","MedianTimingStdNs"]);
end

function tbl = emptyRelativeObservationsLocal()
tbl = table();
end

function tbl = emptyPerPositionTimingLocal()
tbl = table();
end

function tbl = emptyStationaryInstabilityLocal()
tbl = table( ...
    strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),strings(0,1),strings(0,1), ...
    'VariableNames',["CellKey","GSCN","PCI","NumObservations", ...
    "NumPositions","NumCaptures","StationaryTimingInstabilityRmsNs", ...
    "StationaryTimingInstabilityPeakToPeakNs", ...
    "StationaryTimingInstabilityMaxAbsNs", ...
    "ExpectedMeasurementNoiseRmsNs","ExcessStationaryRmsNs", ...
    "StationaryAssessmentStatus","StationaryAssessmentReason"]);
end

function tbl = emptyGeometryMeasurementsLocal()
tbl = table( ...
    strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1), ...
    'VariableNames',["CaptureID","CellKey","GSCN","PCI", ...
    "FramePhaseNs","TimingStdNs","CenterFrequencyHz", ...
    "SelectedSSBOffsetHz","HypothesizedSSBCenterHz", ...
    "MeanPSSMetric","MeanSSSMetric","MeanCFOHz"]);
end

function path = resolvePathLocal(path,baseDir)
path = string(path);
if isempty(regexp(char(path),"^[A-Za-z]:[\\/]","once")) && ...
        ~startsWith(path,filesep)
    path = fullfile(baseDir,path);
end
end

function [x,y] = localOffsetsLocal(lat,lon,lat0,lon0)
earthRadiusM = 6371000;
x = deg2rad(lon-lon0).*earthRadiusM.*cos(deg2rad(lat0));
y = deg2rad(lat-lat0).*earthRadiusM;
end

function out = mergeStructsLocal(base,override)
out = base;
if isempty(override) || isempty(fieldnames(override))
    return;
end
names = fieldnames(override);
for k = 1:numel(names)
    name = names{k};
    if isfield(base,name) && isstruct(base.(name)) && ...
            isstruct(override.(name))
        out.(name) = mergeStructsLocal(base.(name),override.(name));
    else
        out.(name) = override.(name);
    end
end
end
