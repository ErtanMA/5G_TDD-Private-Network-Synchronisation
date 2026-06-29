function result = analyzeCaptureAcrossSSBRaster(iq, meta, opts)
%analyzeCaptureAcrossSSBRaster Detect all credible n78 SSB centres in IQ.
%
% One wideband IQ capture is searched at every legal n78 GSCN whose full
% 7.2 MHz SSB fits inside the recorded bandwidth. A short PSS-only screen is
% followed by full SSS/PCI/timing analysis at promising raster centres.
% Digital frequency translation does not move samples in time, so every
% returned cell remains on the original capture's sample-time axis.

arguments
    iq (:,1) {mustBeNumeric}
    meta struct
    opts struct = struct()
end

[meta, metaReport] = validateMetadata(meta);
searchOpts = defaultSearchOptionsLocal(meta);
if isfield(opts,"SSBRasterSearch")
    searchOpts = mergeStructsLocal(searchOpts,opts.SSBRasterSearch);
end
analysisOpts = removeSearchOptionsLocal(opts);

if ~searchOpts.Enable
    result = analyzeCapture(iq,meta,analysisOpts);
    result.SSBRasterSearch = struct( ...
        "Enable",false, ...
        "Statement","Multi-GSCN raster search disabled.");
    return;
end

rasterOpts = struct( ...
    "SCSkHz",searchOpts.SCSkHz, ...
    "CaptureBandwidthHz",searchOpts.CaptureBandwidthHz, ...
    "EdgeGuardHz",searchOpts.EdgeGuardHz);
raster = enumerateN78SSBRaster(meta,rasterOpts);
if isempty(raster)
    error("analyzeCaptureAcrossSSBRaster:NoRasterCenters", ...
        "No complete n78 SSB raster centre fits inside this IQ capture.");
end

[screenIQ, screenOriginSample0] = screenWindowLocal(iq,meta,searchOpts);
screenRows = cell(height(raster),1);
for k = 1:height(raster)
    if searchOpts.Verbose
        fprintf("SSB raster screen %d/%d: GSCN %d, %.3f MHz\n", ...
            k,height(raster),raster.GSCN(k),raster.SSBCenterHz(k)/1e6);
    end
    shifted = digitalShiftLocal( ...
        screenIQ,meta.SampleRateHz,raster.DigitalShiftHz(k),screenOriginSample0);
    [pss,debug] = estimatePSSCandidates( ...
        shifted,meta,searchOpts.ScreenPSS);
    screenRows{k} = summarizePSSScreenLocal( ...
        raster(k,:),pss,debug,searchOpts);
end

screenTable = vertcat(screenRows{:});
screenTable = sortrows(screenTable, ...
    ["HasPeriodicPSS","PeriodicPairScore","MaxPSSMetric","MaxPeakToMedian"], ...
    ["descend","descend","descend","descend"]);
selected = selectFullAnalysisCentersLocal(screenTable,searchOpts);
screenTable.SelectedForFullAnalysis = ...
    ismember(screenTable.GSCN,selected.GSCN);

analysisOpts.TrimStartSeconds = searchOpts.FullAnalysisStartSeconds;
analysisOpts.AnalysisDurationSeconds = searchOpts.FullAnalysisDurationSeconds;
fullResults = repmat(struct(),height(selected),1);
fullRows = cell(height(selected),1);
pssParts = cell(height(selected),1);
sssParts = cell(height(selected),1);
cellTimingParts = cell(height(selected),1);
detectionTimingParts = cell(height(selected),1);
powerParts = cell(height(selected),1);

for k = 1:height(selected)
    if searchOpts.Verbose
        fprintf("Full SSB analysis %d/%d: GSCN %d, %.3f MHz\n", ...
            k,height(selected),selected.GSCN(k),selected.SSBCenterHz(k)/1e6);
    end
    shifted = digitalShiftLocal( ...
        iq,meta.SampleRateHz,selected.DigitalShiftHz(k),0);
    candidate = analyzeCapture(shifted,meta,analysisOpts);
    fullResults(k).GSCN = selected.GSCN(k);
    fullResults(k).SSBCenterHz = selected.SSBCenterHz(k);
    fullResults(k).DigitalShiftHz = selected.DigitalShiftHz(k);
    fullResults(k).Result = candidate;
    fullRows{k} = summarizeFullResultLocal(candidate,selected(k,:));

    pssParts{k} = addFrequencyIdentityLocal( ...
        candidate.PSSCandidates,selected(k,:),false);
    sssParts{k} = addFrequencyIdentityLocal( ...
        candidate.SSSDetections,selected(k,:),true);
    cellTimingParts{k} = addFrequencyIdentityLocal( ...
        candidate.CellTiming,selected(k,:),true);
    detectionTimingParts{k} = addFrequencyIdentityLocal( ...
        candidate.DetectionTiming,selected(k,:),true);
    if isfield(candidate,"PowerMeasurements")
        powerParts{k} = addFrequencyIdentityLocal( ...
            candidate.PowerMeasurements,selected(k,:),true);
    end
end

fullTable = vertcat(fullRows{:});
fullTable = sortrows(fullTable, ...
    ["IsConfirmed","NumCellTimingRows","MaxDetectionsPerCell", ...
    "BestSSSMetric","MaxPSSMetric"], ...
    ["descend","descend","descend","descend","descend"]);

allPSS = vertcatNonemptyLocal(pssParts);
allSSS = vertcatNonemptyLocal(sssParts);
allCellTiming = vertcatNonemptyLocal(cellTimingParts);
allDetectionTiming = vertcatNonemptyLocal(detectionTimingParts);
allPower = vertcatNonemptyLocal(powerParts);

confirmedGSCN = fullTable.GSCN(fullTable.IsConfirmed);
cellTiming = filterConfirmedCentersLocal(allCellTiming,confirmedGSCN);
detectionTiming = filterConfirmedCentersLocal(allDetectionTiming,confirmedGSCN);
sssDetections = filterConfirmedCentersLocal(allSSS,confirmedGSCN);
powerMeasurements = filterConfirmedCentersLocal(allPower,confirmedGSCN);
powerMeasurements = filterConfirmedCellKeysLocal(powerMeasurements,cellTiming);
detectedCells = detectedCellSummaryLocal(cellTiming);
pciReuse = pciReuseSummaryLocal(cellTiming);

result = struct();
result.Mode = "single_capture_multi_gscn_analysis";
result.Statement = "One wideband IQ capture was digitally searched at every complete legal n78 SSB raster centre. Confirmed cells share the original capture sample-time axis.";
result.Metadata = meta;
result.MetadataReport = metaReport;
result.PSSCandidates = allPSS;
result.SSSDetections = sssDetections;
result.CellTiming = cellTiming;
result.DetectionTiming = detectionTiming;
result.DetectedCells = detectedCells;
result.PowerMeasurements = powerMeasurements;
result.PerFrequencyResults = fullResults;
result.SSBRasterSearch = struct();
result.SSBRasterSearch.Enable = true;
result.SSBRasterSearch.AllRasterCenters = raster;
result.SSBRasterSearch.ScreenTable = screenTable;
result.SSBRasterSearch.FullAnalysisTable = fullTable;
result.SSBRasterSearch.ConfirmedCenters = ...
    fullTable(fullTable.IsConfirmed,:);
result.SSBRasterSearch.NumRasterCenters = height(raster);
result.SSBRasterSearch.NumFullAnalysisCenters = height(selected);
result.SSBRasterSearch.NumConfirmedCenters = nnz(fullTable.IsConfirmed);
result.SSBRasterSearch.SingleCaptureCommonTimeBase = true;
result.SSBRasterSearch.PCIReuseAcrossCenters = pciReuse;
result.SSBRasterSearch.HasAmbiguousPCIReuse = ~isempty(pciReuse) && ...
    height(pciReuse) > 0;
result.SSBRasterSearch.Statement = ...
    "Frequency translation changes carrier phase but not sample arrival time. Cell identity is retained as GSCN plus PCI.";
end

function [window,originSample0] = screenWindowLocal(iq,meta,opts)
originSample0 = max(0,round(opts.ScreenStartSeconds*meta.SampleRateHz));
if originSample0 >= numel(iq)
    error("analyzeCaptureAcrossSSBRaster:ScreenStartOutsideCapture", ...
        "ScreenStartSeconds is outside the capture.");
end
count = min(numel(iq)-originSample0, ...
    max(1,round(opts.ScreenDurationSeconds*meta.SampleRateHz)));
window = iq(originSample0+(1:count));
end

function shifted = digitalShiftLocal(iq,fs,offsetHz,originSample0)
n = originSample0 + (0:numel(iq)-1).';
rotation = exp(-1i*2*pi*offsetHz/fs.*n);
if isa(iq,"single")
    rotation = single(rotation);
end
shifted = iq(:).*rotation;
end

function row = summarizePSSScreenLocal(rasterRow,pss,debug,opts)
maxMetric = max(debug.MetricMax,[],"omitnan");
maxPeakToMedian = NaN;
if ~isempty(pss) && height(pss)>0
    maxPeakToMedian = max(pss.PeakToMedian,[],"omitnan");
end

bestFound = false;
bestScore = 0;
bestNID2 = NaN;
bestPeriodMs = NaN;
bestErrorUs = NaN;
for nid2 = 0:2
    rows = pss(pss.NID2==nid2,:);
    for a = 1:height(rows)
        for b = a+1:height(rows)
            deltaUs = abs(rows.StartTimeUs(b)-rows.StartTimeUs(a));
            [errorUs,periodIndex] = min(abs( ...
                deltaUs-opts.ValidPSSPeriodsMs*1e3));
            if errorUs <= opts.PeriodToleranceUs
                pairStrength = sqrt(rows.PeakMetric(a)*rows.PeakMetric(b));
                score = pairStrength/(1+errorUs/opts.PeriodToleranceUs);
                if score > bestScore
                    bestFound = true;
                    bestScore = score;
                    bestNID2 = nid2;
                    bestPeriodMs = opts.ValidPSSPeriodsMs(periodIndex);
                    bestErrorUs = errorUs;
                end
            end
        end
    end
end

row = table( ...
    rasterRow.GSCN,rasterRow.SSBCenterHz,rasterRow.DigitalShiftHz, ...
    height(pss),maxMetric,maxPeakToMedian,bestFound,bestNID2, ...
    bestPeriodMs,bestErrorUs,bestScore, ...
    'VariableNames',["GSCN","SSBCenterHz","DigitalShiftHz", ...
    "NumPSSCandidates","MaxPSSMetric","MaxPeakToMedian", ...
    "HasPeriodicPSS","PeriodicNID2","PeriodicMs","PeriodErrorUs", ...
    "PeriodicPairScore"]);
end

function selected = selectFullAnalysisCentersLocal(screenTable,opts)
eligible = screenTable(screenTable.HasPeriodicPSS,:);
if height(eligible)>opts.MaxFullAnalysisCenters
    eligible = eligible(1:opts.MaxFullAnalysisCenters,:);
end
if height(eligible)<opts.MinFullAnalysisCenters
    needed = opts.MinFullAnalysisCenters-height(eligible);
    fallback = screenTable(~ismember(screenTable.GSCN,eligible.GSCN),:);
    fallback = fallback(1:min(needed,height(fallback)),:);
    eligible = [eligible;fallback];
end
selected = eligible;
end

function row = summarizeFullResultLocal(candidate,selectedRow)
numTiming = height(candidate.CellTiming);
usable = candidate.SSSDetections(candidate.SSSDetections.IsUsable,:);
maxPSS = NaN;
if ~isempty(candidate.PSSCandidates)
    maxPSS = max(candidate.PSSCandidates.PeakMetric,[],"omitnan");
end
bestSSS = NaN;
maxDetections = 0;
cellList = "";
if ~isempty(usable)
    bestSSS = max(usable.SSSMetric,[],"omitnan");
end
if numTiming>0
    maxDetections = max(candidate.CellTiming.NumDetections);
    cellList = strjoin( ...
        "GSCN"+string(selectedRow.GSCN)+"_PCI"+ ...
        string(candidate.CellTiming.PCI),"|");
end
isConfirmed = numTiming>0 && maxDetections>=2;
if isConfirmed
    reason = "CONFIRMED_REPEATED_PCI_TIMING";
elseif isempty(usable)
    reason = "NO_USABLE_PCI";
else
    reason = "NO_REPEATED_CELL_TIMING";
end
row = table( ...
    selectedRow.GSCN,selectedRow.SSBCenterHz,selectedRow.DigitalShiftHz, ...
    height(candidate.PSSCandidates),height(usable),numTiming,maxDetections, ...
    maxPSS,bestSSS,cellList,isConfirmed,reason, ...
    'VariableNames',["GSCN","SSBCenterHz","DigitalShiftHz", ...
    "NumPSSCandidates","NumUsableDetections","NumCellTimingRows", ...
    "MaxDetectionsPerCell","MaxPSSMetric","BestSSSMetric", ...
    "DetectedCellKeys","IsConfirmed","Reason"]);
end

function tbl = addFrequencyIdentityLocal(tbl,rasterRow,withCellKey)
if isempty(tbl) || height(tbl)==0
    tbl = table();
    return;
end
tbl.GSCN = repmat(rasterRow.GSCN,height(tbl),1);
tbl.SSBCenterHz = repmat(rasterRow.SSBCenterHz,height(tbl),1);
tbl.DigitalShiftHz = repmat(rasterRow.DigitalShiftHz,height(tbl),1);
if withCellKey && ismember("PCI",string(tbl.Properties.VariableNames))
    tbl.CellKey = "GSCN"+string(tbl.GSCN)+"_PCI"+string(tbl.PCI);
end
front = ["GSCN","SSBCenterHz","DigitalShiftHz"];
if ismember("CellKey",string(tbl.Properties.VariableNames))
    front(end+1) = "CellKey";
end
tbl = movevars(tbl,front,"Before",1);
end

function out = vertcatNonemptyLocal(parts)
keep = ~cellfun(@isempty,parts);
if any(keep)
    out = vertcat(parts{keep});
else
    out = table();
end
end

function tbl = filterConfirmedCentersLocal(tbl,confirmedGSCN)
if isempty(tbl) || height(tbl)==0
    return;
end
tbl = tbl(ismember(tbl.GSCN,confirmedGSCN),:);
end

function tbl = filterConfirmedCellKeysLocal(tbl,cellTiming)
if isempty(tbl) || height(tbl)==0
    return;
end
if isempty(cellTiming) || height(cellTiming)==0 || ...
        ~ismember("CellKey",string(tbl.Properties.VariableNames))
    tbl = tbl([],:);
    return;
end
tbl = tbl(ismember(tbl.CellKey,cellTiming.CellKey),:);
end

function summary = detectedCellSummaryLocal(cellTiming)
if isempty(cellTiming) || height(cellTiming)==0
    summary = table();
    return;
end
summary = cellTiming(:,["CellKey","GSCN","SSBCenterHz","PCI", ...
    "NumDetections","FramePhaseNs","TimingStdNs","MeanCFOHz", ...
    "MeanPSSMetric","MeanSSSMetric","AnchorMode"]);
summary = sortrows(summary, ...
    ["GSCN","MeanSSSMetric"],["ascend","descend"]);
end

function reuse = pciReuseSummaryLocal(cellTiming)
if isempty(cellTiming) || height(cellTiming)==0
    reuse = table();
    return;
end
pcis = unique(cellTiming.PCI);
parts = cell(numel(pcis),1);
count = 0;
for k = 1:numel(pcis)
    rows = cellTiming(cellTiming.PCI==pcis(k),:);
    centers = unique(rows.GSCN);
    if numel(centers)>1
        count = count+1;
        parts{count} = table(pcis(k),numel(centers), ...
            strjoin(string(centers.'),"|"), ...
            'VariableNames',["PCI","NumGSCNCenters","GSCNList"]);
    end
end
if count==0
    reuse = table();
else
    reuse = vertcat(parts{1:count});
end
end

function opts = defaultSearchOptionsLocal(meta)
opts = struct();
opts.Enable = true;
opts.SCSkHz = 30;
if isfield(meta,"BandwidthHz") && isfinite(meta.BandwidthHz) && meta.BandwidthHz>0
    opts.CaptureBandwidthHz = meta.BandwidthHz;
else
    opts.CaptureBandwidthHz = meta.SampleRateHz;
end
opts.EdgeGuardHz = 0;
opts.ScreenStartSeconds = 0.003;
opts.ScreenDurationSeconds = 0.06;
opts.FullAnalysisStartSeconds = 0.003;
opts.FullAnalysisDurationSeconds = 0.12;
opts.ValidPSSPeriodsMs = [5 10 20 40];
opts.PeriodToleranceUs = 3;
opts.MinFullAnalysisCenters = 1;
opts.MaxFullAnalysisCenters = Inf;
opts.Verbose = false;
opts.ScreenPSS = struct( ...
    "MaxCandidatesPerNID2",16, ...
    "MaxTotalCandidates",48, ...
    "MinPeakMetric",0.003, ...
    "MinPeakToMedianRatio",12);
end

function opts = removeSearchOptionsLocal(opts)
if isfield(opts,"SSBRasterSearch")
    opts = rmfield(opts,"SSBRasterSearch");
end
if isfield(opts,"SSBOffsetSearch")
    opts = rmfield(opts,"SSBOffsetSearch");
end
end

function out = mergeStructsLocal(base,override)
out = base;
if isempty(override) || isempty(fieldnames(override))
    return;
end
names = fieldnames(override);
for k = 1:numel(names)
    name = names{k};
    if isstruct(override.(name)) && isfield(base,name) && ...
            isstruct(base.(name))
        out.(name) = mergeStructsLocal(base.(name),override.(name));
    else
        out.(name) = override.(name);
    end
end
end
