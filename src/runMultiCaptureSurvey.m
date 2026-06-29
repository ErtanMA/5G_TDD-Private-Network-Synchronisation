function result = runMultiCaptureSurvey(manifest, opts)
%runMultiCaptureSurvey Run the full real-data survey workflow.
%
%   result = runMultiCaptureSurvey(manifest, opts) loads each IQ capture in a
%   manifest, runs analyzeCapture, converts detected timing into survey
%   measurement rows, and runs noLocationSurveyCheck.
%
%   manifest can be a table or a CSV path. Required columns:
%     CaptureID, IQPath, MetadataPath, and receiver position fields
%     (RxX_m/RxY_m or RxLat/RxLon).

arguments
    manifest
    opts struct = struct()
end

[manifestTable, manifestBaseDir] = readManifestLocal(manifest);
manifestTable = canonicalizeManifestLocal(manifestTable, manifestBaseDir);
if height(manifestTable) == 0
    error("runMultiCaptureSurvey:EmptyManifest", ...
        "The survey manifest has no accepted captures.");
end

if ~isfield(opts, "Analyze")
    opts.Analyze = struct();
end
if ~isfield(opts, "Survey")
    opts.Survey = struct();
end
if ~isfield(opts, "Loader")
    opts.Loader = struct();
end
if ~isfield(opts, "RSS")
    opts.RSS = struct();
end
if ~isfield(opts, "EnableRSSPlanning")
    opts.EnableRSSPlanning = false;
end
if ~isfield(opts, "WriteReport")
    opts.WriteReport = false;
end
if ~isfield(opts, "OutputDir")
    opts.OutputDir = "";
end
if ~isfield(opts, "ReportBaseName")
    opts.ReportBaseName = "multi_capture_survey";
end

numCaptures = height(manifestTable);
captureInfo = manifestToCaptureInfoLocal(manifestTable);
captureInfoLocal = receiverXYLocal(captureInfo);
validateDistinctReceiverPositionsLocal(captureInfoLocal);
captureResults = repmat(struct(), numCaptures, 1);
loadRows = cell(numCaptures, 1);
measurementRows = cell(numCaptures, 1);
powerRows = cell(numCaptures, 1);
collectPower = shouldCollectPowerLocal(opts);

for k = 1:numCaptures
    capID = string(manifestTable.CaptureID(k));
    meta = readCaptureMetadata(string(manifestTable.MetadataPath(k)));
    [iq, loadInfo] = loadIQWithOptionsLocal(string(manifestTable.IQPath(k)), opts.Loader);
    validateSavedCaptureLocal(meta, loadInfo, capID);
    localPosition = captureInfoLocal(k,:);
    analyzeOpts = analysisOptionsForCaptureLocal( ...
        opts, localPosition, capID, collectPower);
    captureResult = analyzeCaptureWithOptionalOffsetSearchLocal(iq, meta, analyzeOpts);
    rows = timingMeasurementsFromCaptureResult(captureResult, capID);
    powerRows{k} = powerRowsFromCaptureLocal( ...
        captureResult, localPosition, capID, collectPower);

    captureResults(k).CaptureID = capID;
    captureResults(k).Metadata = meta;
    captureResults(k).LoadInfo = loadInfo;
    captureResults(k).Result = captureResult;

    loadRows{k} = table(capID, string(loadInfo.CapturePath), loadInfo.NumSamples, ...
        height(rows), 'VariableNames', ["CaptureID","IQPath","NumSamples", ...
        "NumMeasurementRows"]);
    measurementRows{k} = rows;
end

if isempty(measurementRows)
    measurementTable = table();
else
    measurementTable = vertcat(measurementRows{:});
end

if collectPower
    powerTable = vertcat(powerRows{:});
else
    powerTable = table();
end

[measurementTable, powerTable, ambiguousRasterPCIs] = ...
    excludeAmbiguousRasterPCIReuseLocal(measurementTable, powerTable);

surveyResult = noLocationSurveyCheck(measurementTable, captureInfo, opts.Survey);

if collectPower && ~isempty(powerTable) && height(powerTable) > 0
    candidateStops = getCandidateStopsLocal(opts.RSS);
    rssPlan = estimateRssPriorsAndPlanStops("all", powerTable, captureInfo, candidateStops, opts.RSS);
    surveyResult.PowerTable = powerTable;
    surveyResult.RssPriors = rssPlan.RssPriors;
    surveyResult.SuggestedReceiverStops = rssPlan.SuggestedReceiverStops;
    surveyResult.RssPlanningStatement = rssPlan.Statement;
else
    rssPlan = struct();
end

result = struct();
result.Mode = "multi_capture_survey_from_iq";
result.Statement = "Loaded IQ captures, analyzed detected PCIs, and ran no-gNB-location survey fitting.";
result.Manifest = manifestTable;
result.CaptureInfo = captureInfo;
result.CaptureInfoLocal = captureInfoLocal;
result.LoadSummary = vertcat(loadRows{:});
result.MeasurementTable = measurementTable;
result.PowerTable = powerTable;
result.ExcludedAmbiguousRasterPCIs = ambiguousRasterPCIs;
result.CaptureResults = captureResults;
result.SurveyResult = surveyResult;
if isfield(rssPlan, "RssPriors")
    result.RssPriors = rssPlan.RssPriors;
    result.SuggestedReceiverStops = rssPlan.SuggestedReceiverStops;
end

if opts.WriteReport
    if strlength(string(opts.OutputDir)) == 0
        error("runMultiCaptureSurvey:MissingOutputDir", ...
            "OutputDir must be set when WriteReport is true.");
    end
    result.ReportPaths = generateSurveyReport(surveyResult, string(opts.OutputDir), string(opts.ReportBaseName));
end

end

function captureResult = analyzeCaptureWithOptionalOffsetSearchLocal(iq, meta, analyzeOpts)
if searchEnabledLocal(analyzeOpts,"SSBRasterSearch")
    captureResult = analyzeCaptureAcrossSSBRaster(iq, meta, analyzeOpts);
elseif searchEnabledLocal(analyzeOpts,"SSBOffsetSearch")
    captureResult = analyzeCaptureWithSSBOffsetSearch(iq, meta, analyzeOpts);
else
    captureResult = analyzeCapture(iq, meta, analyzeOpts);
end
end

function tf = searchEnabledLocal(opts,name)
tf = isfield(opts,name);
if ~tf
    return;
end
searchOpts = opts.(name);
if isfield(searchOpts,"Enable")
    tf = logical(searchOpts.Enable);
end
end

function tf = shouldCollectPowerLocal(opts)
tf = false;
if isfield(opts, "EnableRSSPlanning") && opts.EnableRSSPlanning
    tf = true;
end
if isfield(opts, "Analyze") && isfield(opts.Analyze, "EnablePowerMeasurements") && ...
        opts.Analyze.EnablePowerMeasurements
    tf = true;
end
if isfield(opts, "RSS") && ~isempty(fieldnames(opts.RSS))
    tf = true;
end
end

function analyzeOpts = analysisOptionsForCaptureLocal(opts, localPosition, capID, collectPower)
analyzeOpts = opts.Analyze;
if ~collectPower
    return;
end

analyzeOpts.EnablePowerMeasurements = true;
analyzeOpts.CaptureID = capID;
if isfield(opts, "RSS")
    if isfield(analyzeOpts, "Power")
        analyzeOpts.Power = mergeStructsLocal(opts.RSS, analyzeOpts.Power);
    else
        analyzeOpts.Power = opts.RSS;
    end
end

for name = ["RxX_m","RxY_m","RxZ_m"]
    analyzeOpts.(name) = localPosition.(name);
end
end

function rows = powerRowsFromCaptureLocal(captureResult, localPosition, capID, collectPower)
if ~collectPower || ~isfield(captureResult, "PowerMeasurements") || ...
        isempty(captureResult.PowerMeasurements)
    rows = emptyPowerRowsLocal();
    return;
end

rows = captureResult.PowerMeasurements;
if height(rows) == 0
    rows = emptyPowerRowsLocal();
    return;
end
rows.CaptureID = repmat(capID, height(rows), 1);
names = string(rows.Properties.VariableNames);
if ~ismember("CellKey",names)
    rows.CellKey = "PCI"+string(rows.PCI);
end
if ~ismember("GSCN",names)
    rows.GSCN = NaN(height(rows),1);
end
if ~ismember("SSBCenterHz",names)
    rows.SSBCenterHz = NaN(height(rows),1);
end
if ~ismember("DigitalShiftHz",names)
    rows.DigitalShiftHz = zeros(height(rows),1);
end

rows.RxX_m(:) = localPosition.RxX_m;
rows.RxY_m(:) = localPosition.RxY_m;
rows.RxZ_m(:) = localPosition.RxZ_m;
rows = rows(:,["CaptureID","CellKey","GSCN","SSBCenterHz", ...
    "DigitalShiftHz","PCI","RxX_m","RxY_m","RxZ_m", ...
    "PowerDbFS","NoiseDbFS","SNRdB","NumDetections","QualityFlag"]);
end

function tbl = emptyPowerRowsLocal()
tbl = table( ...
    strings(0,1), strings(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), strings(0,1), ...
    'VariableNames', ["CaptureID","CellKey","GSCN","SSBCenterHz", ...
    "DigitalShiftHz","PCI","RxX_m","RxY_m","RxZ_m", ...
    "PowerDbFS","NoiseDbFS","SNRdB","NumDetections","QualityFlag"]);
end

function [measurements,power,ambiguousPCIs] = ...
    excludeAmbiguousRasterPCIReuseLocal(measurements,power)
ambiguousPCIs = zeros(0,1);
if isempty(measurements) || height(measurements)==0 || ...
        ~ismember("CellKey",string(measurements.Properties.VariableNames))
    return;
end

pcis = unique(measurements.PCI);
for k = 1:numel(pcis)
    keys = unique(measurements.CellKey(measurements.PCI==pcis(k)));
    keys = keys(strlength(keys)>0);
    if numel(keys)>1
        ambiguousPCIs(end+1,1) = pcis(k); %#ok<AGROW>
    end
end
if isempty(ambiguousPCIs)
    return;
end

warning("runMultiCaptureSurvey:PCIReuseAcrossRasterCenters", ...
    "Excluded PCI(s) %s from survey fitting because each appeared on multiple GSCN centres. The current geometry estimator is PCI-indexed and must not merge those cells.", ...
    strjoin(string(ambiguousPCIs.'),", "));
measurements = measurements(~ismember(measurements.PCI,ambiguousPCIs),:);
if ~isempty(power) && height(power)>0
    power = power(~ismember(power.PCI,ambiguousPCIs),:);
end
end

function candidateStops = getCandidateStopsLocal(rssOpts)
if isfield(rssOpts, "CandidateStops")
    candidateStops = rssOpts.CandidateStops;
else
    candidateStops = table();
end
end

function [tbl, baseDir] = readManifestLocal(manifest)
if istable(manifest)
    tbl = manifest;
    baseDir = string(pwd);
elseif ischar(manifest) || isstring(manifest)
    manifestPath = string(manifest);
    if ~isfile(manifestPath)
        error("runMultiCaptureSurvey:ManifestNotFound", ...
            "Manifest file not found: %s", manifestPath);
    end
    tbl = readtable(manifestPath, "TextType", "string", ...
        "VariableNamingRule", "preserve", ...
        "Delimiter", ",", "ReadVariableNames", true);
    baseDir = string(fileparts(manifestPath));
    if strlength(baseDir) == 0
        baseDir = string(pwd);
    end
else
    error("runMultiCaptureSurvey:InvalidManifest", ...
        "Manifest must be a table or CSV path.");
end
end

function [iq, loadInfo] = loadIQWithOptionsLocal(iqPath, loaderOpts)
if isempty(fieldnames(loaderOpts))
    [iq, loadInfo] = loadIQCapture(iqPath);
    return;
end

names = fieldnames(loaderOpts);
args = cell(1, 2*numel(names));
for k = 1:numel(names)
    args{2*k-1} = names{k};
    args{2*k} = loaderOpts.(names{k});
end
[iq, loadInfo] = loadIQCapture(iqPath, args{:});
end

function tbl = canonicalizeManifestLocal(tbl, baseDir)
names = string(tbl.Properties.VariableNames);
captureID = getVarLocal(tbl, ["CaptureID","capture_id","id"], "cap_" + string((1:height(tbl)).'));
iqPath = string(getVarLocal(tbl, ["IQPath","iq_path","CapturePath","capture_path","iq"], strings(height(tbl),1)));
metaPath = string(getVarLocal(tbl, ["MetadataPath","metadata_path","MetaPath","meta_path","metadata"], strings(height(tbl),1)));

if any(strlength(iqPath) == 0) || any(strlength(metaPath) == 0)
    error("runMultiCaptureSurvey:MissingPaths", ...
        "Manifest must include IQPath and MetadataPath columns. Imported columns: %s", ...
        strjoin(names, ", "));
end

for k = 1:height(tbl)
    iqPath(k) = resolvePathLocal(iqPath(k), baseDir);
    metaPath(k) = resolvePathLocal(metaPath(k), baseDir);
end

out = table(string(captureID(:)), iqPath(:), metaPath(:), ...
    'VariableNames', ["CaptureID","IQPath","MetadataPath"]);

copyIfPresent = ["RxX_m","rx_x_m","RxY_m","rx_y_m","RxZ_m","rx_z_m", ...
    "RxLat","rx_lat","RxLon","rx_lon","RxAltM","rx_alt_m", ...
    "RxPositionUncertaintyM","rx_position_uncertainty_m"];
for k = 1:numel(copyIfPresent)
    name = copyIfPresent(k);
    if ismember(name, names)
        canonical = canonicalManifestNameLocal(name);
        out.(canonical) = tbl.(name);
    end
end

tbl = out;
end

function captureInfo = manifestToCaptureInfoLocal(manifestTable)
captureInfo = table();
captureInfo.CaptureID = string(manifestTable.CaptureID);
names = string(manifestTable.Properties.VariableNames);

if all(ismember(["RxX_m","RxY_m"], names))
    captureInfo.RxX_m = manifestTable.RxX_m;
    captureInfo.RxY_m = manifestTable.RxY_m;
    if ismember("RxZ_m", names)
        captureInfo.RxZ_m = manifestTable.RxZ_m;
    else
        captureInfo.RxZ_m = zeros(height(manifestTable), 1);
    end
elseif all(ismember(["RxLat","RxLon"], names))
    captureInfo.RxLat = manifestTable.RxLat;
    captureInfo.RxLon = manifestTable.RxLon;
    if ismember("RxAltM", names)
        captureInfo.RxAltM = manifestTable.RxAltM;
    else
        captureInfo.RxAltM = zeros(height(manifestTable), 1);
    end
else
    error("runMultiCaptureSurvey:MissingReceiverPosition", ...
        "Manifest must include RxX_m/RxY_m or RxLat/RxLon.");
end

if ismember("RxPositionUncertaintyM", names)
    captureInfo.RxPositionUncertaintyM = manifestTable.RxPositionUncertaintyM;
end
end

function captureXY = receiverXYLocal(captureInfo)
cap = captureInfo;
names = string(cap.Properties.VariableNames);
if all(ismember(["RxX_m","RxY_m"], names))
    x = double(cap.RxX_m);
    y = double(cap.RxY_m);
    if ismember("RxZ_m", names)
        z = double(cap.RxZ_m);
    else
        z = zeros(height(cap),1);
    end
elseif all(ismember(["RxLat","RxLon"], names))
    lat = double(cap.RxLat);
    lon = double(cap.RxLon);
    lat0 = lat(1);
    lon0 = lon(1);
    earthRadiusM = 6371000;
    x = deg2rad(lon-lon0) .* earthRadiusM .* cos(deg2rad(lat0));
    y = deg2rad(lat-lat0) .* earthRadiusM;
    if ismember("RxAltM", names)
        z = double(cap.RxAltM)-double(cap.RxAltM(1));
    else
        z = zeros(height(cap),1);
    end
else
    error("runMultiCaptureSurvey:MissingReceiverPosition", ...
        "Manifest must include RxX_m/RxY_m or RxLat/RxLon.");
end

if ismember("RxPositionUncertaintyM", names)
    uncertainty = double(cap.RxPositionUncertaintyM);
else
    uncertainty = 5*ones(height(cap),1);
end
captureXY = table(string(cap.CaptureID),x(:),y(:),z(:),uncertainty(:), ...
    'VariableNames',["CaptureID","RxX_m","RxY_m","RxZ_m", ...
    "RxPositionUncertaintyM"]);
end

function validateDistinctReceiverPositionsLocal(captureXY)
for a = 1:height(captureXY)
    for b = a+1:height(captureXY)
        separationM = hypot( ...
            captureXY.RxX_m(a)-captureXY.RxX_m(b), ...
            captureXY.RxY_m(a)-captureXY.RxY_m(b));
        indistinguishableM = max( ...
            captureXY.RxPositionUncertaintyM([a b]));
        if separationM <= indistinguishableM
            error("runMultiCaptureSurvey:DuplicateReceiverPosition", ...
                "Captures %s and %s are only %.2f m apart, within the %.2f m recorded position uncertainty. Keep repeatability captures out of the geometry manifest.", ...
                captureXY.CaptureID(a),captureXY.CaptureID(b), ...
                separationM,indistinguishableM);
        end
    end
end
end

function validateSavedCaptureLocal(meta, loadInfo, capID)
if ~isfield(meta,"Raw") || ~isstruct(meta.Raw)
    return;
end
raw = meta.Raw;
if isfield(raw,"capture_validated") && ~logical(raw.capture_validated)
    error("runMultiCaptureSurvey:UnvalidatedCapture", ...
        "Capture %s is explicitly marked invalid.",capID);
end
if isfield(raw,"overrun_count") && double(raw.overrun_count) ~= 0
    error("runMultiCaptureSurvey:CaptureOverrun", ...
        "Capture %s reports %.0f overruns and cannot enter the timing fit.", ...
        capID,double(raw.overrun_count));
end
if isfield(raw,"actual_sample_count") && ...
        double(raw.actual_sample_count) ~= double(loadInfo.NumSamples)
    error("runMultiCaptureSurvey:SampleCountMismatch", ...
        "Capture %s metadata reports %.0f samples, but the IQ loader found %.0f.", ...
        capID,double(raw.actual_sample_count),double(loadInfo.NumSamples));
end
if isfield(raw,"expected_sample_count") && ...
        double(raw.expected_sample_count) ~= double(loadInfo.NumSamples)
    error("runMultiCaptureSurvey:IncompleteCapture", ...
        "Capture %s contains %.0f samples; expected %.0f.", ...
        capID,double(loadInfo.NumSamples),double(raw.expected_sample_count));
end
end

function value = getVarLocal(tbl, names, defaultValue)
value = defaultValue;
tableNames = string(tbl.Properties.VariableNames);
for k = 1:numel(names)
    name = names(k);
    if ismember(name, tableNames)
        value = tbl.(name);
        return;
    end
end
end

function path = resolvePathLocal(path, baseDir)
path = string(path);
if isAbsolutePathLocal(path)
    return;
end
path = fullfile(baseDir, path);
end

function tf = isAbsolutePathLocal(path)
path = char(path);
tf = startsWith(path, filesep) || ...
    (~isempty(regexp(path, "^[A-Za-z]:[\\/]", "once")));
end

function canonical = canonicalManifestNameLocal(name)
switch string(name)
    case {"rx_x_m"}
        canonical = "RxX_m";
    case {"rx_y_m"}
        canonical = "RxY_m";
    case {"rx_z_m"}
        canonical = "RxZ_m";
    case {"rx_lat"}
        canonical = "RxLat";
    case {"rx_lon"}
        canonical = "RxLon";
    case {"rx_alt_m"}
        canonical = "RxAltM";
    case {"rx_position_uncertainty_m"}
        canonical = "RxPositionUncertaintyM";
    otherwise
        canonical = string(name);
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
