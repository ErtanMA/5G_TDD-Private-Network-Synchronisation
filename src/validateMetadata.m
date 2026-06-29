function [meta, report] = validateMetadata(meta)
%validateMetadata Validate and canonicalize capture metadata.
%
%   [meta, report] = validateMetadata(meta) checks the fields required for
%   Phase 2 data handling. It intentionally does not require GNSS timing,
%   because the current project mode is internal-clock relative timing.

meta = canonicalizeMetadata(meta);

errors = strings(0,1);
warnings = strings(0,1);

[errors, meta.SampleRateHz] = requirePositiveScalar(errors, meta.SampleRateHz, "SampleRateHz");
[errors, meta.CenterFrequencyHz] = requirePositiveScalar(errors, meta.CenterFrequencyHz, "CenterFrequencyHz");
[errors, meta.BandwidthHz] = requirePositiveScalar(errors, meta.BandwidthHz, "BandwidthHz");

if isempty(errors)
    [errors, warnings, meta.SampleRateHz] = requireIntegerHzLocal( ...
        errors, warnings, meta.SampleRateHz, "SampleRateHz");
end

if ~isnan(meta.RxLat) && (meta.RxLat < -90 || meta.RxLat > 90)
    errors(end+1,1) = "RxLat must be between -90 and 90 degrees.";
end

if ~isnan(meta.RxLon) && (meta.RxLon < -180 || meta.RxLon > 180)
    errors(end+1,1) = "RxLon must be between -180 and 180 degrees.";
end

if isempty(meta.GnssLocked)
    warnings(end+1,1) = "GNSS lock status is absent. This is acceptable for relative mode.";
elseif islogical(meta.GnssLocked) && meta.GnssLocked
    warnings(end+1,1) = "GNSS lock was reported, but absolute mode is not implemented in Phase 2.";
end

if strlength(string(meta.CaptureStartTime)) == 0
    warnings(end+1,1) = "Capture start time is absent. Relative single-capture mode can continue.";
end

report = struct();
report.IsValid = isempty(errors);
report.Errors = errors;
report.Warnings = warnings;
report.Mode = "relative_internal_clock";

if ~report.IsValid
    error("validateMetadata:InvalidMetadata", ...
        "Invalid metadata:\n%s", strjoin(errors, newline));
end

end

function [errors, warnings, value] = requireIntegerHzLocal(errors, warnings, value, fieldName)
roundedValue = round(value);
if abs(value - roundedValue) <= 1
    if value ~= roundedValue
        warnings(end+1,1) = fieldName + " was rounded from " + ...
            string(sprintf("%.12g", value)) + " Hz to " + ...
            string(sprintf("%.0f", roundedValue)) + ...
            " Hz to satisfy integer-valued waveform generation.";
    end
    value = roundedValue;
else
    errors(end+1,1) = fieldName + ...
        " must be integer-valued in Hz for 5G Toolbox waveform generation.";
end
end

function [errors, value] = requirePositiveScalar(errors, value, fieldName)
if isempty(value)
    errors(end+1,1) = fieldName + " is required.";
    value = NaN;
    return;
end

if ischar(value) || isstring(value)
    value = str2double(value);
end

if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value <= 0
    errors(end+1,1) = fieldName + " must be a positive finite scalar.";
end
end

function meta = canonicalizeMetadata(rawMeta)
% Normalize common metadata field spellings from JSON sidecars or MATLAB
% structs into the field names used by the receiver.
if ~isstruct(rawMeta)
    error("validateMetadata:InvalidInput", ...
        "Metadata must be a scalar struct decoded from JSON or built in MATLAB.");
end

if numel(rawMeta) ~= 1
    error("validateMetadata:InvalidInput", ...
        "Metadata must be a scalar struct, not an array.");
end

meta = struct();
meta.SampleRateHz = getFirstPresent(rawMeta, ...
    ["SampleRateHz","sample_rate_hz","sampleRateHz","sampleRate","fs","Fs"]);
meta.CenterFrequencyHz = getFirstPresent(rawMeta, ...
    ["CenterFrequencyHz","center_frequency_hz","centerFrequencyHz","centerFrequency","fc","Fc"]);
meta.BandwidthHz = getFirstPresent(rawMeta, ...
    ["BandwidthHz","bandwidth_hz","bandwidth","capture_bandwidth_hz"]);
meta.RxLat = getFirstPresent(rawMeta, ...
    ["RxLat","rx_lat","rxLat","receiver_lat","receiverLatitude","latitude"], NaN);
meta.RxLon = getFirstPresent(rawMeta, ...
    ["RxLon","rx_lon","rxLon","receiver_lon","receiverLongitude","longitude"], NaN);
meta.RxAltM = getFirstPresent(rawMeta, ...
    ["RxAltM","rx_alt_m","rxAltM","receiver_alt_m","receiverAltitudeM","altitude_m"], NaN);
meta.CaptureStartTime = getFirstPresent(rawMeta, ...
    ["CaptureStartTime","capture_start_time","first_sample_time","first_sample_utc","timestamp_utc"], "");
meta.ClockSource = getFirstPresent(rawMeta, ...
    ["ClockSource","clock_source","clockSource"], "internal");
meta.GnssLocked = getFirstPresent(rawMeta, ...
    ["GnssLocked","gnss_locked","gps_locked","GpsLocked"], []);
meta.Raw = rawMeta;
end

function value = getFirstPresent(s, names, defaultValue)
if nargin < 3
    defaultValue = [];
end

value = defaultValue;
for k = 1:numel(names)
    name = char(names(k));
    if isfield(s, name)
        value = s.(name);
        return;
    end
end
end
