function [iq, info] = loadIQCapture(capturePath, opts)
%loadIQCapture Load complex IQ samples from CSV, MAT, or binary files.
%
%   [iq, info] = loadIQCapture(capturePath) returns a complex column vector.
%
%   Supported formats:
%     .csv/.txt/.tsv : numeric I and Q columns
%     .mat           : variable named iq, waveform, rxWaveform, samples, or rx
%     .bin/.iq/.dat  : interleaved I,Q samples using opts.BinaryPrecision

arguments
    capturePath (1,1) string
    opts.BinaryPrecision (1,1) string = "float32=>double"
    opts.MachineFormat (1,1) string = "ieee-le"
    opts.Scale (1,1) double = 1
    opts.AllowRealOnly (1,1) logical = false
end

if ~isfile(capturePath)
    error("loadIQCapture:FileNotFound", ...
        "Capture file not found: %s", capturePath);
end

[~,~,ext] = fileparts(capturePath);
ext = lower(string(ext));

switch ext
    case {".csv",".txt",".tsv"}
        [iq, sourceInfo] = loadTextIQ(capturePath, opts.AllowRealOnly);
    case ".mat"
        [iq, sourceInfo] = loadMatIQ(capturePath, opts.AllowRealOnly);
    case {".bin",".iq",".dat"}
        [iq, sourceInfo] = loadBinaryIQ(capturePath, opts.BinaryPrecision, opts.MachineFormat);
    otherwise
        error("loadIQCapture:UnsupportedExtension", ...
            "Unsupported capture extension '%s'.", ext);
end

iq = iq(:) .* opts.Scale;

info = struct();
info.CapturePath = char(capturePath);
info.Extension = char(ext);
info.NumSamples = numel(iq);
info.IsComplex = ~isreal(iq);
info.Source = sourceInfo;

if isempty(iq)
    error("loadIQCapture:NoSamples", "Capture contained no IQ samples.");
end

if any(~isfinite(real(iq))) || any(~isfinite(imag(iq)))
    error("loadIQCapture:NonFiniteSamples", ...
        "Capture contains NaN or Inf samples.");
end

end

function [iq, sourceInfo] = loadTextIQ(path, allowRealOnly)
data = readmatrix(path);
data = data(:, all(~isnan(data), 1));

if isempty(data)
    error("loadIQCapture:EmptyTextFile", ...
        "Text capture has no numeric columns.");
end

if size(data,2) >= 2
    iq = complex(data(:,1), data(:,2));
    layout = "I,Q numeric columns";
elseif allowRealOnly
    iq = complex(data(:,1), zeros(size(data,1),1));
    layout = "real-only numeric column";
else
    error("loadIQCapture:TextNeedsIQColumns", ...
        "Text captures must contain at least two numeric columns: I and Q.");
end

sourceInfo = struct("Format","text","Layout",layout);
end

function [iq, sourceInfo] = loadMatIQ(path, allowRealOnly)
s = load(path);
candidateNames = ["iq","waveform","rxWaveform","samples","rx"];

found = "";
raw = [];
for k = 1:numel(candidateNames)
    name = char(candidateNames(k));
    if isfield(s, name)
        raw = s.(name);
        found = candidateNames(k);
        break;
    end
end

if strlength(found) == 0
    names = string(fieldnames(s));
    for k = 1:numel(names)
        value = s.(char(names(k)));
        if isnumeric(value) && ~isempty(value)
            raw = value;
            found = names(k);
            break;
        end
    end
end

if strlength(found) == 0
    error("loadIQCapture:NoMatIQVariable", ...
        "MAT file contains no numeric IQ-like variable.");
end

iq = numericToIQ(raw, allowRealOnly);
sourceInfo = struct("Format","mat","Variable",found);
end

function [iq, sourceInfo] = loadBinaryIQ(path, precision, machineFormat)
fid = fopen(path, "rb", char(machineFormat));
if fid < 0
    error("loadIQCapture:OpenFailed", ...
        "Could not open binary capture: %s", path);
end
cleanup = onCleanup(@() fclose(fid));

data = fread(fid, Inf, char(precision));
if mod(numel(data), 2) ~= 0
    error("loadIQCapture:OddBinarySampleCount", ...
        "Binary interleaved IQ file contains an odd number of scalar samples.");
end

iq = complex(data(1:2:end), data(2:2:end));
sourceInfo = struct("Format","binary","Layout","interleaved I,Q", ...
    "Precision",precision,"MachineFormat",machineFormat);
end

function iq = numericToIQ(raw, allowRealOnly)
if ~isnumeric(raw)
    error("loadIQCapture:NonNumericData", "IQ data must be numeric.");
end

if ~isreal(raw)
    iq = raw(:);
    return;
end

if ismatrix(raw) && size(raw,2) >= 2
    iq = complex(raw(:,1), raw(:,2));
elseif allowRealOnly
    iq = complex(raw(:), zeros(numel(raw),1));
else
    error("loadIQCapture:RealDataNeedsIQColumns", ...
        "Real-valued IQ data must have I and Q columns.");
end
end

