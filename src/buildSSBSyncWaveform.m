function [waveform, info] = buildSSBSyncWaveform(pci, meta, opts)
%buildSSBSyncWaveform Build a locally aligned SS/PBCH-block waveform.
%
% The four SSB symbols do not generally begin at slot symbol zero. For the
% n78/30 kHz Case-C path used by this project, candidate SSBs begin at slot
% symbol 2 or 8. Both starts use the same four short-CP symbol lengths. The
% default start symbol is therefore 2.

arguments
    pci (1,1) double {mustBeInteger,mustBeNonnegative}
    meta struct
    opts struct = struct()
end

[meta, ~] = validateMetadata(meta);
opts = defaultsLocal(opts);

carrier = nrCarrierConfig;
carrier.NSizeGrid = opts.NSizeGrid;
carrier.SubcarrierSpacing = opts.SCSkHz;
carrier.NCellID = pci;

ssbGrid = complex(zeros(12*opts.NSizeGrid, 4));
if opts.IncludePSS
    ssbGrid(nrPSSIndices) = nrPSS(mod(pci, 3));
end
if opts.IncludeSSS
    ssbGrid(nrSSSIndices) = nrSSS(pci);
end
if opts.IncludePBCHDMRS
    ssbGrid(nrPBCHDMRSIndices(pci)) = nrPBCHDMRS(pci, opts.SSBIndex);
end

numSymbolsPerSlot = 14;
startSymbol = opts.SSBStartSymbol;
numOutputSymbols = opts.NumOutputSymbols;
if startSymbol < 0 || startSymbol + 4 > numSymbolsPerSlot
    error("buildSSBSyncWaveform:InvalidStartSymbol", ...
        "SSBStartSymbol must allow four SSB symbols inside one slot.");
end
if numOutputSymbols < 1 || numOutputSymbols > 4
    error("buildSSBSyncWaveform:InvalidOutputSymbols", ...
        "NumOutputSymbols must be an integer from 1 to 4.");
end

slotGrid = complex(zeros(12*opts.NSizeGrid, numSymbolsPerSlot));
ssbColumns = startSymbol + (1:4);
slotGrid(:, ssbColumns) = ssbGrid;

slotWaveform = nrOFDMModulate(carrier, slotGrid, ...
    "SampleRate", meta.SampleRateHz);
ofdmInfo = nrOFDMInfo(carrier, "SampleRate", meta.SampleRateHz);
symbolLengths = double(ofdmInfo.SymbolLengths(:));

startSample0 = sum(symbolLengths(1:startSymbol));
referenceLength = sum(symbolLengths(startSymbol + (1:numOutputSymbols)));
waveform = slotWaveform(startSample0 + (1:referenceLength));
waveform = waveform(:);

if opts.RemoveMean
    waveform = waveform - mean(waveform);
end
if opts.NormalizeReferencePower
    waveform = waveform ./ sqrt(mean(abs(waveform).^2) + eps);
end

info = struct();
info.Carrier = carrier;
info.SSBStartSymbol = startSymbol;
info.SSBColumns = ssbColumns;
info.SlotSymbolLengths = symbolLengths;
info.StartSample0 = startSample0;
info.ReferenceLength = referenceLength;
info.NumOutputSymbols = numOutputSymbols;
info.SCSkHz = opts.SCSkHz;
info.NSizeGrid = opts.NSizeGrid;
end

function opts = defaultsLocal(opts)
defaults = struct();
defaults.SCSkHz = 30;
defaults.NSizeGrid = 20;
defaults.SSBStartSymbol = 2;
defaults.NumOutputSymbols = 4;
defaults.IncludePSS = true;
defaults.IncludeSSS = false;
defaults.IncludePBCHDMRS = false;
defaults.SSBIndex = 0;
defaults.RemoveMean = true;
defaults.NormalizeReferencePower = true;

names = fieldnames(defaults);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(opts, name)
        opts.(name) = defaults.(name);
    end
end
end
