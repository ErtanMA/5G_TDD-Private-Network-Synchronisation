function [ssbGrid, info] = demodulateSSBWaveform(segment, meta, opts)
%demodulateSSBWaveform OFDM-demodulate a local four-symbol SSB segment.
%
% The input segment begins at the SSB boundary. It is inserted at the
% configured absolute slot-symbol position before nrOFDMDemodulate so the
% correct cyclic-prefix lengths are used.

arguments
    segment (:,1) {mustBeNumeric}
    meta struct
    opts struct = struct()
end

[meta, ~] = validateMetadata(meta);
opts = defaultsLocal(opts);

carrier = nrCarrierConfig;
carrier.NSizeGrid = opts.NSizeGrid;
carrier.SubcarrierSpacing = opts.SCSkHz;
carrier.NCellID = opts.NCellID;

ofdmInfo = nrOFDMInfo(carrier, "SampleRate", meta.SampleRateHz);
symbolLengths = double(ofdmInfo.SymbolLengths(:));
startSymbol = opts.SSBStartSymbol;
if startSymbol < 0 || startSymbol + 4 > numel(symbolLengths)
    error("demodulateSSBWaveform:InvalidStartSymbol", ...
        "SSBStartSymbol must allow four SSB symbols inside one slot.");
end

startSample0 = sum(symbolLengths(1:startSymbol));
segmentLength = sum(symbolLengths(startSymbol + (1:4)));
if numel(segment) ~= segmentLength
    error("demodulateSSBWaveform:LengthMismatch", ...
        "Expected %d SSB samples for start symbol %d, received %d.", ...
        segmentLength, startSymbol, numel(segment));
end

slotWaveform = complex(zeros(sum(symbolLengths), 1, "like", segment));
slotWaveform(startSample0 + (1:segmentLength)) = segment;
slotGrid = nrOFDMDemodulate(carrier, slotWaveform, ...
    "SampleRate", meta.SampleRateHz);
ssbColumns = startSymbol + (1:4);
ssbGrid = slotGrid(:, ssbColumns);

info = struct();
info.Carrier = carrier;
info.SSBStartSymbol = startSymbol;
info.SSBColumns = ssbColumns;
info.SlotSymbolLengths = symbolLengths;
info.StartSample0 = startSample0;
info.SegmentLength = segmentLength;
end

function opts = defaultsLocal(opts)
defaults = struct();
defaults.SCSkHz = 30;
defaults.NSizeGrid = 20;
defaults.SSBStartSymbol = 2;
defaults.NCellID = 0;

names = fieldnames(defaults);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(opts, name)
        opts.(name) = defaults.(name);
    end
end
end
