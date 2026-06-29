function raster = enumerateN78SSBRaster(meta, opts)
%enumerateN78SSBRaster Legal n78 SSB centres fully contained in one capture.
%
% The n78/FR1 synchronization raster above 3 GHz is spaced by 1.44 MHz:
%   F_SSB_MHz = 3000 + 1.44 * (GSCN - 7499).
% Only centres whose complete 240-subcarrier SSB fits inside the usable
% complex-IQ bandwidth are returned.

arguments
    meta struct
    opts struct = struct()
end

[meta, ~] = validateMetadata(meta);
opts = defaultsLocal(opts, meta);

usableBandwidthHz = min(meta.SampleRateHz, opts.CaptureBandwidthHz);
ssbBandwidthHz = 240 * opts.SCSkHz * 1e3;
halfUsableHz = usableBandwidthHz / 2;
halfSSBHz = ssbBandwidthHz / 2;

lowestCenterHz = meta.CenterFrequencyHz - halfUsableHz + ...
    halfSSBHz + opts.EdgeGuardHz;
highestCenterHz = meta.CenterFrequencyHz + halfUsableHz - ...
    halfSSBHz - opts.EdgeGuardHz;

if lowestCenterHz > highestCenterHz
    raster = emptyRasterLocal();
    return;
end

firstGSCN = ceil(7499 + (lowestCenterHz/1e6 - 3000)/1.44);
lastGSCN = floor(7499 + (highestCenterHz/1e6 - 3000)/1.44);
firstGSCN = max(firstGSCN, opts.MinGSCN);
lastGSCN = min(lastGSCN, opts.MaxGSCN);
if firstGSCN > lastGSCN
    raster = emptyRasterLocal();
    return;
end

gscn = (firstGSCN:lastGSCN).';
ssbCenterHz = (3000 + 1.44*(double(gscn)-7499))*1e6;
digitalShiftHz = ssbCenterHz - meta.CenterFrequencyHz;
lowerEdgeMarginHz = ssbCenterHz - halfSSBHz - ...
    (meta.CenterFrequencyHz-halfUsableHz);
upperEdgeMarginHz = (meta.CenterFrequencyHz+halfUsableHz) - ...
    (ssbCenterHz+halfSSBHz);

raster = table( ...
    gscn, ssbCenterHz, digitalShiftHz, ...
    repmat(ssbBandwidthHz,numel(gscn),1), ...
    lowerEdgeMarginHz, upperEdgeMarginHz, ...
    'VariableNames', ["GSCN","SSBCenterHz","DigitalShiftHz", ...
    "SSBOccupiedBandwidthHz","LowerEdgeMarginHz","UpperEdgeMarginHz"]);
end

function opts = defaultsLocal(opts, meta)
defaults = struct();
defaults.SCSkHz = 30;
if isfield(meta,"BandwidthHz") && isfinite(meta.BandwidthHz) && meta.BandwidthHz > 0
    defaults.CaptureBandwidthHz = meta.BandwidthHz;
else
    defaults.CaptureBandwidthHz = meta.SampleRateHz;
end
defaults.EdgeGuardHz = 0;
defaults.MinGSCN = 7499;
defaults.MaxGSCN = floor(7499 + (24250-3000)/1.44);

names = fieldnames(defaults);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(opts,name)
        opts.(name) = defaults.(name);
    end
end
end

function tbl = emptyRasterLocal()
tbl = table( ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames', ["GSCN","SSBCenterHz","DigitalShiftHz", ...
    "SSBOccupiedBandwidthHz","LowerEdgeMarginHz","UpperEdgeMarginHz"]);
end
