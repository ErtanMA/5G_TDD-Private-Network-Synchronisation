% test_phase2_loaders
% Lightweight Phase 2 smoke tests for IQ and metadata loading.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot, "src")));

tmpRoot = tempname;
mkdir(tmpRoot);
cleanup = onCleanup(@() rmdir(tmpRoot, "s"));

metaPath = fullfile(tmpRoot, "meta.json");
fid = fopen(metaPath, "w");
fprintf(fid, ['{"sample_rate_hz":30720000,' ...
    '"center_frequency_hz":3770000000,' ...
    '"bandwidth_hz":15000000,' ...
    '"clock_source":"internal"}']);
fclose(fid);

meta = readCaptureMetadata(metaPath);
[meta, report] = validateMetadata(meta);
assert(report.IsValid);
assert(meta.SampleRateHz == 30720000);

csvPath = fullfile(tmpRoot, "iq.csv");
expected = complex((1:8).', -(1:8).');
writematrix([real(expected), imag(expected)], csvPath);
[iqCsv, infoCsv] = loadIQCapture(csvPath);
assert(isequal(iqCsv, expected));
assert(infoCsv.NumSamples == numel(expected));

binPath = fullfile(tmpRoot, "iq.bin");
fid = fopen(binPath, "wb");
interleaved = zeros(2*numel(expected),1);
interleaved(1:2:end) = real(expected);
interleaved(2:2:end) = imag(expected);
fwrite(fid, single(interleaved), "float32");
fclose(fid);
[iqBin, infoBin] = loadIQCapture(binPath);
assert(max(abs(iqBin - expected)) == 0);
assert(infoBin.NumSamples == numel(expected));

matPath = fullfile(tmpRoot, "iq.mat");
iq = expected; %#ok<NASGU>
save(matPath, "iq");
[iqMat, infoMat] = loadIQCapture(matPath);
assert(isequal(iqMat, expected));
assert(infoMat.NumSamples == numel(expected));

disp("Phase 2 loader tests passed.");

