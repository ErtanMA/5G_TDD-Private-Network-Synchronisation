function handles = plotCaptureSanity(iq, meta, opts)
%plotCaptureSanity Plot basic time and spectrum checks for a capture.
%
%   This is a Phase 2 diagnostic helper. It does not perform NR detection.

arguments
    iq (:,1) {mustBeNumeric}
    meta struct
    opts.MaxTimeSamples (1,1) double = 20000
    opts.Nfft (1,1) double = 4096
end

[meta, ~] = validateMetadata(meta);
fs = meta.SampleRateHz;

iq = iq(:);
numTime = min(numel(iq), opts.MaxTimeSamples);
t = (0:numTime-1).' / fs * 1e3;

handles = struct();
handles.TimeFigure = figure("Name","Capture Time Sanity");
subplot(2,1,1);
plot(t, real(iq(1:numTime)));
grid on;
xlabel("Time (ms)");
ylabel("I");
title("In-phase samples");

subplot(2,1,2);
plot(t, abs(iq(1:numTime)).^2);
grid on;
xlabel("Time (ms)");
ylabel("|IQ|^2");
title("Power envelope");

handles.SpectrumFigure = figure("Name","Capture Spectrum Sanity");
if exist("pwelch", "file") == 2
    nfft = min(opts.Nfft, 2^floor(log2(numel(iq))));
    nfft = max(nfft, 256);
    win = hann(nfft, "periodic");
    overlap = floor(0.5*nfft);
    [pxx, f] = pwelch(iq, win, overlap, nfft, fs, "centered");
    plot((f + meta.CenterFrequencyHz)/1e6, 10*log10(pxx + eps));
    ylabel("PSD (dB/Hz)");
else
    nfft = min(opts.Nfft, numel(iq));
    x = iq(1:nfft) .* hann(nfft, "periodic");
    spec = fftshift(abs(fft(x, nfft)).^2);
    f = (-nfft/2:nfft/2-1).' / nfft * fs;
    plot((f + meta.CenterFrequencyHz)/1e6, 10*log10(spec + eps));
    ylabel("Power (dB)");
end
grid on;
xlabel("Frequency (MHz)");
title("Capture spectrum sanity check");

end

