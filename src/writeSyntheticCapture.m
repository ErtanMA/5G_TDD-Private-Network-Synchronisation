function paths = writeSyntheticCapture(outputDir, scenarioName, iq, meta, truth)
%writeSyntheticCapture Write synthetic IQ, metadata, and truth files.
%
%   paths = writeSyntheticCapture(outputDir, scenarioName, iq, meta, truth)
%   writes MAT, CSV metadata, and CSV truth artifacts for repeatable tests.

arguments
    outputDir (1,1) string
    scenarioName (1,1) string
    iq (:,1) {mustBeNumeric}
    meta struct
    truth table
end

if ~isfolder(outputDir)
    mkdir(outputDir);
end

baseName = matlab.lang.makeValidName(scenarioName);
iqPath = fullfile(outputDir, baseName + "_iq.mat");
metaPath = fullfile(outputDir, baseName + "_meta.json");
truthPath = fullfile(outputDir, baseName + "_truth.csv");

save(iqPath, "iq", "-v7.3");

jsonText = jsonencode(meta, "PrettyPrint", true);
fid = fopen(metaPath, "w");
if fid < 0
    error("writeSyntheticCapture:OpenFailed", ...
        "Could not write metadata file: %s", metaPath);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s", jsonText);
clear cleanup;

writetable(truth, truthPath);

paths = struct();
paths.IQ = iqPath;
paths.Metadata = metaPath;
paths.Truth = truthPath;

end

