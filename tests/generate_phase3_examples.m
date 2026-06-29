% generate_phase3_examples
% Optional helper to write all synthetic scenarios under data_examples.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot, "src")));

scenarioTable = listSyntheticScenarios();
outDir = fullfile(projectRoot, "data_examples", "synthetic_phase3");
if ~isfolder(outDir)
    mkdir(outDir);
end

for k = 1:height(scenarioTable)
    scenarioName = scenarioTable.Name(k);
    fprintf("Generating %s...\n", scenarioName);
    [iq, meta, truth] = generateSyntheticScenario(scenarioName);
    writeSyntheticCapture(outDir, scenarioName, iq, meta, truth);
end

fprintf("Wrote synthetic scenarios to:\n  %s\n", outDir);

