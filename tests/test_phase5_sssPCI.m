% test_phase5_sssPCI
% Phase 5 smoke tests for CFO estimation and SSS/PCI recovery.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot, "src")));

% Isolated cell: recover PCI from PSS candidates via SSS.
[iqSingle, metaSingle, truthSingle] = generateSyntheticScenario("single_gnb_not_assessable");
[pssSingle, ~] = estimatePSSCandidates(iqSingle, metaSingle);
[sssSingle, dbgSingle] = detectSSSAndPCI(iqSingle, pssSingle, metaSingle);

assert(~isempty(sssSingle), "Expected SSS/PCI detections for isolated cell.");
assert(dbgSingle.NumProcessed > 0);
usableSingle = sssSingle(sssSingle.IsUsable, :);
assert(~isempty(usableSingle), "Expected at least one usable SSS detection.");
assert(usableSingle.PCI(1) == truthSingle.PCI(1), ...
    "Expected strongest usable PCI to match truth.");

% CFO case: recover PCI and estimate CFO direction/magnitude on an isolated signal.
custom = struct();
custom.NumGNBs = 1;
custom.PCIs = 257;
custom.FrameOffsetsNs = 0;
custom.CFOHz = 250;
custom.GNBPowerdB = 0;
custom.LocationUncertaintyM = 30;
custom.SiteDistanceM = 100;
custom.SNRdB = 35;

[iqCFO, metaCFO, truthCFO] = generateSyntheticScenario("aligned_5gnb", custom);
[pssCFO, ~] = estimatePSSCandidates(iqCFO, metaCFO);
[sssCFO, ~] = detectSSSAndPCI(iqCFO, pssCFO, metaCFO);
usableCFO = sssCFO(sssCFO.IsUsable, :);

assert(~isempty(usableCFO), "Expected usable SSS detection in CFO scenario.");
assert(usableCFO.PCI(1) == truthCFO.PCI(1), ...
    "Expected PCI recovery in CFO scenario.");
matchingCFO = usableCFO(usableCFO.PCI == truthCFO.PCI(1) & ...
    usableCFO.IsSSBIndexUsable, :);
assert(~isempty(matchingCFO), ...
    "Expected PCI-matched detections with usable PBCH DM-RS SSB index.");
assert(abs(median(matchingCFO.CFOHz) - truthCFO.CFOHz(1)) < 80, ...
    "Expected median per-cell CFO estimate to be within 80 Hz of truth.");
assert(all(matchingCFO.SSBIndex == truthCFO.SSBIndex(1)), ...
    "Expected recovered SSB index to match truth.");

% Candidate count guard: processing fewer candidates should still recover the
% strongest isolated cell when the strongest PSS candidate is processed.
[sssLimited, dbgLimited] = detectSSSAndPCI(iqSingle, pssSingle, metaSingle, ...
    struct("MaxCandidatesToProcess", 4));
usableLimited = sssLimited(sssLimited.IsUsable, :);
assert(dbgLimited.NumProcessed <= 4);
assert(~isempty(usableLimited));
assert(any(usableLimited.PCI == truthSingle.PCI(1)));

disp("Phase 5 SSS/PCI tests passed.");
