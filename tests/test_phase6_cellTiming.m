% test_phase6_cellTiming
% Phase 6 smoke tests for per-PCI frame timing estimation.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot, "src")));

% Isolated cell: repeated SSB detections fold to approximately 0 ns phase.
[iqSingle, metaSingle, truthSingle] = generateSyntheticScenario("single_gnb_not_assessable");
pssSingle = estimatePSSCandidates(iqSingle, metaSingle);
sssSingle = detectSSSAndPCI(iqSingle, pssSingle, metaSingle);
[timingSingle, detTimingSingle] = estimateCellTiming(sssSingle, metaSingle);

assert(height(timingSingle) == 1, "Expected one cell timing row.");
assert(timingSingle.PCI(1) == truthSingle.PCI(1));
assert(abs(timingSingle.RelativeArrivalOffsetNs(1)) < 1e-6);
assert(timingSingle.FramePhaseNs(1) < 120, ...
    "Expected isolated synthetic cell phase near 0 ns.");
assert(~isempty(detTimingSingle));

% Two isolated captures combined at detection-table level prove the timing
% combiner can resolve relative offsets without multi-cell SSS contamination.
[detA, metaA, truthA] = isolatedDetectionForPCI(11, 0);
[detB, ~, truthB] = isolatedDetectionForPCI(12, 3000);
combined = [detA; detB];
[timingTwo, ~] = estimateCellTiming(combined, metaA);

assert(height(timingTwo) == 2, "Expected two PCI timing rows.");
rowA = timingTwo(timingTwo.PCI == truthA.PCI(1), :);
rowB = timingTwo(timingTwo.PCI == truthB.PCI(1), :);
assert(~isempty(rowA) && ~isempty(rowB));
relativeDifferenceNs = rowB.RelativeArrivalOffsetNs - rowA.RelativeArrivalOffsetNs;
assert(abs(relativeDifferenceNs - 3000) < 80, ...
    "Expected approximately 3 us relative timing difference.");

% Common-mode case at timing-table level: equal offsets should disappear in
% relative mode.
[detC, metaC, truthC] = isolatedDetectionForPCI(257, 3000);
[detD, ~, truthD] = isolatedDetectionForPCI(503, 3000);
timingCommon = estimateCellTiming([detC; detD], metaC);
assert(all(abs(timingCommon.RelativeArrivalOffsetNs) < 80), ...
    "Common-mode timing shift should not be visible in relative mode.");
assert(all([truthC.FrameOffsetNs; truthD.FrameOffsetNs] == 3000));

disp("Phase 6 cell timing tests passed.");

function [sssDetections, meta, truth] = isolatedDetectionForPCI(pci, offsetNs)
custom = struct();
custom.NumGNBs = 1;
custom.PCIs = pci;
custom.FrameOffsetsNs = offsetNs;
custom.CFOHz = 0;
custom.GNBPowerdB = 0;
custom.LocationUncertaintyM = 30;
custom.SiteDistanceM = 100;
custom.SNRdB = 35;

[iq, meta, truth] = generateSyntheticScenario("aligned_5gnb", custom);
pss = estimatePSSCandidates(iq, meta);
sssDetections = detectSSSAndPCI(iq, pss, meta);
sssDetections = sssDetections(sssDetections.IsUsable & sssDetections.PCI == pci, :);
assert(~isempty(sssDetections), "Expected usable isolated detection for PCI %d.", pci);
end

