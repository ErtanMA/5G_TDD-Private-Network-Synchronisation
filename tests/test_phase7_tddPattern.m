% test_phase7_tddPattern
% Phase 7 smoke tests for TDD pattern and special-slot envelope checking.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot, "src")));

[iqAligned, metaAligned, timingAligned] = timingForScenario("aligned_5gnb");
[tddAligned, ~, specialAligned, ~] = checkTDDPattern(iqAligned, metaAligned, timingAligned);
assert(all(tddAligned.TDDPatternStatus == "PASS"), ...
    "Expected aligned synthetic TDD pattern to pass.");
assert(all(specialAligned.SpecialSlotStatus == "PASS"), ...
    "Expected aligned synthetic special slot to pass.");

[iqWrongTDD, metaWrongTDD, timingWrongTDD] = timingForScenario("wrong_tdd_pattern");
[tddWrong, ~, specialWrongTDD, ~] = checkTDDPattern(iqWrongTDD, metaWrongTDD, timingWrongTDD);
assert(all(tddWrong.TDDPatternStatus == "FAIL"), ...
    "Expected wrong TDD pattern to fail.");
assert(all(specialWrongTDD.SpecialSlotStatus == "PASS"), ...
    "Wrong U-slot energy should not by itself fail the special-slot split.");

[iqWrongSpecial, metaWrongSpecial, timingWrongSpecial] = timingForScenario("wrong_special_slot");
[tddSpecial, ~, specialWrong, specialMeas] = checkTDDPattern(iqWrongSpecial, metaWrongSpecial, timingWrongSpecial);
assert(all(specialWrong.SpecialSlotStatus == "FAIL"), ...
    "Expected wrong special-slot energy placement to fail.");
assert(~isempty(specialMeas), "Expected special-slot measurements.");
assert(all(tddSpecial.TDDPatternStatus == "PASS" | tddSpecial.TDDPatternStatus == "SUSPECT"), ...
    "Wrong special slot should be represented primarily in special-slot status.");

[tddEmpty, slotEmpty, specialEmpty, specialMeasEmpty] = checkTDDPattern(iqAligned, metaAligned, table());
assert(isempty(tddEmpty) && isempty(slotEmpty) && isempty(specialEmpty) && isempty(specialMeasEmpty));

disp("Phase 7 TDD pattern tests passed.");

function [iq, meta, timing] = timingForScenario(scenarioName)
custom = struct();
custom.NumGNBs = 1;
custom.PCIs = 11;
custom.FrameOffsetsNs = 0;
custom.CFOHz = 0;
custom.GNBPowerdB = 0;
custom.LocationUncertaintyM = 30;
custom.SiteDistanceM = 100;
custom.SNRdB = 35;

[iq, meta, ~] = generateSyntheticScenario(scenarioName, custom);
pss = estimatePSSCandidates(iq, meta);
sss = detectSSSAndPCI(iq, pss, meta);
timing = estimateCellTiming(sss, meta);
assert(~isempty(timing), "Expected timing for scenario %s.", scenarioName);
end

