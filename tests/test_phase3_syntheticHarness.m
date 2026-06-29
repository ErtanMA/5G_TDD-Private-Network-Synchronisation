% test_phase3_syntheticHarness
% Phase 3 smoke tests for synthetic validation scenario generation.

clear; clc;

projectRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(projectRoot, "src")));

scenarios = listSyntheticScenarios();
assert(height(scenarios) >= 12, "Expected a broad scenario set.");
assert(any(scenarios.Name == "offset_250ns_one"));
assert(any(scenarios.Name == "common_mode_3us"));
assert(any(scenarios.Name == "wrong_tdd_pattern"));

[iqAligned, metaAligned, truthAligned] = generateSyntheticScenario("aligned_5gnb");
assert(isvector(iqAligned) && ~isreal(iqAligned));
assert(metaAligned.SampleRateHz == 30.72e6);
assert(height(truthAligned) == 5);
assert(all(truthAligned.ExpectedRelativeOffsetNs == 0));
assert(all(truthAligned.ExpectedTDDStatus == "PASS"));

[~, ~, truth250] = generateSyntheticScenario("offset_250ns_one");
assert(any(truth250.ExpectedRelativeOffsetNs == 250));

[~, ~, truth100] = generateSyntheticScenario("offset_100ns_one");
assert(any(truth100.ExpectedRelativeOffsetNs == 100));

[~, ~, truth3us] = generateSyntheticScenario("offset_3us_two");
assert(any(truth3us.ExpectedRelativeOffsetNs == 3000));
assert(any(truth3us.ExpectedRelativeOffsetNs == -3000));

[~, ~, truthCommon] = generateSyntheticScenario("common_mode_3us");
assert(all(truthCommon.FrameOffsetNs == 3000));
assert(all(truthCommon.ExpectedRelativeOffsetNs == 0));

[~, ~, truthTDD] = generateSyntheticScenario("wrong_tdd_pattern");
assert(all(truthTDD.ExpectedTDDStatus == "FAIL"));

[~, ~, truthSpecial] = generateSyntheticScenario("wrong_special_slot");
assert(all(truthSpecial.ExpectedSpecialSlotStatus == "FAIL"));

[~, ~, truthSingle] = generateSyntheticScenario("single_gnb_not_assessable");
assert(height(truthSingle) == 1);

tmpRoot = tempname;
mkdir(tmpRoot);
cleanup = onCleanup(@() rmdir(tmpRoot, "s"));
paths = writeSyntheticCapture(tmpRoot, "offset_250ns_one", iqAligned, metaAligned, truthAligned);
assert(isfile(paths.IQ));
assert(isfile(paths.Metadata));
assert(isfile(paths.Truth));

disp("Phase 3 synthetic harness tests passed.");

