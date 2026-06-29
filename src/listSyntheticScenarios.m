function scenarios = listSyntheticScenarios()
%listSyntheticScenarios Return Phase 3 synthetic validation scenarios.
%
%   scenarios = listSyntheticScenarios() returns a table describing the
%   controlled cases used to validate later receiver stages.

names = [
    "aligned_5gnb"
    "offset_100ns_one"
    "offset_250ns_one"
    "offset_500ns_one"
    "offset_1us_one"
    "offset_1p5us_one"
    "offset_3us_two"
    "offset_negative3us_one"
    "mixed_offsets"
    "common_mode_3us"
    "low_snr_250ns"
    "cfo_anomaly_one"
    "wrong_tdd_pattern"
    "wrong_special_slot"
    "single_gnb_not_assessable"
    "large_geometry_uncertainty"
    "rf_impairment_moderate"
    "rf_impairment_harsh"
    ];

purposes = [
    "Nominal aligned five-gNB relative PASS case"
    "Small sub-sample timing sensitivity case"
    "Warning-threshold timing sensitivity case"
    "Intermediate sub-microsecond timing offset"
    "Fail-threshold timing offset"
    "Strict private-network timing stress case"
    "Two gNBs hard-failed in opposite directions"
    "One gNB appears early by hard-fail threshold"
    "Several offsets in one capture"
    "Common-mode shift invisible to relative mode"
    "250 ns offset under low SNR"
    "One gNB has CFO anomaly"
    "Downlink energy appears in expected uplink slots"
    "Special slot has energy in guard/uplink portion"
    "Only one gNB visible, relative sync cannot be assessed"
    "Timing is measurable but geometry uncertainty is too large"
    "Moderate multipath, IQ imbalance, DC offset, phase noise, tone interference, and clipping"
    "Harsh RF impairment stress case for robustness characterization"
    ];

expected = [
    "PASS"
    "PASS_OR_WARNING"
    "WARNING_OR_SUSPECT"
    "SUSPECT"
    "FAIL"
    "FAIL"
    "FAIL"
    "FAIL"
    "MIXED"
    "RELATIVE_PASS_ABSOLUTE_UNKNOWN"
    "SUSPECT"
    "SUSPECT"
    "FAIL"
    "FAIL"
    "NOT_ASSESSABLE"
    "SUSPECT_OR_NOT_ASSESSABLE"
    "SHOULD_DETECT_MAIN_CELLS"
    "STRESS_CHARACTERIZATION"
    ];

scenarios = table(names, purposes, expected, ...
    'VariableNames', ["Name","Purpose","ExpectedReceiverOutcome"]);

end
