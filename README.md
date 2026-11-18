# T4 permutation entropy for rTMS remission (TDBRAIN challenge)

Target: **remission**. One feature: permutation entropy (order 4, lag 1) of the eyes-closed
resting EEG at T4. Higher values predict remission after prefrontal rTMS.

Processing per recording:
- 19 channels from the .bdf file, demeaned, average reference
- 1-40 Hz band-pass (Butterworth, order 4, zero phase), resampled to 250 Hz
- 8 s non-overlapping segments; a segment is skipped if T4 goes beyond 150 uV or 250 uV peak-to-peak
- PE per segment, median over the kept segments

Held-out prediction: cases with T4 PE at or above the 65th percentile of the discovery MDD rTMS
patients (the top 35% in discovery) are flagged (Remitter = 1). The cutoff is fixed from discovery,
so a case's flag does not depend on the other held-out cases. The Responder column repeats the same
flags only so the template is complete.

## Files
- `run_t4pe.m` runs everything
- `t4pe.m` feature for one recording
- `perm_entropy.m` permutation entropy
- `read_bdf19.m` minimal BDF reader, no FieldTrip needed
- `check_reader.m` compares the reader with FieldTrip-converted copies of a few recordings
- `sha256file.m` checksums of the inputs, written to `output/input_sha256.csv`
- `METHODS.md` methods write-up (challenge format)

Needs MATLAB with the Signal Processing and Statistics toolboxes (tested on R2026a).
Point `TDBRAIN_ROOT` at the folder with `TDBRAIN_Dataset_V3_1`, `TDBRAIN_Tx_Response_ReplicationV3_1`,
`TDBRAIN_participants_V3.1.xlsx` and `TDBRAIN_Tx_Prediction_replicationV3.xlsx` (default: `./data`).
`check_reader.m` also needs `TDBRAIN_FT_DIR`, a folder of FieldTrip-converted recordings.
Output goes to `output/` (feature tables, prediction sheet, `run_info.txt`, checksums).
