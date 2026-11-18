---
title: "Methods: T4 permutation entropy as a predictor of remission after rTMS (TDBRAIN Prediction Challenge, prognostic track)"
author: "Asif Jamil, PhD. Massachusetts General Hospital and Harvard Medical School, Boston, MA, USA"
---

*Sections 1-9 and the list of non-recoverable details were generated with the standard methods prompt (Rules v1.1, §5.5) using Claude (Anthropic), from the submitted MATLAB code: `run_t4pe.m`, `t4pe.m`, `perm_entropy.m`, `read_bdf19.m`, `sha256file.m`, `check_reader.m` and `README.md` (final version). The author's note at the end was written separately and is not derived from the code.*

## Summary of key methodological choices

- **Target:** remission (rTMS, MDD), stated in the code header and README. The prediction is a binary flag per held-out case. The code writes the same flags to both the `Remitter` and `Responder` columns so the template is complete; the submission targets **remission**.
- **Predictor:** one EEG feature, the median normalised permutation entropy (PE; order 4, delay 1 sample) of the eyes-closed resting EEG at channel T4, after a 19-channel average reference, a 1-40 Hz band-pass and resampling to 250 Hz, computed over non-overlapping 8 s segments that pass amplitude limits.
- **No model is trained.** There are no fitted parameters, no feature selection, no scaling or normalisation against training data, no cross-validation and no random numbers.
- **Decision rule:** a fixed cutoff equal to the 65th percentile (`pct = 65`, MATLAB `prctile`) of T4 PE in the discovery MDD rTMS patients (0.629213 in the submitted run). Held-out cases with T4 PE at or above the cutoff are flagged as predicted remitters. Higher PE = predicted remission.
- **Discovery data** provide the cutoff (from feature values only, no outcome labels) and in-sample performance metrics for the rule. The cutoff is the only quantity carried from discovery to the held-out predictions.
- **No demographic or clinical variable enters the prediction.** The participants table is used only to select the discovery evaluation cohort (session 1, `formal_status` = MDD, `rTMS PROTOCOL` in {1, 2, 3}) and to read its outcome labels.

## 1. Prediction task and target definition

Binary prognostic prediction of remission after rTMS treatment for MDD in the 104 cases of the TDBRAIN held-out Treatment Prediction dataset V3.1. For each case the output is 1 (predicted remitter) or 0 (predicted non-remitter); cases without a feature value would receive NaN (none occurred). The `Responder` column is filled with the identical flags; the code makes no separate response prediction. Outcome labels used for the discovery check are the `Remitter` and `Responder` columns of `TDBRAIN_participants_V3.1.xlsx` (definitions as provided by TDBRAIN; the code does not derive them).

## 2. Input data and samples used

**Recordings.** Eyes-closed resting-state EEG, session 1, BDF format: `<dataset>/<ID>/ses-1/eeg/<ID>_ses-1_task-restEC_eeg.bdf`. Eyes-open recordings and other sessions are not used.

**Discovery evaluation cohort** (`run_t4pe.m`): rows of `TDBRAIN_participants_V3.1.xlsx` with `sessID == 1`; duplicate `TDBRAIN_ID` rows reduced to the first occurrence (`unique(..., 'stable')`); then `formal_status == 'MDD'` AND `rTMS PROTOCOL` ∈ {1, 2, 3} AND `Responder` not missing. Recordings are read from `TDBRAIN_Dataset_V3_1`. If a file is missing the case is skipped (message printed, feature left NaN). For each outcome, metrics use the cases with a non-missing label and a non-missing feature. In the run that produced the submission this gave n = 160 for both outcomes (printed output).

**Held-out cases:** all rows of `TDBRAIN_Tx_Prediction_replicationV3.xlsx` (n = 104), recording `TDBRAIN_Tx_Response_ReplicationV3_1/sub-<ID>/ses-1/eeg/sub-<ID>_ses-1_task-restEC_eeg.bdf`.

**Splitting.** There is no train/validation/test split. Discovery and held-out recordings are processed independently with the same fixed function (`t4pe.m`); the discovery T4 PE values of all 160 cohort members (regardless of outcome) define the cutoff.

## 3. Preprocessing

All steps are per recording; none uses information from other recordings.

1. **Reading** (`read_bdf19.m`). BDF header parsed directly (header length, number of records, record duration, number of signals; per signal: label, physical dimension, physical min/max, digital min/max, samples per record). Data read as little-endian signed 24-bit integers (`fread(..., 'bit24=>double')`). The number of complete records is `floor(total samples / samples per record)`, capped at the header's record count if that is positive (an incomplete final record is dropped). Values are scaled as `(v - digital_min) * (physical_max - physical_min) / (digital_max - digital_min) + physical_min`; a unit of `uV`/`µV` is kept, `mV` is multiplied by 1e3 and `V` by 1e6, and any other unit stops the run with an error. Labels T7, T8, P7, P8 are renamed T3, T4, T5, T6; channels are matched case-insensitively (first match); a missing channel is an error. Sampling rate = samples per record of the first selected channel ÷ record duration (500 Hz for TDBRAIN); a selected channel with a different number of samples per record is an error.
2. **Channels:** the 19 channels Fp1, Fp2, F7, F3, Fz, F4, F8, T3, C3, Cz, C4, T4, T5, P3, Pz, P4, T6, O1, O2. Other EEG and all auxiliary channels are ignored.
3. **Demeaning:** each channel's mean over the whole recording is subtracted.
4. **Re-referencing:** common average of the 19 channels (subtract the across-channel mean at every sample).
5. **Filtering:** `butter(4, [1 40]/(fs/2), 'bandpass')` (zero-pole-gain form, converted with `zp2sos`), applied forward and backward with `filtfilt(sos, g, ·)` (zero phase). Only the re-referenced T4 signal is filtered and kept from this point on.
6. **Resampling:** `resample(y, 1, q)` with `q = round(fs/250)` (q = 2 for 500 Hz data; skipped if q = 1), MATLAB's default anti-aliasing filter. Before filtering the code checks that fs/q equals 250 Hz and stops with an error otherwise.
7. **Segmentation:** non-overlapping segments of L = 2000 samples (8 s at 250 Hz) from the start of the recording; `nseg = floor(samples / 2000)`; the remainder is discarded.
8. **Artifact handling:** a segment is kept only if `max(|T4|) ≤ 150` and `max(T4) − min(T4) ≤ 250` (µV); otherwise it is excluded. No ICA, channel interpolation, EOG/EMG regression or other cleaning is done.
9. **Missing data:** a recording with no kept segments gets NaN (median with `'omitnan'` over an all-NaN vector). A warning (no exclusion) is printed if fewer than 10 segments (80 s) are kept. In the submitted run every recording had at least 12 kept segments in discovery and 13 in the held-out set (column `nok` of the output CSVs), and no held-out value was missing.
10. **Normalisation/scaling:** none. The only quantity derived from discovery data is the cutoff (step in §5).

## 4. Features/predictors

One feature per recording: **T4 permutation entropy** (`perm_entropy.m`, `t4pe.m`).

- For each kept segment s (2000 samples): embedding dimension m = 4, delay τ = 1, giving n = 2000 − 3 = 1997 overlapping vectors `[s(i), s(i+1), s(i+2), s(i+3)]`.
- Each vector is sorted ascending (`sort(E, 2)`; MATLAB's sort keeps equal values in their original order). The number of vectors containing tied values is counted, and `t4pe.m` prints a warning for any recording with ties (none were reported in the submitted run). The permutation indices are coded as `code = Σ_k (idx_k − 1) · 4^(k−1)`, k = 1..4 (256 possible codes, of which the 24 permutations occur).
- Pattern probabilities p = counts / n (codes with zero count dropped); PE = −Σ p·ln(p) / ln(4!) (normalised to 0..1).
- Recording value = median of the segment PEs over kept segments.

There is no feature selection in the code; the single feature is fixed.

## 5. Model architecture/algorithm

No statistical or machine-learning model. The prediction is a threshold rule on the single feature: T4 PE ≥ cutoff = predicted remitter, where cutoff = `prctile(discovery T4 PE, 65)` over the 160 discovery cohort members with a value. Fixed parameters in the code: direction (≥), percentile `pct = 65`.

## 6. Training procedure

None beyond computing the cutoff as a percentile of discovery feature values (no outcome labels involved): no loss, optimisation, class weighting, resampling, early stopping, random seeds or ensembling. The procedure is deterministic.

## 7. Hyperparameter tuning and model selection

The code contains no tuning, model selection or cross-validation. All parameters are hard-coded: channel T4; 19-channel average reference; 1-40 Hz, order-4 Butterworth; 250 Hz; 8 s segments; amplitude limits 150 µV and 250 µV peak-to-peak; PE order 4, delay 1; median over segments; direction (higher = remission); `pct = 65`. How these values were chosen is not recorded in the code (see the final section).

**Data leakage:** no preprocessing, feature selection, dimensionality reduction or tuning step is fitted inside or outside a cross-validation loop in the code. The cutoff is computed from discovery feature values without outcome labels; held-out outcome labels are never available. Each held-out flag depends only on that case's recording and the fixed cutoff, not on the other held-out cases.

## 8. Predictions and evaluation

**Held-out predictions.** T4 PE is computed for each of the 104 cases; a case gets 1 if T4 PE ≥ cutoff, otherwise 0; a case with NaN would be set to NaN. The flags are written to both `Responder` and `Remitter`; all other template columns are unchanged. Output: `output/TDBRAIN_Tx_Prediction_replicationV3_T4PE.xlsx` (36 of 104 flagged, 0 missing). The sheet is then read back and the code asserts that column names and every non-prediction cell are identical to the organiser template and that all flags are 0 or 1.

**Discovery metrics** (`score_check` in `run_t4pe.m`, in sample, same rule):

- AUC = (Σ ranks of positives − n₊(n₊+1)/2) / (n₊ n₋), ranks by `tiedrank` (average ranks for ties).
- Flagged = T4 PE ≥ cutoff (the same rule as for the held-out cases). nPPV = 100 × (TP / number flagged) / base rate, TP = positives among the flagged, base rate = proportion of positives.
- Sensitivity = 100 × TP / n₊.

Printed results of the submitted run (MDD rTMS discovery cohort):

| Outcome | n | AUC | Flagged | nPPV | Sensitivity |
|---|---|---|---|---|---|
| Remitter | 160 | 0.624 | 35.0% | 135.1% | 47.3% |
| Responder | 160 | 0.603 | 35.0% | 122.9% | 43.0% |

Intermediate outputs: `output/discovery_T4PE.csv` (TDBRAIN_ID, Responder, Remitter, T4PE, number of kept segments) and `output/heldout_T4PE.csv` (ID, T4PE, number of kept segments). Provenance: `output/input_sha256.csv` (SHA-256 of the participants file, the template and every EEG file read) and `output/run_info.txt` (date, MATLAB and toolbox versions, cutoff value, platform, SHA-256 of the prediction sheet).

## 9. Software/packages and versions

MATLAB (README: tested on R2026a; exact versions written to `run_info.txt` at run time) with the Signal Processing Toolbox (`butter`, `zp2sos`, `filtfilt`, `resample`) and the Statistics and Machine Learning Toolbox (`tiedrank`); base MATLAB `readtable`/`writetable` for the spreadsheets; Java `MessageDigest` (via MATLAB) for checksums. No EEG toolbox is used in the pipeline (own BDF reader; a code comment notes this avoids FieldTrip MEX files on Apple Silicon). `check_reader.m` is a separate check that compares `read_bdf19` with FieldTrip-converted copies (`ft_preprocessing`, no filtering) of five discovery recordings and prints sampling rate, sample count and maximum absolute difference. Data locations are read from the environment variable `TDBRAIN_ROOT` (folder containing `TDBRAIN_Dataset_V3_1`, `TDBRAIN_Tx_Response_ReplicationV3_1`, the participants file and the template; default `./data`) and, for `check_reader.m` only, `TDBRAIN_FT_DIR` (FieldTrip-converted copies; default `./fieldtrip_copies`); outputs go to `./output`. No parallel processing, GPU or random seeds.

## Replication-critical details not recoverable from code

- How the fixed parameters (channel T4, PE order 4 and delay 1, 1-40 Hz band, 8 s segments, 150/250 µV limits, median, direction "higher = remission", 65th-percentile cutoff) were chosen: **Not specified / cannot be determined from submitted code.**
- The code states the target (remission) only in comments and the README; it writes identical flags to both outcome columns.
- Exact TDBRAIN file versions beyond the names used in the code (participants file V3.1, dataset folders `..._V3_1`): **Not specified / cannot be determined from submitted code** (checksums of the files used are written to `output/input_sha256.csv` at run time).
- Exact MATLAB/toolbox versions and platform are not in the code; they are recorded in `output/run_info.txt` when the code runs.
- Details of MATLAB's default `resample` anti-aliasing filter depend on the MATLAB version: **Not specified in code.**

---

## Author's note: how the feature and cutoff were chosen, and supporting checks (not generated from code)

**Selection history.** The feature was not fixed before looking at outcomes, and the discovery numbers above should be read as exploratory. Right temporal entropy measures were already among the best entries on the public leaderboard when we started. We implemented a family of entropy measures (permutation entropy of orders 3-5 and delays 1, 2 and 4 samples at all 19 channels, plus ordinal-pattern complexity, attention entropy, sample entropy and Lempel-Ziv complexity) and looked at their association with outcome in the labelled discovery subjects. T4 permutation entropy of order 4, delay 1 was kept from that comparison. We also ran a much larger exploratory search over about 11,000 spectral and entropy features; nothing in it beat the same search run on shuffled labels, and the submitted feature was not taken from it. Remission was chosen as the target because the response numbers were weaker.

The cutoff changed during preparation. A first version flagged the top 40% within the held-out set. After two external reviews we tested, on the discovery cohort, an age/sex normative version of the feature, an average over delays 1, 2 and 4, a combination with ordinal-pattern statistical complexity, flagged fractions of 33, 35 and 40%, and a fixed versus a within-batch cutoff. The normative and multi-delay versions performed the same as the plain feature (AUC 0.630 and 0.623 vs 0.624); adding statistical complexity made it worse (AUC 0.618). We kept the plain feature and moved to a fixed cutoff at the discovery 65th percentile, which gives a higher nPPV than 40% (135% vs 132%) while keeping sensitivity near 47%, and makes each case's flag independent of the rest of the held-out batch.

**Validation of the selection.** To estimate how much of the in-sample result comes from these choices, we repeated them inside cross-validation on the 160 MDD rTMS discovery patients (50 repeats of 5-fold CV, Python). Choosing the best of 342 candidates (PE order 3-5 × delay 1/2/4 × 19 channels × either direction) and the best flag fraction (30-50%, sensitivity at least 41%) inside each training fold gave an out-of-fold remission nPPV of 121.7% (5th-95th percentile over repeats 107-136%), sensitivity 42%; a T4 variant was chosen in at least 84% of training folds. Choosing among the 21 options of the second round (the six variants above × 33/35/40%) gave 128.5% (121-137%), sensitivity 47%. Our honest expectation for the held-out set is therefore about 120-130% rather than the in-sample 135%. For response the first procedure gave 103.5%, so we make no claim for response.

**Artifact checks** (discovery MDD rTMS patients, remission):

- *Reference:* the same feature with the recorded linked-mastoid reference instead of the average reference gives AUC 0.617 (vs 0.624), so contamination spread through the average reference does not explain the result.
- *Eye movement:* frontal (Fp1/Fp2) 1-4 Hz relative power, a proxy for ocular activity, does not predict remission (AUC 0.48), and T4 PE keeps its accuracy after regressing it out (AUC 0.63).
- *Amount of usable data:* the fraction of segments kept is unrelated to T4 PE (r = −0.03) and to outcome (AUC 0.52).
- *Muscle:* this is the main open question. Recordings with high T4 PE have much more 30-45 Hz power at T4 and T3 (top vs bottom quartile, +0.76 and +0.54 log10 units) but not at O1 or Cz, and less 4-13 Hz power everywhere. The outcome association sits in the temporal 30-45 Hz part (AUC 0.61), not in the slow-wave part (0.53). Ordinal-pattern statistical complexity at T4, which is lower for noise-like signals, is also lower in remitters (AUC 0.40). Activity of the temporalis muscle is therefore a plausible contributor. In earlier analyses of the full labelled set, PE from only the quietest half of each recording's segments correlated 0.97 with the full value and predicted nearly as well (AUC 0.59 vs 0.60), which argues against intermittent muscle bursts but not against steady muscle tone. Scalp EEG cannot fully separate the two, and we welcome the organisers' artifact-dependence checks.

**Cutoff on the held-out set.** The cutoff (0.629213) is fixed from discovery and flags 36 of the 104 held-out cases. A within-batch rule at the same level (top 35%) would flag 37, almost the same cases. Because only MDD cases are scored, we checked, without changing the rule, how the flags fall across the template's rTMS protocol column: about 42% of cases with protocol 1 or 2 (n = 53) are flagged.

**Reader check.** `read_bdf19` reproduces FieldTrip's reading of five discovery recordings to within 4 × 10⁻¹² µV, with identical sampling rates and sample counts. An independent Python implementation of the whole pipeline gives the same T4 PE values (maximum difference < 0.001, rank correlation > 0.9999).
