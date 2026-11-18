%% compare read_bdf19 with the FieldTrip-converted copies of the same recordings
% (analysis_v3/sub-*/ses_1_EC_data_raw.mat: ft_preprocessing on the .bdf, no filtering)

root = getenv('TDBRAIN_ROOT');                 % same as in run_t4pe.m
if isempty(root), root = fullfile(pwd, 'data'); end
root = fullfile(root, 'TDBRAIN_Dataset_V3_1');
ftdir = getenv('TDBRAIN_FT_DIR');              % sub-*/ses_1_EC_data_raw.mat made with FieldTrip
if isempty(ftdir), ftdir = fullfile(pwd, 'fieldtrip_copies'); end
ids = {'sub-87999321', 'sub-88000181', 'sub-88000313', 'sub-88020917', 'sub-88056017'};
chans = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T3','C3','Cz','C4','T4', ...
         'T5','P3','Pz','P4','T6','O1','O2'};

for i = 1:numel(ids)
    fn = fullfile(root, ids{i}, 'ses-1', 'eeg', [ids{i} '_ses-1_task-restEC_eeg.bdf']);
    mf = fullfile(ftdir, ids{i}, 'ses_1_EC_data_raw.mat');
    if ~exist(fn, 'file') || ~exist(mf, 'file'), fprintf('%s: skipped\n', ids{i}); continue; end
    [x, fs] = read_bdf19(fn, chans);
    S = load(mf); d = S.data;
    [~, loc] = ismember(lower(chans), lower(d.label));
    xf = d.trial{1}(loc, :);
    n = min(size(x, 2), size(xf, 2));
    fprintf('%s: %g vs %g Hz, %d vs %d samples, max abs diff %.3g uV\n', ids{i}, fs, d.fsample, ...
        size(x, 2), size(xf, 2), max(abs(x(:, 1:n) - xf(:, 1:n)), [], 'all'));
end
