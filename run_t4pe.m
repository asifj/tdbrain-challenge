%% T4 permutation entropy: TDBRAIN rTMS remission prediction
% Computes T4 PE for the MDD rTMS discovery patients and the 104 held-out cases. The cutoff is
% the 65th percentile of T4 PE in the discovery patients (top 35%); held-out cases at or above
% it are flagged. Also checks the rule on discovery and fills in the prediction template.
% Target: remission. The Responder column gets the same flags so the template is complete.

% folder holding the TDBRAIN v3.1 data, participants file and template:
% set TDBRAIN_ROOT or edit the default below
root  = getenv('TDBRAIN_ROOT');
if isempty(root), root = fullfile(pwd, 'data'); end
disc  = fullfile(root, 'TDBRAIN_Dataset_V3_1');
held  = fullfile(root, 'TDBRAIN_Tx_Response_ReplicationV3_1');
parts = fullfile(root, 'TDBRAIN_participants_V3.1.xlsx');
tmpl  = fullfile(root, 'TDBRAIN_Tx_Prediction_replicationV3.xlsx');
out   = fullfile(pwd, 'output');
pct   = 65;                                 % cutoff percentile of discovery T4 PE

if ~exist(out, 'dir'), mkdir(out); end
ecfile = @(d, id) fullfile(d, id, 'ses-1', 'eeg', [id '_ses-1_task-restEC_eeg.bdf']);

%% discovery: MDD patients treated with prefrontal rTMS (protocol 1-3)
P = readtable(parts, 'VariableNamingRule', 'preserve');
for v = {'sessID', 'rTMS PROTOCOL', 'Responder', 'Remitter'}
    if iscell(P.(v{1})), P.(v{1}) = str2double(P.(v{1})); end
end
P = P(P.sessID == 1, :);
[~, iu] = unique(P.TDBRAIN_ID, 'stable'); P = P(iu, :);
prot = P.('rTMS PROTOCOL');
D = P(strcmp(P.formal_status, 'MDD') & ismember(prot, 1:3) & ~isnan(P.Responder), :);

D.T4PE = nan(height(D), 1); D.nok = nan(height(D), 1); files = {};
for i = 1:height(D)
    fn = ecfile(disc, D.TDBRAIN_ID{i});
    if ~exist(fn, 'file'), fprintf('missing %s\n', fn); continue; end
    [D.T4PE(i), D.nok(i)] = t4pe(fn);
    files{end+1, 1} = fn; %#ok<SAGROW>
end
writetable(D(:, {'TDBRAIN_ID', 'Responder', 'Remitter', 'T4PE', 'nok'}), fullfile(out, 'discovery_T4PE.csv'));

cutoff = prctile(D.T4PE(~isnan(D.T4PE)), pct);
fprintf('\ncutoff (discovery %dth percentile) %.6f\n', pct, cutoff);
fprintf('discovery, MDD rTMS (T4 PE >= cutoff = predicted remitter/responder)\n');
for v = {'Remitter', 'Responder'}
    ok = ~isnan(D.(v{1})) & ~isnan(D.T4PE);
    [a, nppv, sens, fl] = score_check(D.T4PE(ok), D.(v{1})(ok), cutoff);
    fprintf('  %-9s n=%d  AUC %.3f  flagged %.1f%%  nPPV %.1f%%  sensitivity %.1f%%\n', v{1}, nnz(ok), a, fl, nppv, sens);
end

%% held-out cases
T = readtable(tmpl, 'VariableNamingRule', 'preserve');
pe = nan(height(T), 1); nok = pe;
for i = 1:height(T)
    fn = ecfile(held, sprintf('sub-%d', T.ID(i)));
    if ~exist(fn, 'file'), fprintf('missing %s\n', fn); continue; end
    [pe(i), nok(i)] = t4pe(fn);
    files{end+1, 1} = fn; %#ok<SAGROW>
end
writetable(table(T.ID, pe, nok, 'VariableNames', {'ID', 'T4PE', 'nok'}), fullfile(out, 'heldout_T4PE.csv'));

% flag cases at or above the discovery cutoff; same flags for both outcome columns
flag = double(pe >= cutoff);
flag(isnan(pe)) = NaN;
T.Responder = flag; T.Remitter = flag;
xfn = fullfile(out, 'TDBRAIN_Tx_Prediction_replicationV3_T4PE.xlsx');
writetable(T, xfn);
fprintf('\nheld-out: %d of %d flagged, %d without a value\n', nnz(flag == 1), height(T), nnz(isnan(pe)));

% the sheet must match the template apart from the two prediction columns
T0 = readtable(tmpl, 'VariableNamingRule', 'preserve');
T1 = readtable(xfn, 'VariableNamingRule', 'preserve');
keep = setdiff(T0.Properties.VariableNames, {'Responder', 'Remitter'}, 'stable');
assert(isequal(T0.Properties.VariableNames, T1.Properties.VariableNames), 'column names changed');
assert(isequaln(T0(:, keep), T1(:, keep)), 'template cells changed');
assert(all(ismember(T1.Remitter(~isnan(T1.Remitter)), [0 1])), 'flags must be 0/1');
fprintf('template check passed\n');

%% record what was run
files = [{parts; tmpl}; files];
h = cellfun(@sha256file, files, 'UniformOutput', false);
writetable(table(files, h, 'VariableNames', {'file', 'sha256'}), fullfile(out, 'input_sha256.csv'));
fid = fopen(fullfile(out, 'run_info.txt'), 'w');
fprintf(fid, 'run %s\nMATLAB %s\ncutoff %.10f (discovery %dth percentile)\n', datestr(now, 'yyyy-mm-dd HH:MM'), version, cutoff, pct);
for tb = {'MATLAB', 'Signal Processing Toolbox', 'Statistics and Machine Learning Toolbox'}
    vv = ver; vv = vv(strcmp({vv.Name}, tb{1}));
    if ~isempty(vv), fprintf(fid, '%s %s\n', vv.Name, vv.Version); end
end
fprintf(fid, '%s\nprediction sheet sha256 %s\n', computer, sha256file(xfn));
fclose(fid);


function [auc, nppv, sens, fl] = score_check(s, y, cutoff)
y = y(:) == 1; f = s(:) >= cutoff;
r = tiedrank(s);
auc = (sum(r(y)) - nnz(y) * (nnz(y) + 1) / 2) / (nnz(y) * nnz(~y));
tp = nnz(y & f);
nppv = 100 * (tp / nnz(f)) / mean(y);
sens = 100 * tp / nnz(y);
fl = 100 * mean(f);
end
