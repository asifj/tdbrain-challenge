function [x, fs] = read_bdf19(fn, chans)
% Read a TDBRAIN .bdf file and return the requested channels in uV (channels x samples).
% Plain fread so it works without FieldTrip mex files (Apple Silicon).

alias = {'T7','T3'; 'T8','T4'; 'P7','T5'; 'P8','T6'};

fid = fopen(fn, 'r', 'ieee-le');
if fid < 0, error('cannot open %s', fn); end
c = onCleanup(@() fclose(fid));

hdr  = fread(fid, 256, '*char')';
hlen = str2double(hdr(185:192));
nrec = str2double(hdr(237:244));
dur  = str2double(hdr(245:252));
ns   = str2double(hdr(253:256));

rd = @(w) strtrim(cellstr(fread(fid, [w ns], '*char')'));
lab  = rd(16);
rd(80);                                    % transducer
dim  = rd(8);
pmin = str2double(rd(8)); pmax = str2double(rd(8));
dmin = str2double(rd(8)); dmax = str2double(rd(8));
rd(80);                                    % prefilter
nsmp = str2double(rd(8));

for k = 1:size(alias, 1)
    lab(strcmpi(lab, alias{k,1})) = alias(k,2);
end

fseek(fid, hlen, 'bof');
raw = fread(fid, inf, 'bit24=>double');
recLen = sum(nsmp);
nr = floor(numel(raw) / recLen);
if nrec > 0, nr = min(nr, nrec); end
raw = reshape(raw(1:recLen*nr), recLen, nr);
off = [0; cumsum(nsmp(:))];

x = [];
for j = 1:numel(chans)
    i = find(strcmpi(lab, chans{j}), 1);
    if isempty(i), error('%s: no channel %s', fn, chans{j}); end
    if j == 1
        n0 = nsmp(i); fs = n0 / dur; x = zeros(numel(chans), n0 * nr);
    elseif nsmp(i) ~= n0
        error('%s: %s has a different sampling rate', fn, chans{j});
    end
    v = reshape(raw(off(i)+1:off(i+1), :), 1, []);
    v = (v - dmin(i)) * (pmax(i) - pmin(i)) / (dmax(i) - dmin(i)) + pmin(i);
    switch lower(strrep(dim{i}, char(181), 'u'))
        case 'uv'
        case 'mv', v = v * 1e3;
        case 'v',  v = v * 1e6;
        otherwise, error('%s: unexpected unit "%s" for %s', fn, dim{i}, chans{j});
    end
    x(j, :) = v;
end
end
