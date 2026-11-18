function [pe, nok, nseg] = t4pe(fn)
% Median permutation entropy at T4 for one eyes-closed recording.
% 19 channels, average reference, 1-40 Hz (4th order Butterworth, zero phase), 250 Hz,
% 8 s segments; a segment is skipped if T4 exceeds 150 uV or 250 uV peak-to-peak.

chans = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T3','C3','Cz','C4','T4', ...
         'T5','P3','Pz','P4','T6','O1','O2'};
[x, fs] = read_bdf19(fn, chans);

q = round(fs / 250);
if abs(fs / q - 250) > 1e-9, error('%s: %g Hz cannot be brought to 250 Hz', fn, fs); end

x = x - mean(x, 2);
x = x - mean(x, 1);                        % average reference

[z, p, k] = butter(4, [1 40] / (fs/2), 'bandpass');
[sos, g] = zp2sos(z, p, k);
y = filtfilt(sos, g, x(strcmp(chans, 'T4'), :)')';
if q > 1, y = resample(y, 1, q); end

L = 8 * 250;
nseg = floor(numel(y) / L);
h = nan(nseg, 1); ties = 0;
for s = 1:nseg
    seg = y((s-1)*L + (1:L));
    if max(abs(seg)) <= 150 && (max(seg) - min(seg)) <= 250
        [h(s), nt] = perm_entropy(seg, 4, 1);
        ties = ties + nt;
    end
end
nok = sum(~isnan(h));
pe = median(h, 'omitnan');
if nok < 10
    warning('t4pe:short', '%s: only %d usable 8 s segments', fn, nok);
end
if ties > 0
    warning('t4pe:ties', '%s: %d embedding vectors with tied values', fn, ties);
end
end
