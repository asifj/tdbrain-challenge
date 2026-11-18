function h = sha256file(fn)
md = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(fn, 'r');
while true
    b = fread(fid, 2^22, '*uint8');
    if isempty(b), break; end
    md.update(typecast(b, 'int8'));
end
fclose(fid);
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'))', 1, []));
end
