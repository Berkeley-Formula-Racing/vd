function fingerprint = fingerprintFile(path)
%FINGERPRINTFILE Return a lowercase SHA-256 fingerprint for one file.
path = string(path);
if ~isscalar(path) || strlength(strtrim(path)) == 0 || ~isfile(path)
    error('rampSpeed:assetNotFound', 'Asset file does not exist: %s', path);
end

fid = fopen(char(path),'r');
if fid < 0
    error('rampSpeed:assetReadFailed', 'Unable to read asset file: %s', path);
end
cleanup = onCleanup(@()fclose(fid));
bytes = fread(fid,Inf,'*uint8');

digest = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
digest.update(bytes);
digestBytes = typecast(digest.digest(),'uint8');
fingerprint = lower(string(reshape(dec2hex(digestBytes,2).',1,[])));
end
