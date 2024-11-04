%% Author Nicholas DiGregorio, 1220871392
function num = nee_num(B, Bhat, N)
    freqB    = freqz(B, 1, N);
    freqBhat = freqz(Bhat, 1, N);
    num = abs(freqB - freqBhat).^2;
end