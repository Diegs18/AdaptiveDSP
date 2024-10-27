function denom = nee_denom(B, N)
    freqB    = freqz(B, 1, N);
    denom = abs(freqB).^2;
end