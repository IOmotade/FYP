function V = fVoltageSourceSignals(N, vs_mag, nsamp)

extra = log2(nsamp)-N;

H = fOFDMTransformMtrx(nsamp);
V = vs_mag*H(2:end-extra, :);

end