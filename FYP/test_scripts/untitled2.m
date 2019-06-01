%% Generate accurate frequency components
tmpVs = [5; zeros(N-1, 1)];
tmpMemR = storedMemR;
[corr_fmag, corr_fval] = deal(zeros(N));
corr_Vbias = zeros(N);
J = fShiftingMatrix(N);
for i = 1:N
    for j = 1:N
        tmpCircuit = fMacSpiceSim(N, (J^(i-1))*tmpVs, tmpMemR, LRowR, LColR);
        im = tmpCircuit.IO.value(end, :);
        corr_fmag(i, j) = im(j);
        corr_fval(i, j) = fsource(i);
    end
    corr_Vbias(i, :) = tmpCircuit.VI.value(i, :);
end

%%
Corr_Frequency_Components = fUnits(corr_fmag , 'A')
round(mean(round((fmag./corr_fmag).^(1/2), 1, 'significant'), 2), 1, 'significant')
%%
freqVal = fGetFieldValues(Circuit.FreqDom, 'IO', 'value');
[meas_fmag, fval] = deal(zeros(N));
filt_tol = 10;
for colIdx = 1:N
    figure;
    plot(f, freqSpectrum, 'r'); hold on;
    freqSpectrum(abs(freqSpectrum)<max(abs(freqSpectrum))/filt_tol) = 0; %filter
    plot(f, freqSpectrum); hold off;
    title(sprintf("hft for f=%2.1fHz & col:%d", fsource(rowIdx), colIdx))
    for rowIdx = 1:N
        freqSpectrum = freqVal(:, end, colIdx);
        
        minIdx = abs(norm_freq-(fsource(rowIdx)/fsamp));
        meas_fmag(rowIdx, colIdx) = freqSpectrum(minIdx == min(minIdx, [], 2));
        fval(rowIdx, colIdx) = fsamp * norm_freq(minIdx == min(minIdx, [], 2));
    end
end
Measured_Frequency_Components = fUnits(meas_fmag, 'A')