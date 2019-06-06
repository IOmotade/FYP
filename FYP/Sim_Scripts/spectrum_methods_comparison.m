rst
internal_msg_len = 0;

%% Circuit Setup
N = 4;
[MemR, LRowR, LColR] = deal(10e3*ones(N), 1e3*ones(N), 1e3*ones(N));

%% Carry Out Circuit Simulations
vs_mag = 5;

% Setup time samples
basefreq = 10;
fsource = basefreq*(2.^((1:N)-1))';

oversampfactor = 2.^(0:3);
fsamp = oversampfactor(1)*2*max(fsource);
nsamp = 1*(fsamp/basefreq);
norm_freq = ((1:nsamp)-1)/(nsamp-1);

internal_msg_len = ...
    fDisplayInternalMessage('Starting Simulation', internal_msg_len);

circuitIdx = 1;
circuitIdxLength = length(oversampfactor) + 1;

%% Haar-like method
vs = fVoltageSourceSignals(N, vs_mag, nsamp);
Circuit = fMacSpiceSim(N, vs(:, 1), MemR, LRowR, LColR);
Circuit = repmat(Circuit, [nsamp, 1]);

tmp_prog_txtlen = 0;
for idx=1:nsamp
    Circuit(idx) = fMacSpiceSim(N, vs(:, idx), MemR, LRowR, LColR);
    
    tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
    tmp_prog_txtlen = fDisplayInternalMessage(...
        sprintf('Simulation Progress: %2.2f percent', 100*(idx/nsamp)),...
        tmp_prog_txtlen);
end
fClearInternalMessages(tmp_prog_txtlen);
spectrum_method_time{circuitIdx} = Circuit;
f{circuitIdx} = fsamp * norm_freq;
circuitIdx = circuitIdx + 1;

%% Other methods
tmp_prog_txtlen_1 = 0;
for o_idx = 1:length(oversampfactor)
    %     fsource = basefreq*(1:N)';
    fsource = basefreq*(2.^((1:N)-1))';
    fsamp = oversampfactor(o_idx)*2*max(fsource);
    nsamp = 1*(fsamp/basefreq);
    norm_freq = ((1:nsamp)-1)/(nsamp-1);
    
    tsamp = 1/fsamp;
    t = 0:tsamp:(nsamp-1)*tsamp;
    vs = vs_mag*square(2*pi*fsource*t);
    
    Circuit = fMacSpiceSim(N, vs(:, 1), MemR, LRowR, LColR);
    Circuit = repmat(Circuit, [nsamp, 1]);
    
    %%
    tmp_prog_txtlen = 0;
    for idx=1:nsamp
        Circuit(idx) = fMacSpiceSim(N, vs(:, idx), MemR, LRowR, LColR);
        
        tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
        tmp_prog_txtlen = fDisplayInternalMessage(...
            sprintf('Simulation Progress: %2.2f percent', 100*(idx/nsamp)),...
            tmp_prog_txtlen);
    end
    fClearInternalMessages(tmp_prog_txtlen);
    spectrum_method_time{circuitIdx} = Circuit;
    f{circuitIdx} = fsamp * norm_freq;
    
    %%
    tmp_prog_txtlen_1 = fClearInternalMessages(tmp_prog_txtlen_1);
    tmp_prog_txtlen_1 = fDisplayInternalMessage(...
        sprintf('Simulation Progress: %2.2f percent', 100*(circuitIdx/circuitIdxLength)),...
        tmp_prog_txtlen_1);
    
    circuitIdx = circuitIdx + 1;
end
fClearInternalMessages(tmp_prog_txtlen_1);

internal_msg_len = ...
    fDisplayInternalMessage('Simulation Complete', internal_msg_len);

%% Get frequency spectrum
spectrum_method = spectrum_method_time;
spectrum_method{1} = fFrequencyDomain(spectrum_method_time{1});
for circuitIdx = 2:circuitIdxLength
    spectrum_method{circuitIdx} = fFrequencyDomain(spectrum_method_time{circuitIdx}, '1');
end

save('SimData/spec_mthd_comp')

%% Plot Spectrums
load('SimData/spec_mthd_comp');
[plotIdx, plotVar] = deal(1, 'VS');
figure;
spec = fGetFieldValues(spectrum_method{1}.FreqDom, plotVar, 'value');
subplot(2, 4, [1, 2, 5, 6]), fig = stem(f{1}, spec(:, plotIdx, 1));
fsamp = oversampfactor(1)*2*max(fsource);
title(sprintf('f_s_a_m_p=%d*f_o ', fsamp/basefreq));
% title(sprintf('Sampling Frequency: %d*f_o ', fsamp/basefreq))
hold on
fsource = basefreq*(2.^((1:N)-1))';
for i=1:1
    [~, j] = min(abs(f{1}-fsource(i)));
    plot(f{1}(j), spec(j, plotIdx, 1), 'rx')
end
hold off
set(fig,'MarkerFaceColor','blue','Marker','none')
xlim([0, 200]); ylim([0, 5.5]);

% figure;
subplotIdx = [3, 4, 7, 8];
for circuitIdx = 2:circuitIdxLength
    spec = fGetFieldValues(spectrum_method{circuitIdx}.FreqDom, plotVar, 'value');
    subplot(2, 4, subplotIdx(circuitIdx-1)), fig = stem(f{circuitIdx}, spec(:, plotIdx, 1));
    set(fig,'MarkerFaceColor','blue','Marker','none')
    xlim([0, 200]); ylim([0, 5.5]);
    %%
    hold on
    %     fsource = basefreq*(1:N)';
    fsource = basefreq*(2.^((1:N)-1))';
    for i=1:1
        [~, j] = min(abs(f{1}-fsource(i)));
        plot(f{1}(j), spec(j, plotIdx, 1), 'rx');
    end
    hold off
    
    fsamp = oversampfactor(circuitIdx-1)*2*max(fsource);
    title(sprintf('f_s_a_m_p=%d*f_o ', fsamp/basefreq));
end
set(fig,'MarkerFaceColor','blue','Marker','none')
sgtitle(sprintf('Frequency Spectrum for Base Frequency Bias Signal f_o = %dHz, N = %d', basefreq, N));
set(gcf, 'Position',  [100, 100, 1000, 400])