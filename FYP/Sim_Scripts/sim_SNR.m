%% Simulates for the effect of increasing N i.e. size of array.
rst
addpath(genpath('../FYP/Sim_Scripts/'))
addpath(genpath('../FYP/Functions/'))

%% Instatiate baseline case
if ~exist('algo_case', 'var')
    algo_case = 1;
end
base_case;

%% Simulation Variables
SNR_dB = -5:10;
time_stamp = fTimeStamp;
idx_SNR_dB_Length = length(SNR_dB);
BER = zeros(1, idx_SNR_dB_Length);

%% Run simulation for Algorithm X
fDisplayInternalMessage('Starting sim_SNR_dB Simulation');
tmp_prog_txtlen = 0;
for idx_SNR_dB=1:idx_SNR_dB_Length
    setup.SNR = SNR_dB(idx_SNR_dB);
    BER(idx_SNR_dB) = fSimulation(algo, setup, minNumBits);
    
    tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
    tmp_prog_txtlen = fDisplayInternalMessage(...
        sprintf('sim_SNR_dB: Simulation Progress: %2.2f percent', 100*(idx_SNR_dB/idx_SNR_dB_Length)),...
        tmp_prog_txtlen);
end
fClearInternalMessages(tmp_prog_txtlen);
fDisplayInternalMessage('sim_SNR_dB Simulation Complete');

%% Save Data
try
    foldername = 'Sim_Scripts/Results/';
    filename = sprintf('%sBER_SNR_dB_%d_%d_%s',...
        foldername, min(SNR_dB), max(SNR_dB), time_stamp);
    save(filename);
catch
    filename = sprintf('BER_SNR_dB_%d_%d_%s', min(SNR_dB), max(SNR_dB), time_stamp);
    save(filename);
end

plot(SNR_dB, BER);
title('Plot of BER against SNR')
xlabel('SNR/dB'); ylabel('BER')
ylim([0 0.5])
saveas(gcf, filename, 'png')