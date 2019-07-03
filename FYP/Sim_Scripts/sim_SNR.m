%% Simulates for the effect of increasing SNR of Bias Vector Signals
rst
addpath(genpath('../FYP/Sim_Scripts/'))
addpath(genpath('../FYP/Functions/'))

%% Instatiate baseline case
if ~exist('algo_case', 'var')
    algo_case = 1;
end
base_case;
loadPrevious = 1;
minNumBits = 1000;

%% Simulation Variables
SNR_dB = 0:5:40;
time_stamp = fTimeStamp;
idx_SNR_dB_Length = length(SNR_dB);
BER = zeros(length(algo), idx_SNR_dB_Length);

%% Run simulation for Algorithm X
fDisplayInternalMessage('Starting sim_SNR_dB Simulation');
tmp_prog_txtlen = 0;
for idx_SNR_dB=1:idx_SNR_dB_Length
    setup.SNR = SNR_dB(idx_SNR_dB);
    BER(:, idx_SNR_dB) = fSimulation(algo, setup, minNumBits, loadPrevious);
    
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
    filename = sprintf('BER_SNR_dB_%d_%d_%s',...
        min(SNR_dB), max(SNR_dB), time_stamp);
    save([foldername filename]);
catch
    filename = sprintf('BER_SNR_dB_%d_%d_%s', min(SNR_dB), max(SNR_dB), time_stamp);
    save(filename);
end

%%
plot(SNR_dB, BER', 'x-');
title('Plot of BER against SNR')
xlabel('SNR/dB'); ylabel('BER')
ylim([0 0.5]); legend(strcat('Algorithm ', num2str(algo(:))))
grid on
grid minor
foldername = 'Figures/';
saveas(gcf, [foldername filename], 'png')