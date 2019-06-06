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
LRMean = 10.^(-1:8);
time_stamp = fTimeStamp;
idx_LRMean_Length = length(LRMean);
BER = zeros(1, idx_LRMean_Length);

%% Run simulation for Algorithm X
fDisplayInternalMessage('Starting sim_LRMean Simulation');
tmp_prog_txtlen = 0;
for idx_LRMean=1:idx_LRMean_Length
    setup.LRdef(1) = LRMean(idx_LRMean);
    BER(idx_LRMean) = fSimulation(algo, setup, minNumBits);
    
    tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
    tmp_prog_txtlen = fDisplayInternalMessage(...
        sprintf('sim_LRMean: Simulation Progress: %2.2f percent', 100*(idx_LRMean/idx_LRMean_Length)),...
        tmp_prog_txtlen);
end
fClearInternalMessages(tmp_prog_txtlen);
fDisplayInternalMessage('sim_LRMean Simulation Complete');

%% Save Data
try
    foldername = 'Sim_Scripts/Results/';
    filename = sprintf('%sBER_LRMean_%d_%d_%s', min(LRMean), max(LRMean), time_stamp);
    save(filename);
catch
    filename = sprintf('BER_LRMean_%d_%d_%s', min(LRMean), max(LRMean), time_stamp);
    save(filename);
end

semilogx(LRMean, BER);
title('Plot of BER against Mean of Line Resistance')
xlabel('Mean of Line Resistance'); ylabel('BER')
saveas(gcf, filename, 'png')