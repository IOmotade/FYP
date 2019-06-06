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
MemRRange = 1;
MemRMean = 10.^(2:5);
time_stamp = fTimeStamp;
idx_MemRMean_Length = length(MemRMean);
BER = zeros(1, idx_MemRMean_Length);

%% Run simulation for Algorithm X
fDisplayInternalMessage('Starting sim_MemRMean Simulation');
tmp_prog_txtlen = 0;
for idx_MemRMean=1:idx_MemRMean_Length
    setup.var = {[MemRMean(idx_MemRMean)/(10^MemRRange)], [MemRMean(idx_MemRMean)*(10^MemRRange)]};
    BER(idx_MemRMean) = fSimulation(algo, setup, minNumBits);
    
    tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
    tmp_prog_txtlen = fDisplayInternalMessage(...
        sprintf('sim_LRMean: Simulation Progress: %2.2f percent', 100*(idx_MemRMean/idx_MemRMean_Length)),...
        tmp_prog_txtlen);
end
fClearInternalMessages(tmp_prog_txtlen);
fDisplayInternalMessage('sim_MemRMean Simulation Complete');

%% Save Data
try
    foldername = 'Sim_Scripts/Results/';
    filename = sprintf('%sBER_MemRMean_%d_%d_%s', min(MemRMean), max(MemRMean), time_stamp);
    save(filename);
catch
    filename = sprintf('BER_MemRMean_%d_%d_%s', min(MemRMean), max(MemRMean), time_stamp);
    save(filename);
end

semilogx(MemRMean, BER);
title('Plot of BER against Mean of Line Resistance')
xlabel('Log Mean of Memristor Resistance'); ylabel('BER')
saveas(gcf, filename, 'png')