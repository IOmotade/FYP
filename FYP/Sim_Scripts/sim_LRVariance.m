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
LRVariance = 0:0.2:2;
time_stamp = fTimeStamp;
idx_LRVariance_Length = length(LRVariance);
BER = zeros(1, idx_LRVariance_Length);

%% Run simulation for Algorithm X
fDisplayInternalMessage('Starting sim_LRVariance Simulation');
tmp_prog_txtlen = 0;
for idx_LRVariance=1:idx_LRVariance_Length
    setup.LRdef(2) = sqrt(LRVariance(idx_LRVariance));
    BER(idx_LRVariance) = fSimulation(algo, setup, minNumBits);
    
    tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
    tmp_prog_txtlen = fDisplayInternalMessage(...
        sprintf('sim_LRVariance: Simulation Progress: %2.2f percent', 100*(idx_LRVariance/idx_LRVariance_Length)),...
        tmp_prog_txtlen);
end
fClearInternalMessages(tmp_prog_txtlen);
fDisplayInternalMessage('sim_LRVariance Simulation Complete');

%% Save Data
try
    foldername = 'Sim_Scripts/Results/';
    filename = sprintf('%sBER_LRVariance_%d_%d_%s', ...
        foldername, min(LRVariance), max(LRVariance), time_stamp);
    save(filename);
catch
    filename = sprintf('BER_LRVariance_%d_%d_%s', min(LRVariance), max(LRVariance), time_stamp);
    save(filename);
end

semilogx(LRVariance, BER);
title('Plot of BER against Variance of Line Resistance')
xlabel('Variance of Line Resistance'); ylabel('BER')
saveas(gcf, filename, 'png')