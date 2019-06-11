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
numBitspRes = 1:10;
time_stamp = fTimeStamp;
idx_numBitspRes_Length = length(numBitspRes);
BER{1} = zeros(1, idx_numBitspRes_Length);
BER{2} = zeros(1, idx_numBitspRes_Length);

%% Run simulation for Algorithm X

%Linear Rule
setup.ruleNum = 1;

fDisplayInternalMessage('Starting sim_rule_lin_log Simulation');
tmp_prog_txtlen = 0;
for idx_numBitspRes=1:idx_numBitspRes_Length
    setup.numBitspRes = numBitspRes(idx_numBitspRes);
    BER{1}(idx_numBitspRes) = fSimulation(algo, setup, minNumBits);
    
    tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
    tmp_prog_txtlen = fDisplayInternalMessage(...
        sprintf('sim_rule_lin_log1: Simulation Progress: %2.2f percent', 100*(idx_numBitspRes/idx_numBitspRes_Length)),...
        tmp_prog_txtlen);
end
fClearInternalMessages(tmp_prog_txtlen);

% Logarithmic Rule
setup.ruleNum = 2;

tmp_prog_txtlen = 0;
for idx_numBitspRes=1:idx_numBitspRes_Length
    setup.numBitspRes = numBitspRes(idx_numBitspRes);
    BER{2}(idx_numBitspRes) = fSimulation(algo, setup, minNumBits);
    
    tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
    tmp_prog_txtlen = fDisplayInternalMessage(...
        sprintf('sim_rule_lin_log2: Simulation Progress: %2.2f percent', 100*(idx_numBitspRes/idx_numBitspRes_Length)),...
        tmp_prog_txtlen);
end
fClearInternalMessages(tmp_prog_txtlen);

fDisplayInternalMessage('sim_rule_lin_log Simulation Complete');

%% Save Data
try
    foldername = 'Sim_Scripts/Results/';
    filename = sprintf('%sBER_rule_lin_log_numBitspRes_%d_%d_%s'...
        , foldername, min(numBitspRes), max(numBitspRes), time_stamp);
    save(filename);
catch
    filename = sprintf('BER_rule_lin_log_numBitspRes_%d_%d_%s', min(numBitspRes), max(numBitspRes), time_stamp);
    save(filename);
end

%% Plot Image and Save
plot(numBitspRes, BER{1}, numBitspRes, BER{2});
title('Plot of BER against Number of Bits per Memristor')
xlabel('Number of Bits per Memristor'); ylabel('BER');
legend('Linear Rule', 'Logarithmic Rule');
ylim([0 0.5])
saveas(gcf, filename, 'png')