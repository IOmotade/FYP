%% Simulates for the effect of increasing b i.e. number of bits per memristor
rst
addpath(genpath('../FYP/Sim_Scripts/'))
addpath(genpath('../FYP/Functions/'))

%% Instatiate baseline case
if ~exist('algo_case', 'var')
    algo_case = 1;
end
base_case;

loadPrevious = 2;
minNumBits = 1000;

%% Simulation Variables
numBitspRes = 2.^(0:4);
time_stamp = fTimeStamp;
idx_numBitspRes_Length = length(numBitspRes);
BER = zeros(length(algo), idx_numBitspRes_Length);

%% Run simulation for Algorithm X
fDisplayInternalMessage('Starting sim_numBitspRes Simulation');
tmp_prog_txtlen = 0;
for idx_numBitspRes=1:idx_numBitspRes_Length
    setup.numBitspRes = numBitspRes(idx_numBitspRes);
    BER(:, idx_numBitspRes) = fSimulation(algo, setup, minNumBits, loadPrevious);
    
    tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
    tmp_prog_txtlen = fDisplayInternalMessage(...
        sprintf('sim_numBitspRes: Simulation Progress: %2.2f percent', 100*(idx_numBitspRes/idx_numBitspRes_Length)),...
        tmp_prog_txtlen);
end
fClearInternalMessages(tmp_prog_txtlen);
fDisplayInternalMessage('sim_numBitspRes Simulation Complete');

%% Save Data
try
    foldername = 'Sim_Scripts/Results/';
    filename = sprintf('BER_numBitspRes_%d_%d_%s',...
        min(numBitspRes), max(numBitspRes), time_stamp);
    save([foldername filename]);
catch
    filename = sprintf('BER_numBitspRes_%d_%d_%s',...
        min(numBitspRes), max(numBitspRes), time_stamp);
    save(filename);
end

%% Plot Image and Save
plot(numBitspRes, BER', 'x-');
title('Plot of BER against Number of Bits per Memristor')
xlabel('Number of Bits per Memristor'); ylabel('BER')
ylim([0 0.5]); legend(strcat('Algorithm ', num2str(algo(:))))
grid on
grid minor
foldername = 'Figures/';
saveas(gcf, [foldername filename], 'png')