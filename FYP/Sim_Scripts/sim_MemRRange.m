%% Simulates for the effect of increasing Range(M) i.e. Log Range of Possible Memristor Values
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
MemRMean = 10^3;
MemRRange = 1:3;
time_stamp = fTimeStamp;
idx_MemRRange_Length = length(MemRRange);
BER = zeros(length(algo), idx_MemRRange_Length);

%% Run simulation for Algorithm X
fDisplayInternalMessage('Starting sim_MemRMean Simulation');
tmp_prog_txtlen = 0;
for idx_MemRRange=1:idx_MemRRange_Length
    setup.var = {[MemRMean/(10^MemRRange(idx_MemRRange))], [MemRMean*(10^MemRRange(idx_MemRRange))]};
    BER(:, idx_MemRRange) = fSimulation(algo, setup, minNumBits, loadPrevious);
    
    tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
    tmp_prog_txtlen = fDisplayInternalMessage(...
        sprintf('sim_MemRRange: Simulation Progress: %2.2f percent', 100*(idx_MemRRange/idx_MemRRange_Length)),...
        tmp_prog_txtlen);
end
fClearInternalMessages(tmp_prog_txtlen);
fDisplayInternalMessage('sim_MemRRange Simulation Complete');

%% Save Data
try
    foldername = 'Sim_Scripts/Results/';
    filename = sprintf('BER_MemRRange_%1.0f_%1.0f_%s',...
        min(MemRRange), max(MemRRange), time_stamp);
    save([foldername filename]);
catch
    filename = sprintf('BER_MemRRange_%1.0f_%1.0f_%s', min(MemRRange), max(MemRRange), time_stamp);
    save(filename);
end

%%
semilogx(MemRRange, BER', 'x-');
title('Plot of BER against Log Range of Memristor Resistance')
xlabel('Log Range of Memristor Resistance'); ylabel('BER')
ylim([0 0.5]); legend(strcat('Algorithm ', num2str(algo(:))))
grid on
grid minor
foldername = 'Figures/';
saveas(gcf, [foldername filename], 'png')