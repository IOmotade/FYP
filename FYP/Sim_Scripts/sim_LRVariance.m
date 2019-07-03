%% Simulates for the effect of increasing Variance of Line Ressistance Values
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
LRVariance = 0:.01:.1;
time_stamp = fTimeStamp;
idx_LRVariance_Length = length(LRVariance);
BER = zeros(length(algo), idx_LRVariance_Length);

%% Run simulation for Algorithm X
fDisplayInternalMessage('Starting sim_LRVariance Simulation');
tmp_prog_txtlen = 0;
for idx_LRVariance=1:idx_LRVariance_Length
    setup.LRdef(2) = sqrt(LRVariance(idx_LRVariance));
    BER(:, idx_LRVariance) = fSimulation(algo, setup, minNumBits, loadPrevious);
    
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
    filename = sprintf('BER_LRVariance_%d_%d_%s', ...
        min(LRVariance), max(LRVariance), time_stamp);
    save([foldername filename]);
catch
    filename = sprintf('BER_LRVariance_%d_%d_%s', min(LRVariance)*100, max(LRVariance)*100, time_stamp);
    save(filename);
end

%%
plot(LRVariance, BER', 'x-');
title('Plot of BER against Variance of Line Resistance')
xlabel('Variance of Line Resistance'); ylabel('BER')
ylim([0 0.5]); legend(strcat('Algorithm ', num2str(algo(:))))
grid on
grid minor
foldername = 'Figures/';
saveas(gcf, [foldername filename], 'png')