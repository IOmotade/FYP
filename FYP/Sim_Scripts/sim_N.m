%% Simulates for the effect of increasing N i.e. size of array.
rst
addpath(genpath('../FYP/Sim_Scripts/'))
addpath(genpath('../FYP/Functions/'))

%% Instatiate baseline case and other set up variables
if ~exist('algo_case', 'var')
    algo_case = 1;
end
base_case;

loadPrevious = 0;
minNumBits = 1000;

%% Simulation Variables
N = 2.^(0:5);
time_stamp = fTimeStamp;
idx_NLength = length(N);
BER = zeros(length(algo), idx_NLength);

%% Run simulation for Algorithm X
fDisplayInternalMessage('Starting sim_N Simulation');
tmp_prog_txtlen = 0;
for idx_N=1:idx_NLength
    setup.N = N(idx_N);
    BER(:, idx_N) = fSimulation(algo, setup, minNumBits, loadPrevious);
    
    tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
    tmp_prog_txtlen = fDisplayInternalMessage(...
        sprintf('sim_N: Simulation Progress: %2.2f percent', 100*(idx_N/idx_NLength)),...
        tmp_prog_txtlen);
end
fClearInternalMessages(tmp_prog_txtlen);
fDisplayInternalMessage('sim_N Simulation Complete');

%% Save Data
try
    foldername = 'Sim_Scripts/Results/';
    filename = sprintf('BER_N_%d_%d_%s',...
        min(N), max(N), time_stamp);
    save([foldername filename]);
catch
    filename = sprintf('BER_N_%d_%d_%s', min(N), max(N), time_stamp);
    save(filename);
end

%% Plot Image and Save
plot(N, BER', 'x-');
grid on
grid minor
title('Plot of BER against N')
xlabel('N'); ylabel('BER'); legend(strcat('Algorithm ', num2str(algo(:))))
ylim([-.09 0.6])
foldername = 'Figures/';
saveas(gcf, [foldername filename], 'png')