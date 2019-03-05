%% Circuit System vs Time
rst
addpath(genpath('../FYP/Functions/'))

%% Setup Circuit System
N = 3;
MemR = 10e3*ones(N); LRowR = 1e3*ones(N); LColR = 1e3*ones(N);
MemR(1, 1) = 1e3; MemR(1, 2) = 3e3; MemR(1, 3) = 5e3;
MemR(2, 1) = 1e3; MemR(2, 2) = 10e3; MemR(2, 3) = 1e3;
MemR(3, 1) = 1e3; MemR(3, 2) = 1e2; MemR(3, 3) = 60e3;
vs_mag = 5;%Source Voltage Magnitude
% Circuit = fMacSpiceSim(N, vs, MemR, LRowR, LColR);
% fUnits(Circuit.VO.value, 'V')
% fUnits(Circuit.IO.value, 'A')

%% Setup time samples
fsamp = 10e3; tsamp = 1/fsamp;
base_freq = 100;
nsamp = 10*(fsamp/base_freq);
t = 0:tsamp:(nsamp-1)*tsamp;

fsource = base_freq*(2.^((1:N)-1))';
%% Sim Setup
% vs =
vs = vs_mag*square(2*pi*fsource*t);
Circuit = fMacSpiceSim(N, vs(:, 1), MemR, LRowR, LColR);
Circuit = repmat(Circuit, [nsamp, 1]);
freqVal = zeros(size(t));

%% Run Sims
disp("Started")
tic
wBar = waitbar(0, 'Starting Simulation');
complete = zeros(nsamp, 1);
for idx=1:nsamp
    Circuit(idx) = fMacSpiceSim(N, vs(:, idx), MemR, LRowR, LColR, '2');
    complete(idx) = 1;
    progress = sum(complete)/nsamp;
    waitbar(progress, wBar, sprintf('Simulation Progress: %2.2f percent', 100*progress))
end
waitbar(1, wBar, 'Simulation Complete');
toc
delete(wBar)
disp("Finished")
play_alarm()
timeCircuit = Circuit;

%% Save Result
% save('Circuitfsamp10kHz_para_non_uniform_1_scheme_2')

%% Load Result
% load('SimData/Circuitfsamp10kHz_nopara.mat')
% load('Circuitfsamp100kHz_para_non_uniform_1')
timeCircuit = Circuit;

%% Extract Data & Plot Spectrum
Circuit = fFrequencyDomain(timeCircuit);
timeVal = fGetFieldValues(Circuit.TimeDom, 'VO');
freqVal = fGetFieldValues(Circuit.FreqDom, 'VO');
norm_freq = ((1:nsamp)-1)/(nsamp-1);

figure;
plot(norm_freq, freqVal(:, end, 1))
title("Haar Transform of Something")

figure;
plot(t, timeVal(:, end, 1))
title("Time Signal of Something")

%%
function play_alarm
load handel
sound(y,Fs)
end