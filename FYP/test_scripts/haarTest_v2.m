%% Circuit System vs Time
rst
addpath(genpath('../FYP/Functions/'))

%% Setup Circuit System
N = 3;
MemR = 10e3*ones(N); LRowR = 1e3*ones(N); LColR = 1e3*ones(N);
% MemR(1, 1) = 1e3; MemR(1, 2) = 3e3; MemR(1, 3) = 5e3;
% MemR(2, 1) = 1e3; MemR(2, 2) = 10e3; MemR(2, 3) = 1e3;
% MemR(3, 1) = 1e3; MemR(3, 2) = 1e2; MemR(3, 3) = 60e3;
vs_mag = 5;%Source Voltage Magnitude
% Circuit = fMacSpiceSim(N, vs, MemR, LRowR, LColR);
% fUnits(Circuit.VO.value, 'V')
% fUnits(Circuit.IO.value, 'A')

%% Setup time samples
base_freq = 20;

fsource = base_freq*(2.^((1:N)-1))';

fsamp = max(fsource*2); tsamp = 1/fsamp;
nsamp = (fsamp/base_freq);
t = 0:tsamp:(nsamp-1)*tsamp;

%% Sim Setup
% vs =
% vs = vs_mag*square(2*pi*fsource*t);
vs = fVoltageSourceSignals(N, vs_mag, nsamp);
Circuit = fMacSpiceSim(N, vs(:, 1), MemR, LRowR, LColR);
Circuit = repmat(Circuit, [nsamp, 1]);
% freqVal = zeros(size(t));

%% Run Sims
disp("Started")
tic
wBar = waitbar(0, 'Starting Simulation');
complete = zeros(nsamp, 1);
for idx=1:nsamp
    Circuit(idx) = fMacSpiceSim(N, vs(:, idx), MemR, LRowR, LColR);
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
% timeCircuit = Circuit;

%% Extract Data & Plot Spectrum
timeVal = fGetFieldValues(timeCircuit, 'IO', 'value');
fHaarT(timeVal(:, 3, 3))

FullCircuit = fFrequencyDomain(timeCircuit);
timeVal = fGetFieldValues(FullCircuit.TimeDom, 'VS', 'value');
freqVal = fGetFieldValues(FullCircuit.FreqDom, 'VS', 'value');
norm_freq = ((1:nsamp)-1)/(nsamp-1);

idx = 3;

figure;
plot(norm_freq*fsamp, freqVal(:, idx, 1))
title("Haar Transform of Something: Version 1")

figure;
plot(t, timeVal(:, idx, 1))
title("Time Signal of Something")

% figure;
% plot(t, fInvHaarT(transpose(freqVal(:, idx, 1))))
% title("Inverse Haar Transform of Something: Version 1")

%%
timeVal = zeros([1 length(Circuit)]);
for rowIdx = idx:idx
    for colIdx = 1:1
        for tIdx=1:nsamp
            timeVal(tIdx) = Circuit(tIdx).VS.value(rowIdx, colIdx);
        end
        tmpFreqVal = fHaarT(timeVal);
%         for fIdx=1:nSamp
%             freqCircuit(fIdx).VS.value(rowIdx, colIdx) = freqVal(fIdx);
%         end
    end
end
figure;
plot(norm_freq*fsamp, tmpFreqVal)
title("Haar Transform of Something: Version 2")

figure;
plot(t, fInvHaarT(tmpFreqVal))
title("Inverse Haar Transform of Something: Version 2")

%%
function play_alarm
load handel
sound(y,Fs)
end