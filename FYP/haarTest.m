%% Circuit System vs Time
rst
addpath(genpath('../FYP/Functions/'))

%% Setup Circuit System
N = 3;
MemR = 10e3*ones(N); LRowR = 1e0*ones(N); LColR = 1e0*ones(N);
MemR(1, 1) = 1e3; MemR(1, 2) = 3e3; MemR(1, 3) = 5e3;
MemR(2, 1) = 1e3; MemR(2, 2) = 10e3; MemR(2, 3) = 1e3;
MemR(3, 1) = 1e3; MemR(3, 2) = 1e2; MemR(3, 3) = 60e3;
vs_mag = 5;%Source Voltage Magnitude
% Circuit = fMacSpiceSim(N, vs, MemR, LRowR, LColR);
% fUnits(Circuit.VO.value, 'V')
% fUnits(Circuit.IO.value, 'A')

%% Setup time samples
fsamp = 100e3; tsamp = 1/fsamp;
base_freq = 100;
nsamp = 10*(fsamp/base_freq);
t = 0:tsamp:(nsamp-1)*tsamp;

fsource = base_freq*(2.^((1:N)-1))';
%% Sim Setup
% vs =
vs = vs_mag*square(2*pi*fsource*t);
Circuit = fMacSpiceSim(N, vs(:, 1), MemR, LRowR, LColR);
Circuit = repmat(Circuit, [nsamp, 1]);
val = zeros(size(t));

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

%% Save Result
% save('Circuitfsamp100kHz_para_non_uniform_1')

%% Load Result
% load('SimData/Circuitfsamp10kHz_nopara.mat')
load('Circuitfsamp100kHz_para_non_uniform_1')

%% Extract Data
% close all
for idx=1:nsamp
    %     val(idx) = Circuit(idx).VO.value(end, 1);
    %         val(idx) = Circuit(idx).VI.value(end, 1);
%     val(idx) = Circuit(idx).IO.value(end, 3);
    %     val(idx) = Circuit(idx).II.value(end, 1);
        val(idx) = Circuit(idx).VS.value(2, 1);
    %     val(idx) = Circuit(idx).IS.value(end, 1);
    %     val(idx) = Circuit(idx).X.value(end, 1);
end

hft = fHaarT(val);
norm_freq = (((1:length(hft))-1)/length(hft));%*fsamp;
norm_freq(abs(hft)==max(abs(hft)))*fsamp
% hft(abs(hft)<max(abs(hft))/5.5) = 0; %filter
hift = fInvHaarT(hft);

for i = 1:N
    minIdx = abs(norm_freq-(fsource(i)/fsamp));
    val_freqcomp(i) = hft(minIdx == min(minIdx, [], 2));
    f_read(i) = norm_freq(minIdx == min(minIdx, [], 2));
end
val_freqcomp
% f_read
f_read*fsamp

% %% Plots
% close all
figure;
plot(t, val);
xlabel("Time"); ylabel("Value")
title("measured signal")

figure;
plot(norm_freq, hft)
xlabel("Normalized Frequency (f_s_i_g/f_s_a_m_p)"); ylabel("Magnitude")
title("haar transform")

% figure
% plot(norm_freq, angle(hft)/pi)
% xlabel("Normalized Frequency (f_s_i_g/f_s_a_m_p)"); ylabel("pi")
% title("haar transform phase")

figure;
plot(t, hift)
xlabel("Time"); ylabel("Value")
title("recovered signal")

function play_alarm
load handel
sound(y,Fs)
end