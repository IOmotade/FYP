%% Circuit System vs Time
rst
addpath(genpath('../FYP/Functions/'))
[fig_az, fig_el] = deal(145, 30);

%% Setup Circuit System
N = 3;
MemR = abs(10e3*ones(N) + 0e3*randn(N));
LRowR = abs(1e3*ones(N) + 0e3*randn(N));
LColR = abs(1e3*ones(N) + 0e3*randn(N));

vs_mag = [5, zeros(1, N-1)]';%Source Voltage Magnitude
vs_idx = 3;
J = fShiftingMatrix(N);
vs_mag = (J^(vs_idx-1))*vs_mag;
vs_mag(:) = 5;

is_mag = zeros(size(vs_mag(:, 1)));%Source Current Magnitude
is_mag(end) = 1;%0.0014;
is_mag(:) = 1;

C = ones(2*N, 1);%Connection Matrix
C = [0*C(1:N);
    zeros(N, 1)];

vs_ = [0 0 0 0 0 0].';
is_ = [0 0 0 0 0 1].';
C_  = [0 0 0 0 0 0].';

%% Setup time samples
fsamp = 10e3; tsamp = 1/fsamp;
base_freq = 100;
nsamp = 1;%0*(fsamp/base_freq);
t = 0:tsamp:(nsamp-1)*tsamp;

fsource = base_freq*(2.^((1:N)-1))';

%% Sim Setup
% vs =
vs = vs_mag.*square(2*pi*fsource*t);
Circuit = fMacSpiceSim(N, vs(:, 1), is_mag, C, MemR, LRowR, LColR);
Circuit = repmat(Circuit, [nsamp, 1]);
freqVal = zeros(size(t));

%% Run Sims
disp("Started")
tic
wBar = waitbar(0, 'Starting Simulation');
complete = zeros(nsamp, 1);
for idx=1:nsamp
%     Circuit(idx) = fMacSpiceSim(N, vs_, is_, C_, MemR, LRowR, LColR);
%     Circuit(idx) = fMacSpiceSim(N, vs_mag, MemR, LRowR, LColR);
    complete(idx) = 1;
    progress = sum(complete)/nsamp;
    waitbar(progress, wBar, sprintf('Simulation Progress: %2.2f percent', 100*progress))
end
waitbar(1, wBar, 'Simulation Complete');
toc
delete(wBar)
disp("Finished")
% play_alarm()

timeCircuit = Circuit;
Circuit = fFrequencyDomain(timeCircuit);
fUnits(timeCircuit.IO.value, 'A')
return
%% Extract Data & Plot Spectrum
timeVal = fGetFieldValues(Circuit.TimeDom, 'X');

figure;
lgnd = repmat(" ", N, 1);
lgdIdx = 1;
for uIdx = N:N
    plot(1:N, timeVal(uIdx, :));
    xlim([1 N]);
    lgnd(lgdIdx) = sprintf("Plot of U(%s,%s)", 'x', num2str(uIdx));
    hold on
    lgdIdx = lgdIdx + 1;
end
hold off
legend(lgnd)
title("Plot of U");
figure;
surf(timeVal(1:N, 1:N)); view([fig_az fig_el]);
set(gca,'xdir','reverse')
xlabel('col'); ylabel('row'); zlabel('mag');
title("Plot of U");

%% V end
timeVal = fGetFieldValues(Circuit.TimeDom, 'VO');

figure;
lgnd = repmat(" ", N, 1);
lgdIdx = 1;
for vIdx = N:N
    plot(1:N, timeVal(:, vIdx));
    xlim([1 N]);
    lgnd(lgdIdx) = sprintf("Plot of V(%s,%s)", num2str(vIdx), 'y');
    hold on
    lgdIdx = lgdIdx + 1;
end
hold off
legend(lgnd)
title("Plot of V");
figure;
surf(timeVal); view([fig_az fig_el]);
set(gca,'xdir','reverse')
xlabel('col'); ylabel('row'); zlabel('mag');
title("Plot of V");

%% IO end
timeVal = fGetFieldValues(Circuit.TimeDom, 'IO');

figure;
lgnd = repmat(" ", N, 1);
lgdIdx = 1;
for uIdx = N:N
    plot(1:N, timeVal(uIdx, :));
    xlim([1 N]);
    lgnd(lgdIdx) = sprintf("Plot of IO(%s,%s)", 'x', num2str(uIdx));
    hold on
    lgdIdx = lgdIdx + 1;
end
hold off
legend(lgnd)
title("Plot of IO");
figure;
surf(timeVal); view([fig_az fig_el]);
set(gca,'xdir','reverse')
xlabel('col'); ylabel('row'); zlabel('mag');
title("Plot of IO");

%% Plot of Difference Relations
figure;
val = diff(timeVal, 1, 2);
surf(val); view([fig_az fig_el]);
set(gca,'xdir','reverse')
xlabel('col'); ylabel('row'); zlabel('mag');
title("Plot of Difference in X-Dir(Col-Dir)");

%%
function play_alarm
load gong
sound(y(1:10000),Fs)
end