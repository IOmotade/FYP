%% Simple Circuit
% Testing simple Line Resistance Estimation
addpath(genpath('../FYP/Functions/'))
rst

%% Setup Circuit System
N = 4;
% MemR = 10e3*ones(N);
LRowR = abs(1e-3*ones(N)  + 0e3*randn(N));
LColR = abs(1e3*ones(N)  + 0e0*randn(N));

%% Write to Memory
tic
numBitspRes = 4;
[storedMemR, storedBits] = fWriteToMemArray(N, numBitspRes);

MemR = storedMemR;

%% Carry Circuit Simulations
% Setup time samples
base_freq = 10;
fsource = base_freq*(2.^((1:N)-1))';

oversampfactor = 32;
fsamp = oversampfactor*2*max(fsource); tsamp = 1/fsamp;
nsamp = 1*(fsamp/base_freq);
t = 0:tsamp:(nsamp-1)*tsamp;

%% Sim Setup
vs_mag = 5;
vs = vs_mag*square(2*pi*fsource*t);
Circuit = fMacSpiceSim(N, vs(:, 1), MemR, LRowR, LColR);
Circuit = repmat(Circuit, [nsamp, 1]);
val = zeros(size(t));

%% Run Simulation
disp("Started")
% tic
wBar = waitbar(0, 'Starting Simulation');
complete = zeros(nsamp, 1);
for idx=1:nsamp
    Circuit(idx) = fMacSpiceSim(N, vs(:, idx), MemR, LRowR, LColR);
    complete(idx) = 1;
    progress = sum(complete)/nsamp;
    waitbar(progress, wBar, sprintf('Simulation Progress: %2.2f percent', 100*progress))
end
waitbar(1, wBar, 'Simulation Complete');
% toc
delete(wBar)
disp("Finished")

%% Save Result
% save('SimData/Write_Read_Test_4_N_3_LColR_10k_wparaRC.mat')
save('SimData/tmp_var')

%% Load Result
% load('SimData/Write_Read_Test_16_N_4_wpara.mat')
load('SimData/tmp_var')

% Get Haar Transform
timeCircuit = Circuit;

Circuit = fFrequencyDomain(timeCircuit);
% timeVal = fGetFieldValues(Circuit.TimeDom, 'VS', 'value');
% freqVal = fGetFieldValues(Circuit.FreqDom, 'VS', 'value');
norm_freq = ((1:nsamp)-1)/(nsamp-1);
f = fsamp * norm_freq;

% [rowIdx, colIdx] = deal(1, 1);
%
% figure;
% plot(f, freqVal(:, rowIdx, colIdx));
% title(sprintf("Haar Transform of Val(%d, %d)", rowIdx, colIdx));
%
% figure;
% plot(t, timeVal(:, rowIdx, colIdx));
% title(sprintf("Time Signal of Val(%d, %d)", rowIdx, colIdx));
% return
% Organise into frequency components
freqVal = fGetFieldValues(Circuit.FreqDom, 'IO', 'value');
[fmag, fval] = deal(zeros(N));
% filt_tol = 10;
for rowIdx = 1:N
    for colIdx = 1:N
        freqSpectrum = freqVal(:, end, colIdx);
        
        %         figure;
        %         plot(f, freqSpectrum, 'r'); hold on;
        %         freqSpectrum(abs(freqSpectrum)<max(abs(freqSpectrum))/filt_tol) = 0; %filter
        %         plot(f, freqSpectrum); hold off;
        %         title(sprintf("hft for f=%2.1fHz & col:%d", fsource(rowIdx), colIdx))
        
        minIdx = abs(norm_freq-(fsource(rowIdx)/fsamp));
        fmag(rowIdx, colIdx) = freqSpectrum(minIdx == min(minIdx, [], 2));
        fval(rowIdx, colIdx) = fsamp * norm_freq(minIdx == min(minIdx, [], 2));
    end
end
Frequency_Components = fUnits(fmag, 'A')

% %% Memristor Value Finding Algorithm
[~, CalMemR] = fReadFromMemArray(N, numBitspRes, diag(vs_mag*ones(N, 1))/fmag);
prevCalMemRinner = zeros(N);
prevCalMemRouter = prevCalMemRinner;

% %% Calculate Memristor Values
[iterinner, iterouter, Imem] = deal(0, 0, zeros(N));
LColREst = mean(LColR(:))*ones(size(LColR));
while (~fHasConverged(CalMemR, prevCalMemRouter) && (iterouter<=N))
    prevCalMemRouter = CalMemR;
    
    iterinner = 0;
    while (~fHasConverged(CalMemR, prevCalMemRinner) && (iterinner<=N))
        for i = 1:N
            for j = 1:N
                prevCalMemRinner = CalMemR;
                nCF = fNodeCurrentFactors([], N, i, j, 0, abs(CalMemR), LRowR, LColREst);
                totalNCF = prod(nCF(i:end, 1));
                Imem(i, j) = fmag(i, j)/totalNCF;
                CalMemR(i, j) = vs_mag/Imem(i, j) - nCF(i,1)*fEquivalentResistance("Down", N, i, j, abs(CalMemR), LRowR, LColREst);
            end
        end
        iterinner = iterinner + 1;
    end
    
    [~, CalMemR] = fReadFromMemArray(N, numBitspRes, abs(CalMemR));
    iterouter = iterouter + 1;
end
Possible_Memristor_Values = fUnits(unique(sort(MemR))', 'Ohm');
% fUnits([MemR CalMemR], 'Ohm')
toc
% %% Read from Memory
[readBits, readMemR] = fReadFromMemArray(N, numBitspRes, abs(CalMemR));
BER = sum(storedBits~=readBits)/length(readBits)
fprintf("Number Of Bits In Error is %d out of %d bits.\n", (BER * N^2 * numBitspRes), N^2 * numBitspRes)