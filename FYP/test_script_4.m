%% Simple Circuit
% Testing simple Line Resistance Estimation
addpath(genpath('../FYP/Functions/'))
rst

%% Setup Circuit System
N = 4;
% MemR = 10e3*ones(N);
LR = 10e3;
LRowR = abs(LR*ones(N)  + 0e3*randn(N));
LColR = abs(LR*ones(N)  + 0e3*randn(N));

%% Write to Memory
numBitspRes = 4;
[storedMemR, storedBits] = fWriteToMemArray(N, numBitspRes);

%% Carry Circuit Simulations
% Setup time samples
base_freq = 10;
fsource = base_freq*(2.^((1:N)-1))';

oversampfactor = 8;
fsamp = oversampfactor*2*max(fsource); tsamp = 1/fsamp;
nsamp = 1*(fsamp/base_freq);
t = 0:tsamp:(nsamp-1)*tsamp;

%% Sim Setup
MemR = storedMemR;

vs_mag = 5;
vs = vs_mag*square(2*pi*fsource*t);
Circuit = fMacSpiceSim(N, vs(:, 1), MemR, LRowR, LColR);
Circuit = repmat(Circuit, [nsamp, 1]);
val = zeros(size(t));

%% Run Simulation
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

%% Save Result
% save('SimData/Write_Read_Test_4_N_3_bitspM_4_LColR_10k_wparaRC.mat')
save('SimData/tmp_var')

%% Load Result
% load('SimData/Write_Read_Test_4_N_3_bitspM_4_LColR_10k_wparaRC.mat')
load('SimData/tmp_var')

% Get Haar Transform
timeCircuit = Circuit;

Circuit = fFrequencyDomain(timeCircuit);
timeVal = fGetFieldValues(Circuit.TimeDom, 'VS', 'value');
freqVal = fGetFieldValues(Circuit.FreqDom, 'VS', 'value');
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

%% With perfect values
if exist('corr_fmag', 'var')
    fmag = corr_fmag;
end

%% Memristor Value Finding Algorithm
[~, CalMemR] = fReadFromMemArray(N, numBitspRes, diag(vs_mag*ones(N, 1))/fmag);
prevCalMemRinner = zeros(N);
prevCalMemRouter = prevCalMemRinner;

% %% Calculate Memristor Values
[iter, Imem, eqResF] = deal(0, zeros(N), zeros(N));
[Vbias, Ibias] = deal(vs_mag*ones(N), zeros(N, 1));
[iterinner, iterouter] = deal(0, 0);
LColREst = mean(LColR(:))*ones(size(LColR));
exit_cond = false;

while ~exit_cond
    prevCalMemRouter = CalMemR;
%     [Vbias, Ibias] = deal(vs_mag*ones(N), zeros(N, 1));
    iterinner = 0;
    
    % First, estimate Memristor values using the source voltage magnitude
    while ~exit_cond
        for i = 1:N
            for j = 1:N
                if sum(isnan(CalMemR(:)))~=0
                   CalMemR(isnan(CalMemR)) = abs(prevCalMemRinner(isnan(CalMemR)));
                end
                
                prevCalMemRinner = CalMemR;
                
                nCF = fNodeCurrentFactors([], N, i, j, 0, abs(CalMemR), LRowR, LColR)
                eqResF(i, j) = nCF(i,1);
                [i, j] 
                totalNCF = prod(nCF(i:end, 1))
                if totalNCF == 0 ||  (sum(nCF(:)<0)~=0)
                    stop = true;
                end
                Imem(i, j) = fmag(i, j)/totalNCF
%                 Vbias
                CalMemR(i, j) = Vbias(i, j)/Imem(i, j) - nCF(i,1)*fEquivalentResistance("Down", N, i, j, abs(CalMemR), LRowR, LColREst)
            end
        end
        iterinner = iterinner + 1;
        exit_cond = ~(~fHasConverged(CalMemR, prevCalMemRinner) && (iterinner<=N^2));
    end
    [~, CalMemR] = fReadFromMemArray(N, numBitspRes, CalMemR);
    iterouter = iterouter + 1;
    
    Ibias(i) = sum(Imem(i, :), 2);
    for i = 1:N
        for j = 2:N
            tmp = Vbias(i, j-1) - abs((Ibias(i)-sum(Imem(i, 1:j-1), 2))*LRowR(i, j));
            
            % Check if bias voltage is valid value
            is_valid = (tmp>0)&&(tmp<vs_mag);
            Vbias(i, j) = is_valid*(tmp) + (~is_valid)*Vbias(i, j);
            CalMemR(i, j) = Vbias(i, j)/Imem(i, j) - eqResF(i, j)*fEquivalentResistance("Down", N, i, j, abs(CalMemR), LRowR, LColREst);
        end
    end
    Vbias
    exit_cond = ~(~fHasConverged(CalMemR, prevCalMemRouter) && (iterouter<=N^2)) && exit_cond;
end

Possible_Memristor_Values = fUnits(unique(sort(MemR))', 'Ohm');
% % Vbias
% fUnits([MemR CalMemR], 'Ohm')

% %% Work out Voltage adjustments
% i = 1;
% Vbias = vs_mag*ones(N);
% Ibias = sum(Imem, 2);
% for j = 2:N
%     Vbias(i, j) = Vbias(i, j-1) - (Ibias(i)-sum(Imem(i, 1:j-1), 2))*LRowR(i, j)
% end

% %% Read from Memory
[readBits, readMemR] = fReadFromMemArray(N, numBitspRes, CalMemR);
BER = sum(storedBits~=readBits)/length(readBits)
fprintf("Number Of Bits In Error is %d out of %d bits.\n", (BER * N^2 * numBitspRes), N^2 * numBitspRes)