function [BER, LREstimate] = fAlgorithm1(setup, N, numBitspRes, LRdef, ruleNum, varargin)
addpath(genpath('../FYP/Functions/'))

script = false;
%% Test as script
if script
    setup = [];
    [N, numBitspRes] = deal(4, 8);
    LRdef = [100, 0.1];
    ruleNum = 1;
    varargin = {[10], [10e6]};
end

%% Set-up additional variables
if ~isempty(varargin)
    var = squeeze(varargin);
else
    var = [];
end

%% Function Setup
if length(LRdef)==2
    [LRmean, LRsigma] = deal(LRdef(1), LRdef(2)*LRdef(1));
else
    [LRmean, LRsigma] = deal(LRdef(1), 0);
end

if ~exist('setup', 'var')
    default_setup = true;
else
    if isempty(setup)
        default_setup = true;
    else
        default_setup = false;
    end
end

if default_setup
    setup.perfect = false;
    setup.base_freq = 10;
    setup.oversampfactor = 8;
end

%% Setup Circuit System
LRowR = abs(LRmean*ones(N)  + LRsigma*randn(N));
LColR = abs(LRmean*ones(N)  + LRsigma*randn(N));

%% Write to Memory
[storedMemR, storedBits] = fWriteToMemArray(N, numBitspRes, ruleNum, [], var);

%% Carry Out Circuit Simulations
% Setup time samples
base_freq = setup.base_freq;
fsource = base_freq*(2.^((1:N)-1))';

oversampfactor = setup.oversampfactor;
fsamp = oversampfactor*2*max(fsource); tsamp = 1/fsamp;
nsamp = 1*(fsamp/base_freq);
t = 0:tsamp:(nsamp-1)*tsamp;

%% Sim Setup
MemR = storedMemR;

vs_mag = 5;
% vs = vs_mag*square(2*pi*fsource*t);
vs = fVoltageSourceSignals(N, vs_mag, nsamp);
Circuit = fMacSpiceSim(N, vs(:, 1), MemR, LRowR, LColR);
Circuit = repmat(Circuit, [nsamp, 1]);

%% Run Simulation
if ~setup.perfect
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
end

%% Save Result
% save('SimData/tmp_var')
filename = sprintf('SimData/Circuit_N_%d_oversampfactor_%d_numBitspRes_%d_ruleNum_%d_perfect_%d'...
    , N, oversampfactor, numBitspRes, ruleNum, setup.perfect);
try
    save(filename)
catch
    warning('Could not find SimData folder... Saving as tmp file');
    filename = 'tmp';
    save(filename)
end
%% Load Result
% load('SimData/tmp_var')
load(filename)

%% Get Haar Transform
timeCircuit = Circuit;
Circuit = fFrequencyDomain(timeCircuit);
norm_freq = ((1:nsamp)-1)/(nsamp-1);
% f = fsamp * norm_freq;

%% Extract Results from Haar Transform
if ~setup.perfect
    %% Extract measurement values
    freqVal = abs(fGetFieldValues(Circuit.FreqDom, 'IO', 'value'));
    [fmag, fval] = deal(zeros(N));
    for rowIdx = 1:N
        for colIdx = 1:N
            freqSpectrum = freqVal(:, end, colIdx);
            minIdx = abs(norm_freq-(fsource(rowIdx)/fsamp));
            fmag(rowIdx, colIdx) = freqSpectrum(minIdx == min(minIdx, [], 2));
            fval(rowIdx, colIdx) = fsamp * norm_freq(minIdx == min(minIdx, [], 2));
        end
    end
else
    %% Generate Accurate Frequency Components
    tmpVs = [vs_mag; zeros(N-1, 1)];
    tmpMemR = storedMemR;
    [corr_fmag, corr_fval] = deal(zeros(N));
    corr_Vbias = zeros(N);
    J = fShiftingMatrix(N);
    for i = 1:N
        for j = 1:N
            tmpCircuit = fMacSpiceSim(N, (J^(i-1))*tmpVs, tmpMemR, LRowR, LColR);
            im = tmpCircuit.IO.value(end, :);
            corr_fmag(i, j) = im(j);
            corr_fval(i, j) = fsource(i);
        end
        corr_Vbias(i, :) = tmpCircuit.VI.value(i, :);
    end
    fmag = corr_fmag;
end
% Frequency_Components = fUnits(fmag, 'A')

%% Memristor Value Finding Algorithm
% Initial guess for memristor values
[~, CalMemR] = fReadFromMemArray(N, numBitspRes, diag(vs_mag*ones(N, 1))/fmag, ruleNum, var);

% Setup
prevCalMemRinner = zeros(N);
% prevCalMemRouter = prevCalMemRinner;
[Imem, eqResF] = deal(zeros(N), zeros(N));
[~, iterouter] = deal(0, 0);
LColREst = LRmean*ones(size(LColR));
[Vbias, Ibias] = deal(vs_mag*ones(N), zeros(N, 1));
exit_cond = false;

% Calculate Memristor Values
while ~exit_cond
    prevCalMemRouter = CalMemR;
    %[Vbias, Ibias] = deal(vs_mag*ones(N), zeros(N, 1));
    
    % First, estimate Memristor values using the source voltage magnitude
    iterinner = 0;
    while ~exit_cond
        for i = 1:N
            for j = 1:N
                if sum(isnan(CalMemR(:)))~=0
                    CalMemR(isnan(CalMemR)) = abs(prevCalMemRinner(isnan(CalMemR)));
                end
                prevCalMemRinner = CalMemR;
                
                nCF = fNodeCurrentFactors([], N, i, j, 0, abs(CalMemR), LRowR, LColREst);
                totalNCF = prod(nCF(i:end, 1));
                eqResF(i, j) = nCF(i,1);
                
                %if totalNCF == 0 ||  (sum(nCF(:)<0)~=0)
                %   stop = true;
                %end
                
                Imem(i, j) = fmag(i, j)/totalNCF;
                
                CalMemR(i, j) = Vbias(i, j)/Imem(i, j) - ...
                    nCF(i,1)*fEquivalentResistance("Down", N, i, j,...
                    abs(CalMemR), LRowR, LColREst);
            end
        end
        iterinner = iterinner + 1;
        exit_cond = ~(~fHasConverged(CalMemR, prevCalMemRinner) && (iterinner<=N^2));
    end
    [~, CalMemR] = fReadFromMemArray(N, numBitspRes, CalMemR, ruleNum, var);
    iterouter = iterouter + 1;
    
    % Then, re-estimate Memristor values using estimated memristor currents
    Ibias(i) = sum(Imem(i, :), 2);
    for i = 1:N
        for j = 2:N
            %Calculate node bias voltage value and check if valid
            tmp = Vbias(i, j-1) - ...
                abs((Ibias(i)-sum(Imem(i, 1:j-1), 2))*LRowR(i, j));
            is_valid = (tmp>0)&&(tmp<vs_mag);
            
            Vbias(i, j) = is_valid*(tmp) + (~is_valid)*Vbias(i, j);
            
            CalMemR(i, j) = Vbias(i, j)/Imem(i, j) - ...
                eqResF(i, j)*fEquivalentResistance("Down", N, i, j,...
                abs(CalMemR), LRowR, LColREst);
        end
    end
    
    % Repeat till circuit state converges or MAX_ITERATIONS reached
    exit_cond = ~(~fHasConverged(CalMemR, prevCalMemRouter) && (iterouter<=N^2)) && exit_cond;
end

%% Estimate Line Resistances
LREstimate = zeros(N);
for i = 1:N
    for j = 1:N
        nCF = fNodeCurrentFactors([], N, i, j, 0, abs(CalMemR), LRowR, LColR);
        sumnCF = sum(nCF(i:end, 1));
        
        LREstimate(i, j) = (Vbias(i, j)/Imem(i, j) - MemR(i, j))/sumnCF;
        %         LREstimate(i, j) = Restimate(i, j)/eqResF(i, j);
    end
end

% Get rid of invalid values before averaging out
LREstimate = LREstimate(:);
LREstimate(LREstimate<=0) = [];
LREstimate = mean(LREstimate(:));

%% Read from Memory
[readBits, ~] = fReadFromMemArray(N, numBitspRes, CalMemR, ruleNum, var);
BER = sum(storedBits~=readBits)/length(readBits);
% fprintf("Number Of Bits In Error is %d out of %d bits.\n", (BER * N^2 * numBitspRes), N^2 * numBitspRes)

end