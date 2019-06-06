function [BER, LREstimate] = fAlgorithm1(setup, N, numBitspRes, LRdef, ruleNum, varargin)
addpath(genpath('../FYP/Functions/'))
internal_msg_len = 0;
LR_maxdev = 0.1;
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
    [LRmean, LRsigma, LRmeansigma] = deal(LRdef(1), LRdef(2)*LRdef(1), 0);
elseif length(LRdef)==3
    [LRmean, LRsigma, LRmeansigma] = deal(LRdef(1), LRdef(2)*LRdef(1),...
        LRdef(3)*LRdef(1));
else
    [LRmean, LRsigma, LRmeansigma] = deal(LRdef(1), 0, 0);
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
    setup.basefreq = 10;
    setup.oversampfactor = 8;
    setup.SNR = Inf;
end

if ~fFieldExist(setup, 'SNR')
    setup.SNR = Inf;
end

%% Setup Circuit System
mean_offset = LRmeansigma*randn(1);
LRowR = abs(LRmean*ones(N)  + LRsigma*randn(N) + mean_offset);
LColR = abs(LRmean*ones(N)  + LRsigma*randn(N) + mean_offset);

%% Write to Memory
[storedMemR, storedBits] = fWriteToMemArray(N, numBitspRes, ruleNum, [], var);

%% Carry Out Circuit Simulations
% Setup time samples
basefreq = setup.basefreq;
fsource = basefreq*(2.^((1:N)-1))';

oversampfactor = setup.oversampfactor;
fsamp = oversampfactor*2*max(fsource);
nsamp = 1*(fsamp/basefreq);
% tsamp = 1/fsamp;
% t = 0:tsamp:(nsamp-1)*tsamp;

%% Sim Setup
MemR = storedMemR;

vs_mag = 5;
% vs = vs_mag*square(2*pi*fsource*t);
vs = fVoltageSourceSignals(N, vs_mag, nsamp);
vs = vs + (vs_mag/db2mag(setup.SNR))*randn(size(vs));
Circuit = fMacSpiceSim(N, vs(:, 1), MemR, LRowR, LColR);
Circuit = repmat(Circuit, [nsamp, 1]);

%% Run Simulation
if ~setup.perfect
    %     internal_msg_len = ...
    %         fDisplayInternalMessage('Starting Simulation', internal_msg_len);
    
    tmp_prog_txtlen = 0;
    for idx=1:nsamp
        Circuit(idx) = fMacSpiceSim(N, vs(:, idx), MemR, LRowR, LColR);
        
        tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
        tmp_prog_txtlen = fDisplayInternalMessage(...
            sprintf('fAlgorithm1: Simulation Progress: %2.2f percent', 100*(idx/nsamp)),...
            tmp_prog_txtlen);
    end
    fClearInternalMessages(tmp_prog_txtlen);
    
    
    %     internal_msg_len = ...
    %         fDisplayInternalMessage('Simulation Complete', internal_msg_len);
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
% load(filename)

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
            spectrumIdx = (minIdx == min(minIdx, [], 2));
            try
                fmag(rowIdx, colIdx) = freqSpectrum(spectrumIdx);
                fval(rowIdx, colIdx) = fsamp * norm_freq(spectrumIdx);
            catch
                [tmpVal, tmpIdx] = max(freqSpectrum(spectrumIdx));
                fmag(rowIdx, colIdx) = tmpVal;
                tmpVal = norm_freq(spectrumIdx);
                fval(rowIdx, colIdx) = fsamp * tmpVal(tmpIdx);
            end
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
                if sum(isnan(CalMemR(:)))~=0 || sum(CalMemR(:)<=0)~=0
                    CalMemR(isnan(CalMemR)) = abs(prevCalMemRinner(isnan(CalMemR)));
                    CalMemR(CalMemR<0) = prevCalMemRinner(CalMemR<0);
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

%% Read Memory Array
[readBits, CalMemR] = fReadFromMemArray(N, numBitspRes, CalMemR, ruleNum, var);

%% Estimate Line Resistances
allLREstimate = zeros(N);
for i = 1:N
    for j = 1:N
        nCF = fNodeCurrentFactors([], N, i, j, 0, abs(CalMemR), LRowR, LColR);
        sumnCF = sum(nCF(i:end, 1));
        
        allLREstimate(i, j) = (Vbias(i, j)/Imem(i, j) - MemR(i, j))/sumnCF;
    end
end
invalid = false(N, N);
% Get rid of invalid values before averaging out
invalid(allLREstimate<=0) = true;
invalid(allLREstimate>(1+LR_maxdev)*LRmean) = true;
invalid(allLREstimate<(1-LR_maxdev)*LRmean) = true;

validLREstimate = allLREstimate(~invalid);
LREstimate = mean(validLREstimate(:));

%% Calculate Bit Error Rate
BER = sum(storedBits~=readBits)/length(readBits);
% fprintf("Number Of Bits In Error is %d out of %d bits.\n", (BER * N^2 * numBitspRes), N^2 * numBitspRes)

%% Clear all internal messages
% fClearInternalMessages(internal_msg_len);

end