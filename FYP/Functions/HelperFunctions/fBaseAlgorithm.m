function BER = fBaseAlgorithm(setup, N, numBitspRes, LRdef, ruleNum, varargin)
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
fsamp = oversampfactor*2*max(fsource);
nsamp = 1*(fsamp/base_freq);
% tsamp = 1/fsamp;
% t = 0:tsamp:(nsamp-1)*tsamp;

%% Sim Setup
MemR = storedMemR;

vs_mag = 5;
% vs = vs_mag*square(2*pi*fsource*t);
vs = fVoltageSourceSignals(N, vs_mag, nsamp);
Circuit = fMacSpiceSim(N, vs(:, 1), MemR, LRowR, LColR);
Circuit = repmat(Circuit, [nsamp, 1]);

%% Run Simulation
if ~setup.perfect
    internal_msg_len = ...
        fDisplayInternalMessage('Starting Simulation\n', internal_msg_len);
    
    tmp_prog_txtlen = 0;
    for idx=1:nsamp
        Circuit(idx) = fMacSpiceSim(N, vs(:, idx), MemR, LRowR, LColR);
        
        tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
        tmp_prog_txtlen = fDisplayInternalMessage(...
            sprintf('Simulation Progress: %2.2f percent\n', 100*(idx/nsamp)),...
            tmp_prog_txtlen);
    end
    fClearInternalMessages(tmp_prog_txtlen);
    
    
    internal_msg_len = ...
        fDisplayInternalMessage('Simulation Complete\n', internal_msg_len);
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
I = fmag;
% Frequency_Components = fUnits(fmag, 'A')

%% Memristor Value Finding Algorithm
V = diag(vs_mag*ones(N, 1));
G = I/V;
CalMemR = 1./G;

%% Read Memory Array
[readBits, ~] = fReadFromMemArray(N, numBitspRes, CalMemR, ruleNum, var);

%% Calculate Bit Error Rate
BER = sum(storedBits~=readBits)/length(readBits);
% fprintf("Number Of Bits In Error is %d out of %d bits.\n", (BER * N^2 * numBitspRes), N^2 * numBitspRes)

%% Clear all internal messages
fClearInternalMessages(internal_msg_len);

end