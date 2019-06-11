function BER = fAlgorithm2(setup, N, numBitspRes, LRdef, ruleNum, varargin)
addpath(genpath('../FYP/Functions/'))
internal_msg_len = 0;

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
    setup.basefreq = 10;
    setup.oversampfactor = 8;
    setup.SNR = Inf;
end

if ~fFieldExist(setup, 'SNR')
    setup.SNR = Inf;
end

%% Setup Circuit System
LRowR = abs(LRmean*ones(N)  + LRsigma*randn(N));
LColR = abs(LRmean*ones(N)  + LRsigma*randn(N));

%% Fingerprint Set-up
fingerprint_filename = sprintf('Fingerprints/fingerprint_N_%d_LRmean_%2.0f_numBitspRes_%d_MemRange_%2.0f_%2.0f',...
    N, LRmean, numBitspRes, var{1}, var{2});
vs_mag = 5;
memLength = N^2 * numBitspRes;
idxLength = 2^memLength;

try
    load(fingerprint_filename, 'Fingerprint');
catch
    internal_msg_len = fDisplayInternalMessage...
        ('Generating Fingerprint Database', internal_msg_len);
    
    Fingerprint = zeros(idxLength, N, N);
    
    tmp_prog_txtlen = 0;
    for idx = 1:idxLength
        bitsToStore = dec2bin(idx-1, memLength);
        bitsToStore = bitsToStore(:);
        bitsToStore = str2num(bitsToStore);
        
        [fingerprintMemR, ~] = fWriteToMemArray(N, numBitspRes, ruleNum, bitsToStore, var);
        Fingerprint(idx, :, :) = fIdealOFDMSolution(N, vs_mag, fingerprintMemR, LRmean*ones(N), LRmean*ones(N));
        
        tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
        tmp_prog_txtlen = fDisplayInternalMessage(...
            sprintf('Completed %d out of %d Fingerprints', idx, idxLength),...
            tmp_prog_txtlen);
    end
    fClearInternalMessages(tmp_prog_txtlen);
    
    internal_msg_len = fDisplayInternalMessage...
        ('Generation Complete', internal_msg_len);
    
    save(fingerprint_filename, 'Fingerprint');
    
    internal_msg_len = fDisplayInternalMessage...
        ('Fingerprint Database Saved', internal_msg_len);
end

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

% vs = vs_mag*square(2*pi*fsource*t);
vs = fVoltageSourceSignals(N, vs_mag, nsamp);
vs = vs + (vs_mag/db2mag(setup.SNR))*randn(size(vs));
Circuit = fSpiceSim(N, vs(:, 1), MemR, LRowR, LColR);
Circuit = repmat(Circuit, [nsamp, 1]);

%% Run Simulation
if ~setup.perfect
%     internal_msg_len = ...
%         fDisplayInternalMessage('Starting Simulation', internal_msg_len);
    
    tmp_prog_txtlen = 0;
    for idx=1:nsamp
        Circuit(idx) = fSpiceSim(N, vs(:, idx), MemR, LRowR, LColR);
        
        tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
        tmp_prog_txtlen = fDisplayInternalMessage(...
            sprintf('fAlgorithm2: Simulation Progress: %2.2f percent', 100*(idx/nsamp)),...
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

%% Work out best fit to circuit
distance = zeros(1, idxLength);
for idx = 1:idxLength
    val = fmag - squeeze(Fingerprint(idx, :, :));
    distance(idx) = norm(val, 'fro')^2;
end
[~, solutionIdx] = min(distance);

%% Read/Rework Bits Stored
readBits = dec2bin(solutionIdx-1, memLength);
readBits = readBits(:);
readBits = str2num(readBits);

%% Actual stored state
% pwr = (2.^(length(storedBits)-1:-1:0));
% correctIdx = sum(storedBits(:).*pwr(:))+1;

%% Read from Memory
% [readBits, ~] = fReadFromMemArray(N, numBitspRes, CalMemR, ruleNum, var);
BER = sum(storedBits~=readBits)/length(readBits);

%% Clear all internal messages
fClearInternalMessages(internal_msg_len);

end
