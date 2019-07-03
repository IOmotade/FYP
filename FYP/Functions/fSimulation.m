% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 28/06/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Carry out array circuit time simulations using OFDM method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% algo (1xA Double) = Algorithms to be evaluated
% setup (Struct) = Circuit Simulation Setup 
% minNumBits (Integer) = Minimum number of bits to be evaluated per system
% and algorithm
% readFromFile (Double) = Option to read from file (0/1/2 => New
% Simulations/ Read from Time Simulations/ Read from Ideal Simulation Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% meanBER (1xA Double) = Mean Bit Error Rate for all algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function meanBER = fSimulation(algo, setup, minNumBits, readFromFile)
%% Unpack Setup Structure
N = setup.N;
numBitspRes = setup.numBitspRes;
LRdef = setup.LRdef;
ruleNum = setup.ruleNum;
var = setup.var;
vs_mag = 5;

%% Array Simulation Setup Variable
simsetup.perfect = false;   %Get frequency components without using time sims
simsetup.basefreq = 10;     %Signal Vector Base Frequency
simsetup.oversampfactor = 1;    %Oversampling factor
if fFieldExist(setup, 'SNR')
    simsetup.SNR = setup.SNR;
else
    simsetup.SNR = Inf;
end

if ~isinf(simsetup.SNR)
    readFromFile = 1;
end

% Mean, Standard Deviation of and of Mean of Line Resistance
if length(LRdef)==2
    [LRmean, LRsigma, LRmeansigma] = deal(LRdef(1), LRdef(2)*LRdef(1), 0);
elseif length(LRdef)==3
    [LRmean, LRsigma, LRmeansigma] = deal(LRdef(1), LRdef(2)*LRdef(1),...
        LRdef(3)*LRdef(1));
else
    [LRmean, LRsigma, LRmeansigma] = deal(LRdef(1), 0, 0);
end

%% Decide if to read from saved file if possible
if ~exist('readFromFile', 'var')
    readFromFile = 1;
end

% To cut down on simulation time
if N>8
    readFromFile = 2;
end

%% Setup Circuit System
mean_offset = LRmeansigma*randn(1);
LRowR = abs(LRmean*ones(N)  + LRsigma*randn(N) + mean_offset);
LColR = abs(LRmean*ones(N)  + LRsigma*randn(N) + mean_offset);

%% General Setup
bitsPerIter = N^2 * setup.numBitspRes;
idxLength = ceil(minNumBits/bitsPerIter);
algoIdxLength = length(algo);

%% Create necessary variables
BER = zeros(algoIdxLength, idxLength);

%% Run Simulations
tmp_prog_txtlen = 0;
for idx = 1:idxLength
    %% Generate Sim Results Filename
    fileLoc = fSimResultsFilename(N, idx, numBitspRes, var, LRmean, LRsigma, simsetup.SNR, ruleNum);
    
    %% Decide if to run time simulation or read from stored files
    if (readFromFile~=2) || ~isinf(simsetup.SNR)
        if exist(char([fileLoc, '.mat']), 'file') && readFromFile==1
            %% Load required variables from file
            load(fileLoc,...
                'storedMemR', 'storedBits',...
                'CircuitSim', 'nsamp', 'fsource', 'fsamp');
        else
            %% Write to Memory
            [storedMemR, storedBits] = fWriteToMemArray(N, numBitspRes, ruleNum, [], var);
            
            %% Run Array Simulation
            [CircuitSim, nsamp, fsource, fsamp] = fArraySim(N, vs_mag, storedMemR, LRowR, LColR, ...
                simsetup.basefreq, simsetup.oversampfactor, simsetup.SNR, 1);
            
            %% Save Time Simulation Results
            save(fileLoc,...
                'storedMemR', 'storedBits',...
                'CircuitSim', 'nsamp', 'fsource', 'fsamp');
        end
        
        %% Get 'Frequency' Transform
        timeCircuitSim = CircuitSim;
        CircuitSim = fFrequencyDomain(timeCircuitSim);
        
        %% Extract Frequency Components from Transform
        [fmag, ~] = fExtractFrequencyComponents(N, CircuitSim, nsamp, fsource, fsamp);
    else
        try
            %%
            load(fileLoc,...
                'storedMemR', 'storedBits',...
                'CircuitSim', 'nsamp', 'fsource', 'fsamp');
            
            %% Get 'Frequency' Transform
            timeCircuitSim = CircuitSim;
            CircuitSim = fFrequencyDomain(timeCircuitSim);
            
            [fmag, ~] = fExtractFrequencyComponents(N, CircuitSim, nsamp, fsource, fsamp);
        catch
            try
                load([fileLoc '_fmag'],...
                    'storedMemR', 'storedBits', 'fmag');
            catch
                %% Write to Memory
                [storedMemR, storedBits] = fWriteToMemArray(N, numBitspRes, ruleNum, [], var);
                
                fmag = fIdealOFDMSolution(N, vs_mag, storedMemR, LRowR, LColR);
                
                save([fileLoc '_fmag'],...
                    'storedMemR', 'storedBits', 'fmag');
            end
        end
    end
    for algoIdx = 1:algoIdxLength
        switch algo(algoIdx)
            case 0
                CalMemR = fBaseAlgorithm(N, vs_mag, fmag);
            case 1
                [CalMemR, ~] = fAlgorithm1(N, vs_mag, fmag, LRmean, numBitspRes, ruleNum, var);
            case 1.1
                [CalMemR, ~] = fAlgorithm1_1(N, vs_mag, fmag, LRmean, numBitspRes, ruleNum, var);
            case 1.2
                [CalMemR, ~] = fAlgorithm1_2(N, vs_mag, fmag, LRmean, numBitspRes, ruleNum, var);
            case 2
                %[CalMemR, ~] =
        end
        [readBits, ~] = fReadFromMemArray(N, numBitspRes, CalMemR, ruleNum, var);
        BER(algoIdx, idx) = sum(storedBits~=readBits)/length(readBits);
        
        [prog_idx, prog_len] = deal((algoIdx+(idx-1)*algoIdxLength), idxLength*algoIdxLength);
        tmp_prog_txtlen = fUpdateProgress(prog_idx, prog_len, tmp_prog_txtlen);
    end
end

%% Calculate mean BER
meanBER = mean(BER, 2);

%% Clear Remaining Text
fClearInternalMessages(tmp_prog_txtlen);

end

%% Generate unique filenames to store results
function fileLoc = fSimResultsFilename(N, idx, numBitspRes, var, LRmean, LRsigma, SNR, ruleNum)
%  fileLoc = sprintf('SimData/fSimulation/CircuitSim_N_%d_idx_%d_LRmean_%d_numBitspRes_%d_minNumBits_%d',...
%         N, idx, LRmean, numBitspRes, minNumBits);
fileLoc = sprintf('SimData/fSimulation/CircuitSim_N_%d_idx_%d_LRmean_%d_numBitspRes_%d',...
    N, idx, LRmean, numBitspRes);

%% Check for special file cases
if ruleNum~=1
    fileLoc = char([fileLoc, sprintf('_ruleNum_%d', ruleNum)]);
end

if ~isinf(SNR)
    fileLoc = char([fileLoc, sprintf('_SNR_%1.0f', simsetup.SNR)]);
    readFromFile = 1;
end

if var{1}~=1e3
    fileLoc = char([fileLoc, sprintf('_MemRMin_%1.0f', var{1})]);
end

if var{2}~=1e6
    fileLoc = char([fileLoc, sprintf('_MemRMax_%1.0f', var{2})]);
end

if LRmean~=100
    fileLoc = char([fileLoc, sprintf('_LRmean_%1.0f', LRmean)]);
end

if LRsigma~=0
    fileLoc = char([fileLoc, sprintf('_LRsigma_%dpc', LRsigma*100)]);
end
end

%% Push progress text to screen
function tmp_prog_txtlen = fUpdateProgress(idx, idxLength, tmp_prog_txtlen)
tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
tmp_prog_txtlen = fDisplayInternalMessage(...
    sprintf('fSimulation: Simulation Progress: %2.2f percent', 100*(idx/idxLength)),...
    tmp_prog_txtlen);
end