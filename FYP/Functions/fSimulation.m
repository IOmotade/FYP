function meanBER = fSimulation(algo, setup, minNumBits)

%% Unpack Setup Structure
N = setup.N;
numBitspRes = setup.numBitspRes;
LRdef = setup.LRdef;
ruleNum = setup.ruleNum;
var = setup.var;

%% Algorithm Setup Variable
algosetup.perfect = false;
algosetup.basefreq = 10;
algosetup.oversampfactor = 1;

%% General Setup
bitsPerIter = N^2 * setup.numBitspRes;
idxLength = ceil(minNumBits/bitsPerIter);

%% Create necessary variables
BER = zeros(1, idxLength);

%% Run Simulations
if algo==0
    for idx = 1:idxLength
        %fAlgorithm1(setup, N, numBitspRes, LRdef, ruleNum, varargin)
        BER(idx) = fBaseAlgorithm(algosetup, N, numBitspRes, LRdef, ruleNum, var{1}, var{2});
    end
elseif algo==1
    for idx = 1:idxLength
        %fAlgorithm1(setup, N, numBitspRes, LRdef, ruleNum, varargin)
        BER(idx) = fAlgorithm1(algosetup, N, numBitspRes, LRdef, ruleNum, var{1}, var{2});
    end
elseif algo==2
    for idx = 1:idxLength
        %fAlgorithm2(setup, N, numBitspRes, LRdef, ruleNum, varargin)
        BER(idx) = fAlgorithm2(algosetup, N, numBitspRes, LRdef, ruleNum, var{1}, var{2});
    end
% elseif algo==3
end

%% Calculate mean BER
meanBER = mean(BER);

end