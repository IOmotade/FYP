% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 26/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts a simulated circuit from the time domain to frequency domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% N (Integer) = Order of memristor array i.e. NxN array
% numBitspRes(Integer) = Number of bits represented by a single memristor
% in the array
% rulenum (Integer) = Rule used for writing/reading to/from memristor array
% default value is 1
% bitsToStore (Rx1 Integers) = Bits(0/1s) to store in the memristor array
% (randomly generated otherwise)
% varargin = input variables for the rule to be used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% MemR (NxN Double) = Memristor array values implemented for stored bits
% storedBits (Circuit) = Bits read from memristor array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MemR, storedBits] = fWriteToMemArray(N, numBitspRes, ruleNum, bitsToStore, varargin)
%% Calculate Total Memory Size
memLength = N^2 * numBitspRes;
if ~isempty(varargin)
    var = squeeze(varargin{:});
else
    var = [];
end

%% Generate bit combinations and mapping
posBits = (1:2^numBitspRes)-1; %Stored Bit Words
posBits = fBin2Gray(dec2bin(posBits));
posSymbs = bin2dec(posBits);

%% Save possible bits to cut down on processing time in fReadFromMemArray
filename = sprintf('posBits_numBitspRes_%d.txt', numBitspRes);
fid = fopen(filename, 'w');
fwrite(fid, posBits);

%% Depending on Encoding Rule Used
if ~exist('ruleNum', 'var')
    %Rule 1
    %rule1(numBitspRes, MemRLowLim, MemRUpLim);
    ruleNum = 1;
    if length(var)==2
        posMemRValues = rule1(numBitspRes, var{1}, var{2});
    else
        posMemRValues = rule1(numBitspRes);
    end
else
    if ruleNum==1
        %Rule 1
        %rule1(numBitspRes, MemRLowLim, MemRUpLim);
        if length(var)==2
            posMemRValues = rule1(numBitspRes, var{1}, var{2});
        else
            posMemRValues = rule1(numBitspRes);
        end
    elseif ruleNum==2
        %Rule 1
        %rule1(numBitspRes, MemRLowLim, MemRUpLim);
        if length(var)==2
            posMemRValues = rule2(numBitspRes, var{1}, var{2});
        else
            posMemRValues = rule2(numBitspRes);
        end
    end
end

%% Instatiate Stored Values
if ~exist('bitsToStore', 'var')
    storedBits = round(rand(memLength, 1));
else
    if isempty(bitsToStore)
        storedBits = round(rand(memLength, 1));
    else
        storedBits = bitsToStore;
    end
end
tmp = reshape(storedBits, memLength/numBitspRes, numBitspRes);
storedSymbs = bin2dec(num2str(tmp));

%% Assign Storage Values to Memristor Values
[rowIdx, ~] = find(posSymbs==storedSymbs.');
MemR = reshape(posMemRValues(rowIdx), N, N);

end