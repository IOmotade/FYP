% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 26/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts a simulated circuit from the time domain to frequency domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% N (Integer) = Order of memristor array i.e. NxN array
% numBitspRes(Integer) = Number of bits represented by a single memristor
% in the array
% MemR (NxN Double) = Memristor array values
% rulenum (Integer) = Rule used for writing/reading to/from memristor array
% MimMemRValue (Double) = Minimum Possible Memristor Resistance Value
% MaxMemRValue (Double) = Maximum Possible Memristor Resistance Value
% bitsStored (1xR Integer) = Bits originally stored in array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% readBits (Circuit) = Bits read from memristor array
% closestMemR (NxN Double) = Best memristor array approximation for stored
% bits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [readBits, closestMemR] = fReadFromMemArray(N, numBitspRes, MemR, ruleNum, varargin)
internal_msg_len = 0;

%% Calculate Total Memory Size
% memLength = N^2 * numBitspRes;
if ~isempty(varargin)
    var = squeeze(varargin{:});
else
    var = [];
end

%% Read possible bits to cut down on processing time in fReadFromMemArray/fWriteFromMemArray
filePath = 'Functions/Circuit_Sim_Functions/Encoding_Rules/';
fileLoc = sprintf('%sposBits_numBitspRes_%d.txt', filePath, numBitspRes);
fid = fopen(fileLoc, 'r');
if fid~=-1
    posBits = reshape(fread(fid, 'uint8=>char'), [], numBitspRes);
    fclose(fid);
else
    internal_msg_len = fDisplayInternalMessage...
        ('Generating Gray Code Encoding', internal_msg_len);
    
    %% Generate bit combinations and mapping
    posBits = (1:2^numBitspRes)-1; %Stored Bit Words
    posBits = fBin2Gray(dec2bin(posBits));
    
    internal_msg_len = fDisplayInternalMessage...
        ('Gray Code Encoding Generation Complete', internal_msg_len);
    
    %% Save possible bits to cut down on processing time in fReadFromMemArray/fWriteFromMemArray
    fid = fopen(fileLoc, 'w');
    fwrite(fid, posBits);
    fclose(fid);
    internal_msg_len = fDisplayInternalMessage...
        ('Gray Code Encoding Saved To File', internal_msg_len);
    
end
posSymbs = bin2dec(posBits);

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
        posMemRValues = log10(posMemRValues);
        MemR = log10(MemR);
    elseif ruleNum==3
        %Rule 1
        %rule1(numBitspRes, MemRLowLim, MemRUpLim);
        if length(var)==2
            posMemRValues = rule2(numBitspRes, var{1}, var{2});
        else
            posMemRValues = rule2(numBitspRes);
        end
    elseif ruleNum==4
        if length(var)==2
            posMemRValues = rule4(numBitspRes, var{1}, var{2});
        else
            posMemRValues = rule4(numBitspRes);
        end
    end
end

%% Estimate Memristor Array State and Stored Bits
closestMemR = MemR(:);
closestMemR = abs(closestMemR - posMemRValues);
[~, idx] = min(closestMemR, [], 2);
closestMemR = posMemRValues(idx);

[idx, ~] = find((closestMemR(:)==posMemRValues)');
readSymbs = posSymbs(idx);
readBits = dec2bin(readSymbs, numBitspRes);
readBits = str2num(readBits(:));

closestMemR = reshape(closestMemR, N, N);

%% Final adjustments to values
if ruleNum==2
    closestMemR = 10.^closestMemR;
end

fClearInternalMessages(internal_msg_len);
end