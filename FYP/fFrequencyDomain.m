function Circuit = fFrequencyDomain(timeCircuit)

%% Set-up
freqCircuit = timeCircuit;
nSamp = length(timeCircuit);
N = size(timeCircuit(1).VO.value, 1);

%% VO
timeVal = zeros([1 length(timeCircuit)]);
for rowIdx = 1:N
    for colIdx = 1:N
        for tIdx=1:nSamp
            timeVal(tIdx) = timeCircuit(tIdx).VO.value(rowIdx, colIdx);
        end
        freqVal = fHaarT(timeVal);
        for fIdx=1:nSamp
            freqCircuit(fIdx).VO.value(rowIdx, colIdx) = freqVal(fIdx);
        end
    end
end

%% VI
timeVal = zeros([1 length(timeCircuit)]);
for rowIdx = 1:N
    for colIdx = 1:N
        for tIdx=1:nSamp
            timeVal(tIdx) = timeCircuit(tIdx).VI.value(rowIdx, colIdx);
        end
        freqVal = fHaarT(timeVal);
        for fIdx=1:nSamp
            freqCircuit(fIdx).VI.value(rowIdx, colIdx) = freqVal(fIdx);
        end
    end
end

%% IO
timeVal = zeros([1 length(timeCircuit)]);
for rowIdx = 1:N
    for colIdx = 1:N
        for tIdx=1:nSamp
            timeVal(tIdx) = timeCircuit(tIdx).IO.value(rowIdx, colIdx);
        end
        freqVal = fHaarT(timeVal);
        for fIdx=1:nSamp
            freqCircuit(fIdx).IO.value(rowIdx, colIdx) = freqVal(fIdx);
        end
    end
end

%% II
timeVal = zeros([1 length(timeCircuit)]);
for rowIdx = 1:N
    for colIdx = 1:N
        for tIdx=1:nSamp
            timeVal(tIdx) = timeCircuit(tIdx).II.value(rowIdx, colIdx);
        end
        freqVal = fHaarT(timeVal);
        for fIdx=1:nSamp
            freqCircuit(fIdx).II.value(rowIdx, colIdx) = freqVal(fIdx);
        end
    end
end

%% VS
timeVal = zeros([1 length(timeCircuit)]);
for rowIdx = 1:N
    for colIdx = 1:N
        for tIdx=1:nSamp
            timeVal(tIdx) = timeCircuit(tIdx).VS.value(rowIdx, colIdx);
        end
        freqVal = fHaarT(timeVal);
        for fIdx=1:nSamp
            freqCircuit(fIdx).VS.value(rowIdx, colIdx) = freqVal(fIdx);
        end
    end
end

%% IS
timeVal = zeros([1 length(timeCircuit)]);
for rowIdx = 1:N
    for colIdx = 1:N
        for tIdx=1:nSamp
            timeVal(tIdx) = timeCircuit(tIdx).IS.value(rowIdx, colIdx);
        end
        freqVal = fHaarT(timeVal);
        for fIdx=1:nSamp
            freqCircuit(fIdx).IS.value(rowIdx, colIdx) = freqVal(fIdx);
        end
    end
end

%% X
timeVal = zeros([1 length(timeCircuit)]);
for rowIdx = 1:N
    for colIdx = 1:N
        for tIdx=1:nSamp
            timeVal(tIdx) = timeCircuit(tIdx).X.value(rowIdx, colIdx);
        end
        freqVal = fHaarT(timeVal);
        for fIdx=1:nSamp
            freqCircuit(fIdx).X.value(rowIdx, colIdx) = freqVal(fIdx);
        end
    end
end

%% Delete Redundant Information
freqCircuit = rmfield(freqCircuit,'Memristor');
freqCircuit = rmfield(freqCircuit,'RowLineResistance');
freqCircuit = rmfield(freqCircuit,'ColLineResistance');

%% Finalize Structure
Circuit.TimeDom = timeCircuit;
Circuit.FreqDom = freqCircuit;

end