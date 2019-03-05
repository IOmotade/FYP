function Circuit = fGenerateSpiceFile(N, vs, MemR, LRowR, LColR, fileLoc, scheme, spiceDir)
if N~=length(vs)
    fprintf('Dimensions do not match');
    return
end

if(~exist('spiceDir', 'var'))
    spiceDir = '.op\n';
end

if ~exist('scheme', 'var')
    scheme = '1';
end

%% Generate Matrices
[Vo, Io, Vs, Is, X] = deal(zeros(N));
Vs(:, 1) = vs;

%% Write-up Spice File
Circuit = fWriteToSpiceFile(N, MemR, LRowR, LColR, Vo, Io, Vs, Is, X, fileLoc, scheme, spiceDir);

end

function Circuit = fWriteToSpiceFile(N, MemR, LRowR, LColR, Vo, Io, Vs, Is, X_, fileLoc, scheme, spiceDir)
R_Amm = 1e-9; %Ohms Used as Ammeter

%% Definition of Component Base Names
vCompName = "V";
memRCompName = "MemR";
lRRCompName = "LRR";
lCRCompName = "LCR";

%Resistances used as Ammetters
rarCompName = "RAR_V";  %Row
carCompName = "CAR_I";  %Column

rowTag = "Row";
colTag = "Col";

%% Generate component aliases
res.tag = [memRCompName lRRCompName lCRCompName rarCompName carCompName];
res.alias = res.tag;
for idx = 1:length(res.tag)
    res.alias(idx) = sprintf('R%d', idx);
end

%% Variables storing list of all spice node names
%%Measurement Values
%Output Voltage Matrix
VO.value = Vo;
VO.nodename = repmat(" ", size(Vo));

%Input Voltage Matrix
VI.value = Vo;
VI.nodename = repmat(" ", size(Vo));

%Output Current Matrix
IO.value = Io;
IO.devname = repmat(" ", size(Io));

%Output Current Matrix
II.value = Io;
II.devname = repmat(" ", size(Io));

%Source Voltage Matrix
VS.value = Vs;
VS.nodename = repmat(" ", size(Vs));

%Source Current Matrix
IS.value = Is;
IS.devname = repmat(" ", size(Is));

%X Voltage Matrix
X.value = X_;
X.nodename = repmat(" ", size(Is));

%%Device Properties
%Memristor
Memristor.value = MemR;
Memristor.devname = repmat(" ", size(MemR));

%Row Line Resistance
RowLineResistance.value = LRowR;
RowLineResistance.devname = repmat(" ", size(LRowR));

%Column Line Resistance
ColLineResistance.value = LColR;
ColLineResistance.devname = repmat(" ", size(LColR));

%% Generation of spice script
fileID = fopen(fileLoc, 'w+');
%%Comment: File Location & Name
fprintf(fileID, "*%s\n", fileLoc);

%% Write Voltage Sources into file
%%Comment: Voltage Sources
fprintf(fileID, "\n\n*Voltage Sources\n");

%%Generate Voltage Source Names
for rowIdx = 1:N
    elemName = sprintf('%s%d', vCompName, rowIdx);
    posNode = sprintf('%s_%d', vCompName, rowIdx);
    negNode = '0';
    value = Vs(rowIdx, 1); %vs(rowIdx);
    text = sprintf('%s %s %s %f\n', elemName, posNode, negNode, value);
    fprintf(fileID, text);
    
    VS.nodename(rowIdx, 1) = lower(string(posNode));
end
VS.nodename(:, 2:end) = lower(string(negNode));

%% Write Row Ammeter Resistances
%%Comment: Row Ammeter Resistances
fprintf(fileID, "\n\n*Row 'Ammeter' Resistances\n");

%%Generate Row Ammeter Resistances
colIdx = 1;
compName = res.alias(strcmp(res.tag, rarCompName));
for rowIdx = 1:N
    elemName = sprintf("%s%d", compName, rowIdx);
    posNode = sprintf('%s_%d', vCompName, rowIdx);
    negNode = sprintf('%s%s%d', memRCompName, rowTag, fPairFunction(rowIdx, colIdx));
    text = sprintf('%s %s %s %.12f\n', elemName, posNode, negNode, R_Amm);
    fprintf(fileID, text);
    
    IS.devname(rowIdx, 1) = lower(string(elemName));
end

%% Write Memristor Resistance values into file
%%Comment: Memristor Resistances
fprintf(fileID, "\n\n*Memristor Resistances\n");

%%Generate Memristor Resistances
compName = res.alias(strcmp(res.tag, memRCompName));
for rowIdx = 1:N
    for colIdx = 1:N
        elemName = sprintf('%s%d', compName, fPairFunction(rowIdx, colIdx));
        posNode = sprintf('%s%s%d', memRCompName, rowTag, fPairFunction(rowIdx, colIdx));
        negNode = sprintf('%s%s%d', memRCompName, colTag, fPairFunction(rowIdx, colIdx));
        value = MemR(rowIdx, colIdx);
        text = sprintf('%s %s %s %f\n', elemName, posNode, negNode, value);
        fprintf(fileID, text);
        
        Memristor.devname(rowIdx, colIdx) = lower(string(elemName));
        VI.nodename(rowIdx, colIdx) = lower(string(posNode));
    end
end

%% Write Line Row Resistance values into file
%%Comment: Line Row Resistances
fprintf(fileID, "\n\n*Line Row Resistances\n");

%%Generate Line Row Resistances
compName = res.alias(strcmp(res.tag, lRRCompName));
for rowIdx = 1:N
    for colIdx = 1:N
        elemName = sprintf('%s%d', compName, fPairFunction(rowIdx, colIdx));
        posNode = sprintf('%s%s%d', memRCompName, rowTag, fPairFunction(rowIdx, colIdx));
        negNode = sprintf('%s%s%d', memRCompName, rowTag, fPairFunction(rowIdx, colIdx+1));
        value = LRowR(rowIdx, colIdx);
        text = sprintf('%s %s %s %f\n', elemName, posNode, negNode, value);
        fprintf(fileID, text);
        
        RowLineResistance.devname(rowIdx, colIdx) = lower(string(elemName));
        if colIdx==N
            VO.nodename(rowIdx, N) = lower(string(negNode));
        end
    end
end
VO.nodename(:, 1:N-1) = VI.nodename(:, 2:N);
IS.devname(:, 2:N) = RowLineResistance.devname(:, 1:N-1);

%% Write Line Column Resistance values into file
%%Comment: Line Column Resistances
fprintf(fileID, "\n\n*Line Column Resistances\n");

%%Generate Line Column Resistances
compName = res.alias(strcmp(res.tag, lCRCompName));
for rowIdx = 1:N
    for colIdx = 1:N
        elemName = sprintf('%s%d', compName, fPairFunction(rowIdx, colIdx));
        posNode = sprintf('%s%s%d', memRCompName, colTag, fPairFunction(rowIdx, colIdx));
        negNode = sprintf('%s%s%d', memRCompName, colTag, fPairFunction(rowIdx+1, colIdx));
        value = LColR(rowIdx, colIdx);
        text = sprintf('%s %s %s %f\n', elemName, posNode, negNode, value);
        fprintf(fileID, text);
        
        ColLineResistance.devname(rowIdx, colIdx) = lower(string(elemName));
        X.nodename(rowIdx, colIdx) = lower(string(posNode));
    end
end
IO.devname = ColLineResistance.devname;
II.devname(2:N, :) = IO.devname(1:N-1, :);

%% Write Column Ammeter Resistances
%%Comment: Column Ammeter Resistances
fprintf(fileID, "\n\n*Column 'Ammeter' Resistances\n");

%%Generate Column Ammeter Resistances
compName = res.alias(strcmp(res.tag, carCompName));
rowIdx = N+1;
for colIdx = 1:N
    elemName = sprintf('%s%d', compName, colIdx);
    posNode = sprintf('%s%s%d', memRCompName, colTag, fPairFunction(rowIdx, colIdx));
    if strcmp(scheme, '1')
        negNode = '0';
    elseif strcmp(scheme, '2')
        negNode = sprintf('%s%s%d', carCompName, colTag, fPairFunction(rowIdx+1, colIdx));
    end
    text = sprintf('%s %s %s %.12f\n', elemName, posNode, negNode, R_Amm);
    fprintf(fileID, text);
end

%% Write Spice Directives
%%Comment: Spice Directives
fprintf(fileID, "\n\n*Spice Directives\n");
fprintf(fileID, spiceDir);
fprintf(fileID, '*.backanno\n');
% fprintf(fileID, '.run\n');
fprintf(fileID, '.end\n');

%% Close File
fclose(fileID);

%% Collate all circuit information into one structure
Circuit.VO = VO;
Circuit.VI = VI;
Circuit.IO = IO;
Circuit.II = II;
Circuit.VS = VS;
Circuit.IS = IS;
Circuit.X = X;
Circuit.Memristor = Memristor;
Circuit.RowLineResistance = RowLineResistance;
Circuit.ColLineResistance = ColLineResistance;

end

function f = fPairFunction(a, b)
%Cantor Pairing Function
f = (1/2)*(a+b)*(a+b+1) + b;
end