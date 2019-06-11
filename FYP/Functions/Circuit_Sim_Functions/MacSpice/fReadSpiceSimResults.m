% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 26/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read MacSpice result files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% fileLoc (String) = Full Absolute/Relative Filepath and Name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% Circuit (Circuit) = Simulated/Resulting Circuit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Circuit = fReadSpiceSimResults(fileLoc, Circuit)
%% Create temporary file
[filePath, name, ext] = fileparts(char(fileLoc));
tmpFileLoc = fullfile(filePath, strcat('tmp', name, ext));
copyfile(char(fileLoc), char(tmpFileLoc))

%% Pre-process File
pre_process(fileLoc);

%% Extract Data
SimData = extractData(fileLoc);

%% Correlate Simulation Data with Circuit
Circuit = fCorrelateCD(Circuit, SimData);

end

%% Main Sub-Functions
%% Pre-Process Simulation Results File for further processing
function pre_process(fileLoc)
if ismac
    %% MacSpice
    unwantedLineStarts = [
        "Cir";
        "Date";
        "Oper";
        "Node";
        "---";
        "Sou";
        "Sou";
        "Isourc";
        "Vsource";
        "Resist";
        "mode";
        "rsh";
        "narr";
        "tc";
        "def";
        "resistan";
        "p";
        "acm";
        "dc";
        "CPU";
        "Tota";
        "MacSp";
        ];
    
elseif ispc
    %% ngspice
    unwantedLineStarts = [
        "Doing,";
        "No.";
        "short";
        "l, ";
        "kf";
        "af";
        "lf";
        "wf";
        "ef";
        "m, ";
        "sin";
        "sffm";
        "am";
        "trnoise";
        "trrandom";
        "current";
        "ac";
        "bv_max";
        "noisy";
        ];
end

fRemoveExtraWhitespace(fileLoc);
fRemoveLines_v1(fileLoc, unwantedLineStarts)

end

%% Extract Data from Simulation Results File
function SimData = extractData(fileLoc)
%% Open file
fileID = fopen(fileLoc);

%% Read in all node voltages
data = textscan(fileID, '%s %f');
nodenames = strrep(string(data{:, 1}), ',', '');
nodevalues = data{:, 2};
nodenames = nodenames(1:length(nodevalues));

%% Read in all devices & their respective currents
%reload file
fclose(fileID);
fileID = fopen(fileLoc);

%remove already processed lines
numOfLines = length(nodevalues);
fRemoveLines_v2(fileID, numOfLines);

%further extract data
% data = textscan(fileID, '%s %s %s %s');
% ProcessedData = string(zeros(numel(data{:,1}), numel(data)));
% for idx = 1:numel(data)
%     ProcessedData(:,idx) = strrep(string(data{:,idx}), ',', '');
% end

numCol = 2;
lineIdx = 1;
pttrn = '%s ';
ProcessedData = repmat(" ", 1, numCol);
text_in = fgetl(fileID);
while ischar(text_in)
    colLength = length(strfind(text_in,',')) + 1;
    %     linePttrn = repmat(pttrn, 1, colLength);
    data = textscan(text_in, pttrn);
    if mod(lineIdx, numCol) == 1
        tmpProcessedData = repmat(" ", colLength, numCol);
    end
    tmpProcessedData(:, (mod(lineIdx-1, numCol)+1)) = ...
        strrep(string(data{:}), ',', '');
    if mod(lineIdx, numCol) == 0
        ProcessedData = [ProcessedData; tmpProcessedData(2:end, :)];
    end
    text_in = fgetl(fileID);
    lineIdx = lineIdx + 1;
end

%Process extracted data
ProcessedData = ProcessedData(2:end, :);
devnames = ProcessedData(:, 1);
devvalues_str = ProcessedData(:, 2);

devnames = devnames(:);
devvalues_str = devvalues_str(:);

devvalues = zeros(size(devvalues_str));
for idx = 1:length(devvalues_str)
    devvalues(idx) = sscanf(devvalues_str(idx), '%f');
end

%% Add in ground node
nodenames(end+1) = "0";
nodevalues(end+1) = 0;

%% Add in non-existent devices i.e. no current
devnames(end+1) = " ";
devvalues(end+1) = 0;

%% Collate Data
SimData.paramname = [nodenames; devnames];
SimData.value = [nodevalues; devvalues];

%% Close File
fclose(fileID);

end

%% Put Extracted Data into Circuit Structure
function Circuit = fCorrelateCD(Circuit, SimData)
%% Output Voltage Matrix
[numRows, numCols] = size(Circuit.VO.nodename);
for rowIdx = 1:numRows
    for colIdx = 1:numCols
        Circuit.VO.value(rowIdx, colIdx) = SimData.value(...
            strcmp(SimData.paramname, Circuit.VO.nodename(rowIdx, colIdx)));
    end
end

%% Input Voltage Matrix
[numRows, numCols] = size(Circuit.VI.nodename);
for rowIdx = 1:numRows
    for colIdx = 1:numCols
        Circuit.VI.value(rowIdx, colIdx) = SimData.value(...
            strcmp(SimData.paramname, Circuit.VI.nodename(rowIdx, colIdx)));
    end
end

%% Output Current Matrix
[numRows, numCols] = size(Circuit.IO.devname);
for rowIdx = 1:numRows
    for colIdx = 1:numCols
        Circuit.IO.value(rowIdx, colIdx) = SimData.value(...
            strcmp(SimData.paramname, Circuit.IO.devname(rowIdx, colIdx)));
    end
end

%% Input Current Matrix
[numRows, numCols] = size(Circuit.II.devname);
for rowIdx = 1:numRows
    for colIdx = 1:numCols
        Circuit.II.value(rowIdx, colIdx) = SimData.value(...
            strcmp(SimData.paramname, Circuit.II.devname(rowIdx, colIdx)));
    end
end

%% Source Voltage Matrix
[numRows, numCols] = size(Circuit.VS.nodename);
for rowIdx = 1:numRows
    for colIdx = 1:numCols
        Circuit.VS.value(rowIdx, colIdx) = SimData.value(...
            strcmp(SimData.paramname, Circuit.VS.nodename(rowIdx, colIdx)));
    end
end

%% Source Current Matrix
[numRows, numCols] = size(Circuit.IS.devname);
for rowIdx = 1:numRows
    for colIdx = 1:numCols
        Circuit.IS.value(rowIdx, colIdx) = SimData.value(...
            strcmp(SimData.paramname, Circuit.IS.devname(rowIdx, colIdx)));
    end
end

%% X Voltage Matrix
[numRows, numCols] = size(Circuit.X.nodename);
for rowIdx = 1:numRows
    for colIdx = 1:numCols
        Circuit.X.value(rowIdx, colIdx) = SimData.value(...
            strcmp(SimData.paramname, Circuit.X.nodename(rowIdx, colIdx)));
    end
end

end

%% Auxiliary Functions
function fRemoveLines_v1(fileLoc, lineChars)
if size(fileLoc, 1) ~= size(lineChars, 1)
    fileLoc = repmat(fileLoc, size(lineChars));
end
arrayfun(@fRemoveLine, fileLoc, lineChars);
end

function fRemoveLine(fileLoc, lineChars)

[filePath, name, ext] = fileparts(char(fileLoc));
tmpFileLoc = fullfile(filePath, strcat('tmp', name, ext));
copyfile(char(fileLoc), char(tmpFileLoc))

expr = char(strcat(lineChars, '(\w*)'));

%Make temporary copy of file
copyfile(char(fileLoc), char(tmpFileLoc))

%Open Files
fid_in  = fopen(tmpFileLoc);
fid_out = fopen(fileLoc,'w+');

text_in = fgetl(fid_in);
while ischar(text_in)
    %Check if an empty line
    valid = ~isempty(text_in);
    
    %Remove unneccessary lines
    lenChars = length(char(lineChars)) + 1;
    if length(text_in)>=lenChars
        strtIdx = regexp(text_in(1:lenChars),expr, 'ONCE');
    else
        strtIdx = regexp(text_in, expr, 'ONCE');
    end
    
    valid = valid & isempty(strtIdx);
    if  valid
        fprintf(fid_out,'%s\n',text_in);
    end
    text_in = fgetl(fid_in);
end

fclose(fid_in);
fclose(fid_out);

delete(tmpFileLoc)

end

function fRemoveLines_v2(fileID, N)
for i = 1:N
    fgetl(fileID);
end
end

function fRemoveExtraWhitespace(fileLoc)

[filePath, name, ext] = fileparts(char(fileLoc));
tmpFileLoc = fullfile(filePath, strcat('tmp', name, ext));
copyfile(char(fileLoc), char(tmpFileLoc))

%Make temporary copy of file
copyfile(char(fileLoc), char(tmpFileLoc))

%Open Files
fid_in  = fopen(tmpFileLoc);
fid_out = fopen(fileLoc,'w+');

text_in = fgetl(fid_in);
while ischar(text_in)
    % Clear leading whitespace
    text_out = strtrim(text_in);
    
    % Remove other whitespaces
    text_out = regexprep(text_out,'(\t|\s+)',', ');
    fprintf(fid_out,'%s\n',text_out);
    text_in = fgetl(fid_in);
end

fclose(fid_in);
fclose(fid_out);

delete(tmpFileLoc)

end