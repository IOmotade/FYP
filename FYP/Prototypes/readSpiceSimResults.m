%% Test line removal function
% rst
fileDir = '/Users/IOmotade/Desktop/School/4th_Year/FYP/Sims/';
fileName = string([fileDir 'out.txt']);
tmpFileName = string(['../FYP/' 'tmp.txt']);
unwantedLineStarts = ["Cir";
    'Date';
    'Oper';
    'Node';
    '---';
    'Sou';
    'Resist';
    'mode';
    'rsh';
    'narr';
    'tc';
    'def';
    'resistan';
    'p';
    'acm';
    'dc';
    'CPU';
    'Tota';
    ];

copyfile(char(fileName), char(tmpFileName))
fRemoveLines(tmpFileName, unwantedLineStarts)
fRemoveExtraWhitespace(tmpFileName);
tmpFileID = fopen(tmpFileName);

%% Memristor Nodes
data = textscan(tmpFileID, '%s %f');
nodenames = strrep(string(data{:, 1}), ',', '');
nodevalues = data{:, 2};
nodenames = nodenames(1:length(nodevalues));

%% Devices & Currents
%reload file
fclose(tmpFileID);
tmpFileID = fopen(tmpFileName);

%remove already processed lines
numOfLines = length(nodevalues);
fRemoveLines_v2(tmpFileID, numOfLines);

%further extract data
data = textscan(tmpFileID, '%s %s %s %s');
% devnames = strrep(string(data{:, 1}), ',', '');
ProcessedData = string(zeros(numel(data{:,1}), numel(data)));
for idx = 1:numel(data)
    ProcessedData(:,idx) = strrep(string(data{:,idx}), ',', '');
end

%Process extracted data
ProcessedData = ProcessedData(:, 2:end);
devnames = ProcessedData(1:2:end, :);
devvalues_str = ProcessedData(2:2:end, :);

devnames = devnames(:);
devvalues_str = devvalues_str(:);

devvalues = zeros(size(devvalues_str));
for idx = 1:length(devvalues_str)
   devvalues(idx) = sscanf(devvalues_str(idx), '%f');
end

SimResults.paramname = [nodenames; devnames];
SimResults.values = [nodevalues; devvalues];

fclose(tmpFileID);