%% Full-Simulation
% function Circuit = fMacSpiceSim(N, vs, is, connVec, MemR, LRowR, LColR)
function Circuit = fMacSpiceSim(N, varargin)
%% Set-up
[vs, is, connVec, MemR, LRowR, LColR] = fSetup(N, varargin);

%% Definition of different types of files
spiceFileName = 'mat_array.cir';
simOutputFile = 'out.txt';
tmpSimOutputFile = string(['tmp' simOutputFile]);

%% Generate Spice File
fileLoc = spiceFileName;
Circuit = fGenerateSpiceFile(N, vs, is, connVec, MemR, LRowR, LColR, fileLoc);

%% Perform Simulation Using Macspice
command = sprintf('/Applications/MacSpice.app/Contents/MacOS/MacSpice -b mat_array.cir > %s',...
    simOutputFile);
status = system(command);

%% Read Simulation Results
if status ==0
    copyfile(char(simOutputFile), char(tmpSimOutputFile))
    Circuit = fReadSpiceSimResults(tmpSimOutputFile, Circuit);
else
    fprintf('Spice simulation could not run');
end

%% Delete Temporary Files
delete(tmpSimOutputFile)

end

function [vs, is, connVec, MemR, LRowR, LColR] = fSetup(N, var)
%% Unpack Input Arguments
%Condition 1: vs, MemR, LRowR, LColR, scheme
if length(var)==4
    % Passed Arguments
    vs = var{1};
    MemR = var{2};
    LRowR = var{3};
    LColR = var{4};
    
    %Unpassed Arguments
    is = zeros(N, 1);
    connVec = ones(N, 1);
end

%Condition 2: vs, is, connVec, MemR, LRowR, LColR
if length(var)==6
    vs = var{1};
    is = var{2};
    connVec = var{3};
    MemR = var{4};
    LRowR = var{5};
    LColR = var{6};
end

%% Further Setup
% Default Column Source/Sink Values: 0
if length(vs)==N
    vs = [vs;
        zeros(N, 1)];
end

if length(is)==N
    is = [is;
        zeros(N, 1)];
end

% Default Connection: To Voltage Source/Sink
lnght = length(connVec);
if lnght<=(2*N)
    connVec = [connVec;
        ones(((2*N)-lnght), 1)];
end
end