% Call cases:
% -----------
% function Circuit = fMacSpiceSim(N, vs, is, connVec, MemR, LRowR, LColR)
% function Circuit = fMacSpiceSim(N, s, connVec, MemR, LRowR, LColR)
%
% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 26/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Carries out a macspice simulation for a memristor array with parasitic
% line resistances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% N (Integer) = Order of memristor array i.e. NxN array
% s (Nx1 or 2Nx1 Double) = General Bias Source/Sink values
% vs (Nx1 or 2Nx1 Double) = Voltage Source/Sink values
% is (Nx1 or 2Nx1 Double) = Current Source/Sink values
% C (Nx1 or 2Nx1 Integer) = Connection Vector (describing whether nodes are connected to is or vs(0/1))
% MemR (NxN Double) = Memristor array values
% LRowR (NxN Double) = Line Row Resistances array values
% LColR (NxN Double) = Line Column Resistances array values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% Circuit (Circuit) = Simulation circuit set-up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Circuit = fSpiceSim(N, varargin)
%% Set-up
[vs, is, connVec, MemR, LRowR, LColR] = fSetup(N, varargin);
id = fTimeStamp();

%Programs
mac_exec = '/Applications/MacSpice.app/Contents/MacOS/MacSpice';
win_exec = '..\Spice64\bin\ngspice_con';

%% Definition of different types of files
foldername = 'MacSpiceSimFiles/';
spiceFileName = sprintf('mat_array_%s.cir', id);
simOutputFile = sprintf('out_%s.txt', id);
tmpSimOutputFile = string(['tmp' simOutputFile]);

%% Generate Spice File
fileLoc = spiceFileName;
Circuit = fGenerateSpiceFile(N, vs, is, connVec, MemR, LRowR, LColR, fileLoc);

%% Perform Simulation Using Spice
try
    if ismac
        command = sprintf('%s -b %s > %s', mac_exec, spiceFileName, simOutputFile);
        [status] = system(command);
    elseif ispc
        command = sprintf('%s -b %s > %s', win_exec, spiceFileName, simOutputFile);
        [status] = system(command);
%     elseif isunix
    end
catch
    error("Spice Simulation could not be carried out. Ensure paths and relative paths are well defined")
end
%% Read Simulation Results
if status==0
    copyfile(char(simOutputFile), char(tmpSimOutputFile))
    Circuit = fReadSpiceSimResults(tmpSimOutputFile, Circuit);
else
    fprintf('Spice simulation could not run');
end

%% Delete Temporary Files
delete(char(tmpSimOutputFile))

%% Move files to defined destination folder
try
    if ~strcmp(foldername, '')
        movefile(char(spiceFileName), char([foldername, 'mat_array.cir']));
        movefile(char(simOutputFile), char([foldername, 'out.txt']));
    end
catch
    movefile(char(spiceFileName), char(foldername));
    movefile(char(simOutputFile), char(foldername));
end

%% Clear previous files that might exist in root folder
delete('mat_array*.cir')
delete('out*.txt')
delete('tmp*mat_array*.cir')
delete('tmp*out*.txt')

end

function [vs, is, connVec, MemR, LRowR, LColR] = fSetup(N, var)
%% Unpack Input Arguments
%Condition 1: vs, MemR, LRowR, LColR
if length(var)==4
    % Passed Arguments
    vs = var{1};
    MemR = var{2};
    LRowR = var{3};
    LColR = var{4};
    
    %Unpassed Arguments
    is = zeros(N, 1);
    connVec = zeros(N, 1);
end

%Default Condition
%Condition 2: s, connVec, MemR, LRowR, LColR
if length(var)==5
    % Passed Arguments
    vs = var{1};
    is = var{1};
    connVec = var{2};
    MemR = var{3};
    LRowR = var{4};
    LColR = var{5};
end

% %Condition 3: vs, is, connVec, MemR, LRowR, LColR
% if length(var)==6
%     vs = var{1};
%     is = var{2};
%     connVec = var{3};
%     MemR = var{4};
%     LRowR = var{5};
%     LColR = var{6};
% end

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
        zeros(((2*N)-lnght), 1)];
end

% Set Unconnected Sources to zero
vs = vs.*(connVec==0);
is = is.*(connVec==1);

end