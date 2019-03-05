%% Full-Simulation
function Circuit = fMacSpiceSim(N, vs, MemR, LRowR, LColR, scheme)
if N~=length(vs)
    fprintf('Dimensions do not match');
    return
end

if ~exist('scheme', 'var')
   scheme = '1'; 
end

%% Definition of different types of files
spiceFileName = 'mat_array.cir';
simOutputFile = 'out.txt';
tmpSimOutputFile = string(['tmp' simOutputFile]);

%% Generate Spice File
fileLoc = spiceFileName;
Circuit = fGenerateSpiceFile(N, vs, MemR, LRowR, LColR, fileLoc, scheme);

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