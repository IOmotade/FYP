addpath(genpath('../FYP/Functions/'))
internal_msg_len = 0;

%% Quick Analsis: Standard Deviation
sd  = 0:1:10;
BER = zeros(size(sd));
idxLength = length(BER);

%% Set-up
algo = 1;
minNumBits = 100;
setup.N = 1;
setup.numBitspRes = 1;
setup.ruleNum = 1;
setup.var = {[10], [100e3]};

%% Trial Script
disp('Starting Simulation');

tmp_prog_txtlen = 0;
for idx =1:idxLength
    setup.LRdef = [10e3 sd(idx)];
    BER(idx) = fSimulation(algo, setup, minNumBits);
    
    tmp_prog_txtlen = fClearInternalMessages(tmp_prog_txtlen);
    tmp_prog_txtlen = fDisplayInternalMessage(...
        sprintf('%2.2f Percent Complete', 100*(idx/idxLength)),...
        tmp_prog_txtlen);
end

disp('Simulation Complete!!!');

%% Plot Results
loglog(sd, BER);
semilogx(sd, BER)
save('sd_BER_tmp')

%% Clear all internal messages
% fClearInternalMessages(internal_msg_len);