addpath(genpath('../FYP/Sim_Scripts/'))
addpath(genpath('../FYP/Functions/'))


%% Base Case for Algorithm 1
algo = [0 1 1.1 1.2];

setup.N = 8;
setup.numBitspRes = 4;
setup.ruleNum = 1;  %Memristor Encoding Rule
setup.var = {[1e3], [1e6]}; %Memristor upper and lower bounds for resistance values
setup.LRdef = [10 0 0]; %[mean, variance, variance of mean] of line resistance
