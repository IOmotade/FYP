%% Add Relevant Paths
addpath(genpath('../FYP/Sim_Scripts/'))
addpath(genpath('../FYP/Functions/'))

%% Run Scripts
sim_N;
sim_numBitspRes;
sim_MemRMean;
sim_MemRRange;
sim_rule_lin_log;
sim_LRMean;
sim_LRVariance;
sim_SNR;

%%
exit