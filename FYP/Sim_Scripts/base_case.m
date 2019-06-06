addpath(genpath('../FYP/Sim_Scripts/'))
addpath(genpath('../FYP/Functions/'))

if algo_case==1
    %% Base Case for Algorithm 1
    algo = 1;
    
    setup.N = 8;
    setup.numBitspRes = 4;
    setup.ruleNum = 1;
    setup.var = {[100], [1e6]};
    setup.LRdef = [10e3 0 0];
    
    minNumBits = 100;
end