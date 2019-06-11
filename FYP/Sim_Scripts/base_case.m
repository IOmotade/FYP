addpath(genpath('../FYP/Sim_Scripts/'))
addpath(genpath('../FYP/Functions/'))

if algo_case==1
    %% Base Case for Algorithm 1
    algo = 1;
    
    setup.N = 8;
    setup.numBitspRes = 4;
    setup.ruleNum = 1;
    setup.var = {[1e3], [1e6]};
    setup.LRdef = [100 0 0];
    
    minNumBits = 1000;
end