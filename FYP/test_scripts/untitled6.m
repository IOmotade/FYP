%%
setup.basefreq = 5;
setup.perfect = false;
setup.oversampfactor = 1;%64*2;
tic
% BER = fBaseAlgorithmTest(setup, 4, 8, [100, 0], 1, 1e3, 1e6);
% [BER1, LREstimate1] = fAlgorithm1Test(setup, 4, 8, [100, 0], 1, 1e3, 1e6);
[BER2, LREstimate2] = fAlgorithm1_1Test(setup, 4, 2, [100, 0], 1, 1e3, 1e6);
% [BER3, LREstimate3] = fAlgorithm1_2Test(setup, 4, 8, [100, 0], 1, 1e3, 1e6);
% [BER4, LREstimate4] = fAlgorithm1_2(setup, 4, 16, [100, 0], 1, 1e3, 1e6);
% BER = fAlgorithm2Test(setup, 2, 1, [100, 0], 1, 10, 10e6);
toc

%%
setup.base_freq = 10;
setup.perfect = false;
setup.oversampfactor = 2;%64*2;
% [BER, LREstimate] = fAlgorithm2(setup, 4, 4, [1e3, 0], 1, 10, 100e3);
% [BER, LREstimate]
BER= fBaseAlgorithm(setup, 4, 8, [1e0, 0], 1, 10, 100e3);
% BER= fAlgorithm2(setup, 1, 3, [1e3, 0], 1, 10, 100e3);

% %%
% N = 4; numBitspRes = 4;
% memLength = N^2 * numBitspRes;
% idxLength = memLength^2;
% Fingerprint = zeros(idxLength, N, N);
% testBits = zeros(memLength, idxLength);
% for idx = 1:idxLength
%     bitsToStore = dec2bin(idx-1, memLength);
%     bitsToStore = bitsToStore(:);
%     testBits(:, idx) = str2num(bitsToStore);
% end

%%
% N = 4; numBitspRes = 4;
% memLength = N^2 * numBitspRes;
% idxLength = memLength^2;
% Fingerprint = zeros(idxLength, N, N);
% testBits = zeros(memLength, idxLength);
% disp("Generating Fingerprint Database");
% Fingerprint = zeros(idxLength, N, N);
% prog_txtlen = 0;
% for idx = 1:idxLength
%     prog_txt = sprintf("Index %d out of %d\n", idx, idxLength);
%     pause(0.1)
%     fprintf(strcat(repmat('\b', 1, prog_txtlen), prog_txt))
%     prog_txtlen = strlength(prog_txt);
% end
% disp("Generation Complete");


%%
rst
algo = [0 1 1.1 1.2];
minNumBits = 100;
setup.N = 4;
setup.numBitspRes = 4;
setup.var = {[1e3], [1e6]};
setup.LRdef = [100 0 0];
setup.SNR = 100;

setup.ruleNum = 1;
BER = fSimulation(algo, setup, minNumBits, false);
r1 = BER;

setup.ruleNum = 2;
BER = fSimulation(algo, setup, minNumBits, false);
r2 = BER;

setup.ruleNum = 3;
BER = fSimulation(algo, setup, minNumBits, false);
r3 = BER;

setup.ruleNum = 4;
BER = fSimulation(algo, setup, minNumBits, false);
r4 = BER;

[r1 r2 r3 r4]