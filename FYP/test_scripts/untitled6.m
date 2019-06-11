%%
setup.basefreq = 5;
setup.perfect = false;
setup.oversampfactor = 1;%64*2;
tic
% BER = fBaseAlgorithm(setup, 1, 1, [1e2, 0], 1, 10, 10e6);
[BER, LREstimate] = fAlgorithm1(setup, 4, 4, [100, 0], 1, 10, 10e6);
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
algo = 1;
minNumBits = 1000;
setup.N = 8;
setup.numBitspRes = 4;
setup.ruleNum = 1;
setup.var = {[1e3], [1e6]};
setup.LRdef = [100 0 0];
BER = fSimulation(algo, setup, minNumBits)