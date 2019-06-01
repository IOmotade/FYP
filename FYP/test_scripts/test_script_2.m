%% Simple Circuit
% Testing simple Line Resistance Estimation
% rst
addpath(genpath('../FYP/Functions/'))
[fig_az, fig_el] = deal(145, 30);

%% Setup Circuit System
N = 3;
p = 1e-9;
Rp_estimate = zeros(size(p));

MemR  = abs(1e3*ones(N) + 0e3*randn(N));
LRowR = abs(1e4*ones(N)  + 0e3*randn(N));
LColR = abs(1e4*ones(N)  + 0e3*randn(N));

s_orig = [5, zeros(1, N-1)]';       %Source Magnitude
c  = zeros(2*N, 1);                 %Connection Matrix
% c(1:N-1) = 1;
%     c(N+1:end) = 1;
R = zeros(N);

idx = 2;
s = fShiftingMatrix(length(s_orig))^(idx-1) * s_orig;

%% Simulation
Circuit = fMacSpiceSim(N, s, c, MemR, LRowR, LColR);

%% Results & Calculation of Memristor values
% nCF = fNodeCurrentFactors([], N, 1, 1, 0, MemR, LRowR, LColR);
clc;
IO = Circuit.IO.value;
IMem = IO-Circuit.II.value;
im = IO(end, :);
CalMemR = zeros(N);
fprintf("----------\n")
for i = idx:idx%1:N
    for j = 1:N
        fprintf("Row %d, Col %d;" , i, j);
        %         [i, j] = deal(3, 3);
        nCF = fNodeCurrentFactors([], N, i, j, 0, MemR, LRowR, LColR)
        totalNCF = prod(nCF(i:end, 1));
        actualNCF = im(i)/IMem(i, j);
        fprintf("Total NCF: %f, Actual NCF: %f\n----------\n", totalNCF, actualNCF)
        actualIO = im(i)/totalNCF;
        %         fprintf("The actual current is: %s\n", fUnits(Circuit.IO.value(i, j), 'A'));
        %         fprintf("The measured current is: %s\n", fUnits(im(i), 'A'));
        %         fprintf("The corrected current is: %s\n", fUnits(actualIO, 'A'));
        CalMemR(i, j) = s(i)/actualIO - nCF(i,1)*fEquivalentResistance("Down", N, i, j, MemR, LRowR, LColR);
    end
end
fUnits(CalMemR, 'Ohm')

%%
% .48e-6/4.87e-3