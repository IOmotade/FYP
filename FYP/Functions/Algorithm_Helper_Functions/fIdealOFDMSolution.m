% Omotade Iluromi, GROUP (EE4), 2019, Imperial College.
% 26/05/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the correct frequency spectrum matrix expected for the
% simulation setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% N (Integer) = Order of memristor array i.e. NxN array
% Rx, Ry (Double) = Resistance values with Rx and Ry in parallel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% Req (Double) = Equivalent resistance of circuit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [corr_fmag] = fIdealOFDMSolution(N, vs_mag, MemR, LRowR, LColR)
%% Generate Accurate Frequency Components
Vs = [vs_mag; zeros(N-1, 1)];
[corr_fmag, ~] = deal(zeros(N));
J = fShiftingMatrix(N);
for i = 1:N
    for j = 1:N
        Circuit = fMacSpiceSim(N, (J^(i-1))*Vs, MemR, LRowR, LColR);
        im = Circuit.IO.value(end, :);
        corr_fmag(i, j) = im(j);
        %             corr_fval(i, j) = fsource(i);
    end
    %         corr_Vbias(i, :) = Circuit.VI.value(i, :);
end
% fmag = corr_fmag;

end